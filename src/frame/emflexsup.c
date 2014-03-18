
/********************************************************************************
 *										*
 *				 M A C E K					*
 *		("Matroids Also Computed Efficiently" Kit)			*
 *										*
 * A set of tools and routines for computations with representable matroids.	*
 * Copyright (C) 2001,2002  Petr Hlineny.					*
 *										*
 *   This program is free software; you can redistribute it and/or modify	*
 *   it under the terms of the GNU General Public License as published by	*
 *   the Free Software Foundation; either version 2 of the License, or		*
 *   (at your option) any later version.					*
 *   You should have received a copy of the GNU General Public License		*
 *   along with this program; if not, write to the Free Software		*
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.	*
 * 										*
 * See the top-level "ABOUT" and "doc/ *" files for more information.		*
 *										*
 ********************************************************************************/













/**
 * These are supplementary functions directly related to the lexical scanner
 * in the file "emflex.l".
 * The file is NOT processed separately, but included in "emflex.c/l" instead!
 *            ~~~~~~~~~~~~~~~~~~~~~~~~~~
 * See more comments on the top of "emflex.l".
 * 
**/










/******************	Calling the analyzer	*******************/
/******************************************************************/


#undef CURDLEV
#define CURDLEV		6


/**
 * This function is called from outside to start scanning of the given string ins or list insl.
 * (It is the only public function of the analyzer.)
 * Before scanning starts, the character replacements requested by SPECREPLACE are applied.
 * Only one of ins, insl should be given.
 * The parameter name is an optional name for the first frame.
 * If name==NULL is given, then the first included filename is substituted.
 * The rec parameter is used internally to indicate that we call recursively for
 * the input list, it should not be used by outside calls.
 * 
 * The function returns the first scanned frame pointer (others are its subframes).
 * 
**/

static framestruc	*curframe = NULL;	/* (keeps the current frame while scanning) */
static ematrix		*savefrm = NULL;	/* (keeps the previous matrix when appending frame) */

framestruc*	frame_doinput_ext(char *ins, char **insl, framestruc *infr, char *name, int rec) {
	char		*repl[2] = SPECREPLACE, buf[50],*bnm, *inp, **xin;
	int		i;
	framestruc	*frs=NULL, *f;
	
	if ((ins && insl) || (!ins && !insl) || (insl?!insl[0]:0))
		{PROGERROREXIT("Exactly one of ins, insl must be given! %p %p",ins,insl);}
	if (!rec && curframe)
		{PROGERROREXIT("The lexical scanner is not re-entrant - current frame is %p [%s]",curframe,FRNAME(curframe));}
	lerror_occured = 0;	/* clear the local error-flag */
	
	if (insl) {
			/* scans the input list of strings one by one, indexing the given name */
		DEBUG(CURDLEV-1,"Calling to scan a list of input frames...\n");  SDEBUG(CURDLEV+1,"\t-----------\n");
		for (xin=insl,i=0; *xin; xin++,i++) {
			if (name) {
				snprintf(buf,48,name,i);  buf[48] = 0;  bnm = buf;
			} else  bnm = NULL;
			f = frame_doinput_ext(*xin,NULL,NULL,bnm,1);
			if (i==0)  frs = f;	/* (only the first frame is returned!) */
		}
	} else {
			/* scans one given string, applies the replacements requested by SPECREPLACE */
		DEBUG(CURDLEV-1,"Input frame - lexical scanning \"%s\":\n",ins);
		inp = yy_specrepl(ins,repl);
			/* here the new frame (or subframe) is created for the current input */
		yy_startframe(name,infr);
		yyincl(inp,0);
		fullineno = 1;  framelev = 0;
		frame_flex();		/* this calls the lexical analyzer preprocessed by flex */
		if (framelev>0)  {LERROR("Unclosed subframe(s) on the input!");}
		FREE(inp);
		frs = yy_endframe();	/* (the frame is finished and returned here) */
	}
	if (!rec)  curframe = NULL;	/* must clear the current frame reference for the next call */
	if (lerror_occured<0 && !infr)  return NULL;	/* nothing is returned on an error (and no cleaning here) */
	return frs;
}


/**
 * This function is called to read input matrix(ces) from the string (or file name) ins.
 * The commands given with the matrix are processed before giving the matrix.
 * If one==1, then exactly one frame with one matrix is expected, otherwise an error is issued.
 * The new matrix(ces) are returned in a list (with names inherited from the frames).
 * Do not forget to free the list later.
**/

ematrix**	frame_inputmatrix_ext(char *ins, int one) {
	framestruc	*fr, **fwa, **x;
	ematrix		*e, **el=NULL;
	
	printlev -= 2;		/* (to lower debug printing inside here) */
	fr = frame_doinput(ins);
	if (fr)  frame_processcommands(fr);
	printlev += 2;
	if (!fr)  return NULL;
	
	fwa = frame_gettree_dep(fr);
	for (x=fwa; x?*x:0; x++) {
		if (FRMATRIX(*x) && FRPFINDEX(*x)!=pfield_curindex())
			{USERERROR("The frame %p on input has matrix over a different pfield %d!=%d.",*x,FRPFINDEX(*x),pfield_curindex()); break;}
		e = FRMATRIX(*x)? ematrix_copy(FRMATRIX(*x)):NULL;
		if (e && FRNAME(*x))  EMSETNAME(e,FRNAME(*x));    
		if (e)  el = alist_append(el,e);  
	}
	if (one && alist_getlength(el)!=1)  {USERERROR("Exactly one matrix is expected from \'%s\'.",ins);}
	alist_free(fwa);
	dispose_frame_recur(fr);
	return el;
}



/**
 * These functions handle creation and finishing of frames (and subframes) during
 * scanning the input:
 * yy_startframe() starts a new frame as a son of the current frame, with an optional name.
 * (If no name is given, then the first included filename is used, or a variation of
 * the parent's name.)
 * Moreover, one may optionally give a frame to append all input to (use carefully!).
 * yy_endframe() ends the current frame - the matrix is reduced to the proper size,
 * and the parent becomes the next current frame (it stays the same if no parent).
 * 
 * These functions are directly called from the lexical analyzer.
 * The current frame pointer is stored in a static variable curframe here, and it may be
 * accessed only within the functions in this file.
**/

void	yy_startframe(char *name, framestruc *frst) {
	ematrix	*ee;
	
	if (frst && curframe)  {PROGERROR("frst may be given only at the start of scanning!");}
	errhint = -1;
	if (frst)  DEBUG(CURDLEV,"%*s# Appending input frame \"%s\" (fr-lev %d) ...\n",framelev,"",FRNAME(frst),framelev);
	else  DEBUG(CURDLEV,"%*s# New input frame \"%s\" (fr-lev %d) starting...\n",framelev,"",name?name:"?",framelev);
	
	if (frst)  curframe = frst;		/* the given frame frst is used for input */
	else curframe = new_frame(curframe);	/* a new frame structure is created here */
	if (name)  frame_setname(curframe,name);
	if (!FRMATRIX(curframe)) {		/* a new matrix for the frame (stores also cur. pfield) */
		frame_setmatrix(curframe, new_ematrix(16,16,16,16));
		frst = NULL;	/* (not to be used as current matrix dimensions next) */
		savefrm = NULL;
	} else {
		savefrm = ematrix_copy(FRMATRIX(curframe));
	}			/* (we may set the previous matrix back to the frame on error) */
	if (frst) {
		ee = FRMATRIX(frst);		/* the current (real) matrix dimensions */
		FRMATRIXROWS(curframe) = mxline = ROWSM(ee);
		FRMATRIXCOLS(curframe) = COLSM(ee);
		ROWSM(ee) = ROWSMAX(ee);  COLSM(ee) = COLSMAX(ee);
	} else {
		mxline = FRMATRIXROWS(curframe) = FRMATRIXCOLS(curframe) = 0;
	}
}

framestruc*	yy_endframe(void) {
	framestruc	*fr;
	ematrix		*ee,*en;
	char		*nm, buf[30], **noval;
	int		sz;
	
	ee = FRMATRIX(curframe);
	ROWSM(ee) = FRMATRIXROWS(curframe);	/* the actual dimensions of the matrix */
	COLSM(ee) = FRMATRIXCOLS(curframe);
	sz = ROWSM(ee)+COLSM(ee);
	if (lerror_occured<0 && savefrm) {	/* replace the old matrix on syntax errors */
		en = savefrm;
	} else if (sz>0) {			/* the proper-size copy of the matrix, with tr==0 */
		en = ematrix_copy_ext(ee,ROWSM(ee)+1,COLSM(ee)+1,0);
		ematrix_resetid(en);
	} else  en = NULL;
	frame_setmatrix(curframe,en);
						/* fills in some name if not determined yet */
	if (!FRNAME(curframe) && FRPARENT(curframe)) {
		nm = FRNAME(FRPARENT(curframe));
		if (nm)  snprintf(buf,25,FR_DEFSUBNAME,FRNUMSONS(FRPARENT(curframe)));
		else  snprintf(buf,25,FR_DEFSUBNONAME,framelev,FRNUMSONS(FRPARENT(curframe)));
		frame_setname2(curframe,nm,buf);
	}
	if (!FRNAME(curframe))  frame_setname(curframe,FR_DEFNONAME);
	noval = NULL;
#ifndef FASTPROG
	noval = frame_getoptionval_last(curframe,"nopfcheck");  if (noval)  dispose_alist(noval);
	if (noval==NULL)
		if (sz>2 && sz<20) if (ematrix_inpfield_rand(FRMATRIX(curframe))<0) {
			EMATOUTPUT(FRMATRIX(curframe),"!\t");
			ematrix_inpfield_printed(FRMATRIX(curframe));
			LERROR("The matrix is not represented over the pfield!");
		}
	DEBUG(CURDLEV,"%*s=> Input of frame %p [%s] (fr-lev %d) finished.\n",framelev,"",curframe,FRNAME(curframe),framelev);
	if (CURDLEV+1<=DEBUGLEV)  frame_fprint_all(debugout,curframe,"\t");
	SDEBUG(CURDLEV+1,"\t-----------\n");
#endif
	if (FRPFINDEX_SAVE(curframe)>=0 && FRPFINDEX_SAVE(curframe)!=pfield_curindex()) {
			/* back to the original pfield of this frame - <0 or stored at '@inputpf XX' */
		pfield_switchto_fast(FRPFINDEX_SAVE(curframe));
		DEBUG(CURDLEV-2,"Switched input arithmetics back to pfield \"%s\"\n",pfield_curname());
	}
	fr = curframe;
	if (FRPARENT(curframe))  curframe = FRPARENT(curframe);	/* back to the parent frame */
	errhint = -1;
	mxline = FRMATRIXROWS(curframe);	/* (stored number of matrix lines for the parent frame) */
	return fr;
}














/******************	Modification functions	*******************/
/******************************************************************/


#undef CURDLEV
#define CURDLEV		5


/**
 * This function performs the special line-input replacements given by repl[] --
 * characters of repl[0] are changed to corresponding characters of repl[1].
 * The replacement is not applied inside ".."-enclosed sections.
 * Morever, the leading and ending "s are stripped, and the internal \" are
 * replaced by single ", \..\"s are replaced by one less of \s.
 * The function returns back a copy of the modified string.
**/

char*	yy_specrepl(char *st, char *repl[]) {
	char	*inp;
	int	i,j,l, q;
	
	if (repl)  DEBUG(CURDLEV+1,"Line-replacing characters (\"%s\" -> ...) in \'%s\'\n",repl?repl[0]:"",st);
	
	l = strlen(st);		/* stripping the leading and ending "s here: */
	if (st[0]=='"' && st[l-1]=='"') {
		inp = MSTRDUP(st+1);  inp[l-2] = 0;
	} else  inp = MSTRDUP(st);
				/* handling \..\" here (one less active backslash): */
	for (i=j=0; inp[i]; i++,j++) {
		if (inp[i]=='\\' && inp[i+1]=='"')  i++;
		inp[j] = inp[i];
	}
	inp[j] = 0;
				/* replacing characters required by repl[] here: */
	if (repl)  for (i=0; repl[0][i]; i++) {
		if (repl[0][i]=='"' || repl[1][i]=='"')  {PROGERROR("Do not mix quotas in the replace string.");}
		for (j=q=0; inp[j]; j++) {
			if (inp[j]=='"' && (j<=0 || inp[j-1]!='\\'))  q = !q;
			if (q==0 && repl[0][i]==inp[j])  inp[j] = repl[1][i];
		}
		if (i==0 && q>0)  {LERROR("Unmatched double-quotas in the input string? \'%s\'",st);}
	}
	if (repl)  DEBUG(CURDLEV+2,"Line-replaced \'%s\' -> \'%s\'\n",st,inp);
	return inp;
}


/**
 * This function modifies the given string st in the following way:
 * The function yyb_processvalue() (in a separate supplementary lexical scanner emflexb.l)
 * is called to process possible "$"-macros in the value - text-replacing them.
 * It returns back a copy of the modified string.
**/

char*	yy_ovaluesubst(char *st) {
	char	*sr;
	
	if (!st)  return NULL;
	sr = yyb_processvalue(curframe,st);	/* (call "yyb_"-scanner from emflexb.l) */
	return sr;
}


/**
 * This function is called from the lexical analyzer to handle requested symbol
 * replacements in the input expressions.
 * The current string is given in st, and the replacement string is returned
 * via the pointer *rep.
 * The function returns how long prefix of st was replaced (or 0).
 * 
 * The requested replacements are obtained from the option "replace" stored
 * in the current frame, or from "repl-PF" for the current pfield PF.
 * (The options are read only if something has changed in the frame or the options.)
 * 
 * The string value returned via *rep must not be changed !!!
 * 
**/

int	yy_entrysubst(char *st, int wh, char **rep) {
	static framestruc	*lcurf = NULL;
	static int		loi = -1;
	static char		*lpf = NULL;
	static char		**rlist[2] = {NULL,NULL};
	char		buf[50], **x;
	int		l=0;
	
	if (wh<0 || wh>=1)  {PROGERROR("Wrong value wh=%d",wh); return 0;}
	DEBUG(CURDLEV+2,"calling with st=%s, wh=%d\n",st,wh);
				/* read the options again only if they might have changed... */
	if (lcurf!=curframe || loi!=FRLASTOPTI(curframe) || lpf!=pfield_curname() || !rlist) {
		lcurf = curframe;  loi = FRLASTOPTI(curframe);  lpf = pfield_curname();
		if (rlist[0])  dispose_alist(rlist[0]);
		rlist[0] = frame_getoptionval_exval(curframe,"replace",2);
		sprintf(buf,"repl-%.20s",lpf);
		x = frame_getoptionval_exval(curframe,buf,2);
		rlist[0] = alist_applist(rlist[0],x);
	}
	*rep = NULL;		/* search for prefix match in the replacement values */
	if (rlist[wh]) for (x=rlist[wh]; x[0] && x[1]; x+=2)
		if (strncmp(st,x[0],strlen(x[0]))==0) {
			*rep = x[1];  l = strlen(x[0]);
	}
	return *rep? l:0;
}



/**
 * These two functions handle the stack of include files, and also properly preserve
 * the (debugging) filenames and linenumbers.
 * (This is not the place where subframes are handled!)
 * If a subframe is to be included, it is already created and this file is just included
 * into the subframe.
 * If fil==0 is given in yyincl(), then the given string is included itself, not as
 * a filename.
 * Filenames are searched according to the general frame-reading path.
 * The return value of yywrap() is defined by flex, the return value of yyincl()
 * indicates a success (0) or an error (-1).
**/

#define	MAX_INCLUDE_DEPTH	300
static YY_BUFFER_STATE	include_stack[MAX_INCLUDE_DEPTH];
static FILE		*include_stack_file[MAX_INCLUDE_DEPTH+1];
static char		*include_stack_fn[MAX_INCLUDE_DEPTH];
static int		include_stack_ln[MAX_INCLUDE_DEPTH];
//static int		include_stack_ptr = -1;		moved to emflex.l above!!!
static int		include_app = 0;


int	yywrap() {
	
	if (!include_app) {	/* all input files and strings get a newline appended to them */
		include_app = 1;
		yy_scan_string("\n");
		return 0;
	} else {
		include_app = 0;
		if (inpfile) {	/* (this is just an informal input text for printing-out...) */
			FREE(inpfile);  inpfile = NULL;
		}
			/* delete the input buffer, and possibly close the input file */
		yy_delete_buffer(YY_CURRENT_BUFFER);
#if CLOSE_AFTER_FLEX>0
		if (include_stack_file[include_stack_ptr])
			fclose(include_stack_file[include_stack_ptr]);
#endif
			/* pop the previous saved state from the stack */
		if (--include_stack_ptr>=0) {
			yy_switch_to_buffer(include_stack[include_stack_ptr]);
			inpfile = include_stack_fn[include_stack_ptr];
			inpline = include_stack_ln[include_stack_ptr];
			yy_pop_state();
			return 0;
		} /*** else  yyterminate();	does not work ??? ***/
	}
	return 1;	/* (means the end of current input) */
}


int	yyincl(char *fn, int fil) {
	char		*s;

	DEBUG(CURDLEV+1,"Going to include (ptr=%d) %s \'%s\'.\n",include_stack_ptr,fil?"file":"string",fn);
	if (include_stack_ptr>=MAX_INCLUDE_DEPTH || include_stack_ptr<-MAX_INCLUDE_DEPTH) {
		PROGERROR("Includes are nested too deeply! %d",include_stack_ptr);
		exit(-1);
	}
	if (include_stack_ptr<0) {
		include_stack_ptr = -1;	 include_app = 0;
		inpfile = NULL;  include_stack_file[0] = NULL;
	}
			/* opens the given include file for reading */
	if (fil) {
		if (include_stack_ptr>0)	/* (tries also the subdir of the previous include) */
			s = fdir_extract(include_stack_fn[include_stack_ptr-1]);
		else  s = NULL;
		yyin = fopen_path(fn,"r",frame_path_read,s,frame_extension);
		//********* add stdin as an input "-" (interactive input) ........
		if (!yyin) {
			LERROR("Cannot open include file \'%s\'!",fn);
			return -1;
		}
		if (!FRNAME(curframe)) {
			frame_setname(curframe,s=fname_extract(fn));
			FREE(s);	/* (this s is a private copy of the string) */
		}
	}
	if (include_stack_ptr>=0) {	/* stores the current buffer, and opens a new buffer */
		include_stack[include_stack_ptr] = YY_CURRENT_BUFFER;
	}
	if (fil) {
		yy_switch_to_buffer(yy_create_buffer(yyin,YY_BUF_SIZE));
		if (include_stack_ptr+1>=0)  include_stack_file[include_stack_ptr+1] = yyin;
	} else {
		yy_scan_string(fn);
		if (include_stack_ptr+1>=0)  include_stack_file[include_stack_ptr+1] = NULL;
	}
	if (include_stack_ptr>=0) {	/* stores the whole current state on the stack */
		yy_push_state(INITIAL);
		include_stack_fn[include_stack_ptr] = inpfile;
		include_stack_ln[include_stack_ptr] = inpline;
	}
	inpfile = MSTRDUP(fn);	/* (freed in yywrap()) */
	inpline = 1;
	include_stack_ptr++;
	inpspecial = !fil;	/* whether special-lines should apply to this include - see SPECMATRIXLINE */
	return 0;
}















/******************	Evaluating expressions	*******************/
/******************************************************************/


#undef CURDLEV
#define CURDLEV		5	/*5*/


/**
 * This handles an arithmetic-expression stack used when evaluating one complex matrix entry.
 * In general, we allow an arbitrary symbolic expression (defined in the pfield) to be
 * on the program input.
 * Since we need to handle bracketing and usual operator precedence, we need to store
 * the atoms and operators of the expression on this stack.
 * 
 * The function yy_entryop() is called either for an expression atom - given as a raw
 * string in nst (with op==YYOP_NOP),
 * or for an operator - given in op by a value from enum YY_ENTRYOPS.
 * Additionally, the operator YYOP_POWER gets its integral exponent in pow.
 * (Operators include brackets.)
 * 
 * How the stack evaluation works:
 *  - Given atoms are stored on the stack.
 *  - If an operator arrives, all immediately previous operations with the same or higher
 *    priority are evaluated, except a possible left bracket.
 *  - A power operator is always immediately applied to the last atom.
 *  - A left bracket is always simply stored on the stack.
 *  - A right bracket evaluates the whole preceeding expression to its matching left bracket.
 *    (Like if it has the lowest priority...)
 *  - If no operator is given between atoms, then multiplication is used.
 *  - If no operator preceeds a left bracket, then multiplication is used, again.
 *  - If no atom is given before an operator, then 0 is used (+ -) or 1 is used (* /).
 *    (This correctly evaluates expressions like (-x) or (/x).
 * 
 * The return value for an atom entry is the scanned length of the atom (from pfield_scanvalue_one).
 * The return value is 0 otherwise.
 * See also yy_entryeval() below.
**/

#define MAX_OP_DEPTH	200
static struct {
	exp_t	x;
	sign_t	g;
	int	op;
}			op_stack[MAX_OP_DEPTH];
static int		op_stack_ptr = -1,
			op_error = 0;

int	yy_entryop(char *nst, int op, int pow) {
	int	k=0, am,sx,sg;
	exp_t	xx;
	sign_t	gg;
	
	if (op_stack_ptr<0) {
		op_stack[op_stack_ptr = 0].g = -11;
		op_error = 0;
	}
	if ((op==YYOP_LBRACE || op==YYOP_NOP) && op_stack[op_stack_ptr].g>=-1) {
		yy_entryop(NULL,YYOP_TIMES,0);
	}
	DEBUG(CURDLEV+2,"  %*s. op stack called:  nst=\"%s\", op=%d, pow=%d\n",op_stack_ptr,"",nst?nst:"",op,pow);
		/**
		 * An atom is given here - stored on the stack.
		 * (If a second atom is to be stored, a YYOP_TIMES operation is inserted.)
		 * The return value is copied from pfield_scanvalue_one().
		**/
	if (op==YYOP_NOP) {
		k = pfield_scanvalue_one(1,nst,&xx,&gg);
		if (k<=0)  return 0;
		op_stack[op_stack_ptr].x = xx;
		op_stack[op_stack_ptr].g = gg;
		DEBUG(CURDLEV+2,"  %*s(%d) .  stored atom %s len %d\n",op_stack_ptr,"",op_stack_ptr,pfield_pvalue(9,xx,gg),k);
		return k;
	}
		/**
		 * An operation without an atom - puts 0 or 1 in front of it.
		**/
	if (op_stack[op_stack_ptr].g<-1) {
		pfield_setzeroexp(&op_stack[op_stack_ptr].x);
		op_stack[op_stack_ptr].g = (op==YYOP_TIMES || op==YYOP_OVER)? 1:0;
	}
		/**
		 * This handles the given operations (YYOP_POWER is applied immediately):
		 * All preceeding aperations of the same or higher priority on the stack
		 * are evaluated first (up to a left bracket).
		 * Brackets are handled specially...
		**/
	if (op==YYOP_POWER) {
		pfield_power(pow,op_stack[op_stack_ptr].x,op_stack[op_stack_ptr].g,
				&op_stack[op_stack_ptr].x,&op_stack[op_stack_ptr].g);
	} else {
		while (op_stack_ptr>0) {
			if (op==YYOP_NOP)  break;
			if (op>op_stack[op_stack_ptr-1].op)  break;
			if (op>YYOP_RBRACE && op_stack[op_stack_ptr-1].op==YYOP_LBRACE)  break;
			
			sx = sg = 1;  am = -1;
			switch (op_stack[op_stack_ptr-1].op) {
			
			case YYOP_PLUS:		am = 0;
					break;
			case YYOP_MINUS:	am = 0;  sg = -1;
					break;
			case YYOP_TIMES:	am = 1;
					break;
			case YYOP_OVER:		am = 1;  sx = -1;
					break;
			case YYOP_LBRACE:	am = 3;
						if (op!=YYOP_EVAL)  op = YYOP_NOP;
					break;
			default:	PROGERROREXIT("Corrupted op stack!?!");
			}
				/* this part evaluates the prepared operation */
			if (am==3) {
				op_stack[op_stack_ptr-1].x = op_stack[op_stack_ptr].x;
				op_stack[op_stack_ptr-1].g = op_stack[op_stack_ptr].g;
			} else {
				k = pfield_arithmetic_ext(am,0,1,op_stack[op_stack_ptr-1].x,op_stack[op_stack_ptr-1].g,
					sx,op_stack[op_stack_ptr].x,sg*op_stack[op_stack_ptr].g,
					&op_stack[op_stack_ptr-1].x,&op_stack[op_stack_ptr-1].g);
				if (k<0) {
					op_error = 1;
					/****if (YY_START==Imxnumbers)  {LERROR("Undefined operation in the matrix expression!");}*/
				}
			}
			op_stack_ptr--;
		}
		/**
		 * Finally, after the preceeding higher operation were evaluated,
		 * the current operation is stored on the stack, and the function returns.
		**/
		op_stack[op_stack_ptr].op = op;
		DEBUG(CURDLEV+2,"  %*s(%d) .  stored operation %d\n",op_stack_ptr,"",op_stack_ptr,op);
		if (op!=YYOP_NOP) {
			if (op_stack_ptr>=MAX_OP_DEPTH-1) {
				PROGERROR("Expression nested too deeply!");
				op_stack_ptr /= 2;  op_error = 1;
			}
			op_stack_ptr++;
			op_stack[op_stack_ptr].g = -11;
		}
	}
	return op_error?-1: 0;
}



/**
 * Here we have functions called to evaluate a whole entry expression (already fed in by
 * the function yy_entryop()).
 * The first function yy_entryevalcheck() actually does the evaluation, checks for errors
 * (return -1 if an error occured), and stores the result into curframe->xx,gg.
 * The second function yy_mxentryeval() actually stores the previously evaluated entry
 * to the matrix in curframe, and possibly enlarges the matrix if neccessary.
 * The function yy_mxentryline() finishes the whole matrix line.
 * 
 * ....
**/

int	yy_entryevalcheck(int cd) {
	
	if (cd!=Imxnumbers && cd!=Ivalnumbers)  {PROGERROREXIT("Invalid condition occured! %d",cd);}
	if (op_stack[op_stack_ptr].g<-1)  {LERROR("Probably unfinished expression..."); return -1;}
	
	if (op_stack_ptr>0)  yy_entryop(NULL,YYOP_EVAL,0);
	if (op_error || op_stack[0].g<-1) {
		/***if (cd==Imxnumbers)  {LERROR("Computation error or no value to evaluate?!");}*/
		op_stack[0].g = -111;
	}
	op_stack_ptr = -1;	/* clear the stack for the next expression */
	if (curframe) {
		FRVALUEEXP(curframe) = op_stack[0].x;  FRVALUESIG(curframe) = op_stack[0].g;
	}			/* stores the evaluated value into the current frame */
	return op_stack[0].g<-1? -1:1;
}


int	yy_mxentryeval(int r, int c) {
	ematrix	*ee;
	
	if (!curframe || !FRMATRIX(curframe))  {PROGERROREXIT("The frame of the matrix must be initialized here!");}
	if (FRVALUESIG(curframe)<-1)  {PROGERROREXIT("Undefined matrix entry over %s to store ?!",pfield_curname());}
		
	ee = FRMATRIX(curframe);	/* automatically enlarging the matrix if neccessary */
	if (c>=COLSMAX(ee)-1 || r>=ROWSMAX(ee)-1) {
		ee = ematrix_copy_to(ee,2*ROWSMAX(ee)+2,2*COLSMAX(ee)+2);
		ROWSM(ee) = ROWSMAX(ee);  COLSM(ee) = COLSMAX(ee);
		frame_setmatrix(curframe,ee);
	}
			/* setting the matrix entry by the result evaluated in yy_entryevalcheck() */
	SETEXSIGM(ee,r,c,FRVALUEEXP(curframe),FRVALUESIG(curframe));
	if (c+1>FRMATRIXCOLS(curframe))  FRMATRIXCOLS(curframe) = c+1;
	//DEBUG(CURDLEV-2,"  - col %d, max %d\n",c,FRMATRIXCOLS(curframe));
	//*************** some error when a second line is appended to existing frame!!!!!!!!
	
	DEBUG(CURDLEV+2,"  - matrix entry (%d,%d): %s\n",r,c,pfield_printvalue_ext(NULL,25,op_stack[0].x,op_stack[0].g,0));
	return 0;
}

int	yy_mxentryline(int r, int c, int ln) {
	
	FRMATRIXROWS(curframe) = r+1;	/* just remembers the number of scanned matrix lines */
	DEBUG(CURDLEV+2,"== close matrix line (%d) at %d\n",r,ln);
	return c=ln=0;
}


/**
 * This function is called to close an evaluation of an input expression other than
 * the matrix entries.
 * There is no functionality here now, just storing the result for later retrieval...
**/

static exp_t	valxxe;
static sign_t	valgge;

int	yy_valentryeval(int cd) {
	
	valxxe = FRVALUEEXP(curframe);
	valgge = FRVALUESIG(curframe);
	DEBUG(CURDLEV+1," - Value param  (%d):  %s\n",cd,valgge>=-1?pfield_pvalue(25,valxxe,valgge):"undefined");
	return cd=0;
}

void	yy_valentrylast(exp_t *x, sign_t *g) {
	if (x)  *x = valxxe;
	if (g)  *g = valgge;
}




















/******************	Handling options and commands	*******************/
/**************************************************************************/


#undef CURDLEV
#define CURDLEV		5


/**
 * These functions are called to start and to end recording of one input option/command.
 * The start function gets the option name, the current lex condition, the logical line
 * of this option, and a "continuation flag" (see below).
 * The end function gets the total number of values, and the logical line again.
 * The option values (command parameters) are stored by the next two functions...
 * 
 * The current option pointer is kept in curoption just below.
 * When an option is ended, it is given to yy_applyoption() for possible immediate
 * application in the lexical scanner.
 * When a command is ended, it is only stored in the current frame.
 * However, when a command-procedure is ended, it is disposed, since it was previously
 * text-translated into include of a file of other commands (with "@sub-xx" parameters).
 * 
**/

static optionstruc	*curoption=NULL;	/* the current option pointer */
static int		curoptvalix=0;
static int		curoptcd=0;		/* the current option condition from lex */

void	yy_optionstart(char *nm, int cd, int ln, int ct) {
					/* (ct==1 starts a "continuing" option, not from lex) */
	
	if (curoption)  {PROGERROR("A current option is already set %p!",curoption);}
	if (nm?!nm[0]:1)  {PROGERROR("No name is given for the new option!"); nm="";}
	if (!ct)  DEBUG(CURDLEV+2,"A new %s \"%s\" is created (ln %d)...\n",cd!=Ioption?"command":"option",nm,ln);
	
	curoption = new_comoption_ext(nm,(cd!=Ioption),ct);
	OPTINDEX(curoption) = ln;
	curoptvalix = 0;  curoptcd = cd;
}

void	yy_optionend(int ix, int ln) {
	optionstruc	*curo;
	
	if (ln>=0 && ln!=OPTINDEX(curoption))  {PROGERROR("Option ends with different line number ?!? %d!=%d",ln,OPTINDEX(curoption));}
	if (ix>=0 && curoptcd!=Icomproc)  DEBUG(CURDLEV+2," -  %s \"%s\" is stored (ln %d) [ %s %s .. ]\n",
			curoptcd!=Ioption?"command":"option",OPTNAME(curoption),ln,OPTPARAM(curoption,0),OPTPARAM(curoption,1));
	if (ix>=0 && curoptcd==Icomproc)  DEBUG(CURDLEV+1," -  getting procedure %s (ln %d) [ %s %s .. ]\n",
			OPTNAME(curoption),ln,OPTPARAM(curoption,0),OPTPARAM(curoption,1));
	
	curo = curoption;  curoption = NULL;
	if (curoptcd==Icommand) {
		frame_addcommand(curframe,curo,ln);
	} else if (curoptcd==Ioption) {
		frame_addoption(curframe,curo,ln);
		yy_applyoption(curo,ln);
	} else if (curoptcd==Icomproc) {
		dispose_comoption(curo);
	}
}

/**
 * Supplements for handling "&procedure"s with their "@sub-paramI xxx" parameters...
 * 
**/

optionstruc*	yy_optionproc(void) {
	return (curoptcd!=Icomproc?NULL: curoption);
}

void	yy_optionprocpar(int ix, char *val) {
	char		buf[50];
	optionstruc	*sub;
	
	if ((val?val[0]:0)==FRGET_RESULT) {
		strncpy(buf,FR_SUBPARAMRES,45);  buf[45] = 0;
		val++;
	} else {
		strncpy(buf,FR_SUBPARAMNAME,45);  buf[45] = 0;
		buf[strlen(buf)-1] = '1'+ix;
	}
	sub = new_option_onestr(buf,val?val:"",OPTNOPARAM);
	frame_addoption(curframe,sub,fullineno);
	DEBUG(CURDLEV+2,"   -  proc parameter \'%s\' \'%s\'\n",buf,val);
}


/**
 * These functions are called to store one option value to the current option.
 * The parameters are the index of the value (starting with 0), the text value,
 * and the numeric integer value (or OPTNOPARAM).
 * Optionally, a pfield value may be stored for an option, but only as the first value.
 * 
 * The number of possible values for one option structure is fixed, and so a subsequent
 * given value forces us to create a "continuing" copy of the option (with the same name,
 * but marked as a continuation) that keeps more values, and so on...
 * 
 * ......
 * 
**/

void	yy_optionval(int ix, char *val, long dv) {
	char	*s;
	int	l;
	
	if (!curoption || ix<curoptvalix)  {PROGERROREXIT("No option has started %p, or invalid index %d.",curoption,ix);}
	if (!val && dv==OPTNOPARAM)  return;
	
	if (ix-curoptvalix>=MAXOPTPARAMS) {
		s = OPTNAME(curoption);  l = OPTINDEX(curoption);
		yy_optionend(-1,-1);
		yy_optionstart(s,curoptcd,l,1);	/* a "continuing" option to accommodate all values */
		curoptvalix += MAXOPTPARAMS;
	}
	option_setparam(curoption,ix-curoptvalix,val,dv);	/* setting the value */
}

void	yy_optionpfval(int ix, exp_t x, sign_t g) {
	
	if (ix!=0)  {PROGERROR("Only the first value may be a pfield element! %d",ix);}
	OPTPAREX(curoption) = x;
	OPTPARSG(curoption) = g;
}



/******************	Immediate option application	*******************/
/**************************************************************************/


#undef CURDLEV
#define CURDLEV		5


/**
 * This function is called after each recorder input option.
 * If the option op should take immediate effect (like @require ...), then the appropriate code
 * is executed here.
 * The function so far does not work with continued options!
 * (Consider not writing such immediate options... fropts.c)
 * 
**/

void	yy_applyoption(optionstruc *op, int ln) {
	int	r,x;
	ematrix	*ex;
	
	if (OPTISOPTION(op) && OPTISNAME(op,"require")) {
		if (OPTPARSG(op)<-1)  {LERROR("Undefined required value over %s at l. %d!",pfield_curname(),ln);}
		if (OPTPARNUM(op)>=2 && OPTNPARAM(op,1)<10)
			if ((OPTPARSG(op)==0)!=(OPTNPARAM(op,1)==0))
				{LERROR("Wrong required value over %s at l. %d:  %s ! %s !",
					pfield_curname(),ln,pfield_pvalue(20,OPTPAREX(op),OPTPARSG(op)),OPTNPARAM(op,1)==0?"==0":"!=0");}
		if (OPTNPARAM(op,1)<0 || (OPTNPARAM(op,1)>1 && OPTNPARAM(op,1)!=OPTNOPARAM))
				{LERROR("Expecting no value, or 0, or 1 for \"@require expr X\"!");}
		
	} else if (OPTISOPTION(op) && OPTISNAME(op,"name")) {
		if (OPTPARAM(op,0))  frame_setname(curframe,OPTPARAM(op,0));
	
	} else if (OPTISOPTION(op) && OPTISNAME(op,"comment")) {
		if (OPTPARAM(op,0))  FRCOMMENT(curframe) = MSTRDUP(OPTPARAM(op,0));
	
	} else if (OPTISOPTION(op) && OPTISNAME(op,"inputpf")) {
		x = pfield_curindex();
		r = pfield_switchto_ext(0,OPTPARAM(op,0));
		if (r<0)  {LERROR("Failed to switch input pfield to \"%s\", now in \"%s\"\n",OPTPARAM(op,0),pfield_curname());}
		else  DEBUG(CURDLEV-2,"Switched input arithmetics to a new pfield \"%s\"\n",pfield_curname());
		if (mxline>0 && r>=0)  {LERROR("Do not switch input pfield inside a matrix!");}
		else  FRPFINDEX(curframe) = pfield_curindex();
		if (r>=0 && FRPFINDEX_SAVE(curframe)<0)  FRPFINDEX_SAVE(curframe) = x;
	
	} else if (OPTISOPTION(op) && OPTISNAME(op,"transpose")) {
		if ((ex=FRMATRIX(curframe))!=NULL) {
			ematrix_transpose(ex);
			mxline = FRMATRIXCOLS(curframe);
			FRMATRIXCOLS(curframe) = FRMATRIXROWS(curframe);
			FRMATRIXROWS(curframe) = mxline;
		}
	
	//} else if (OPTISOPTION(op) && OPTISNAME(op,"")) {
		
	} else  return;
	
	DEBUG(CURDLEV+1,"Option \"%s\" applied (l. %d).\n",OPTNAME(op),ln);
	if (OPTCONTINUED(op) || curoptvalix>=MAXOPTPARAMS)
		DEBUG(CURDLEV-2,"** Warning ** - a continued option \"%s\" is not correctly processed here!\n",OPTNAME(op));
}







































