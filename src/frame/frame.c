
/********************************************************************************
 *										*
 *				 M A C E K					*
 *		("Matroids Also Computed Efficiently" Kit)			*
 *										*
 * A set of tools and routines for computations with representable matroids.	*
 * Copyright (C) 2001--2011  Petr Hlineny.					*
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
 * A "frame" is the basic general data entity in the program.
 * Read more about a frame and its data in ../include/frame.h.
 * 
 * This file contains basic handling functions for frame options/commands, and for
 * the frame itself.
 * These include creating, disposing, printing, and basic data-setting...
 * 
 * 
**/





#include "macek.h"  
#include "fr.h"









/******************	Basic option handling	*******************/
/******************************************************************/


#undef CURDLEV
#define CURDLEV         6

static int	fr_aloptions = 0, fr_alframes = 0;


/**
 * These functions are provided for creating and disposing option/command structures
 * (oc distinguishes option 0 or command 1, cont is for the "continued" flag).
 * The name nm must be given and it is stored as a private copy of the given string.
 * A newly created option/command has no values/parameters.
 * A copy of an option/command gets copies of the name and of all the values.
**/

optionstruc*	new_comoption_ext(char *nm, int oc, int cont) {
	optionstruc	*op;
	int		i;
	
	if (!nm)  {PROGERROR("The option name must be given!");}
	if (oc<0 || oc>1)  {PROGERROR("Invalid option/command flag %d!",oc);}
	
	op = MMALLOC(sizeof(op[0]));		/* creates the option with its name */
	OPTNAME(op) = (nm?MSTRDUP(nm):NULL);
	OPTCOM(op) = oc;
	OPTCONTINUED(op) = (cont>0? cont:0);	/* a "continued" option flag - for continuing a long list of option values */
	OPTINDEX(op) = -1;
	OPTPARSG(op) = 0;  pfield_setzeroexp(&OPTPAREX(op));
	
	for (i=0; i<MAXOPTPARAMS; i++) {	/* nulls the option parameters (values) */
		OPTPARAM(op,i) = NULL;
		OPTNPARAM(op,i) = OPTNOPARAM;
	}
	OPTPARNUM(op) = 0;
	fr_aloptions++;
	return op;
}

void	dispose_comoption(optionstruc *op) {
	int	i;
	
	if (OPTCONTINUED(op)<0)  {PROGERROREXIT("Disposing an option %p second time ?!",op);}
	OPTCONTINUED(op) = -1;
	if (OPTNAME(op))  FREE(OPTNAME(op));
	for (i=0; i<MAXOPTPARAMS; i++)
		if (OPTPARAM(op,i))  FREE(OPTPARAM(op,i));
	FREE(op);
	fr_aloptions--;
}

optionstruc*	option_copy(optionstruc *op) {
	optionstruc	*opn;
	int		i;
	
	opn = new_comoption_ext(OPTNAME(op),OPTCOM(op),OPTCONTINUED(op));
	OPTINDEX(opn) = OPTINDEX(op);
	OPTPARSG(opn) = OPTPARSG(op);
	OPTPAREX(opn) = OPTPAREX(op);
	for (i=0; i<MAXOPTPARAMS; i++)
		if (OPTPARAM(op,i))  option_setparam(opn,i,OPTPARAM(op,i),OPTNPARAM(op,i));
	return opn;
}


/**
 * The function option_setparam() sets the ix-th value/parameter of op to the string
 * (a copy of) val, and to the numeric value dv (OPTNOPARAM when number undefined).
 * The values are indexed from 0.
 * The number of values in op is updated as well.
 * There is op as the return value, the given option is modified instead.
 * 
**/

optionstruc*	option_setparam(optionstruc *op, int ix, char *val, long dv) {
	char	buf[20];
	
	if (OPTCONTINUED(op)<0)  {PROGERROR("Accession invalid option %p !",op);}
	if (OPTPARAM(op,ix))  FREE(OPTPARAM(op,ix));
	OPTNPARAM(op,ix) = dv;
	if (!val && dv!=OPTNOPARAM) {
		val = buf;  snprintf(val,15,"%ld",dv);  val[15] = 0;
	}
	OPTPARAM(op,ix) = MSTRDUP(val?val:"");
	if (OPTPARNUM(op)<=ix)  OPTPARNUM(op) = ix+1;
	return op;
}


/**
 * This function prints the current option from the list ol to the buffer buf (of length max).
 * A continued options are concatenated to one string.
 * It returns the next position in the option list.
**/

optionstruc**	option_print_to(optionstruc **ol, char *buf, int max) {
	optionstruc	**x=ol;
	int		i,l;
	char		*sv, *s,*t;
	extern char	*STROPTIONSEP;	/* (from emflex.l - what characters need to be quoted) */
	
	buf[0] = 0;
	if (ol)  for (x=ol; ; x++) {		/* prints the requested option values */
		if (*x? (x!=ol && !OPTCONTINUED(*x)): 1)  break;
		if (!buf[0])  sprintf(buf,"%s%s ",FR_OPTION,OPTNAME(*x));
		for (i=0; i<OPTPARNUM(*x); i++) {
			l = strlen(buf);
			s = OPTPARAM(*x,i);  sv = MMALLOC(2*strlen(s)*sizeof(s[0])+8);
			for (t=sv; *s; *(t++) = *(s++))
				if (*s=='"')  *(t++) = '\\';
			t[0] = 0;	/* escape quotas as \" */
			if (l<max-1) {
				if (strpbrk(sv,STROPTIONSEP)!=NULL)
					snprintf(buf+l,max-l-1,"\"%s\" ",sv);
				else
					snprintf(buf+l,max-l-1,"%s ",sv);
			}		/* (use double-quotes in the print only when really neccessary) */
			FREE(sv);
		}
	}
	if ((int)strlen(buf)>max-5)  {PROGERROR("Possibly truncated option print-out \'%s\'.\n",buf);}
	
	return x;	/* returns the next position in the option list */
}











/******************	Basic frame handling	*******************/
/******************************************************************/


#undef CURDLEV
#define CURDLEV         6


/**
 * This function creates a new empty frame that is a son of the given parent frame.
 * The new frame has no sons, no matrix, no name, no nothing...,
 * but it is arranged as one of the sons of the given parent frame.
 * 
 * The next function disposes the whole frame, including names and all options/commands.
 * If rec==1 is given, then also all sons are recursively disposed, and the frame is
 * taken from the son-list of its parent.
 * (Be careful not to refer the name or options/commands anywhere else!)
 * 
 * At any time the whole frame structure must form a simple tree!!!
**/

framestruc*	new_frame(framestruc *parent) {
	framestruc	*fr;
	
	fr = MMALLOC(sizeof(fr[0]));
	FRMAGICSET(fr);		/* (a magic number is used to test correct access to frames) */
	FRNAME(fr) = FREXTENS(fr) = FRCOMMENT(fr) = NULL;
	FRSONS(fr) = NULL;
	FRNUMSONS(fr) = 0;
	FRPARENT(fr) = NULL;
	frame_setson(fr,parent);	/* adds to the parent's list of sons */
	
	FRMATRIX(fr) = NULL;	/* clears all important values in this frame */
	FRPFINDEX(fr) = FRPFINDEX_SAVE(fr) = -1;
	FRMATRIXROWS(fr) = FRMATRIXCOLS(fr) = 0;
	FRVALUESIG(fr) = 0;  pfield_setzeroexp(&FRVALUEEXP(fr));
	FRNUMBER(fr) = 0;
	FROPTIONS(fr) = NULL;
	FRCOMMANDS(fr) = NULL;
	FRLASTOPTI(fr) = -1;	/* the "index" of the last added option */
	FRDELMARK(fr) = FRTREEMARK(fr) = 0;
	fr_alframes++;
	return fr;
}

void	dispose_frame_ext(framestruc *fr, int rec, int par) {
	framestruc	**rr;
	int		nums;
	
	DEBUG(CURDLEV+1,"Call to dispose frame %p, rec=%d, par=%d\n",fr,rec,par);
	FRMAGIC(fr);	/* (a magic number test is used to check correct access to frames) */
	nums = FRNUMSONS(fr);
	if (nums==-12345)  {PROGERROR("Probably cycled recursion when disposing frames!"); return;}
	FRNUMSONS(fr) = -12345;
				/* recursive disposal of all sons of this frame */
	if (rec && FRSONS(fr))
		for (rr=FRSONS(fr); *rr; rr++) {
			dispose_frame_ext(*rr,1,0);	/* (must not call with par==1 !!!) */
			*rr = (void*)1;
			nums--;
		}
	if (nums!=0)  {PROGERROR("Wrong number of sons %d when disposing.",nums);}
	
				/* removing from the list of sons of the parent frame */
	if (par)  frame_setson(fr,NULL);
				/* disposing the matrix, option/command lists, the the son list */
	if (FRMATRIX(fr)!=NULL)
		dispose_ematrix(FRMATRIX(fr));
	if (FROPTIONS(fr))
		dispose_alist_copts(FROPTIONS(fr));
	if (FRCOMMANDS(fr))
		dispose_alist_copts(FRCOMMANDS(fr));
	if (FRSONS(fr))  alist_free(FRSONS(fr));
	
	if (FRNAME(fr))  FREE(FRNAME(fr));
	if (FREXTENS(fr))  FREE(FREXTENS(fr));
	if (FRCOMMENT(fr))  FREE(FRCOMMENT(fr));
	FRMAGICCLR(fr);		/* finally, trashing this frame structure */
	FREE(fr);
	fr_alframes--;
}

	/* this is for use in disposing whole lists, be careful... */
void	dispose_frame_ls(framestruc *fr) {
	if (FRPARENT(fr) || FRSONS(fr))  {PROGERROR("does not work with parent or sons!");}
	dispose_frame_ext(fr,0,0);
}

	/* this is for final reporting of lost frames and options/commands: */
void	frame_finishall(void) {
	if (fr_aloptions>33)  DEBUG(CURDLEV-2,"Frame_ finished with %d options/commands not freed!\n",fr_aloptions);
	if (fr_alframes>10)  DEBUG(CURDLEV-2,"Frame_ finished with %d frames not freed!\n",fr_alframes);
}


/**
 * This function adds a new son fr to the parent frame parent.
 * If where is given, then the son fr is inserted into the exact position of a son where
 * among the other sons.
 * If there is a current parent, then the frame fr is deleted from the list of sons
 * of the current parent.
**/

void	frame_setson_ext(framestruc *fr, framestruc *parent, framestruc *where) {
	framestruc	*pold, **rr;
	
	if (!fr)  return;
	pold = FRPARENT(fr);
	if (pold) {		/* removing from the son-list of its old parent */
		if (!FRSONS(pold))  rr = NULL;
		else  for (rr=FRSONS(pold); *rr && *rr!=fr; rr++) ;
		if (rr?(*rr==fr):0) {
			alist_delete(rr);
			FRNUMSONS(pold)--;
		} else  {PROGERROR("Probably broken list of sons (%p) in the parent %p.",fr,pold);}
	}
	FRPARENT(fr) = parent;
	if (parent) {		/* the new parent */
		if (!where) {
			FRSONS_XD(parent) = alist_append(FRSONS(parent),fr);
		} else {
			if (!FRSONS(parent))  rr = NULL;
			else  for (rr=FRSONS(parent); *rr && *rr!=where; rr++) ;
			if (rr?(*rr!=where):1)  {PROGERROR("Nonexistent \"where\" position given.");}
			FRSONS_XD(parent) = alist_insert(FRSONS(parent),rr,fr);
		}
		FRNUMSONS(parent)++;
	}
}


/**
 * This function sets a new name nm1nm2 and a new extension ext for the given frame.
 * The given strings are internally copied before setting.
**/

void	frame_setname_ext(framestruc *fr, char *nm1, char *nm2, char *ext) {
	char	*buf;
	
	if (nm1 || nm2) {
		buf = MMALLOC(((nm1?strlen(nm1):0)+(nm2?strlen(nm2):0)+3)*sizeof(buf[0]));
		if (nm1)  strcpy(buf,nm1);  else  buf[0] = 0;
		if (nm2)  strcat(buf,nm2);
		if (FRNAME(fr))  FREE(FRNAME(fr));
		FRNAME(fr) = buf;
	}
	if (ext) {
		if (FREXTENS(fr))  FREE(FREXTENS(fr));
		FREXTENS(fr) = MSTRDUP(ext);
	}
}


/**
 * This function sets a new matrix for the frame.
 * The given new matrix en must be a private copy for this frame, and the old matrix is disposed!
**/

void	frame_setmatrix(framestruc *fr, ematrix *en) {
	
	if (FRMATRIX(fr))  dispose_ematrix(FRMATRIX(fr));
	FRMATRIX(fr) = en;
	FRMATRIXROWS(fr) = en?ROWSM(en):0;  FRMATRIXCOLS(fr) = en?COLSM(en):0;
	FRPFINDEX(fr) = pfield_curindex();
}


/**
 * This function extracts a copy of the matrix of the given frame.
 * (With a name inherited from the frame...)
**/

ematrix*	frame_extractmatrix(framestruc *fr) {
	ematrix		*em;
	
	if (fr?!FRMATRIX(fr):1)  return NULL;
	em = ematrix_copy(FRMATRIX(fr));
	if (!EMNAME(em))  EMSETNAME(em,FRNAME(fr));
	return em;  
}


/**
 * This function creates a copy of the given frame fr, with a copy of a matrix.
 * If nnm is given, then the new frame has this name.
 * If copt==1, then the options of fr are copied as well.
 * If ccom==1, then the commands of fr are copied as well.
 * The return value is the new frame, which has no parent and no sons.
**/

framestruc*	frame_copy_ext(framestruc *fr, char *nnm, int copt, int ccom) {
	framestruc	*frn;
	optionstruc	**ol,*op;
	
	frn = new_frame(NULL);
	if (FRMATRIX(fr)) {
		frame_setmatrix(frn,ematrix_copy(FRMATRIX(fr)));
		FRPFINDEX(frn) = FRPFINDEX(fr);
	}
	if (nnm)  frame_setname(frn,nnm);
	else  frame_setname_ext(frn,FRNAME(fr),NULL,FREXTENS(fr));
	if (FRCOMMENT(fr) && !nnm)  FRCOMMENT(frn) = MSTRDUP(FRCOMMENT(fr));
	FRNUMBER(frn) = FRNUMBER(fr);
	
	if (copt && FROPTIONS(fr))
		for (ol=FROPTIONS(fr); *ol; ol++) {
			op = option_copy(*ol);
			frame_addoption(frn,op,OPTINDEX(op));
		}
	if (ccom && FRCOMMANDS(fr))
		for (ol=FRCOMMANDS(fr); *ol; ol++) {
			op = option_copy(*ol);
			frame_addcommand(frn,op,OPTINDEX(op));
		}
	FRLASTOPTI(frn) = FRLASTOPTI(fr);
	return frn;
}


/**
 * This function is called to collect all frames in the tree of the given frame fr
 * to one output list (including fr itself).
 * The returned list points directly to the frames in the tree, not to copies.
 * The parameter ord determines the order of collection:
 *  ord==0 means breadth-first, ord==1 depth-first, and ord==2 reversed depth-first.
 * 
 * Do not forget to free the list later.
**/

framestruc**	frame_gettree_ext(framestruc *fr, int ord) {
	framestruc	**fin, **fout;
	
	if (!fr)  return NULL;
	fin = alist_append(NULL,fr);
	fout = NULL;
	frame_gettree_recur(fin,ord,&fout);
	alist_free(fin);
	return fout;
}

void	frame_gettree_recur(framestruc **fin, int ord, framestruc ***fout) {
	framestruc	**x;
	
	if (!fin)  return;
	if (ord==0)  for (x=fin; *x; x++)
		*fout = alist_append(*fout,*x);
	for (x=fin; *x; x++) {
		if (ord==1)  *fout = alist_append(*fout,*x);
		if (FRTREEMARK(*x)==11)  {PROGERROR("Probably cycled recursion in frames!"); return;}
		FRTREEMARK(*x) = 11;
		frame_gettree_recur(FRSONS(*x),ord,fout);
		FRTREEMARK(*x) = 0;
		if (ord==2)  *fout = alist_append(*fout,*x);
	}
}




/**
 * General printing function for the frame - use for debug or informal output.
 * (Not for formal writing to files!)
 * The function prints the frame fr to the output fout, using prefix pref.
 * If verb>0+, then more is printed, and the frame matrix is printed if mat>0 (>1).
 * 
 * The next function prints the whole frame tree of fr (possibly to a string buffer).
 * It has few more "limiting" parameters:
 * - maxbuf limits the number of characters printed to a string buffer buf[] (nothing for file).
 * - maxson limits the number of sons printed for a frame (first and last few are printed out).
 * - maxdepth limits the depth of the frame tree printed out.
 * Values of <0 are substituted by the default limits.
**/

void	frame_print_ext(FILE *fout, framestruc *fr, int verb, int mat, char *pref) {
	int		i;
	optionstruc	**op;
	char		pref2[50];
	
	if (!fr) {
		PROGERROR("Not printing NULL frame.");  return;
	}
	fprintf(fout,"%sFrame %p [%.18s%s%.4s], parent [%.9s], %s%d\n",pref,fr,
			FRNAME(fr)?FRNAME(fr):"",FREXTENS(fr)?" ":"",FREXTENS(fr)?FREXTENS(fr):"",
			(FRPARENT(fr)?(FRNAME(FRPARENT_XD(fr))?FRNAME(FRPARENT_XD(fr)):""):""),
			FRSONS_XD(fr)?"sons ":"",FRNUMSONS(fr));
	if (verb>0 && FRCOMMENT(fr))  fprintf(fout,"%s  - comment \"%s\"\n",pref,FRCOMMENT(fr));
	if (verb>0)  fprintf(fout,"%s  - matrix %d x %d,  number %d,  pf value %s,\n",pref,
			FRMATRIXROWS(fr),FRMATRIXCOLS(fr),FRNUMBER(fr),pfield_printvalue(25,FRVALUEEXP(fr),FRVALUESIG(fr)));
	if (FRMATRIX(fr)) if (ROWSM(FRMATRIX(fr))+COLSM(FRMATRIX(fr))>2) {
		strncpy(pref2,pref,40);  pref2[39] = 0;  strcat(pref2,"  :\t");
		if (mat>1)  ematrix_fprint_pref(fout,FRMATRIX(fr),pref2);
		else if (mat>0)  ematrix_fprint_nofr(fout,FRMATRIX(fr),pref2);
	}
	if (verb>1 && FROPTIONS(fr)) {
		fprintf(fout,"%s  - options:  ",pref);
		for (op=FROPTIONS(fr),i=0; *op; op++,i++)
			if (i<4)  fprintf(fout,"%s %.8s %.8s ; ",OPTNAME(*op),OPTPARAM(*op,0)?OPTPARAM(*op,0):"",OPTPARAM(*op,1)?OPTPARAM(*op,1):"");
		fprintf(fout,"..(%d)\n",i);
	}
	if (verb>1 && FRCOMMANDS(fr)) {
		fprintf(fout,"%s  - commands:  ",pref);
		for (op=FRCOMMANDS(fr),i=0; *op; op++,i++)
			if (i<4)  fprintf(fout,"%s %.8s %.8s ; ",OPTNAME(*op),OPTPARAM(*op,0)?OPTPARAM(*op,0):"",OPTPARAM(*op,1)?OPTPARAM(*op,1):"");
		fprintf(fout,"..(%d)\n",i);
	}
}


#define	PTREE_MAXDEPTH	20
#define	PTREE_MAXSONS	99999

void	frame_printtree_ext(FILE *fout, char *buf, framestruc *fr, int verb, char *pref,
						 int maxbuf, int maxson, int maxdepth) {
	framestruc	*ff, ***yr;
	ematrix		*ffem;
	int		*ir,*lr, ii,k,l, maxsl,maxsl2;
	char		bpr[100];
	
	if (!fout && !buf)  {PROGERROR("a file or a buffer must be given!"); return;}
	if (maxbuf<0 && buf)  {PROGERROR("maxbuf limit must be given when printing to a buffer!"); maxbuf=10;}
	if (maxdepth<0)  maxdepth = PTREE_MAXDEPTH;
	if (maxson<0)  maxson = PTREE_MAXSONS;
	maxsl = maxson/2;  maxsl2 = maxson/3;
	if (maxsl<2)  maxsl = 2;
	if (buf)  buf[0] = 0;
	yr = MMALLOC((maxdepth+2)*sizeof(yr[0]));
	ir = MMALLOC(2*(maxdepth+2)*sizeof(ir[0]));  lr = ir+maxdepth+2;
	
			/**
			 * Here we print the subtree of given fr:
			 * We start a depth-first search in the tree of fr,
			 * and we print each frame when it is reached the first time.
			 * The array yr[] is used to keep the current list of sons
			 * at each level of the search - incremented along the way.
			 * The array ir[] counts the printed sons, and lr[] keeps
			 * the lengths of the whole son lists.
			**/
	yr[0] = NULL;
	for (k=0; k>=0; ) {
		if (yr[k]==NULL) {
			yr[k] = (k<=0?&fr: FRSONS(*yr[k-1]));
			ir[k] = 0;
			lr[k] = (k<=0?1: alist_getlength(yr[k]));
			if (lr[k]<=0) { k--; continue; }
		} else {
			if (k>0)  yr[k]++;
			if (k<=0 || *yr[k]==NULL) { k--; continue; }
		}
			/**
			 * We print the actual frame *yr[k]=ff, according to buf and verb.
			 * However, if the list of sons is longer than the given limit maxson,
			 * then only first half and last third of them is printed which
			 * is indicated by a printed message.
			 * Printed format depends on buf/fout and on verb.
			**/
		ff = *yr[k];  ii = ++ir[k];
		if (lr[k]>maxson && maxsl<lr[k]-maxsl2 && verb<=1) {
			if (ii==maxsl) {
			  if (!buf)
				fprintf(fout,"%s%s%*s(%d.%d-%d) \t...  skipping sons in a long list  ...\n%s",
					(verb>0?"\n":""),pref,2*k,"",k+1,maxsl,lr[k]-maxsl2,(verb>0?"\n":""));
			  else if ((int)strlen(buf)<maxbuf-99-2*k)
			  	snprintf(buf+strlen(buf),90+2*k,"%s%s%*s(%d.%d-%d) ...  skipping sons in a long list  ...\n%s",
					(verb>0?"\n":""),pref,2*k,"",k+1,maxsl,lr[k]-maxsl2,(verb>0?"\n":""));
			}
			if (ii>=maxsl && ii<=lr[k]-maxsl2)  continue;
		}
		ffem = FRMATRIX(ff);
		if (!ffem)  bpr[0] = 0;
		else  snprintf(bpr,80,"m%dx%d",ROWSM(ffem),COLSM(ffem));
		if (buf) {
			l = strlen(buf);
			if (l<maxbuf-5)  snprintf(buf+l,maxbuf-l-4,"%*s(%d.%d)fr [%.15s] %s  \"%.*s\"\n",
				2*k,"",k+1,ii,FRNAME(ff),bpr,44-3*k,FRCOMMENT(ff)?FRCOMMENT(ff):"");
			buf[maxbuf-1] = 0;
		} else if (verb<=1) {
			fprintf(fout,"%s%*s(%d.%d)fr [%.15s] %s  \"%.*s\"\n",pref,2*k,"",
				k+1,ii,FRNAME(ff),bpr,55-3*k,FRCOMMENT(ff)?FRCOMMENT(ff):"");
		} else {
			if (verb>2)  fprintf(fout,"\n");
			snprintf(bpr,80,"%s%*s(%d.%d)  ",pref,3*k,"",k+1,ii);
			frame_print_ext(fout,ff,verb-1,(verb>3),bpr);
		}
			/**
			 * We dive deeper in the tree here (after printing the frame),
			 * unless the maximal given bound maxdepth is reached.
			 * Deeper frames in the tree are simply skipped.
			**/
		if (k<maxdepth)  yr[++k] = NULL;
		else if (FRSONS(*yr[k])!=NULL) {
			if (!buf)
				fprintf(fout,"%s%s%*s(%d.x) \t...  skipping too deeply nested frames  ...\n%s",
					(verb>0?"\n":""),pref,2*k+2,"",k+2,(verb>0?"\n":""));
			else if ((int)strlen(buf)<maxbuf-99-2*k)
			  	snprintf(buf+strlen(buf),90+2*k,"%s%s%*s(%d.x) \t...  skipping too deeply nested frames  ...\n%s",
					(verb>0?"\n":""),pref,2*k+2,"",k+2,(verb>0?"\n":""));
		}
	}
	if (fout && verb>2)  fprintf(fout,"\n");
	FREE(yr);  FREE(ir);
}













/******************	Frame options (commands)	*******************/
/**************************************************************************/


#undef CURDLEV
#define CURDLEV         6



/**
 * These functions add the given option or command to the frame fr.
 * The option/command must already be properly created using the "option" functions below.
 * (It may not be used in any other frame!)
 * Optionally, the given index li indicates that a new option was added - so that
 * functions depending on the current options need not  look through all options
 * each time.
 * 
 * If an option is added, it is first checked against the list of all known options
 * optdescdefs[], including the number of parameters (values).
 * (The option name may be also matched by an exact prefix, without respect to the values.)
**/

void	frame_addoption(framestruc *fr, optionstruc *op, int li) {
	int	i,k;
	char	*s;
	
	FROPTIONS(fr) = alist_append(FROPTIONS(fr),op);
	FRLASTOPTI(fr) = li;
	for (i=k=0; k<4 && optdescdefs[i].name; i++) {
		s = optdescdefs[i].name;
		if (k>=0 && strncmp(OPTNAME(op),s,strlen(s))==0)  k = 1;
		if (OPTISNAME(op,s))  k = 2;
		if (k>=2)  k = (optdescdefs[i].numpar<0 || OPTPARNUM(op)==optdescdefs[i].numpar)?4:-1;
	}
	if (k<=0)  {USERERROR("Unknown option \'%s\' or wrong number of option values (%d).",OPTNAME(op),OPTPARNUM(op));}
}

void	frame_addcommand(framestruc *fr, optionstruc *op, int li) {
	
	FRCOMMANDS(fr) = alist_append(FRCOMMANDS(fr),op);
	junk = li;	/* (should we use this as above???) */
}


/**
 * Here we get the options of given name onm from the given frame fr (plus ancestors).
 * If all==1, then all options of the given name are extracted from the frame and
 * its ancestors.
 * If all==0, only the last instance of the option is returned (possibly as a "continued"
 * chain of options).
 * If anc==1, then the options are taken also from the ancestors of fr.
 * The return value of a list of the extracted options, or NULL.
 * 
 * Additionally, special options "erase" and "eraseall" can erase the last one
 * or all previous option occurences from the frame (and ancestors).
 * 
 * Do not alter the returned options!!!
 * It is suggested not to use this function directly, but to extract required information
 * in a separate function, and then use only these (copied!) information.
 * 
**/

optionstruc**	frame_getoption_ext(framestruc *fr, char *onm, int anc, int all) {
	optionstruc	**ol, **x,**y;
	int		i;
	
	if (!onm || !fr)  {PROGERROREXIT("Must give fr and onm !");}
	ol = NULL;
	if (anc && all && FRPARENT(fr)) {	/* gets options from the parent if anc requested */
		ol = frame_getoption_ext(FRPARENT(fr),onm,1,all);
	}
					/* gets options named onm from this frame */
	for (x=FROPTIONS(fr); x?*x:0; x++) {
		if (OPTISNAME(*x,onm)) {
			ol = alist_append(ol,*x);
					/* or, erases the previously got options */
		} else if (OPTISNAME(*x,FR_ERASEOPT) || OPTISNAME(*x,FR_ERASEALLOPT)) {
			for (i=0; ol && i<MAXOPTPARAMS; i++)
			  if (OPTPARAM(*x,i)? strcasecmp(OPTPARAM(*x,i),onm)==0:0) {
				if (OPTISNAME(*x,FR_ERASEALLOPT) || alist_getlength(ol)<=1) {
					alist_free(ol);  ol = NULL;
				} else  alist_delete_last(ol);
			}
		}
	}				/* get the list from ancestors if nothing found in fr */
	if (!all && anc && !ol && FRPARENT(fr)) {
		ol = frame_getoption_ext(fr,onm,anc,1);
	}
	if (ol && !all) {		/* finally, extracting only the last option if all=0 */
		for (i=alist_getlength(ol)-1; i>=0; i--)
			if (!OPTCONTINUED(ol[i]))  break;
		y = ol;  ol = NULL;
		for (x=y+i; *x; x++)  ol = alist_append(ol,*x);
		alist_free(y);
	}
	return ol;
}


/**
 * This function extracts the option parameters from (all or the last) options of name onm
 * of the frame fr and its ancestors.
 * If nv>0, then exactly nv parameters are taken from each option instance (substituting
 * "" for nonexistent parameters), otherwise, all parameters are taken from each instance.
 * The parameters are returned as a list of string copies.
 * (The list and the strings should be freed after use.)
**/

char**	frame_getoptionval_ext(framestruc *fr, char *onm, int all, int nv) {
	optionstruc	**ol, **x;
	char		*s, **cl;
	int		i;
	
	if (nv<0 || nv>MAXOPTPARAMS)  {PROGERROR("Wrong value of nv = %d given.\n",nv);}
	cl = NULL;
	ol = frame_getoption_ext(fr,onm,1,all);	/* gets all options of the given name */
	
	for (x=ol; x?*x:0; x++) {		/* extracts the requested option values */
		if (nv>0 && OPTCONTINUED(*x))  continue;
		for (i=0; i<MAXOPTPARAMS; i++) {
			s = OPTPARAM(*x,i);
			if (nv>0 && s==NULL)  s = "";
			if (nv>0 && i>=nv)  s = NULL;	/* (the strings are copied for the list...) */
			if (s)  cl = alist_append(cl,MSTRDUP(s));
			else if (!cl)  cl = new_alist(4);
		}			/* (at least an empty list is returned when no values) */
	}
	if (ol)  alist_free(ol);
	return cl;	/* returns the list of values */
}


/**
 * This function copies the options of names given in onml from fr (and ancestors if anc==1).
 * All the option copies are returned in an option list (ready to be given to another frame).
 * 
**/

optionstruc**	frame_getoptionlist_ext(framestruc *fr, char **onml, int anc, int all) {
	optionstruc	**ol,**x, **out;
	char		**ss;
	
	out = NULL;
	for (ss=onml; ss?*ss:0; ss++) {
		ol = frame_getoption_ext(fr,*ss,anc,all);	/* gets the options of the given name */
		for (x=ol; x?*x:0; x++)
			out = alist_append(out,option_copy(*x));
		if (ol)  alist_free(ol);
	}
	return out;
}


/**
 * This function extracts the option numeric parameters from the last option of name onm
 * of the frame fr (and its ancestors if anc==1).
 * The parameters are returned in the given array nout[] (of size nsz), by default value def.
**/

int	frame_getoptionnum_ext(framestruc *fr, char *onm, int anc, long nout[], int nsz, long def) {
	optionstruc	**ol;
	int		i,r;
	
	ol = frame_getoption_ext(fr,onm,anc,0);	/* gets the last option of the given name */
	r = (ol?(ol[0]!=NULL):0);
	if (nout && r) {
	        for (i=0; i<nsz; i++)  nout[i] = def;	
	  					/* extracts the option numeric parameters */
                for (i=0; i<nsz && i<MAXOPTPARAMS; i++)
			if (OPTNPARAM(ol[0],i)!=OPTNOPARAM)  nout[i] = OPTNPARAM(ol[0],i);
	}
	if (ol)  alist_free(ol);
	return r;
}


/**
 * This function extracts the option parameters from (all or the last) options of name onm
 * of the frame fr and its ancestors, and prints them into strings as on the input.
 * The strings are returned as a list of string copies.
 * (The list and the strings should be freed after use.)
**/

char**	frame_getoptionprint_ext(framestruc *fr, char *onm, int anc, int all) {
	optionstruc	**ol, **x;
	char		*buf, **cl;
	
	buf = MMALLOC(5100*sizeof(buf[0]));
	cl = NULL;				/* gets all options of the given name */
	ol = frame_getoption_ext(fr,onm,anc,all);
	
	if (ol)  for (x=ol; *x; ) {		/* prints the requested option values */
		x = option_print_to(x,buf,5000);
		cl = alist_append(cl,MSTRDUP(buf));
	}
	if (ol)  alist_free(ol);
	FREE(buf);
	return cl;	/* returns the list of strings */
}




































