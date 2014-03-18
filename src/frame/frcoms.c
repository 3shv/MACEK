
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
 * This source file contains declarations of command-handling functions.
 * The first part consists of definitions of the handling functions
 * (they have unified parameter-format).
 * Then there is a definition of the array comhdefs[] that holds pointers
 * to these handling functions, in association with names, comments,
 * and parameter descriptions of the commands.
 * 
 * 
 * When writing new command handles, follow these RULES:
 * 
 * make similar to existing handles.....
 * 
 * getting parameters.....
 * 	- input matrices are already copied, frames are not(!!!),
 * 	- output should be copied here (ready for putting into the frame tree),
 * 	  output frames must not have parents, but they may have descendants,
 * 	- even the output matrices must be copied again, since the input ones are disposed later!
 * 	- set the names and options here for the output frames if you want...
 * 
 * 
 * It would be great to make a common printing interface for command handles...
 * 
**/





#include "macek.h"  
#include "fr.h"

#include <sys/wait.h>









/******************	Frame i/o and handling	*******************/
/******************************************************************/


#undef CURDLEV
#define CURDLEV         5


int	ret_yes=1, ret_no=0, ret_def=2;	/* used for yes/no lists - pointers point to these variables */


/**
 * One may control output printing verbosity within commands by a command and by an option...
 * The value of frcom_printbrief is set in frame_applycommand() based on the option @prbrief,
 * and the value of frcom_verbose is adjusted here by special commands.
**/

int	frcom_verbose = 1;
int	frcom_printbrief = 0;

int	frcom_verbose_out(void) {
	return frcom_verbose+frcom_printbrief;
}

int	COMFUNCTION(frcom_setverbose) {
	int	adj = inum[0];
	frcom_verbose += adj;
	return 0;	/* no return values */
}



/**
 * Here we defined frame-, matrix- and tree-printing functions, simple text printing, and setting a new name.
 * They do not change the frame tree.
**/

int	COMFUNCTION(frcom_prframe) {
	framestruc	**flin = (framestruc**)inls[0];
	framestruc	**x;
	
	if (flin) for (x=flin; *x; x++) {
		if (frcom_verbose_out()>0)  FROUTPUT(*x,"    ~\t");
		else  FROUTPUTS(*x,"    ~\t");
	}
	return 0;	/* no return values */
}

int	COMFUNCTION(frcom_print) {
	framestruc	**flin = (framestruc**)inls[0];
	framestruc	**x;
	ematrix		*em;
	int		k,xpf;
	
	for (k=0,x=flin; x?*x:0; x++,k++) {
		em = frame_extractmatrix(*x);
		if (em) {
			xpf = pfield_curindex();
			if (FRPFINDEX(*x)!=xpf)  pfield_switchto_fast(FRPFINDEX(*x));
			else  xpf = -1;
			if (frcom_verbose_out()>0)  OUTPUT("Matrix of the frame %p [%s] in %s: \"%.35s\"\n",
					*x,FRNAME(*x),pfield_curname(),FRCOMMENT(*x)?FRCOMMENT(*x):"");
			if (frcom_verbose_out()>0)  EMATOUTPUT(em,"    ~\t");
			else  EMATOUTPUTS(em,"    ~\t");
			dispose_ematrix(em);
			if (xpf>=0)  pfield_switchto_fast(xpf);
		} else  SOUTPUT("\tThere is no matrix in the frame [%s].\n",FRNAME(*x));
	}
	if (frcom_verbose_out()>0 && k>10)  OUTPUT("Printed out total %d matrices.\n",k);
	return 0;	/* no return values */
}

int	COMFUNCTION(frcom_printmore) {
	ematrix		**elin = (ematrix**)inls[0];
	ematrix		**el;
	long		i,j, *rem;
	
	rem = MMALLOC(sizeof(rem[0])*(elin?alist_getlength(elin)+2:2));
	for (el=elin,i=0; el?*el:0; el++,i++) {
		if (frcom_verbose_out()>0) {
			if (el!=elin) SOUTPUT("\n\t\t========  [%.25s]  ========\n\n",EMNAME(*el)?EMNAME(*el):"");
			EMATOUTPUT(*el,"    ~\t");
		}
		rem[i] = ematrix_printmore(*el,frcom_verbose_out());
		if (frcom_verbose_out()>=0)
		  for (j=i-1; j>=0; j--)  if (rem[j]==rem[i]) {
			if (frcom_verbose_out()>2)  EMATOUTPUT(elin[j],"    ~\t\t");
			OUTPUT(" ## Repeated hash value %ld - same as for matrix #%ld [%s].\n",
				rem[j],j+1,EMNAME(elin[j]));
			break;
		}
	}
	FREE(rem);
	return 0;	/* no return values */
}

int	COMFUNCTION(frcom_printbasecirc) {
	ematrix		**elin = (ematrix**)inls[0];
	int		ix1 = (inum?inum[0]:0), ix2 = (inum?inum[1]:0), ix3 = (inum?inum[2]:0);
	int		whp = (strcmp("circ",instr?instr[0]:"")? 1:11);
	ematrix		**el;
	int		i, cix[10];
	
	cix[0] = ix1;  cix[1] = ix2;
	cix[2] = ix3;  cix[3] = 0;
	for (el=elin,i=0; el?*el:0; el++,i++) {
		if (frcom_verbose_out()>0) {
			if (el!=elin) SOUTPUT("\n\t\t========  [%.25s]  ========\n\n",EMNAME(*el)?EMNAME(*el):"");
			EMATOUTPUT(*el,"    ~\t");
		}
		ematrix_printbasecirc(*el,whp,frcom_verbose_out(),cix);
	}
	return 0;	/* no return values */
}


int	COMFUNCTION(frcom_printtree) {
	framestruc	**flin = (framestruc**)inls[0];
	framestruc	**x;
	
	if (flin) for (x=flin; *x; x++) {
		if (frcom_verbose_out()>0)  OUTPUT("Printing the subtree of the frame %p [%s]:\n",*x,FRNAME(*x));
		if (frcom_verbose_out()<3)
			frame_printtree_norep(printout,*x,0,frcom_verbose_out()>1?40:16,"    ~\t");
		else	frame_printtree(printout,*x,frcom_verbose_out()-2,"    ~\t");
		if (frcom_verbose_out()>0)  OUTPUT("------------------------------------\n");
	}
	return 0;	/* no return values */
}

int	COMFUNCTION(frcom_printtext) {
	char	*text = instr?instr[0]:"";
	OUTPUT("%s\n\n",text);
	return 0;	/* no return values */
}


int	COMFUNCTION(frcom_setname) {
	framestruc	**flin = (framestruc**)inls[0];
	framestruc	**x;
	char		buf[100];
	
	for (x=flin; (x?*x:0) && (instr?instr[0]:0); x++) {
		snprintf(buf,90,instr[0],(int)(x-flin)+1);
		buf[90] = 0;  frame_setname(*x,buf);
	}
	return 0;	/* no return values */
}


/**
 * Here we define frame and frame-tree file writing functions.
 * The filename is taken from the name of the current frame, or "noname", or
 * can be set by the instr[0] string parameter (only for the first frame!).
 * Moreover, we may request writing additional matroid comments to the file
 * by given value of morecom>0.
**/

int	COMFUNCTION(frcom_write) {
	framestruc	**flin = (framestruc**)inls[0];
	int		morecom = (inum?inum[0]:-1);
	char		*dir, *fn = instr?instr[0]:NULL;
	extern char	*fopen_fname_used;	/* (records the complete name found in the path during open) */
	framestruc	**x;
	int		k;
	
	dir = (fn? (fn[0]&&fn[strlen(fn)-1]=='/'? fn:NULL):NULL);
	for (x=flin; x?*x:0; x++) {
		if (x==flin && fn && !dir)  frame_setname(*x,fn);
		if (frcom_verbose_out()>0)  OUTPUT("Writing the frame [%s]%s  \"%.35s\".\n",
				FRNAME(*x),morecom>0?" (commented)":"",FRCOMMENT(*x)?FRCOMMENT(*x):"");
		if (frcom_verbose_out()>2)  FROUTPUT(*x,"    ~\t");
		else if (frcom_verbose_out()>1)  FROUTPUTS(*x,"    ~\t");
		
		if (morecom>0)  k = frame_writed_onemc(*x,dir,FRCOMMENT(*x));
		else  k = frame_writed_one(*x,dir,FRCOMMENT(*x));
		
		if (k>=0)  OUTPUT("Written to \"%s\"; OK.\n",fopen_fname_used);
		else  OUTPUT("Write \"%s\" ERROR!\n",FRNAME(*x));
	}
	return 0;	/* no return values */
}

int	COMFUNCTION(frcom_writetree) {
	framestruc	**flin = (framestruc**)inls[0];
	char		*dir, *fn = instr?instr[0]:NULL;
	extern char	*fopen_fname_used;	/* (records the complete name found in the path during open) */
	framestruc	**x;
	int		k;
	char		*buf;
	
	dir = (fn? (fn[0]&&fn[strlen(fn)-1]=='/'? fn:NULL):NULL);
	buf = MMALLOC(50600);
	if (flin)  for (x=flin; *x; x++) {
		if (x==flin && fn && !dir)  frame_setname(*x,fn);
		if (frcom_verbose_out()>0)  OUTPUT("Writing the subtree of the frame %p [%s]  \"%.35s\".\n",
				*x,FRNAME(*x),FRCOMMENT(*x)?FRCOMMENT(*x):"");
		if (frcom_verbose_out()>3)  frame_printtree_norep(printout,*x,frcom_verbose_out()/2-2,16,"    ~\t");
		else if (frcom_verbose_out()>2)  FROUTPUT(*x,"    ~\t");
		else if (frcom_verbose_out()>1)  FROUTPUTS(*x,"    ~\t");
		
		strcpy(buf,"\nFrame tree:\n");
		frame_printtree_to(buf+strlen(buf),50000,*x,12);  buf[50000] = 0;
		if (strlen(buf)>50000-40)  sprintf(buf+strlen(buf),"\n\t.......(possibly unfinished printout).......\n");
		k = frame_writed_tree(*x,dir,buf);
		
		if (k>=0)  OUTPUT("Written to \"%s\"; OK.\n",fopen_fname_used);
		else  OUTPUT("Write \"%s\" ERROR!\n",FRNAME(*x));
	}
	FREE(buf);
	return 0;	/* no return values */
}


/**
 * This reads an additional input frame (with a possible subtree) from the given
 * string parameter, and arranges the new frame within the frame tree.
 * Commands in the new frame are processed before returning.
 * The second function acts similarly, with the diff that only all matrices
 * from the scanned frame are returned (and stored).
**/

int	COMFUNCTION(frcom_read) {
	char		*rdnm = instr?instr[0]:"NO";
	framestruc	*in;
	int		r;
	
	in = frame_doinput(rdnm);
	r = in? frame_processcommands(in):-111;
	if (r<0)  {USERERROR("An error %d occured when reading \"%s\", cancelling.",r,rdnm);}
	else  outls[0] = alist_append(outls[0],in);
	return FRES_STOREFR|FRES_STOREREM;	/* of type FCRET_ONE */
}

int	COMFUNCTION(frcom_mread) {
	char		*rdnm = instr?instr[0]:"NO";
	ematrix		**el;
	
	el = frame_inputmatrices(rdnm);
	if (!el)  {USERERROR("No matrix was read from \"%s\"...",rdnm);}
	else  outls[0] = alist_applist(outls[0],el);
	return FRES_STOREMX|FRES_STOREMXNM|FRES_STOREREM;	/* of type FCRET_LIST */
}


/**
 * Here we define function appending given (additional) input string to the given
 * frame(s), appending commands/options, and openning the matrix again for input.
 * Nothing is returned, the frame itself is modified.
**/

int	COMFUNCTION(frcom_append) {
	framestruc	**flin = (framestruc**)inls[0];
	char		*apstr = instr?instr[0]:NULL;
	framestruc	**x, *f;
	
	for (x=flin; x&&apstr?*x:0; x++) {
		OUTPUT("Appending \'%.20s\' to frame \"%s\"...\n",apstr,FRNAME(*x));
		f = frame_doinput_ap(apstr,*x);
		if (f!=*x)  {PROGERROR("The frame address has changed after an append command %p!=%p!",*x,f);}
	}
	return 0;	/* no return values */
}


/**
 * This is a general manipulating function for the frame tree.
 * It can move, copy or delete frames depending on the input and output parameters.
 * (The functionality is actually defined by the parameter processing in frame_applycommand(),
 * this function only copies the given frames to the output list.)
 * The next one is a similar function just for matrices (copied to new frames).
**/

int	COMFUNCTION(frcom_move) {
	framestruc	**flin = (framestruc**)inls[0];
	framestruc	**x;
	
	outls[0] = new_alist(10);
	if (flin) for (x=flin; *x; x++) {
		outls[0] = alist_append(outls[0],frame_copy_all(*x));
		//********** renaming copied frames........ ???
	}
	return FRES_STOREFR|FRES_DELETEMARK|FRES_STOREREM;	/* of type FCRET_PARTLIST */
}

int	COMFUNCTION(frcom_mmove) {
	ematrix		**elin = (ematrix**)inls[0];
	ematrix		**el, *em;
	
	outls[0] = new_alist(10);
	for (el=elin; *el; el++) {
		em = ematrix_copy(*el);
		if (EMNAME(*el))  EMSETNAME(em,EMNAME(*el));
		outls[0] = alist_append(outls[0],em);
	}
	return FRES_STOREMX|FRES_STOREMXNM|FRES_DELETEMARK|FRES_STOREREM;/* of type FCRET_PARTLIST */
}


/**
 * This is a general manipulating function for subtrees in the frame tree.
 * All frames in the subtrees of the given frames are copied to one list
 * which is stored on the output.
**/

int	COMFUNCTION(frcom_flatten) {
	framestruc	**flin = (framestruc**)inls[0];
	framestruc	**x,**y, **tr;
	
	outls[0] = new_alist(10);
	if (flin) for (x=flin; *x; x++) {
		tr = frame_gettree_dep(*x);
		for (y=tr; y?*y:0; y++)
			outls[0] = alist_append(outls[0],frame_copy_all(*y));
	}
	return FRES_STOREFR|FRES_DELETEMARK|FRES_STOREREM;	/* of type FCRET_LIST */
}












/******************	Basic matrix commands	*******************/
/******************************************************************/



/**
 * This function is used as a general finishing to single-matrix manipulating functions.
 * Depending on rettype, it prints or stores the resulting matrix em.
**/

int	frcom_decide_printmx(ematrix *em, void **outls[], int rettype) {
	
	if (FCRET_PRINTOWN(rettype)) {	/* printed-only version of the command */
		if (frcom_verbose_out()>0)  EMATOUTPUT(em,"    ~\t");
		else  EMATOUTPUTS(em,"    ~\t");
		dispose_ematrix(em);
	} else {
		outls[0] = alist_append(outls[0],em);
	}
	return (FCRET_PRINTOWN(rettype)?0: FRES_STOREMX|FRES_STOREREM);
}


/**
 * This is a handle for matrix pivoting; "!pivot row col matrix".
 * Only one input matrix is taken, it is copied and then the pivot is applied.
 * See frcom_decide_printmx() for final finishing...
**/

int	COMFUNCTION(frcom_pivot) {
	ematrix		*ein = inls[0][0];
	int		r = inum[0], c = inum[1];
	
	if (!ein)  return 0;
	if (r>0 && r<=ROWSM(ein) && c>0 && c<=COLSM(ein) && SIGNM(ein,r-1,c-1)!=0) {
		ein = ematrix_copy(ein);
		ematrix_pivot(ein,r-1,c-1);
	} else {
		USERERROR("Invalid pivot on an entry %d,%d.",r,c);  EMATDEBUG(0,ein,"!\t");
		return -1;
	}
	return frcom_decide_printmx(ein,outls,rettype);	/* of type FCRET_REPLMX1 */
}


/**
 * These functions delete (contract) the line given by its id (not row or column!).
 * The result is pronted or stored by frcom_decide_printmx().
**/

int	COMFUNCTION(frcom_delete) {
	ematrix		*ein = inls[0][0];
	int		id = inum[0];
	
	if (!ein)  return 0;
	ein = ematrix_removematid_del(ein,id);
	return frcom_decide_printmx(ein,outls,rettype);	/* of type FCRET_REPLMX1 */
}

int	COMFUNCTION(frcom_contract) {
	ematrix		*ein = inls[0][0];
	int		id = inum[0];
	
	if (!ein)  return 0;
	ein = ematrix_removematid_contr(ein,id);
	return frcom_decide_printmx(ein,outls,rettype);	/* of type FCRET_REPLMX1 */
}


/**
 * These functions delete (contract) all elements of the given matrix one by one,
 * creating an output list of minors...
**/

int	COMFUNCTION(frcom_removeeach) {
	ematrix		**elin = (ematrix**)inls[0];
	int		wh = (inum?inum[0]:0);
	ematrix		**el, *ee;
	int		i,z;
	
	for (el=elin; el?*el:0; el++)	/* all given input matrices are processed here */
		for (i=0; i<ROWSM(*el)+COLSM(*el); i++)
		  for (z=0; z<2; z++)  if (wh==z || wh>=2) {
			//*** ee = ematrix_removemat(*el,z,i>=ROWSM(*el),i<ROWSM(*el)?i:i-ROWSM(*el));
			ee = ematrix_removemat_rcsum(*el,z,i);
			if (EMNAME(ee))  FREE(EMNAME(ee));
			EMSETNAME(ee,EMNAME(*el));
			outls[0] = alist_append(outls[0],ee);
		  }
	return FRES_STOREMX|FRES_STOREREM|FRES_DELETEMARK|FRES_STOREMXNM;
}				/* of type FCRET_LISTMX */


/**
 * This function simply transposes the given matrix(ces) in the list.
**/

int	COMFUNCTION(frcom_dual) {
	ematrix		**elin = (ematrix**)inls[0];
	ematrix		**el, *ee;
	
	for (el=elin; el?*el:0; el++) {		/* all given input matrices are transposed here */
		ee = ematrix_copydual(*el);
		outls[0] = alist_append(outls[0],ee);	/* the output list is then stored back to the frames */
	}
	return FRES_STOREMX|FRES_STOREREM;	/* of type FCRET_REPLMXL */
}


/**
 * This function is for switching the current pfield to the new by given name.
 * 
**/

int	COMFUNCTION(frcom_pfield) {
	char	*pftx = instr?instr[0]:"";
	int	r;
	
	r = pfield_switchto_ext(0,pftx);
	if (r<0)  {USERERROR("Failed to switch pfield to \"%s\", now in \"%s\"\n",pftx,pfield_curname());}
	else  OUTPUT("Switched arithmetics to a new pfield \"%s\"\n",pfield_curname());
	return 0;	/* no return values */
}

/**
 * This function is provided for importing matrices of the given frames by the
 * given pf translation.
 * The imported copies of matrices are returned back for storing in the frames.
**/

int	COMFUNCTION(frcom_import) {
	char		*trtx = instr?instr[0]:"";
	framestruc	**flin = (framestruc**)inls[0];
	framestruc	**x;
	ematrix		*ee;
	
	if (pfield_switchtrans(trtx)<0) {	/* first must select the given pf translation */
		USERERROR("Setting pfield translation \"%s\" failed!",trtx);
		return -1;
	}
	if (flin) for (x=flin; *x; x++) {	/* imports copies of matrices of all given frames */
		pfield_translation(&FRVALUEEXP(*x),&FRVALUESIG(*x));
		ee = frame_extractmatrix(*x);
		if (!ee)  continue;
		ematrix_import(ee);	/* (calls pfield_translation() for all entries inside) */
		outls[0] = alist_append(outls[0],ee);	/* the output list is then stored back to the frames */
		FRPFINDEX(*x) = pfield_curindex();
	}
	return FRES_STOREMX|FRES_STOREREM;	/* of type FCRET_REPLMXL */
}






/**
 * 
 * 
 * sample...
**/

int	COMFUNCTION(frcom_) {
	ematrix		*ein = inls[0][0];
	
	outls[0] = alist_append(outls[0],ein);
	return FRES_STOREMX;	/* of type FCRET_REPLMX1 */
}














/******************	Structural matrix tests	*******************/
/******************************************************************/



/**
 * This function is used as a general finishing to frame-list returning functions.
 * Depending on rettype, it ...
**/

int	frcom_decide_outlist(void **outls[], long outnum[], int rettype) {
	
	if (outnum[0]<0)	/* must give the default list value to process the yes/no list */
		outnum[0] = ((FCRET_ISYESNOY(rettype) || rettype==FCRET_YESNOPRINT)? ret_yes:ret_no);
	
	if (FCRET_ISYESNO(rettype) && !FCRET_PRINTREQ(rettype)) {    /* of type FCRET_YESNO? */
		if ((rettype)==FCRET_YESNOFR || (rettype)==FCRET_NOYESFR)
			return FRES_STOREFR|FRES_DELETEMARK|FRES_STOREREM;
		else  return FRES_STOREFR|FRES_STOREREM;
	} else  return 0;		/* of type FCRET_PRINTED? */
}


/**
 * This function is used to mark all frames in the given list for later use.
**/

int	COMFUNCTION(frcom_marklist) {
	framestruc	**flin = (framestruc**)inls[0];
	framestruc	**x;
	
	for (x=flin; x?*x:0; x++)
		outls[0] = alist_append(outls[0],&ret_yes);
	outnum[0] = ret_yes;
	return FRES_STOREFR|FRES_STOREREM;	/* of type FCRET_LREMFR */
}


/**
 * This function handles the command '!msize mat-list R S cmp [S]'.
 * The matrices in the input frame list are tested for their size and number of elements,
 * compared to the input numbers according to the "cmp" string "[<=>.]+".
 * Frames woth no matrices can be matched only by size equal to 0x0, not otherwise.
**/

int	COMFUNCTION(frcom_hasmsize) {
	framestruc	**flin = (framestruc**)inls[0];
	int		row = (inumn>0&&inum? inum[0]: 0),
			col = (inumn>1&&inum? inum[1]: 0),
			sum = (inumn>2&&inum? inum[2]: -1);
	char		*cmpi = instr?instr[0]:"==";
	char		cmp[5];
	framestruc	**x;
	ematrix		*xem;
	int		i,k, rm,cm, ret;
	
	cmp[0] = cmp[1] = cmp[2] = cmp[3] = 0;
	strncpy(cmp,cmpi,3);
	if (!cmp[0] || row<0)  cmp[0] = '.';
	if (!cmp[1] || col<0)  cmp[1] = '.';
	if (!cmp[2] || sum<0 || inumn<3)  cmp[2] = '.';
	
	for (x=flin,i=k=0; x?*x:0; x++,i++) {
		xem = FRMATRIX(*x);
		rm = (xem? ROWSM(xem):0);
		cm = (xem? COLSM(xem):0);
		ret = ((cmp[0]=='=' && rm==row) || (xem && (cmp[0]=='.' ||
			(cmp[0]=='<' && rm<row) || (cmp[0]=='>' && rm>row))));
		ret &= ((cmp[1]=='=' && cm==col) || (xem && (cmp[1]=='.' ||
			(cmp[1]=='<' && cm<col) || (cmp[1]=='>' && cm>col))));
		ret &= ((cmp[2]=='=' && rm+cm==sum) || cmp[2]=='.' || (xem &&
			((cmp[2]=='<' && rm+cm<sum) || (cmp[2]=='>' && rm+cm>sum))));
		if (ret)  k++;
		if (FCRET_ISYESNOPR(rettype))
			outls[0] = alist_append(outls[0],(ret?&ret_yes:&ret_no));
	}
	if (FCRET_PRINTOWN(rettype) || frcom_verbose_out()>1) if (i>1)
		OUTPUT("Total %d out of %d frames have matrices of required size (%c%d, %c%d, %c%d).\n",
				k,i,cmp[0],row,cmp[1],col,cmp[2],sum);
	return frcom_decide_outlist(outls,outnum,rettype);
}


/**
 * This function is handling fans in the given matrices (frames).
 * If inum[] param is given, then we ask for a fan of size inum[0], otherwise for
 * the maximal fan.
 * The question about a fan of size inum[0] may generate a yes/no list.
**/

int	COMFUNCTION(frcom_hasfan) {
	framestruc	**flin = (framestruc**)inls[0];
	int		fan = ((inumn>0&&inum?(inum[0]!=0):0)? inum[0]: -1);
	framestruc	**x;
	ematrix		*ee;
	int		i,k,r;
	
	for (x=flin,i=k=0; x?*x:0; x++,i++) {
		ee = frame_extractmatrix(*x);
		if (!ee && FCRET_ISYESNOPR(rettype))  outls[0] = alist_append(outls[0],&ret_no);
		if (!ee)  continue;
		if (!FCRET_PRINTOWN(rettype) || frcom_verbose_out()<=1) {
			if (fan<0)  r = struct_hasfan_max(ee);
			else  r = struct_hasfan_full(ee,fan);
		} else if (frcom_verbose_out()>2) {
			r = struct_hasfan_printall(ee,fan);
		} else {
			r = struct_hasfan_print(ee,fan);
		}
		k += (r>0);
		if (FCRET_ISYESNOPR(rettype))
			outls[0] = alist_append(outls[0],(r>0?&ret_yes:&ret_no));
		if (FCRET_PRINTOWN(rettype) || frcom_verbose_out()>2)
			OUTPUT("The #%d matroid [%.20s] %s fan (>=3) of size %d.\n",
					i+1,FRNAME(*x),r>0?"+HAS+":"has -NO-",fan>0?fan:r);
		dispose_ematrix(ee);
	}
	if (fan>0 && i>2) if (FCRET_PRINTOWN(rettype) || frcom_verbose_out()>1)
		OUTPUT("Total %d out of %d matroids have fans of size %d.\n",k,i,fan);
	return frcom_decide_outlist(outls,outnum,rettype);
}


/**
 * This function is computing connectivity of the given matrices (frames).
 * If inum[] param is given, then we ask for connectivity >=inum[0], otherwise for
 * the full connectivity value.
 * The question about connectivity >=inum[0] may generate a yes/no list.
 * Newly, the function also tests simple/cosimple matroids - by icon[0]="s"/"c".
**/

int	COMFUNCTION(frcom_isconnected) {
	framestruc	**flin = (framestruc**)inls[0];
	int		con = (inumn>0&&inum? (inum[0]>=2? inum[0]:2): -1);
	char		*icon = instr?instr[0]:"";
	framestruc	**x;
	ematrix		*ee;
	int		i,k,r, rc,ri;
	
	for (x=flin,i=k=0; x?*x:0; x++,i++) {
		ee = frame_extractmatrix(*x);
		if (!ee && FCRET_ISYESNOPR(rettype))  outls[0] = alist_append(outls[0],&ret_no);
		if (!ee)  continue;
		if (icon[0]=='s')  r = struct_issimple(ee);
		else if (icon[0]=='c')  r = struct_iscosimple(ee);
		else if (con<0)  r = struct_connectivity(ee);
		else if (icon[0]!='i')  r = struct_isconnected(ee,con);
		else { rc = struct_iconnectivity(ee,&ri);  r = (rc+(!!ri)>con); }
		k += (r>0);
		if (FCRET_ISYESNOPR(rettype))
			outls[0] = alist_append(outls[0],(r>0?&ret_yes:&ret_no));
		if (FCRET_PRINTOWN(rettype) || frcom_verbose_out()>2) {
		  if (icon[0]=='s' || icon[0]=='c')
			OUTPUT("The #%d matroid [%.20s] is%s (%d) %ssimple.\n",i+1,FRNAME(*x),
				r>0?"":" -NOT-",r>0?0:(icon[0]=='c'?ROWSID(ee,-r-1):COLSID(ee,-r-1)),
				(icon[0]=='s'?"":"co"));
		  else
			OUTPUT("The #%d matroid [%.20s] has%s connectivity %s%s%d.\n",
				i+1,FRNAME(*x),r>0?"":" -NOT-",con>=0?"at least ":"exactly ",
				(icon[0]=='i'?"int-":""),(con>=0? con+(icon[0]=='i'):r));
		}
		dispose_ematrix(ee);
	}
	if (con>0 && i>2) if (FCRET_PRINTOWN(rettype) || frcom_verbose_out()>1) {
	  if (icon[0]=='s' || icon[0]=='c')
		OUTPUT("Total %d out of %d matroids are %ssimple.\n",k,i,(icon[0]=='s'?"":"co"));
	  else
		OUTPUT("Total %d out of %d matroids have connectivity at least %s%d.\n",
				k,i,(icon[0]=='i'?"int-":""),con+(icon[0]=='i'));
	}
	return frcom_decide_outlist(outls,outnum,rettype);
}


/**
 * This command handle computes the girth of given matroids (inum[0]<0),
 * or it decides whether the girth is at least inum[0].
 * All short cycles are possibly printed out (regardless whether the asked
 * girth has been achieved).
 * The question about the girth at least inum[0] may generate a yes/no list.
 * With the special value inum[0]==11111, the function tests paving matroids
 * (girth>=rank).
**/

int	COMFUNCTION(frcom_hasgirth) {
	framestruc	**flin = (framestruc**)inls[0];
	int		gt = ((inumn>0&&inum?(inum[0]!=0):0)? inum[0]: -1);
	framestruc	**x;
	ematrix		*ee, **gl=NULL,**y;
	int		i,j,k,r, gtx=gt;
	
	if (frcom_verbose_out()>1+(!FCRET_PRINTOWN(rettype)) && alist_getlength(flin)>1) {
		if (gtx!=11111)  OUTPUT("Computing girth (shortest cycle) >=%d in %d matroids.\n",gt,alist_getlength(flin));
		else  OUTPUT("Testing paving matroid (girth>=rank) in %d matroids.\n",alist_getlength(flin));
	}
	for (x=flin,i=k=0; x?*x:0; x++,i++) {
		ee = frame_extractmatrix(*x);
		if (!ee && FCRET_ISYESNOPR(rettype))  outls[0] = alist_append(outls[0],&ret_no);
		if (!ee)  continue;
		if (gtx==11111)  gt = ROWSM(ee);	/* (asking for paving matroid) */
		if (!FCRET_PRINTOWN(rettype) || frcom_verbose_out()<=1) {
			if (gtx<0)  r = struct_matgirth(ee);
			else  r = struct_hasmatgirth(ee,gt);
		} else {
			OUTPUT("Testing girth (shortest cycle) >=%d in #%d matroid [%.20s]...\n",gt,i+1,FRNAME(*x));
			gl = NULL;
			r = struct_hasmatgirth_list(ee,gt,&gl);
			for (y=gl; (y?*y:0) && (r>0 || frcom_verbose_out()>2); y++) {
				SOUTPUT("    ~\t - %d-cycle \t(",ROWSM(*y)+COLSM(*y));
				for (j=0; j<ROWSM(*y); j++)  SOUTPUT(" %d,",ROWSID(*y,j));
				SOUTPUT(" ");
				for (j=0; j<COLSM(*y); j++)  SOUTPUT(" %d,",COLSID(*y,j));
				SOUTPUT(" )\n");
				if (frcom_verbose_out()<=2) break;
			}
			if (gl)  dispose_alist_mats(gl);
		}
		k += (r>0);
		if (FCRET_PRINTOWN(rettype) || frcom_verbose_out()>2) {
			if (gtx!=11111)  OUTPUT("The #%d matroid [%.20s] %s girth %s%d.\n",
						i+1,FRNAME(*x),r>0?"+HAS+":"has -NOT- at least",gtx>=0&&r>0?">=":"",r>0?r:gt);
			else  OUTPUT("The #%d matroid [%.20s] %s paving (girth>=rank).\n",
						i+1,FRNAME(*x),r>0?"+IS+":"is -NOT-");
		}
		if (FCRET_ISYESNOPR(rettype))
			outls[0] = alist_append(outls[0],(r>0?&ret_yes:&ret_no));
		dispose_ematrix(ee);
	}
	if (gt>0 && i>2) if (FCRET_PRINTOWN(rettype) || frcom_verbose_out()>1)
		OUTPUT("Total %d out of %d matroids have girth at least %s%d.\n",
				k,i,gtx==11111?"rank (paving)":"",gtx==11111?0:gtx);
	return frcom_decide_outlist(outls,outnum,rettype);
}


/**
 * This is a handling function for matroid minors.
 * It looks whether the given frame(s) has a minor in the given minor-list inls[1].
 * Moreover, if the minor-list contains an equal matrix, during printing(!) command only,
 * this is reported and the minor removed.
**/

int	COMFUNCTION(frcom_minorlist) {
	framestruc	**flin = (framestruc**)inls[0];
	ematrix		**elmin = inls[1]? alist_copy((ematrix**)inls[1]):NULL;
	ematrix		*ee,*em, **el, **y;
	framestruc	**x;
	int		i,k,p,r;
	
	if (alist_getlength(flin)>2) if (FCRET_PRINTOWN(rettype) || frcom_verbose_out()>2)
		OUTPUT("Testing %d matroids for a list of {%.7s %.7s%s}%d minors.\n",alist_getlength(flin),
				(elmin?elmin[0]:0)?EMNAME(elmin[0]):"",(elmin?elmin[0]&&elmin[1]:0)?EMNAME(elmin[1]):"",
				(elmin?elmin[0]&&elmin[1]&&elmin[2]:0)?"..":"",alist_getlength(elmin));
	for (x=flin,i=k=0; x?*x:0; x++) {
		ee = frame_extractmatrix(*x);
		if (!ee && FCRET_ISYESNOPR(rettype))  outls[0] = alist_append(outls[0],&ret_no);
		if (!ee)  continue;
		i++;  el = elmin;
		if (FCRET_PRINTOWN(rettype))	/* looking if there is one equal in the minors */
		  for (y=elmin; y?*y:0; y++)
			if (ematrix_isequal(*y,ee)) {
				OUTPUT("The #%d matroid [%.20s] is +EQUAL+ to the minor #%d in the minor list, removing.\n",
					i,FRNAME(*x),(int)(y-elmin)+1);
				el = alist_getstart(alist_delete(alist_copy(y)));
				break;	/* an equal matrix is removed from the minor-list for printing */
			}
						/* testing for the minors is called here: */
		if (el?!el[0]:1)  r = 0;
		else if (!FCRET_PRINTOWN(rettype) || frcom_verbose_out()<=1) {
			r = struct_hasminorlist(ee,el);
		} else {
			em = frame_extractmatrix(*x);
			for (y=el, r=0; (y?*y:0); y++) {
				if (frcom_verbose_out()>2)  p = struct_allminors_printall(em,*y,1);
				else  p = struct_allminors_print(em,*y,1);
				if (p>=0)  r = (y-el)+1;
				if (r>0 && frcom_verbose_out()<=2)  break;
			}
			dispose_ematrix(em);
		}
		if (FCRET_ISYESNOPR(rettype))
			outls[0] = alist_append(outls[0],(r?&ret_yes:&ret_no));
		k += (r>0);
		if ((FCRET_PRINTOWN(rettype) || frcom_verbose_out()>2) && (r || el==elmin))
			OUTPUT("%sThe #%d matroid [%.20s] %s minor #%d [%.20s] in the %slist {%.7s %.7s%s}.\n",el!=elmin?"  ":"",
				i,FRNAME(*x),r?"+HAS+":"has -NO-",r,r?EMNAME(el[r-1]):"",el!=elmin?"reduced ":"",
				(el?el[0]:0)?EMNAME(el[0]):"",(el?el[0]&&el[1]:0)?EMNAME(el[1]):"",(el?el[0]&&el[1]&&el[2]:0)?"..":"");
		if (el!=elmin)  alist_free(el);
					/* a speed-up - move the found minor towards list start: */
		if (!FCRET_PRINTOWN(rettype) || frcom_verbose_out()<=0)
			if (el==elmin && r>=2) {
				em = el[r-1];  el[r-1] = el[2*r/3-1];  el[2*r/3-1] = em;
		}
	}
	if (i>2) if (FCRET_PRINTOWN(rettype) || frcom_verbose_out()>1)
		OUTPUT("Total %d out of %d matroids have minors in the list of length %d.\n",k,i,alist_getlength(elmin));
	if (elmin)  alist_free(elmin);	/* (a copy of the input list is used here!) */
	return frcom_decide_outlist(outls,outnum,rettype);
}


/**
 * This is a handling function for matroid minors using a specified element (by lineid).
 * It looks whether the given frame(s) has a minor in the given minor-list inls[1],
 * with the element id in 
**/

int	COMFUNCTION(frcom_minorusing) {
	framestruc	**flin = (framestruc**)inls[0];
	ematrix		*ein = inls[1]? ((ematrix**)inls[1])[0]:NULL;
	int		id = inum[0];
	ematrix		*ee, **el,**emout;
	framestruc	**x;
	int		i,k,p,r,j, rom,com;
	
	if (!ein)  return 0;
	if (alist_getlength(flin)>2) if (FCRET_PRINTOWN(rettype) || frcom_verbose_out()>2)
		OUTPUT("Testing %d matroids for a minor [%.12s] using element id%d.\n",alist_getlength(flin),EMNAME(ein),id);
	rom = ROWSM(ein);  com = COLSM(ein);
	
	for (x=flin,i=k=0; x?*x:0; x++) {
		ee = frame_extractmatrix(*x);
		for (j=0; j<ROWSM(ee)+COLSM(ee); j++)
			if (id==(j<ROWSM(ee)?ROWSID(ee,j):COLSID(ee,j-ROWSM(ee))))  break;
		if (j>=ROWSM(ee)+COLSM(ee)) {
			ee = NULL;
			if (FCRET_PRINTOWN(rettype) || frcom_verbose_out()>1)
				OUTPUT("The matroid [%.20s] has -NO- element of id%d!\n",FRNAME(*x),id);
		}
		if (!ee && FCRET_ISYESNOPR(rettype))  outls[0] = alist_append(outls[0],&ret_no);
		if (!ee)  continue;
			/**
			 * Here we find all all diplayed ein minors in ee (with self-maps
			 * factored out), and look for those using the element id.
			 * 
			**/
		i++;  emout = NULL;
		p = struct_allminors(ee,ein,1,&emout);
		if (p<0 || !emout) {
			if (FCRET_PRINTOWN(rettype) || frcom_verbose_out()>2)
				OUTPUT("The #%d matroid [%.20s] has -NO- [%.7s] minor at all.\n",i,FRNAME(*x),EMNAME(ein));
			continue;
		}
			/* the list contains also copies of bases, so look only at refmatrices! */
		for (el=emout, r=0; el?*el:0; el++) if (ISREFMAT(*el)) {
			for (j=0; j<rom+com; j++)
				if (id==(j<rom?ROWSID(*el,j):COLSID(*el,j-rom)))  break;
			if (j>=rom+com)  continue;
			r = 1;  break;
		}
		if (FCRET_PRINTOWN(rettype) || frcom_verbose_out()>2)
			OUTPUT("The #%d matroid [%.20s] %s [%.7s] minor using id%d.\n",
					i,FRNAME(*x),r?"+HAS+":"has -NO-",EMNAME(ein),id);
		if (FCRET_ISYESNOPR(rettype))
			outls[0] = alist_append(outls[0],(r?&ret_yes:&ret_no));
		k += (r>0);
		dispose_alist_mats(emout);  emout = NULL;
	}
	if (i>2) if (FCRET_PRINTOWN(rettype) || frcom_verbose_out()>1)
		OUTPUT("Total %d out of %d matroids have %.7s minors using id%d.\n",k,i,EMNAME(ein),id);
	return frcom_decide_outlist(outls,outnum,rettype);
}


/**
 * This function tests equivalence of the matroids in the first list to some members
 * of the second list.
 * If only one list is given, then equivalent pairs in this list are found, i.e.
 * each one compared against all further ones.
 * 
 * In the first case, if two equal matrices are found, then this is printed out separately.
 * (However, equality is considered as equivalence for purposes of filtering.)
**/

int	COMFUNCTION(frcom_equivalent) {
	framestruc	**flin = (framestruc**)inls[0];
	ematrix		**eleq = (ematrix**)inls[1];
	framestruc	**x;
	ematrix		**eleq2, **y, *ee,*ef;
	int		i,j,k, r,s,q, *outq=NULL;
	
				/* this is called for one given matrix list: */
	if (!eleq || alist_getlength(eleq)==0) {
		eleq2 = NULL;	/* extracting the matrices from the frame list, and preparing printlevel */
		for (x=flin; x?*x:0; x++) if (FRMATRIX(*x))
			eleq2 = alist_append(eleq2,frame_extractmatrix(*x));
				/* determining output verbosity (see inside struct_equivlistself_()) */
		i =  frcom_verbose_out();  j = FCRET_PRINTOWN(rettype)?1:0;
		s = (i>3? 5: (i>3-2*j? 3: (i>2-2*j? 2: (i>1-2*j? 1:0))));
		if (eleq2)  struct_equivlistself_ch(s,eleq2,&outq);
		
				/* converting the computation results */
		if (FCRET_ISYESNOPR(rettype))
			for (i=0; i<alist_getlength(eleq2); i++)
				outls[0] = alist_append(outls[0],(outq[i]<i?&ret_yes:&ret_no));
		if (outq)  FREE(outq);
		if (eleq2)  dispose_alist_mats(eleq2);
		
				/* this is called for two given matrix lists: */
	} else { for (x=flin,i=k=0; x?*x:0; x++,i++) {
		if (!FRMATRIX(*x) && FCRET_ISYESNOPR(rettype))  outls[0] = alist_append(outls[0],&ret_no);
		if (!FRMATRIX(*x))  continue;
		r = q = 0;  ee = NULL;
		for (y=eleq; (y?*y:0); y++) {
				/* looking if there is one equal in the second list for printing out */
			if (FCRET_PRINTOWN(rettype)) if (ematrix_isequal(*y,FRMATRIX(*x))) {
				if (frcom_verbose_out()>3)  EMATOUTPUTS(*y,"    ~\t\t");
				if (frcom_verbose_out()>1 || !q)
					OUTPUT("%sThe #%d matroid [%.12s] is +EQUAL+ to the #%d matroid [%.12s] in the second list.\n",
							q?"  ":"",i+1,FRNAME(*x),(int)(y-eleq)+1,EMNAME(*y));
				q = 1;  continue;
			}
			j = 0;		/* looking for equivalence here, one matrix by one */
			ef = frame_extractmatrix(*x);
			s = (frcom_verbose_out()>3? frcom_verbose_out()-1:0);
			if (struct_isequivalent_ch(s,ef,*y)>0) {
				j = r = (y-eleq)+1;  ee = *y;
			}
			if (ef)  dispose_ematrix(ef);
			if (r && frcom_verbose_out()<=2)  break;
			if (j)  OUTPUT("  The #%d matroid [%.12s] IS equivalent to the next #%d [%.12s].\n",
						i+1,FRNAME(*x),j,ee?EMNAME(ee):"");
		}
		if (FCRET_ISYESNOPR(rettype))
			outls[0] = alist_append(outls[0], ((r||q)?&ret_yes:&ret_no));
		k += (r||q);
		if ((FCRET_PRINTOWN(rettype) || frcom_verbose_out()>2) && (r || !q))
			OUTPUT("%sThe #%d matroid [%.20s] %s equivalent to #%d [%.20s] in the %slist {l %d}.\n",
				(q&&!r)?"  ":"",i+1,FRNAME(*x),r?"+IS+":"is -NOT-",r,r&&ee?EMNAME(ee):"",q?"reduced ":"",alist_getlength(eleq));
	}
		if (i>2) if (FCRET_PRINTOWN(rettype) || frcom_verbose_out()>1)
			OUTPUT("Total %d out of %d matroids have equivalents in the list {%.8s %.8s}%d.\n",k,i,
					(eleq?eleq[0]:0)?EMNAME(eleq[0]):"",(eleq?eleq[0]&&eleq[1]:0)?EMNAME(eleq[1]):"",alist_getlength(eleq));
	}
	return frcom_decide_outlist(outls,outnum,rettype);
}


/**
 * This is a handling function for testing the input frames as tight majors of the next
 * given minor list.
**/

int	COMFUNCTION(frcom_tmajorlist) {
	framestruc	**flin = (framestruc**)inls[0];
	ematrix		**elmin = (ematrix**)inls[1];
	framestruc	**x;
	ematrix		*ee;
	int		i,k,r,rr, id;
	
	if (alist_getlength(flin)>2) if (FCRET_PRINTOWN(rettype) || frcom_verbose_out()>2)
		OUTPUT("Testing %d matroids as tight majors of a list {%.8s %.8s%s}%d.\n",alist_getlength(flin),
				(elmin?elmin[0]:0)?EMNAME(elmin[0]):"",(elmin?elmin[0]&&elmin[1]:0)?EMNAME(elmin[1]):"",
				(elmin?elmin[0]&&elmin[1]&&elmin[2]:0)?"..":"",alist_getlength(elmin));
	for (x=flin,i=k=0; x?*x:0; x++,i++) {
		ee = frame_extractmatrix(*x);
		if (!ee && FCRET_ISYESNOPR(rettype))  outls[0] = alist_append(outls[0],&ret_no);
		if (!ee)  continue;
		if (frcom_verbose_out()>=2) {
			rr = struct_hasminorlist(ee,elmin);
			if (!rr)  OUTPUT("The #%d matroid [%.20s] has !NO! minor in the given \"tight-major\" list.\n",i+1,FRNAME(*x));
		} else  rr = 1;
		if (!FCRET_PRINTOWN(rettype) || frcom_verbose_out()<=1) {
			r = struct_hasdelcontr(ee,elmin);
		} else if (frcom_verbose_out()==2) {
			r = struct_hasdelcontr_print(ee,elmin);
		} else {
			r = struct_hasdelcontr_printall(ee,elmin);
		}
		id = !r?0: (r-1<ROWSM(ee)? ROWSID(ee,r-1): COLSID(ee,r-1-ROWSM(ee)));
		
		if (FCRET_ISYESNOPR(rettype))
			outls[0] = alist_append(outls[0],(!r?&ret_yes:&ret_no));
		k += (!r);
		if (FCRET_PRINTOWN(rettype) || frcom_verbose_out()>2)
			OUTPUT("%sThe #%d matroid [%.20s] %s a tight major (%s%d) of the given list {len %d}.\n",
					rr?"":"  ",i+1,FRNAME(*x),!r?"+IS+":"is -NOT-",r?"rem id ":"",id,alist_getlength(elmin));
		dispose_ematrix(ee);
	}
	if (i>2) if (FCRET_PRINTOWN(rettype) || frcom_verbose_out()>1)
		OUTPUT("Total %d out of %d matroids are tight majors of the list of length %d.\n",k,i,alist_getlength(elmin));
	return frcom_decide_outlist(outls,outnum,rettype);
}


/**
 * This is a handling function for testing proper matrix representation over the pfield.
 * (If some matrices are not properly represented, use option @nopfcheck...)
**/

int	COMFUNCTION(frcom_pfreprlist) {
	framestruc	**flin = (framestruc**)inls[0];
	framestruc	**x;
	ematrix		*em;
	int		i,r;
	
	for (x=flin,i=0; x?*x:0; x++,i++) {
		em = frame_extractmatrix(*x);
		r = ematrix_inpfield(em)>=0;
		if (!r && FCRET_PRINTOWN(rettype) && frcom_verbose_out()>1)
			ematrix_inpfield_printed(em);
		
		if (FCRET_ISYESNOPR(rettype))
			outls[0] = alist_append(outls[0],(r?&ret_yes:&ret_no));
		if (FCRET_PRINTOWN(rettype) || frcom_verbose_out()>1)
			OUTPUT("The #%d matroid [%.20s] %s properly represented over %s.\n",
				i+1,FRNAME(*x),r?"+IS+":"is -NOT-",pfield_curname());
		dispose_ematrix(em);
	}
	return frcom_decide_outlist(outls,outnum,rettype);
}


/**
 * This is a handling function for testing branchwidth3.
 * So far it works only for 3-connected matroids!!!
 * 
**/

int	COMFUNCTION(frcom_bwidth3test) {
	framestruc	**flin = (framestruc**)inls[0];
	framestruc	**x;
	ematrix		*em;
	int		i,k,b;
	char		*bo,**bop;
	
	bop = (frcom_verbose_out()>1 ?&bo:NULL);
	for (x=flin,i=k=0; x?*x:0; x++,i++) {
		em = frame_extractmatrix(*x);
		if (!em && FCRET_ISYESNOPR(rettype))  outls[0] = alist_append(outls[0],&ret_no);
		if (!em)  continue;
		if (bop)  *bop = NULL;
		if (ROWSM(em)+COLSM(em)<10 || (ROWSM(em)+COLSM(em)<14 && IFRANDDEBUGLESS(22)) || IFRANDDEBUGLESS(222))
			if (!struct_isconnected(em,3))  {USERERROR("We can correctly test branch-width 3 only for 3-connected matroids! [%s]",FRNAME(*x)); continue;}
		
		b = struct_hasbwidth3_sprint(em,bop);
		k += (b!=0);
		if (FCRET_ISYESNOPR(rettype))
			outls[0] = alist_append(outls[0],(b?&ret_yes:&ret_no));
		if (FCRET_PRINTOWN(rettype) && frcom_verbose_out()>2) {
			if (x!=flin) SOUTPUT("\n");  EMATOUTPUTS(em,"    ~\t");
		}
		if (FCRET_PRINTOWN(rettype) || frcom_verbose_out()>1)
			OUTPUT("The #%d matroid [%.20s] %s brach-width <=3.\n",i+1,FRNAME(*x),b?"+HAS+":"has -NOT-");
		if (bop?*bop:0) {
			OUTPUT(" Decomposition:  %s\n",*bop);  FREE(*bop);
		}
		dispose_ematrix(em);
	}
	if (i>2) if (FCRET_PRINTOWN(rettype) || frcom_verbose_out()>1)
		OUTPUT("Total %d out of %d matroids have brach-width <=3.\n",k,i);
	return frcom_decide_outlist(outls,outnum,rettype);
}


/**
 * This special handling function is provided for "cleaning" lists of matroids from
 * isomorphic (non-equivalent) pairs.
 * The idea is that the first given list is that one to clean, and the second given list
 * contains alternative representations of the matroids from the first list (in one-to-one
 * correspondence).
 * Only those matroids of the first list are kept that have no equivalent one further
 * down in any of the two lists.
**/

int	COMFUNCTION(frcom_makeunique) {
	framestruc	**flin = (framestruc**)inls[0];
	ematrix		**ealt = (ematrix**)inls[1];
	framestruc	**x;
	ematrix		*ee, **y,**z,**le;
	int		i,j, r,er;
	
	if (!flin || !ealt)  {USERERROR("Both lists must be given here!"); return 0;}
	if (alist_getlength(flin)!=alist_getlength(ealt))  {USERERROR("The given lists are supposed to have equal length!"); return 0;}
	er = 0;  le = NULL;
			/* preparing the corresponding pairs into one list */
	for (x=flin,i=0; x?*x:0; x++,i++) {
		ee = frame_extractmatrix(*x);
		if (!ee)  continue;
		if (ROWSM(ee)+COLSM(ee)<14 && IFRANDDEBUGLESS(11))
		  if (ematrix_getmhash(ee)!=ematrix_getmhash(ealt[i]))
			if (!er)  {er = 1; USERERROR("The given lists are supposed to contain pairs of isomorphic matroids!");}
		le = alist_append(le,ee);
		ee = ematrix_copy(ealt[i]);
		if (EMNAME(ealt[i]))  EMSETNAME(ee,EMNAME(ealt[i]));
		le = alist_append(le,ee);
	}
	if (!er && frcom_verbose_out()>(FCRET_PRINTOWN(rettype)?0:2))
		OUTPUT("Looking for isomorphic pairs in a double list of length %d:\n",alist_getlength(flin));
	
			/* processing the prepared list - first of each pair against all remaining pairs */
	for (y=le,i=0; (y?y[0]&&y[1]:0) && !er; y+=2,i++) {
		ee = *y;
		for (z=y+2,r=0; *z && !r; z++)
			r = struct_isequivalent(ee,*z)? 1:r;
		
		if (FCRET_ISYESNOPR(rettype))
			outls[0] = alist_append(outls[0],(r?&ret_no:&ret_yes));
		j = (int)(z-y-1)/2+i;
		if (r && FCRET_PRINTOWN(rettype) && frcom_verbose_out()>1) {
			EMATOUTPUTS(ee,"    ~\t\t ");  EMATOUTPUTS(z[-1],"    ~\t\t~");
		}
		if (r && (FCRET_PRINTOWN(rettype) || frcom_verbose_out()>1))
			OUTPUT("The #%d matroid [%s] is isomorphic to the #%d matroid [%s]%s.\n",
					i+1,FRNAME(flin[i]),j+1,FRNAME(flin[j]),FCRET_ISYESNOPR(rettype)?", REMOVED":"");
		if (!r && frcom_verbose_out()>(FCRET_PRINTOWN(rettype)?0:2))
			OUTPUT("The #%d matroid [%s] is unique in the lists.\n",i+1,FRNAME(flin[i]));
	}
	dispose_alist_mats(le);
	return frcom_decide_outlist(outls,outnum,rettype);
}


/**
 * This function is looking for matroids of given hash value in the list.
 * (A matroid hash value is computed from the basis-structure of the matrix,
 * and hence it is matroid-invariant and may be used to informally compare
 * two matroids over different pfields.)
 * The function may generate a yes/no list for filtering.
**/

int	COMFUNCTION(frcom_matroidhash) {
	framestruc	**flin = (framestruc**)inls[0];
	long		hash = ((inumn>0&&inum?(inum[0]!=0):0)? inum[0]: -1);
	framestruc	**x;
	ematrix		*ee;
	int		i,k,r;	long	hh;
	
	for (x=flin,i=k=0; x?*x:0; x++,i++) {
		ee = frame_extractmatrix(*x);
		if (!ee && FCRET_ISYESNOPR(rettype))  outls[0] = alist_append(outls[0],&ret_no);
		if (!ee)  continue;
		hh = ematrix_getmhash(ee);
		r = (hh==hash? 1:-1);  k += (r>0);
		
		if (FCRET_ISYESNOPR(rettype))
			outls[0] = alist_append(outls[0],(r>0?&ret_yes:&ret_no));
		if (FCRET_PRINTOWN(rettype)) {
			if (frcom_verbose_out()>1)  EMATOUTPUTS(ee,"    ~\t\t");
			if (hash==-1)  OUTPUT("The #%d matroid [%.20s] has hash-value \t%ld.\n",i+1,FRNAME(*x),hh);
		}
		if (hash!=-1 && frcom_verbose_out()>(FCRET_PRINTOWN(rettype)?-1:(r>0?1:2)))
			OUTPUT("The #%d matroid [%.20s] %s matroid hash-value \t%ld %s %ld.\n",
					i+1,FRNAME(*x),r>0?" +HAS+ ":"has NOT",hh,r>0?"==":"!=",hash);
	}
	if (i>2 && hash!=-1) if (FCRET_PRINTOWN(rettype) || frcom_verbose_out()>1)
		OUTPUT("Total %d out of %d matroids have hash-value %ld.\n",k,i,hash);
	return frcom_decide_outlist(outls,outnum,rettype);
}


/**
 * This function tests for matroid isomorphism in one or two given matroid
 * lists -- similarly as in frcom_equivalent...
 * Not much is here, the whole testing is implemented in a separate
 * function strmag_isomorphlist_().
 * Notice that the input frames to this function may be over different
 * pfields, and this function can handle that.
**/

int	COMFUNCTION(frcom_isomorph) {
	framestruc	**flin = (framestruc**)inls[0],
			**flin2 = (framestruc**)inls[1];
	framestruc	**x;
	ematrix		*ee, **el[2]={NULL,NULL};
	int		i,s, *xpf[2]={NULL,NULL}, *outq=NULL;
	
	xpf[0] = MMALLOC((alist_getlength(flin)+alist_getlength(flin2)+2)*sizeof(xpf[0][0]));
	xpf[1] = xpf[0]+(alist_getlength(flin)+1);
		/* extracting the input matrix lists and their pfields */
	for (i=0; i<2; i++)
	  for (x=i?flin2:flin; x?*x:0; x++) {
		ee = frame_extractmatrix(*x);
		if (!ee)  continue;
		el[i] = alist_append(el[i],ee);
		xpf[i][alist_getlength(el[i])-1] = FRPFINDEX(*x);
	}
		/* preparing verbosity, and calling the isomorphism test over el[0],el[1] */
	if (FCRET_PRINTOWN(rettype))  s = frcom_verbose_out()+1+(!flin2);
	else  s = frcom_verbose_out()+ (frcom_verbose_out()>1? (!flin2):0);
	if (alist_getlength(flin)==1 && s>1)  s++;
		/* (all output is printed out inside strmag_isomorphlist_()...) */
	if (el[0])  strmag_isomorphlist_ch(s,el[0],xpf[0],el[1],xpf[1],&outq);
	if (FCRET_ISYESNOPR(rettype))
		for (i=0; i<alist_getlength(el[0]); i++)
			outls[0] = alist_append(outls[0], \
				((el[1]?outq[i]>=0:outq[i]<i)? &ret_yes:&ret_no));
	
	if (outq)  FREE(outq);  if (xpf[0])  FREE(xpf[0]);
	if (el[0])  dispose_alist_mats(el[0]);
	if (el[1])  dispose_alist_mats(el[1]);
	return frcom_decide_outlist(outls,outnum,rettype);
}














/******************	Matrix constructions	*******************/
/******************************************************************/



/**
 * This is a handle for matrix self- line maps -- the maps are printed on the screen.
 * Nothing is returned.
**/

int	COMFUNCTION(frcom_automorph) {
	ematrix		*ein = inls[0][0];
	int		p;
	
	if (!ein)  return 0;
	if (frcom_verbose_out()>0)  EMATOUTPUT(ein,printoutpref);
	if (frcom_verbose_out()>0)  p = struct_matrixselfmaps_print(ein);
	else  p = struct_matrixselfmaps_number(ein);
	OUTPUT("\tFound %d (scaled) self- line maps (\"automorphisms\") of the matrix [%s].\n",p,EMNAME(ein));
	return 0;	/* of type FCRET_PRINTED */
}



/**
 * This function generates 3-connected extensions of the matrix(ces) in the frame(s) flin.
 * The string parameter what tells which extensions should be generated - r for a row,
 * c for a column, b for both.
 * One may give a sequence like "rcbb..." to generate several steps at once.
 * Or, one may give, instead of what, two number parameters rrow,rcol to generate all
 * extensions up to the dimensions rrow,rcol of the resulting matrix.
 * 
 * If using any form of multi-step extension, the input may contain only one matrix.
 * The result is returned as a list of all generated matrices (not containing input ones),
 * or only those of dimensions rrow,rcol if this is requested.
 * Read ../gener/gener.c.
 * 
 * IMPORTANT: The generating process uses internal options @bsize and @signature
 * that are initialized at the first step of generating.
 * They should never be set by hand, as doing so may compromise the whole process.
 * All steps of generating extensions must be 3-connected, including the starting one!
 * Read carefully the theory in ../include/gener.h.
 * 
**/

int	COMFUNCTION(frcom_extend) {
	framestruc	**flin = (framestruc**)inls[0];
	char		*what = instr?instr[0]:NULL;
	int		rrow = inumn>=2?inum[0]:0,  rcol = inumn>=2?inum[1]:0;
	char		*whx, *upto="==";
	framestruc	**x,**y,**z, **out;
	ematrix		*xem;
	int		i,j,k,r, ext=0, cox=0;
	
	if (!what)  what = (rrow>0? "==":"b");
	if (rrow>0 && rcol>0) {
		upto = what;
		what = "nnnnnnnnnnnnnnnnnnnnnnnnn";
	}
	out = flin;  whx = what;
	for (; (out?*out:0) && what[0]; what++) {
		if (out!=flin && flin[1])  {USERERROR("Cannot correctly do multiple-steps with more than one input matrix!");}
		if (what[0]!='r' && what[0]!='c' && what[0]!='b' && what[0]!='n')  {USERERROR("Wrong extension description string \"%s\".",what);}
		if (frcom_verbose_out()>0 && (alist_getlength(out)>1||what[1]))
			OUTPUT("This generation step extends \"%c\" for %d sequences {%.7s %.7s%s}:\n",
					what[0],alist_getlength(out),FRNAME(out[0]),out[1]?FRNAME(out[1]):"",out[1]&&out[2]?"..":"");
		y = NULL;
		ext = (what[0]=='c' || what[0]=='b');	/* (extend column or both) */
		cox = (what[0]=='r' || what[0]=='b');	/* (extend row or both) */
		
		for (x=out,i=0; x?*x:0; x++,i++) {
			if (!FRMATRIX(*x))  continue;
			if (what[0]=='n') {	/* special - extending all to given size rrow x rcol */
				xem = FRMATRIX(*x);
				ext = (COLSM(xem)<rcol);  cox = (ROWSM(xem)<rrow);
				if ((!ext || upto[1]=='<') && (!cox || upto[0]=='<'))
					outls[0] = alist_append(outls[0],frame_copy_all(*x));
			}
				/* must (temporarily) give the parent for inheriting options! */
			if (out!=flin)  frame_setson(*x,flin[0]);
			for (j=0, k=cox; j<2; j++, k=ext) if (k) {
				z = NULL;
				if (frcom_verbose_out()<=0)  r = gener_extframe_nopr(*x,j,&z);
				else if (frcom_verbose_out()<=2)  r = gener_extframe(*x,j,&z);
				else  r = gener_extframe_print(*x,j,&z);
				if (frcom_verbose_out()>2)  SOUTPUT("\n");
				y = alist_applist(y,z);
			}
			if (out!=flin)  frame_setson(*x,NULL);
		}
				/* dispose intermediate list if generating to size, or store it if generating by symbols */
		if (out!=flin && what[0]=='n')  dispose_alist_frams(out);
		if (out!=flin && what[0]!='n')  outls[0] = alist_applist(outls[0],out);
		out = y;	/* (we cannot append the list y here since appending would destroy it!!!) */
	}
	
	if (out!=flin)  outls[0] = alist_applist(outls[0],out);
	if (frcom_verbose_out()>0) if (alist_getlength(flin)>1 || what[1])
		OUTPUT("In total %d (co-)extensions of %d matrix-sequences generated for \"%s\" over %s.\n",
				alist_getlength(outls[0]),alist_getlength(flin),rrow>0?"size":whx,pfield_curname());
	if (!outls[0])  outls[0] = new_alist(4);
	return FRES_DELETEMARK|FRES_STOREFR|FRES_STOREREM;
}


/**
 * Here we test pfield representability / generate pfield representations of the given
 * matroids in the list (the pfield is given as the first parameter in the first case).
 * When generating representations, the first parameter tells to generate only one
 * (default), all ("all"), or all nonequivalent ("allq") representations for each.
 * Read theory in gener/grepres.c ...
 * For nonequivalent representations, struct_equivlistself_ch() is used as a filter.
 * 
 * The first function returns a yes-no printed filter of input frames,
 * the second function returns a list of generated matrix representations.
**/

int	COMFUNCTION(frcom_reprtest) {
	framestruc	**frin = (framestruc**)inls[0];
	char		*whatpf = instr?instr[0]:"";
	framestruc	**x;
	ematrix		*ee;
	int		k,p, xp;
	
	xp = pfield_curindex();
	if (whatpf[0])  pfield_switchto(whatpf);
	else  xp = -1;
	for (x=frin, k=0; x?*x:0; x++) {
		ee = frame_extractmatrix(*x);
		if (!ee && FCRET_ISYESNOPR(rettype))  outls[0] = alist_append(outls[0],&ret_no);
		if (!ee)  continue;
		p = grepr_isrepresented(ee,FRPFINDEX(*x));
		if (p>0)  k++;
		if (FCRET_PRINTOWN(rettype) || frcom_verbose_out()>1) {
			OUTPUT("There %s %s-representation of the #%d matroid [%.18s] (%dx%d, %s).\n",
				(p>0?"+IS+ a":"is -NO-"),pfield_curname(),(int)(x-frin)+1,FRNAME(*x),
				ROWSM(ee),COLSM(ee),pfield_ixname(FRPFINDEX(*x)));
		}
		if (FCRET_ISYESNOPR(rettype))
			outls[0] = alist_append(outls[0],(p>0?&ret_yes:&ret_no));
		dispose_ematrix(ee);
	}
	if (frcom_verbose_out()>0) if (alist_getlength(frin)>2)
		OUTPUT("In total %d of %d given matroids have %s-representations.\n",k,alist_getlength(frin),pfield_curname());
	if (xp>=0)  pfield_switchto_fast(xp);
	return frcom_decide_outlist(outls,outnum,rettype);
}

int	COMFUNCTION(frcom_reprgener) {
	framestruc	**frin = (framestruc**)inls[0];
	char		*how = instr?instr[0]:"";
	framestruc	**x;
	ematrix		*ee, **lout,**y;
	int		i,k,p, *eqo=NULL;
	
	for (x=frin; x?*x:0; x++) {
		ee = frame_extractmatrix(*x);
		if (!ee)  continue;
		lout = NULL;	/* generating representations first: */
		if (frcom_verbose_out()>3)
			p = grepr_generate_prall(ee,FRPFINDEX(*x),&lout);
		else if (strncmp(how,"all",3)==0)
			p = grepr_generate_all(ee,FRPFINDEX(*x),&lout);
		else
			p = grepr_generate(ee,FRPFINDEX(*x),&lout);
		
				/* processing generated list - possible equivalence testing: */
		if (strncmp(how,"allq",4)==0 && alist_getlength(lout)>0)
			struct_equivlistself_ch(0,lout,&eqo);
		else if (strncmp(how,"all",3)==0)  eqo = NULL;
		else  eqo = (void*)1;
		for (i=k=0, y=lout; y?*y:0; y++,i++) {
			if (EMNAME(*y))  FREE(EMNAME(*y));
			EMSETNAME(*y,EMNAME(ee));
			if (!eqo || (eqo==(void*)1 && i==0) || ((unsigned long)eqo>1? eqo[i]==i:0)) {
				k++;
				outls[0] = alist_append(outls[0],*y);
			} else  dispose_ematrix(*y);
		}
		if (k==0 && frcom_verbose_out()>0)
			OUTPUT("There is -NO- %s-representation of the #%d matroid [%.18s] (%dx%d, %s).\n",
					pfield_curname(),(int)(x-frin)+1,FRNAME(*x),ROWSM(ee),COLSM(ee),pfield_ixname(FRPFINDEX(*x)));
		if (frcom_verbose_out()>1 && eqo!=(void*)1 && k>0)
			OUTPUT("There are %d %s%s-representations of #%d matroid [%.18s] (%dx%d, %s).\n",
					k,eqo?"nonequiv ":"",pfield_curname(),
					(int)(x-frin)+1,FRNAME(*x),ROWSM(ee),COLSM(ee),pfield_ixname(FRPFINDEX(*x)));
		dispose_ematrix(ee);
		if (lout)  alist_free(lout);
		if (eqo && eqo!=(void*)1)  FREE(eqo);
	}
	if (frcom_verbose_out()>0) if (alist_getlength(frin)>1)
		OUTPUT("In total %d %s-representations (%s) of %d given matroids generated.\n",
				alist_getlength(outls[0]),pfield_curname(),how,alist_getlength(frin));
	if (!outls[0])  outls[0] = new_alist(4);
	return FRES_STOREMX|FRES_STOREREM|FRES_DELETEMARK|FRES_STOREMXNM;
}				/* of type FCRET_LISTMX */



















/******************	Command-flow Control	*******************/
/******************************************************************/



/**
 * This function requests to restart command processing in the whole tree.
 * The commands which were already processed are not repeated, but new commands
 * later added with '!append' get chance...
**/

int	COMFUNCTION(frcom_restart) {
	
	OUTPUT("\tRestarting command processing in the whole tree...\n");
	if (frcom_verbose_out()>0)  OUTPUT("\t (Commands already processed are not repeated.)\n");
	return FRES_RESTART;
}


/**
 * These command handles are for skipping and exitting from command processing.
 * 
**/

int	COMFUNCTION(frcom_comexit) {
	int	num = inumn>=1?inum[0]:0;
	
	if (frcom_verbose_out()>0)  OUTPUT("\tExit %d from the program command processing.\n",num);
	if (num<0)  DEBUG(CURDLEV-2,"Negative return value %d is returned as %d !\n",num,FRES_PROCSTRIP2(num));
	frame_skipcommands = -1;
	return FRES_PROCSTRIP2(num)|FRES_SKIPCOM;
}

int	COMFUNCTION(frcom_comskip) {
	int	num = inumn>=1?inum[0]:0;
	
	if (frcom_verbose_out()>0)  OUTPUT("\tSkipping the next %d commands in processing.\n",num);
	frame_skipcommands = num>0?num:0;
	return FRES_SKIPCOM;
}


/**
 * This is a command testing length of the given frame list.
 * If the test is true, then the next command is executed, otherwise it is skipped.
 * 
**/

int	COMFUNCTION(frcom_iflist) {
	framestruc	**flin = (framestruc**)inls[0];
	int		len = inumn>=1?inum[0]:0;
	char		*what = instr?instr[0]:"";
	int		r,ll;
	
	ll = alist_getlength(flin);
	r = 0;
	if (what[0]=='!' || strcmp(what,"<>")==0) {
		r = (len!=ll);
	} else if (what[0]=='=' && strcmp(what,"=!")!=0) {
		r = (len==ll);
	} else if (what[0]=='<') {
		r = what[1]=='='? (len<=ll): (len<ll);
	} else if (what[0]=='>') {
		r = what[1]=='='? (len>=ll): (len>ll);
	} else {
		USERERROR("Invalid comparing \"%s\" in !iflist.",what);
	}
	if (frcom_verbose_out()>0)  OUTPUT("If-list length test %d %s %d, result %s.\n",len,what,ll,r?"YES":"NO");
	frame_skipcommands = r?0:1;
	return (r?0:FRES_SKIPCOM)|FRES_STOREREM;
}

/**
 * This is a command testing the return status of the given shell command.
 * If the test is true, then the next command is executed, otherwise it is skipped.
 * 
**/

int	COMFUNCTION(frcom_ifshell) {
	char		*shx = instr?instr[0]:"false";
	int		rs = inumn>=1?inum[0]:0;
	int		f,w=0,r=0;
	
	if (safexec>0) {	/* this is not allowed in a safe mode... */
		USERERROR("NO shell command execution is allowed in a safe mode (-s or -S)!");
		return -1;
	}
	if (frcom_verbose_out()>1)  OUTPUT("Executing  \'%.60s\' ...\n",shx);
	f = fork();
	if (f==0) {
		if (execlp("/bin/sh","/bin/sh","-c",shx,NULL)<0)  exit(-1);
	} else {
		wait(&w);
		w = WIFEXITED(w)?WEXITSTATUS(w):-1;
	}
	r = (w==rs);
	if (frcom_verbose_out()>0)  OUTPUT("If-shell test \'%.28s\', return status %d %s.\n",shx,w,r?"YES":"NO");
	frame_skipcommands = r?0:1;
	return (r?0:FRES_SKIPCOM)|FRES_STOREREM;
}

/**
 * This is a command testing readability of the given file
 * If the test is true, then the next command is executed, otherwise it is skipped.
**/

int	COMFUNCTION(frcom_iffile) {
	char		*fnm = instr?instr[0]:NULL;
	int		r=0;
	FILE		*fp=NULL;
	
	if (fnm)  fp = fopen_path(fnm,"r",frame_path_read,NULL,frame_extension);
	r = (fp!=NULL);
	if (fp)  fclose(fp);
	if (frcom_verbose_out()>0 && fnm)
		OUTPUT("Given file \'%.60s\' is%s in the read-search path.\n",fnm,r?"":" -NOT-");
	frame_skipcommands = r?0:1;
	return (r?0:FRES_SKIPCOM)|FRES_STOREREM;
}






















int	COMFUNCTION(frcom_test) {
	//ematrix	*ein = inls[0][0];	ematrix		**exl;
	framestruc	**flin = (framestruc**)inls[0];
	framestruc	**x,**y=NULL;
	int	i;
	
	//exl = gener_matextens_ext(2,ein);
	//dispose_alist_mats(exl);
	
	if (flin) for (x=flin,i=0; *x; x++,i++) {
		//y = gener_extframe_ext(*x,0,0);
	}
	
	outls[0] = (void**)y;
	//return 0;
	//return FRES_STOREMX;
	return FRES_DELETEMARK|FRES_STOREFR;
}














#include "frcomst.c"



















