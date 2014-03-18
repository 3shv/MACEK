
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
 * Some theory about testing matroid representabilty for start...
 * We use a brute-force algorithm here that tries to fill nonzero entries of
 * a new matrix with pfield elements so that all subdeterminants agree with
 * the given matroid.
 * The implementation is short, using the gmatext.c code, but not much efficient.
 * (It would be good to write a separate code tailored to this problems from scratch.)
 * 
 * 
 * See also gmatext.c for related information and routines...
 * 
**/


#include "macek.h"
#include "gen.h"













/******************	Generating matroid representations	*******************/
/**********************************************************************************/


#define	CURDLEV		6		/*7 (use <=DEBUGLEV+2 to perform extra checks) */


/**
 * This function adds one line to the matrix ex so that the matroid of ex+ is isomorphic
 * to the matroid of given ed (represented over xpfd).
 * The matrix ex may be given NULL -- when ed has one row or one column.
 * All distinct extensions are stored to the output list *lout (if given).
 * The value of tr determines whether we are adding a row (0) or a column(1);
 * when tr=-1, then this is guessed automatically from the sizes of ex,ed.
 * The return value is the number of generated extensions, or 0, or <0 on an error.
**/

int	grepr_addline(ematrix *ex, ematrix *ed, int tr, int xpfd, ematrix ***lout) {
	ematrix		*ee,*ef, **el=NULL,**elo=NULL, **x;
	int		i,j,k,p,r, ro,co, xp=-1;
	
	if (!ed)  {PROGERROR("the matrix ed must be given!"); return -1;}
	ro = ROWSM(ed);  co = COLSM(ed);
	if (tr<0)  tr = ((ex? ROWSM(ex)<ro: ro<=1)? 0:1);	/* (tr need not be given - determining from matrix sizes) */
	if (ex? (ro!=ROWSM(ex)+(!tr) || co!=COLSM(ex)+(!!tr)): (ro!=1 && co!=1))
		{PROGERROR("the matrices ex, ed have wrong sizes! tr=%d; ro=%d, co=%d - %d, %d",tr,ro,co,ROWSM(ex),COLSM(ex)); return -1;}
	if ((ex?ISREFMAT(ex):0) || ISREFMAT(ed))  {PROGERROR("do not give refering matrices here"); return -1;}
	if (lout?*lout:0)  {PROGERROR("the list lout must be given empty (null)"); return -1;}
	DEBUG(CURDLEV+1,"Adding a %s to matrix %p to (%dx%d), pf=%d\n",tr?"column":"row",ex,ro,co,xpfd);
	if (ex && lout)  EMATDEBUGS(CURDLEV+2,ex,"\t\t-\t");
	
			/**
			 * We first try the easy extensions:
			 * If ed is 1xC or Rx1, then ex has entries 0 or 1 depending on ed.
			 * If some other line of ed is parallel (and nonzero) to the added line,
			 * then we make ex by copying the corresponding line to append.
			 * (We have to use the proper pfield xpfd for ed!)
			**/
	xp = pfield_curindex();
	if (xpfd>=0) pfield_switchto_fast(xpfd);
	if (lout)  EMATDEBUGS(CURDLEV+2,ed,"\t\t~\t");
	p = tr? ematrix_parallel_col(ed,co-1): ematrix_parallel_row(ed,ro-1);
	if (p>=0) if (ematrix_linezero(ed,tr,p)>0)  p = -1;
	if (xpfd>=0 && xp>=0) pfield_switchto_fast(xp);
	if (!ex) {
		ee = new_ematrix(ro,co,ro,co);
		for (i=0; i<ro; i++) for (j=0; j<co; j++) {
			SETEXSIGMZERO(ee,i,j);
			if (SIGNM(ed,i,j)!=0)  SIGNM(ee,i,j) = 1;
		}
		elo = alist_append(elo,ee);  k = 1;
		DEBUG(CURDLEV+1,"Creating a %s matrix with %s.\n",tr?"1xC":"Rx1",gener_extprintline(0,ee,tr));
	} else if (p>=0 && ro+co>=3) {
		ef = ematrix_copy_to(ex,ro,co);
		ef = ematrix_append_rc(ef,tr);
		if (tr) { for (i=0; i<ro; i++)
			if (SIGNM(ed,i,co-1))  COPYEXSIGM(ef,i,co-1,ef,i,p);
		} else { for (i=0; i<co; i++)
			if (SIGNM(ed,ro-1,i))  COPYEXSIGM(ef,ro-1,i,ef,p,i);
		}
		elo = alist_append(elo,ef);  k = 1;
		DEBUG(CURDLEV+1,"Appending a parallel %s %s.\n",tr?"column":"row",gener_extprintline(0,ef,tr));
	} else {
			/**
			 * In all other cases we use the next routine:
			 * We generate all one-line extensions of ex which have the same
			 * zero pattern in the appended line as ed.
			 * Then we check the subdeterminants of each extension against
			 * the matrix ed, and store copies of isomorphic ones.
			**/
		r = gener_matextens_zeros(ex,tr,ed,&el);
		DEBUG(CURDLEV+1,"Filtering %d (zero-respecting) extension lines...\n",alist_getlength(el));
		if (r<0)  return r;
		for (x=el, k=0; x?*x:0; x++)
			if (ematrix_havesamedets_repres(*x,ed,xpfd,tr)) {
				k++;
				elo = alist_append(elo,ematrix_copy(*x));
				DEBUG(CURDLEV+1+2*(k>2),"Appending a #%d(%d) %s %s.\n",k,alist_getlength(elo),tr?"column":"row",gener_extprintline(0,*x,tr));
			}
		if (el)  dispose_alist_mats(el);
	}
#ifndef FASTPROG	
	if (IFRANDDEBUGLESS(222) && k>0) {
		x = elo+RANDOM()%k;
		if (!ematrix_havesamedets(*x,ed,xpfd))  {PROGERROR("the extended matrix has not same subdeterminants!");}
	}
	if (IFRANDDEBUGLESS(222) && ro>2 && co>2 && lout) {
		ee = ematrix_copydual(ex);  ef = ematrix_copydual(ed);
		ematrix_swaprows(ee,0,co-2);  ematrix_swaprows(ef,0,co-2);
		DEBUG(CURDLEV+1,"Recursive testing the number of added lines...\n");
		r = grepr_addline(ee,ef,(tr<0?-1:!tr),xpfd,NULL);
		if (r!=k && (r>=0 || k>=0))  {PROGERROR("wrong result of a recursive call %d!=%d",r,k);}
	}
#endif
	if (lout)  *lout = (*lout? alist_applist(*lout,elo): elo);
	else if (elo)  dispose_alist_mats(elo);
	return k;
}



/**
 * Here we generate (labelled-)distinct unit-scaled representations of the
 * matroid eg (which is over xpfd) over the current pfield.
 * Or, we simply test representability if ch<=1.
 * The representations are returned through the give matrix list *lout.
 * (*lout must be initialized empty or NULL.)
 * 
 * The generated representations are collected in the list *lout (if given, ch>=2).
 * The function returns the number of generated representations, or -1 on an error.
 * If ch>=1, then the results are printed out.
 * If ch>=2, then the generated representations are all collected (and printed out ch>=4).
 * The return value is the number of generated representations (>0), or 0,
 * or -1 for an error.
**/

int	grepr_generate_ext(int ch, ematrix *eg, int xpfd, ematrix ***lout) {
	ematrix		*ee,*re,**rel=NULL, **elo=NULL,**x, ***epl=NULL;
	int		i,k,p,br, ro,co, step, xp=-1, ret=0;
	char		*s, buf[100];
	
	if (!eg)  {PROGERROR("the matrix eg must be given!"); return -1;}
	ro = ROWSM(eg);  co = COLSM(eg);
	if (ro<=0 || co<=0)  {PROGERROR("the matrix eg has wrong size!"); return -1;}
	if (lout?*lout:0)  {PROGERROR("the list lout must be given empty (null)"); return -1;}
	xp = pfield_curindex();
	s = pfield_curname();
	if (xpfd>=0) pfield_switchto_fast(xpfd);
	DEBUG(CURDLEV-1,"Generating representation(s) of a %s matrix %p[%.12s] (%dx%d) over %s...\n",
			pfield_curname(),eg,EMNAME(eg)?EMNAME(eg):"",ro,co,s);
	EMATDEBUG(CURDLEV+0,eg,"\t:\t");
	if (xpfd>=0 && xp>=0) pfield_switchto_fast(xp);
			/**
			 * We allocate working data, and we prepare the submatrices rel[]
			 * of the matrix eg which display how the generated representation
			 * grows line by line (total step=ro steps).
			 * We should find a more clever order of added lines for faster
			 * computation...
			**/
	step = (ro>co? ro: co);
	rel = MMALLOC((step+2)*sizeof(rel[0]));
	if (step==ro)  re = ematrix_refer(eg,-1,-1,0,co);
	else  re = ematrix_refer(eg,0,ro,-1,-1);
	for (i=0; i<step; i++) {
		if (step==ro)  ematrix_refadd_row(re,i);
		else  ematrix_refadd_col(re,i);
		rel[i] = ematrix_copy(re);
	}
	dispose_ematrix(re);
	epl = MMALLOC((step+2)*sizeof(epl[0]));
	for (i=0; i<step+1; i++)  epl[i] = NULL;
	
			/**
			 * Here is the main backtracking code of this function.
			 * At each level k (0..step-1), we keep in epl[k] the list of
			 * partial matrices generated at this level (matching rel[k]);
			 * and we collect the full representations found in elo.
			 * At each step, we call grepr_addline() to add one more line
			 * (row) to the previous partial matrices.
			**/
	k = br = 0;
	while (k>=0) {
		if (k>=step)  {PROGERROR("cannot get here!");}
		if (alist_getlength(epl[k])>0)  {PROGERROR("cannot get nonempty list here!");}
		DEBUG(CURDLEV+2,"%*s calling GREPR at k=%d, left %d\n",2*k,"",k,(k>0?alist_getlength_part(epl[k-1]):-1));
		if (k>0) {
			ee = *epl[k-1];
			if (!ee || br) {
				dispose_alist_mats(epl[k-1]);  epl[k-1] = NULL;
				if (--k<=0)  break;  else  continue;
			} else  epl[k-1]++;
		} else  ee = NULL;
		p = grepr_addline(ee,rel[k],-1,xpfd,epl+k);
		if (p<=0)  continue;
			/**
			 * When some partial representations are generated in grepr_addline(),
			 * we store them for use in the next step; or in the result list elo.
			 * We are finished with the first representation(s) if ch<=1.
			**/
		if (k<step-1) { ++k;  continue; }
		ret += p;
		if (ch<=1)  br = 1;
		DEBUG(CURDLEV,"Found %d (%d tot) representations of %p[%.12s] (%dx%d) over %s here.\n",
				alist_getlength(epl[k]),ret,eg,EMNAME(eg)?EMNAME(eg):"",ro,co,pfield_curname());
		EMATDEBUGS(CURDLEV+2,*epl[k],"\t\t.=.\t");
		elo = alist_applist(elo,epl[k]);
		epl[k]= NULL;
	}
	
#ifndef FASTPROG	
	if (xp!=pfield_curindex())  {PROGERROR("The current pfield has changed!");}
	if (ret>0 && alist_getlength(elo)<=0)  {PROGERROR("Where is the generated representation?!");}
	DEBUG(CURDLEV-1-(ret>0),"Generated total %d representations of a matrix %p[%.12s] (%dx%d) over %s...\n",
				ret,eg,EMNAME(eg)?EMNAME(eg):"",ro,co,pfield_curname());
	if (IFRANDDEBUGLESS(222) && ch>=0) {
		if (xpfd>=0) pfield_switchto_fast(xpfd);
		ee = ematrix_copydual(eg);
		for (i=0; i<4 && i<ro && i<co; i++)
			if (SIGNM(ee,i,i)!=0)  ematrix_pivot(ee,i,i);
		if (xpfd>=0 && xp>=0) pfield_switchto_fast(xp);
		p = grepr_generate_ext(-1,ee,xpfd,NULL);
		if ((p>0)!=(ret>0))  {PROGERROR("Wrong result of a recursive call! %d!=%d",p,ret);}
	}
#endif
	if (ch==1 || ch>=3)  OUTPUT("There %s %s-representation of the matroid [%.18s] (%dx%d).\n",
				(ret>0?"+IS+ a":"is -NO-"),pfield_curname(),EMNAME(eg)?EMNAME(eg):"",ro,co);
	if (ret>0 && ch>=3)  OUTPUT(" Generated total %d representations of the matroid [%.18s] (%dx%d) over %s.\n",
				ret,EMNAME(eg)?EMNAME(eg):"",ro,co,pfield_curname());
	if (ret>0 && ch>=4)  for (x=elo; x?*x:0; x++) {
		snprintf(buf,80,"%sr%d\t",printoutpref,(int)(x-elo)+1);
		EMATOUTPUTS(*x,buf);
	}
	if (rel)  FREE(rel);  if (epl)  FREE(epl);
	if (lout)  *lout = (*lout? alist_applist(*lout,elo): elo);
	else if (elo)  dispose_alist_mats(elo);
	return ret;
}



/**
 * This function tests representability of the given matroid eg (viewed over pfield xpfd)
 * in the current pfield.
 * The result is >=1 for representable, and 0 for not.
 * 	.......
**/

int	grepr_isrepresented(ematrix *eg, int xpfd) {
	
	if (xpfd==pfield_curindex())  return 5;
	
	//******** decide "included" representability, like regular, binary, etc...
	
	//******** consider also excluded minors? like U24...
	
	//******** make more clever..... (basis with the most zeros, special cases) !!
	
	return (grepr_generate_ext(0,eg,xpfd,NULL)>0);
}























