
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
 * Some theory for start... about matroid connectivity.
 * (Read general description of "structural" functions on the top of ../include/struct.h.)
 * 
 * See the function ematrix_whatsep() in ../ematrix/ematrixop.c for computing the
 * connectivity function of a separation.
 * We use the formula for the connectivity of (A,B) as lambda(A) = r(A)+r(B)-r(M).
 * So the value of this function is 0 for 1-separation, 1 for exact 2-separation, etc...
 * A matroid is called k-connected (k>=2) if it has no separation of lambda <k-1.
 * (A connected matroid has a connected matrix, and it corresponds to a simple
 * vertex-2-connected graph. We mostly consider 3-connected matroids.)
 * 
 * 
 * 
**/






#include "macek.h"
#include "str.h"










/******************	Matroid connectivity	*******************/
/******************************************************************/


#undef CURDLEV
#define CURDLEV		7		/* (need <=DEBUGLEV+2 to perform paranoic tests) */


/**
 * This function computes a "connected order" of columns/rows (tr==1/0) of the matrix e.
 * Connected order means the order in which the columns of e should be taken (starting with 0)
 * so that the partial matrix stays connected at each step (up to possible zero cross-lines!).
 * If the matrix is disconnected, a complete order is still computed so that it is
 * stepwise connected on each component.
 * The return value is 1 for a connected order, and -1 for disconnected.
 * 
 * The order of cols is stored to the given array cor[], or nowhere if cor==NULL is given.
 * (cor[] must be given allocated to proper size!)
 * If cp is given, then the indices of components are stored to cp[] for the columns (1..).
 * (These are indexed the same as cor[], or by absolute matrix indices if cor=NULL.)
 * If cxp is given, then the same component indices are recorded also for (perpendicular)
 * rows/columns of the matrix, and zero lines get index 0.
 * The function may also compute an adjacency matrix cadj[][] for the columns/rows of e,
 * where -1 means no adjacency, and p>=0 means adjacency in a row/column p.
 * The adjacency matrix is indexed by pairs of columns/rows of e (not by cor[]).
**/

int	struct_connorder_ext(ematrix *e, int tr, int *cor, int *cp, int *cxp, int **cadj) {
	int	i,ii,j,k, r,c, cornul=0, *ud,udstack[30];
	
	if (!tr)  ematrix_transpose(e);
	if (COLSM(e)<14)  ud = udstack;
	else  ud = MMALLOC((2*COLSM(e)+2)*sizeof(ud[0]));
	if (cxp && cor)  {PROGERROR("not allowed to use cxp[] with cor[] !!"); cxp=NULL;}
	if (cor==NULL) { cornul = 1;  cor = ud+COLSM(e); }
	for (j=0; j<COLSM(e); j++)  ud[j] = cor[j] = 0;
	if (cxp)  for (j=0; j<ROWSM(e); j++)  cxp[j] = 0;
	if (cadj)  for (i=0; i<COLSM(e); i++)
		for (ii=0; ii<COLSM(e); ii++)  cadj[i][ii] = -1;
		/**
		 * We simply search the columns of e, for each one looking at the rows it
		 * has nonzero, and adding other columns which have nonzeros in this row.
		 * ud[ii] marks columns that were already reached, starting from ud[0].
		 * If k>j happens, that means we ran out of connected columns to continue
		 * the search, and we start with the first next component.
		 * (This happens always in the first run.)
		**/
	j = -1;  c = 0;  r = 1;
	for (k=0; k<COLSM(e); k++) {
		if (k>j) {
			for (ii=0; ii<COLSM(e) && ud[ii]; ii++) ;
			if (ii>=COLSM(e) || k>j+1)  {PROGERROR("must find another component here!"); break;}
			cor[++j] = ii;  ud[ii] = ++c;
			if (k>0)  r = 0;		/* next component reached in the graph -> disconnected */
		}
		for (i=0; i<ROWSM(e); i++)
		  if (SIGNM(e,i,cor[k])!=0)
			for (ii=0; ii<COLSM(e); ii++)
			  if (SIGNM(e,i,ii)!=0 && ii!=cor[k]) {
				if (!ud[ii]) {		/* next column reached in the graph */
					cor[++j] = ii;  ud[ii] = c;
				}
				if (cadj)		/* record adjacency between the columns */
					cadj[cor[k]][ii] = cadj[ii][cor[k]] = i;
			}
		if (j+1>=COLSM(e) && !cadj)  break;
	}
		/**
		 * Additional data: Store component indices for the columns and the rows
		 * of (transposed?) matrix e.
		 * These are indexed in the same way as cor[] if given (only for cp[]),
		 * or by absolute matrix line indices otherwise.
		**/
	if (cp)  for (j=0; j<COLSM(e); j++) {
		if (!cornul)  cp[j] = ud[cor[j]];
		else  cp[j] = ud[j];
	}
	if (cxp)  for (i=0; i<ROWSM(e); i++) {
		for (j=0; j<COLSM(e) && !SIGNM(e,i,j); j++) ;
		if (j<COLSM(e))  cxp[i] = ud[j];
	}
#ifndef FASTPROG
	if (DEBUGLEV>=CURDLEV+2) {
		DEBUG(CURDLEV+2,"Found%s connected col order: ",r?"":" -NO-");
		for (ii=0; ii<COLSM(e); ii++)  SDEBUG(CURDLEV+2,"%d, ",cor[ii]);
		SDEBUG(CURDLEV+2,"\n");  EMATDEBUG(CURDLEV+4,e,"\t\t\t<\t");
	}
	if (cadj && DEBUGLEV>=CURDLEV+2) {
		EMATDEBUGS(CURDLEV+2,e,"\t\tconn\t");
		DEBUG(CURDLEV+2," - cadj[][] adjacencies (-1 for no adj):\n");
		for (i=0; i<ROWSM(e); i++)  for (j=0; j<COLSM(e); j++)
			SDEBUG(CURDLEV+2,"%s% d%c",j==0?"\t\t":"",cadj[i][j],j==COLSM(e)-1?'\n':'\t');
	}
	if (j+1>=COLSM(e) && IFRANDDEBUG(33))  for (ii=0; ii<COLSM(e); ii++)
		if (!ud[ii])  {PROGERROR("wrong computation of a connected order!");}
	if (j+1<COLSM(e) && IFRANDDEBUGLESS(333)) {
		for (i=0; i<COLSM(e); i++)  for (ii=i+1; ii<COLSM(e); ii++)
			if (ud[i]!=ud[ii])  for (j=0; j<ROWSM(e); j++)
				if (SIGNM(e,j,i)!=0 && SIGNM(e,j,ii)!=0)  {PROGERROR("wrong computation of a (dis-)connected order!");}
	}
#endif
	if (ud && ud!=udstack)  FREE(ud);
	if (!tr)  ematrix_transpose(e);
	return (r?1:-1);
}



/**
 * This function tests/determines connectivity of the given matroid.
 * (See above and ematrix_whatsep() for the definition of matroid connectivity.)
 * If cn==-1, then the function returns the actual connectivity of the matroid (>=1).
 * If cn>0, then the function tests whether the matroid is (at least) cn-connected,
 * and the return value is 1 for yes and -1 for no.
 * If sepout is given (must be initialized to *sepout==NULL!), then *sepout collects
 * the list of all interesting small separations.
 * If the pointer ixcon is given, then the function moreover tests for internal
 * (cn+1)-connectivity (all cn-separations have one side <=cn), and sets ixco==1 if so.
 * Here U_2,n is not considered internally 4-connected for any n.
**/

int	struct_connectivity_ext(int ch, ematrix *e, int cn, int *ixcon, ematrix ***sepout) {
	ematrix		**mat,**x,*y;
	int		i,j,k, r,s, mr,con;
	
	if (!e)  {PROGERROR("The matrix e must be given here!"); return -1;}
	if (sepout?*sepout:0)  {PROGERROR("The list sepout must be given initialized to NULL!");}
	if (cn==0)  {PROGERROR("Do not ask about connectivity 0 - disconnected matroids have connectivity 1!");}
	if (ixcon)  *ixcon = 1;	/* one may ask for internal (cn+1)-connectivity as an addition */
	DEBUG(CURDLEV+1,"Computing matroid connectivity of %p [%s], cn=%d ...\n",e,EMNAME(e)?EMNAME(e):"",cn);
	mr = (ROWSM(e)<COLSM(e)? ROWSM(e):COLSM(e));
	if (cn>mr+1)  return -1;
	con = mr+1;
		/**
		 * We separately test connectivity and 3-connectivity by faster routines.
		 * A disconnected matroid either has a loop or a coloop, or it has no
		 * "connected order" of columns.
		 * A matroid that is not 3-connected is either not simple, or it has
		 * a nonempty all-zero submatrix with the complement of rank 1.
		 * These separations here are nontrivial automatically.
		**/
	if (con>=mr && (cn<0 || cn>1)) {
		if (ematrix_linezero_any(e)>=0)  con = 1;
		else if (struct_connorder_col(e,NULL)<0)  con = 1;
	}
	if (con>=mr && (cn<0 || cn>2)) {
		if (ematrix_parallel_any(e)>=0)  con = 2;
		else {
			mat = ematrix_submatrices_zero(e);	/* (rank 0 nonempty submatrices) */
			for (x=mat; (x?*x:0) && con>2; x++) {
				y = ematrix_refextract_xall(e,*x);
				if (ematrix_ismatrank_less(y,2))  con = 2;
				dispose_ematrix(y);
			}
			if (mat)  dispose_alist_mats(mat);
		}
	}
		/**
		 * If there is an exact k-separation in the matroid (whatsep==k-1),
		 * then one of the two complementary submatrices has rank at most (k-1)/2,
		 * so we look for it by brute force...
		 * Only non-trivial (>=k both sides) separations are considered.
		 * We also collect the small separations in *sepout if requested.
		 * If ixcon!=NULL is given, then we have to look for one larger separations
		 * for testing also internal (cn+1)-connectivity.
		 * The meaning of variables is the following:
		 *  cn is the asked connectivity (cn==-1 to determine ot).
		 *  con is the best upper bound for connectivity found so far (starts at rank+1).
		**/
	if (sepout && con<3) {
		con = 3;  k = 1;
		if (ixcon)  *ixcon = 0;	/* (no internal connectivity tested for <3) */
	} else  k = 3;
	for (; k<con; k++)  if (con>=mr && (cn<0 || k<cn+(ixcon!=NULL))) {
		mat = ematrix_submatrices_conn(e,(k-1)/2);
		for (x=mat; x?*x:0; x++)
		  if (ROWSM(e)-ROWSM(*x)+COLSM(*x)>=k && COLSM(e)+ROWSM(*x)-COLSM(*x)>=k) {
			y = ematrix_refextract_xrow(e,*x);
			s = ematrix_whatsep_bound(e,y,k);
			if (s<k)  con = k;	/* (we have found a smaller separation) */
			if (s<k-1)  {PROGERROR("How is it that we have found a smaller separation then before? %d<%d",s,k-1);}
			if (sepout && s<k)  *sepout = alist_append(*sepout,y);
			else  dispose_ematrix(y);
			if (ixcon && s<k) if (ROWSM(e)-ROWSM(*x)+COLSM(*x)>k && COLSM(e)+ROWSM(*x)-COLSM(*x)>k)
				*ixcon = 0;
			if (s<k && !sepout && !ixcon)  break;
		}
		if (mat)  dispose_alist_mats(mat);
	}
	r = (cn<0? con: (con>=cn?1:-1));
	
#ifndef FASTPROG		/* paranoic testing of the connectivity computation: */
	if (IFRANDDEBUGLESS(222) && ch>=0) {
		y = ematrix_copydual(e);
		for (k=0; k<4; k++) {	/* (extra testing with pivoted matrix) */
			i = RANDOM()%ROWSM(y);  j = RANDOM()%COLSM(y);
			if (SIGNM(y,i,j)!=0)  ematrix_pivot(y,i,j);
		}
		if (r!=struct_connectivity_ext(-1,y,cn,NULL,NULL))  {PROGERROR("incorrect computation of matroid connectivity, ret=%d",r);}
		dispose_ematrix(y);
	}
	DEBUG(CURDLEV+1," - connectivity computed %s %d\n",(cn<0?"exactly":(r>=0?"at least":"less")),con);
#endif
	ch = i = j = k;
	return r;	/* the return value  (cn<0? con: (con>=cn?1:-1))  computed above */
}



/**
 * This function computes induced paths between vertices of the graph in cadj[][].
 * We define a graph on the indices of cadj where edges connect pairs
 * that have cadj[i][j]>=0.
 * For each pair of vertices in this graph we find an induced path joining them,
 * and using only other vertices before (=not later than the second one in order by cor[])
 * if the vertex order cor[0...sz-1] is given.
 * We record the internal vertices of these paths in the bitfields ppat[][],
 * and the cadj-values needed for the paths in the bitfields mpat[][].
 * (These paths need not be the shortest ones, just induced.)
 * The paths are all indexed in order of cor[] (not as cadj[][]!) if it is given.
**/

int	struct_connpaths_ext(int sz, int **cadj, int *cor, long **ppat, long **mpat) {
	int	i,j,k, ii,jj,kk,ll;
	
	if (!cadj || !ppat)  {PROGERROR("both cadj, ppat must be given here"); return -1;}
	for (i=k=0; i<sz; i++)  for (j=0; j<sz; j++) {
		if (cadj[i][j]>k)  k = cadj[i][j];
		if (cadj[i][j]!=cadj[j][i])  {PROGERROR("the array cadj[][] must be symmetric");}
	}
	if (sz>=NUM_LONG_BITS-2 || (k>=NUM_LONG_BITS-2 && mpat))
		{PROGERROR("connecting paths do not fit into long bits %d %d",sz,k); return -1;}
	
		/**
		 * We use classical "matrix" algorithm for computing a graph metric, adapted to
		 * restricted induced paths.
		 * At round kk, we find all paths using only internal vertices from
		 * cor[0],...,cor[kk]==k (indices with respect to cadj[][]).
		**/
	for (kk=-1; kk<sz; kk++) {
		k = (cor&&kk>=0? cor[kk]:kk);
		ll = 1;
		for (ii=0; ii<sz; ii++)
		  for (jj=ii+1; jj<sz; jj++) {
		  	i = (cor? cor[ii]:ii);
		  	j = (cor? cor[jj]:jj);
			if (kk==-1) {	/* the first run records edges between indices */
				if (cadj[i][j]>=0) {
					ppat[ii][jj] = 0;
					if (mpat)  mpat[ii][jj] = (1l<<cadj[i][j]);
				} else  ppat[ii][jj] = -1;
					/* the next runs look for better paths between them */
			} else if (!cor || kk<=ii || kk<=jj) {
			  		/* (only rows before ii,jj in cor[] may be used inside paths) */
				if (ppat[ii][jj]<0)
				  if (ppat[ii][kk]>=0 && ppat[kk][jj]>=0) {
					ppat[ii][jj] = (ppat[ii][kk] | ppat[kk][jj] | (1l<<kk));
					if (mpat)  mpat[ii][jj] = (mpat[ii][kk] | mpat[kk][jj]);
				}
			}
			ppat[jj][ii] = ppat[ii][jj];
			if (mpat)  mpat[jj][ii] = mpat[ii][jj];
			if (ppat[ii][jj]<0)  ll = 0;
		}
		if (ll==1)  break;	/* when all paths (connections) are already found */
	}
#ifndef FASTPROG
	//
#endif
	return 1;
}


//********* function to decompose a matroid into 3-connected components......???



/**
 * This function tests (co)simplicity of the given matroid (tr=1 simple / tr=0 cosimple).
 * Actually, the function looks for zero, parallel, or unit columns / rows in the matrix.
 * The return value is 1 for (co)simple, and -lo-1 if line lo is parallel to something else.
**/

int	struct_simplicity_ext(int ch, ematrix *e, int tr) {
	int	j, ret;

	if (!e)  {PROGERROR("The matrix e must be given here!"); return -1;}
	DEBUG(CURDLEV+1,"Computing matroid simplicity of %p [%s], tr=%d ...\n",e,EMNAME(e)?EMNAME(e):"",tr);
	if (!tr)  ematrix_transpose(e);
	junk = ch;  ret = 1;
	for (j=0; j<COLSM(e) && ret>0; j++)
		if (ematrix_parallel_col(e,j)>=0 || ematrix_lineunit(e,1,j)>=0)
			ret = -j-1;
	if (!tr)  ematrix_transpose(e);
	DEBUG(CURDLEV+1," - %ssimplicity computed %d\n",tr?"":"co",ret);
	return ret;
}



















/******************	Finding short cycles	*******************/
/******************************************************************/


#undef CURDLEV
#define CURDLEV		7		/* (need <=DEBUGLEV+2 to perform paranoic tests) */


/**
 * This function tests/computes the shortest cycle of the given matroid.
 * If cc==-1, then the function returns the length of the shortest cycle (>=1).
 * (If the matroid has no cycle, then the return value is 9999.)
 * If cc>=0, then the function returns number cc if the length of the shortest cycle is >=cc,
 * but it returns -1 if there is a cycle shorter than cc.
 * If depout is given (must be initialized to *depout==NULL!), then *depout collects
 * the list of all shortest cycles in e (each cycle by a submatrix on the cycle elements).
**/

int	struct_matgirth_ext(int ch, ematrix *e, int cc, ematrix ***depout) {
	ematrix		**mat,**x,*y;
	int		i,j,k,r;
	
	if (!e)  {PROGERROR("The matrix e must be given here!"); return -1;}
	if (depout?*depout:0)  {PROGERROR("The list depout must be given initialized to NULL!");}
	DEBUG(CURDLEV+1,"Computing the shortest cycle (girth) of %p(%dx%d) [%s], cc=%d ...\n",e,ROWSM(e),COLSM(e),EMNAME(e)?EMNAME(e):"",cc);
	if (ROWSM(e)<=0)  return 0;
	if (COLSM(e)<=0)  return 9999;
	if (cc>ROWSM(e)+1)  return -1;
		/**
		 * We look for the shortest cycle(s) by brute force - trying all
		 * k-element subsets of e and computing their rank in the matroid.
		 * Whenever a cycle is found, we stop further search -- no bigger k.
		 * We also collect the shortest cycles in *depout if requested.
		**/
	r = j = -1;
	for (k=1; (k<=cc || cc<0); k++) {
		if (k>ROWSM(e)+1)  {PROGERROR("something wrong - no cycle found?! k=%d, cc=%d, rows=%d",k,cc,ROWSM(e)); break;}
		mat = ematrix_submatrices_sub(e,k);
		DEBUG(CURDLEV+2,"Trying list of %d subsets of size %d...\n",alist_getlength(mat),k);
		for (x=mat; x?*x:0; x++) {
			i = ematrix_setrank(e,*x);
			if (i>=k)  continue;
			if (i<k-1)  {PROGERROR("something wrong - missed a shorter cycle?! k=%d,i=%d",k,i);}
			if (r<=0)  r = k;
			if (depout && r==k)  *depout = alist_append(*depout,ematrix_copy(*x));
			else  break;
		}
		if (mat)  dispose_alist_mats(mat);
		if (r>0)  break;
	}
	if (cc>0 && r<0)  r = cc;	/* (if no cycle shorter than cc found above) */
	
#ifndef FASTPROG		/* paranoic testing of the cycle computation: */
	if (IFRANDDEBUGLESS(222) && ch>=0) {
		y = ematrix_copy(e);
		for (k=0; k<6; k++) {	/* (extra testing with pivoted matrix) */
			i = RANDOM()%ROWSM(y);  j = RANDOM()%COLSM(y);
			if (SIGNM(y,i,j)!=0)  ematrix_pivot(y,i,j);
		}
		if (r!=struct_matgirth_ext(-1,y,cc,NULL))  {PROGERROR("incorrect computation of shortest cycle, ret=%d",r);}
		dispose_ematrix(y);
	}
	DEBUG(CURDLEV+1," - shortest cycle (girth) computed %s %d\n",(cc<0?"exactly":(r>=0?"":"greater than")),r);
#endif
	return r;	/* the return value  k, cc, or -1  computed above */
}

















/******************	Matroid branch-width	*******************/
/******************************************************************/


#undef CURDLEV
#define CURDLEV		6		/* (need <=DEBUGLEV+2 to perform paranoic tests) */



/**
 * We provide a simple greedy (but precise) algorithm to test for branch-width 3 of a matroid.
 * We use the fact that any 3-separating set of size >=4 can be displayed in some
 * bw-3 decomposition, even when considering an already partitioned matroid.
 * On the other hand, if a bw-3 decomposition exists, then there always is 3-separating
 * subset of <=6 elements (using two pairs or triangles or triads),
 * or <=4 blocks if using bigger than 1-element part(s).
 * Hence we just look for such possible 3-separating set.
 * If we find one, we make its union as a new block in the matroid.
 * If we succeed all the way to just two blocks at the end, we have a bw-3 decomposition.
 * Otherwise, the matroid has no such decomposition.
 * 
 * The return value is 1 for branch-width <=3 (printed via *decomp), and -1 if not.
 * (ldecomp is not yet implemented.)
 * If decomp is given, then a copy of a string describing the resulting decomposition
 * is stored into *decomp; user has to free(!) this string after use.
 * If ch>1, then the decomposition is printed here.
 * The function requires the matroid to be 3-connected !!!
**/

#define	GETREFMLINE(e,r,i)	((i)<ROWSM(r)? GETREFMROW(r,i):GETREFMCOL(r,(i)-ROWSM(r))+ROWSM(e))

int     struct_hasbwidth3_ext(int ch, ematrix *em, char **decomp, int *ldecomp) {
	ematrix		*ee, *ex,*e1,*e2, *eu, *ef,*eff;
	ematrix		**list, **lpair, **lln,*llnst[30], **x,**y;
	int		i,j,k, a,b, r, msz;
	char		buf[400];
	
	if (ldecomp!=NULL)  {PROGERROR("return ldecomp[] is not implemented yet!"); ldecomp = NULL;}
	if (!em)  {PROGERROR("the matrix must be given here"); return -1;}
	ee = ematrix_copy(em);
	msz = ROWSM(ee)+COLSM(ee);
#ifndef	FASTPROG
	if (msz<10 || (msz<14 && IFRANDDEBUGLESS(22)) || IFRANDDEBUGLESS(222)) {
		if (!struct_isconnected(ee,3))  {PROGERROR("Can correctly test branch-width 3 only for 3-connected matrices!");}
	}
	DEBUG(CURDLEV,"Testing branch-width 3 for the matroid %p [%s]...\n",ee,EMNAME(em)?EMNAME(em):"");
#endif
	lln = (msz<28?llnst: MMALLOC((msz+2)*sizeof(lln[0])));
	list = ematrix_submatrices_sub(ee,1);
	if (alist_getlength(list)!=msz)  {PROGERROR("wrong list of singletons generated here ?!?");}
	for (x=list; x?*x:0; x++) {
		snprintf(buf,5," %d ",(ROWSM(*x)>0?ROWSID(*x,0):COLSID(*x,0)));
		EMSETNAME(*x,buf);
		lln[GETREFMLINE(ee,*x,0)] = *x;
	}
			/**
			 * We start with an initial partition  list  of the matroid into singletons.
			 * Then we cycle the next procedure until the list has one part only.
			 * Trivial cases are sorted first - a matroid on <=6 elements always
			 * has bw <=3, and if only <=3 singletons remain besides one larger
			 * part, then we are done as well.
			 * The supplementary array refers to those parts in  list  that are
			 * still singletons (by their row/column indices in ee).
			**/
	while (alist_getlength(list)>1) {
		DEBUG(CURDLEV+2,"  new cycle for list %d\n",alist_getlength(list));
		lpair = NULL;  eu = e1 = e2 = NULL;
		r = 0;
		for (x=list, a=b=0,y=NULL; *x && a<2 && b<7; x++,b++)
			 if (ROWSM(*x)+COLSM(*x)>1) { a++;  y = x; }
		if (a==0 && b<=6) {
			e1 = ematrix_refer_all(ee);  eu = ematrix_refer_all(ee);  r = 1;
		} else if (a==1 && b<=4) {
			e1 = *y;  e2 = ematrix_refextract_xall(ee,*y);
			eu = ematrix_refer_all(ee);  r = 2;
		}
			/**
			 * Here we try all pairs of parts in  list.
			 * If their union is 3-separating, and not both are singletons,
			 * then we have a new part for replacement.
			**/
		for (x=list; *x && !r; x++) for (y=x+1; *y && !r; y++) {
		  if ((a=ROWSM(*x)+COLSM(*x))>1 || (b=ROWSM(*y)+COLSM(*y))>1) {
			eu = ematrix_union(*x,*y);
			if (ematrix_whatsep(ee,eu)<3) {
				if (a>=b) { e1 = *x;  e2 = *y; }
				else { e1 = *y;  e2 = *x; }
				r = 3;
			} else {
				dispose_ematrix(eu);
			}
			/**
			 * If both in the pair are singletons, then we save this pair for
			 * later testing, and we try our luck with the line-(co)closure.
			 * If the closure has >=4 elements among singletons, then it forms
			 * a new part.
			 * We save a possible closure triple for later testing, but only
			 * if the third element is in the middle of the pair.
			**/
		  } else {
			ef = ematrix_union(*x,*y);
			lpair = alist_append(lpair,ef);
			for (k=0; k<2 && !r; k++) {
				eff = ematrix_closure_tr(ee,ef,k);
				a = 1;
				for (i=ROWSM(eff)+COLSM(eff)-1; i>=0; i--) {
					j = GETREFMLINE(ee,eff,i);
					if (lln[j]==NULL) {
						if (i<ROWSM(eff))  ematrix_remove_row(eff,i);
						else  ematrix_remove_col(eff,i-ROWSM(eff));
					} else {
						a = a && (j==GETREFMLINE(ee,*x,0) || j==GETREFMLINE(ee,*y,0)
							|| (j>GETREFMLINE(ee,*x,0) && j>GETREFMLINE(ee,*y,0)));
					}
				}
				if (ROWSM(eff)+COLSM(eff)>3) {
					e1 = eff;  eu = ematrix_refer_all(eff);  r = 1;
				} else if (a && ROWSM(eff)+COLSM(eff)==3) {
					lpair = alist_append(lpair,eff);
				} else  dispose_ematrix(eff);
			}
		  }
		  DEBUG(CURDLEV+3,"   cycle x=%d, y=%d, r=%d, pair %d\n",(int)(x-list),(int)(y-list),r,alist_getlength(lpair));
		}
			/**
			 * Then we have to try all pairs that take one non-singleton part
			 * and the other as one of the above pairs or triples,
			 * or two disjoint above pairs or triples.
			 * Again, we look at their union whether it is 3-separating.
			**/
		for (x=list; *x && !r; x++) if (ROWSM(*x)+COLSM(*x)>1) {
			for (y=lpair; (y?*y:0) && !r; y++) {
				eu = ematrix_union(*x,*y);
				if (ematrix_whatsep(ee,eu)<3) {
					e1 = *x;  e2 = ematrix_refer_all(*y);  r = 2;
				} else  dispose_ematrix(eu);
			}
		}
		for (x=lpair; (x?*x:0) && !r; x++) {
			if (ROWSM(*x)+COLSM(*x)>3)  {PROGERROR("only pairs and triples should be here");}
			for (y=x+1; *y && !r; y++)
			  if (ematrix_isdisjoint(*x,*y)) {
				eu = ematrix_union(*x,*y);
				if (ematrix_whatsep(ee,eu)<3) {
					e1 = ematrix_refer_all(*x);  e2 = ematrix_refer_all(*y);
					r = 1;
				} else  dispose_ematrix(eu);
			}
		}
		if (lpair)  dispose_alist_mats(lpair);
			/**
			 * Finally, if we have found a union (in eu) for the new part (r>0),
			 * then we have to remove the current parts of the union,
			 * print this new branch for records, and store the union into the list.
			**/
		if (r<=0)  break;
		buf[0] = 0;
		for (k=0, ex=e1; ex && k<2; k++, ex=e2) {
		  if (r>k+1) {
			snprintf(buf+strlen(buf),180,"(%s)",EMNAME(ex));
			list = alist_delete_val(list,ex);
		  } else {
			for (i=ROWSM(ex)+COLSM(ex)-1, a=0; i>=0; i--) {
				j = GETREFMLINE(ee,ex,i);
				ef = lln[j];  lln[j] = NULL;
				if (!ef)  {PROGERROR("should find only singleton elements here"); continue;}
				if (!a++)  snprintf(buf+strlen(buf),3,"(");
				if (strlen(buf)<380)  snprintf(buf+strlen(buf),5,"%s",EMNAME(ef));
				list = alist_delete_val(list,ef);
				dispose_ematrix(ef);
			}
			if (a>0)  snprintf(buf+strlen(buf),3,")");
		  }
		  dispose_ematrix(ex);
		}
		DEBUG(CURDLEV+1," - new branch found   %s\n",buf);
		EMSETNAME(eu,buf);
		list = alist_append(list,eu);
	}
			/**
			 * Here we prepare the return values and possible debug printing.
			 * We also try the result recursively on an equivalent dual matrix.
			**/
	if (alist_getlength(list)<=1) {
		r = 1;
		if (decomp && list[0])  *decomp = MSTRDUP(EMNAME(list[0]));
		if (ch>1)  OUTPUT("A width-3 branch decomposition of [%s] is:  %s\n",EMNAME(em)?EMNAME(em):"",EMNAME(list[0]));
	} else {
		r = -1;
		if (ch>1)  OUTPUT("There is NO width-3 branch decomposition of [%s].\n",EMNAME(em)?EMNAME(em):"");
	}
#ifndef	FASTPROG
	if (r>=0 && ch>=0) {
		DEBUG(CURDLEV-1,"Found a width 3 branch-decomposition of %p [%s] here\n\t\t\t\t\t%s.\n",
					ee,EMNAME(em)?EMNAME(em):"",list[0]?EMNAME(list[0]):"");
		if (msz>1 && (list[0]?ROWSM(list[0])+COLSM(list[0])!=msz:1))  {PROGERROR("wrong final partition in computation of branch-width 3");}
	} else if (ch>=0) {
		DEBUG(CURDLEV-1,"NO width 3 branch-decomposition of %p [%s] exists.\n",ee,EMNAME(em)?EMNAME(em):"");
		for (x=list,a=0; x?*x:0; x++)  a += ROWSM(*x)+COLSM(*x);
		if (a!=msz)  {PROGERROR("lost or extra elements in the partition list! %d!=%d",msz,a);}
	}
	if (IFRANDDEBUGLESS(222) && ch>=0) {
		ematrix_transpose(ee);	/* extra debug testing with pivoted dual matrix */
		for (k=0; k<4; k++) {
			i = RANDOM()%ROWSM(ee);  j = RANDOM()%COLSM(ee);
			if (SIGNM(ee,i,j)!=0)  ematrix_pivot(ee,i,j);
		}
		k = struct_hasbwidth3_ext(-1,ee,NULL,NULL);
		if (k!=r)  {PROGERROR("wrong result when recursivly calling r%d o%d !",k,r);
				EMATDEBUG(0,ee," !r!\t"); EMATDEBUG(0,em," !o!\t");}
	}	/* (ee is modified here!!!) */
#endif
	dispose_alist_mats(list);
	dispose_ematrix(ee);
	if (lln && lln!=llnst)  FREE(lln);
	return r;
}





/**
 * 
 * 
**/













/******************	Checking for fans	*******************/
/******************************************************************/


#undef CURDLEV
#define CURDLEV		7


/**
 * This function looks whether the matrix e has a "fan" of given length, moreover
 * ending with the line lo if tr>=0 (row for tr==0, column for tr==1),
 * and using only the lines of the refering submtrix ec if given.
 * A fan of length d is a sequence of d element where consecutive triples switch between
 * being triangles and triads (parallel or serial pairs are not counted here).
 * 
 * If one asks about fans ending with the given line (tr>=0), then the matroid e-lo
 * without that line must be 3-connected (!otherwise possibly incorrect answer!).
 * This assumption of 3-connectivity implies that the last added line must be the end of
 * a possible fan, and that the end triple must be a triangle (tr==1) or triad (tr==0).
 * 
 * If tr>=0, then only fans ending with the given line lo (row for tr==0, column for tr==1)
 * are searched.
 * If lo==-1, then the last line is considered (determined row/col by tr).
 * If tr==-1, then the the whole matrix is searched without respect to the (last) line.
 * If fan==-1, then the length of the longest fan is returned (so all of them are searched).
 * If ec!=NULL, then only the lines of the refering matrix ec are considered in the search.
 * If ch<-3, then the triple-row/triple-column tests are skipped, as they require proper
 * representability (which may not be correct yet in the generating process).
 * Hence ch<-3 gives only a partial answer, a fan may be missed (this really happens in generating!).
 * 
 * The function returns 0 if there is NO such fan, and 1 if there is some;
 * or the length of the longest fan if fan==-1 was given.
**/

int	struct_hasfan_ext(int ch, ematrix *e, ematrix *ec, int tr, int lo, int fan) {
	ematrix		*ee, *eec;
	char		buf[1000];
	int		i,k,kk, st,tt1,tt2,tt,ttm=0, br, ret, lf,
			*fn=NULL,*fnm, fnstack[104];
	
	if (ec)  if (!ISREFMAT(ec) || REFMAT(ec)!=e || ISTRANSPM(e)!=ISTRANSPM(ec))
		{ PROGERROR("When using ec, it must refer to e in the same transpose state"); ec = NULL; }
	if (ec && tr>=0) { PROGERROR("When using ec, do not use tr>=0"); tr = -1; }
	if (fan>=0 && fan<3) { PROGERROR("What is fan of length %d <3 ??",fan); fan = 3; }
	if (lo<0)  lo = tr? COLSM(e)-1:ROWSM(e)-1;
	buf[0] = 0;  kk = 0;
#ifndef FASTPROG
	DEBUG(CURDLEV+(ch<0),"Looking for %d-fans in %p [%s] %dx%d (tr=%d, lo=%d%s)...\n",
			fan,e,EMNAME(e)?EMNAME(e):"",ROWSM(e),COLSM(e),tr,lo,ec?", ec":"");
	EMATDEBUGS(CURDLEV+1,e,"\t\tf\t");  if (ec) EMATDEBUGS(CURDLEV+1+(ch<0),ec,"\t\t\t-\t");
	if (ch>=0 && tr>=0 && IFRANDDEBUGLESS(555)) {	/* (debug-check for 3-connectivity without lo) */
		ee = ematrix_copy(e);  ematrix_remove_rc(ee,tr,lo);
		if (!struct_isconnected(ee,3))  {PROGERROR("wrong - not 3-connected without lo=%d",lo); EMATDEBUG(0,ee,"!!\t");}
		dispose_ematrix(ee); }
#endif
	if (fan>ROWSM(e)+COLSM(e))  return 0;
	if (ch>2) {
		OUTPUT("Looking for %d%s in the matroid [%s] %dx%d...\n",
			fan,fan<0?" the longest fan":"-fans",EMNAME(e)?EMNAME(e):"",ROWSM(e),COLSM(e));
		if (tr>=0 || ec)  DEBUG(0,"Do not call fan printing with tr=%d>=0 or ec=%p.\n",tr,ec);
		if (fan<0)  DEBUG(CURDLEV-3,"Do not call all fan printing with fan=%d<0 (max fan).\n",fan);
	}
	if (ROWSM(e)+COLSM(e)<50)  fn = fnstack;
	else  fn = MMALLOC((ROWSM(e)+COLSM(e)+2)*2*sizeof(fn[0]));
	fnm = fn+ROWSM(e)+COLSM(e)+2;
	if (tr<0) {		/* no last line for the fan given - trying everything */
		st = 0;  tt1 = 0;  tt2 = 1;
	} else {		/* the last fan line is given, also determines tt by 3-connectivity */
		st = 1;  tt1 = tt2 = tr? 1:0;
		fn[0] = tr? ROWSM(e)+lo:lo;
		DEBUG(CURDLEV+2,".. hasfan pre-determined %s fn[0] = %d, tt=%d\n",tr?"col":"row",lo,tt1);
	}
		/**
		 * We cycle all choices of k lines and starting tt values (triad 0/triangle 1),
		 * beginning with a choice of fn[st], where st and tt1,tt2 depend on the input.
		 * One choice is a row for <ROWSM and a column +ROWSM otherwise.
		 * We record the longest fan length in lf, and the fan itself in fnm[].
		 * If we look only for a fan longer than the given value, we stop immediately
		 * after we find it, otherwise we search through all fans to find the longest one.
		 * As we generate our choice, we test the consecutive triples for being
		 * triangles/triads (by tt).
		**/
	ret = lf = 0;
	for (tt=tt1; tt<=tt2 && !ret; tt+=(ret?0:1)) {
	  k = st;  fn[k] = -1;
	  while (k>=st && !ret) {
			/* the value of fn[k] is the current choice - row for <ROWSM, column+ROWSM otherwise */
		if (++fn[k]>=ROWSM(e)+COLSM(e)) {
			--k;  continue;		/* (no more choices at this level) */
		}
		for (i=0; i<k; i++)  if (fn[i]==fn[k])  break;
		if (i<k)  continue;		/* (if fn[k] was already chosen previously) */
		if (ec) {
			if (fn[k]<ROWSM(e))  for (i=ROWSM(ec)-1; i>=0 && GETREFMROW(ec,i)!=fn[k]; i--) ;
			else  for (i=COLSM(ec)-1; i>=0 && GETREFMCOL(ec,i)!=fn[k]-ROWSM(e); i--) ;
			if (i<0)  continue;	/* (if fn[k] does not belong to ec if given) */
		}
		DEBUG(CURDLEV+2,"%*s. hasfan(%d) trying choice %s fn[%d] = %d (%d) for %s\n",2*k,"",fan,
				fn[k]<ROWSM(e)?"row":"col",k,fn[k]<ROWSM(e)?fn[k]:fn[k]-ROWSM(e),fn[k]<ROWSM(e)?ROWSID(e,fn[k]):COLSID(e,fn[k]-ROWSM(e)),(k+tt)%2?"triang":"triad");
		
		if (k>=2) {	/* looking for a triangle ((k+tt)%2==1) or triad ((k+tt)%2==0) at this level */
			br = struct_hasfan_triax((ch<-3?-1:0),e,(k+tt)%2,fn[k-2],fn[k-1],fn[k],1);
			if (br<=0)  continue;	/* if there was no triaxx, then cycle other choices */
		}
		if (++k>lf) {			/* records the longest fan found so far */
			lf = k;  ttm = tt;
			for (i=0; i<k; i++)  fnm[i] = fn[i];
			DEBUG(CURDLEV+2,"-%*s. hasfan(%d) found longer = %d\n",2*k,"",fan,lf);
		}
		if (k<fan || fan<0) {		/* to the next choice of line */
			fn[k] = -1;  continue;
		}
		if (k>ROWSM(e)+COLSM(e))  {PROGERROREXIT("Out of index range looking for a fan");}
		fn[k] = 11111;			/* fan of given length found here, plus printing */
		if (ch>=2 || CURDLEV-1+(ch<0)<=printlev) {
			sprintf(buf,"fan: (%s+..) ",ttm?"triangle":"triad");
			for (i=0; i<lf; i++)  snprintf((buf[800]=0,buf+strlen(buf)),50,"%s%d(%d), ",
					fn[i]<ROWSM(e)?"r":"c",fn[i]<ROWSM(e)?fn[i]:fn[i]-ROWSM(e),fn[i]<ROWSM(e)?ROWSID(e,fn[i]):COLSID(e,fn[i]-ROWSM(e)));
		}
		if (ch>2)  OUTPUT("\t[%s]#%d %s\n",EMNAME(e)?EMNAME(e):"",++kk,buf);
		else  ret = 1;		/* (this breaks the cycle, unless all printing is required) */
	}}
		/**
		 * If fan>0, then the fans of this length are recorded above, and they are
		 * all printed out if ch>2 is requested.
		 * On the other hand, fan<0 the longest fan is printed out below, separately.
		 * The return value is the longest fan, or the requested fan length, or 0.
		**/
	if (fan<=0)  ret = (lf>=3? lf:0);
	else  ret = (lf>=fan? fan:0);
	if (fan<0) if (ch>=2 || CURDLEV-1+(ch<0)<=printlev) {
		sprintf(buf,"longest fan: (%s+..) ",ttm?"triangle":"triad");
		for (i=0; i<lf; i++)  snprintf((buf[800]=0,buf+strlen(buf)),50,"%s%d(%d), ",
				fnm[i]<ROWSM(e)?"r":"c",fnm[i]<ROWSM(e)?fnm[i]:fnm[i]-ROWSM(e),fnm[i]<ROWSM(e)?ROWSID(e,fnm[i]):COLSID(e,fnm[i]-ROWSM(e)));
	}
	if (ch==2 || (ch>=2 && fan<0))  OUTPUT("\t[%s] %s\n",EMNAME(e)?EMNAME(e):"",buf);
	ee = eec = NULL;
#ifndef FASTPROG
	DEBUG(CURDLEV-1+(!ret)+(ch<0),"- hasfan(%d) test for %p [%s] (l %d, tr=%d%s) found max %d  %s.\n",
			fan,e,EMNAME(e)?EMNAME(e):"",lo,tr,ec?", ec":"",lf,fan>0?"":(ret?"+HAS+":"NO fan"));
	SDEBUG(CURDLEV-1+(!ret)+(ch<0),"\t\t\t ->  %s\n",buf);
			/* (debug-testing the full-test for diff basis of e) */
	if (ch>=0 && tr<0 && IFRANDDEBUG(222)) {
		ee = ematrix_copy(e);  EMSETNAME(ee,"rec-test");
		for (eec=ec, tt=0; tt<5; tt++) {
			i = random()%ROWSM(ee);  k = random()%COLSM(ee);
			if (SIGNM(ee,i,k)==0)  continue;
			if (ec!=NULL) {	/* (not to interfere with the matrix ec if given) */
				for (tt1=ROWSM(ec)-1; tt1>=0 && GETREFMROW(ec,tt1)!=i; tt1--) ;
				for (tt2=COLSM(ec)-1; tt2>=0 && GETREFMCOL(ec,tt2)!=k; tt2--) ;
				if (tt1>=0 && tt2>=0)  ematrix_pivot(ee,i,k);
			} else  ematrix_pivot(ee,i,k);
		}
		if (ec) { eec = ematrix_refer_all(ec);  ematrix_rerefer(eec,ee); }
		else if (random()%2==1)  ematrix_transpose(ee);
		if (ret!=(kk=struct_hasfan_ext(-1,ee,eec,-1,-1,fan)))  {PROGERROR("The result must be the same after pivoting (transpose)! ret=%d != %d",ret,kk);}
		if (ec) dispose_ematrix(eec);  dispose_ematrix(ee);
	}
#endif
	if (fn && fn!=fnstack)  FREE(fn);
	return ret;
}


/**
 * This function is for use in estruct_hasfan_ext() and estruct_hastail_ext() only...
 * It tests a triangle (trt==1) or triad (trt==0) on the given lines of e, taking lX for a row
 * when <ROWSM(e), and for a column otherwise.
 * If nps==1, then only proper trangles and triads are allowed, no parallel or series pairs.
 * If ch<0, then the triple-row/triple-column tests are skipped, as they require proper
 * representability (which may not be correct yet in the generating process).
 * The return value is 1 for a triaxx, 0 for no one, and -1 for parallel/series pair.
**/

int	struct_hasfan_ppair(ematrix *e, int l1, int l2, int t1, int t2) {
	switch (t1+2*t2) {
	case 0:	return 0;
	case 1:	return ematrix_isparallel_rc(e,l2,l1);
	case 2:	return ematrix_isparallel_rc(e,l1,l2);
	case 3:	return ematrix_isparallel_cc(e,l1,l2);
	}
	return 0;
}

int	struct_hasfan_triax(int ch, ematrix *e, int trt, int l1, int l2, int l3, int nps) {
	int	r,t1,t2,t3, ret;
	
	trt = (trt?1:0);  r = ROWSM(e);
	if (!trt)  ematrix_transpose(e);
	t1 = t2 = t3 = !trt;	/* deciding row/column and the test to perform... */
	if (l1>=r) { t1 = trt;  l1 -= r; }
	if (l2>=r) { t2 = trt;  l2 -= r; }
	if (l3>=r) { t3 = trt;  l3 -= r; }
	
	ret = 0;		/* rejecting parallel pairs in the triple if nps is given */
	if (nps)  if (struct_hasfan_ppair(e,l1,l2,t1,t2) ||
			struct_hasfan_ppair(e,l1,l3,t1,t3) || struct_hasfan_ppair(e,l3,l2,t3,t2))
		ret = -1;
				/* testing for the triangle (triad is now transposed): */
	if (ret>=0)  switch (t1+2*t2+4*t3) {
	case 1:	ret = ematrix_istriangle_rrc(e,l2,l3,l1);
		break;
	case 2:	ret = ematrix_istriangle_rrc(e,l1,l3,l2);
		break;
	case 4:	ret = ematrix_istriangle_rrc(e,l1,l2,l3);
		break;
	case 3:	ret = ematrix_istriangle_rcc(e,l3,l1,l2);
		break;
	case 5:	ret = ematrix_istriangle_rcc(e,l2,l1,l3);
		break;
	case 6:	ret = ematrix_istriangle_rcc(e,l1,l2,l3);
		break;
	case 7:	ret = (ch>=0? ematrix_istriangle_ccc(e,l1,l2,l3):0);
		/* (the matrix may not be representable for ch<0, so cannot call _ccc !) */
		break;
	case 0:	ret = 0;	/* (no triangle on three rows) */
		break;
	default:PROGERROR("Cannot get here!");	ret = -1;  break;
	}
	if (!trt)  ematrix_transpose(e);
	DEBUG(CURDLEV+2," Test returning %s %s.\n",trt?"triangle":"triad",ret>0?"YES":(ret<0?"parallel":"NO"));
	return ret;	/* (the above tests return 1 for triaxx, and 0 or -1 for no one) */
}























































