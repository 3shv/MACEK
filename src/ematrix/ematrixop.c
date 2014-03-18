
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
 * A matrix is a structure holding an array (considered 2-dimensional) of pfield elements,
 * together with its dimensions and maximal dimensions.
 * The basic manipulating functions for matrices (creating, freeing, accessing, refering,
 * printing, etc) are defined in ematrix.c.
 * 
 * This file defines more-involved matrix functions, like computations (pivoting, determinants)
 * and various simple matroid-related tests (rank, dependency, etc).
 * (More complicated structural functions are defined in struct/ * ...)
 * 
 * Look for more information in ../include/ematrix.h.
 * 
**/




#include "macek.h"
#include "emx.h"








/******************	Matrix computations	*******************/
/******************************************************************/


#undef CURDLEV
#define CURDLEV		7

static int	inrecur = 0;	/* (used for some recursive debug tests) */



/**
 * Pivoting operation on the given element ro,so of the matrix e.
 * (Remember that the matrix is supposed to be in a standard form, i.e. without
 * the identity part, and it remains in the form afterwards.)
 * If the computation gets out of in pfield numbers, an error is issued if ch>0.
 * A variant that does not print error messages if ch<=0 is provided as well.
 * Return value 0 means success, and -1 failure of computation.
 * (The modified matrix is stored in the same e.)
**/

int	ematrix_pivot_ext(int ch, ematrix *e, int ro, int co) {
	int	i,j, r;
	exp_t	xx;
	sign_t	gg;
	
	DEBUG(CURDLEV+2,"pivoting on %d,%d in a matrix %p\n",ro,co,e);
	r = 0;
	if (SIGNM(e,ro,co)==0) { r = -1;  PROGERROR("cannot pivot on a 0 element at %d,%d",ro,co); }
	
			/* classical matrix pivoting -- replace the column by a unit vector, but store the other vector instead */
	for (i=0; i<ROWSM(e) && r>=0; i++)
		if (i!=ro) {
			r = pfield_mul_ch(ch,1,EXPM(e,i,co),SIGNM(e,i,co),
				-1,EXPM(e,ro,co),-SIGNM(e,ro,co),&(EXPM(e,i,co)),&(SIGNM(e,i,co)));
		}
	for (i=0; i<ROWSM(e) && r>=0; i++)  for (j=0; j<COLSM(e) && r>=0; j++)
		if (i!=ro && j!=co) {
			r = pfield_mul_ch(ch,1,EXPM(e,i,co),SIGNM(e,i,co),
					1,EXPM(e,ro,j),SIGNM(e,ro,j),&xx,&gg);
			if (r>=0)  r = pfield_sum_ch(ch,1,EXPM(e,i,j),SIGNM(e,i,j),
					1,xx,gg,&(EXPM(e,i,j)),&(SIGNM(e,i,j)));
		}		/* (arithmetic errors are reported for ch>0) */
	pfield_setzeroexp(&(EXPM(e,ro,co)));  SIGNM(e,ro,co) = 1;
			/* the pivoting exchanges the row and column ro,co id's */
	i = COLSID(e,co);  COLSID(e,co) = ROWSID(e,ro);  ROWSID(e,ro) = i;
	
	if (r<0 && ch>0)  {PROGERROR("pivoting %d,%d in %p failed!",ro,co,e); EMATDEBUG(0,e,"  !!\t");}
	return r<0?-1:0;
}
  

/**
 * This function multiplies a submatrix on rows r1,...,r2-1 and columns c1,...,c2-1 of
 * the given matrix e by x,g.
 * (The modified matrix is stored in the same e.)
**/

void	ematrix_multiply_rc(ematrix *e, int r1, int r2, int c1, int c2, int sx, exp_t x, sign_t g) {
	int	i,j, ret=0;
	
	if (g==1) if (pfield_iszeroexp(x))  return;
	for (i=r1; i<r2; i++)  for (j=c1; j<c2; j++) {
		ret = pfield_mul(1,EXPM(e,i,j),SIGNM(e,i,j),sx,x,g,
				&EXPM(e,i,j),&SIGNM(e,i,j)) <0?-1:ret;
	}
	if (ret<0)  {PROGERROR("multiplying %d_%d x %d_%d submatrix of %p failed!",r1,r2,c1,c2,e);}
}



/**
 * This function computes a new matrix obtained from em by deleting (dc==1) or contracting
 * (dc==0) the element so (row if tr==0, column if tr==1).
 * (Rows and columns are indexed from 0.)
 * If necessary, the element for removing is first pivoted in the matrix.
 * Contracting/deleting a loop/coloop results in the oposite operation.
 * The new matrix is returned as a new copy.
 * 
 * The next function does the same for the element of the given id label (not line-number).
**/

ematrix*	ematrix_removemat_ext(ematrix *em, int dc, int tr, int so) {
	ematrix	*ee;
	int	j;
	
	ee = ematrix_copy(em);
	if (tr)  ematrix_transpose(ee);
	if (so<0 || so>=ROWSM(ee))  {PROGERROR("line number out of matrix size %d !!",so);}
	if (tr==dc) {		/* the easy removal with the right line */
		ematrix_remove_row(ee,so);
	} else {		/* must pivot before removing from the other line dir. */
		for (j=0; j<COLSM(ee) && SIGNM(ee,so,j)==0; j++) ;
		if (j<COLSM(ee)) {
			ematrix_pivot(ee,so,j);
			ematrix_remove_col(ee,j);
		} else  ematrix_remove_row(ee,so);
	}
	if (tr)  ematrix_transpose(ee);
	return ee;
}

ematrix*	ematrix_removematid_ext(ematrix *em, int dc, int id) {
	ematrix		*ee;
	int		i,k,t;
	
	for (i=0,k=t=-1; i<ROWSM(em)+COLSM(em); i++)
		if ((i<ROWSM(em)? ROWSID(em,i):COLSID(em,i-ROWSM(em)))==id) {
			if (k>=0)  {PROGERROR("Duplicated id's %d in lines %d, %d!",id,k,i); EMATDEBUG(1,em,"\t!\t"); break;}
			t = (i>=ROWSM(em));  k = (t? i-ROWSM(em):i);
		}
	ee = ematrix_removemat_ext(em,dc,t,k);
	return ee;
}



/**
 * This is for computation of a determinant of the given matrix.
 * The computation creates a local copy of the matrix not to change the original.
 * The result is stored into given pointers x,g;  and 0, or -1 for an error, is returned.
 * (Errors include non-square matrix, out of in pfield numbers, or invalid entries.)
 * A variant is provided that does not print error messages for ch<=0.
**/

int	ematrix_determinant_ch(int ch, ematrix *e, exp_t *x, sign_t *g) {
	ematrix		*e1;
	int		i,j, ret = 0;
	exp_t		xx, xm;
	sign_t		gg, gm;
	
	if (ROWSM(e)!=COLSM(e) || ROWSM(e)==0) {
		if (ch>0)  {PROGERROR("Not computing determinant of a non-square matrix  %d x %d",ROWSM(e),COLSM(e));}
		ret = -1;
	} else {
		pfield_setzeroexp(&xx);  gg = 1;
		e1 = ematrix_copy_to(e,ROWSM(e),COLSM(e));
		/**
		 * How do we compute the determinant by the Gauss elimination:
		 * We find a nonzero entry i,j in the last row (if it exists),
		 * and then we pivot on that entry i,j.
		 * In this way the (original) column j becomes a unit row, not
		 * affecting the determinant value, and so the lines i,j are removed.
		 * The smaller determinant is then computed by recursion.
		 * Finally, we must multiply the recurrent value by the original
		 * entry at i,j since pivoting divided the determinant by this entry,
		 * and we must also add a proper sign by (i+j)%2.
		**/
		while (ROWSM(e1)>0 && ret>=0) {
			i = ROWSM(e1)-1;
			for (j=COLSM(e1)-1; j>=0; j--)
				if (SIGNM(e1,i,j)!=0)  break;
			if (j<0) {		/* zero line means 0 determinant value */
				gg = 0;  break;
			}
			xm = EXPM(e1,i,j);	/* multiply by the entry value before pivoting */
			gm = SIGNM(e1,i,j) * ((i+j)%2?-1:1);
			ret = pfield_mul_ch(ch,1,xm,gm,1,xx,gg,&xx,&gg);
			if (ROWSM(e1)>1 && ret>=0) {
				ret = ematrix_pivot_check(ch,e1,i,j);
				ematrix_remove(e1,i,j);
			} else  break;
		}
		dispose_ematrix(e1);
#ifndef FASTPROG
		if (ret>=0 && IFRANDDEBUGLESS(222) && ch>0) {	/* testing correct computation (only if supposed in pfield) */
			e1 = ematrix_copy(e);
			if (ROWSM(e1)==2)  ematrix_determinant_2x2(e1,0,1,0,1,&xm,&gm);
			else  if (ROWSM(e1)>2) {
				ematrix_swaprows(e1,1,2);  ematrix_swapcols(e1,COLSM(e1)-1,COLSM(e1)-2);
				ematrix_transpose(e1);  ematrix_determinant_ch(-1,e1,&xm,&gm);
			}		/* (call with ch==-1 to prevent inf recursion) */
			if (ROWSM(e1)>=2)  if (!pfield_isequal(xx,gg,xm,gm)) {
				PROGERROR("incorrect determinant computation!");  EMATDEBUG(0,e1,"!\t");
				DEBUG(0," - det = %s ,  detch = %s\n",pfield_printvalue(20,xx,gg),pfield_printvalue(20,xm,gm)); }
			dispose_ematrix(e1);
		}
#endif
	}
	if (ret<0)  gg = 0;	/* returns 0 after an error */
	if (x || g)  pfield_tostandform(&xx,&gg);
	if (x)  *x = xx;  if (g)  *g = gg;
	if (ch>0 && ret<0)  {PROGERROR("Computing determinant of %p [%dx%d] has failed!",e,ROWSM(e),COLSM(e));}
	return ret;
}


/**
 * This is a simple computation of a 2 x 2 subdeterminant on r1,r2 x c1,c2 lines of e.
 * See above for the full determinant computation.
**/

int	ematrix_determinant_2x2ext(int ch, ematrix *e, int r1, int r2, int c1, int c2,
								exp_t *x, sign_t *g) {
	exp_t	x1,x2,xx;
	sign_t	g1,g2,gg;
	int	ret = 0;
	
	ret = pfield_mul_ch(ch,1,EXPM(e,r1,c1),SIGNM(e,r1,c1),
					1,EXPM(e,r2,c2),SIGNM(e,r2,c2),&x1,&g1);
	if (ret>=0)  ret = pfield_mul_ch(ch,1,EXPM(e,r1,c2),SIGNM(e,r1,c2),
					1,EXPM(e,r2,c1),SIGNM(e,r2,c1),&x2,&g2);
	if (ret>=0)  ret = pfield_sum_ch(ch,1,x1,g1,1,x2,-g2,&xx,&gg);
	
	if (x)  *x = xx;  if (g)  *g = gg;
	if (ret<0 && ch>0)  {PROGERROR("computing 2x2 determinant %d,%d x %d,%d of %p failed! %s %s; %s %s",r1,r2,c1,c2,e,
				pfield_pvalue(8,EXPM(e,r1,c1),SIGNM(e,r1,c1)),pfield_pvalue(8,EXPM(e,r1,c2),SIGNM(e,r1,c2)),pfield_pvalue(8,EXPM(e,r2,c1),SIGNM(e,r2,c1)),pfield_pvalue(8,EXPM(e,r2,c2),SIGNM(e,r2,c2)));
				EMATDEBUG(0,e," !D!\t");}
	return ret;
}



/**
 * This function is called to "import" the matrix ee into the current pfield
 * by the current selected translation (see pfield_translation_ext() in pfield.c).
 * Importing means translating the exponents of the entries of ee into powers
 * of the images of generators, as defined by the translation.
 * 
 * The imported entries are stored back to the given matrix ee.
 * The return value is >=0 for success, and <0 for failure.
**/

int	ematrix_import_ext(int ch, ematrix *ee) {
	int	i,j, g, r = 0;
	
	for (i=0; i<ROWSM(ee) && r>=0; i++)
		for (j=0; j<COLSM(ee) && r>=0; j++) {
			g = SIGNM(ee,i,j);
			r = pfield_translation_ch(ch,&EXPM(ee,i,j),&SIGNM(ee,i,j))<0? -1:r;
			if (g && !SIGNM(ee,i,j)) { r = -1;  if (ch>=0) USERERROR("Pfield translation to zero 0."); }
		}
	if (ch>0 && r<0)  {PROGERROR("Translation of the matrix %p has failed!",ee);}
	return r; 
}


/**
 * This function is called to apply the pfield endomorphism #le to the matrix ee.
 * (See pfield_endomorph_ext() in ../pfield/pfmore.c).
 * The translated entries are stored back to the given matrix ee.
 * The return value is >=0 for success, and <0 for failure.
**/

int	ematrix_pfendomorph_ext(int ch, int le, ematrix *ee) {
	int	i,j, r = 0;
	
	if (le<0 || le>=pfield_endomorph_number())  {PROGERROR("Pfield endomorphism index out of bounds! %d",le); r = -1;}
	for (i=0; i<ROWSM(ee) && r>=0; i++)
		for (j=0; j<COLSM(ee) && r>=0; j++) {
			r = pfield_endomorph(le,&EXPM(ee,i,j),&SIGNM(ee,i,j))<0? -1:r;
		}
	if (ch>0 && r<0)  {PROGERROR("Pfield endomorphism of the matrix %p has failed!",ee);}
	return r; 
}

















/******************	Simple matrix tests	*******************/
/******************************************************************/


#undef CURDLEV
#define CURDLEV		7


/**
 * This function tests whether the two lines s1,s2 of the matrices e1,e2 are parallel up to
 * the excluded two rows no1,no2 (use -1 to not exclude them).
 * The return value is either -1 for not parallel, or r>=0 if s1,s2 are parallel up to no1,no2.
 * (In case when one of the lines is 0, the answer is parallel as well.)
**/

int	ematrix_twoparallel_rc(ematrix *e1, int tr, int s1, ematrix *e2, int s2, int no1, int no2) {
	int	i,f,z, ret;
	exp_t	de, ae;
	sign_t	dg, ag;
	
	pfield_setzeroexp(&de);  dg = 0;
	if (!tr)  ematrix_transpose(e1);
	if (!tr && e1!=e2)  ematrix_transpose(e2);
	if (ROWSM(e1)!=ROWSM(e2))
		{PROGERROR("cannot test parallel lines that have different lengths %d!=%d",ROWSM(e1),ROWSM(e2));}
	
	ret = -1;  f = 0;		/* tests all entries of s1,s2 except excluded no1,no2 */
	for (i=0; i<ROWSM(e1); i++)  if (i!=no1 && i!=no2) {
		if (SIGNM(e1,i,s1)==0 && SIGNM(e2,i,s2)==0)  continue;
					/* computes the fraction - OK to "divide by 0" which results in 0 */
		pfield_mul(1,EXPM(e1,i,s1),SIGNM(e1,i,s1),-1,EXPM(e2,i,s2),SIGNM(e2,i,s2), &ae,&ag);
		z = (SIGNM(e2,i,s2)==0? -1:1);		/* (indicates division by 0: z==-1) */
		if (f==0) {
			de = ae;  dg = ag;  f = z;	/* remembers the difference between the first nonzeros */
		} else {
			if (f!=z || !pfield_isequal(ae,ag,de,dg))  break;
		}			/* (must distinguish when s2 is a 0-multiple of s1 as well) */
	}
	if (i>=ROWSM(e1))  ret = 0;	/* if the above cycle was not broken for j, then parallel */
	if (!tr)  ematrix_transpose(e1);
	if (!tr && e1!=e2)  ematrix_transpose(e2);
	
#ifndef FASTPROG			/* an extra test for symmetric result... */
	if (IFRANDDEBUGLESS(222) && !inrecur) {
		inrecur = 1;		/* (to avoid inf recursion) */
		if ((ret<0)!=(ematrix_twoparallel_rc(e2,tr,s2,e1,s1,no2,no1)<0))
			{PROGERROR("Non-symmetric twoparallel computation");}
		inrecur = 0; }
#endif
	return ret;
}


/**
 * This function tests whether the line so of the matrix e is all-zero up to
 * the excluded three rows no1,no2,no3 (use -1 to not exclude them).
 * If nz==1, then the excluded rows must be nonzero to pass.
 * If ct>0, then up to ct arbitrary nonzero entries (in addition to no?) pass.
 * The return value is either -i-1 for other nonzero, or 1 for all but no?,ct zero entries.
**/

int	ematrix_linezero_ext(ematrix *e, int tr, int so, int no1, int no2, int no3, int nz, int ct) {
	int	i,ret;
	
	DEBUG(CURDLEV+1,"tr=%d, so=%d, no1=%d, no2=%d, no3=%d, nz=%d, ct=%d\n",tr,so,no1,no2,no3,nz,ct);
	if (!tr)  ematrix_transpose(e);
	if (ct<0)  ct = 0;
	for (i=0; i<ROWSM(e) && ct>=0; i++)
		if (i!=no1 && i!=no2 && i!=no3 && SIGNM(e,i,so)!=0)  ct--;
	ret = (ct<0? -i-1:1);
	if (nz>0) if ((no1>=0? SIGNM(e,no1,so)==0:0) || (no2>=0? SIGNM(e,no2,so)==0:0) ||
		(no3>=0? SIGNM(e,no3,so)==0:0))  ret = -1;
	if (!tr)  ematrix_transpose(e);
	DEBUG(CURDLEV+1,"ct=%d, ret=%d\n",ct,ret);
	return ret;
}


/**
 * This function computes whether the given matrix e has a column parallel to the given
 * column (tr==1) so of a matrix eo, or similarly for the rows (tr==0).
 * If e==eo, then the column so is not tested against itself (of course).
 * The parameter ch is not used now.
 * The return value is either -1 for no parallel, or r>=0 if a column r is parallel to so.
**/

int	ematrix_lineparallel_ext(int ch, ematrix *e, int tr, ematrix *eo, int so) {
	int	j, ret;
	
	if (!tr)  ematrix_transpose(e);
	if (!tr && e!=eo)  ematrix_transpose(eo);
	junk = ch;  ret = -1;
	for (j=0; j<COLSM(e) && ret<0; j++)
		if (j!=so || eo!=e) {
			ret = ematrix_twoparallel_rc(e,1,j,eo,so,-1,-1);
			if (ret>=0)  ret = j;
		}
	if (!tr)  ematrix_transpose(e);
	if (!tr && e!=eo)  ematrix_transpose(eo);
	return ret;
}


/**
 * This function computes whether the given matrix e has two parallel columns or rows.
 * The return value is either -1 for no parallel, or r=0 if two rows are parallel,
 * or r=1 if two columns are parallel, or r=2 if both happen.
 * 
 * An analogous function is provided for testing any zero line in the matrix.
**/

int	ematrix_parallel_any(ematrix *e) {
	int	i,t, ret;
	
	ret = -1;
	for (t=0; t<2; t++)  for (i=0; i<(t?COLSM(e):ROWSM(e)); i++)
		if (ematrix_lineparallel(e,t,e,i)>=0) {
			ret = (ret<0? t:2);
			break;
		}
	return ret;
}

int	ematrix_linezero_any(ematrix *e) {
	int	i,t, ret;
	
	ret = -1;
	for (t=0; t<2; t++)  for (i=0; i<(t?COLSM(e):ROWSM(e)); i++)
		if (ematrix_linezero(e,t,i)>=0) {
			ret = (ret<0? t:2);
			break;
		}
	return ret;
}


/**
 * This function returns the rank of the given matrix e (as an integer).
 * If lrk>=0 is given, then only rank up to lrk is computed, even if it is higher.
 * If the matrix is not representable and this is found here, then the result is -1.
 * (The matrix e is not modified here.)
**/

int	ematrix_matrank_ext(ematrix *e, int lrk) {
	int	i,j,c, rk;
	ematrix	*ee;
	
	ee = ematrix_copy(e);
	rk = 0;
	while (ROWSM(ee)>0 && COLSM(ee)>0 && rk>=0 && (lrk<0 || rk<lrk)) {
		j = COLSM(ee)-1;
		for (i=ROWSM(ee)-1; i>=0 && SIGNM(ee,i,j)==0; i--) ;
		if (i<0) {		/* zero column is ignored */
			ematrix_remove_col(ee,j);
		} else {
			rk++;		/* nonzero column is pivoted and the row-column removed */
			c = ematrix_pivot_check(0,ee,i,j);
			ematrix_remove(ee,i,j);
			if (c<0)  rk = -1;	/* for the case that the matrix is not representable */
		}
	}
	dispose_ematrix(ee);
#ifndef FASTPROG		/* recursive debug test of the result... */
	if (IFRANDDEBUGLESS(333) && rk>=0 && lrk<0 && ROWSM(e)>1 && COLSM(e)>1) {
		ee = ematrix_copydual(e);
		i = RANDOM()%ROWSM(ee);  j = RANDOM()%COLSM(ee);
		if (SIGNM(ee,i,j)==0)  c = 0;
		else { c = 1;  ematrix_pivot(ee,i,j);  ematrix_remove(ee,i,j); }
		c += ematrix_matrank_ext(ee,rk+1);
		if (c!=rk)  {PROGERROR("Wrong matrank computation, ret=%d check=%d",rk,c); EMATDEBUG(0,e,"!\t");}
		dispose_ematrix(ee);
	}
#endif
	return rk;
}


/**
 * This function returns whether the three lines so1,so2,so3 (row/col by tr) of the matrix e
 * are dependent (1 for dep, -1 for not).
 * If the matrix is not representable, then the return value is like for dependent.
 * (This is for testing triangles and triads in the matrix, but only without parallel/serial elements.)
**/

int	ematrix_triax_rc(ematrix *e, int tr, int so1, int so2, int so3) {
	ematrix	*re;
	int	r;
	
	if (!tr)  ematrix_transpose(e);
	re = ematrix_refer(e,0,ROWSM(e),-1,-1);
	ematrix_refadd_col(re,so1);
	ematrix_refadd_col(re,so2);
	ematrix_refadd_col(re,so3);
	r = ematrix_matrank(re);
	dispose_ematrix(re);
	if (!tr)  ematrix_transpose(e);
	return r<3? 1:-1;
}


/**
 * This function returns the matroid-set rank of the element subset given by
 * lines of refering re in the matroid given by e.
**/

int	ematrix_setrank(ematrix *e, ematrix *re) {
	ematrix		*rr;
	int		k;
	
	if (re?REFMAT(re)!=e:1)  {PROGERROR("Can compute matroid set rank only for a refering matrix re->e.");}
	rr = ematrix_refextract_xrow(e,re);
	k = ROWSM(re)+ematrix_matrank(rr);
	dispose_ematrix(rr);
	return k;
}














/******************	Closures and Separations	*******************/
/**************************************************************************/


#undef CURDLEV
#define CURDLEV		7


/**
 * This function returns the (co-)closure of the given submatrix re in the matrix e.
 * "Closure" means all matroid elements (rows and columns of e) that belong to the subspace
 * spanned by the subset corresponding to the lines of re.
 * "Co-closure" is the closure in the dual matroid (transposed matrix).
 * Closure is computed for tr==1, and co-closure for tr==0.
 * The matrix must be properly representable.
 * The return value is a new refering matrix (to e) containing all lines of re and of the
 * closure.
**/

ematrix*	ematrix_closure_ext(ematrix *em, ematrix *re, int tr) {
	ematrix		*e, *rre, *e2,*re2,*rre2;
	int		i,ii,j,k, *cm,*rm, cmstack[30];
	
	if (ISREFMAT(em) || !ISREFMAT(re) || REFMAT(re)!=em)  {PROGERROREXIT("call only for a refering matrix re -> e");}
	if (ISTRANSPM(em)!=ISTRANSPM(re))  {PROGERROR("call only for matrices e,re of the same transpose");}
	if (!tr) { ematrix_transpose(em);  ematrix_transpose(re); }
	e = ematrix_copy(em);	/* (the matrix e is pivoted below!) */
	
		/**
		 * We prepare two arrays cm[],rm[] for marking matrix lines that belong
		 * to re (xm=-2), and for marking columns that were pivoted from rows
		 * not in re (cm[ii]==i iff the row i was pivoted into the column ii).
		**/
	if (ROWSM(e)+COLSM(e)<25)  cm = cmstack;
	else  cm = MMALLOC((ROWSM(e)+COLSM(e)+5)*sizeof(cm[0]));
	rm = cm+COLSM(e)+2;
	for (i=0; i<ROWSM(e); i++)  rm[i] = 0;
	for (i=0; i<ROWSM(re); i++)  rm[GETREFMROW(re,i)] = -2;
	for (j=0; j<COLSM(e); j++)  cm[j] = -1;
	for (j=0; j<COLSM(re); j++)  cm[GETREFMCOL(re,j)] = -2;
		/**
		 * We try to pivot the columns of re to get as much as possible of them 
		 * into rows of the matrix.
		 * The remaining columns must have all zeros in the rows outside re.
		**/
	for (j=0; j<COLSM(re); j++) {
		ii = GETREFMCOL(re,j);
		for (i=0; i<ROWSM(e) && cm[ii]==-2; i++)
		  if (rm[i]==0 && SIGNM(e,i,ii)!=0) {
			ematrix_pivot(e,i,ii);
			rm[i] = -2;  cm[ii] = i;
		}
	}
		/**
		 * Finally, the closure contains all columns (of the new pivoted e)
		 * that are dependent on the pivoted rows of new re.
		 * We use cm[] to recover the original positions of these columns in em.
		**/
	rre = ematrix_refer_all(re);
	for (ii=0; ii<COLSM(e); ii++)  if (cm[ii]>=-1) {
		for (i=k=0; i<ROWSM(e); i++)
			if (rm[i]>=-1 && SIGNM(e,i,ii)!=0)  k = 1;
		if (k)  continue;
		if (cm[ii]==-1)  ematrix_drefadd_col(rre,ii);
		else  ematrix_drefadd_row(rre,cm[ii]);
	}
	
	rre2 = re2 = e2 = NULL;
#ifndef FASTPROG		/* some extra paranoic tests - double closure and pivoted closure */
	if ((IFRANDDEBUGLESS(33)) && !inrecur && ROWSM(e)>1 && COLSM(e)>1) {
		inrecur = 1;		/* (to avoid inf recursion) */
		rre2 = ematrix_closure_ext(em,rre,1);
		if (!ematrix_isequal(rre2,rre))  {PROGERROR("wrong double-closure computation");}
		dispose_ematrix(rre2);
		e2 = ematrix_copy(em);	/* trying to pivot the matrix not touching the re-lines */
		re2 = ematrix_refer_all(re);  ematrix_rerefer(re2,e2);
		i = ROWSM(rre)>ROWSM(re)? GETREFMROW(rre,ROWSM(rre)-1): RANDOM()%ROWSM(e2);
		j = RANDOM()%COLSM(e2);
		for (ii=0; ii<ROWSM(re2); ii++)  if (GETREFMROW(re2,ii)==i) { i = -1;  break; }
		for (ii=0; ii<COLSM(re2); ii++)  if (GETREFMCOL(re2,ii)==j) { j = -1;  break; }
		if (i>=0 && j>=0 && SIGNM(e2,i,j)!=0) {
			ematrix_pivot(e2,i,j);
			rre2 = ematrix_closure_ext(e2,re2,1);
			if (ROWSM(rre2)+COLSM(rre2)!=ROWSM(rre)+COLSM(rre)) {
				PROGERROR("wrong closure computation after pivoting");
				EMATDEBUG(0,e," e\t"); EMATDEBUG(0,rre," +r\t");
				EMATDEBUG(0,e2," e2\t"); EMATDEBUG(0,re2," in\t"); EMATDEBUG(0,rre2," !r!\t"); }
			dispose_ematrix(rre2);
		}
		dispose_ematrix(re2);  dispose_ematrix(e2);
		inrecur = 0;
	}
#endif
	if (cm && cm!=cmstack)  FREE(cm);
	dispose_ematrix(e);
	if (!tr) {
		ematrix_transpose(em);  ematrix_transpose(re);
		ematrix_transpose(rre);	/* (the output matrix must be transposed back as well!) */
	}
	return rre;		/* the closure is returned as a new refering matrix */
}


#if	0	/* this is an old implementation no longer used... */

ematrix*	ematrix_closure_ext(ematrix *e, ematrix *re, int tr) {
	ematrix		*ex, *ee,*ee2,*re2, *rre;
	int		i,ii,j, r0,r1;
	
	if (ISREFMAT(e) || !ISREFMAT(re) || REFMAT(re)!=e)  {PROGERROREXIT("call only for a refering matrix re -> e");}
	if (ISTRANSPM(e)!=ISTRANSPM(re))  {PROGERROR("call only for matrices e,re of the same transpose");}
	if (!tr) {
		ematrix_transpose(e);  ematrix_transpose(re);
	}
	rre = ematrix_refer_all(re);
	ex = ematrix_refextract_xrow(e,re);
	r0 = ematrix_matrank(ex);
		/**
		 * A column of e belongs to the closure of re iff it is dependent on the submatrix
		 * of e on the columns as re and on the rows not in re (matrix ex here).
		 * So we add corresponding part of a column to the matrix ex and compare the rank.
		**/
	if (COLSM(re)<COLSM(e))  for (j=0; j<COLSM(e); j++) {
		for (ii=0; ii<COLSM(ex) && GETREFMCOL(ex,ii)!=j; ii++) ;
		if (ii<COLSM(ex))  continue;		/* (column already in re) */
		ematrix_drefadd_col(ex,j);
		r1 = ematrix_matrank(ex);
		ematrix_remove_col(ex,COLSM(ex)-1);
		if (r1<0 || r0<0)  {PROGERROR("the matrix must be representable here");}
		if (r1<r0)  {PROGERROR("wrong rank computation of added column");}
		if (r1==r0)  ematrix_drefadd_col(rre,j);
	}
		/**
		 * We then find all rows not in re but belonging to the closure.
		 * We actually compute with a row like with an extra column parallel to that row
		 * (so we are appending unit columns to a copy of the matrix ex).
		 * The rest is the same as previously.
		**/
	ee = ematrix_copy(ex);
	ematrix_resetid(ee);		/* (must prevent duplicate line id's for append) */
	ee = ematrix_append_col(ee);
	if (ROWSM(re)<ROWSM(e))  for (i=0; i<ROWSM(e); i++) {
		for (ii=0; ii<ROWSM(ex) && GETREFMROW(ex,ii)!=i; ii++) ;
		if (ii>=ROWSM(ex))  continue;		/* (row already in re - not in ex) */
		for (j=0; j<ROWSM(ee); j++)  SETEXSIGMZERO(ee,j,COLSM(ee)-1);
		SIGNM(ee,ii,COLSM(ee)-1) = 1;		/* must index rows according to ex ! */
		r1 = ematrix_matrank(ee);
		if (r1<0 || r0<0)  {PROGERROR("the matrix must be representable here");}
		if (r1<r0)  {PROGERROR("wrong rank computation of added row");}
		if (r1==r0)  ematrix_drefadd_row(rre,i);
	}
	dispose_ematrix(ee);
	dispose_ematrix(ex);
	re2 = ee2 = NULL;
#ifndef FASTPROG		/* some extra paranoic tests - double closure and pivoted closure */
	if ((IFRANDDEBUGLESS(111)) && !inrecur && ROWSM(e)>1 && COLSM(e)>1) {
		inrecur = 1;		/* (to avoid inf recursion) */
		ee = ematrix_closure_ext(e,rre,1);
		if (!ematrix_isequal(ee,rre))  {PROGERROR("wrong double-closure computation");}
		dispose_ematrix(ee);
		ee2 = ematrix_copy(e);	/* trying to pivot the matrix not touching the re-lines */
		re2 = ematrix_refer_all(re);  ematrix_rerefer(re2,ee2);
		i = ROWSM(rre)>ROWSM(re)? GETREFMROW(rre,ROWSM(rre)-1): RANDOM()%ROWSM(ee2);
		j = RANDOM()%COLSM(ee2);
		for (ii=0; ii<ROWSM(re2); ii++)  if (GETREFMROW(re2,ii)==i) { i = -1;  break; }
		for (ii=0; ii<COLSM(re2); ii++)  if (GETREFMCOL(re2,ii)==j) { j = -1;  break; }
		if (i>=0 && j>=0 && SIGNM(ee2,i,j)!=0) {
			ematrix_pivot(ee2,i,j);
			ee = ematrix_closure_ext(ee2,re2,1);
			if (ROWSM(ee)+COLSM(ee)!=ROWSM(rre)+COLSM(rre)) {
				PROGERROR("wrong closure computation after pivoting");
				EMATDEBUG(0,e," e\t"); EMATDEBUG(0,rre," +r\t");
				EMATDEBUG(0,ee2," e2\t"); EMATDEBUG(0,re2," in\t"); EMATDEBUG(0,ee," !r!\t"); }
			dispose_ematrix(ee);
		}
		dispose_ematrix(re2);  dispose_ematrix(ee2);
		inrecur = 0;
	}
	//DEBUG(1,"## closure results ===================================\n");
	//EMATDEBUG(1,e,"\t\tin\t");EMATDEBUG(1,re,"\t\tre\t");EMATDEBUG(1,ex,"\t\tex\t");EMATDEBUG(1,rre,"\t\trre\t");
#endif
	if (!tr) {
		ematrix_transpose(e);  ematrix_transpose(re);
		ematrix_transpose(rre);	/* (the output matrix must be transposed back as well!) */
	}
	return rre;		/* the closure is returned as a new refering matrix */
}
#endif		/* this is an old implementation no longer used... */




/**
 * This function computes the connectivity of the given submatrix re in e,                                       
 * i.e. the separation value for the separation (re,e-re) minus one.
 * We use the formula for (A,B) as lambda(A) = r(A)+r(B)-r(M).
 * So the value of this function is 0 for 1-separation, 1 for exact 2-separation, etc...
 * The above formula can be simply implemented via computing the ranks of the submatrices
 * opposite to those of A and of B (rxr and rxc below).
 * 
 * If lsp>=0 is given, then the function computes the separation value only up to lsp
 * (even when it may be higher).
 * The separation value is returned back.
**/

int	ematrix_whatsep_ext(ematrix *e, ematrix *re, int lsp) {
	ematrix		*rxr,*rxc,*rr;
	int		sp;
	
	if (ISREFMAT(e) || !ISREFMAT(re) || REFMAT(re)!=e)  {PROGERROREXIT("call only for a refering matrix re -> e");}
	if (ISTRANSPM(e)!=ISTRANSPM(re))  {PROGERROR("call only for matrices e,re of the same transpose");}
	rr = rxc = NULL;
	rxr = ematrix_refextract_xrow(e,re);
	sp = ematrix_matrank_ext(rxr,(lsp<0?-1:lsp));
	if (lsp<0 || sp<lsp) {
		rxc = ematrix_refextract_xcol(e,re);
		sp += ematrix_matrank_ext(rxc,(lsp<0?-1:lsp-sp));
	}
#ifndef FASTPROG
	DEBUG(CURDLEV+1,"Computed separation value %d (bound %d) for:\n",sp,lsp);
	EMATDEBUGS(CURDLEV+1,rxr,"\t\t\t"); EMATDEBUGS(CURDLEV+1,rxc,"\t\t\t");
	if (IFRANDDEBUGLESS(222) && lsp<0 && sp>=0) {
		rr = ematrix_refextract_xall(e,re);
		if (sp!=ematrix_whatsep_ext(e,rr,sp+1))  {PROGERROR("wrong separation value for the complement! %d!=%d",sp,ematrix_whatsep(e,rr));
				EMATDEBUG(1,re,"\t\tre\t");EMATDEBUG(1,rr,"\t\trr\t");}
		dispose_ematrix(rr); }
#endif
	dispose_ematrix(rxr);  if (rxc) dispose_ematrix(rxc);
	return sp;
}












































