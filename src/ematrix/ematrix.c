
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
 * Specific functions are provided for allocating and freeing such matrices.
 * Macros EXPM, SIGNM, ROWSM, COLSM, etc should be used to access elements of a matrix.
 * (More access macros in matrix.h ...)
 * 
 * Additionally, a matrix may just refer to a submatrix of another existing matrix.
 * In such case it does not have its own entries, but reference arrays for rows/columns.
 * (This gets little tricky when you try to refer an already refering matrix...)
 * 
 * Look for more information and for the easy-access macros in ../include/ematrix.h.
 * 
**/





#include "macek.h"
#include "emx.h"








/******************	Matrix allocation and clearing	*******************/
/**************************************************************************/


#define CURDLEV		7


/**
 * These are allocation a freeing functions for our matrices.
 * Each matrix uses just one piece of memory, the dynamic arrays are referred to the mem
 * just after the base structure.
 * Notice that ematrix_alloc already sets the m->dlength field(!).
 * 
 * Moreover, we implement "saving" of allocated matrices for later use, which is supposed
 * to make allocation faster... ("recycling"?)
 * The saved matrices are stored in two-way lists starting at hash-indices of their sizes.
**/

#ifndef EM_ALCMAX
#define EM_ALCMAX	70000
#endif
#define	SIZEMOD	(11+EM_ALCMAX/10)
struct	alchain { int size;  void *prev;  void *next; } ;
#define ALCSIZE(p)	((struct alchain*)(p))->size
#define ALCPREV(p)	((struct alchain*)(p))->prev
#define ALCNEXT(p)	((struct alchain*)(p))->next
static void	*alcpoints[SIZEMOD];
static int	alcstart = 1, alcnum = 0;
static int	alctotal = 0;

#define	EM_MAXMXSIZE	500000	/* (no matrix may be bigger than this size) */

void*	ematrix_alloc(int sz, int clr) {
	int	i;
	void	*dt, *x;
	
	dt = NULL;
	if (sz>EM_MAXMXSIZE)  {PROGERROREXIT("this matrix is too big for me %d>%d",sz,EM_MAXMXSIZE);}
	if (sz<(int)sizeof(ematrix))  {PROGERROREXIT("probably too small matrix %d<%d",sz,(int)sizeof(ematrix));}
	if (alcstart==1) {		/* an initial clearing of indices */
		for (i=0; i<SIZEMOD; i++)  alcpoints[i] = NULL;
		alcstart = 0;
	}
	for (i=0, x=alcpoints[sz%SIZEMOD]; x && i<5 && !dt; i++, x=ALCNEXT(x))
	  if (ALCSIZE(x)==sz) {		/* a saved matrix of the right size found here */
		dt = x;  alcnum--;
		if (i==0)  alcpoints[sz%SIZEMOD] = ALCNEXT(x);
		else  ALCNEXT(ALCPREV(x)) = ALCNEXT(x);
		if (ALCNEXT(x))  ALCPREV(ALCNEXT(x)) = ALCPREV(x);
		if (clr)  MEMSET(x,0,sz);
	}
	if (!dt) {			/* if no saved matrix was found above, then allocate a new one */
		if (clr)  dt = CALLOC(sz+2+sizeof(long),1);
		else  dt = MALLOC(sz+2+sizeof(long));
	}
	if (!dt)  {PROGERROREXIT("cannot allocate new ematrix structure");}
	((ematrix*)dt)->dlength = sz;
	EMMAGIC(dt) = EMMAGICNUM;	/* ("magic number" kept against mem overflows in the matrix - depends on sz !) */
	((ematrix*)dt)->name = NULL;  ((ematrix*)dt)->data = NULL;
	alctotal++;
	return dt;
}

void	ematrix_free(void *dt) {
	int	sz;
	
	if (!EMMAGICTEST(dt))  {PROGERROREXIT("damaged matrix memory or incorrect pointer %p",dt);}
	sz = ((ematrix*)dt)->dlength;
	((ematrix*)dt)->dlength = 0;		/* clearing the matrix memory - magic value refered by dlength */
	alctotal--;
	if (alcnum<=EM_ALCMAX) {		/* saving the matrix for later use, unless too many saved */
		ALCSIZE(dt) = sz;  alcnum++;
		ALCNEXT(dt) = alcpoints[sz%SIZEMOD];
		alcpoints[sz%SIZEMOD] = dt;
		if (ALCNEXT(dt))  ALCPREV(ALCNEXT(dt)) = dt;
		ALCPREV(dt) = NULL;
	} else	FREE(dt);
}

/**
 * This function flushes all saved matrix allocations.
**/

void	ematrix_finishall(void) {
	int	i;
	void	*x,*xx;
	
	DEBUG(CURDLEV-1,"Freeing %d saved matrix allocations...\n",alcnum);
	if (alctotal>10)  DEBUG(CURDLEV-3,"There are %d matrices in the program not disposed yet!\n",alctotal);
	for (i=0; i<SIZEMOD; i++) {
		for (x=alcpoints[i]; x; x=xx) {
			xx = ALCNEXT(x);
			FREE(x);
		}
		alcpoints[i] = NULL;
	}
	alcnum = 0;
}

int	ematrix_alctotal(void) {
	return alctotal;
}


/**
 * We must use specific "memcpy()" implementation to copy matrices - to preserve the dlength
 * field of the destination matrix correctly, and to clear the matrix name.
**/

void	ematrix_memcpy(void *ds, void *sr, int len) {
	int	sz;
	
	sz = ((ematrix*)ds)->dlength;
	if (len>sz)  {PROGERROREXIT("Cannot copy longer than the allocated size %d>%d",len,sz);}
	MEMCPY(ds,sr,len);
	((ematrix*)ds)->dlength = sz;
	((ematrix*)ds)->name = NULL;  ((ematrix*)ds)->data = NULL;
}




/**
 * Here we allocate (stand-alone matrix only, no referring), and dispose a matrix.
 * A new matrix has an initial dimension and a maximal (allocated) dimension.
 * After creating, the matrix is filled with zeros.
**/

ematrix*	new_ematrix(int ro, int co, int mr, int mc) {
	ematrix	*e;
	long	i, sz;

	if (mr<ro)  mr = ro;
	if (mc<co)  mc = co;
	mr += mr%2;  mc += mc%2;	/* (we make the dimensions even for mem alignment) */
	sz = MEMALIGNSIZE(sizeof(*e));
	sz += MEMALIGNSIZE((mr*mc)*sizeof(e->exp[0])) + MEMALIGNSIZE((mr*mc)*sizeof(e->sgn[0]));
	sz += MEMALIGNSIZE(mr*sizeof(e->lid[0][0])) + MEMALIGNSIZE(mc*sizeof(e->lid[0][0]));
	e = ematrix_alloc(sz,1);	/* allocates cleared memory of size sz */
	
	if (e!=NULL) {
		e->transp = e->isref = e->numrefm = 0;
		ROWSMAX(e) = mr;  COLSMAX(e) = mc;
		ROWSM(e) = ro;  COLSM(e) = co;
		e->refm = NULL;  REFROW(e) = REFCOL(e) = NULL;
					/* the arrays of the matrix are stored in the subsequent memory */
		e->exp = (exp_t*)MEMALIGNPOINT(e+1);
		e->sgn = (sign_t*)MEMALIGNPOINT(e->exp+mr*mc);
		e->lid[0] = (short*)MEMALIGNPOINT(e->sgn+mr*mc);
		e->lid[1] = (short*)MEMALIGNPOINT(e->lid[0]+mr);
		for (i=0; i<ro; i++)  ROWSIDS(e)[i] = i+1;
		for (i=0; i<co; i++)  COLSIDS(e)[i] = -i-1;
		
		DEBUG(CURDLEV+1,"new matrix %p created of dims %d,%d and size %ld\n",e,ro,co,sz);
		if ((byte1*)(e->lid[1]+mc)-(byte1*)e >sz)  {PROGERROREXIT("wrong mem allocation %ld > %ld",(long)((byte1*)(e->lid[1]+mc)-(byte1*)e),sz);}
	}
	return e;
}

void	dispose_ematrix(ematrix *e) {
	
	if (!e)  return;
	if (e->dlength<=(long)sizeof(*e) || e->dlength>EM_MAXMXSIZE)  {PROGERROREXIT("disposing invalid matrix %p or disposing second time",e);}
	if (e->numrefm>0)  {PROGERROREXIT("disposing a matrix %p that is referred elsewhere %d",e,e->numrefm);}
#ifndef FASTPROG
	if ((!ISREFMAT(e))!=(!e->refm))  {PROGERROREXIT("disposing a matrix %p that refers incorectly",e);}
	if (e->numrefm<0 || e->numrefm>99999 || e->transp<0 || e->transp>1 || e->isref<0 || e->isref>1)
		{PROGERROREXIT("disposing a matrix %p that is in an incorrect state",e);}
	if (EMNAME(e))  if (strlen(EMNAME(e))>300)  {PROGERROR("wrong name in a matrix %p? : %s",e,EMNAME(e));}
#endif
	if (ISREFMAT(e))  REFMAT(e)->numrefm--;
	if (EMNAME(e))  FREE(EMNAME(e));
	if (EMDATA(e))  FREE(EMDATA(e));
	DEBUG(CURDLEV+2,"matrix %p disposed\n",e);
	ematrix_free(e);	/* (the e->dlength value is zeroed here inside) */
}













/******************	Basic matrix manipulation	*******************/
/**************************************************************************/


exp_t	macexp;		/* (for use in macros...) */


/**
 * This function computes matrix index for the entry r,c.
 * It contains some bound-checking, and supports evaluation for referred matrices.
 * The index is returned as an integer (to be used in the *exp and *sgn arrays).
 * For a really fast access the macros should be used instead !
**/

int	ematrix_function_indexm(ematrix *m, int r, int c) {
	static int	inerr=0;
	
	if (m->dlength<(long)sizeof(*m) || m->dlength>EM_MAXMXSIZE) {
		PROGERROREXIT("accessing invalid or freed matrix %p",m);
	}
#if CURDLEV<=DEBUGLEV
	if (!EMMAGICTEST(m)) {
		PROGERROREXIT("accessing invalid or freed matrix %p",m);
	}
	if (ROWSM(m)>ROWSMAX(m) || COLSM(m)>COLSMAX(m))  if (!inerr) {
		inerr = 1;  EMATDEBUG(0,m,"#!\t");  inerr = 0;
		PROGERROREXIT("dimensions out of bounds rmax=%d, cmax=%d in matrix %p",ROWSMAX(m),COLSMAX(m),m);
	}
	if (ISREFMAT(m)) if ( !REFMAT(m) || !REFROW(m) || !REFCOL(m) ||
			ISREFMAT(REFMAT(m)) || REFMAT(m)->numrefm<=0 ) {
		PROGERROREXIT("incorrect submatrix reference inside %p",m);
	}
	if (m->transp<0 || m->transp>1) {
		PROGERROR("incorrect transp %d value in %p", m->transp,m);
		return 0;
	}
#endif
	if (r<0 || c<0 || r>=ROWSM(m) || c>=COLSM(m))  if (!inerr) {
		inerr = 1;  EMATDEBUG(0,m,"#!\t");  inerr = 0;
		PROGERROR("index out of bounds r=%d, c=%d in matrix %p",r,c,m);
		return 0;
	}
	if (ISREFMAT(m))  return INDEXR(m,r,c);
	return INDEXM(m,r,c);
}


/**
 * This function checks distinct (ch==0), or "resets" (ch==1) the matrix line
 * identification labels to the standard form (like returned from new_ematrix()).
 * The return value is the largest absolute label,
 * or -1 if two same or zero labels are found.
**/

int	ematrix_lineids(int ch, ematrix *e) {
	int	i,j,x,y, *us, ustack[30];
	
	if (ISREFMAT(e)) { PROGERROR("Do not call for refering matrices %p",e); return -1; }
	if (ch<=0) {
		for (j=x=0; j<2; j++) {
			for (i=0; i<ROWSM(e); i++) { y = abs(ROWSID(e,i));  if (y>x)  x = y; }
			ematrix_transpose(e);
		}
		if (2*x<28)  us = ustack;
		else  us = MMALLOC((2*x+2)*sizeof(us[0]));
		for (i=0; i<=2*x; i++)  us[i] = 0;
		us[x] = 1;
		for (j=y=0; j<2; j++) {
			for (i=0; i<ROWSM(e); i++)
				if (++us[x+ROWSID(e,i)]>1)  y = 1;
			ematrix_transpose(e);
		}
		if (us!=ustack)  FREE(us);
		if (y)  x = -1;
	} else {
		for (j=0; j<2; j++) {
			for (i=0; i<ROWSM(e); i++)  ROWSID(e,i) = ISTRANSPM(e)? -i-1:i+1;
			ematrix_transpose(e);
		}
		x = (ROWSM(e)>COLSM(e)? ROWSM(e):COLSM(e));
	}
	return x;
}


/**
 * This creates an identical copy (by memory) of the given matrix - returns a pointer to it.
 * It works both for stand-alone and refering matrices, and keeps the transpose state.
**/

ematrix*	ematrix_copy(ematrix *e) {
	ematrix *e2;
	long	df;
	
	if (e->dlength<=0 || e->dlength>EM_MAXMXSIZE || e->transp<0 || e->isref<0)
		{PROGERROREXIT("copying a matrix %p that is in an incorrect state",e);}
	
	if (ISREFMAT(e)) {
		e2 = ematrix_copy_to(e,ROWSMAX(e),COLSMAX(e));
	
	} else if ((e2 = ematrix_alloc(e->dlength,0))!=NULL) {
		ematrix_memcpy(e2,e,e->dlength);
		e2->refm = NULL;  e2->isref = e2->numrefm = 0;
		df = (byte1*)e2-(byte1*)e;		/* the byte difference between e2 and e */
		e2->exp = (exp_t*)((byte1*)e->exp+df);
		e2->sgn = (sign_t*)((byte1*)e->sgn+df);
		e2->lid[0] = (short*)((byte1*)e->lid[0]+df);
		e2->lid[1] = (short*)((byte1*)e->lid[1]+df);
		DEBUG(CURDLEV+2,"copy matrix %p created from %p\n",e2,e);
	}
	return e2;
}

/**
 * This creates a different-sized copy (by entries) of the given matrix - the copy has exactly
 * the same entries and the actual size and transposition state,
 * but it has the given maximal dimensions (not from e).
 * The new matrix is returned as a pointer.
 * If trn==1, then the copy has the same transposition state as e, otherwise it has 0-transp.
 * It works both for stand-alone and refering matrices (but always creates a stand-alone matrix).
**/

ematrix*	ematrix_copy_ext(ematrix *e, int maxr, int maxc, int trn) {
	ematrix *e2;
	int	i,j;
	
	if (e->dlength<=0 || e->dlength>EM_MAXMXSIZE || e->transp<0 || e->isref<0)
		{PROGERROREXIT("copying a matrix %p that is in an incorrect state",e);}
	
	if (trn && ISTRANSPM(e)) {
		e2 = new_ematrix(COLSM(e),ROWSM(e),maxc,maxr);
		if (e2!=NULL)  ematrix_transpose(e2);
	} else {
		e2 = new_ematrix(ROWSM(e),COLSM(e),maxr,maxc);
	}
	if (e2!=NULL) {
		for (i=0; i<ROWSM(e2); i++) for (j=0; j<COLSM(e2); j++)
			COPYEXSIGM(e2,i,j,e,i,j);
		for (i=0; i<ROWSM(e2); i++)  ROWSID(e2,i) = ROWSID(e,i);
		for (i=0; i<COLSM(e2); i++)  COLSID(e2,i) = COLSID(e,i);
		DEBUG(CURDLEV+2,"copy matrix %p of maxdims %d,%d created from %p\n",e2,ROWSMAX(e2),COLSMAX(e2),e);
	}
	return e2;
}


/**
 * This creates a "reference" matrix - a matrix that refers to a submatrix of the matrix e.
 * In other words, the matrix created here has no own entris, it just points to
 * selected entries (submatrix) of e from rows ro1,ro1+1,...,ro2-1 and columns co1,...,co2-1.
 * All changes made to the entries of the new matrix appear in the matrix e.
 * It returns the newly created matrix (having the same transposition state as e).
 * 
 * If it is called for a refering matrix e, then the new matrix actually refers to the same
 * matrix as e referred to (behaves more like a copy).
 * However, a reference to the matrix e is remembered as well to be used later in ematrix_refadd_rc().
**/

ematrix*	ematrix_refer(ematrix *e, int ro1, int ro2, int co1, int co2) {
	ematrix *e2;
	int	i, sz,lsz;
	
	if ((co1<0 && co2!=co1) || co2>COLSM(e) || (ro1<0 && ro2!=ro1) || ro2>ROWSM(e))
		{PROGERROR("refering lines of %p out of range %d,%d %d,%d!",e,ro1,ro2,co1,co2);}
	lsz = MEMALIGNSIZE(ROWSMAX(e)*sizeof(e->refs[0][0]));
	lsz += MEMALIGNSIZE(COLSMAX(e)*sizeof(e->refs[0][0]));
	sz = MEMALIGNSIZE(sizeof(*e2))+lsz;
	e2 = ematrix_alloc(sz,0);
	if (e2!=NULL) {
		ematrix_memcpy(e2,e,sizeof(*e));	/* (the exp and sign references are copied here as well) */
		REFROW(e2) = (int*)MEMALIGNPOINT(e2+1);	/* reference arrays start right after the ematrix structure */
		REFCOL(e2) = (int*)MEMALIGNPOINT(REFROW(e2)+ROWSMAX(e));
		ROWSM(e2) = ro2-ro1;  COLSM(e2) = co2-co1;
		e2->numrefm = 0;  e2->isref = 1;
		e2->refm = ISREFMAT(e)? REFMAT(e):e;
		e2->refma = e;	/* (must remember the matrix e as well for later use in refadd!) */
		REFMAT(e2)->numrefm++;
		/*if (ISREFMAT(e) && ro1==0 && co1==0 && ro2>1 && co2>1) {
			MEMCPY(REFROW(e2),REFROW(e),lsz);
		} else ******** why this does not work??? */{
			for (i=ro1; i<ro2; i++)		/* creates the references for e2 to e */
				REFROW(e2)[i-ro1] = ISREFMAT(e)? REFROW(e)[i]:INDEXM(e,i,0);
			for (i=co1; i<co2; i++)
				REFCOL(e2)[i-co1] = ISREFMAT(e)? REFCOL(e)[i]:INDEXM(e,0,i);
		}
		DEBUG(CURDLEV+1,"new refering matrix %p created for %p: %d-%d,%d-%d\n",e2,e,ro1,ro2,co1,co2);
	}
	return e2;
}


/**
 * This changes the current reference to another matrix in a refering matrix e to point to
 * the given matrix eto.
 * It is supposed that both matrices have the same dimensions.
 * The function does not work well with "id tags" of rows and columns (x transpose)!!!
 * Use with care...
**/

void	ematrix_rerefer(ematrix *e, ematrix *eto) {
	
	if (!ISREFMAT(e) || ISREFMAT(eto))
		{PROGERROREXIT("rerefer works only for a refering matrix (%p) to a nonrefering matrix (%p)",e,eto);}
	if (ROWSM(REFMAT(e))!=ROWSM(eto) || COLSM(REFMAT(e))!=COLSM(eto) ||\
			ROWSMAX(REFMAT(e))!=ROWSMAX(eto) || COLSMAX(REFMAT(e))!=COLSMAX(eto))
		{PROGERROR("re-refering to a different size matrix (%p -> %p)",e,eto);}
	if (ISTRANSPM(REFMAT(e))!=ISTRANSPM(eto))
		{PROGERROR("re-refering to a different transpose-state matrix (%p -> %p)",e,eto);}
	DEBUG(CURDLEV+2,"re-refering matrix %p to a new matrix %p\n",e,eto);
	e->refm->numrefm--;
	e->refm = eto;
	e->refm->numrefm++;
	e->exp = eto->exp;  e->sgn = eto->sgn;
	ROWSIDS(e) = ROWSIDS(eto);  COLSIDS(e) = COLSIDS(eto);
}


/**
 * This transposes the given matrix - simply reverses the meaning of rows and columns.
 * (The modified matrix is stored in the same e - must be so since the return value is usually not used.)
 * Transpose is now implemented by a separate macro !!!
**/

/*	ematrix*	ematrix_transpose_set(ematrix *e, int tr) {
	e->transp = tr?1:0;
	return e;
}*/
ematrix*	ematrix_copydual(ematrix *e) {
	ematrix *e2;
	e2 = ematrix_copy(e);
	ematrix_transpose(e2);
	return e2;
}


/**
 * Given row and column ro,co are removed from the matrix e.
 * (Rows and columns are indexed from 0.)
 * If ro,co are not valid indices, they are not removed (no message printed).
 * If e is a refering matrix, then the references in the lines are forgotten.
 * The line id's in the matrix are moved accordingly.
 * (The modified matrix is stored in the same e.)
**/

void	ematrix_remove(ematrix *e, int ro, int co) {
	int	i,j, i2,j2, r,c;
	
	if (e->numrefm>0)  {PROGERROR("removing lines from a matrix %p that is referred from elsewhere!",e);}
	if (co>=COLSM(e) || ro>=ROWSM(e))  {PROGERROR("removing lines from %p out of range %d,%d!",e,ro,co);}
	
	if (co<0 || co>=COLSM(e))  c = co = COLSM(e);
	else  c = COLSM(e)-1;
	if (ro<0 || ro>=ROWSM(e))  r = ro = ROWSM(e);
	else  r = ROWSM(e)-1;
	if (co<c || ro<r) {
		if (!ISREFMAT(e)) {
			for (i=0; i<r; i++)  for (j=0; j<c; j++) {
				if (i<ro && j<co)  continue;
				i2 = i+(i>=ro);  j2 = j+(j>=co);
				COPYEXSIGM(e,i,j,e,i2,j2);
			}
			for (i=ro; i<r; i++)  ROWSID(e,i) = ROWSID(e,i+1);
			for (j=co; j<c; j++)  COLSID(e,j) = COLSID(e,j+1);
		} else {
			for (i=ro; i<r; i++)  REFROW(e)[i] = REFROW(e)[i+1];
			for (j=co; j<c; j++)  REFCOL(e)[j] = REFCOL(e)[j+1];
	}	}
	COLSM(e) = c;  ROWSM(e) = r;
}


/**
 * Given row and column ro,co of the referred matrix REFMAT(e) are added as the last lines
 * to the refering matrix e (provided dr==0).
 * However, if e was created as refering to another refering matrix er, then the lines ro,co
 * are actually taken from their positions in e->refma provided that dr>0
 * (but finally pointing to e, of course).
 * This is done to make refering a completely transparent process unless dr==0 is given.
 * Special care is taken for a refered matrix of a different transpose state.
 * The row or column is not added if it is negative (-1).
 * (Rows and columns are indexed from 0.)
 * The modified matrix is stored in the same e.
**/

void	ematrix_refadd_ext(ematrix *e, int dr, int ro, int co) {
	int	tr, i;
	ematrix	*er;
	
	if (!ISREFMAT(e))  {PROGERROREXIT("refadd works only for a reference matrix (%p)",e);}
	er = (dr>0? e->refma:REFMAT(e));
		/* e->refma is the matrix remembered in ematrix_refer() when creating e - defined only there */
	tr = ISTRANSPM(er);
	ematrix_transpose_set(er,ISTRANSPM(e));	/* must be transposed same for refering to work */
	
	if (ro>=ROWSM(er) || co>=COLSM(er))
		{PROGERROREXIT("refadd (%p) out of dest range %d,%d >= %d,%d",e,ro,co,ROWSM(er),COLSM(er));}
	if (ro>=0 && ro<ROWSM(er)) {
		if (++ROWSM(e)>ROWSMAX(e))
			{PROGERROREXIT("got out of matrix size when adding %d > %d",ROWSM(e),ROWSMAX(e));}
		else
			REFROW(e)[ROWSM(e)-1] = ISREFMAT(er)? REFROW(er)[ro]:INDEXM(er,ro,0);
	}				/* must use the same full indexing as in ematrix_refer() */
	if (co>=0 && co<COLSM(er)) {
		if (++COLSM(e)>COLSMAX(e))
			{PROGERROREXIT("got out of matrix size when adding %d > %d",COLSM(e),COLSMAX(e));}
		else
			REFCOL(e)[COLSM(e)-1] = ISREFMAT(er)? REFCOL(er)[co]:INDEXM(er,0,co);
	}
#ifndef FASTPROG
	if (IFRANDDEBUGLESS(111)) {
		for (i=0; i<ROWSM(e)-1; i++)
			if (REFROW(e)[i]==REFROW(e)[ROWSM(e)-1])  {PROGERROR("duplicated row reference %d in %p (%dx%d) -> (%dx%d)",GETREFMROW(e,i),e,ROWSM(e),COLSM(e),ROWSM(er),COLSM(er)); EMATDEBUG(0,e," !\t");}
		for (i=0; i<COLSM(e)-1; i++)
			if (REFCOL(e)[i]==REFCOL(e)[COLSM(e)-1])  {PROGERROR("duplicated col reference %d in %p (%dx%d) -> (%dx%d)",GETREFMCOL(e,i),e,ROWSM(e),COLSM(e),ROWSM(er),COLSM(er));}
	}
#endif
	ematrix_transpose_set(er,i=tr);
}


/**
 * This function appends a new column to the matrix e (as the last one),
 * and fills the new column with all zeros.
 * If the maximum dimension of e is reached, then a new larger copy of e is created.
 * The new line id is given according to the transposition state and used id's.
 * The modified matrix is returned (usually just stored in the same e).
**/

ematrix*	ematrix_append_rc(ematrix *e, int tr) {
	int		i,j,k;
	ematrix		*e2;
	
	if (ISREFMAT(e))  {PROGERROR("cannot append row/column to a reference matrix %p",e); return e;}
	if (tr)  ematrix_transpose(e);
	if (ROWSM(e)>=ROWSMAX(e)) {
		e2 = ematrix_copy_to(e,ROWSMAX(e)+1,COLSMAX(e));
		dispose_ematrix(e);  e = e2;
		DEBUG(CURDLEV-2,"enlarging the matrix %p and disposing the old matrix %p !\n",e2,e);
	}
	ROWSM(e)++;
	for (i=0; i<COLSM(e); i++)  SETEXSIGMZERO(e,ROWSM(e)-1,i);
	
	for (i=k=0; i<ROWSM(e)+COLSM(e); i++) {
		j = (i<ROWSM(e)? ROWSID(e,i):COLSID(e,i-ROWSM(e)));
		j = (ISTRANSPM(e)? -j:j);
		if (k<j)  k = j;	/* (using the next available line id +-k) */
	}
	ROWSID(e,ROWSM(e)-1) = ISTRANSPM(e)? -k-1:k+1;
	if (tr)  ematrix_transpose(e);
	return e;
}


/**
 * These functions swap two rows/columns of the matrix e.
 * (The modified matrix is stored in the same e.)
**/

void	ematrix_swap_rc(ematrix *e, int tr, int s1, int s2) {
	int	i;
	exp_t	xx;
	sign_t	gg;
	
	if (s1==s2)  return;
	if (!tr)  ematrix_transpose(e);
	if (!ISREFMAT(e)) {
		for (i=0; i<ROWSM(e); i++) {
			xx = EXPM(e,i,s1);  EXPM(e,i,s1) = EXPM(e,i,s2);  EXPM(e,i,s2) = xx;
			gg = SIGNM(e,i,s1);  SIGNM(e,i,s1) = SIGNM(e,i,s2);  SIGNM(e,i,s2) = gg;
		}
		i = COLSID(e,s1);  COLSID(e,s1) = COLSID(e,s2);  COLSID(e,s2) = i;
	} else {
		i = REFCOL(e)[s1];  REFCOL(e)[s1] = REFCOL(e)[s2];  REFCOL(e)[s2] = i;
	}
	if (!tr)  ematrix_transpose(e);
}


/**
 * This function extracts a submatrix of the matrix e determined by the refering matrix re.
 * For md==0, the referred submatrix is extracted, for md==1 the other rows (not of re) are
 * taken together with the columns of re, for md==2 similar happens for the other columns,
 * and for md==3 the rows and columns not in re are extracted.
 * The extracted matrix is allocated as refering to e (with the same max dimensions as e)
 * and returned (in the same transposition state as e).
**/

ematrix*	ematrix_refextract_ext(ematrix *e, ematrix *re, int md) {
	int	i,j, *ar,*arc, arstack[30];
	ematrix	*ex;
	
	if (ISREFMAT(e) || !ISREFMAT(re) || REFMAT(re)!=e)  {PROGERROREXIT("call only for a refering matrix re -> e");}
	
	if (ROWSM(e)+COLSM(e)<28)  ar = arstack;	/* (a faster stack allocation for small) */
	else  ar = MMALLOC((ROWSM(e)+COLSM(e)+2)*sizeof(ar[0]));
	for (i=0; i<ROWSM(e)+COLSM(e); i++)  ar[i] = 0;
	arc = ar+ROWSM(e);	/* mark what lines of e are present in re */
	for (i=0; i<ROWSM(re); i++)  ar[GETREFMROW(re,i)] = 1;
	for (i=0; i<COLSM(re); i++)  arc[GETREFMCOL(re,i)] = 1;
				/* construct the output refering matrix here */
	ex = ematrix_refer_empty(e);
	for (i=0; i<ROWSM(e); i++)
		if (ar[i]==((md&1)==0))  ematrix_refadd_row(ex,i);
	for (j=0; j<COLSM(e); j++)
		if (arc[j]==((md&2)==0))  ematrix_refadd_col(ex,j);
#ifndef FASTPROG
	if (ROWSM(ex)!=((md&1)?ROWSM(e)-ROWSM(re):ROWSM(re)) || COLSM(ex)!=((md&2)?COLSM(e)-COLSM(re):COLSM(re)))
		{PROGERROR("wrong extraction %d of a matrix, why ???",md); EMATDEBUG(0,e,"\t"); EMATDEBUGS(0,re,"\t <\t"); EMATDEBUGS(0,ex,"\t->\t");}
#endif
	if (ar && ar!=arstack)  FREE(ar);
	return ex;
}


/**
 * This function returns a "union" of the two refering matrices re1,re2.
 * The matrices re1,re2 must refer to a common matrix, and they may overlap in some lines.
 * The returned matrix refers to the same matrix as re1,re2, and it contains a union
 * of the lines of re1 and of re2 (no repeated elements).
**/

ematrix*	ematrix_union_ext(ematrix *re1, ematrix *re2, int dj) {
	ematrix		*e, *x, *ex;
	int 		i,j,k, r=0, *ar,*arc, arstack[30];
	
	if (!ISREFMAT(re1) || !ISREFMAT(re2) || REFMAT(re1)!=REFMAT(re2))  {PROGERROREXIT("call only for refering matrices re1,re2 to the same matrix");}
	
	e = REFMAT(re1);
	if (ROWSM(e)+COLSM(e)<28)  ar = arstack;	/* (a faster stack allocation for small) */
	else  ar = MMALLOC((ROWSM(e)+COLSM(e)+2)*sizeof(ar[0]));
	for (i=0; i<ROWSM(e)+COLSM(e); i++)  ar[i] = 0;
	arc = ar+ROWSM(e);	/* mark what lines of e are present in re1,re2 */
	for (x=re1,k=0; k<2; x=re2,k++) {
		for (i=0; i<ROWSM(x); i++)  ar[GETREFMROW(x,i)]++;
		for (j=0; j<COLSM(x); j++)  arc[GETREFMCOL(x,j)]++;
	}
	if (dj>0) {		/* asking only for disjoint submatrices */
		ex = NULL;
		for (i=r=0; i<ROWSM(e)+COLSM(e) && !r; i++)
			if (ar[i]>1)  r = 1;
	} else {		/* construct the output refering matrix here */
		ex = ematrix_refer_empty(e);
		for (i=0; i<ROWSM(e); i++)
			if (ar[i]>0)  ematrix_refadd_row(ex,i);
		for (j=0; j<COLSM(e); j++)
			if (arc[j]>0)  ematrix_refadd_col(ex,j);
	}
#ifndef FASTPROG
	if (dj==0 && IFRANDDEBUGLESS(33)) {
		x = ematrix_union_ext(ex,RANDOM()%2?re1:re2,-1);
		if (ROWSM(x)!=ROWSM(ex) || COLSM(x)!=COLSM(ex))  {PROGERROR("wrong union of matrices, why???");}
		dispose_ematrix(x); }
#endif
	if (ar && ar!=arstack)  FREE(ar);
	return (dj>0? (r?NULL:(void*)1): ex);
}


/**
 * This function compares the two given matrices, first for their dimensions,
 * then for entry-wise equality over the pfield.
 * (No matter whether the matrix is stand-alone or refering...)
 * Moreover, the second variant compares also the transposition states which
 * must be identical.
**/

int	ematrix_isequal(ematrix *e1, ematrix *e2) {
	int	i,j;
	
	if (!e1 || !e2)
		return (!e1 && !e2);
	if (ROWSM(e1)!=ROWSM(e2) || COLSM(e1)!=COLSM(e2))
		return 0;
	for (i=0; i<ROWSM(e2); i++)  for (j=0; j<COLSM(e2); j++)
		if (!ematrix_isequal_entry(e1,i,j,e2,i,j))
			return 0;	/* (macro ematrix_isequal_entry() defined in the header) */
	return 1;
}

int	ematrix_isequaltr(ematrix *e1, ematrix *e2) {
	if (e1 && e2) if (ISTRANSPM(e1)!=ISTRANSPM(e2))  return 0;
	return ematrix_isequal(e1,e2);
}
                                                











/******************	Matrix printing		*******************/
/******************************************************************/


#undef CURDLEV
#define CURDLEV		7


/**
 * This prints the given matrix to the file fout (on the screen),
 * possibly highlights specified element kr,kc, and prefixes each line with pref.
 * The function is not intended for any kind of formal output, just for a quick
 * view of the matrix on screen (like debug messages).
 * For kr==-9 and kc==-9, no frame-lines are printed around.
**/

#define	EMPRFR	"--------------------------------------------------------------\n"

void	ematrix_print_ext(FILE *fout, ematrix *e, int kr, int kc, char *pref) {
	int	i,j;
	
	if (!e) {
		PROGERROR("not printing NULL matrix");  return;
	}
	fprintf(fout,"%s%s%s"
		"matrix %p [%s], r=%d, c=%d, tr=%d, ref=%p\n%s ",
		(kr==-9&&kc==-9)?"":pref,(kr==-9&&kc==-9)?"":EMPRFR,
		pref,e,EMNAME(e)?EMNAME(e):"",ROWSM(e),COLSM(e),ISTRANSPM(e),
		ISREFMAT(e)?REFMAT(e):NULL,pref);
	if (kr>-9 && ISREFMAT(e)) {
		for (i=0; i<ROWSM(e); i++)
			fprintf(fout," r%d=%d(%d)",i,GETREFMROW(e,i),REFROW(e)[i]);
		if (ROWSM(e)>0)  fprintf(fout,"\n%s ",pref);
		for (i=0; i<COLSM(e); i++)
			fprintf(fout," c%d=%d(%d)",i,GETREFMCOL(e,i),REFCOL(e)[i]);
		if (COLSM(e)>0)  fprintf(fout,"\n%s ",pref);
	}
	for (j=0; j<COLSM(e); j++)
		fprintf(fout,"\t '%d')",COLSID(e,j));
	if (!(kr==-9&&kc==-9))  fprintf(fout,"\n%s ",pref);
	
	for (i=0; i<ROWSM(e); i++) {
		if (COLSM(e)>0)  fprintf(fout,"\n%s ",pref);
		for (j=0; j<COLSM(e); j++) {
			if (j==0)  fprintf(fout,"'%d')\t",ROWSID(e,i));
			if (i==kr && j==kc)  fprintf(fout,"* ");
			if (!(kr==-9&&kc==-9))  fprintf(fout,"%5s\t",pfield_printvalue(16,EXPM(e,i,j),SIGNM(e,i,j)));
			else  fprintf(fout,"%5s\t",pfield_pvalue(7,EXPM(e,i,j),SIGNM(e,i,j)));
		}
	}
	fprintf(fout,"\n%s%s",(kr==-9&&kc==-9)?"":pref,(kr==-9&&kc==-9)?"":EMPRFR);
}





















































