
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
 * A matrix is a structure holding an array (considered 2-dimensional) of pfield elements
 * (see pfield.h for the format of the pfield numbers and the arithmetic),
 * together with its dimensions and maximal dimensions.
 * Specific functions are provided for allocating and freeing such matrices.
 * Macros EXPM, SIGNM, ROWSM, COLSM should be used to access elements of a matrix.
 * It is easy to switch between reading the matrix by rows or by columns ("transpose").
 * There are "id tags" (labels) attached to each row and column that follow them during matrix
 * operations (swap, pivot, etc.) - rows have initial tags 1,2,..., columns -1,-2,...
 * 
 * Additionally, a matrix may just refer to a submatrix of another existing matrix.
 * In such case it does not have its own entries, but reference arrays for rows/columns.
 * The id tags are referred as well, respecting the transpose state.
 * Nothing outside ematrix.c should manipulate refering matrices directly (!).
 * (This gets liitle tricky when you try to refer an already refering matrix...)
 * 
 * Practical implementation of a matrix is in one piece of memory that holds the ematrix
 * structure at the beginning and referred data after that.
 * The field dlength keeps the allocated length, it is also cleared before freeing so
 * that access functions can recognize access to an invalid matrix.
 * The allocated matrices are not really freed, but recycled for new matrices of the same
 * size to speed-up the allocation process.
 * 
**/



#ifndef EMATRIX_H 
#define EMATRIX_H 





/******************	The matrix structure	*******************/
/******************************************************************/


	/**
	 * This is the matrix structure.
	 * Remember that, in addition to the structure members, the matrix uses variable-
	 * length arrays (exp,sign,lid, or refs) that are stored after the structure in memory.
	 * Because of this complex structure, special functions MUST be used to create,
	 * copy or dispose an ematrix.
	 * 
	 * How a stand-alone (i.e. not refering) matrix looks like:
	 * The entries are stored in a separate exponent and sign arrays (exp,sgn).
	 * The meaning of row/column for matrix access depends on the value of transp.
	 * (Also concerns max dimensions and ids.)
	 * 
	 * How a refering matrix looks like:
	 * There are no entries, but instead there are reference indices in refs[]
	 * that point to the entry arrays in the refered matrix.
	 * One may easily elnlarge refering matrix by adding other refered lines.
	 * Moreover, the refered matrix remembers how many matrices refer to it,
	 * for checking correct deallocation.
	 * 
	 * An optional name can be set for a matrix for easier identification in debugging.
	 * Each line of a matrix has its "id" label - a positive or negative number.
	 * These ids move with the matrix line when swapping or pivoting, etc.
	 * The ids are used in some functions to keep track of a line through a
	 * complex computation.
	 * They should be kept distinct for all lines, but they may be changed in general.
	 * Be careful when using id labels during computation (and check them!).
	 * 
	**/

typedef struct struct_ematrix {
	
	int	maxdims[2];	/* maximal dimensions for the current matrix */
	exp_t	*exp;		/* the arrays of exponents and signs of the matrix */
	sign_t	*sgn;
	int	dims[2];	/* current dimensions of the matrix (0-rows, 1-cols) */
	short	*lid[2];	/* "id tags" for rows and columns of the matrix */
	int	transp;		/* a 0/1 flag whether the matrix is currently transposed */
	
	int	isref;		/* a flag whether this matrix refers to entries of another matrix */
	int	*refs[2];	/* reference indices (0-rows, 1-cols) and the referred matrix */
	struct struct_ematrix	*refm;
	struct struct_ematrix	*refma;	/* "intermediate" ref matrix for use in refadd - do not use anywhere else */
	int	numrefm;	/* the number of matrices referring to this one */
	
	char	*name;		/* the name of this matrix, or NULL by default (for output use) */
	void	*data;		/* a place to refer arbitrary data (free()'d when disposing!) */
	long	dlength;	/* keeps this matrix data length, see also exp,sgn or rowref */
	
}	ematrix;

#define	EMMAGICNUM	111222
#define	EMMAGIC(m)	(*((int*)((byte1*)(m)+((ematrix*)(m))->dlength)))
#define	EMMAGICTEST(m)	(((m)?EMPREMAGIC((ematrix*)(m)):0)? EMMAGIC(m)==EMMAGICNUM:0)
#define	EMPREMAGIC(m)	((m)->maxdims[0]>=0&&(m)->maxdims[1]>=0 && (m)->maxdims[0]<10000&&(m)->maxdims[1]<10000 && \
			 (m)->dlength>0 && (m)->transp>=0&&(m)->transp<100)


	/**
	 * These are the macros provided to access the numbers of rows and columns,
	 * the exponents and signs of entries, and other matrix attributes.
	 * 
	 * There are two indexing schemes for access to matrix entries --
	 * one (INDEXM) applies to a stand-alone matrix (i.e. having its own data),
	 * the other one (INDEXR) applies to a matrix which refers to another submatrix.
	 * In the first scheme we compute usual 2-dim array index.
	 * In the second scheme, we remember the offsets (in the original matrix)
	 * of the rows and columns that are referred from our refering matrix, and
	 * then we simply sum the offsets to get an entry position.
	 * The macros GETREFMROW/COL provide "inverse to indexing", i.e. they return
	 * the row/column index in the original matrix for our refering matrix entry.
	 * (Note that, during use of GETREFMROW/COL, the original matrix must not be transposed!)
	 * 
	**/

#define TRANSPM(m)	((m)->transp)
#define ISTRANSPM(m)	((m)->transp?1:0)
#define ISREFMAT(m)	((m)->isref?1:0)

#define	COLSM(m)	(m)->dims[!(m)->transp]
#define	COLSMAX(m)	(m)->maxdims[!(m)->transp]
#ifdef		FASTPROG
#define	ROWSM(m)	(m)->dims[(m)->transp]
#define	ROWSMAX(m)	(m)->maxdims[(m)->transp]
#define	REFMAT(m)	((m)->refm)
#define	REFROW(m)	(m)->refs[(m)->transp]
#else
#define	ROWSM(m)	(m)->dims[!!(m)->transp]
#define	ROWSMAX(m)	(m)->maxdims[!!(m)->transp]
#define	REFMAT(m)	((m)->isref?0:PROGERROR("this matrix %p is not refering!",m),(m)->refm)
#define	REFROW(m)	(m)->refs[!!(m)->transp]
#endif
#define	REFCOL(m)	(m)->refs[!(m)->transp]

#define INDEXM(m,r,c)	((m)->transp?(m)->maxdims[1]*(c)+(r):(m)->maxdims[1]*(r)+(c))
#define INDEXR(m,r,c)	(REFROW(m)[r]+REFCOL(m)[c])

#define	GETLINIX(m,x)	((x)<(m)->maxdims[1]? (x):(x)/(m)->maxdims[1])
#define	GETREFMROW(m,r)	GETLINIX((m)->refm,REFROW(m)[r])
#define	GETREFMCOL(m,c)	GETLINIX((m)->refm,REFCOL(m)[c])
		/* GETREFMROW/COL(m,i) are "inverse functions" to INDEXM(e,i,0)/(e,0,i) */

#define	ROWSIDS(m)	(m)->lid[!!(m)->transp]
#define	ROWSID(m,r)	ROWSIDS(m)[ISREFMAT(m)? GETREFMROW(m,r):(r)]
#define	COLSIDS(m)	(m)->lid[!(m)->transp]
#define	COLSID(m,c)	COLSIDS(m)[ISREFMAT(m)? GETREFMCOL(m,c):(c)]
#define	NOLINEID	-999999	/* (supposed not to be a line id - out of short range anyway) */

#ifdef		FASTPROG
#define	EXPM(m,r,c)	(m)->exp[ISREFMAT(m)? INDEXR(m,r,c):INDEXM(m,r,c)]
#define	SIGNM(m,r,c)	(m)->sgn[ISREFMAT(m)? INDEXR(m,r,c):INDEXM(m,r,c)]
#else		/* (the only difference between these two is in index-bound checking) */
#define	EXPM(m,r,c)	(m)->exp[ematrix_function_indexm(m,r,c)]
#define	SIGNM(m,r,c)	(m)->sgn[ematrix_function_indexm(m,r,c)]
#endif

#define	SETEXSIGM(m,r,c,x,g)	(EXPM(m,r,c)=(x),SIGNM(m,r,c)=(g))
#define	SETEXSIGMZERO(m,r,c)	(pfield_setzeroexp(&(EXPM(m,r,c))),SIGNM(m,r,c)=0)
#define	COPYEXSIGM(m1,r1,c1,m2,r2,c2)	(EXPM(m1,r1,c1)=EXPM(m2,r2,c2),SIGNM(m1,r1,c1)=SIGNM(m2,r2,c2))

#define	EMSETNAME(m,s)	(m)->name = MSTRDUP(s)
#define	EMNAME(m)	(m)->name
#define	EMDATA(m)	(m)->data

extern exp_t	macexp;
int     ematrix_function_indexm(ematrix *m, int r, int c) ;	/* like macro indexing, but with checking */







/******************	Functions from  ematrix.c	*******************/
/**************************************************************************/


/**
 * These are the basic manipulating functions for matrices.
 * Their meaning is pretty obvious... (see more in ematrix.c)
 * ....
 * 
 * One feature of ematrix is an easy fast transposition operation (which just
 * changes an internal flag in ematrix).
 * However, keep in mind that some operations (in particular those involving
 * refering matrices) require you to keep a particular transposition state
 * of the matrix.
 * Also some additional information about the matrix may be dependent on
 * the particular transposition state (or better on the 0-transposition state).
 * 
 * Read carefuly the documentation about "refering matrices" (those that
 * do not have entries, but refer to a submatrix of another matrix instead)
 * before playing tricks with them!!! - ematrix_refer..(), ematrix_refadd..()
 * In general, refering a matrix is a transparent action - you may work with
 * the refering matrix in the same way as with a stand-alone matrix, including
 * further refering and copies.
 * Of course, all entry-changes made to the refering matrix appear in the referred matrix.
 * (For example, pivoting would be quite damaging to the referred matrix...)
 * But, the actions stay transparent only as long as you do not use functions
 * like ematrix_drefadd..() or ematrix_rerefer() !
 * 
**/


ematrix*        new_ematrix(int ro, int co, int mr, int mc) ;
void		dispose_ematrix(ematrix *e) ;

int     ematrix_lineids(int ch, ematrix *e) ;
#define	ematrix_resetid(e)	ematrix_lineids(1,e)
#define	ematrix_checkid(e)	ematrix_lineids(0,e)

#define	ematrix_transpose_set(e,tr)	(TRANSPM(e)=(tr))
#define	ematrix_transpose(e)		ematrix_transpose_set(e,!ISTRANSPM(e))
ematrix*        ematrix_copydual(ematrix *e) ;

ematrix*        ematrix_copy(ematrix *e) ;
ematrix*        ematrix_copy_ext(ematrix *e, int maxr, int maxc, int trn) ;
#define		ematrix_copy_to(e,maxr,maxc)	ematrix_copy_ext(e,maxr,maxc,1)
#define		ematrix_copy_notr(e)		ematrix_copy_ext(e,ROWSMAX(e),COLSMAX(e),0)

ematrix*        ematrix_refer(ematrix *e, int ro1, int ro2, int co1, int co2) ;
#define		ematrix_refer_all(e)	ematrix_refer(e,0,ROWSM(e),0,COLSM(e))
#define		ematrix_refer_empty(e)	ematrix_refer(e,-1,-1,-1,-1)
void    	ematrix_rerefer(ematrix *e, ematrix *eto) ;

void    ematrix_remove(ematrix *e, int ro, int co) ;
#define	ematrix_remove_rc(e,tr,so)	ematrix_remove(e,((tr)?-1:(so)),((tr)?(so):-1))
#define	ematrix_remove_row(e,ro)	ematrix_remove(e,ro,-1)
#define	ematrix_remove_col(e,co)	ematrix_remove(e,-1,co)

void    ematrix_refadd_ext(ematrix *e, int dr, int ro, int co) ;
#define ematrix_refadd_rc(e,ro,co)	ematrix_refadd_ext(e,1,ro,co)
#define	ematrix_refadd_row(e,ro)	ematrix_refadd_rc(e,ro,-1)
#define	ematrix_refadd_col(e,co)	ematrix_refadd_rc(e,-1,co)
	/* (the drefadd_... versions specifically request lines from the final refered matrix) */
#define ematrix_drefadd_rc(e,ro,co)	ematrix_refadd_ext(e,0,ro,co)
#define	ematrix_drefadd_row(e,ro)	ematrix_drefadd_rc(e,ro,-1)
#define	ematrix_drefadd_col(e,co)	ematrix_drefadd_rc(e,-1,co)

ematrix*        ematrix_append_rc(ematrix *e, int tr) ;
#define	ematrix_append_row(e)	ematrix_append_rc(e,0)
#define	ematrix_append_col(e)	ematrix_append_rc(e,1)

void    ematrix_swap_rc(ematrix *e, int tr, int s1, int s2) ;
#define	ematrix_swaprows(e,r1,r2)	ematrix_swap_rc(e,0,r1,r2)
#define	ematrix_swapcols(e,c1,c2)	ematrix_swap_rc(e,1,c1,c2)

ematrix*        ematrix_refextract_ext(ematrix *e, ematrix *re, int md) ;
#define	ematrix_refextract_xrow(e,re)	ematrix_refextract_ext(e,re,1)
#define	ematrix_refextract_xcol(e,re)	ematrix_refextract_ext(e,re,2)
#define	ematrix_refextract_xall(e,re)	ematrix_refextract_ext(e,re,3)

ematrix*        ematrix_union_ext(ematrix *re1, ematrix *re2, int dj) ;
#define ematrix_union(e1,e2)		ematrix_union_ext(e1,e2,0)
#define ematrix_isdisjoint(e1,e2)	(ematrix_union_ext(e1,e2,1)!=NULL)

#define	ematrix_isequal_entry(e1,r1,c1,e2,r2,c2)	\
		pfield_isequal(EXPM(e1,r1,c1),SIGNM(e1,r1,c1),EXPM(e2,r2,c2),SIGNM(e2,r2,c2))
int     ematrix_isequal(ematrix *e1, ematrix *e2) ;
int     ematrix_isequaltr(ematrix *e1, ematrix *e2) ;

void    ematrix_print_ext(FILE *fout, ematrix *e, int kr, int kc, char *pref) ;
#define	ematrix_print_mark(e,kr,kc)	ematrix_print_ext(stdout,e,kr,kc,"")
#define	ematrix_print_pref(e,pr)	ematrix_print_ext(stdout,e,-1,-1,pr)
#define	ematrix_print(e)		ematrix_print_pref(e,"")
#define	ematrix_print_nofr(e,pr)	ematrix_print_ext(stdout,e,-9,-9,pr)

#define	ematrix_fprint_pref(fo,e,pr)	ematrix_print_ext(fo,e,-1,-1,pr)
#define	ematrix_fprint(fo,e)		ematrix_fprint_pref(fo,e,"")
#define	ematrix_fprint_nofr(fo,e,pr)	ematrix_print_ext(fo,e,-9,-9,pr)

#define EMATOUTPUT(e,pf)	ematrix_fprint_pref(printout,e,pf)
#define EMATOUTPUTS(e,pf)	ematrix_fprint_nofr(printout,e,pf)
#define	EMATDEBUG(l,e,pf)	(junk = TESTDLEV(l)? (ematrix_fprint_pref(debugout,e,pf),0):0)
#define	EMATDEBUGS(l,e,pf)	(junk = TESTDLEV(l)? (ematrix_fprint_nofr(debugout,e,pf),0):0)







/******************	Functions from  ematrixop.c	*******************/
/**************************************************************************/


/**
 * These functions implement basic matrix operations - pivoting, scaling,
 * removing a line, computing determinant, and importing from another pfield.
 * Read their source descriptions...
**/


int     ematrix_pivot_ext(int ch, ematrix *e, int ro, int co) ;
#define	ematrix_pivot(e,ro,co)		ematrix_pivot_ext(1,e,ro,co)
	/* (validity of computation: ch>=1 errors are printed, ch==0 no errors but returns -1 for invalid) */
#define	ematrix_pivot_check(ch,e,ro,co)	ematrix_pivot_ext(ch,e,ro,co)

void    ematrix_multiply_rc(ematrix *e, int r1, int r2, int c1, int c2, int sx, exp_t x, sign_t g) ;
#define ematrix_multiply_row(e,r,s,x,g)	ematrix_multiply_rc(e,r,(r)+1,0,COLSM(e),s,x,g)
#define ematrix_multiply_col(e,c,s,x,g)	ematrix_multiply_rc(e,0,ROWSM(e),c,(c)+1,s,x,g)
#define ematrix_multiply_entry(e,r,c,s,x,g)	ematrix_multiply_rc(e,r,(r)+1,c,(c)+1,s,x,g)

ematrix*        ematrix_removemat_ext(ematrix *em, int dc, int tr, int so) ;
#define	ematrix_removemat(em,dc,tr,so)		ematrix_removemat_ext(em,dc,tr,so)
#define	ematrix_removemat_rcsum(em,dc,sum)	ematrix_removemat_ext(em,dc,\
						 (sum)>=ROWSM(em),(sum)<ROWSM(em)?(sum):(sum)-ROWSM(em))
#define	ematrix_removemat_rowdel(em,ro)		ematrix_removemat_ext(em,1,0,ro)
#define	ematrix_removemat_colcon(em,co)		ematrix_removemat_ext(em,0,1,co)
ematrix*        ematrix_removematid_ext(ematrix *em, int dc, int id) ;
#define	ematrix_removematid_dc(em,dc,id)	ematrix_removematid_ext(em,dc,id)
#define	ematrix_removematid_contr(em,id)	ematrix_removematid_dc(em,0,id)
#define	ematrix_removematid_del(em,id)		ematrix_removematid_dc(em,1,id)

int     ematrix_determinant_ch(int ch, ematrix *e, exp_t *x, sign_t *g) ;
	/* (validity of computation: ch>=1 errors are printed, ch==0 no errors but returns -1 for invalid) */
#define ematrix_determinant(e,x,g)	ematrix_determinant_ch(1,e,x,g)
#define ematrix_determinant_check(e)	ematrix_determinant_ch(0,e,NULL,NULL)

int     ematrix_determinant_2x2ext(int ch, ematrix *e, int r1, int r2, int c1, int c2, exp_t *x, sign_t *g) ;
#define	ematrix_determinant_2x2ch(ch,e,r1,r2,c1,c2,x,g)	ematrix_determinant_2x2ext(ch,e,r1,r2,c1,c2,x,g)
#define	ematrix_determinant_2x2check(e,r1,r2,c1,c2)	ematrix_determinant_2x2ch(0,e,r1,r2,c1,c2,NULL,NULL)
#define	ematrix_determinant_2x2(e,r1,r2,c1,c2,x,g)	ematrix_determinant_2x2ext(1,e,r1,r2,c1,c2,x,g)

int     ematrix_import_ext(int ch, ematrix *ee) ;
#define	ematrix_import(ee)	ematrix_import_ext(1,ee)

int     ematrix_pfendomorph_ext(int ch, int le, ematrix *ee) ;
#define	ematrix_pfendomorph(le,ee)	ematrix_pfendomorph_ext(1,le,ee)



/**
 * The following functions are about linear dependencies in the matrix.
 * They test two parallel lines (ignoring up to two entries) or zero lines.
 * The return value is either -1 for not parallel, or r>=0 if parallel??.
 * Next ematrix_lineparallel_..() look for other lines parallel to the given line
 * (-1 if no such found),
 * and ematrix_parallel_any() tests whether there are any two parallel lines
 * (value -1 if not, 0 if two rows are parallel, or r=1 if two columns, or r=2 if both).
 * In a parallel test when one of the lines is 0, the answer is parallel as well.
**/

int     ematrix_twoparallel_rc(ematrix *e1, int tr, int s1, ematrix *e2, int s2, int no1, int no2) ;

int     ematrix_linezero_ext(ematrix *e, int tr, int so, int no1, int no2, int no3, int nz, int ct) ;
#define	ematrix_linezero_rc(e,tr,so,no1,no2)	ematrix_linezero_ext(e,tr,so,no1,no2,-1,0,0)
#define	ematrix_linezero(e,tr,so)		ematrix_linezero_rc(e,tr,so,-1,-1)
#define	ematrix_lineunit(e,tr,so)		ematrix_linezero_ext(e,tr,so,-1,-1,-1,0,1)
	/* (the skipped entries must be nonzero if nz==1 below!) */
#define	ematrix_linezero_col3(e,co,n1,n2,n3)	ematrix_linezero_ext(e,1,co,n1,n2,n3,1,0)

int     ematrix_lineparallel_ext(int ch, ematrix *e, int tr, ematrix *eo, int so) ;
#define ematrix_lineparallel(e,tr,eo,so)	ematrix_lineparallel_ext(0,e,tr,eo,so)
#define ematrix_parallel_row(e,ro)		ematrix_lineparallel(e,0,e,ro)
#define ematrix_parallel_trow(e,et,ro)		ematrix_lineparallel(e,0,et,ro)
#define ematrix_parallel_col(e,co)		ematrix_lineparallel(e,1,e,co)
#define ematrix_parallel_tcol(e,et,co)		ematrix_lineparallel(e,1,et,co)

/**
 * The next two look for any zero or parallel lines in the matrix.
 * The return value is either -1 for no parallel, or r=0 if two rows are parallel,
 * or r=1 if two columns are parallel, or r=2 if both happen.
**/

int     ematrix_parallel_any(ematrix *e) ;
int     ematrix_linezero_any(ematrix *e) ;


/**
 * These macros test for parallel pairs in the matrix - either row and column,
 * or two columns.
 * The next are macros testing for triangles/triads of a particular shape in
 * the matrix, assuming that there are no parallel/serial pairs(!).
**/

#define	ematrix_isparallel_rc(e,r,c)	(ematrix_linezero_rc(e,1,c,r,-1)>=0)
#define	ematrix_isparallel_cc(e,c1,c2)	(ematrix_twoparallel_rc(e,1,c1,e,c2,-1,-1)>=0)

#define	ematrix_istriangle_rcc(e,r,c1,c2)	(ematrix_twoparallel_rc(e,1,c1,e,c2,r,-1)>=0)
#define	ematrix_istriangle_rrc(e,r1,r2,c)	(ematrix_linezero_rc(e,1,c,r1,r2)>=0)
#define	ematrix_istriangle_ccc(e,c1,c2,c3)	(ematrix_triax_rc(e,1,c1,c2,c3)>=0)

#define	ematrix_istriad_rrc(e,r1,r2,c)		(ematrix_twoparallel_rc(e,0,r1,e,r2,c,-1)>=0)
#define	ematrix_istriad_rcc(e,r,c1,c2)		(ematrix_linezero_rc(e,0,r,c1,c2)>=0)
#define	ematrix_istriad_rrr(e,r1,r2,r3)		(ematrix_triax_rc(e,0,r1,r2,r3)>=0)



/**
 * The function ematrix_matrank() returns the rank of the given matrix,
 * and ematrix_triax_rc() tests whether the three lines are dependent (-1 if not).
 * (Both can work for non-pfield matrices, but the result is strange.)
 * If lrk>=0 is given, then only rank up to lrk is computed, even if it is higher.
 * The function ematrix_setrank() returns the matroid-set rank of the element subset
 * given by lines of refering re in the matroid given by e.
 * 
 * The function ematrix_closure..() computes the linear (co)closure of the given
 * submatrix in the whole matrix, and
 * ematrix_whatsep() returns the connectivity of the given submatrix re in e,
 * i.e. the separation value (s-1) for the separation (re,e-re).
 * We use the formula for the connectivity of (A,B) as lambda(A) = r(A)+r(B)-r(M).
 * 
**/

int     ematrix_matrank_ext(ematrix *e, int lrk) ;
#define	ematrix_matrank(e)		ematrix_matrank_ext(e,-1)
#define	ematrix_ismatrank_less(e,lr)	((ROWSM(e)<(lr)||COLSM(e)<(lr))? 1: \
					 (ematrix_matrank_ext(e,(lr))<(lr)))

int     ematrix_triax_rc(ematrix *e, int tr, int so1, int so2, int so3) ;

int     ematrix_setrank(ematrix *e, ematrix *re) ;

ematrix*        ematrix_closure_ext(ematrix *e, ematrix *re, int tr) ;
#define	ematrix_closure_tr(e,re,tr)	ematrix_closure_ext(e,re,tr)
#define	ematrix_closure(e,re)		ematrix_closure_ext(e,re,1)
#define	ematrix_coclosure(e,re)		ematrix_closure_ext(e,re,0)

int     ematrix_whatsep_ext(ematrix *e, ematrix *re, int lsp) ;
#define	ematrix_whatsep(e,re)		ematrix_whatsep_ext(e,re,-1)
#define	ematrix_whatsep_bound(e,re,bs)	ematrix_whatsep_ext(e,re,bs)








/******************	Functions from  ematexp.c	*******************/
/**************************************************************************/


/**
 * The following function prepares a list of all submatrices of the given matrix e
 * according to several given parameters.
 * If rnd>0, then only random partial list is generated:
 *  - all submatrices of size up to 3 (up to 2) are generated if rnd==1 (rnd==2);
 *  - the list gets sparser with higher values of rnd>2.
 * The particular application macros follow...
 * 
 * All bases are generated by the next function.
 * Only square submatrices are returned if piv==0, or pivoted new matrices if piv==1.
 * If hh[],hhb[] are given, then only those bases of e for which the element codes in hh[]
 * give the same collection (multiset) as all the codes in hhb[] are generated.
**/

ematrix**       ematrix_submatrices_ext(ematrix *e, int sq, int nz,
                              int rnd, int tot, int mrk, int mins, int maxs, int kr, int kc) ;

#define	ematrix_submatrices_inpf(e,rnd,min,kr,kc) \
		ematrix_submatrices_ext(e,1,1,rnd,-1,-1,min,-1,((kr)?ROWSM(e)-1:-1),((kc)?COLSM(e)-1:-1))
#define	ematrix_submatrices_bases(e,rnd)	ematrix_submatrices_ext(e,1,1,rnd,-1,-1,-1,-1,-1,-1)
#define	ematrix_submatrices_dets(e,rnd,kr,kc)	ematrix_submatrices_ext(e,1,0,rnd,-1,-1,1,-1,kr,kc)
#define	ematrix_submatrices_sub(e,s)	ematrix_submatrices_ext(e,0,0,0,s,-1,-1,-1,-1,-1)
#define	ematrix_submatrices_all(e)	ematrix_submatrices_ext(e,0,0,0,-1,-1,-1,-1,-1,-1)
#define	ematrix_submatrices_zero(e)	ematrix_submatrices_ext(e,0,0,0,-1,0,1,-1,-1,-1)
#define	ematrix_submatrices_conn(e,rk)	ematrix_submatrices_ext(e,0,0,0,-1,rk,-1,-1,-1,-1)

ematrix**       ematrix_getbases_ext(int ch, ematrix *e, int piv, int rnd, long hh[], long hhb[]) ;
#define	ematrix_getbases_sq(e)	ematrix_getbases_ext(0,e,0,0,NULL,NULL)
#define	ematrix_getbases(e)	ematrix_getbases_ext(0,e,1,0,NULL,NULL)
#define	ematrix_getbases_rand(e)	ematrix_getbases_ext(0,e,1,3,NULL,NULL)
#define	ematrix_getbases_randl(e)	ematrix_getbases_ext(0,e,1,4,NULL,NULL)
	/* (this function rnd==4 returns only bases more than 2 lines from the current one) */
#define	ematrix_getbases_hh(e,hh,hhb)	ematrix_getbases_ext(0,e,1,0,hh,hhb)


/**
 * The next function tests whether the given matrix em is in the current pfield, i.e.
 * whether all of its subdeterminants (of any size) are defined in the pfield numbers.
 * (A matrix is properly represented over the pfield if all subdeterminants are
 * defined in the pfield. Then also all subsequent pivots are defined.)
 * Since the ematrix_inpfield() test is very slow in general, a randomized
 * (incomplete, but fast) test is provided as well.
 * The function returns 0 for a pfield matrix, and -1 otherwise.
**/

int     ematrix_inpfield_ext(int ch, ematrix *em, int klr, int klc) ;

#define ematrix_inpfield_check(ch,e)	ematrix_inpfield_ext(ch,e,0,0)
		/* (ch>=0 full test; ch>0 bad printed with error; ch==-7 randomized test) */
#define	ematrix_inpfield(e)		ematrix_inpfield_check(0,e)
#define	ematrix_inpfield_printed(e)	ematrix_inpfield_check(1,e)
#define	ematrix_inpfield_rand(e)	ematrix_inpfield_check(-7,e)
		/* (request to check only subdeterminants containing the last row or/and column) */
#define ematrix_inpfield_last(e,kr,kc)	ematrix_inpfield_ext(0,e,kr,kc)
#define	ematrix_inpfield_randlast(e,kr,kc)	ematrix_inpfield_ext(-7,e,kr,kc)


/**
 * This function tests whether the two given matrices define the same matroid.
 * (That means whether all subdeterminants have simultaneously zero/nonzero values.)
 * When xpf2>=0 is given, then the second matrix e2 is represented over another pfield xpf2.
 * The return value is 1 for identical dets, 0 for different size matrices, and -d
 * when the smallest distinct subdeterminat pair is of size d.
 * Parameter undf>0 prints error for undefined, undf==0 returns false for undefined,
 * and undf<0 treats undefined subdeterminants as nonzero.
**/

int     ematrix_samedets_ext(int ch, ematrix *e1, ematrix *e2, int xpf2, int undf, int klr, int klc) ;
#define	ematrix_havesamedets(e1,e2,xp2)	(ematrix_samedets_ext(0,e1,e2,xp2,1,-1,-1)>0)
#define	ematrix_havesamedets_repres(e1,e2,xp2,trq) \
		(ematrix_samedets_ext(0,e1,e2,xp2,0,((trq)?-1:ROWSM(e1)-1),((trq)?COLSM(e1)-1:-1))>0)


/**
 * This is a function for an extended structural printing of a matroid.
 * Many interesting matroid-invariant numbers are printed here, derived from the
 * list of all bases, from the closure operator (flats), separations, etc...
 * 
 * The function also computes a matroid hash-value from the basis-structure of the matrix,
 * which is matroid-invariant and may be used to informally compare two matroids.
 * When the hash-value changes, the version of it must be updated, as well as the
 * value EM_HASHR10 for the R10 matroid (for debug checking).
 * The hash-value should not be used anywhere inside the program routines!
 * (Better see element "magic-values" for a matroid defined in struct/strmagic.c ...)
**/


long    ematrix_printmore_ext(ematrix *ee, int lev, char *bto, int mx) ;
#define	ematrix_printmore(ee,lev)		ematrix_printmore_ext(ee,lev,NULL,0)
#define	ematrix_printmore_to(ee,lev,bf,mx)	ematrix_printmore_ext(ee,lev,bf,mx)

#define	EM_HASHVER	"1.0"
#define	EM_HASHR10	26051682l
#define	ematrix_getmhash(ee)		ematrix_printmore(ee,-1)

/**
 * This function is used to print out all matroid bases and all circuits.
 * The matroid is given in e, parameter whp determines the function - <=10 for printing bases,
 * >=10 for printing circuits, and lev is the verbosity printing level.
 * One may possibly specify elements (via their id's) that must be contained in the printed
 * bases/circuits, in the array cix[], where NULL or 0 value mean no element.
**/

void    ematrix_printbasecirc_ext(ematrix *e, int whp, int lev, int *cix, char *bto, int mx) ;
#define	ematrix_printbasecirc(e,w,l,cx)	ematrix_printbasecirc_ext(e,w,l,cx,NULL,0)
#define	ematrix_printbases(e,l,cx)	ematrix_printbasecirc_ext(e,1,l,cx,NULL,0)
#define	ematrix_printcircuits(e,l,cx)	ematrix_printbasecirc_ext(e,11,l,cx,NULL,0)


















#endif	/* (of #ifndef EMATRIX_H) */

























