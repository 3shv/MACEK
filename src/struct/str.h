
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
 * This part contains fuctions for "structural work" with represented matroids.
 * (We do not attempt to work with matroids without a matrix representation - too complicated.)
 * Read ematrix.h to learn about our matrices, and pfield.h about our pfield arithmetic.
 * Read ../include/struct.h for more information about these structural functions.
 * 
 * The functions here are only for internal use in other struct/ * functions...
**/



#ifndef STR_H 
#define STR_H 










/******************     Magic numbers (strmagic.c)	*******************/
/**************************************************************************/


/**
 * When working with matrices or with the underlying matroids, it is often useful to
 * consider some "magic numbers" associated with these objects.
 * The magic numbers are, by default, invariant on the operations that we consider
 * (that means scaling and permuting lines for matrices, or taking various representations
 * for an abstract matroid in the second case).
 * We may then apply various tricks using injections or symmetric functions to these numbers...
 * 
 * To preserve compatibility between functions, we store all magic numbers into the long type.
 * 
**/


/**
 * The first function decides whether there is an injection from the array fr[] to the array to[]
 * such that elements of value other than nomag (==0) are mapped to elements of the same value.
 * The return value is 0 for no injection, or d+1 where d is the number of distinct values.
 * 
 * Next we compute a simple symmetric function of the values in the array ar[] of length arn.
 * 
**/

int     strmag_isinjection_ext(long *fr, int frn, long *to, int ton, long nomag) ;
#define	strmag_isinjection(fr,frn,to,ton)	strmag_isinjection_ext(fr,frn,to,ton,0l)

long    strmag_symmfunction_ext(int lev, long *ar, int arn) ;
#define	strmag_symmfunction(ar,arn)	strmag_symmfunction_ext(2,ar,arn)




/**
 * This function computes the zero-pattern of the lines (rows tr==0, columns tr==1) of
 * the matrix e; storing the zeros of each line into bits of an element of zp[],
 * and the number of zeros +1 in the line to zn[].
 * The total number of zeros +1 is returned from the function.
 * 
 * Next we find the clusters of zero 2x2 subdeterminants of the two given rows r1,r2 of ee.
 * (The zero subdeterminants form cliques if viewed as edges of a graph.)
 * The same is computed for columns if tr==1.
 * The result is stored in ar[]; ar[i] contains the index of the first column j<i for
 * which the subdeterminant of columns i,j is zero, or ar[i]=i if there is no zero.  
 * The columns of both-0 entries are treated specially - they form their own cluster.
 * Subdeterminants may be undefined, then they are treated as nonzeros.
 * 
 * The next is a magic number computed from all entries of the given matrix.
 * 
**/

long    strmag_zerolines_ext(ematrix *e, int tr, int s1, int s2, long *zp, long *zn) ;
#define	strmag_zeropattern_cols(e,zp)	strmag_zerolines_ext(e,1,0,COLSM(e),zp,NULL)
#define	strmag_numzero_cols(e,nz)	strmag_zerolines_ext(e,1,0,COLSM(e),NULL,nz)
#define	strmag_numzero_rows(e,nz)	strmag_zerolines_ext(e,0,0,ROWSM(e),NULL,nz)
#define	strmag_numzeros(e)		strmag_zerolines_ext(e,1,0,COLSM(e),NULL,NULL)

long    strmag_twolinedet_ext(ematrix *e, int tr, int r1, int r2, int *ar) ;
#define	strmag_twolinedet_rows(e,r1,r2,ar)	strmag_twolinedet_ext(e,0,r1,r2,ar)
#define	strmag_twolinedet_numz(e,tr,r1,r2)	strmag_twolinedet_ext(e,tr,r1,r2,NULL)

long    strmag_matrixentries_ext(ematrix *e) ;
#define	strmag_matrixentries(e)	strmag_matrixentries_ext(e)




/**
 * This function computes magic numbers for matroid elements based on small flats in
 * the matroid.  These numbers are matroid-invariant!  
 * The numbers are returned in the array flt[] (indexed first by rows, then by columns).
**/
 
long    strmag_flatlines_ext(ematrix *e, int minr, int maxr, int tf, long flt[]) ;
#define	strmag_flatlines(e,mi,flt)	strmag_flatlines_ext(e,mi,(mi)+3,ROWSM(e)*COLSM(e)/3+1,flt)
#define	strmag_flatlines_more(e,mi,flt)	strmag_flatlines_ext(e,mi,(mi)+4,9*ROWSM(e)*COLSM(e),flt)
#define	strmag_flatnumber(e,mi)		strmag_flatlines_ext(e,mi,(mi)+3,2*ROWSM(e)*COLSM(e)+1,NULL)










/******************	Matrix line maps (strlmap.c)	*******************/
/**************************************************************************/


/**
 * Let us consider all injective line maps of the given matrix E to E'.
 * (that means injective maps from rows/columns of E to rows/columns of E').
 * To state that VERY clearly, an entry Eij at (i,j) is mapped by a line map to an entry
 * Eij*sr[pr[i]]*sc[pc[j]] at E'(pr[i],pc[j]), where pr,pc are the row and column maps
 * and sr,sc are the row and column scales (scales indexed by the destination matrix indices!).
 * Hence Eij*sr[pr[i]]*sc[pc[j]]==E'(pr[i],pc[j]) for all choices of i,j in E.
 * Lines of E' not covered by the map may be arbitrary, but sr[],sc[] is 0 on uncovered lines.
 * 
 * These maps are useful in describing matrix "automorphism", displayed minors,
 * and similar structural properties that are tied with a matroid representation
 * against a particular basis.
**/


/**
 * The following structure is supplied for description of a matrix line maps.
 * A matrix line map consists of arbitrary permutations of rows and columns and of
 * nonzero scaling of rows and columns (indexed by the destination matrix).
 * The associated row and column permutations are stored in prc[],
 * and the row and column scalings are stored in xrc[],grc[].
 * The interpretation of rows/columns depends on the current transposition state
 * of the source matrix, and the difference in transposition of the destination
 * matrix against the source one is stored in xtr.
 * The (maximal) dimensions of the source and destination matrices are in ns[], nd[].
 * 
 * Use the next macros for access to the structure fields, based on the source matrix
 * transposition.
 * Functions are provided for allocating a new line map structure (dipose it with "free()"),
 * and for checking and printing the structure.
 * There is a possibility of "recycling" the line maps, but this should be used with caution.
 * (The "recycling" definitions are in include/struct.h.)
**/

typedef struct {
	
	short	*prc[2];
	exp_t	*xrc[2];
	sign_t	*grc[2];
	int	ns[2], nd[2];
	int	xtr;
	int	sz;		/* (the size of this memory piece - for easy duplication) */
	
}	emlinemap;

#define	EMLMMAGICNUM	456321
#define	EMLMMAGIC(a)	(*((int*)((a)->sz>=0&&(a)->sz<9999? (byte1*)(a)+(a)->sz-sizeof(int):(byte1*)(a))))

#define	EMLMPR(a,e,i)	((a)->prc[ISTRANSPM(e)][i])
#define	EMLMPC(a,e,i)	((a)->prc[!ISTRANSPM(e)][i])
#define	EMLMXR(a,e,i)	((a)->xrc[ISTRANSPM(e)][i])
#define	EMLMXC(a,e,i)	((a)->xrc[!ISTRANSPM(e)][i])
#define	EMLMGR(a,e,i)	((a)->grc[ISTRANSPM(e)][i])
#define	EMLMGC(a,e,i)	((a)->grc[!ISTRANSPM(e)][i])

#define	EMLMSR(a,e)	((a)->ns[ISTRANSPM(e)])
#define	EMLMSC(a,e)	((a)->ns[!ISTRANSPM(e)])
#define	EMLMDR(a,e)	((a)->nd[ISTRANSPM(e)])
#define	EMLMDC(a,e)	((a)->nd[!ISTRANSPM(e)])


emlinemap*	strmap_new_ext(ematrix *es, ematrix *ed, int init) ;
#define	strmap_new(es,ed)	strmap_new_ext(es,ed,0)
#define	strmap_new_id(es,ed)	strmap_new_ext(es,ed,1)
/* (it is enough to "free(a)" to dispose a line map) */

emlinemap*      strmap_copy(emlinemap *a) ;

int     strmap_check_ext(int ch, emlinemap *a, ematrix *es, ematrix *ed) ;
	/* (errors are printed out only if ch>0, return value of -1 indicates an error) */
#define	strmap_check_size(a,es,ed)	strmap_check_ext(-1,a,es,ed)
#define	strmap_check(a,es,ed)		strmap_check_ext(1,a,es,ed)
#define	strmap_check_comp(a,es,ed)	strmap_check_ext(3,a,es,ed)

void    strmap_fprint(FILE *fout, emlinemap *a, ematrix *e, char *pref) ;
#define	strmap_print(a,e,p)	strmap_fprint(stdout,a,e,p)


/**
 * This function translates the given line map a (es->ed) into a refering submatrix to ed.
 * The map a is stored in the resulting refering matrix, so it must NOT be freed separately.
 * The new matrix is returned.
**/

ematrix*        strmap_tosubmatrix(ematrix *es, ematrix *ed, emlinemap *a) ;


/**
 * This function is called to find injective line maps of columns of the matrix efr to columns
 * of the matrix eto(x) (unscaled mappings if scal==0) that respect fixed rows of these matrices.
 * The matrices to map are in efr, eto - they need not have the same number of columns,
 * but they must have the same number of rows.
 * Alternatively, etox is a larger matrix from which we select the rows given in the row
 * map ax if given (and this is reflected in the output maps in *al).
 * The arrays hscfr[],hscto[] may give column hash codes - if a column j of efr can be mapped
 * to a column jj of eto, then must be hcfr[j]==hcto[jj].
 * If the list afo is given, then its col map is supposed to be a permutation group  
 * on the columns of efr, and these permutations are factored out from all col maps here.
 * 
 * The return value is -1 for no mapping found, or 0 for existing mappings but al==NULL,
 * or the number of mappings found and stored in the list *al.
 * If ch>=4, then all mappings are found and printed out.
 * The parameter ll gives the number of previously printed maps for printing results.
 * 
 * One may "transpose" the generated maps simply by transposing the matrices efr,eto.
 * To conserve memory, we never generate more than SM_MAXLINEMAPS.
 * (Matrix connectivity is no longer required here.)
**/

#ifndef SM_MAXLINEMAPS
#define SM_MAXLINEMAPS	400000
#endif

int     strmap_fixedrows_ext(int ch, ematrix *efr, ematrix *etox, int scal, emlinemap **afo,
                                long hcfr[], long hcto[], emlinemap *ax, emlinemap ***al, int ll) ;

/**
 * The next functions are similar to strmap_fixedrows_() above, but they are less restrictive
 * on generated maps - just the same number of columns, or all maps in different ways.
 * To make the computation faster, call always the most restrictive of these functions
 * you can, or the most suitable one according to the function description in strlmap.c.
 * 
 * Possible line hash codes given in hrxx[],hcxx[] are always fully respected in all
 * generated maps.
 * If the list afo is given, then it is supposed to form permutation groups on the lines
 * (rows/columns) of efr, and these permutations are factored out from all maps here.
 * The return values are the same as above, in particular -1 (<0) for no mapping found.
 * 
**/

int     strmap_samecols_ext(int ch, ematrix *efr, ematrix *etox, int scal, emlinemap **afo,
                                long hrfr[], long hrto[], long hcfr[], long hctox[],
                                emlinemap *ax, emlinemap ***al) ;
#define	strmap_samelines_ch(ch,efr,eto,sc,al)	strmap_samecols_ext(ch,efr,eto,sc,NULL,NULL,NULL,NULL,NULL,NULL,al)
#define	strmap_samecols(efr,eto,afo,h1,h2,h3,h4,al)	strmap_samecols_ext(0,efr,eto,1,afo,h1,h2,h3,h4,NULL,al)
#define	strmap_samelines_hh(efr,eto,sc,h1,h2,h3,h4,al)	strmap_samecols_ext(0,efr,eto,sc,NULL,h1,h2,h3,h4,NULL,al)

int     strmap_allbysubs_ext(int ch, ematrix *efr, ematrix *eto, int scal, emlinemap **afo,
                        long hrfr[], long hrto[], long hcfr[], long hcto[], emlinemap ***al) ;
	/* (All (unordered) subsets of columns are tried here, so consider this for a transposition.) */
#define	strmap_allbysubs(efr,eto,afo,h1,h2,h3,h4,al)	strmap_allbysubs_ext(0,efr,eto,1,afo,h1,h2,h3,h4,al)

//************ strmap_allbyorder_ext() transferred later....









/******************	Matrix minors (strminor.c)	*******************/
/**************************************************************************/


/**
 * These are supplementary functions for matroid minor tests.
 * The first one (also used globally) asks for a displayed minor in a matrix.
 * The second one looks for a represented minor in all bases of a matrix,
 * while the third one looks for a represented minor by trying all line removals.
 * If rnd>0, then only a randomized test is performed (for a sublist of bases).
 * 
 * The return values are similar as above, in particular <0 for no mapping found.
 * Call struct_hasminor() or struct_allminors() from global functions instead.
 * 
**/

/*int     struct_dispminor_ext(int ch, ematrix *e, ematrix *emin, void **afox, ematrix ***mls,
                                long hrem[], long hre[], long hcem[], long hce[], int rnd) ;*/
#define	struct_dispminor_number(e,em,afo)	struct_dispminor_ch(1,e,em,afo,NULL)
#define	struct_dispminor_eqv(es,ed,hes,hed)	struct_dispminor_ext(0,es,ed,NULL,NULL,hes,(hes)+ROWSM(es),hed,(hed)+ROWSM(ed),0)


int     struct_minorbybas_ext(int ch, ematrix *e, ematrix *emin, ematrix **basl, void **afox,
                                ematrix ***mls, long hhe[], long hhem[], int rnd) ;
#define	struct_minorbybas(ch,e,emin,mls,rnd)	struct_minorbybas_ext(ch,e,emin,NULL,(void**)1,mls,NULL,NULL,rnd)
#define	struct_minorbybas_all(ch,e,emin,af,mls)	struct_minorbybas_ext(ch,e,emin,NULL,af,mls,NULL,NULL,0)
#define	struct_minorbybas_print(ch,e,emin,af)	struct_minorbybas_all(ch,e,emin,af,NULL)
#define	struct_minorbybas_bl(ch,e,emin,bl,mls,hhe,hhem)	struct_minorbybas_ext(ch,e,emin,bl,NULL,mls,hhe,hhem,0)

#ifndef STR_MAXMINREM
#define STR_MAXMINREM   6	/* (cannot have bigger size difference for struct_minorbyrem) */
#endif
int     struct_minorbyrem_ext(int ch, ematrix *e, ematrix *emin, int rnd) ;
#define	struct_minorbyrem(ch,e,em)	struct_minorbyrem_ext(ch,e,em,0)


/**
 * Testing matrix equivalence is implemented separately from minor testing, and
 * it is much faster in general (for nonsymmetric matroids).
 * The point is that we may use equivalence for minor testing when the minor is large
 * (see struct_minorbyrem_ext() in strminor.c).
 * !! There is one important difference here from struct_dispminor_ext() and others --
 *   the line hash-codes hhs[],hhd[] (indexed by rows and then by columns) are not
 *   guaranteed to be reflected in the maps, since they are locally modified inside the
 *   function !!
 * 
 * The return values are similar as above, in particular <0 for no mapping found.
 * If rnd>0, then only a randomized test is performed - for a sublist of bases.
 * So rnd>0 may return NO even if equivalent, but a YES answer is always sure.
 * The random return value is -10 if the NO answer is already sure.
**/

int     struct_matequivalence_ext(int ch, ematrix *es, ematrix *ed, ematrix ***mls,
                                        long hhs[], long hhd[], long fls[], long fld[], int rnd) ;
#define	struct_matequivalence_mbr(ch,es,ed,fls,fld)	struct_matequivalence_ext(ch,es,ed,NULL,NULL,NULL,fls,fld,0)
#define	struct_matequivalence_rn(ch,es,ed,rn)		struct_matequivalence_ext(ch,es,ed,NULL,NULL,NULL,NULL,NULL,rn)



/**
 * This is a fast special check for a triad/triangle used when looking for fans in strconn.c.
 * 
**/

int     struct_hasfan_triax(int ch, ematrix *e, int trt, int l1, int l2, int l3, int nps) ;















#endif	/* (of #ifndef STR_H) */























