
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
 * Some theory for start... about structure of matroid representations:
 * 
 * This part contains fuctions for "structural work" with represented matroids.
 * (We do not attempt to work with matroids without a matrix representation - too complicated.)
 * Read ematrix.h to learn about our matrices, and pfield.h for the pfield arithmetic.
 * Next you find a description of our "represented matroids" and of representation "equivalence".
 * 
 * 
 * We consider matrix representations in the standard form - i.e. without the leading unit matrix.
 * The matrix lines are labeled with the matroid elements, here given initially
 * as 1...r and -1...-r*.
 * A particular matrix representation is determined by a choice of a basis in it
 * - forming the rows of the matrix.
 * The representation is invariant on scaling of matrix lines, and we may as well consider
 * it as invariant on line-permutations (provided we permute also the line-labels).
 * However, it is often good to stick with certain line-orders in a representation.
 * Sometimes we do not even care about the line-labels ("unlabeled" matroids).
 * 
 * There is one BIG PROBLEM with matroid representations - caused by the fact that
 * non-equivalent representations exist for the same matroid.
 * Two matroid representations are "strongly equivalent" if one can be obtained from the other
 * one by pivoting (changing basis), row/column permutations, and nonzero scales.
 * (In this definition, we stick the line-labels with the lines, and we do not allow
 * pfield automorphisms.)
 * However, the nature of this program needs also a slightly weaker definition of an
 * "unlabeled (strong) equivalence", when a matrix strongly equivalent to the first matrix
 * is matroid-isomorphic to the second matrix.
 * (That means we forget the line-labels from the representations.)
 * 
 * For example, the matroid U24 has two strongly-non-equivalent representations over
 * GF(4), but these two are unlabeled-equivalent.
 * These problems must be taken seriously especially when looking for matroid minors
 * or matroid isomorphism.
 * Similar (but somehow different as well) problems exist when generating matrix
 * extensions, as described in gener.h ...
 * 
 * 
 * More description of the "structural" functions is given in the sections below...
 * Look also at struct/str.h for some local definitions and routines.
 * 
**/



#ifndef STRUCT_H 
#define STRUCT_H 











/******************     Matroid connectivity (strconn.c)	*******************/
/**********************************************************************************/


/**
 * See the function ematrix_whatsep() in ../ematrix/ematrixop.c for computing the
 * connectivity function of a separation.
 * We use the formula for the connectivity of (A,B) as lambda(A) = r(A)+r(B)-r(M).
 * So the value of this function is 0 for 1-separation, 1 for exact 2-separation, etc...
 * 
 * A matroid is called k-connected (k>=2) if it has no (k-1)-separation of lambda <k-1.
 * (A connected matroid has a connected matrix, and it corresponds to a simple
 * vertex-2-connected graph. We mostly consider 3-connected matroids.)
 * 
**/


/**
 * The first function computes the "connected order" of columns/rows (tr==1/0) of the matrix e.
 * Connected order means the order in which the columns of e should be taken (starting with 0)
 * so that the partial matroid stays connected up to possible zero lines(!).
 * If the matrix is disconnected, the full order is still computed so that it is
 * connected on each component.
 * The return value is 1 for connected order, and -1 for disconnected.
 * The order is stored to the given array co[], or nowhere if co==NULL is given.
 * If cp is given, then the indices (1,2..) of components are stored there for the columns.
 * (These are indexed the same as cor[], or by absolute col indices if cor=NULL.)
 * If cxp is given, then the same component indices are recorded also for (perpendicular)
 * rows/columns of the matrix, and zero lines get index 0.
 * The function may also compute an adjacency matrix for the columns/rows of e,
 * where -1 means no adjacency, and p>=0 means adjacency in a row/column p.
**/

int     struct_connorder_ext(ematrix *e, int tr, int *cor, int *cp, int *cxp, int **cadj) ;
#define	struct_connorder_col(e,co)		struct_connorder_ext(e,1,co,NULL,NULL,NULL)
#define	struct_connorder_radj(e,co,cp,cad)	struct_connorder_ext(e,0,co,cp,NULL,cad)
#define	struct_connorder_rcp(e,cp)		struct_connorder_ext(e,0,NULL,cp,NULL,NULL)
#define	struct_connorder_rcxp(e,cp,cx)		struct_connorder_ext(e,0,NULL,cp,cx,NULL)


/**
 * The next function tests/determines the connectivity of the given matroid.
 * If cn==-1, then the function returns the actual connectivity of the matroid (>=1).
 * If cn>0, then the function tests whether the matroid is (at least) cn-connected.
 * If sepout is given (must be initialized to *sepout==NULL!), then *sepout collects
 * the list of all interesting small separations.
**/

int     struct_connectivity_ext(int ch, ematrix *e, int cn, int *ixcon, ematrix ***sepout) ;
#define	struct_connectivity(e)	struct_connectivity_ext(0,e,-1,NULL,NULL)
#define	struct_iconnectivity(e,ixc)	struct_connectivity_ext(0,e,-1,ixc,NULL)
#define	struct_isconnected(e,c)	(struct_connectivity_ext(0,e,((c)>2?(c):2),NULL,NULL)>=0)


/**
 * This function computes induced paths between vertices of the graph in cadj[][].
 * For each pair of vertices in this graph we find an induced path joining them,
 * and using only other verticess before (=not larger than the second one in cor[])
 * if the vertex order cor[] is given.
 * (These paths need not be the shortest, just induced.)
 * We record the internal vertices of these paths in the bitfields ppat[][],
 * and the cadj-values needed for the paths in the bitfields mpat[][].
 * The paths are all indexed in order of cor[] if it is given.
**/

int	struct_connpaths_ext(int sz, int **cadj, int *cor, long **ppat, long **mpat) ;
#define	struct_connpaths(sz,cad,cor,pp,mp)	struct_connpaths_ext(sz,cad,cor,pp,mp)


/**
 * The function tests (co)simplicity of the given matroid (tr=1 simple / tr=0 cosimple).
 * The return value is 1 for (co)simple, and -lo-1 if line lo is parallel to something else.
**/

int	struct_simplicity_ext(int ch, ematrix *e, int tr) ;
#define	struct_issimple(e)	struct_simplicity_ext(0,e,1)
#define	struct_iscosimple(e)	struct_simplicity_ext(0,e,0)


/**
 * This function tests/computes the shortest cycle of the given matroid.
 * If cc==-1, then the function returns the length of the shortest cycle (>=1).
 * (If the matroid has no cycle, then the return value is 9999.)
 * If cc>=0, then the function returns number cc if the length of the shortest cycle is >=cc,
 * but it returns -1 if there is a cycle shorter than cc.
 * If depout is given (must be initialized to *depout==NULL!), then *depout collects
 * the list of all shortest cycles in e, and the return value is the length.
**/

int	struct_matgirth_ext(int ch, ematrix *e, int cc, ematrix ***depout) ;
#define	struct_matgirth(e)	struct_matgirth_ext(0,e,-1,NULL)
#define	struct_hasmatgirth(e,g)	struct_matgirth_ext(0,e,g,NULL)
#define	struct_hasmatgirth_list(e,g,l)	struct_matgirth_ext(0,e,g,l)




/**
 * We provide a simple greedy (but precise) algorithm to test for branch-width 3 of a matroid.
 * 
 * We use the fact that any 3-separating set of size >=4 can be displayed in some
 * bw-3 decomposition, even when considering an already partitioned matroid
 * [Hall, Oxley, Semple, Whittle].
 * On the other hand, if a bw-3 decomposition exists, then there always is 3-separating
 * subset of <=6 elements, or <=4 blocks if using bigger than 1-element parts.
 * Hence we just look for such possible 3-separating set.
 * If we find one, we make its union as a new block in the matroid.
 * If we succeed all the way to just two blocks at the end, we have a bw-3 decomposition.
 * Otherwise, the matroid has no such decomposition.
 * 
 * We refer to the paper "On the Excluded Minors for Matroids of Branch-Width Three"
 * which is submitted for publication and available on the MACEK web page for viewing.
**/

/**
 * This function tests branch-width 3 of the given matroid ee.
 * The return value is 1 for branch-width <=3 (printed via *decomp), and -1 if not.
 * (ldecomp is not yet implemented.)
 * The function requires the matroid to be 3-connected !!!
**/
  
int     struct_hasbwidth3_ext(int ch, ematrix *ee, char **decomp, int *ldecomp) ;
#define	struct_hasbwidth3_sprint(e,buf)	(struct_hasbwidth3_ext(0,e,buf,NULL)>=0)
	/* (the printed output is returned via *buf, and it must be freed later) */
#define	struct_hasbwidth3(e)		struct_hasbwidth3_sprint(e,NULL)



/**
 * A fan of length d is a sequence of d element where consecutive triples switch between
 * being triangles and triads.
 * The function struct_hasfan..() looks whether the matrix e has a "fan" of the given length,
 * moreover ending with the given line lo if tr>=0 (row for tr==0, column for tr==1).
 * If one asks about fans ending with the given line, then the matroid must be 3-connected
 * after this line is removed(!), otherwise the answer is incorrect.
 * One may also ask for a fan only in the given submatrix (ec) of the whole matrix.
 * 
**/

int     struct_hasfan_ext(int ch, ematrix *e, ematrix *ec, int tr, int lo, int fan) ;
#define struct_hasfan_ch(ch,e,tr,lo,fan)	struct_hasfan_ext(ch,e,NULL,tr,lo,fan)
#define struct_hasfan(e,tr,lo,fan)		struct_hasfan_ch(0,e,tr,lo,fan)
#define struct_hasfan_full(e,fan)		struct_hasfan_ch(0,e,-1,-1,fan)
#define struct_hasfan_max(e)			struct_hasfan_full(e,-1)

#define struct_hasfan_print(e,fan)		struct_hasfan_ch(2,e,-1,-1,fan)
#define struct_hasfan_printall(e,fan)		struct_hasfan_ch(3,e,-1,-1,fan)

#define struct_hasfan_part(ch,e,ec,fan)		struct_hasfan_ext(ch,e,ec,-1,-1,fan)














/******************	Matrix Minors (strminor.c)	*******************/
/**************************************************************************/


/**
 * The word "minor" is used in its matroidal sense, but that has little meaning in
 * matrices alone.
 * Thus we speak about a "displayed minor" which means a submatrix displayed in
 * a suitable representation obtained by a choice of a basis (via pivoting and scaling).
 * (We shortly say a "dispminor".)
 * A "represented minor" is a minor that may be displayed in some basis.
 * So be careful when working with matroids having nonequivalent representations.
 * Read about general matroid representations at the top.
 * 
 * In general, there are two ways how to search for a represented minor --
 * we may look at dispminors among all matrices for all choices of bases for the matroid,
 * or we may try all deletions/contractions in the matrix to get an equivalent representation.
 * For the first way, we implement a function for finding dispminors, using the
 * above line maps of matrices.
 * For the second one, we implement a (fast) equivalence testing function that
 * uses various matroid-invariant properties of a matrix for faster rejections.
 * (However, in some bad cases we have to search through all bases anyway...)
 * 
 * The functions here work with matrix lists (ematrix**) that are variants of "alist"
 * from misc.c; these allows to handle collections of an arbitrary length.
 * These lists, as well as the matrices in them, should be freed by caller when
 * no longer used, by calling dispose_alist_mats(ls).
 * 
 * 
 * We shortly say that we consider "bases" of a matrix when speaking about all
 * pivoted matrices displaying all bases of the underlying matroid.
 * When reporting displayed minors, we usually factor out all self- line maps of
 * the minor matrix (it is faster).
 * 
 * Keep in mind that when asking about an abstract minor in the given matrix,
 * one has to give a list of all representations of the minor up to unlabeled strong
 * equivalence.
 * 
 * Read particular descriptions of the functions in struct/strminor.c and in struct/str.h.
**/



/**
 * The following functions produce a list of (displayed) bases for use in the next functions.
 * To save comp time, a list of all bases may be "recycled", i.e. stored for possible next call
 * for the same list (this must not be used for randomized bases or altered lists).
 * !!! Recycled bases are not guaranteed to have the same id labels for lines !!!
 * 
**/

ematrix**       struct_minorbases_ext(int ch, ematrix *e, int rnd) ;
#define	struct_minorbases_ch(ch,e)	struct_minorbases_ext(ch,e,0)
#define	struct_minorbases(e)		struct_minorbases_ch(0,e)
#define	struct_minorbases_rand(e)	struct_minorbases_ext(0,e,1)
#define	struct_minorbases_randlg(e)	struct_minorbases_ext(0,e,2)
#define	struct_minorbases_rn(e,rnd)	struct_minorbases_ext(0,e,rnd)
#define	struct_minorbases_norecycle(e)	struct_minorbases_rn(e,-1)

	/* one may "recycle" a full non-altered list of bases if STR_CACHEDBASES>0 by calling this: */
void    struct_minorbases_recycle(ematrix **bl) ;
void    struct_minorbases_finish(void) ;
#ifndef	STR_CACHEDBASES
#define	STR_CACHEDBASES	0	/* recycling bases is not allowed now... - to be implemented? */
#endif


/**
 * This function generates all self- line maps (i.e. row and column permutations that,
 * possibly after rescaling, identically preserve all entries) of the matrix e.       
 * If scal==0, then rescaling after permuting is not allowed.
 * 
 * Matrix line maps are described in struct/str.h.
 * General access to line maps is restricted to structural functions only, and so the
 * pointer is given as void* here.
 * However, a public function for printing a line map to a file is provided next.
**/

int     struct_matrixselfmaps_ext(int ch, ematrix *e, int scal, void ***fls) ;
#define	struct_matrixselfmaps(e,fls)		struct_matrixselfmaps_ext(0,e,1,fls)
#define	struct_matrixselfmaps_unscaled(e,fls)	struct_matrixselfmaps_ext(0,e,0,fls)
#define	struct_matrixselfmaps_print(e)		struct_matrixselfmaps_ext(2,e,1,NULL)
#define	struct_matrixselfmaps_number(e)		struct_matrixselfmaps_ext(1,e,1,NULL)

	/* one may "recycle" a full non-altered list of self-maps by this call: */
void    struct_matrixselfmaps_recycle(void *fls) ;
void    struct_matrixselfmaps_finish(void) ;

	/* this is for printing one matrix line map - given as void *a: */
void    struct_lmapfprint(FILE *fout, void *a, ematrix *e, char *pref) ;
#define EMLMAPDEBUGS(l,a,e,pf)	(junk = TESTDLEV(l)? (struct_lmapfprint(debugout,a,e,pf),0):0)
#define EMLMAPOUTPUT(a,e,pf)	struct_lmapfprint(printout,a,e,pf)


/**
 * This function compares the two given matrices for ``equality'' up to scaling (scal==1),
 * or up to matrix row maps (fixr==1), or up to matrix col maps (fixc==1).
 * The generated maps may be stored in *fls.
 * There is no use of ch so far...
 * The return value is >0 for ``equality'', and -1 for not.
**/

int	struct_mapmatrix_ext(int ch, ematrix *e1, ematrix *e2, int fixr, int fixc, int scal, void ***fls) ;
#define	struct_isscaledequal(e1,e2)	struct_mapmatrix_ext(0,e1,e2,1,1,1,NULL)
#define	struct_islmapequal(e1,e2)	struct_mapmatrix_ext(0,e1,e2,0,0,1,NULL)



/**
 * The following functions look for displayed minors in a matrix.
 * Keep in mind that a minor is found only in the given representation class
 * if there are non-equivalent representations!!! - see above.
 * The return value is -1 for no dispminor found, or 0 or the number of minors found here.
**/

int     struct_dispminor_ext(int ch, ematrix *e, ematrix *emin, void **afox, ematrix ***mls,
                                long hrem[], long hre[], long hcem[], long hce[], int rnd) ;
#define	struct_dispminor_ch(ch,e,em,afo,mls)	struct_dispminor_ext(ch,e,em,afo,mls,NULL,NULL,NULL,NULL,0)
#define	struct_dispminor(e,em,afo,mls)		struct_dispminor_ch(0,e,em,afo,mls)
#define	struct_dispminor_print(e,em)		struct_dispminor_ch(3,e,em,(void**)1,NULL)
#define	struct_isdispminor(e,em)		(struct_dispminor_ch(0,e,em,(void**)1,NULL)>=0)


/**
 * The following functions look for matroid-minors in general (displayed over all bases).
 * Keep in mind that a minor is found only in the given representation class
 * if there are non-equivalent representations!!! - see above.
 * No more than STR_MAXMINORS may be generated to save memory.
 * 
 * In struct_hasminorlist_(), the return value is i>0 if lmin[i-1] is a minor in e
 * (the first one found), or 0 if no minor.
 * No other output is collected or printed,
 * and the function is generaly faster than the next one.
 * If rnd>0, then only randomized test is performed -- against a random sublist of bases.
 * So rnd>0 may return NO even if the minor exists, but a YES answer is always sure.
 * 
 * In struct_allminors_(), return value is -1 for no minor found,
 * or 0 or the number of minors found here.
 * The collected minors are stored to the list *mls (which must be given empty or NULL)
 * as refering matrices together with copies of bases of e that display them.
 * (The copies of bases are stored in the same list, so watch for them!
 *  At the end, dispose_alist_mats() is enough...)
 * If af!=0, then the self-maps of emin are factored out from all (disp)minors.
 * 
**/

#ifndef STR_MAXMINORS
#define STR_MAXMINORS   400000
#endif

int     struct_hasminorlist_ext(int ch, ematrix *e, ematrix *emin, ematrix **lmin, int rnd) ;
#define	struct_hasminorlist_ch(ch,e,lmin)	struct_hasminorlist_ext(ch,e,NULL,lmin,0)
#define	struct_hasminorlist(e,lmin)		struct_hasminorlist_ch(0,e,lmin)
#define	struct_hasminorlist_rand(e,lmin)	struct_hasminorlist_ext(0,e,NULL,lmin,1)
#define	struct_hasminor_ch(ch,e,emin)		struct_hasminorlist_ext(ch,e,emin,NULL,0)
#define	struct_hasminor(e,emin)			struct_hasminor_ch(0,e,emin)

int     struct_allminors_ext(int ch, ematrix *e, ematrix *emin, int af, ematrix ***mls) ;
#define	struct_allminors(e,em,af,mls)	struct_allminors_ext(0,e,em,af,mls)
#define	struct_allminors_print(e,em,af)	struct_allminors_ext(3,e,em,af,NULL)
#define	struct_allminors_printall(e,em,af)	struct_allminors_ext(4,e,em,af,NULL)


/**
 * Testing matrix equivalence is implemented separately from minor testing, and
 * it is much faster in general (for nonsymmetric matroids).
 * The point is that we may use equivalence for minor testing when the minor is large
 * (see struct_minorbyrem_ext() in strminor.c).
 * !! There is one important difference here from struct_dispminor_ext() and others --
 *   the line hash-codes hhs[],hhd[] (indexed by rows and then by columns) are not
 *   guaranteed to be reflected in the maps, since they are locally modified inside the
 *   function !!
 * The return value is i>0 if ld[i-1] is equivalent to es (the first one found), or 0 if no one.
 * If ch>=3(4) is given, then we ask for all equivalence mappings from es to ed or to ld[]
 * (with no self-maps factored out), and these mappings are printed out.
 * 
 * The next function struct_equivlistself_ext() works similarly, but it tests for
 * equivalence all pairs of matrices from the given one list le.
 * The results are stored to an array *out (allocated inside function, free later!),
 * where (*out)[i]=j means that #j is the first matrix in the list equivalent to #i matrix.
**/

int     struct_isequivlist_ext(int ch, ematrix *es, ematrix *ed, ematrix **ld) ;
#define	struct_isequivalent_ch(ch,e1,e2)	struct_isequivlist_ext(ch,e1,e2,NULL)
#define	struct_isequivalent(e1,e2)		struct_isequivalent_ch(0,e1,e2)
#define	struct_isequivlist_ch(ch,e1,el)		struct_isequivlist_ext(ch,e1,NULL,el)
#define	struct_isequivlist(e1,el)		struct_isequivlist_ch(0,e1,el)

int     struct_equivlistself_ext(int ch, ematrix **le, int **out) ;
#define	struct_equivlistself_ch(ch,le,out)	struct_equivlistself_ext(ch,le,out)


/**
 * A "tight major" N of F is defined with respect to a (finite) family F of matroids.
 * N is a tight major of F if N has a minor from F and no element of N may be both contracted
 * or deleted keeping some minor isomorphic to one of F.
 * The function struct_hasdelcontr...() is provided for testing a tight major against
 * the list in tghl (in the opposite meaning of the return value).
 * The return value is 0 if e is(!) a tight major, and r>0 if the row r-1 / column r-1-ROWS
 * of e is removable in both ways.
**/

int     struct_hasdelcontrlist_ext(int ch, ematrix *e, ematrix *tgh, ematrix **tghl, int rnd) ;
#define	struct_hasdelcontr(e,tghl)		struct_hasdelcontrlist_ext(0,e,NULL,tghl,0)
#define	struct_hasdelcontr_one(e,tgh)		struct_hasdelcontrlist_ext(0,e,tgh,NULL,0)
#define	struct_hasdelcontr_print(e,tghl)	struct_hasdelcontrlist_ext(2,e,NULL,tghl,0)
#define	struct_hasdelcontr_printall(e,tghl)	struct_hasdelcontrlist_ext(3,e,NULL,tghl,0)
#define	struct_hasdelcontr_rand(e,tghl)		struct_hasdelcontrlist_ext(0,e,NULL,tghl,2)















/******************	Matroid Stuff (strmagic.c)	*******************/
/**************************************************************************/


/**
 * The functions here are pure matroidal - they reflect only the properties of the
 * underlying matroids, and not of particular representations.
 * (Some of the "magic" functions are available only locally, see in struct/str.h.)
 * 
 * The most important function here is an isomorphism testing for matroids.
 * We provide three separate functions -- the first one for testing isomorphism with
 * a fixed matroid basis, the second one for general isomorphism on a pair of matroids,
 * and the third one for testing ismorphisms in list(s) of matroids.
 * (The isomorphism functions may take optional "labels" for elements in lab1,lab2[],
 * indexed by rows and then columns, which must be respected by the isomorphisms.
 * Moreover, supplementary matroid-invariant hash codes may be given in h?[].)
 * 
 * 
 * Read particular descriptions of the functions in struct/strmagic.c and in struct/str.h.
**/


int     strmag_isomorphic_fixb(int ch, ematrix *e1, ematrix *e2, long he1[], long he2[],
                                                        long lab1[], long lab2[], int xpf2) ;

int     strmag_isomorphic_ext(int ch, ematrix *e1, ematrix *e2, int xpf2,
                                        long lab1[], long lab2[], long hl1[], long hl2[]) ;
#define	strmag_isomorphic_ch(ch,e1,e2,xp2)	strmag_isomorphic_ext(ch,e1,e2,xp2,NULL,NULL,NULL,NULL)
#define	strmag_isisomorphic(e1,e2,xp2)		(strmag_isomorphic_ch(0,e1,e2,xp2)>0)

int     strmag_isomorphic_map(int ch, ematrix *e1, ematrix *e2, int xpf2,
					int mp1, int mp2, long hl1[], long hl2[]) ;
#define	strmag_isautmap_h(e,mp1,mp2,hl)		(strmag_isomorphic_map(0,e,e,-1,mp1,mp2,hl,hl)>0)
#define	strmag_isautmap(e,mp1,mp2)		strmag_isautmap_h(e,mp1,mp2,NULL)

int	strmag_isomorphlist_ext(int ch, ematrix **el1, int xpf1[],
					ematrix **el2, int xpf2[], int **out) ;
#define	strmag_isomorphlist_ch(ch,el1,x1,el2,x2,out)	strmag_isomorphlist_ext(ch,el1,x1,el2,x2,out)


                                                                














/******************	Handling 3-separations (strtail.c)	*******************/
/**********************************************************************************/


#if 0

/**
 * Next we compute a "full-closure" of a submatrix (like of a subset of elements of a matroid).
 * A full closure is defined by a sequence of closures and coclosures applied to the submatrix,
 * until this process stops (related to the pathwidth of a matroid).
 * 
 * The full-closure of a triangle or a triad in the matroid is called a "tail" of width 3.
 * It is a piece of the matroid that has path-width 3.
 * The function estruct_hastail3_ext() is provided to compute that, but
 * this function works only for 3-connected matroids!!!
 * (It is possible to specify in txl the list of trangles and triads to close in the search.)
 * Additionally, one may request certain modifications to the tail length for special tail types
 * (like a fan or a line) by the parameter tmod - see the macros TMOD_xxx... and the function.
 * 
**/

ematrix*        struct_fclosure_ext(int ch, ematrix *e, ematrix *re) ;
#define	struct_fclosure_ch(ch,e,re)	struct_fclosure_ext(ch,e,re)
#define	struct_fclosure(e,re)		struct_fclosure_ch(0,e,re)

int     struct_hastail3_ext(int ch, ematrix *e, int tr, int lo, ematrix **txl, int tail, int tmod) ;
#define struct_hastail3_ch(ch,e,tl,tm)		struct_hastail3_ext(ch,e,-1,-1,NULL,tl,tm)
#define struct_hastail3(e,tl,tm)		struct_hastail3_ch(0,e,tl,tm)

/**
 * The value of TMOD.. allows to request modifications to the tail size for particular tail
 * types (size up or down); the tail types must be self-dual.
 * Precisely speaking, when a sub-tail of the particular type is found inside the whole tail
 * (or if the whole tail is of that particular type for TMOD_xxxONLY flag),
 * then the value of TMOD_xxx modifier is ADDed to the length of the particular-type tail, 
 * and the resulting tail length is adjusted accordingly up.
 * For negative modifier value, if the whole tail length plus mod falls below the particular
 * type length, then the resulting tail length is adjusted accordingly down.
 * See more inside struct_hastail3_core() ...
**/
#define	TMOD_FAN(md)		(((md)>=0? ((md)&7):8+((-(md))&7))<<5)
#define	TMOD_FANONLY(md)	(TMOD_FAN(md)|(1<<9))
#define	ISTMOD_FAN(tm)		((((tm)>>5)&7)*(((tm)&(1<<8))?-1:1))
#define	ISTMOD_FANONLY(tm)	((tm)&(1<<9))

#define struct_notail3all(e,tl,tm)		(struct_hastail3(e,tl,tm)? -1:1)
#define struct_notail3lastline(e,tr,tl,tm)	struct_notail3all(e,tl,tm)
//********** last-line oriented checking is not yet implemented (would use txl as well...)

#endif













#endif	/* (of #ifndef STRUCT_H) */


































