
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
 * Read ../include/gener.h for more information about the generating process.
 * Most of the structures and functions for generating extensions and testing elimination
 * sequences (as described in ../include/gener.h) are defined privately here.
 * See more in the corresponding source files...
 * 
**/



#ifndef GEN_H 
#define GEN_H 











/******************	Elimination sequence	*******************/
/******************************************************************/


#include "eseq.h"



/**
 * These are supplementary functions used when extending a matrix,
 * and some profiling global variables (not affecting the program computation).
 * In source files gmatext.c and gener.c.
 * Read ../include/gener.h for more information...
**/

int     gener_matextens_cdet(ematrix *ex, int rcon[], long **ppat, long **mpat, ematrix *edet) ;
int     gener_matextens_cycl(int k, long **ppat, ematrix *edet, int selsg[], int selex[], int sel1, int rnf) ;
int     gener_matextens_cpscal(int k, ematrix *edet, int cco, int rcon[], int rcpi[], ematrix *ee, int ro, int co) ;

#ifndef	GEN_USEEXTPREPF
	/* whether (smaller) submatrices are pfield-checked when looking for matrix extensions */
#define	GEN_USEEXTPREPF		1
#endif

char*   gener_extprintline(int ch, ematrix *e, int trq) ;

extern int     stseq3, ststruc3, stcanon3, ststruc2, ststruc1, stcanon1, stseq0, ststruc0, stcanon0;

int     gener_extensions(int ch, int ctrl, void *vq, int trq, void ***qout) ;
int     gener_extensions_pass(int ch, int ctrl, elimseqm *qn, int tr, ematrix *ee, int ix) ;










/******************	Step-tests in generating	*******************/
/**************************************************************************/


	/* (this header is needed for work with the emlinemap structure) */
#include "../struct/str.h"


/**
 * The following functions put together all tests that must be performed at each (co)extension
 * step when reconstructing an elimination sequence.
 * In general, there are 5 test classes distinguished by what are the tested properties
 * invariant on - only line scales, also base-minor automorphisms, all sequences of a fixed
 * basis, all equivalences of the matrix, or specially the sequence signature bits only.
 * More precise description follows:
 * 
 * genstep_linecheck...()
 *  checks that look at the added line exactly up to scaling (include setting to unit-scaled form).
 * genstep_sequencecheck...()
 *  checks that look at the elimination sequence order, but are independent from
 *  base-minor automorphisms and line-scalings.
 * genstep_basischeck...()
 *  checks that look at the chosen basis for the resulting matrix, but independent from the sequence.
 * genstep_structurecheck...()
 *  checks that purely look at (matroidal) structural properties, invariant on equivalence.
 * genstep_seqsigncheck...()
 *  checks that look only at the sequence signature and matrix dimensions, not at entries.
 * 
 * The value of tr tells what was the last added line (row if tr==0, column if tr==1) of the
 * matrix in test - given against the current(!) transposition state of it.
 * If the value tr==-1 is given, then the tests should not use information about the last added line.
 * The value of lev determines the level of test - 0 means to do all tests thoroughly,
 * while values >0 (1,3 for now) request faster but incomplete tests - to speed up rejection.
 * Keep in mind that the matrix may not be fully valid at levels lev>0, and that the matrix
 * may even not be fully representable at lev>1, so design checks carefully.
 * Specially, lev==-3 is used to request a test for test_elimseqm_ext(), which may work differently.
 * The value lev==4 is used to request the mandatory tests even when this class of tests is excluded by ctrl.
 * 
 * The value lev==0 must always perform a full test, even when this means a repetition of
 * some previous faster tests.
 * (See the functions in elimseq.c for more information.)
 * All sequence-dependent properties checked here must be hereditary(!) - that means if
 * something passed at some of the previous steps (of the sequence), then it must be OK
 * also at all subsequent steps (so that the whole sequence need not be re-checked).
 * 
 * The functions return 1 if the matrix passed the test, and -1 otherwise.
**/

int     genstep_linecheck_ext(int set, int lev, elimseqm *q, ematrix *e, int tr) ;
#define	genstep_linecheck_set(q,e,tr)		genstep_linecheck_ext(1,0,q,e,tr)
#define	genstep_linecheck_setref(q,e,tr)	genstep_linecheck_ext(2,0,q,e,tr)
#define	genstep_linecheck(q,e,tr)		genstep_linecheck_ext(0,0,q,e,tr)

int     genstep_sequencecheck_ext(int lev, elimseqm *q, ematrix *e, int tr) ;
#define	genstep_sequencecheck_lev(lev,q,e,tr)	genstep_sequencecheck_ext(lev,q,e,tr)
#define	genstep_sequencecheck(q,e,tr)		genstep_sequencecheck_lev(0,q,e,tr)
#define	genstep_sequencecheck_unord(q,e,tr)	genstep_sequencecheck_ext(-10,q,e,tr)

int     genstep_basischeck_ext(int lev, elimseqm *q, ematrix *e) ;
#define	genstep_basischeck(q,e)		genstep_basischeck_ext(0,q,e)

int     genstep_structurecheck_ext(int lev, elimseqm *q, ematrix *e, int tr) ;
#define	genstep_structurecheck_lev(lev,q,e,tr)	genstep_structurecheck_ext(lev,q,e,tr)
#define	genstep_structurecheck(q)		genstep_structurecheck_lev(0,q,NULL,-1)

#ifndef	GEN_RANDSTRUCTCH
	/* a speedup choice - use randomized structure checks at level 1 for at least this number of elements */
#define	GEN_RANDSTRUCTCH	9
#endif

int     genstep_seqsigncheck_ext(int lev, elimseqm *q, eseqs_t nsq) ;
#define	genstep_seqsigncheck(q,nsq)	genstep_seqsigncheck_ext(0,q,nsq)


/**
 * Additionally, it is possible to use pre-computed data for q in the checks - start with calling
 * genstep_precheck_start() before generating starts, and genstep_precheck_stop() afterwards.
 * These calls may use the matrix in q only without(!) the last added line.
 * (The validity of all checks must not depend on availability of precomputed data.)
 * Call genstep_precheck_isvalid() to check for correctly precomputed data.
 * 
 * An example of a pre-computed quick check is genstep_pretest_tight() below.
 * Data are also precomputed for canonical checks for use in the gcanon_...() functions.
 * This includes matrix line self-maps of the base minor and of the matrix to extend.
**/

int     genstep_precheck_ext(elimseqm *q, int tr, int wh, ematrix *e) ;
#define genstep_precheck_start(q,tr)	genstep_precheck_ext(q,tr,1,NULL)
#define genstep_precheck_stop(q,tr)	genstep_precheck_ext(q,tr,-1,NULL)
#define genstep_precheck_isvalid(q,tr,e)	(genstep_precheck_ext(q,tr,0,e)>=0)

#ifndef	GEN_USEPRECHECK
	/* a speedup choice - use some pre-computed data for tests when filtering extensions */
#define	GEN_USEPRECHECK	1
#endif

int     genstep_precheck_tight(elimseqm *q, ematrix *e, int tr) ;


/**
 * In addition to the above tests performed when generating elim sequences,
 * one may also apply a "special" sort of tests designed to filter-out certain
 * extensions (like to allow only extensions with the apecified prefix).
 * These tests must stand separately since they are NOT applied inside the
 * canonical-minimality test.
 * (That is logical since special tests are here to "kill" certain sequences
 * for good, so these sequences should not appear even in an equivalent form.)
**/

int	genstep_specialcheck_ext(int lev, elimseqm *q, ematrix *e, int tr) ;
#define	genstep_specialcheck_lev(lev,q,e,tr)	genstep_specialcheck_ext(lev,q,e,tr)
#define	genstep_specialcheck(q,e,tr)		genstep_specialcheck_lev(0,q,e,tr)



/**
 * Here we provide functions testing our canonical ordering of elimination sequences.
 * Briefly speaking, an extension sequence must be "minimal" by the canonical ordering among
 * all elimination sequences leading to equivalent (to e) matrices, to pass this test.
 * The ordering is described on the top of this file gener.h ...
 * (Remember that the sequence signatures are compared first, and only then the matrix lines
 * in their natural order and possibly other resulting-matrix properties.)
 * 
 * The value of lev determines the level of test - 0 means to do all tests thoroughly,
 * while values >0 (1,2,3 for now) request faster but incomplete tests - to speed up rejection.
 * Keep in mind that the matrix may not be fully valid at levels lev>0, and that the matrix
 * may even not be fully representable at lev>1, so design checks carefully.
 * 
 * The value of tr determines what was the last added line (0 row/1 column for the current transp).
 * This may be used for fast preliminary tests, but the full canonical test should be
 * caried out regardless of tr, taking the whole resulting matrix into an account.
 * The functions return 1 if the matrix passed the test, and -1 otherwise.
**/

int     genstep_canonicalcheck_ext(int lev, elimseqm *q, ematrix *e, int tr) ;
#define	genstep_canonicalcheck_lev(lev,q,e,tr)	genstep_canonicalcheck_ext(lev,q,e,tr)
#define	genstep_canonicalcheck(q)		genstep_canonicalcheck_lev(0,q,NULL,-1)

#ifndef	GEN_SUBSEQCOMP
	/* a speedup choice - compare (canon) incomplete subsequences when looking for a smaller one */
#define	GEN_SUBSEQCOMP	2
#endif
#ifndef	GEN_SUBSEQUNORD
	/* a speedup choice - added lines of the same type may be considered unordered */
#define	GEN_SUBSEQUNORD	1
#endif
#ifndef	GEN_RANDCANON
	/* a speedup choice - use randomized canonical check - against some random bases */
#define	GEN_RANDCANON	9
#endif




/**
 * These are supplementary functions used when checking the canonical minimality of a sequence.
 *  gcanon_fullselfmaps() uses precomputed self-maps of the previous matrix for a quick check.
 *  gcanon_allseqtest() checks minimality against all equivalent sequences, or randomized,
 *   or only against the same displayed basis.
 *  gcanon_dispseqtest() and gcanon_dispseqtest_scal() are used inside the previous function.
**/

int     gcanon_fullselfmaps(elimseqm *q, ematrix *e, int tr) ;

int     gcanon_allseqtest_ext(elimseqm *q, ematrix *e, int ran) ;
#define	gcanon_allseqtest(q,e)		gcanon_allseqtest_ext(q,e,0)
#define	gcanon_allseqtest_rand(q,e)	gcanon_allseqtest_ext(q,e,2)
#define	gcanon_allseqtest_one(q,e)	gcanon_allseqtest_ext(q,e,11)

int	gcanon_dispseqtest(int tt, elimseqm *qo, elimseqm *qn, ematrix *bs, ematrix *mn, emlinemap **aut) ;
int     gcanon_dispseqtest_scal(int tt, ematrix *e, emlinemap *am, emlinemap *au,
                                                int xc[], elimseqm *qn, elimseqm *qs) ;













/******************	Additions to generating		*******************/
/**************************************************************************/


/**
 * The functions here are provided as an "interface" for adding new attributes and
 * tests for the elimination sequence (read in ../include/gener.h).
 * These functions and additions actually appear in the file gener-more.inc.
 * Be very careful how to add new tests so that the computation is still correct,
 * the extension-generating process is very complex and fragile.
 * First make sure that you fully understand the algorithm outlined here and
 * in the file genstep.c.
 * 
**/


void    elimseqm_fprint_incl(FILE *fo, int ch, elimseqm *q, char *pref) ;

int     genstep_sequencecheck_incl(int lev, elimseqm *q, ematrix *e, int tr) ;
int     genstep_structurecheck_incl(int lev, elimseqm *q, ematrix *e, int tr) ;

int     genstep_pretest_incl(elimseqm *q, int tr, int wh) ;

int     genstep_specialcheck_incl(int lev, elimseqm *q, ematrix *e, int tr) ;

elimseqm*       gener_extframe_incl(framestruc *fr, int trq, int ctrl, elimseqm *q,
                          elimseqm **qel, framestruc **frc, framestruc **flist, int wh) ;




























#endif	/* (of #ifndef GEN_H) */








