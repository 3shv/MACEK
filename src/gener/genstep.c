
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
 * Some theory about generating matroids for start...
 * The master description of the generating process is in ../include/gener.h.
 * 
 * This file puts together all tests that must be performed at each (co)extension
 * step when generating elimination sequences -- "admissibility tests" and "canonical tests".
 * These test classes are presented in the next two sections.
 * It is very important to understand(!) and to follow the structure of these test classes
 * when adding new particaluar tests (to gener-more.inc).
 * 
 * Moreover, a new test class - "special", has been added for tests (like line prefix)
 * which are applied to filter-out certain extensions for good.
**/


#include "macek.h"
#include "gen.h"










/******************	Checking generated matrices	*******************/
/**************************************************************************/


#undef CURDLEV
#define	CURDLEV		6		/*6 (use <=DEBUGLEV+2 to perform extra checks) */


/**
 * The following functions put together all tests that must be performed at each (co)extension
 * step when reconstructing an eliminatin sequence - "admissibility tests".
 * (These do not include canonical-order testing, which is written in the next section.)
 * Read more about the general theory in the file gener.h ...
 * 
 * In general, there are 5 test classes distinguished by what are the tested properties
 * invariant on - only line scales, also base-minor automorphisms, all sequences of a fixed
 * basis, all equivalences of the matrix, or specially the sequence signature bits only.
 * Remember that one call usually performs the test only for the last added line (except structural),
 * so one must call the tests repeatedly along the sequence to determine validity of the
 * whole elim sequence from scratch (this is implicite in normal generating process).
 * More precise description follows at the functions below.
 * 
 * The given value of lev determines the level of tests - 0 means to do all tests thoroughly,
 * while values >0 (1,2,3 for now) request faster but incomplete tests - to speed up rejection.
 * The value lev==0 must always perform a full test, even when this means a repetition of
 * some previous faster tests.
 * Specially, lev==-3 is used to request a test for test_elimseqm_ext(), which may work
 * slightly differently.
 * The value lev==4 is used to request the mandatory tests (like ?-connectivity) even
 * when this class of tests is excluded by ctrl.
 * 
 * The given transposition value tr is applied against the current transposition state
 * of the matrix e (1 to col ext), which may not be 0 in general.
 * (However, the attributes inside the sequence q are stored against 0-transp.)
 * If the value tr==-1 is given, then the tests should not use information about the
 * last added line.
 * Keep in mind that the matrix may not be fully valid at levels lev>0, and that the matrix
 * may even not be fully representable at lev>=2, so design additional tests carefully.
 * 
 * 
 * All sequence-dependent properties checked here must be hereditary(!) - that means if
 * something passed at some of the previous steps (of the sequence), then it must be OK
 * also at all subsequent steps (so that the whole sequence need not be re-checked).
 * Moreover, the sequence-dependent checks must be repeated along the sequence for all added
 * lines since they work only on the last line in one call.
 * 
 * The tests should not assume at lev>=2 proper representability of the supplied
 * ext matrix, since no full representability check is performed at lev>2 (would
 * take too long to complete).
 * However, from the highest-level structure test (genstep_structurecheck(), lev==3),
 * all 2x2 subdeterminants are defined, which is enough for our usual tests.
 * 
 * Additionally, it is possible to use pre-computed data for q in the checks - by calling
 * genstep_pretest_start() before generating starts, and genstep_pretest_stop() afterwards.
 * This call may use the matrix in q only without(!) the last added line.
 * (The validity of the checks must not depend on availability of precomputed data.)
 * See  static elimseqm *precheck_seq;  and  genstep_pretest_ext()...
 * 
**/


	/* standard start for subsequent functions - fills e from q if e==NULL: */
#ifndef FASTPROG
#define	GENER_COMMON(q,e,tr)	\
		timeel = TIME();	\
		if (!e) { e = ESMATRIX(q);  if (tr>=0)  tr = (ISTRANSPM(e)!=ESELIMSG(q)%2); }	\
		else if (ROWSM(ESMATRIX(q))!=ROWSM(e) || COLSM(ESMATRIX(q))!=COLSM(e) ||	\
						ISTRANSPM(e)!=ISTRANSPM(ESMATRIX(q)))	\
			{PROGERROR("given different size/state matrix %p for seq %p",e,q);}	\
		if (tr>=0 && tr!=(ISTRANSPM(e)!=ESELIMSG(q)%2))	\
			{PROGERROR("wrong transp %d~%d compared to the seq %p",tr,ISTRANSPM(e),q);}

static int	timeel,tiel;
#define	TIMEELAPSED()	if ((tiel=TIME()-timeel)>3)  DEBUG(CURDLEV-1,"\t  - time elapsed %d.\n",tiel)
#else

#define	GENER_COMMON(q,e,tr)	\
		if (!e) { e = ESMATRIX(q);  if (tr>=0)  tr = (ISTRANSPM(e)!=ESELIMSG(q)%2); }
#define	TIMEELAPSED()	junk = 0
#endif



/**
 * This function checks the matrix e (or q) line-properties - i.e. properties dependent on the
 * particular elimination sequence and order of matrix lines (including the base minor lines).
 * (Remember that all nonzero scalings of a matrix line are "the same" for matroids.)
 * A compulsory property to check is the standard scaling form of a line.
 * If set>0, then the line is multiplied to get to the standard scaling form instead.
 * If e is a refering matrix and set>1, then the refered matrix(!!!) by e is scaled instead of e.
 * (Generaly, set>0 should set the property otherwise checked here, but it may sometimes fail.)
 * 
 * The check is carried out on the last line (row if tr==0, column if tr==1) of the
 * matrix e (which is supposed to fit as the resulting matrix of the sequence q).
 * This means that to perform a full check of the whole sequence, one must call this function
 * repeatedly along the elimination sequence!
 * If e is given, then the matrix entries in q are never accessed, only the shape.
 * These checks must not assume that the matrix is properly represented over the pfield.
 * 
 * Do not try to put here any long test, since this function is called very often.
 * In fact, there should probably be no more tests than the unit-scaled line one.
 * The return value is 1 if the matrix passed the test, and -1 otherwise.
**/

int	genstep_linecheck_ext(int set, int lev, elimseqm *q, ematrix *e, int tr) {
	int	p, r = 1;
	
	if (tr<0)  {PROGERROR("Line-check must be given the last line, it does not work on the whole matrix!");}
	GENER_COMMON(q,e,tr);
	junk = lev;		/* (level not used here) */
	
		/* so far, we only check/set the unit-scale for the matrix line, */
	p = elimchm_unitlastline(set,e,tr,ESREQCONN(q));
	r = (p>=0);		/* (lines with all-zero p==0 are also allowed here) */
	
	DEBUG(CURDLEV+2*r,"   Line-check (set %d) for q %p (e %p, tr %d) found %s\n",set,q,e,tr,r?"pass":"NO");
	TIMEELAPSED();
	return r? 1:-1;
}




/**
 * This function checks the matrix sequence-properties - i.e. properties dependent on the
 * particular elimination sequence and the added lines, but invariant on line scaling and
 * base minor self-maps.
 * (The major difference from genstep_linecheck_ext() is in allowing base minor self-maps.)
 * A sequence-connectivity check should be carried out here (see ESCONNECT).
 * 
 * The check is carried out on the last line (row if tr==0, column if tr==1) of the
 * matrix e (which is supposed to fit as the resulting matrix of the sequence q).
 * This means that to perform a full check of the whole sequence, one must call this function
 * repeatedly along the elimination sequence!
 * If e is given, then the resulting matrix entries in q are never accessed, only the shape.
 * These checks must not assume that the matrix is properly represented over the pfield for lev>1.
 * 
 * A special call at lev==-10 is used to find out whether a subsequence of additions of the
 * same type (all rows or all columns) has to be ordered or not.
 * If unordered (ret>0), then admissibility of a partial sequence after adding this subsequence
 * of lines must not depend on the internal order of these lines.
 * 
 * Do not try to put here any long tests or debug outputs, since this function is called very often
 * (especially during canonical checking).
 * The return value is 1 if the matrix passes the test, and -1 otherwise.
**/

int	genstep_sequencecheck_ext(int lev, elimseqm *q, ematrix *e, int tr) {
	int	i,tm, r = 1;
	
	if (tr<0 && lev>-3)  {PROGERROR("Sequence check must be given the last line, it does not work on the whole matrix!");}
	GENER_COMMON(q,e,tr);		/* (calls at lev==-3 from test_elimseqm() may not set tr) */
	
		/* special - asks whether a subsequence of same-type lines may be considered unordered */
	if (lev==-10) {
		r = (ESNOFAN(q)<=0 && ESNOTAIL(q)<=0);
		/* (if adding more involved tests here, better check in gcanon_dispseqtest()...) */
		if (r)  r = genstep_sequencecheck_incl(lev,q,e,tr)>0;
		return r? 1:-1;
	}
	if (lev==0 && tr<0)  {PROGERROR("Sequence check needs the last added line even at lev 0!");}
	
		/* whether the last line is added so that no parallel/series pair occurs ->
		    3-connectivity or simple/cosimple extension, or just connectivity */
	if (r && (lev>=3 || lev<=0) && tr>=0) {
		DEBUG(CURDLEV+1,"testing %dconn: %d %d, %d %d\n",ESCONNECT(q),tr,ISTRANSPM(e),ESREQSIMPLE(q),ESREQCOSIMPLE(q));
		EMATDEBUGS(CURDLEV+2,e,"\tconn-test\t");
		if ((tr!=ISTRANSPM(e))? ESREQSIMPLE(q): ESREQCOSIMPLE(q))
			r = elimchm_3connlastline(e,tr)>=0;
		else if (ESREQCONN(q))
			r = elimchm_connlastline(e,tr)>=0;
	}
		/* and whether there are no fans exceeding limit ESNOFAN(q) at this step of the sequence */
	if (r && (lev<=2 || lev<=0))  if (ESNOFAN(q)>0) {
		if ((lev==2 || lev<0) && tr>=0) {	/* (may not yet be representable here) */
			r = struct_nofanlastline_nofull(e,tr,ESNOFAN(q))>=0;
		} else if (lev<=1 && tr>=0) {	/* (full check of fans ending with the last line) */
			r = struct_nofanlastline(e,tr,ESNOFAN(q))>=0;
		} else if (lev==0) {		/* (this is not called under normal circumstances!) */
			r = struct_nofanall(e,ESNOFAN(q))>=0;
		}
	}
		/* and whether there are no tails exceeding limit ESNOTAIL(q) at this step of the sequence */
	if (r && lev<=1 && lev>-3)  if (ESNOTAIL(q)>0) {
		tm = 0;  i = ESNOFAN(q)-ESNOTAIL(q);
		PROGERROR("no-tail checking code is not implemented correctly yet!");
#if 0
		if (i>0)  tm = TMOD_FANONLY(-i);
		//************* handle tm in general....... - rewrite tail handling completly ???
		if (lev==1 && tr>=0) {		/* (a full last-line test, needs representability!) */
			r = struct_notail3lastline(e,tr,ESNOTAIL(q),tm)>=0;
		}
		if (lev<=0) {
			r = struct_notail3all(e,ESNOTAIL(q),tm)>=0;
		}
		//************ distinguish different types of tails (related to fans ???)
#endif
	}
		/* finally, here we may call user-defined extra sequence checks: */
	if (r)  r = genstep_sequencecheck_incl(lev,q,e,tr)>0;
	
	DEBUG(CURDLEV+2*r,"%*sSequence-check at lev %d for q %p (e %p, tr %d) found %s\n",lev,"",lev,q,e,tr,r?"pass":"NO");
	TIMEELAPSED();
	return r? 1:-1;
}




/**
 * This function checks the elimination sequence signature nsq against the shape of q.
 * The check has absolutely nothing(!) to do with entries of the resulting matrix in q
 * (not even with the basis), only the matrix shape and other sequence attributes are accessed.
 * 
 * It returns 1 if the signature nsq passed the test, and -1 otherwise.
**/

int	genstep_seqsigncheck_ext(int lev, elimseqm *q, eseqs_t nsq) {
	int	r = 1;
	
	junk = lev;		/* (level not used here) */
	
		/* we look whether the signature nsq has the right number of extensions and coextensions */
	r = elimchm_signbits_poss(q,nsq)>=0;
	
	DEBUG(CURDLEV+1+r,"Signature-check for q %p and nsq %lX found %s\n",q,nsq,r?"pass":"NO");
	return r? 1:-1;
}


/**
 * This function checks the matrix base-properties - i.e. properties dependent on the
 * chosen basis matrix e as if it was the resulting matrix of q (regardless of ESMATRIX(q)).
 * The checks are invariant on the sequence or on the order of lines in the matrix,
 * but they may depend on the sequence signature in q.
 * Nothing is implemented here so far, and if something will be implemented later, then
 * that would probably require to change extension-generating to consider all bases(!)...
 * 
 * Returns 1 if the matrix passed the test, and -1 otherwise.
**/

int	genstep_basischeck_ext(int lev, elimseqm *q, ematrix *e) {
	int	trx=0, r = 1;
	
	junk = lev;		/* (level not used here) */
	if (e)  trx = (ISTRANSPM(e)!=ESELIMSG(q)%2);
	GENER_COMMON(q,e,trx);
	
		/* nothign is tested here so far ... (what about representability here?) */
	
	DEBUG(CURDLEV+2*r,"Basis-check at lev %d for q %p (e %p) found %s\n",lev,q,e,r?"pass":"NO");
	TIMEELAPSED();
	return r? 1:-1;
}




/**
 * This function checks the structral properties of the matrix - i.e. properties not dependent
 * on matrix equivalence (pivoting, permuting, scaling) or the elimination sequence.
 * The checks must include full testing the matrix e in the partial field(!) at levels lev<=2.
 * Other possible tests include forbidden minors and similar stuff...
 * 
 * The matrix e is supposed to fit as the resulting matrix of the sequence q.
 * The check may use an information about the last added line (in tr) to the matrix, but it
 * should NOT entirely depend on such information - remember it is sequence-independent.
 * (For example, it is OK to use tr in a fast incomplete check at lev>0, but not at lev==0...)
 * The checks for lev>1 do not assume that the matrix is properly represented over the pfield.
 * For lev>0, we use randomized structure checks to speed-up rejection.
 * (Specially, lev==-3 is used to request a test for test_elimseqm_ext() - somehow restricted to be faster.)
 * 
 * It returns 1 if the matrix passed the test, and -1 otherwise.
**/

int	genstep_structurecheck_ext(int lev, elimseqm *q, ematrix *e, int tr) {
	int	tri, r = 1;
	
	if ((lev<=0 && tr>=0 && tr<=1) || (lev>0 && tr<0))  {PROGERROR("Wrong use of tr=%d at level %d",tr,lev);}
	GENER_COMMON(q,e,tr);
	
		/* depending on lev, we look whether the matrix e is proper in current pfield */
	if (r && lev>=3 && tr>=0)
		r = elimchm_represented_randlast(q,e,tr)>=0;
	if (r && lev==2 && tr>=0)
		r = elimchm_represented_last(q,e,tr)>=0;
	if (r && lev<=0 && lev>-3)	/* (maybe not necessary, but we formally require a FULL matrix test for lev==0) */
		r = elimchm_represented(q,e)>=0;
	if (!r)  DEBUG(CURDLEV,"     - structure: NOT in pfield found.\n");
	if (!r && lev<=1 && lev>=0)  {PROGERROR("We should find not-represented already at level >1 !");}
	TIMEELAPSED();
	
		/* we additionally look at forbidden fans regardless of the sequence order, just to be sure */
	if (r && lev==0)  if (ESNOFAN(q)>0) {
		r = struct_nofanall(e,ESNOFAN(q))>=0;
		if (!r)  {PROGERROR("We should find no-fan %d already in seq check!",ESNOFAN(q));}
	}
	
		/* for lev==3, we look whether the matrix e is a tight major against displayed minors from ESTIGHT */
	if (r && lev==3 && tr>=0)  if (ESTIGHT(q)) {
		r = genstep_precheck_tight(q,e,tr)>=0;	/* (no representability required here) */
	}
	
		/* (further checks assume properly represented matrix e over the pfield !) */
	tri = ISTRANSPM(e);
	ematrix_transpose_set(e,0);	/* (we must look at the 0-transp state) */
	printlev -= 1;
	
		/* depending on lev, we look whether the matrix e has a forbidden minor from ESFORBID */
	if (r && lev<=1 && lev>-3)
	  if (ESFORBID(q))  if (lev<=0 || (GEN_RANDSTRUCTCH>0 && ROWSM(e)+COLSM(e)>=GEN_RANDSTRUCTCH)) {
		if (lev>0) {
			r = !struct_hasminorlist_rand(e,ESFORBID(q));	/* (only fast randomized test) */
		} else {
			r = !struct_hasminorlist(e,ESFORBID(q));
		}
	}

		/* depending on lev, we look whether the matrix e is a tight major for the list ESTIGHT */
	if (r && lev<=1 && lev>-3)
	  if (ESTIGHT(q))  if (lev<=0 || (GEN_RANDSTRUCTCH>0 && ROWSM(e)+COLSM(e)>=GEN_RANDSTRUCTCH)) {
		if (lev>0) {
			r = !struct_hasdelcontr_rand(e,ESTIGHT(q));	/* (only fast randomized test) */
		} else {
			r = !struct_hasdelcontr(e,ESTIGHT(q));
		}
	}
	ematrix_transpose_set(e,tri);
	printlev += 1;
	
		/* finally, here we may call user-defined extra structure checks: */
	if (r)  r = genstep_structurecheck_incl(lev,q,e,tr)>0;
	
	DEBUG(CURDLEV+2*r,"%*sStructure-check at lev %d for q %p (e %p, tr %d) found %s\n",lev,"",lev,q,e,tr,r?"pass":"NO");
	TIMEELAPSED();
	return r? 1:-1;
}












/******************	Pre-computing for generated sequences	*******************/
/**********************************************************************************/


#undef CURDLEV
#define	CURDLEV		6


	/* the precomputed data are here - see next genstep_pretest_start() for description: */

static elimseqm	*precheck_seq = NULL;	/* mark the sequence - to check correct data */
static ematrix	*precheck_ecopy = NULL;	/* a copy of the matrix before extension, transposed so that a column is to be added */

static long	*disptight;	/* the row patterns of displayed "tight" minors in the matrix to extend */
static int	disptlen = 0;


static emlinemap **fullselfmaps = NULL;	/* self-maps of the whole matrix without the added line */
static int	nofullselfmaps = 1;	/* (>0 if there is no non-identical self-map) */
static ematrix	*fsmaps_ev = NULL;	/* a copy of the extended matrix, with undefined new line */
		
static emlinemap **minselfmaps = NULL;	/* self-maps of the base minor - not used now... */



/**
 * This function is called to precompute some data for the extension generating process.
 * See precheck_seq and related variables above.
 * The data should not be essential to other computing, just use them to speedup some checks.
 * The function is called from gener_extensions() (gener.c) if requested.
 * For wh>0, it prepares the data for the sequence q; and for wh<0, it clears the data.
 * It returns 1 for success or -1 on failure (but this is not much useful anyway).
 * For wh==0, the input q,e is tested against precomputing copies precheck_seq,precheck_ecopy.
 * (e is used only for wh==0.)
 * 
 * The matrix in q is already given extended by the new line, so one must erase that line.
 * This function expects the added line to be given as the last column of ESMATRIX(q) in the
 * current transposition state, subject to tr (otherwise an error is reported).
 * So tr is not an actually active parameter, it is given just for control purposes.
 * 
 * Row-patterns are precomputed for displayed "tight" minors in the matrix - patterns have
 * 1-bits at rows that support the displayed minors.
 * These patterns are used in genstep_structurecheck_ext() to see whether the extended matrix
 * is "displayed" as not a tight major -- see genstep_precheck_tight().
 * 
 * More precomputing may be added in the included function genstep_pretest_incl();
 * 
 * Additionally, data are precomputed for canonical checks for use directly in the
 * gcanon_...() functions below.
 * This includes matrix line self-maps of the base minor and of the matrix to extend.
 * See more in the functions...
**/

int	genstep_precheck_ext(elimseqm *q, int tr, int wh, ematrix *e) {
	int		j,l, r;
	ematrix		**exl,**x;
	
	if (wh || !q) if (!q || (precheck_seq==NULL)!=(wh>0) || (precheck_seq==q)!=(wh<0))
		{PROGERROR("Wrong precheck call for q=%p, e=%p, wh=%d.",q,e,wh);  return -1;}
	if (wh==0 && !e)  return -1;
	if (!tr)  ematrix_transpose(ESMATRIX(q));
	if ((!ISTRANSPM(ESMATRIX(q)))!=(ESELIMSG(q)%2)) {	/* (expect a column to be added!) */
		PROGERROR("Wrong precheck call for matrix transpose.");  return -1;
	}
	r = 1;
		/* clearing pre-viously precomputed data here, _incl() additions first */
	if (wh<0) {
		genstep_pretest_incl(q,tr,wh);
		if (disptlen>0) { FREE(disptight);  disptight = NULL; }
		disptlen = 0;
		if (fullselfmaps)  struct_matrixselfmaps_recycle(fullselfmaps);
		fullselfmaps = NULL;
		if (fsmaps_ev)  dispose_ematrix(fsmaps_ev);  fsmaps_ev = NULL;
		if (minselfmaps)  struct_matrixselfmaps_recycle(minselfmaps);  minselfmaps = NULL;
		if (precheck_ecopy)  dispose_ematrix(precheck_ecopy);
		precheck_ecopy = NULL;  precheck_seq = NULL;
		DEBUG(CURDLEV-1,"Pre-computed data for seq %p cleared.\n",q);
		if (!tr)  ematrix_transpose(ESMATRIX(q));
		return r;
	}
		/* checking correct pre-computation state -- against precheck_seq,precheck_ecopy */
	if (wh==0) {
		if (!tr && e!=ESMATRIX(q))  ematrix_transpose(e);
		if (precheck_seq!=q || ROWSM(e)!=ROWSM(precheck_ecopy) || COLSM(e)-1!=COLSM(precheck_ecopy)
				|| (IFRANDDEBUGLESS(111)?elimchm_compare_col(e,0,precheck_ecopy,0)!=0:0))
			{PROGERROR("The given sequence (matrix) is not the same as the precomputed one!"); r = -1;}
		if (!tr && e!=ESMATRIX(q))  ematrix_transpose(e);
		if (!tr)  ematrix_transpose(ESMATRIX(q));
		return r;
	}
	
		/* finally, pre-computing extension data here... (wh>0): */
	DEBUG(CURDLEV-2,"Pre-computing data for generating seq %p [%s] (%d x %d+1).\n",q,ESNAME(q),ROWSM(ESMATRIX(q)),COLSM(ESMATRIX(q))-1);
	precheck_seq = q;	/* (empty precheck_seq,precheck_ecopy already tested above) */
	precheck_ecopy = ematrix_copy(ESMATRIX(q));
	ematrix_remove_col(precheck_ecopy,COLSM(precheck_ecopy)-1);
	
			/* matrix self-maps used in canonical checks below */
	fullselfmaps = minselfmaps = NULL;
	struct_matrixselfmaps(precheck_ecopy,(void***)&fullselfmaps);
	nofullselfmaps = (alist_getlength(fullselfmaps)<=1);	/* (only has sense when there is a non-identical aut) */
	fsmaps_ev = (nofullselfmaps?NULL: ematrix_copy(ESMATRIX(q)));
	struct_matrixselfmaps(ESMINOR(q),(void***)&minselfmaps);
	DEBUG(CURDLEV-1," - %d minor self-maps and %d full self-maps found\n",alist_getlength(minselfmaps),alist_getlength(fullselfmaps));
	
			/* row-patterns for displayed "tight" minors in the matrix precheck_ecopy */
	disptlen = 0;
	if (ESTIGHT(q) && ROWSM(ESMATRIX(q))<NUM_LONG_BITS-4 \
				&& ROWSM(ESMATRIX(q))-ROWSM(ESMINOR(q))>2) {
		exl = NULL;
		for (x=ESTIGHT(q); x?*x:0; x++)		/* (the dispminors are all collected in exl) */
			struct_dispminor(precheck_ecopy,*x,(void**)1,&exl);
		disptlen = alist_getlength(exl);
		disptight = MMALLOC((disptlen+2)*sizeof(disptight[0]));
		for (x=exl,l=0; x?*x:0; x++,l++) {
			disptight[l] = 0;
			for (j=0; j<ROWSM(*x); j++)  disptight[l] |= (1l<<GETREFMROW(*x,j));
		}
		dispose_alist_mats(exl);
		if (l>0)  DEBUG(CURDLEV-1," - row patterns for %d displayed \"tight\" minors found\n",l);
	}
		/* calling possible third-side pre-test function for starting */
	genstep_pretest_incl(q,tr,wh);
	
	if (!tr)  ematrix_transpose(ESMATRIX(q));
	return r;
}


/**
 * This is a fast tight-major pre-test, based on data precomputed above.
 * For each occurence of a displayed minor from the "tight" list in the original matrix,
 * we remembered the row pattern of the minor.
 * If we now find that the non-zero pattern of the added column is disjoint from some
 * of these row patterns, then we know that the extended matrix is not a tight major;
 * any nonzero may be pivoted not changing the displayed minor from the "tight" list.
**/

int	genstep_precheck_tight(elimseqm *q, ematrix *e, int tr) {
	long	lzp=0;
	int	i,j, r=1;
	
	if (q==precheck_seq && disptlen>0) {
		if (tr<0)  {PROGERROR("Wrong use of tr=%d at level>0",tr); return 1;}
		if (!tr)  ematrix_transpose(e);
		if (!genstep_precheck_isvalid(q,tr,e))  r = 0;	/* (error reported inside) */
		if (r) for (i=lzp=0; i<ROWSM(e); i++)
			if (SIGNM(e,i,COLSM(e)-1)!=0)  lzp |= (1l<<i);
			/* test against rows of "tight" minors - must miss some pattern to be tight: */
		for (j=0; j<disptlen && r; j++)
			if ((disptight[j]&lzp)==0)  r = 0;
		if (!tr)  ematrix_transpose(e);
		DEBUG(CURDLEV+2*r,"     - fast tight-major test for %p found %s\n",e,r?"pass":"NOT");
	}
	return r?1:-1;
}














/******************	"Special" check class	*******************/
/******************************************************************/


#undef CURDLEV
#define	CURDLEV		6

/**
 * In addition to the above tests performed when generating elim sequences,
 * one may also apply a "special" sort of tests designed to filter-out certain
 * extensions (like to allow only extensions with the apecified prefix).
 * These tests must stand separately since they are NOT applied inside the
 * canonical-minimality test.
 * (That is logical since special tests are here to "kill" certain sequences
 * for good, so these sequences should not appear even in an equivalent form.)
**/

int	genstep_specialcheck_ext(int lev, elimseqm *q, ematrix *e, int tr) {
	int	r = 1;
	
	if (tr<0)  {PROGERROR("Special check must be given the last line, it does not work on the whole matrix!");}
	GENER_COMMON(q,e,tr);
	
		/* prefix - see the option '@ext-entries "e n t r i e s"' and gener_matextens_() */
	if (r && (lev==3 || lev<=0)) if (ESPREFIXM(q)) {
		r = (elimchm_prefixlastline(e,tr,ESPREFIXM(q))>=0);
	}
		/* finally, here we may call user-defined extra sequence checks: */
	if (r)  r = genstep_specialcheck_incl(lev,q,e,tr)>0;
	
	DEBUG(CURDLEV-1+3*r,"%*sSpecial-check at lev %d for q %p (e %p, tr %d) found %s\n",lev,"",lev,q,e,tr,r?"pass":"NO");
	TIMEELAPSED();
	return r? 1:-1;
}






























/******************	Checking canonical forms	*******************/
/**************************************************************************/



/**
 * The following functions put together all tests that must be performed at each (co)extension
 * step when reconstructing an elimination sequence in order to refuse repeated occurences of
 * an equivalence class of matrices starting with the (identically) same minor.
 * Briefly speaking, an extension sequence must be "minimal" by the canonical ordering among
 * all elimination sequences leading to equivalent (to e) matrices, to pass this test.
 * Read more about the related theory in the file ../include/gener.h ...
 * (Remember that the sequence signatures are compared first, and only then the matrix lines
 * in their natural order and possibly other resulting-matrix properties.)
 * 
 * It is vital that the functions here use exactly the same genstep checks to test admissibility
 * of a matrix as the generating process in gener.c - otherwise some extensions may be lost(!).
 * Besides those, we use the elimchm_signbitscomp() and elimchm_seqcompare() tests for finding
 * which sequence (matrix) is smaller, to reject sequences that are not canonically-minimal.
 * 
 * General scheme of the canonical testing follows:
 *  - cycle all admissible signatures that are not larger than that one of current sequence;
 *    - cycle all bases of the current matrix, check their admissibility;
 *      - cycle all displayed base-minors for the chosen basis (self-maps factored out);
 *        - cycle all orders of the added lines with respect to the chosen signature,
 *          check admissibility of the constructed sequence (no line-scales yet);
 *          - cycle base-minor automorphisms, apply them to the constructed matrix,
 *            set up the displayed minor identically,
 *            set unit-scales for all added lines in the sequence,
 *            compare the constructed sequence with the current one.
 * To speed-up the computation, we add a pre-comparing check for subsequences during cycling
 * all orders of added lines if SUBSEQCOMP>0 is defined.
 * (For this to work, the ordering of sequences must be hereditary on subsequences!)
 * We also implement randomized pre-check if RANDCANON>0 - using just random subselection of bases.
 * 
 * No third-side additions are possible to the canonical testing, as it would not have much sense
 * and might introduce serious errors...
 * 
**/


#undef CURDLEV
#define	CURDLEV		6


#ifndef FASTPROG
static int	timeelc,tielc;
#define	GENER_COMMON_C(q,e,tr)	timeelc = TIME();  GENER_COMMON(q,e,tr)
#define	TIMEELAPSED_C()		if ((tielc=TIME()-timeelc)>2) DEBUG(CURDLEV-1,"\t  - time canonical elapsed %d.\n",tielc)

#else
#define	GENER_COMMON_C(q,e,tr)	GENER_COMMON(q,e,tr)
#define	TIMEELAPSED_C()		junk = 0
#endif



/**
 * This function checks whether the matrix e is minimal by the canonical ordering among
 * all equivalent (to e) matrices and all their sequences.
 * The matrix e is supposed to fit as the resulting matrix of the sequence q.
 * (However, if lev>0, then the matrix e may not yet be fully valid generated matrix;
 * and if lev>1, then e may not even be properly represented in the pfield.)
 * 
 * The given value of lev determines the level of tests - 0 means to do all tests thoroughly,
 * while values >0 (1,3 for now) request faster but incomplete tests - to speed up rejection.
 * The value lev==0 must always perform a full test, even when this means a repetition of
 * some previous faster tests.
 * The value of tr determines what was the last added line (0 row/column 1).
 * This may be used for fast preliminary tests, but the full canonical test lev==0 should be
 * caried out regardless of tr, taking the whole resulting matrix into an account.
 * 
 * The return value is 1 if the matrix passed the test, and -1 otherwise.
**/

int	genstep_canonicalcheck_ext(int lev, elimseqm *q, ematrix *e, int tr) {
	int	r = 1;
	
	if ((lev<=0 && tr>=0 && tr<=1) || (lev>0 && tr<0))  {PROGERROR("Wrong use of tr=%d at level %d",tr,lev);}
	GENER_COMMON_C(q,e,tr);
	
		/* at lev==3 (highest), we look whether the last line is minimal for full self-maps */
	if (r && lev==3 && !nofullselfmaps) {
		r = gcanon_fullselfmaps(q,e,tr)>=0;
	}
	
	/* further checks (lev<=2) assume properly represented matrix e over the pfield ! */
	
		/* for lev==2, we check equivalent sequences displayed in the same matrix */
	if (r && lev==2) {
		r = gcanon_allseqtest_one(q,e)>=0;
	}
		/* for lev==1, we check randomized equivalent sequences and matrices */
	if (r && lev==1 && GEN_RANDCANON>0 && ROWSM(e)+COLSM(e)>=GEN_RANDCANON) {
		r = gcanon_allseqtest_rand(q,e)>=0;
	}
		/* for lev==0, we check all equivalent sequences and matrices - the full test */
	if (r && lev<=0 && lev>-3) {
		r = gcanon_allseqtest(q,e)>=0;
	}
	DEBUG(CURDLEV-(lev<=1),"%*sCanonical-check at lev %d for q %p (e %p, tr %d) found %s\n",lev,"",lev,q,e,tr,r?"pass":"NO");
	TIMEELAPSED_C();
	return r? 1:-1;
}





/**
 * This function performs a quick check of the added line against all self-maps
 * of the (before-extended part of) matrix e.
 * (These self-maps must change nothing at the previous lines/steps of the sequence,
 * so only the added line is considered, which is quite fast).
 * The check is carried out on the last line (row if tr==0, column if tr==1) of the
 * matrix e, which is supposed to fit as the resulting matrix of the sequence q.
 * The test does not assume proper representability (in pfield) of the matrix.
 * 
 * The function returns 1 if the matrix passed the test, and -1 otherwise.
**/

int	gcanon_fullselfmaps(elimseqm *q, ematrix *e, int tr) {
	int		i,lc, ru, r=1;
	emlinemap	**aut, **a;
	ematrix		*ev;
	
#ifndef	FASTPROG
	if (tr<0)  {PROGERROR("Must give tr=%d>=0",tr);}
	if (IFRANDDEBUGLESS(22))  if (genstep_linecheck(q,e,tr)<0)  {PROGERROR("We expect the added column to be unit-leading scaled!");}
#endif
		/**
		 * We prepare a "working copy" ev of the matrix e for our fast canonical check.
		 * ev is then used for applying the self-maps of ev to the extended column.
		 * If the self-maps are already precomputed in fullselfmaps, then we use them.
		**/
	if (!tr)  ematrix_transpose(e);
	ev = (fsmaps_ev? fsmaps_ev: ematrix_copy(e));
	lc = COLSM(ev)-1;  ematrix_remove_col(ev,lc);
	r = 1;
	if (fullselfmaps!=NULL)  aut = fullselfmaps;
	else { aut = NULL;  struct_matrixselfmaps(ev,(void***)&aut); }
	ev = ematrix_append_col(ev);
		/**
		 * Then we cycle through all self-maps of ev without the added column,
		 * and we apply each of them to the added column (relocation and scaling).
		 * Then we look whether we get a "smaller" line after the map --
		 * which would mean our matrix was not canonical, so it should be rejected.
		**/
	for (a=aut; *a && r; a++) {
		for (i=0; i<ROWSM(ev); i++)
			COPYEXSIGM(ev,EMLMPR(*a,ev,i),lc, e,i,lc);
		for (i=0; i<ROWSM(ev); i++)
			ematrix_multiply_entry(ev,i,lc, 1,EMLMXR(*a,ev,i),EMLMGR(*a,ev,i));
		ru = genstep_linecheck_set(q,ev,1)>=0;
			/* the following test is just a shortcut for the full test in elimchm_seqcompare() !!!*/
		if (ru)  r = elimchm_compare_col(ev,lc,e,lc)>=0;
	}
	
#ifndef FASTPROG	/* for safety, we check whether the above shortcut test was correct */
	if (ROWSM(ev)!=ROWSM(e) || COLSM(ev)!=COLSM(e) || ISTRANSPM(ev)!=ISTRANSPM(e))  {PROGERROR("Wrong working copy ev!");}
	if (aut==fullselfmaps) if (!genstep_precheck_isvalid(q,tr,e))  r = 1;  /* (error reported inside) */
	if (!r && IFRANDDEBUGLESS(222)) {
		elimseqm  *q2 = elimseqm_copyfor(q,ev);
		if (elimchm_seqcompare(q2,q)>=0)  {PROGERROR("Incorrect shortcut for sequence-compare test!");}
		FREE(q2);
		if (!struct_isdispminor(ev,e))  {PROGERROR("The smaller matrix is not equivalent to the tested one!");}
	}
	if (!r)  DEBUG(CURDLEV+0,"  rejecting full self-map #%d found -> col %s\n",(int)(a-aut),gener_extprintline(-1,ev,1));
	else  DEBUG(CURDLEV+1,"  unsuccessfull canonical try-reject with %d full self-maps\n",alist_getlength(aut));
#endif
	if (aut && fullselfmaps!=aut) {
		struct_matrixselfmaps_recycle(aut);  aut = NULL;
	}		/* (frees all the automorphisms if not precomputed in fullselfmaps) */
	if (ev && ev!=fsmaps_ev)  dispose_ematrix(ev);
	if (!tr)  ematrix_transpose(e);
	return r? 1:-1;
}









/******************	Canonically-minimal sequence	*******************/
/**************************************************************************/


#undef CURDLEV
#define	CURDLEV		6		/* (OK with randomized debug tests if >=5) */



/**
 * This function performs a thorough check of the extended matrix e against all elimination
 * sequences of all equivalent (to e) resulting matrices, with respect to the (identically)
 * same base minor as in q.
 * The question is whether some equivalent sequence to q is canonically smaller (-1).
 * More is contained in the general algorithm description at the top of this section.
 * The matrix e must be properly represented over the pfield for this function (ran<10).
 * 
 * If ran>0, then only randomized pre-test is carried out - not complete(!!!).
 * (One may call ran>0 even when the matrix is not yet fully validated, so be careful...)
 * If ran>10, then only the current basis of the matrix is considered.
 * If ran>20, then only the currently displayed base-minor is considered.
 * No attention is paid to the current transposition state or to the last added line
 * in the test - everything is tested against the 0-transp state of matrices.
 * 
 * The function returns 1 if the matrix e passed the test, and -1 otherwise.
 * (A return value 1 pass is only partially valid if ran>0 is given!)
**/

static int	linords, linordp, linordps, idfound = 0;	/* (used for debugging only) */

int	gcanon_allseqtest_ext(elimseqm *q, ematrix *e, int ran) {
	eseqs_t		lq;
	int		i,k,t, nlq=0,nbm=0, ret=1;
	ematrix		*ee,*esm,*re, **bas,**b, **bmin,**m, ***bmlists,***bx;
	elimseqm	*qn;
	emlinemap	**aut;
	extern int	gener_extensions_ix;	/* (used only for debug printing of seq number) */
	
#ifndef FASTPROG
	if ((q&&e)? ROWSM(ESMATRIX(q))!=ROWSM(e):1)  {PROGERROR("Wrong given parameters q,e ! %p,%p",q,e);}
	if (ESLENGTH(q)>GEN_MAXELIMLENGTH)  {PROGERROREXIT("cannot work with elim seq longer than %d > %d",ESLENGTH(q),GEN_MAXELIMLENGTH);}
	linords = linordp = linordps = 0;	/* (used for debug - profiling) */
	idfound = 0;		/* (used for debug - whether the identical sequence was tested) */
	DEBUG(CURDLEV-1+(ran>3),"Going to look for %s equivalent sequences (seq %p, mat %p) ...\n",ran>0?(ran<10?"random":"onedisp"):"all",q,e);
	ELIMQDEBUGS(CURDLEV+2,q,ran<=0?"\t\t|=\t":"\t\t|-\t");
	if (IFRANDDEBUGLESS(22))  test_elimseqm_notough(q);
#endif
		/**
		 * We prepare data for our canonical-minimality test:
		 * ee is a copy of the matrix in q, transposed to 0-tr as required by the
		 *  definition of an elim sequence, and stored in the working seq qn.
		 * bas gets all the bases that we display in our test, also in 0-transp state.
		 * aut gets all self-maps of the base minor esm, possibly using precomputed list.
		 * But first we call the test for the default choice of the base minor
		 * in the current basis, trying all smaller signatures (quite fast).
		**/
	ret = 1;
	esm = ESMINOR(q);  aut = NULL;
	if (minselfmaps!=NULL)  aut = minselfmaps;
	else  struct_matrixselfmaps(esm,(void***)&aut);
	ee = ematrix_copy(e);
	qn = elimseqm_copyfor(q,ee);
	ematrix_transpose_set(ee,0);
	re = ematrix_refer(ee,0,ROWSM(esm),0,COLSM(esm));
	EMDATA(re) = strmap_new_id(re,ee);
	for (lq=0; lq<=ESELIMSG(q) && ret>0; lq++) {
		if (genstep_seqsigncheck(q,lq)<0)  continue;
		if (elimchm_signbitscomp(q,lq)<0)  continue;
		qn = elimseqm_setfor_noch(qn,ee,NULL,lq);
		ret = gcanon_dispseqtest(0,q,qn,ee,re,aut)>=0;
	}
	dispose_ematrix(re);
		/**
		 * Here we cycle through all smaller or equal sequence signatures,
		 * through all admissible bases and all displayed base-minors.
		 * Then we call gcanon_dispseqtest() to do the minimality test for each one
		 * displayed minor.
		 * (Recall that base-minor automorphisms are factored out for dispminors.)
		**/
	bas = NULL;  bmlists = NULL;
	if (ran>10 && ran<20)  bas = alist_append(NULL,ematrix_copy(ee));
	else if (ran>0)  bas = struct_minorbases_rn(ee,ran);
	else if (ran<=0)  bas = struct_minorbases(ee);
	DEBUG(CURDLEV+(ran>3)," - seqtest using %d %sbases...\n",alist_getlength(bas),ran>10?"one ":(ran>0?"random ":""));
	t = 1;
	if (bas)  for (lq=0; lq<ESMAXLSG(q) && ret>0; lq++) {
		if (elimchm_signbitscomp(q,lq)<0)  continue;
		if (genstep_seqsigncheck(q,lq)<0)  continue;
		qn = elimseqm_setfor_noch(qn,ee,NULL,lq);
			/* (qn is not a fully valid elim sequence by this setting!) */
		
		DEBUG(CURDLEV+(ran>3)," ~ considering seq signature %lX (length %d)...\n",lq,ESLENGTH(qn));
		nlq++;
		for (b=bas, i=0; (b?*b:0) && ret>0; b++, i++) {
			if (genstep_basischeck(qn,*b)<0)  continue;
			if (t) {	/* (only the first time we compute all dispminors) */
				bmin = NULL;  printlev -= 2;
				k = struct_dispminor(*b,esm,(void**)aut,&bmin);
				printlev += 2;
				if (bmin && k<=0)  {PROGERROR("minors generated even for <=0 ?!?, %d",alist_getlength(bmin));}
				if (k<=0)  bmin = new_alist(8);
				bmlists = alist_append(bmlists,bmin);
			} else  bmin = bmlists[i];
			
			if (!bmin)  {PROGERROR("where is the bmin list ?!?");}
			for (m=bmin; (m?*m:0) && ret>0; m++)
				ret = gcanon_dispseqtest((ran>0?0:1),q,qn,*b,*m,aut)>=0;
			if (t)  nbm += m-bmin;
		}
		t = 0;
		DEBUG(CURDLEV+1,"\t(used %d auts, %d bases, %d dispmins, %d linorders, %d(%d) linprecmps)\n",
				alist_getlength(aut),alist_getlength(bas),nbm,linords,linordp,linordps);
	}
		/**
		 * After all, we clean up local data and check results.
		 * In particular, in the full test ran==0, the sequence identical to given q
		 * must be found somewhere (this is not trivial in the algorithm).
		 * All lists in bmlists and the list bas are properly disposed.
		**/
	if (ret>0 && idfound==0 && ran<=0)  {PROGERROR("The identical sequence was not checked -- non-admissible sequence ???"); ret=0;}
#ifndef FASTPROG
	if (alist_getlength(bmlists)>alist_getlength(bas))  {PROGERROR("Too long bmlists ???");}
	if (aut==minselfmaps) if (!genstep_precheck_isvalid(q,(ISTRANSPM(e)!=ESELIMSG(q)%2),e))  ret = 1;
	DEBUG(CURDLEV-2+(ret>0)+(ran>3),"Testing (%d) seq [%.9s] {ext %.20s} for %s canonical min found %s.\n",
			gener_extensions_ix,ESNAME(q),gener_extprintline(-1,e,1),ran>0?(ran<10?"random":"onedisp"):"all",ret?"PASS":"-NO-");
	SDEBUG(CURDLEV-1+(ran>3),"\t\t\t\t(used %d smaller sigs, %d auts, %d bases, %d dispmins, %d linorders, %d(%d) linprecmps)\n",
			nlq,alist_getlength(aut),alist_getlength(bas),nbm,linords,linordp,linordps);
#endif
	for (bx=bmlists; bx?*bx:0; bx++)	/* (each minor list in bmlists must be freed separately) */
		if (*bx)  dispose_alist_mats(*bx);
	if (bmlists)  alist_free(bmlists);
	if (aut && minselfmaps!=aut) {		/* (frees all the automorphisms if not precomputed in minselfmaps) */
		struct_matrixselfmaps_recycle(aut);  aut = NULL;
	}
	if (ran==0 && bas)  struct_minorbases_recycle(bas);
	else if (bas)  dispose_alist_mats(bas);
	FREE(qn);  dispose_ematrix(ee);
	return (ret>0? 1:-1);
}



/**
 * This function examines all equivalent elim sequences that use the basis bs with the 
 * displayed base-minor mn.
 * Given is qo the original sequence, qn a new sequence with bs as the resulting matrix,
 * bs the resulting matrix of qn equivalent to the matrix of qo, mn a matrix refering
 * to the base-minor in bs (in correct order, but not scaled), aut the list of self-maps
 * of the base-minor.
 * The parameter tt determines level of debug-checks we may afford here - full for tt>0.
 * (The given matrix/sequence may not yet be fully valid during randomized calls...)
 * 
 * The question of this test is to find out whether some other ordering of lines in bs
 * (with the same displayed minor mn and the same signature) is canonically smaller than qo.
 * The return value is 1 if qo passed the test, and -1 otherwise (i.e. smaller seq found).
 * This is a full test - nothing random or incomplete here.
**/

int	gcanon_dispseqtest(int tt, elimseqm *qo, elimseqm *qn, ematrix *bs, ematrix *mn, emlinemap **aut) {
	ematrix		*rem, *small=NULL;
	int		i,j,k, br,qnl,eqsg, ret=1;
	int		xc[GEN_MAXELIMLENGTH+2], sel[GEN_MAXELIMLENGTH+2], *usr=NULL,*usc, ustack[65];
	elimseqm	*qq, *qqs;
	emlinemap	**a, **aur[GEN_MAXELIMLENGTH+2];
	
	if (ESLENGTH(qn)>GEN_MAXELIMLENGTH)  {PROGERROR("cannot work with elim seq longer than %d > %d",ESLENGTH(qn),GEN_MAXELIMLENGTH); return -1;}
#ifndef	FASTPROG
	if (!ISREFMAT(mn) || REFMAT(mn)!=bs || ROWSM(bs)!=ROWSM(ESMATRIX(qn)) || COLSM(bs)!=COLSM(ESMATRIX(qn)))
		{PROGERROR("Call only for a refering matrix  mn -> bs (bs ~ q-matrix) !"); return -1;}
	if (ISTRANSPM(bs) || ISTRANSPM(mn))  {PROGERROR("Call only for a 0-transp matrix bs !"); return -1;}
	if (ESMINOR(qo)!=ESMINOR(qn) || !aut)  {PROGERROR("Call only for the same base minor and given self-maps!"); return -1;}
	if (IFRANDDEBUGLESS(22)) {
		printlev -= 2;
		if (elimchm_signbitscomp(qo,ESELIMSG(qn))<0)  {PROGERROR("Do not call for a larger qn-signature (than qo)!"); return -1;}
		if (aut[0]?strmap_check_size(aut[0],ESMINOR(qo),ESMINOR(qo))<0:1)  {PROGERROR("The (first) line map from aut is wrong!"); return -1;}
		if (EMDATA(mn)?strmap_check_size(EMDATA(mn),ESMINOR(qo),bs)<0:1)  {PROGERROR("The line map of the minor must be given at mn!"); return -1;}
		if (IFRANDDEBUGLESS(111)) if (!struct_isdispminor(ESMINOR(qn),mn))  {PROGERROR("The matrix mn is not (scaled-)equal to the minor in qn!");}
		printlev += 2; }
	DEBUG(CURDLEV+1,"  - testing one basis %p for canonical, base minor %p.\n",bs,mn);
	EMATDEBUGS(CURDLEV+3,bs,"\t\t\t");  EMATDEBUGS(CURDLEV+3,mn,"\t\t\t>\t");
#endif	
		/**
		 * We consider all orders of lines of bs added to the base-minor mn.
		 * For that we use a refering matrix rem that starts from mn and grows to
		 * the whole bs.
		 * (We cannot use bs directly since bs does not display the base minor
		 * in the top-left corner as a sequence matrix should.)
		 * More info in the general algorithm description at the top of this section.
		 * 
		 * The sequence qq is used in finding acceptable orders of added lines,
		 * and qqs is temporarily used when calling gcanon_dispseqtest_scal() later.
		 * The array xc[] contains separate reversed bits of the signature of qn.
		 * The array sel[] keeps the row selections (at level k) in backtracking.
		 * The arrays usr[],usc[] indicate which lines of bs have already been used in rem.
		 * The array aur[] is (possibly) used to refine the set of self-maps aut
		 * to those that still can produce a smaller sequence.
		**/
	qnl = ESLENGTH(qn);
	eqsg = (elimchm_signbitscomp(qo,ESELIMSG(qn))==0);
	qq = elimseqm_copyfor(qn,ESMATRIX(qn));
	qqs = elimseqm_copyfor(qn,ESMATRIX(qn));
	for (i=0; i<qnl; i++)  xc[qnl-1-i] = (ESELIMSG(qn)&(1l<<i))?1:0;
	rem = ematrix_refer_all(mn);
	if (ROWSM(bs)+COLSM(bs)<60)  usr = ustack;
	else  usr = MMALLOC((ROWSM(bs)+COLSM(bs)+5)*sizeof(usr[0]));
	for (i=0; i<ROWSM(bs)+COLSM(bs); i++)  usr[i] = 0;
	usc = usr+ROWSM(bs);
	for (i=0; i<ROWSM(rem); i++)  usr[GETREFMROW(rem,i)] = 1;
	for (i=0; i<COLSM(rem); i++)  usc[GETREFMCOL(rem,i)] = 1;
	ret = 1;
		/**
		 * Here we try to add lines to the base minor mn in all possible orders
		 * respecting the given signature of qn, using sel[] and usr[],usc[].
		 * With each added line, we check for the sequence test.
		 * Finally, we compare the newly created sequence with the original one.
		 * The base minor mn refered in top-left of rem is not touched.
		 * (Notice that we cannot break the while-cycle at k>=0 since we have
		 *  to free some lists collected at each step k,k-1,...)
		**/
	aur[0] = aut;  aur[1] = NULL;  sel[0] = -1;
	k = 0;
	while (k>=0) {		/* (ret<=0 exits the cycle below by steps k-1,k-2...) */
		j = ROWSM(rem)+COLSM(rem)-(ROWSM(mn)+COLSM(mn));
		while (--j>=k) {
			if (xc[j]==0) {		/* (remove previously added extra lines...) */
				usr[GETREFMROW(rem,ROWSM(rem)-1)] = 0;
				ematrix_remove_row(rem,ROWSM(rem)-1);
			} else {
				usc[GETREFMCOL(rem,COLSM(rem)-1)] = 0;
				ematrix_remove_col(rem,COLSM(rem)-1);
			}
			if (aur[j+1] && aur[j+1]!=aur[j])  alist_free(aur[j+1]);
			aur[j+1] = NULL;	/* (j>=0 here!) */
		}
		if (++sel[k]>=(xc[k]?COLSM(bs):ROWSM(bs)) || ret<=0) {
			--k;  continue;		/* (no more choice at this level) */
		}
		if (xc[k]==0? (usr[sel[k]]>0): (usc[sel[k]]>0))
			continue;		/* (the line sel[k] is already used in rem) */
		
		if (xc[k]==0) {		/* add the line sel[k], and try it in the sequence check */
			usr[sel[k]] = 1;  ematrix_drefadd_row(rem,sel[k]);
		} else {
			usc[sel[k]] = 1;  ematrix_drefadd_col(rem,sel[k]);
		}
		aur[k+1] = aur[k];
		qq = elimseqm_setfor_noch(qq,rem,NULL,(ESELIMSG(qn)>>(qnl-k-1)));
			/* (be careful here - rem is only scaled-equal to the base minor in qq!) */
		printlev -= 1;
		br = genstep_sequencecheck(qq,rem,xc[k])<0;
		printlev += 1;
		if (br)  continue;
		DEBUG(CURDLEV+3,"%*s . passed subsequence %d, last %s %d\n",2*k,"",k,xc[k]?"col":"row",sel[k]);
		if (ESLENGTH(qq)!=k+1)  {PROGERROR("Wrong qq-subsequence length %d!=%d.",ESLENGTH(qq),k+1);}
			/**
			 * In the special case that qn has the same signature as qo,
			 * we may also pre-compare subsequences if they are too long, to save time.
			 * (Remember, however, that a smaller signature beats all other comparings!)
			 * If the current subsequence qqs for rem, after all possible maps aut,
			 * is always larger that the relevant part of qo, then rem may never
			 * produce a smaller sequence than qo since the ordering is hereditary.
			 * Moreover, we may discard those self-maps from aut that are larger.
			**/
		if (eqsg && GEN_SUBSEQCOMP>0 && k<qnl-GEN_SUBSEQCOMP) {
			linordp++;  DEBUG(CURDLEV+3,"   subsequence pre-compare for k=%d, qnl=%d\n",k,qnl);
			aur[k+1] = NULL;
			for (a=aur[k]; a?*a:0; a++) {
				j = gcanon_dispseqtest_scal(tt,rem,EMDATA(mn),*a,xc,qq,qqs);
					/* (qqs is set up to the scaled submatrix inside) */
				if (j<=0)  continue;
				if (IFRANDDEBUGLESS(444))  test_elimseqm_notough(qqs);
				br = elimchm_seqcompare_part(qo,qqs,k+1);
						/* (<0 means qo smaller than qqs) */
				dispose_ematrix(ESMATRIX(qqs));
				elimseqm_setfor(qqs,NULL,0);
				if (br>=0)  aur[k+1] = alist_append(aur[k+1],*a);
			}
			if (!aur[k+1]) { linordps++; DEBUG(CURDLEV+2,"  - no smaller sequence possible by a pre-compare at k=%d\n",k); }
			if (!aur[k+1])  continue;
		}
			/**
			 * Here we proceed to the next level with an initial choice of line.
			 * Normally, the initial choice starts with -1 - all possibilities.
			 * However, to speed up the search, if adding lines of the same type
			 * is unordered (reported from genstep_sequencecheck_unord()), then
			 * we start next with the current line choice.
			**/
		if (k<qnl-1) {
			++k;
			if (GEN_SUBSEQUNORD>0 && !eqsg) {
				br = (genstep_sequencecheck_unord(qq,rem,xc[k-1])>0);
				sel[k] = ((br && xc[k]==xc[k-1])? sel[k-1]: -1);
			} else  sel[k] = -1;
			continue;
		}
			/**
			 * Above, we have found all the lines sel[k] at levels k, for all levels.
			 * Now the matrix rem refers to whole bs in the order chosen by sel[].
			 * If the signature of the new qqs is smaller than that of qo, we are
			 * actually done in the first cycle iteration here.
			 * If the signatures are equal, we have to try all applications of
			 * self-maps of the base minor (**aut) to the matrix rem, each
			 * properly scaled to unit-form inside gcanon_dispseqtest_scal().
			 * (The scaled copy of rem is stored directly into qqs there.)
			 * Then we compare the matrices lexicographically (see pre-comparing above).
			 * Finally, if the new sequence is smaller, then we reject the original one.
			**/
		linords++;
		for (a=aur[k]; (a?*a:0) && ret>0; a++) {
			j = gcanon_dispseqtest_scal(tt,rem,EMDATA(mn),*a,xc,qq,qqs);
			if (j<=0)  continue;		/* (actually, j<=0 should never happen below) */
			br = (eqsg? elimchm_seqcompare(qo,qqs): 1);
			if (br==0)  idfound = 1;	/* (found the equal sequence to qo, OK) */
			if (IFRANDDEBUGLESS(br>0?55:444))  test_elimseqm_notough(qqs);
			if (!eqsg || br>0) {
				ret = -1;		/* (>0 means qqs smaller than qo -> reject qo) */
				small = ESMATRIX(qqs);
			} else  dispose_ematrix(ESMATRIX(qqs));
			elimseqm_setfor(qqs,NULL,0);
		}
	}
#ifndef FASTPROG
	if (ret<0 && small) {
		DEBUG(CURDLEV-1,"  - a smaller elim sequence found using sig %lX (~%lX):\n",ESELIMSG(qn),ESELIMSG(qo));
		EMATDEBUGS(CURDLEV,small,"\t\t^^\t");
		ematrix_transpose_set(small,ISTRANSPM(ESMATRIX(qo)));
		if (tt>0 && IFRANDDEBUGLESS(333))  if (!struct_isequivalent(small,ESMATRIX(qo)))
			{PROGERROR("The matrix of smaller seq is not equivalent to the original one!"); EMATDEBUGS(0,bs,"\t\t");EMATDEBUGS(0,small,"\t\t");EMATDEBUGS(0,ESMATRIX(qo),"\to\t");}
	}
	DEBUG(CURDLEV-1+2*(ret>0),"  - testing one basis %p for canonical minimality found %s\n",bs,ret>0?"pass":"NO");
#endif
	dispose_ematrix(rem);
	if (small)  dispose_ematrix(small);
	FREE(qq);  FREE(qqs);
	if (usr && usr!=ustack)  FREE(usr);
	return ret;
}


/**
 * This function produces a copy of the given matrix e after the scales from am and
 * the map and scales from au are applied, and after unit-scaling the remaining lines.
 * (So that the base-minor of qn is identically displayed in the top-left corner.)
 * Here am is the line map that shows the base minor in the basis bs above (already
 * top-left in e), au is a self-map of the base minor, and
 * xc[] gives the sequence signature separated each bit into one array element from highest bit.
 * Additionally, qn is the new constructed elim sequence (used only for the base minor
 * and signature), and qs is a copy of qn used for temporary data and return value here.
 * 
 * The resulting scaled matrix is returned inside qs as the resulting matrix.
 * The return value is 1 for a successfull scale.
**/

int	gcanon_dispseqtest_scal(int tt, ematrix *e, emlinemap *am, emlinemap *au,
						int xc[], elimseqm *qn, elimseqm *qs) {
	int		i,j,k,l, r, rm,cm,rs,cs;
	ematrix		*emin, *es, *res;
	exp_t		xx;
	sign_t		gg;
	
	if (!qs)  {PROGERROR("must have qs (a sequence copy) given in advance here"); return -1;}
	emin = ESMINOR(qn);
	rm = ROWSM(emin);  cm = COLSM(emin);
	es = ematrix_copy(e);
	rs = ROWSM(es);  cs = COLSM(es);
	r = 1;
		/**
		 * Here we copy the entries of e to es according to the map au.
		 * Then we scale the lines intersecting the base minor by both am,au.
		 * (Recall that line scale are indexed by the destination matrix in maps!)
		**/
	for (i=0; i<rs; i++)  for (j=0; j<cs; j++) {
		COPYEXSIGM(es,(i<rm?EMLMPR(au,emin,i):i),(j<cm?EMLMPC(au,emin,j):j), e,i,j);
		if (i>=rm && j>=cm)  break;
	}
	for (i=0; i<rm; i++) {
		k = EMLMPR(am,emin,i);  l = EMLMPR(au,emin,i);
		pfield_mul(-1,EMLMXR(am,emin,k),EMLMGR(am,emin,k),
				1,EMLMXR(au,emin,l),EMLMGR(au,emin,l), &xx,&gg);
		ematrix_multiply_row(es,l, 1,xx,gg);
		//DEBUG(0,"rmul  %s (%s.%s)\n",pfield_pvalue(8,xx,gg),pfield_pvalue(8,EMLMXR(am,emin,k),EMLMGR(am,emin,k)),pfield_pvalue(8,EMLMXR(au,emin,l),EMLMGR(au,emin,l)));
	}
	for (j=0; j<cm; j++) {
		k = EMLMPC(am,emin,j);  l = EMLMPC(au,emin,j);
		pfield_mul(-1,EMLMXC(am,emin,k),EMLMGC(am,emin,k),
				1,EMLMXC(au,emin,l),EMLMGC(au,emin,l), &xx,&gg);
		ematrix_multiply_col(es,l, 1,xx,gg);
	}
	res = ematrix_refer(es,0,rm,0,cm);
	if (IFRANDDEBUGLESS(111)) if (!ematrix_isequal(res,emin))
		{PROGERROR("Scaled es is not equal to the base minor, why ?!?"); EMATDEBUG(0,es,"\t"); EMATDEBUG(0,emin,"\t>\t");}
		/**
		 * Then we unit-scale the remaining lines in the matrix es through the
		 * refering matrix res (res refers to es in the same order of lines).
		 * If everything goes OK (which should always be), then the scaled matrix es
		 * is set as the resulting one in the given sequence qs.
		**/
	for (k=0, i=rm,j=cm; r>0 && (i<rs||j<cs); k++) {
		if (xc[k]==0)  ematrix_refadd_row(res,i++);
		else  ematrix_refadd_col(res,j++);
		elimseqm_setfor(qs,res, (ESELIMSG(qn)>>(ESLENGTH(qn)-k-1)));
		r = genstep_linecheck_setref(qs,res,xc[k])>=0;
			/* (it should not happen that setref fails, but the test is here for safety) */
	}
	if (r<=0 || i>rs || j>cs || k>ESLENGTH(qn))  {PROGERROR("Whole matrix was not unit-scaled, why??? r=%d, i=%d>%d, j=%d>%d, k=%d>%d",r,i,rs,j,cs,k,ESLENGTH(qn)); r = -1;}
	dispose_ematrix(res);
	if (r>0)  elimseqm_setfor(qs,es,ESELIMSG(qn));
	else { dispose_ematrix(es);  elimseqm_setfor(qs,NULL,0); }
	return r;
}




































