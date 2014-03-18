
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
 * This is just a separate part of the local header gen.h...
 * Read also ../include/gener.h.
**/









/******************	Elimination sequence	*******************/
/******************************************************************/


/**
 * The structure elimseqm keeps an elimination sequence from the base minor min to
 * the whole matrix mat, everything stored for the base transposition (tr==0) state.
 * The lines are taken exactly in the order as they apper in the whole matrix mat,
 * starting with the "base minor" min in the top-left corner.
 * The steps of the elimination sequence are encoded by a bit-field where 1 means that a column
 * is extended, and 0 that a row is co-extended (against the base trasposition state).
 * The lines closest to the minor min are taken for the highest bits.
 * This bitfield is stored simply as a binary number, and it is called the "signature" of seq.
 * (See the definition of ESNEXTSG(qs,tr) to see how the signature is constructed.)
 * 
 * Additionally, the sequence structure keeps several parameters of the construction:
 * A list of "forbidden" minors that may not appear in the constructed matroid.
 * A list of minors that define the term "tight major" - for no element x of the matroid of mat,
 *  both contraction and deletion of x may produce a minor from the list tig[].
 * A size of the smallest "fan" that may not appear anywhere in the elim sequence.
 * ...... and more later - some added also through "gener-more.inc".
 * (See the definitions of minors, fans, etc in struct.h ...)
 * 
 * All the parameter-lists in the sequence must be of type alist (as defined in misc.h).
 * These supplementary parameters should not(!) affect the elim sequence itself,
 * only on particular passing-tests.
 * The lists stored in the elimseqm structure are not changed, and they may be common
 * to all sequences derived from the given one.
 * They are never disposed inside elimseqm functions, nor the sequence name is disposed.
 *                         !!!                        !!!
 * One consequence of this is that the lists are not disposed when freeing the sequence,
 * so one must take care of freeing them separately at the end.
 * 
**/


typedef	long	eseqs_t;			/* the type for a sequence signature bitfield */
#define	GEN_MAXELIMLENGTH	(NUM_LONG_BITS-2)	/* maximal length of an elim sequence */

typedef struct {		/* the sequence structure type: */
	
	char		*name;		/* the name for this sequence (updated as needed when generating) */
	
	ematrix		*mat;		/* the current matrix */
	ematrix		*min;		/* the base (starting) minor matrix - indep. of mat */
	eseqs_t		elsg;		/* bitfield of the elim sequence signature */
	
	ematrix		**forb;		/* a NULL-terminated list of forbidden minors */
	ematrix		**tigh;		/* a NULL-terminated list of "tight" minors */
	
	int		nofan;		/* size of a fan which is forbidden, -1 for nothing */
	int		notail;		/* size of a tail which is forbidden, -1 for nothing */
	int		gconn;		/* what value of connectivity we require in generation:
					    (3 3-conn, 2 conn, 1 simple, 0 cosimple, -1 nothing) */
	ematrix		*eprefix;	/* a matrix with one row whose entries start each extension line */

#define	GENER_ELIMQDEF 1
#include "gener-more.inc"		/* third-party additions to the elimseqm structure... */
#undef	GENER_ELIMQDEF
	
	long		sz;		/* (the size of this structure - used for testing correctness) */
	
}	elimseqm;


elimseqm*       new_elimseqm(void) ;
/* (it is enough to "free()" an elim sequence when disposing) */

#define	ESNAME(q)	(q)->name

#define	ESMATRIX(q)	(q)->mat
#define	ESMINOR(q)	(q)->min
#define	ESELIMSG(q)	(q)->elsg

#define	ESFORBID(q)	(q)->forb
#define	ESTIGHT(q)	(q)->tigh
#define	ESCONNECT(q)	(q)->gconn
#define	ESCONNECT_S(q)	(ESCONNECT(q)>2?"3-conn":(ESCONNECT(q)==2?"connect":(ESCONNECT(q)==1?"simple":"cosimple")))
#define	ESNOFAN(q)	(q)->nofan
#define	ESNOTAIL(q)	(q)->notail
#define	ESPREFIXM(q)	(q)->eprefix

#define	ESLENGTH(q)	((ESMATRIX(q)?(ROWSM(ESMATRIX(q))+COLSM(ESMATRIX(q))):0)-\
				(ESMINOR(q)?(ROWSM(ESMINOR(q))+COLSM(ESMINOR(q))):0))
#define	ESNEXTSG(qs,tr)	(2*(qs)+(tr?1:0))	/* (the next sequence signature for step tr) */
#define	ESMAXLSG(q)	(1l<<ESLENGTH(q))	/* (the max length of this sequence signature) */

#define	ESDEFCONN	3
#define	ESREQCONN(q)	((q)->gconn>=2)		/* whether matroids are gen. connected only */
#define	ESREQ3CONN(q)	((q)->gconn>=3)		/* whether matroids are gen. 3-connected only */
#define	ESREQSIMPLE(q)	((q)->gconn==1||(q)->gconn>=3)	/* whether matroids are gen. simple (no loops, no parallel) only */
#define	ESREQCOSIMPLE(q) ((q)->gconn==0||(q)->gconn>=3)	/* whether matroids are gen. cosimple (no coloops, no serial) only */
		/* (Notice that the requirements to be simple or cosimple are crucial also
		   for 3-connectivity of the sequence - that is actually tested at each step!) */


/**
 * These are functions for basic handling of elimination sequences.
 * Each sequence should be initially created by new_elimseqm(), and it may be disposed by free().
 * Functions for simple changes, printing, etc. are provided next.
 * (Find more in elimseq.c ...)
 * Keep in mind that the elim sequence attributes are never disposed inside elimseqm functions,
 * nor the sequence name is freed.
**/

elimseqm*       elimseqm_setfor_ext(int ch, elimseqm *q, ematrix *e, ematrix *emin, eseqs_t qs) ;
#define	elimseqm_setfor(q,e,qs)		elimseqm_setfor_ext(0,q,e,NULL,qs)
#define	elimseqm_setfor_mat(q,e)	elimseqm_setfor(q,e,ESELIMSG(q))
#define	elimseqm_setfor_noch(q,e,em,qs)	elimseqm_setfor_ext(-1,q,e,em,qs)

elimseqm*       elimseqm_copyfor_ext(int sz, elimseqm *q, ematrix *e, eseqs_t qs) ;
#define	elimseqm_copyfor_sig(q,e,qs)	elimseqm_copyfor_ext(1,q,e,qs)
#define	elimseqm_copyfor(q,e)		elimseqm_copyfor_sig(q,e,ESELIMSG(q))

#define	elimseqm_next(q,e,trq)	elimseqm_copyfor_ext(0,q,e,ESNEXTSG(ESELIMSG(q),trq))


void    elimseqm_fprint_ext(FILE *fo, int ch, elimseqm *q, char *pref) ;
#define	elimseqm_fprint_pref(fo,q,pref)		elimseqm_fprint_ext(fo,1,q,pref)
#define	elimseqm_fprint(fo,q)			elimseqm_fprint_pref(fo,q,printoutpref)
#define	elimseqm_fprint_full(fo,q,pref)		elimseqm_fprint_ext(fo,5,q,pref)
#define	elimseqm_fprint_brief(fo,q,pref)	elimseqm_fprint_ext(fo,0,q,pref)

#ifndef	FASTPROG
#define	ELIMQDEBUG(l,q,pf)	(junk = TESTDLEV(l)? (elimseqm_fprint_ext(debugout,2,q,pf),0):0)
#define	ELIMQDEBUGS(l,q,pf)	(junk = TESTDLEV(l)? (elimseqm_fprint_brief(debugout,q,pf),0):0)
#endif
#define ELIMQOUTPUT(q,pf)	elimseqm_fprint_full(printout,q,pf)
#define ELIMQOUTPUTS(q,pf)	elimseqm_fprint_pref(printout,q,pf)



/**
 * This function debug-tests the given elimination sequence for internal consistency,
 * and it prints errors if found.
 * (There is NO testing for the canonical-minimal form or matroidal-minor properties,
 * only for pfield, connectivity, unit-scale, etc.)
 * This function is intended for debug purposes only, not for testing generated matrices!
 * 
 * The version test_elimseqm_notough() performs a lightweight version of the test - it is
 * faster and may be used even when the matrix is not yet known to be valid.
**/

void    test_elimseqm_ext(int ch, elimseqm *q) ;	/* paranoic consistency check - not for generating, only debugging! */
#define	test_elimseqm_ch(ch,q)		test_elimseqm_ext(ch,q)
#define	test_elimseqm(q)		test_elimseqm_ch(2,q)
#define	test_elimseqm_rand(q)		test_elimseqm_ch(1,q)
#define	test_elimseqm_notough(q)	test_elimseqm_ch(0,q)	/* (skips some tests - for use inside generating process) */













/******************	Basic elimination checks	*******************/
/**************************************************************************/


/**
 * The following functions are basic checks used when generating an elimination sequence.
 * (I.e. checks that should be always performed for the generating process to work...)
 * 
 * The first one tests the bits of sequence signature - whether there are right number of
 * extensions and co-extensions (against the base minor).
 * The second one looks whether the given matrix is properly represented (no matter of q).
 * The next one looks whether the specified matrix line is in unit-scaled form (leading 1 in each comp).
 * The next one checks that the specified entries of epr equal to the first entries of the last line.
 * The next two check whether the specified line is a 3-connected (simple/cosimple)
 * and connected (co)extension of the rest.
 * (The current transpose state of the matrix is used, modified by tr (0 row, 1 col).)
 * 
 * All these funtions return 1 on success and <0 on failure.
**/


int     elimchm_signbits_ext(elimseqm *q, eseqs_t qs) ;
#define	elimchm_signbits_poss(q,qs)	elimchm_signbits_ext(q,qs)
#define	elimchm_signbits_check(q)	elimchm_signbits_poss(q,ESELIMSG(q))

#define	elimchm_represented(q,e)		(ematrix_inpfield((e)?(e):ESMATRIX(q))<0? -1:1)
#define	elimchm_represented_last(q,e,tr)	(ematrix_inpfield_last(e,(tr)?0:1,(tr)?1:0)<0? -1:1)
		/* (randomized test looks only at randomly selected subeterminants, plus all 2x2) */
#define	elimchm_represented_rand(q,e)		(ematrix_inpfield_rand((e)?(e):ESMATRIX(q))<0? -1:1)
#define	elimchm_represented_randlast(q,e,tr)	(ematrix_inpfield_randlast(e,(tr)?0:1,(tr)?1:0)<0? -1:1)

int     elimchm_unitline_ext(int set, ematrix *e, int tr, int lo, int st, int acn) ;
#define	elimchm_unitlastline(set,e,tr,acn)	elimchm_unitline_ext(set,e,tr,-1,0,acn)
#define	elimchm_unitlastcol_check(e,acn)	elimchm_unitlastline(0,e,1,acn)
#define	elimchm_unitlastline_set(e,tr,acn)	elimchm_unitlastline(1,e,tr,acn)
#define	elimchm_unitlastcol_set(e,acn)		elimchm_unitlastline_set(e,1,acn)
#define	elimchm_unitlastline_setref(e,tr)	elimchm_unitlastline(2,e,tr)

int     elimchm_prefixline_ext(ematrix *e, int tr, int lo, ematrix *epr, int len) ;
#define	elimchm_prefixlastline(e,tr,ep)	elimchm_prefixline_ext(e,tr,-1,ep,-1)

int     elimchm_3connline_ext(ematrix *e, int tr, int lo) ;
#define	elimchm_3connlastline(e,tr)	elimchm_3connline_ext(e,tr,-1)
#define	elimchm_3connlastcol(e)		elimchm_3connlastline(e,1)

#define elimchm_connline_ext(e,tr,lo)	(ematrix_linezero(e,tr,lo)>0? -1:1)
#define	elimchm_connlastline(e,tr)	elimchm_connline_ext(e,tr,(tr?COLSM(e):ROWSM(e))-1)
#define	elimchm_connlastcol(e)		elimchm_connlastline(e,1)

//******changed******	int     elimchm_3connsequence(elimseqm *q) ;
int     elimchm_xconnsequence(elimseqm *q, int cn) ;


/* testing for fans, tight major, forbidden minors, etc. is implemented more generaly in struct.h */

#define	struct_nofanall(e,fan)			(struct_hasfan_full(e,(fan))? -1:1)
#define	struct_nofanlastline(e,tr,fan)		(struct_hasfan(e,tr,-1,(fan))? -1:1)
	/* remember that when asking for "lastline" fans, the matrix must be 3-connected the step before! */
#define	struct_nofanlastcol(e,fan)		struct_nofanlastline(e,1,fan)
	/* special ch<-3 -> not a full test since the matrix may not be represented, so cannot compute rank */
#define	struct_nofanlastline_nofull(e,tr,fan)	(struct_hasfan_ch(-4,e,tr,-1,fan)? -1:1)




/**
 * Here we provide functions implementing our canonical order tests for elimination sequences.
 * The ordering is described on the top of this file gener.h and detailed in elimseq.c ...
 * (Remember that the sequence signatures are compared first, and only then the matrix lines
 * in their natural order and possibly other resulting-matrix properties.)
 * 
 * The first function compares the sequence signature of the given seq q against another
 * given seq signature (comparing as binary number for now, but that may change in future).
 * The secod function compares compatible parts of two matrices lexicoraphically entry by
 * entry, using the macro EM_COMPARE().
 * It is important to have a special separate comparism function for use in extension generating
 * to properly define the canonical order of sequences.
 * 
 * The third function compares two elimination sequences of the same dimensions and signatures
 * using line-comparing for the added lines along the sequence.
 * The compared lines of sequences must be unit-scaled as given by genstep_linecheck_...()!
 * Moreover, the lines are (and must be!) compared lexicographically with the first added lines
 * compared first.
 * (In this way, sequence-comparing is hereditary on subsequences...)
 * A version comparing only a subsequence prefix (len>0) is provided as well.
 * 
 * The return value is 0 when equal, -1 when the first is smaller, and 1 when the second is
 * smaller (like if we subtracted them...).
**/


int     elimchm_signbitscomp_ext(elimseqm *q, eseqs_t qs) ;
#define	elimchm_signbitscomp(q,qs)	elimchm_signbitscomp_ext(q,qs)


#define	EM_COMPARE(xa,ga,xb,gb)	(((ga)==0&&(gb)==0)?0: ((ga)!=(gb)? (ga)-(gb): pfield_compareexp(xa,xb)))

int     elimchm_compare_ext(int ch, ematrix *e1, int r1, int c1, ematrix *e2, int r2, int c2) ;
#define	elimchm_compare_entry(e1,r1,c1,e2,r2,c2)	elimchm_compare_ext(1,e1,r1,c1,e2,r2,c2)
#define	elimchm_compare_row(e1,r1,e2,r2)	elimchm_compare_ext(1,e1,r1,-1,e2,r2,-1)
#define	elimchm_compare_col(e1,c1,e2,c2)	elimchm_compare_ext(1,e1,-1,c1,e2,-1,c2)
#define	elimchm_compare_all(e1,e2)		elimchm_compare_ext(1,e1,-1,-1,e2,-1,-1)


int     elimchm_seqcompare_ext(int ch, elimseqm *q1, elimseqm *q2, int len) ;
//#define	elimchm_seqcompare_ch(ch,q1,q2)	elimchm_seqcompare_ext(ch,q1,q2,-1)
#define	elimchm_seqcompare(q1,q2)	elimchm_seqcompare_ext(0,q1,q2,-1)
		/* (only a subsequence of length len is compared when len>=0) */
#define	elimchm_seqcompare_part(q1,q2,len)	elimchm_seqcompare_ext(0,q1,q2,len)
































