
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













/******************************************************************************************/
/******************	Generating (?-connected) represented matroids	*******************/
/******************************************************************************************/



/**
 * Some theory about generating matroids for start...
 * 
 * Our task is to "generate" all matroids having the given matroid (called further
 * the "base minor") as a minor, subject to representability in the pfield and to other
 * conditions.
 * In fact, we generate all non-equivalent representations of all such matroids extending
 * the given base minor representation, but that means no difference when the base minor
 * (strongly) stabilizes the representations.
 * No extension generating functions are implemented for abstract matroids here.
 * 
 * In default interpretation,
 * we require the base minor and all the generated matroids to be 3-connected.
 * Then, by Seymour's splitter theorem, there is a sequence of single steps (extensions /
 * co-extensions) building the matroid from the base minor keeping 3-connectivity,
 * except the case when the base minor is a wheel or a whirl.
 * (The assumption of connectivity or 3-connectivity is important in the algorithms!)
 * From version 1.2.09, we allow to construct also disconnected matroids - simple or
 * cosimple, which generalizes some concepts of the generation scheme.
 * 
 * 
 * Such a sequence of extensions/co-extensions (or conversly of deletions/contractions) is
 * called the "elimination sequence" for the (resulting) matroid M over the base minor S.
 * Formally, an elimination sequence consists of the base minor S, of the resulting matroid M,
 * the order of lines of M as they are added to S, and the "signature" of the sequence telling
 * us which lines of M are extended and which are co-extended.
 * (The sequence signature is taken separately from the order since the order naturaly follows
 * from the matrix representation of M, while the signature does not.)
 * The base minor S is displayed in the matrix representation of M.
 * Read in gen.h and elimseq.c about the elimination sequence implementation.
 * 
 * Another important note concerns the words "non-equivalent representations".
 * We use them in the sense of unlabelled strong matrix equivalence, as discussed
 * in include/struct.h.
 * In particular, one abstract matroid may be generated in several non-equivalent
 * representations in pfields that have non-equivalent representations (above GF3).
 * 
 * The whole extension generating algorithm is theoretically described in the paper
 * [Petr Hlineny: Equivalence-free exhaustive generation of represented matroids].
 * Extensive comments are also contained in the Macek manual.
 * 
**/



/******************	Admissible sequences; Canonical sequence	*******************/
/******************************************************************************************/


/**
 * An elimination sequence may have several more attributes controlling what matrices are
 * generated from the sequence and how (including restrictions on steps of the sequence).
 * Examples include forbidden minors, tight majors, no large fans along, etc.
 * (See the note on the functions genstep_....() described below.)
 * Moreover, each generated matrix is automatically stored in a unit-leading scaled form.
 * (This concerns only the extra lines over the base minor S, not S itself.)
 * 
 * These restrictions define what sequences are "admissible" in generating.
 * (It is often much faster to test admissibility right through the generating
 * process than filter-out wrong sequences later...)
 * Each generated matrix is also produced "canonically minimal", as described next.
 * An extensive description of the tests applied in generating is in genstep.c.
 * 
 * 
 * We call two sequences equivalent if their base minors are identical, they have the same
 * additional attributes, and their resulting matrices are (strongly) equivalent.
 * To avoid generating equivalent sequences many times, we require each eliminations sequence
 * we generate to be minimal with respect to the "canonical order":
 * We compare two sequences first by their signatures (prefering the signature bits corresponding
 * to lines closer to the base minor), and then by their lines as they are added in the sequence
 * (again prefering the lines closer to the base minor).
 * Recall that the two matrices share the same base minor representation, and that their added
 * lines are scaled to unit-leading form, so this is a well defined order.
 * Read in genstep.c about how the canonical order is implemented.
 * 
 * Of course, for every ?-connected representable matroid M having S as a minor and satisfiyng
 * additional conditions which we require in the sequence,
 * there is a canonically-minimal elimination sequence leading to M --
 * just take the smallest one among those leading to M (if any one exists).
 * 
 * To implement the canonical test in the generating process correctly, we have to separate
 * and group all tests used when generating matrices according to their effect on the matrices.
 * It is crucial that exactly the same tests are carried out when extending a matrix as well as
 * when checking a sequence minimality, otherwise some extensions may be lost.
 * This "grouping" of tests is achieved in the functions genstep_....() described in gen.h,
 * and all tests added in the future must(!) be incorporated within these functions as well.
 * (Look at the functions ..._incl() at the end of gen.h, and into gener-more.inc.)
 * 
**/



#ifndef GENER_H 
#define GENER_H 














/******************	Extending a matroid	*******************/
/******************************************************************/


/**
 * This function generates column/row (by tr=1/0) extensions to the matrix ex.
 * All distinct extension lines are generated here, normalized to unit-leading form.
 * When extending over a partial field, then also some pfield-checks of the matrix
 * are used to restrict the number of generated choices.
 * (But the extension is not guarranteed to be properly pfield-represented unless tst>=1!)
 * 
 * The generated extensions are collected in the list *els (if given).
 * The function returns the number of generated extensions, or -1 on an error.
 * If ch>1 (>2), then the generated extensions are all checked and printed out.
 * If rnd>0, then only a random sublist of extensions is generated, about rnd choices per each entry.
 * If tst==0, then all generated extensions are tested in a partial field randomly,
 * and if tst==1, then they are fully tested (all aplies only to partial fields).
 * 
 * If ez is given, then the extended line of ex must have the same zero pattern as
 * that line in ez (the matrix ez may be over a different pfield!).
 * Moreover, if ez and eeq>0 are given, then the first eeq entries of the extension
 * line must match the first eeq entries (in row 0) of ez exactly.
**/

int     gener_matextens_ext(int ch, ematrix *ex, int tr, int rnd, int tst,
                                                ematrix *ez, int eeq, ematrix ***els) ;
#define	gener_matextens_col(ex,els)	gener_matextens_ext(0,ex,1,0,-1,NULL,0,els)
#define	gener_matextens_zeros(ex,tr,ez,els)	gener_matextens_ext(0,ex,tr,0,0,ez,-1,els)
#define	gener_matextens_prefix(ex,ez,els)	gener_matextens_ext(0,ex,1,0,0,ez,COLSM(ez),els)

#ifndef GEN_MAXMATEXTENS
	/* the function will not generate more than this number (estimated) of extensions */
#define	GEN_MAXMATEXTENS	1000000
#endif



/**
 * The function gener_extframe() extends the matrix of the given frame fr by a column
 * if trq==1, or by a row if trq==0.
 * The optional parameter ctrl may be used to control the generating process --
 * which test classes are applied in the extension filtering in gener_extensions_pass().
 * (See the GCTRL_* macros below, and in gener.c ...)
 * 
 * The result value is a list of new frames *frout holding the generated matrices,
 * and having no parents or sons.
 * The number of generated frames, or -1, is returned.
 * 
 * The parameteres of the generating process are given by options in the frame fr.
 * The two main options (that must be present right in the frame fr, no inheritance)
 * are GEN_OPTPREFIX"bsize" and GEN_OPTPREFIX"signature".
 * They determine the base minor dimensions (top-left submatrix of fr) and the sequence
 * signature, as defined in ../include/gener.h.
 * These options are also automatically stored to the new frames (signature updated).
 * The new frames also get names determined by the pattern GEN_??EXTNAME, and comments.
 * The sequence signature should never be set by hand(!), but taken from previous computations.
 * 
 * Another important (new from version 1.2.09) parameter is GEN_OPTPREFIX"connect",
 * which change the connectivity requirements on the generated sequence.
 * The values are: 0 cosimple, 1 simple, 2 connected, 3 default 3-connected.
 * String versions of these values are GEN_OPTPREFIX"cosimple", "simple", "connected".
 * 
 * More attributes of the elimination sequence are obtained next from other options.
 * (These options are already inherited, and they are not written to the new frames.)
 * For example, one may give list of forbidden minors, length of forbidden fan, etc.
 * See the description of options in ../frame/fropts.c.
 * It is possible to implement more options for an elimination sequence in the function
 * gener_extframe_incl() in gener-more.inc.
 * 
**/

int     gener_extframe_ext(int ch, framestruc *fr, int trq, int ctrl, framestruc ***frout) ;
#define	gener_extframe_nopr(fr,trq,fout)	gener_extframe_ext(0,fr,trq,0,fout)
#define	gener_extframe(fr,trq,fout)		gener_extframe_ext(1,fr,trq,0,fout)
#define	gener_extframe_print(fr,trq,fout)	gener_extframe_ext(2,fr,trq,0,fout)

#define	GEN_OPTPREFIX	"ext-"		/* common prefix of extension options in the frame */

#define	GEN_EXTNAME	"%s_c%d"	/* name-patterns for extended frames */
#define	GEN_COEXTNAME	"%s_r%d"

					/* macros for use in the ctrl parameter of generating: */
#define	GCTRL_NOPRETEST		64
#define	GCTRL_CANONLEV(l)	((l)&7)
#define	GCTRL_NOSEQTEST		16
#define	GCTRL_NOSTRUCTEST	32
#define	GCTRL_NOSPECTEST	128


/**
 * This function is called from gener_extframe_ext() to do the actual work with
 * generating extensions.
 * All input information are given in the elimination sequence vq (must be void here
 * since we do not export the declaration of elimseqm).
 * The output is returned in the list *qout of sequences with the generated matrices.
**/

int     gener_extensions(int ch, int ctrl, void *vq, int trq, void ***qout) ;









/******************	Generating representations	*******************/
/**************************************************************************/


/**
 * Here we generate (labelled-)distinct unit-scaled representations of the
 * matroid eg (which is over xpfd) over the current pfield.
 * Or, we simply test representability if ch<=1.
 * The representations are returned through the give matrix list *lout.
 * (*lout must be initialized empty or NULL.)
 * 
 * The generated representations are collected in the list *lout (if given, ch>=2).
 * The function returns the number of generated representations, or -1 on an error.
 * If ch==1 or ch>=3, then the results are printed out.
 * If ch>=2, then all the generated representations are collected (and printed out ch>=4).
 * 
 * We provide a separate function for simply testing representability, without
 * generating anything -- supposed to be faster in specific situations...
 * 
**/


int     grepr_generate_ext(int ch, ematrix *eg, int xpfd, ematrix ***lout) ;
#define	grepr_generate(eg,xp,lo)	grepr_generate_ext(0,eg,xp,lo)
#define	grepr_generate_pr(eg,xp,lo)	grepr_generate_ext(1,eg,xp,lo)
#define	grepr_generate_all(eg,xp,lo)	grepr_generate_ext(2,eg,xp,lo)
#define	grepr_generate_prall(eg,xp,lo)	grepr_generate_ext(4,eg,xp,lo)

int     grepr_isrepresented(ematrix *eg, int xpfd) ;














#endif	/* (of #ifndef GENER_H) */





















