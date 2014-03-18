
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
 * Some theory about fields and partial fields for start...
 * 
 * A field is taught in every algebra course...
 * A partial field (shortly called here the "pfield") is defined similarly as a field,
 * but the addition operation is not required to be always defined.
 * (An easy example is a multiplicative subgroup of a field, with 0.)
 * In general, a partial field is always formed by a multiplicative subgroup of
 * invertible elements in a ring, with 0.
 * 
 * We consider only those pfields in which the multiplicative group is finitely
 * generated, and the equation x-1=y has only finitely many solutions.
 * (However, the pfield itself may contain infinitely many elements.)
 * 
 * Hence, we represent an element of a partial field in the form of a pair (x,g),
 * where g is a "sign" (one of -1,0,1) and x is an "exponent" (a vector of integers
 * indexed by the generators of the multiplicative group) - types sign_t, exp_t.
 * In this way, an element may have more than one representation, so one must be careful.
 * Important role is played by so called "fundamental" elements, i.e. those x for
 * which x-1 is defined (except 0).
 * Briefly speaking, a partial field is given by the subgroups (Z or Z_m) of the
 * generators of its multiplicative group, by the list of fundamental elements and their
 * x-1 values, and by the proper representation of -1.
 * The complete declaration of pf_definition[] can be found in pf.h.
 * 
 * Multiplication is then simply implemented by summing the exponents.
 * Division uses the same routine with negative exponents, and so division by 0
 * always results in 0, not an error.
 * A sum a+b is computed as ((a/-b)-1)*(-b), using the fundamental elements.
 * This may not be defined!
 * Read more in pfield_precompile(), pfield_minusone(), pfield_tostandform_ext(),
 * and pfield_arithmetic_ext().
 * 
 * 
**/


#ifndef PFIELD_H
#define	PFIELD_H






/******************	General (p)field definitions	*******************/
/**************************************************************************/


/**
 * This defines the (maximal) number of exponents that one pfield element has.
 * (If less exponents are really used, then the remaining ones are equal to zero.)
 * Try to keep this number as small as practical, to speed the program up.
 * Some functions are optimized for use of 1 or 2 exponents.
 * 
 * The next "pfnumexp" is the current number of exponents actually used in exp_t.
 * (Do not rely on pfnumexp unless in special situations...)
 * Use the next two functions to get the current pfield name and index,
 * and get the current number of fundamentals (excluding 0) from pfield_numfundamentals().
 * The function pfield_hasnegative() tells whether the pfield contains "negative elements",
 * i.e. elements that have sign -1 in the standard form.
 * 
**/

#ifdef  BINARYONLY
#define	PF_NUMEXPS	0
#else
#ifndef	PF_NUMEXPS
#define	PF_NUMEXPS	2
#endif
#endif

extern int	pfnumexp,	/* the actual number of exponents used */
		pfispartial;	/* whether the current pfield is partial (sum) */

char*   pfield_curname(void) ;
int     pfield_curindex(void) ;
char*   pfield_ixname(int k) ;
int     pfield_numfundamentals(void) ;
int     pfield_hasnegative(void) ;
int     pfield_characteristic(void) ;

void    pfield_pfieldhelp(int ch, FILE *fo) ;


/**
 * Basic data types for a sign and exponents of the pfield elements.
 * A sign may have values only -1,0,1, while the exponents are theoretically unlimited,
 * practically we do not allow them to be bigger than PFMAXEXPONENT (warning when >PFMAXEXPONENT/4).
**/

#if PF_NUMEXPS==0
typedef	struct { short xe[1]; }			exp_t;
#else
typedef	struct { short xe[PF_NUMEXPS]; }	exp_t;
#endif
#define PF_MAXEXPONENT				(SHRT_MAX/2)

typedef signed char				sign_t;


/**
 * This function switches the program into a new pfield given by the name pfnm.
 * A pfield of the given name is searched through all definitions in pfdefins[],
 * each one may have several optional names.
 * Only definitions that have non-NULL first name are considered.
 * The last exact name match wins, case insensitive.
 * If no name match is found, and ch>0, then an error is issued.
 * The function returns 0 on success, and -1 on failure.
 * 
 * For already used pfields, one may quickly switch between them using pfield_switchto_fast(),
 * giving as a parameter the pfield indices from pfield_curindex();
**/

int     pfield_switchto_ext(int ch, char *pfnm) ;
#define	pfield_switchto(pfnm)	(pfield_switchto_ext(1,pfnm)>=0?1:-1)

void    pfield_switchto_fast(int k) ;


/**
 * The program offers a possibility to "translate" pfield values from one pfield
 * to another one.
 * The function pfield_switchtrans() is provided for switching the program to
 * the given pfield translation by a name (see pftrans[] in pfdefs.c).
 * The function pfield_translation() then translates the given exponents and sign
 * according to this translation - which gives images of the generators of the
 * given exponent to the current pfield.
 * 
 * (See also the function ematrix_import() in ematrixop.c ...)
**/

int     pfield_switchtrans_ext(int ch, char *trnm) ;
#define	pfield_switchtrans(trnm)	(pfield_switchtrans_ext(1,trnm)>=0?1:-1)

int     pfield_translation_ext(int ch, exp_t x1, sign_t g1, exp_t *x, sign_t *g) ;
#define	pfield_translation_ch(ch,x,g)	pfield_translation_ext(ch,*(x),*(g),x,g)
#define	pfield_translation(x,g)		pfield_translation_ch(1,x,g)


/**
 * Variables for use in local macros...
**/

extern exp_t	pfmacexp;
extern sign_t	pfmacsgn;
extern int	pfmacint;







/******************	(P)Field computation	*******************/
/******************************************************************/



/**
 * These simple functions set and test all-zero exponent value, compare exponents
 * (as a difference x1-x2 lexicographically), convert or test the standard form
 * of a pfield element, compare elements, and compute hashvalue for an element
 * of the pfield.
 * The functions are considered low-level, and they should be used only in special
 * situations which require their use.
 * The functions are also very fast.
 * 	
**/

void    pfield_setzeroexp_ext(exp_t *x) ;
int     pfield_iszeroexp_ext(exp_t x) ;
int     pfield_compareexp_ext(exp_t x1, exp_t x2) ;

#if PF_NUMEXPS==0
#define pfield_setzeroexp(x)	((x)->xe[0]=0)
#define pfield_iszeroexp(x)		1
#define pfield_compareexp(x1,x2)	0
#elif PF_NUMEXPS==1
#define pfield_setzeroexp(x)	((x)->xe[0]=0)
#define pfield_iszeroexp(x)	((x).xe[0]==0||pfnumexp<=0)
#define pfield_compareexp(x1,x2)	(pfnumexp<=0? 0: (x1).xe[0]-(x2).xe[0])
#elif PF_NUMEXPS==2
#define pfield_setzeroexp(x)	((x)->xe[0]=(x)->xe[1]=0)
#define pfield_iszeroexp(x)	(((x).xe[0]==0||pfnumexp<=0)&&((x).xe[1]==0||pfnumexp<=1))
#define pfield_compareexp(x1,x2) \
				(pfnumexp<=0? 0: \
				 (pfnumexp<=1||(x1).xe[0]!=(x2).xe[0])? (x1).xe[0]-(x2).xe[0]: \
				  (x1).xe[1]-(x2).xe[1])
#else
#define	pfield_setzeroexp(x)	pfield_setzeroexp_ext(x)
#define	pfield_iszeroexp(x)	pfield_iszeroexp_ext(x)
#define	pfield_compareexp(x1,x2)	pfield_compareexp_ext(x1,x2)
#endif

int     pfield_tostandform_ext(int ch, exp_t *x, sign_t *g) ;
#define	pfield_tostandform(x,g)		pfield_tostandform_ext(0,x,g)

#define	pfield_isstandform(x,g)		(pfield_tostandform_ext(-1,&(x),&(g))>=0?1:0)
int     pfield_isequal_ext(int ch, exp_t x1, sign_t g1, exp_t x2, sign_t g2) ;
#define	pfield_isequal(x1,g1,x2,g2)	pfield_isequal_ext(0,x1,g1,x2,g2)
#define	pfield_isequal_fast(x1,g1,x2,g2)	pfield_isequal_ext(-1,x1,g1,x2,g2)

unsigned long    pfield_hashvalue_ext(exp_t x, sign_t g) ;
#define	pfield_hashvalue(x,g)	pfield_hashvalue_ext(x,g)


/**
 * This function computes selected endomorphism id of the multiplicative subgroup of
 * the actual partial field, modifying directly given x,g.
 * (The index id is in range 0...pfield_endomorph_number()-1.)
 * The selection of endomorphisms always include "mapping all to 1" and "absolute
 * value" (this is computed by squaring one generator when no negatives exist).
 * And there may be more endomorphisms supplied for more complicated pfields.
 * (There is no desire to compute all endomorphisms, just selected useful ones!)
 * Zero input value is never changed.
 * 
 * The macro PF_MAXENDOMORPH limits the number of endomorphisms that we consider
 * in cases when static allocations are necessary (but the actual number of
 * endomorphisms may be higher).
 * 
**/

int     pfield_endomorph_ext(int id, exp_t *x, sign_t *g) ;
#define	pfield_endomorph(id,x,g)	pfield_endomorph_ext(id,x,g)
#define	pfield_endomorph_number()	pfield_endomorph_ext(-1,NULL,NULL)

#define	PF_MAXENDOMORPH		4



/**
 * The function pfield_arithmetic_ext() performs all usual arithmetic operations with
 * the pfield elements, depending on the value of am.
 * The result is stored in the standard form into given x,g pointers.
 * 
 * When multiplying, the signs are simply multiplied and the exponents summed together.
 * The power (sx1,sx2) simply multiplies the exponents.
 * Sum a1+a2 is implemented as ((a1/-a2)-1)*(-a2) using the fundamental elements of the pfield.
 * While the multiplication (even /0) and power are always defined, the sum may not exist.
 * An error is reported when the result is undefined and ch>0.
 * 
 * The functions are (and should be) quite fast even for complicated pfields.
 * The return value is 0 for success, or -1 for undefined result (then x,g=0).
 * 
**/

int	pfield_arithmetic_ext(int am, int ch, int sx1, exp_t x1, sign_t g1,
					int sx2, exp_t x2, sign_t g2, exp_t *x, sign_t *g) ;

#define	pfield_sum_ch(ch,sx1,x1,g1,sx2,x2,g2,x,g)	pfield_arithmetic_ext(0,ch,sx1,x1,g1,sx2,x2,g2,x,g)
#define pfield_sum(x1,g1,x2,g2,x,g)	pfield_sum_ch(1,1,x1,g1,1,x2,g2,x,g)
#define pfield_sum_test(x1,g1,x2,g2)	(pfield_sum_ch(0,1,x1,g1,1,x2,g2,NULL,NULL)>=0?1:-1)

#define	pfield_mul_ch(ch,sx1,x1,g1,sx2,x2,g2,x,g)	pfield_arithmetic_ext(1,ch,sx1,x1,g1,sx2,x2,g2,x,g)
#define pfield_mul(sx1,x1,g1,sx2,x2,g2,x,g)	pfield_mul_ch(1,sx1,x1,g1,sx2,x2,g2,x,g)

#define	pfield_power_ch(ch,pw,x1,g1,x,g)	pfield_arithmetic_ext(2,ch,pw,x1,g1,1,pfmacexp,0,x,g)
#define pfield_power(pw,x1,g1,x,g)		pfield_power_ch(1,pw,x1,g1,x,g)


/**
 * These functions are used to get all possible nonzero values for the number xo,go.
 * For a partial field, one must supply an additional number xx,gg,
 * and then only those values xo,go are returned which have defined sum with xx,gg.
 * For a finite field, all possible non-zero elements are returned.
 * 
 * The first call (n==-1) returns the total number of possibilities.
 * The values are then returned for n=1,2,...,max.
 * Return value <0 means that the n-th value is not possible to set.
**/

int     pfield_getpvalue_ext(int n, exp_t xx, sign_t gg, exp_t *xo, sign_t *go) ;
#define	pfield_getpvalue(n,xx,gg,xo,go)	pfield_getpvalue_ext(n,xx,gg,xo,go)
#define	pfield_getpvalue_num(xx,gg)	pfield_getpvalue(-1,xx,gg,NULL,NULL)

int     pfield_getavalue_ext(int n, exp_t *xo, sign_t *go) ;
#define	pfield_getavalue(n,xo,go)	pfield_getavalue_ext(n,xo,go)
#define	pfield_getavalue_num()		pfield_getavalue(-1,NULL,NULL)



/**
 * These special re-definitions are used to speed-up computation in a specific
 * pfields (here for binary)...
 * 
**/

#ifdef  BINARYONLY
#define	pfield_tostandform(x,g)		(*(g) = (*(g))%2)
#define	pfield_isstandform(x,g)		((g)==0||(g)==1)
#define	pfield_isequal(x1,g1,x2,g2)	(!(g1)==!(g2))
#define	pfield_isequal_fast(x1,g1,x2,g2)	pfield_isequal(x1,g1,x2,g2)

#define	pfield_sum_ch(ch,sx1,x1,g1,sx2,x2,g2,x,g)	((g)? (*(g) = ((g1)+(g2)+2222)%2):0)
#define pfield_sum(x1,g1,x2,g2,x,g)		pfield_sum_ch(1,1,x1,g1,1,x2,g2,x,g)
#define pfield_sum_test(x1,g1,x2,g2)			1
#define	pfield_mul_ch(ch,sx1,x1,g1,sx2,x2,g2,x,g)	((g)? (*(g) = ((g1)*(g2))%2):0)
#define pfield_mul(sx1,x1,g1,sx2,x2,g2,x,g)	pfield_mul_ch(1,sx1,x1,g1,sx2,x2,g2,x,g)
#endif







/******************	(P)Field input/output	*******************/
/******************************************************************/


/**
 * These are basic pfield input/output functions - all printing or reading 
 * in the program should call them.
 * The functions are not too fast, as it is not supposed that they are called often.
 * 
 * The format of pfield elements is directed by values stored in the structure pf_definition
 * about the current pfield.
 * The general format of printing (/scanning) pfield numbers is implicitely coded-in
 * the functions, and also in the lexical scanner for matrix files in emflex.l.
 * (See emflex.l for a discussion of these dependencies.)
 * 
 * It should be possible to implement an arbitrary partial field by a new pf_definition:
 * just give the string representations of the pfield generators, and possibly
 * extra strings for small powers of generators (not necessary).
 * 
 * The value is printed into the given buffer tx, or into an internal buffer if tx=NULL.
 * Not more than maxl characters are printed (including \0).
**/

char*   pfield_printvalue_ext(char *tx, int maxl, exp_t x, sign_t g, int form) ;
#define	pfield_printvalue_formal(tx,maxl,x,g)	pfield_printvalue_ext(tx,maxl,x,g,1)
#define	pfield_printvalue_to(tx,maxl,x,g)	pfield_printvalue_ext(tx,maxl,x,g,0)
#define	pfield_printvalue(maxl,x,g)		pfield_printvalue_to(NULL,maxl,x,g)
		/* (be careful - pfield_pvalue() prints quite loose expressions!) */
#define	pfield_pvalue_to(tx,maxl,x,g)		pfield_printvalue_ext(tx,maxl,x,g,-1)
#define	pfield_pvalue(maxl,x,g)			pfield_pvalue_to(NULL,maxl,x,g)


/**
 * These are special constant strings used when printing (or scanning, see above)
 * partial field values.
 * Be careful - these constants are closely related to the definition of ematrix
 * lexical analyzer in emflex.l, so change them with consideration!
 * The *1 versions are used when printing formally, the * versions when printing informally.
 * 
 * The next is the input scanning function, with given control parameters.
 * The special input characters are taken from the previous defines.
 * Specially, NOATOMINPUT is a list of characters that may not appear in scanned informal
 * strings if pfdefins[].noatominp is set -- not to scan expressions like "1-a" as atoms.
 * (Anyway, we do not use this function in its full power, we only call
 *  pfield_scanvalue_one() from the lexical analyzer.)
 * ......
**/

#define PFPRINTZERO	"o"
#define PFPRINTZERO1	"0"
#define PFPRINTONE	"1"
#define PFPRINTONE1	"1"

#define PFPRINTMINUS	"-"
#define PFPRINTMINUS1	"-"
#define PFPRINTTIMES	"."
#define PFPRINTTIMES1	"*"
#define PFPRINTPOWER	"^"
#define PFPRINTPOWER1	"^"

#define PFNOATOMINPUT	PFPRINTMINUS PFPRINTMINUS1 PFPRINTTIMES PFPRINTTIMES1 PFPRINTPOWER PFPRINTPOWER1


int     pfield_scanvalue_ext(int ch, int n, char *inp, char *ign, char *noat, char *one,
                                char *minus, char *zero, char *pow, exp_t *x, sign_t *g) ;
#define	pfield_scanvalue_one(ch,inp,x,g)	\
		pfield_scanvalue_ext(ch,1,inp, PFPRINTTIMES1 PFPRINTTIMES, PFNOATOMINPUT,\
			PFPRINTONE1 PFPRINTONE, PFPRINTMINUS1 PFPRINTMINUS,\
			PFPRINTZERO1 PFPRINTZERO, PFPRINTPOWER1 PFPRINTPOWER, x,g)
#define	pfield_scanvalue_all(ch,inp,x,g)	\
		pfield_scanvalue_ext(ch,0,inp, PFPRINTTIMES1 PFPRINTTIMES "()", "",\
			PFPRINTONE1 PFPRINTONE, PFPRINTMINUS1 PFPRINTMINUS,\
			PFPRINTZERO1 PFPRINTZERO, PFPRINTPOWER1 PFPRINTPOWER, x,g)









#endif	/* (of #ifndef PFIELD_H) */































