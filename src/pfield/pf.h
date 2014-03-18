
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
 * This header is for "private" definitions of the pfield computation.
 * See the main include "pfield.h".
**/


/*#include "pfield.h"*/

#ifndef PF_H
#define	PF_H






/******************	Definitions of partial fields	*******************/
/**************************************************************************/



/**
 * The structure describing one particular partial field as given in the source code:
 * (This information must be precompiled before using in the program computation
 * since it is optimized for easy human input of new (p)fields, not for computing.)
 *  - There are up to PF_NUMNAMES optional names for this pfield.
 *  - Next the chracteristic is stored (-1 for undefined).
 *  - The number of generators (i.e. of actually used exponents) numexp.
 *  - A flag ispart whether the pfield is partial or not.
 *  - The modulos defining multiplicative subgroups of generators, Z_m or Z for m=0.
 *  - The value mnone of the standard form -1.
 *  - An array of all (up to PF_NUMFUNDAMENTALS) fundamental elements excluding 0;
 *    the next element after each fundamental one is the value-1.
 *  - The numeric bases for generators, and their numeric value modulos (usable for Z_p).
 *  - The basic string representatons of generators in exstr[], no longer than 9(!) chars.
 *    Only the first one may be NULL (then exbase and valmod are used to scan numeric input),
 *    all others must be given.
 *  - Optionally, one may provide extra strings reprersenting gen^i (in exstrp[i-1])
 *    or g^-i (in exstrn[i-1]), which are then used when printing depending on .
 *    Again optionally, one may give strings for negative values -gen and -gen^i,-gen^-i
 *    right after the strings for the positive values (like to give "1-a" for -(a-1)).
 *  - Flags telling how to print numbers at various levels of formality (-1,0,1,?),
 *    and how to scan the informal strings (for this see also PFNOATOMINPUT).
 *  
 *  Several important notes about the strings representing generators and their informal powers:
 *  The generator string is always scanned from the input, and so we suggest to enclose
 *  it in square brackets if it looks like an expression -- for example "[a-1]".
 *  If an expression-atom is used anywhere, it should be correctly parsed as an expression.
 *  The optional strings may be strange symbols (like "ai" for a-1), or also expressions.
 *  However, optional informal strings are not scanned from the input if they contain
 *  characters from PFNOATOMINPUT (like "+-*^/") and noatominp>0 is given.
 * 
**/

			/* max number of fundamental elements we can store for one pfield */
#define PF_NUMFUNDAMENTALS	(20+2*PF_NUMEXPS*PF_NUMEXPS*PF_NUMEXPS)
			
#define	PF_NUMNAMES	10	/* max number of string names that one pfield may have */

#define	PF_NUMSPECPOW	10	/* max number of special power-atoms (a^i) one may define */

typedef struct {
	
	char		*pfnames[PF_NUMNAMES];	/* various pfield names */
	char		*pfdescr;		/* short description of the pfield */
	int		characteristic;		/* pfield characteristic, -1 for undef */
	
	int		numexp, ispart;		/* number of generators, whether partial */
	int		exmod[PF_NUMEXPS+1];	/* modulos for exponents of generators */
	
	struct { sign_t	sg;  exp_t ex; }
			mnone,			/* -1 value in standard form (or -1) */
			funds[2*PF_NUMFUNDAMENTALS];	/* fundamental elements and their's-1 (not including 0 !) */
	
	int		exbase[PF_NUMEXPS+1];	/* numeric values of generators - for numeric scanning or hashvalue */
	int		valmod[PF_NUMEXPS+1];	/* modulos for numeric generator values (optional) */
	
	char		*exstr[2*PF_NUMEXPS+2],	/* strings representing the generators (numeric-only if NULL) */
			*exstrp[PF_NUMSPECPOW][2*PF_NUMEXPS+2],	/* optional +power strings */
			*exstrm[PF_NUMSPECPOW][2*PF_NUMEXPS+2];	/* optional -power strings */
						/* - the strings for negative values follow after the positive values */
	
	int		optstrf,		/* optional strings are printed in formats <=optstrf */
			nobracef,		/* no brackets around generators in formats <=nobracef */
			noatominp;		/* >0 means not to input optional strings with chars in PFNOATOMINPUT */
	
}	pf_definition;

extern pf_definition	pfdefins[];

int     pfield_precompile(int k) ;



/**
 * The structure describing one particular partial field in the program computation:
 *  It refers to the original pf_definition,
 *  and it provides extra hashtable of fundamental elements to speed-up computation.
 *  Moreover, a random hash seed is chosen when compiling the pfield to be used
 *  in later computing of hashvalues for pfield elements.
 * 
 * The precompiled definition is kept private in pfield.c, there is no outside access to it.
 * All other program functions should call some of the pfield functions declared in pfield.h.
 * 
**/

#define PF_FUNDHASHSIZE	(8*PF_NUMFUNDAMENTALS+3)	/* the size of the hash table for fundamental elements */

typedef struct {
	
	char		**pfnames;
	pf_definition	pfdef;		/* refers to the original pf_definition */
	
	int		numfund;	/* the actual number of fundamentals, excluding 0 */
	
	unsigned long	hashseed;	/* the hash seed use when computing hash values of pfield elements */
	
	struct { int val;  sign_t sg; exp_t ex;  sign_t sg1; exp_t ex1; }
			fundhash[PF_FUNDHASHSIZE];	/* the hashtable of fundamental elements */
	
}	pf_computing;


extern pf_computing	pfcomputing[], *pfcnow;
extern int		pfnumexp, pfindex, pfispartial;

void    pfield_getfundamental(int ix, exp_t *x, sign_t *g) ;



/**
 * The structure describing translations from another pfield (identified by the translation
 * name) to the current field (which must be recognized by the translation, or "" for any).
 * The field ptnumexp tells the number of exponents in the pfield "from".
 * The array ptimages[] gives the generator images for these exponents.
 * 
 * The translations are defined in pfdefs.c (with additions in pftran-more.inc).
 * They are handled by functions in pfmore.c .
**/

typedef struct {

	char		*tonames[PF_NUMNAMES];	/* names of the "to" pfields (to which we may translate) */
	char		*ptnames[PF_NUMNAMES];	/* various translation names */
	char		*ptdescr;		/* short description of the translation */
	
	int		ptnumexp;		/* the number of exponents to translate */
	struct { sign_t sg;  exp_t ex; }	
			ptimages[PF_NUMEXPS+1];	/* the images of the source generators */

}	pf_translation;

extern pf_translation	pftrans[];











#endif	/* (of #ifndef PF_H) */




























