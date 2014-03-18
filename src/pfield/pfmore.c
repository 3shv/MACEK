
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
 * This file contains additional functions for handling partial fields in the program.
 * In particular, this includes extensions to matrix entries, pfield endomorphisms,
 * and pfield translations.
 * The functions for switching between pfields and between translations are here.
 * 
 * Read about the partial fields in ../include/pfield.h .
**/



#include "macek.h"
#include "pf.h"

#include "pfdefs.c"	/* (the pfield definitions must be included here...) */






/******************	(P)Field definition	*******************/
/******************************************************************/


#define CURDLEV		7


/**
 * Here we keep a "precompiled" pfield description.
 * See the file pf.h for declarations of the structures pf_definition and pf_computing.
 * Generally, pf_computing contains all information as pf_definition, but
 * some of them are better organized for computing (like a hashtable of
 * fundamental elements).
 * 
 * Besides that, we keep the current number of used exponents in a global variable pfnumexp.
**/

pf_computing	pfcomputing[sizeof(pfdefins)/sizeof(pfdefins[0])],
		*pfcnow = NULL;

int		pfnumexp = 0,
		pfindex = -1,
		pfispartial = 0;


/**
 * Printing simple online help on the partial fields in the program.
 * 
**/

void	pfield_pfieldhelp(int ch, FILE *fo) {
	int	i,j;
	
	fprintf(fo,"Short help on %s %s partial fields: (not case sensitive)\n",PROGNAME,PROGVER);
	for (i=0; i<(int)(sizeof(pfdefins)/sizeof(pfdefins[0])); i++) {
		if (!pfdefins[i].pfnames[0])  continue;
		fprintf(fo,"  %s \t[",pfdefins[i].pfnames[0]);
		for (j=1; j<PF_NUMNAMES; j++)
			if (pfdefins[i].pfnames[j])  fprintf(fo," %s",pfdefins[i].pfnames[j]);
		fprintf(fo," ]\n");
		if (ch>0)  fprintf(fo,"\t (%s)\n",pfdefins[i].pfdescr);
	}
	fprintf(fo,"\nShort help on partial field translations:\n");
	for (i=0; i<(int)(sizeof(pftrans)/sizeof(pftrans[0])); i++) {
		if (!pftrans[i].ptnames[0])  continue;
		fprintf(fo,"  %s\t[",pftrans[i].ptnames[0]);
		for (j=1; j<PF_NUMNAMES; j++)
			if (pftrans[i].ptnames[j])  fprintf(fo," %s",pftrans[i].ptnames[j]);
		fprintf(fo," ]\t used in [");
		for (j=0; j<PF_NUMNAMES; j++)
			if (pftrans[i].tonames[j])  fprintf(fo," `%s",pftrans[i].tonames[j]);
		fprintf(fo," ]\n");
		if (ch>0)  fprintf(fo,"\t (%s)\n",pftrans[i].ptdescr);
	}
}







/******************	(P)Field selection	*******************/
/******************************************************************/


/**
 * This function switches the program into a new pfield given by the name pfnm.
 * A pfield of the given name is searched through all definitions in pfdefins[],
 * each one may have several optional names.
 * Only definitions that have non-NULL first name are considered.
 * The last exact name match wins, case insensitive.
 * 
 * If no name match is found, and ch>0, then an error is issued.
 * The function returns 0 on success, and -1 on failure.
 * (See also pfield_precompile() when switching to a completely new pfield...)
 * 
**/

int	pfield_switchto_ext(int ch, char *pfnm) {
	int	i,j,k;
	
	if (pfcnow==NULL) {	/* initial call to this function - clears precompiled pfields */
		for (i=0; i<(int)(sizeof(pfcomputing)/sizeof(pfcomputing[0])); i++)
			pfcomputing[i].pfnames = NULL;
	}
	if (!pfnm)  return -1;
	DEBUG(CURDLEV-2,"Switching to new pfield \"%s\"...\n",pfnm);
	
			/* try to find the given name among all names of all definitions */
	for (i=0,k=-1; i<(int)(sizeof(pfdefins)/sizeof(pfdefins[0])); i++)
		for (j=0; j<PF_NUMNAMES; j++) {
			if (pfdefins[i].pfnames[j]!=NULL) {
				if (!strcasecmp(pfdefins[i].pfnames[j],pfnm))
					k = i;
			} else if (j==0)  break;
		}
	if (k<0) {	/* the name was not found, error */
		if (ch>0)  {USERERROR("Cannot find definition for the pfield \"%s\".",pfnm);}
		return -1;
	}
			/* switching to the new (precompiled) pfield definition */
	pfindex = k;
	pfcnow = &pfcomputing[k];
	pfnumexp = pfdefins[k].numexp;
	pfispartial = pfdefins[k].ispart;
	DEBUG(CURDLEV-2,"  OK, new pfield %s (ix %d).\n",pfdefins[k].pfnames[0],k);
	DEBUG(CURDLEV-2,"   \"%s\"\n",pfdefins[k].pfdescr);
	if (pfcnow->pfnames!=NULL)  return 0;
			/* really new - needs precompiling before use */
	return pfield_precompile(k);
}

void	pfield_switchto_fast(int k) {
	
	if (pfcomputing[k].pfnames==NULL) {
		PROGERROR("Cannot switch to #%d which is not pre-compiled!",k);  return;
	}
	pfindex = k;
	pfcnow = &pfcomputing[k];
	pfnumexp = pfdefins[k].numexp;
	pfispartial = pfdefins[k].ispart;
}



/**
 * This function switches the current pfield translation in the program
 * (as used by pfield_translation_ext() in pfield.c) to the new one of name trnm.
 * If the translation description does not list the current pfield as an acceptable
 * one, then an error occurs (only prefix is compared).
 * The return value is the translation index, or -1.
**/

int	curtranslation = -1;

int	pfield_switchtrans_ext(int ch, char *trnm) {
	int	i,j, k,l;
	
			/* try to find the given name among all names of all definitions */
	DEBUG(CURDLEV-2,"Switching to new pf translation \"%s\"...\n",trnm);
	for (i=0,k=-1; i<(int)(sizeof(pftrans)/sizeof(pftrans[0])); i++)
		for (j=0; j<PF_NUMNAMES; j++) {
			if (pftrans[i].ptnames[j]!=NULL) {
				if (!strcasecmp(pftrans[i].ptnames[j],trnm))
					k = i;
			} else if (j==0)  break;
		}
	if (k<0) {	/* the name was not found, error */
		if (ch>0)  {USERERROR("Cannot find definition for the pf translation \"%s\".",trnm);}
		return -1;
	}
			/* checking that the current pfield is accepted by the translation */
	for (j=0,l=-1; j<PF_NUMNAMES; j++) {
		if (pftrans[k].tonames[j]!=NULL) {
			if (!strncasecmp(pftrans[k].tonames[j],pfield_curname(),
					strlen(pftrans[k].tonames[j])))  l = i;
		} else if (j==0)  break;
	}
	if (l<0) {	/* the translation is not allowed in the current pfield, error */
		if (ch>0)  {USERERROR("The pf translation \"%s\" is not allowed in %s.",trnm,pfield_curname());}
		return -1;
	}
	return curtranslation = k;
}









/******************	Making extension matrices	*******************/
/**************************************************************************/


#undef CURDLEV
#define	CURDLEV		7


/**
 * In order to create a matrix extension, we need a function that provides the possible
 * nonzero matrix entries for the new matrix line.
 * There are two general cases - a finite field, where simply a list of all nonzero
 * elements is given;
 * and a partial field, where possible elements for the given sign are determined
 * such that a sum with another given element is defined.
 * 
 * The functions return the number m of possibilities for n==-1, and then they return these
 * possibilities in subsequent calls for n=1,2,...,m (not n==0 !).
 * The possible values of exponent and sign are stored into *xo, *go.
 * For the pfield case, xx,gg is the number that must have defined sum with *xo,*go.
 * The return value is 0 if the exp value was successfully set (not guaranteed), and -1
 * if the setting of index n is not possible.
 * 
**/

int	pfield_getpvalue_ext(int n, exp_t xx, sign_t gg, exp_t *xo, sign_t *go) {
	int	r = 0;
	
#ifndef FASTPROG
	if (n==0 || n>pfield_numfundamentals())  {PROGERROR("invalid value of index n=%d",n); return -1;}
	if (!pfispartial && !pfield_iszeroexp(xx))  {PROGERROR("it has no sense to ask for a summed-value when not partial");}
	if (pfispartial && gg==0)  {PROGERROR("cannot determine 0 + value for the pfield when partial"); return -1;}
#endif	
	if (n<0)  r = pfield_numfundamentals();
	else {
		pfield_getfundamental(n-1,xo,go);
		r = pfield_mul_ch(0,1,xx,-gg,1,*xo,*go,xo,go)<0 ?-1:0;
	}
#ifndef FASTPROG
	if (n<0)  DEBUG(CURDLEV+3,"   - the number of choices (%s + ?) is %d\n",pfield_printvalue(10,xx,gg),r);
	else  DEBUG(CURDLEV+4,"   - choosing (n=%d) value (%s +) : %s\n",n,pfield_printvalue(10,xx,gg),pfield_printvalue(10,*xo,*go));
#endif	
	return r;
}

int	pfield_getavalue_ext(int n, exp_t *xo, sign_t *go) {
	exp_t	xx;
	
	pfield_setzeroexp(&xx);
	return pfield_getpvalue_ext(n,xx,-1,xo,go);
}







/******************	Pfield endomorphisms	*******************/
/******************************************************************/


#undef CURDLEV
#define	CURDLEV		7


/**
 * This function computes selected endomorphisms of the multiplicative subgroup of
 * the actual partial field.
 * This selection of endomorphisms always include "mapping all to 1" and "absolute
 * value" (this is computed by squaring one generator when no negatives exist).
 * And there may be more endomorphisms supplied for more complicated pfields.
 * (There is no desire to compute all endomorphisms, just selected useful ones!)
 * 
 * If id<0 is given, then the number of all endomorphisms is returned as int.
 * Then, id = 0...max-1 returns the particular endomorphism - it changes directly
 * given pfield elements x,g.
 * (The function return value is not much significant then...)
 * Zero input value is never changed.
 * 
**/

int	pfield_endomorph_ext(int id, exp_t *x, sign_t *g) {
	int	i;
	
	if (id>=0) if (!x || !g) { PROGERROR("Must give both x,g.");  return -1; }
	if (id>=0)  pfield_tostandform(x,g);
	if (id>=0 && *g==0)  return 0;	/* zero value is returned unchanged */
	i = 1;
	if (id==0) {		/* "mapping all to 1" endomorphism */
		pfield_setzeroexp(x);  *g = 1;
	}
	if (pfnumexp>=1 && pfield_characteristic()!=2) {
		i++;
		if (id==1) {	/* "absolute value" endomorphism, only if an exponent is used */
			if (pfield_hasnegative()) {
				*g = 1;		/* (this works if the pfield has negative elements) */
			} else if (pfnumexp>=1) {
				x->xe[0] *= 2;	/* (squaring the first exponent instead of the sign) */
			}
		}
	}
		/**
		 * 
		 * 
		**/
	//*********** add more involved endomorphisms if useful...
	// add  i += ... (!!!)
	
	if (id<0)  return i;	/* returns the number of defined endomorphisms */
	else  pfield_tostandform(x,g);
	return 1;
}




































