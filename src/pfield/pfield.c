
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
 * indexed by the generators of the multiplicative group).
 * In this way, an element may have more than one representation, so one must be careful.
 * Important role is played by so called "fundamental" elements, i.e. those x for
 * which x-1 is defined (except 0).
 * Briefly speaking, a partial field is given by the subgroups (Z or Z_m) of the
 * generators of its multiplicative group, by the list of fundamental elements and their
 * x-1 values, and by the proper representation of -1.
 * 
 * Multiplication is then simply implemented by summing the exponents.
 * Division uses the same routine with negative exponents, and so division by 0
 * always results in 0, not an error.
 * A sum a+b is computed as ((a/-b)-1)*(-b), using the fundamental elements.
 * This may not be defined!
 * 
 * See more in pfield_precompile(), pfield_minusone(), pfield_tostandform_ext(),
 * and pfield_arithmetic_ext().
 * 
**/





#include "macek.h"
#include "pf.h"


	/* global variables used in macros and elsewhere... (be careful with them!) */
exp_t	pfmacexp = {{0}};
sign_t	pfmacsgn = -2;
int	pfmacint = 0;
int	inrecur = 0;







/******************	Compiling (p)field definition	*******************/
/**************************************************************************/


#define CURDLEV		7


/**
 * The definition of pfcomputing[] (precompiled pfield definitions), and other
 * global variables must be in pfmore.c, since pfdefs.c is required as an include there.
**/


/**
 * Few more easy functions...
 * 
**/

char*	pfield_curname(void) {
	return pfcnow? pfcnow->pfnames[0]:NULL;
}
int	pfield_curindex(void) {
	return pfindex;
}
char*	pfield_ixname(int k) {
	return k>=0? pfdefins[k].pfnames[0]:NULL;
}

int	pfield_numfundamentals(void) {
	return pfcnow? pfcnow->numfund:0;
}
int	pfield_hasnegative(void) {
	return pfcnow? (pfcnow->pfdef.mnone.sg==-1): 0;
}
int	pfield_characteristic(void) {
	return pfcnow? pfcnow->pfdef.characteristic: -1;
}


/**
 * This is a private function that "precompiles" the pfield of given index k into
 * the format of pf_computing for faster computation (stored in pfcnow-> ).
 * It is intended to be called only from pfield_switchto_ext() below.
 * 
 * The return value is 0 for success, and -1 for an error.
**/

int	pfield_precompile(int k) {
	short		a[PF_FUNDHASHSIZE];
	unsigned long	hbs,l;
	int		i,j, b,bs, nfu,m,mm=0;
	exp_t		xx;
	sign_t		gg;
	
	DEBUG(CURDLEV-2,"Pre-compiling a pfield definition (%d) %s ...\n",k,pfdefins[k].pfnames[0]);
	DEBUG(CURDLEV-3,"  \"%s\"\n",pfdefins[k].pfdescr);
	pfcnow->pfdef = pfdefins[k];	/* a reference to the original description */
	pfcnow->pfnames = pfcnow->pfdef.pfnames;
	hbs = bs = PF_NUMFUNDAMENTALS;
		/**
		 * We find a random hash seed that best distributes the fundamental
		 * elements in the hashtable of size PF_FUNDHASHSIZE.
		 * The same seed will be then used in all computations involving a
		 * hashvalue of a pfield element.
		 * (See pfield_hashvalue_ext().)
		**/
	for (j=0; j<30; j++) {
		pfcnow->hashseed = RANDOM();
		for (i=0; i<PF_FUNDHASHSIZE; i++)  a[i] = 0;
		for (i=0; i<PF_NUMFUNDAMENTALS; i++)
		  if (pfcnow->pfdef.funds[2*i].sg!=0) {	/* counts repetitions for each table entry */
			l = pfield_hashvalue(pfcnow->pfdef.funds[2*i].ex,
						pfcnow->pfdef.funds[2*i].sg);
			a[l%PF_FUNDHASHSIZE]++;
		}
		for (i=b=0; i<PF_FUNDHASHSIZE; i++)
			if (a[i]>b)  b = a[i];		/* find the biggest repetition */
		if (b<=bs) { bs = b;  hbs = pfcnow->hashseed; }
		if (bs<=1)  break;	/* (no repetition - very good) */
	}
	DEBUG(CURDLEV-2,"  - negative -1 is %s (%s)\n",pfield_pvalue(14,pfcnow->pfdef.mnone.ex,pfcnow->pfdef.mnone.sg),pfcnow->pfdef.mnone.sg<0?"has -":"no -");
	DEBUG(CURDLEV-2,"  - hashtable of fundamental elements of size %d with repetition %d:\n",PF_FUNDHASHSIZE,bs);
	if (bs>3) {
		DEBUG(0,"## No good hashseed %d found when compiling a new pfield %d!",bs,k);
	}
		/**
		 * Here we permanently set the best hash seed found above,
		 * and we use the derived hashvalues to set up the hashtable of
		 * fundamental elements.
		 * Each entry of the table contains the fundamental f and the value f-1.
		**/
	pfcnow->hashseed = hbs;
	for (i=0; i<PF_FUNDHASHSIZE; i++)	/* (marks invalid table entries) */
		pfcnow->fundhash[i].val = 0;
	for (i=nfu=0; i<PF_NUMFUNDAMENTALS; i++)
	  if (pfcnow->pfdef.funds[2*i].sg!=0) {	/* (only nonzero input fundamentals are valid) */
		nfu++;
		l = pfield_hashvalue(pfcnow->pfdef.funds[2*i].ex,pfcnow->pfdef.funds[2*i].sg)%PF_FUNDHASHSIZE;
		pfcnow->fundhash[l].val = 1;	/* (indicates an occupied table entry) */
		pfcnow->fundhash[l].sg = pfcnow->pfdef.funds[2*i].sg;
		pfcnow->fundhash[l].ex = pfcnow->pfdef.funds[2*i].ex;
		pfcnow->fundhash[l].sg1 = pfcnow->pfdef.funds[2*i+1].sg;
		pfcnow->fundhash[l].ex1 = pfcnow->pfdef.funds[2*i+1].ex;
		SDEBUG(CURDLEV-1,"\tfund(%d@%ld): %s -> %s ;",i,l,pfield_pvalue(14,pfcnow->fundhash[l].ex,pfcnow->fundhash[l].sg),pfield_pvalue(14,pfcnow->fundhash[l].ex1,pfcnow->fundhash[l].sg1));
	}
	SDEBUG(CURDLEV-1,"  (%d total)\n",nfu);
	pfcnow->numfund = nfu;
		/**
		 * These are some tests to see whether the partial field is correctly given.
		 * More tests may be added later............
		 * 	...
		 * 
		**/
	DEBUG(CURDLEV-2,"  - checking correct pfield definition...\n");
	if (pfcnow->pfdef.numexp<0 || pfcnow->pfdef.numexp>PF_NUMEXPS)  {PROGERROR("wrong number of generators %d",pfcnow->pfdef.numexp);}
	if (pfcnow->pfdef.mnone.sg!=-1 && pfcnow->pfdef.mnone.sg!=1)  {PROGERROR("wrong -1 value sign %d",pfcnow->pfdef.mnone.sg);}
	if (pfcnow->pfdef.mnone.sg==-1 && pfcnow->pfdef.characteristic>0)  {PROGERROR("cannot have negative elements for characteristic>0");}
	for (i=0; i<pfnumexp; i++)  if (pfcnow->pfdef.exbase[i]<=0)  {PROGERROR("exbase must be positive %d",pfcnow->pfdef.exbase[i]);}
	for (i=0; i<pfnumexp; i++)  if ((pfcnow->pfdef.exstr[i]?strlen(pfcnow->pfdef.exstr[i]):0)>=10)  {PROGERROR("cannot have generator \'%s\' longer than 10",pfcnow->pfdef.exstr[i]);}
	for (i=1; i<pfnumexp; i++)  if (!pfcnow->pfdef.exstr[i])  {PROGERROR("only the first generator may be numeric");}
	if (!pfcnow->pfdef.ispart)  for (i=0,mm=1; i<pfnumexp; i++) {
		m = pfcnow->pfdef.exmod[i];  mm *= m;
		if (m<1)  {PROGERROR("cannot have unbounded exponent when not partial");}
	}
	if (!pfcnow->pfdef.ispart && pfcnow->numfund!=mm)  {PROGERROR("all elements must be fundamental when not partial");}
	
	//******** test endomorphisms ???	pfield_endomorph_ext()
	
	for (i=0; i<PF_NUMFUNDAMENTALS; i++)  if (pfcnow->pfdef.funds[2*i].sg!=0) {
		if (!pfield_isstandform(pfcnow->pfdef.funds[2*i].ex,pfcnow->pfdef.funds[2*i].sg))  {PROGERROR("fundamental #%d not in standard format",i);}
		pfield_setzeroexp(&xx);  gg = 1;
		pfield_sum_ch(1,-1,pfcnow->pfdef.funds[2*i].ex,pfcnow->pfdef.funds[2*i].sg,1,xx,-gg,&xx,&gg);
		pfield_mul_ch(1,1,pfcnow->pfdef.funds[2*i].ex,pfcnow->pfdef.funds[2*i].sg,1,xx,gg,&xx,&gg);
		pfield_sum_ch(1,1,pfcnow->pfdef.funds[2*i].ex,pfcnow->pfdef.funds[2*i].sg,1,xx,gg,&xx,&gg);
		if (gg!=1 || !pfield_iszeroexp(xx))  {PROGERROR("wrong #%d result of (1/x-1)*x+x != 1 but %s",i,pfield_pvalue(14,xx,gg));}
	}
	if (perror_occured)  {PROGERROR("Failed compilation of the pfield \"%s\", exitting!",pfdefins[k].pfnames[0]); exit(-1);}
	DEBUG(CURDLEV-2,"  ... pfield compilation passed OK.\n");
	return 0;
}



/**
 * This is a private function used in pfield_arithmetic_ext() below to compute sum.
 * The function computes (x,g)-1 value for the given x,g, the result stored back
 * to x,g, or 0 result if undefined.
 * A precompiled hashtable of fundamental elements is used in the computation.
 * (See pfield_precompile() for description of the hashtable.)
 * 
 * The parameter ch affects only debug prints - none if ch<0.
 * The return value (in addition to x,g) is 1 if defined, and -1 if undefined.
**/

int	pfield_minusone(int ch, exp_t *x, sign_t *g) {
	int	h,i,j, ret=0;
	
#ifndef	FASTPROG
	if (!x || !g)  {PROGERROREXIT("Must give both x,g!");}
	if (IFRANDDEBUGLESS(111))  if (!pfield_isstandform(*x,*g)) {
		PROGERROR("The given element must be in the standard form!");}
	DEBUG(CURDLEV+2+(ch<0),"Looking for fundamental elem %s -> ",pfield_pvalue(14,*x,*g));
#endif
	if (*g==0) {		/* returning the 0-1 value: */
		*x = pfcnow->pfdef.mnone.ex;
		*g = pfcnow->pfdef.mnone.sg;
		ret = 1;
	} else {		/* looking for x,g in the table of fundamental elements: */
		h = (int)(pfield_hashvalue(*x,*g)%PF_FUNDHASHSIZE);
		for (i=h,j=0; ret==0 && j<PF_NUMFUNDAMENTALS; j++) {
			if (pfcnow->fundhash[i].val==0)
				ret = -1;		/* no valid match found */
			else if (pfcnow->fundhash[i].sg==*g)
				if (pfield_compareexp(pfcnow->fundhash[i].ex,*x)==0) {
					ret = 1;
					break;		/* matched entry i of the table */
				}
			i = (i+1)%PF_FUNDHASHSIZE;
		}
		if (ret<=0) {
			pfield_setzeroexp(x);  *g = 0;	/* error - not found (x,g)-1, undefined */
		} else {
			*x = pfcnow->fundhash[i].ex1;	/* the value of (x,g)-1 is returned */
			*g = pfcnow->fundhash[i].sg1;
		}
	}
	//SDEBUG(CURDLEV+2+(ch<0),"%s%s\n",ret>0?" ":"\n\t",ret>0?pfield_pvalue(14,*x,*g):"XX fundamental NOT found");
	return ret>0?1:-1;
}


/**
 * This function returns the ix-th fundamental element (as taken from the definition
 * structure, no special order) through x,g.
**/

void	pfield_getfundamental(int ix, exp_t *x, sign_t *g) {
	int	i;
	
	if (ix<0 || ix>=pfield_numfundamentals())  {PROGERROR("fund-index out of bounds %d",ix); return;}
	for (i=0; i<PF_NUMFUNDAMENTALS; i++)
		if (pfcnow->pfdef.funds[2*i].sg!=0)
			if (--ix<0)  break;
	*x = pfcnow->pfdef.funds[2*i].ex;
	*g = pfcnow->pfdef.funds[2*i].sg;
}














/******************	(P)Field arithmetic	*******************/
/******************************************************************/


#undef CURDLEV
#define CURDLEV		7


/**
 * These simple functions set and test all-zero exponent value.
 * (For small number of exponents, macros are used instead, see include/pfield.h .)
**/

void	pfield_setzeroexp_ext(exp_t *x) {
	int	i;
	for (i=0; i<pfnumexp; i++)  x->xe[i] = 0;
}

int	pfield_iszeroexp_ext(exp_t x) {
	int	i;
	for (i=0; i<pfnumexp; i++) 
		if (x.xe[i]!=0)  return 0;
	return 1;
}


/**
 * This simple function returns the first (lexicographic) difference between the given
 * exponent types  x1 - x2.
**/

int	pfield_compareexp_ext(exp_t x1, exp_t x2) {
	int	i;
	for (i=0; i<pfnumexp; i++) 
		if (x1.xe[i]!=x2.xe[i])  return (x1.xe[i]-x2.xe[i]);
	return 0;
}



/**
 * This function turns the given pfield element x,g into the standard form (unique!).
 * What a standard form is determined by pf_definition:
 * A possible negative sign is turned into positive via multiplying by the "-1" value.
 * Then all exponents that are from Z_m (i.e. modulo m) are shifted into the
 * range 0...m-1; exponents from Z are left unchanged.
 * Morevoer, if the sign absolute value is bigger than 1, then that many times is the
 * result summed together (like multiplied by a scalar).
 * 
 * Two pfield elements are equal iff their standard forms are bitwise-equal!
 * If ch==-1, then it only tests for the standard form - returns -1 if not, 0 if yes.
 * 
**/

int	pfield_tostandform_ext(int ch, exp_t *x, sign_t *g) {
	int	i,m,r, gg;
	exp_t	xm;
	sign_t	gm;
	
	if (!x || !g) { PROGERROR("Must give both x,g.");  return -1; }
	gg = *g;  *g = (gg>0?1:(gg<0?-1:0));
				/* checking negative value first */
	if (gg<0 && pfcnow->pfdef.mnone.sg>0) {
		if (ch<0)  return -1;
		*g = 1;
		for (i=0; i<pfnumexp; i++)	/* multiply by a standard repr of -1 */
			x->xe[i] += pfcnow->pfdef.mnone.ex.xe[i];
	}
	if (gg==0) {		/* setting unique representation of zero */
		if (ch<0) if (!pfield_iszeroexp(*x))  return -1;
		pfield_setzeroexp(x);
	}			/* then checking exponents against their modulos */
	if (gg!=0) for (i=0; i<pfnumexp; i++) {
		m =  pfcnow->pfdef.exmod[i];
		if (m<1)  continue;	/* (no modulo here - the whole set Z) */
		if (ch<0) if (x->xe[i]<0 || x->xe[i]>=m)  return -1;
		x->xe[i] = (x->xe[i]>0? x->xe[i]%m: ((1-m)*x->xe[i])%m);
	}
	if (gg<0)  gg = -gg;
	if (gg>1) {		/* summing the result together gg-times if gg>1: */
		if (ch<0)  return -1;
		xm = *x;  gm = *g;
		for (r=0,i=1; i<gg && r>=0; i++)
			r = pfield_sum_ch(0,1,xm,gm,1,*x,*g,x,g);
		if (r<0)  {USERERROR("Failed to make standard form from %d-multiple of %s.",gg,pfield_pvalue(12,xm,gm));  *x = xm; *g = gm;}
	}
	return 0;		/* successful finish */
}


/**
 * This function simply compares the two given pfield elements for pfield equality.
 * For ch<0, the function assumes that the numbers are already in their standard forms.
**/

int	pfield_isequal_ext(int ch, exp_t x1, sign_t g1, exp_t x2, sign_t g2) {
	
#ifndef FASTPROG
	DEBUG(CURDLEV+3,"comparing %s ~ %s\n",pfield_pvalue(12,x1,g1),pfield_pvalue(12,x2,g2));
	if (ch<0 && IFRANDDEBUG(333))
		if (!pfield_isstandform(x1,g1) || !pfield_isstandform(x2,g2))
			{PROGERROR("The compared numbers are supposed to be in their standard forms here!");}
#endif
	if (ch>=0)  pfield_tostandform(&x1,&g1);
	if (ch>=0)  pfield_tostandform(&x2,&g2);
	return (g1==g2 && pfield_compareexp(x1,x2)==0);
}


/**
 * Here a "hash value" for the given pfield element is computed.
 * The basic (necessary!!!) properties of this value are:
 *  - non-negative value,
 *  - the same element (whatever representation) gets always the same hash value
 *    through one program run (see also pfield_tostandform_ext()...),
 *  - the hash value (logairtmically) respects multiplication h(x*y) = h(x)+h(y)
 *    and power h(x^r) = h(x)*r for nonzero(!) elements,
 *  - the value h(0) = 0x10000 (optional).
 * The computation depends on pfcnow->hashseed value (set up when compiling the pfield definition).
 * 
 * The resulting value is returned as a long.
 * 
**/

unsigned long	pfield_hashvalue_ext(exp_t x, sign_t g) {
	unsigned long	h, sd;
	int		i;
	
	pfield_tostandform(&x,&g);
	if (g==0)  return 0x10000l;	/* special handling of zero value */
	
	h = 0l;  sd = pfcnow->hashseed;
	if (g<0) {			/* adding the highest bit for negative values */
		g = 1;  h = 1l<<(NUM_LONG_BITS-1);
	}
	for (i=0; i<pfnumexp; i++) {	/* summing hash values over all exponents */
		h += pfcnow->pfdef.exbase[i]*(sd&0xffffl)*x.xe[i];
		sd = 1003*sd+10003l;	/* (with shifting hash seed) */
	}
	return h;
}



/**
 * This function performs all usual arithmetic operations with the pfield elements.
 * It computes (x1,g1)^sx1 + (x2,g2)^sx2 for am=0 and (x1,g1)^sx1 * (x2,g2)^sx2 for am=1.
 * If am=2, only the power (x1,g1)^sx1 is computed.
 * (Use -g2 for subtraction, or sx2=-1 for division.)
 * The result is stored in the standard form into given x,g pointers.
 * 
 * When multiplying, the signs are simply multiplied and the exponents summed together.
 * The power (sx1,sx2) multiplies the exponents.
 * Sum a1+a2 is implemented as ((a1/-a2)-1)*(-a2) using the fundamental elements of the pfield.
 * (See pfield_minusone() description above.)
 * While the multiplication (even /0) and power are always defined, the sum may not exist.
 * An error is reported when undefined and ch>0.
 * Nothing is printed if ch<0, but the computation is complete.
 * 
 * The return value is 0 for success, or -1 for undefined result (then x,g=0).
 * The function can be used just to verify a defined sum, then give x,g as NULL.
 * 
**/

int	pfield_arithmetic_ext(int am, int ch, int sx1, exp_t x1, sign_t g1,
					int sx2, exp_t x2, sign_t g2, exp_t *x, sign_t *g) {
	int	i, ret = 0;
	
#ifndef	FASTPROG	/* debug - used to recursively test the computation... (only if ch>=0) */
	exp_t	tx1,tx2;	sign_t	tg1,tg2=0;
	tx1 = x1;  tg1 = g1;
	if (ch>=0 && am<2 && IFRANDDEBUGLESS(222)) {
		tx2 = x2;  tg2 = g2;
	} else  tg2 = -2;
#endif
	if (!pfcnow)  {PROGERROREXIT("No pfield is defined for computation!");}
#ifndef	FASTPROG
	DEBUG(CURDLEV+2+(ch<0),"Computing  %s [^%d] %c %s [^%d] in %s.\n",pfield_pvalue(14,x1,g1),sx1,
			am==0?'+':(am==1?'*':' '),am>=2?"":pfield_pvalue(14,x2,g2),sx2,pfcnow->pfnames[0]);
	if (ch>=0 && IFRANDDEBUG(66))  for (i=0; i<pfnumexp; i++)
		if (abs(sx1*x1.xe[i])>PF_MAXEXPONENT/4 || (am>=2?0: abs(sx2*x2.xe[i])>PF_MAXEXPONENT/4)) {
			USERERROR("The exponent is very large (%s, %s), expect arithmetic overflow...",
					pfield_pvalue(16,x1,g1),am>=2?"":pfield_pvalue(16,x2,g2));
			break;
		}	
	if (g1<-100 || g1>100 || (am>=2?0: (g2<-100 || g2>100))) {
		PROGERROR("Invalid sign value in arithmetic %d, %d.",g1,g2);
		return -1; }
#endif
	if (g1<-1 || g1>1)  pfield_tostandform(&x1,&g1);
	if (g2<-1 || g2>1)  pfield_tostandform(&x2,&g2);
	
		/**
		 * First, the powers (sx1,sx2) multiply the given exponents.
		 * Then the selected operation is performed:
		 *  If am=2, only the power (x1,g1)^sx1 is computed.
		 *  If am=1, the signs are multiplied and the exponents summed together.
		 *  If am=0, then the sum is computed.
		 * The sum a1+a2 is implemented as ((a1/-a2)-1)*(-a2) using the fundamental
		 * elements of the pfield (see pfield_minusone() above).
		 * Finally, the result is turned into the standard form.
		**/
	if (sx1!=1 || sx2!=1)  for (i=0; i<pfnumexp; i++) {
		x1.xe[i] *= sx1;
		if (g1<0 && (sx1&1)==0)
			g1 = -g1;
		if (am<2) {
			x2.xe[i] *= sx2;
			if (g2<0 && (sx2&1)==0)  g2 = -g2;
		}
	}
	switch(am) {	/* then perform the requested operation : */
	
	case 2:	ret = 0;	/* (power is already computed by sx1...) */
		break;
	
	case 1:	g1 *= g2;	/* mul: multiply signs, and sum exponents */
		for (i=0; i<pfnumexp; i++)
			x1.xe[i] += x2.xe[i];
		ret = 0;
		break;
	
	case 0:			/* sum: compute (a/-b)-1 by recursive use of multiplication */
		if (g2==0)  break;
		if (g1==0) { x1 = x2;  g1 = g2;  break; }
		ret = pfield_mul_ch(-1,1,x1,g1,-1,x2,-g2,&x1,&g1);
		if (ret>=0)  ret = pfield_minusone(ch,&x1,&g1);
		if (ret>=0)  ret = pfield_mul_ch(-1,1,x1,g1,1,x2,-g2,&x1,&g1);
		if (ret<0 && ch>0)  {PROGERROR("Undefined sum  %s (already modified!) + %s in %s.",
				pfield_pvalue(14,x1,g1),pfield_pvalue(14,x2,g2),pfcnow->pfnames[0]);}
#ifndef	FASTPROG
		if (ret<0 && ch==0)  {DEBUG(CURDLEV-1,"Undefined sum  %s + %s in %s.\n",
				pfield_pvalue(14,tx1,tg1),pfield_pvalue(14,x2,g2),pfcnow->pfnames[0]);}
#endif
		break;
	
	default:	PROGERROREXIT("Undefined operation %d.",am);
	}			/* finally, result to standard pfield form */
	if (ret>=0) {
		if (x || g)  pfield_tostandform(&x1,&g1);
	}
	  else {		/* result 0 if not defined properly... */
		pfield_setzeroexp(&x1);  g1 = 0;
		//DEBUG(CURDLEV+1+2*(ch<0),"  XX NO correct result, setting to 0.\n");
	}
#ifndef	FASTPROG
	if (ret>=0) {		/* debug message, and some paranoic tests... */
		DEBUG(CURDLEV+2+(ch<0),"   OK, result is  %s .\n",pfield_pvalue(14,x1,g1));
		if (tg2>-2 && ch>=0) {
			pfield_arithmetic_ext(am,-1,sx2,tx2,am==0?-tg2:tg2,sx1,tx1,-tg1,&tx2,&tg2);
			pfield_sum_ch(-1,1,tx2,tg2,1,x1,g1,&tx2,&tg2);
			if (tg2!=0)  {PROGERROR("Computation error (not comutative?? ) - difference %s.",pfield_pvalue(14,tx2,tg2));
				DEBUG(0,"\t%s (^%d) %c %s =? %s\n",pfield_pvalue(14,tx1,tg1),sx1,am==0?'+':'*',pfield_pvalue(14,x2,g2),pfield_pvalue(14,x1,g1));}
	}	}
#endif
	if (x)  *x = x1;	/* returning the result through given pointers */
	if (g)  *g = g1;
	return ret;
}


/**
 * This function translates the given number x1,g1 to the current pfield according to
 * the current pfield translation (see pfield_switchtrans()).
 * (The translation defines images of the pfield generators.)
**/

int	pfield_translation_ext(int ch, exp_t x1, sign_t g1, exp_t *x, sign_t *g) {
	extern int	curtranslation;
	pf_translation	*ctr = pftrans+curtranslation;
	int		i;
	exp_t		xx;
	sign_t		gg;
	
	if (curtranslation<0)  {PROGERROR("No pf translation was defined!"); return -1;}
	pfield_tostandform(&x1,&g1);
	pfield_setzeroexp(&xx);  gg = g1;
	for (i=0; i<ctr->ptnumexp; i++) {
		if (x1.xe[i]==0)  continue;
		pfield_tostandform(&ctr->ptimages[i].ex,&ctr->ptimages[i].sg);
		pfield_mul(x1.xe[i],ctr->ptimages[i].ex,ctr->ptimages[i].sg,1,xx,gg,&xx,&gg);
	}
	DEBUG(CURDLEV+1,"Translated from %s to %s.\n",pfield_pvalue(12,x1,g1),pfield_pvalue(12,xx,gg));
	if (x)  *x = xx;  if (g)  *g = gg;
	return 1;
}
















/******************	(P)Field input/output	*******************/
/******************************************************************/


#undef CURDLEV
#define CURDLEV		7	/*7*/


/**
 * This function scans the given string inp for a pfield value (stored through x,g).
 * An "atom" is each generator or a power of it.
 * 
 * The function reads first n atoms from inp (or all inp if n=0).
 * Special characters in inp are handled as follows:
 *  chars from ign are totally ignored (used just as separators), but not after last atom if n>0,
 *  chars from minus invert the scanned value,
 *  chars from zero make the result equal to zero, chars from one do not change it,
 *  chars from pow separate generator string from its power (as +-integer number).
 * 
 * The value of x,g is set to 1 at beginning.
 * Then each scanned generator (or power of it) multiplies the current value.
 * (Powers ^ apply only to their immediately preceeding generators.)
 * No complex expression are allowed.
 * If ch>0, then syntax errors are verbosely reported.
 * 
 * The call
 *  pfield_scanvalue_ext(ch,0,inp,PFPRINTTIMES,PFPRINTMINUS,PFPRINTZERO,PFPRINTPOWER,x,g)
 * should be an inverse to pfield_printvalue_ext().
 * (Be careful with differences between formal/informal printing! - see pfield_printvalue_ext().)
 * However, practically we call this function from the lexical analyzer just to scan
 * one atom from an expression.
 * 
 * The scanned value is stored to x,g in the standard form.
 * The function returns the number of characters that were scanned from inp,
 * or -1,-2 if an error happened.
 * 
**/

int	pfield_scanvalue_ext(int ch, int n, char *inp, char *ign, char *noat, char *one,
				char *minus, char *zero, char *pow, exp_t *x, sign_t *g) {
	int	i,j,jj, m,b,bb,d, l,jm,jjm,lm, ret=0;
	char	*s, *t, *tt;
	exp_t	xx;
	sign_t	gg;
	
	gg = 1;  pfield_setzeroexp(&xx);
	DEBUG(CURDLEV+1,"Going to scan %d atoms from \"%s\"...\n",n,inp);
	if (!pfcnow)  {PROGERROREXIT("No pfield is defined for computation!");}
	
			/* cycles through the given number n of atoms to scan */
	for (i=jjm=0, s=inp; *s && ret>=0; ) {
				/* (i is incremented below after each atom, number or 0/1 scanned) */
		if (n<=0 || i<n)  while (*s && index(ign,*s)!=NULL)  s++;
		for (j=0; j<pfnumexp; j++)  if (abs(xx.xe[j])>PF_MAXEXPONENT/2)
			{USERERROR("Too large exponent on input %d.",xx.xe[j]);}
		if (!*s || (n>0 && i>=n))  break;
			/**
			 * This part first looks at generator names, and their special powers.
			 * Except the base generators, the scanned string must not contain
			 * a character from noat for noatominp>0.
			 * The longest matching prefix-substring wins, generators first.
			 * (Globally lm keeps the longest match found so far, jm keeps the
			 *  corresponding power, and jjm the generator index.)
			 * Later, we must also check that the substring found here is not
			 * just a prefix of a longer integer(!).
			**/
		lm = jm= 0;
		for (j=-1; j<2*PF_NUMSPECPOW; j++)  for (jj=0; jj<2*pfnumexp; jj++) {
			if (j<0)  t = pfcnow->pfdef.exstr[jj];
			else if (j<PF_NUMSPECPOW)  t = pfcnow->pfdef.exstrp[j][jj];
			else  t = pfcnow->pfdef.exstrm[j-PF_NUMSPECPOW][jj];
			if (!t)  continue;
			tt = "";	/* (not accepting noat characters in the string) */
			if (pfcnow->pfdef.noatominp>0 && (j>=0 || jj>=pfnumexp))
				for (tt=noat; *tt && index(t,*tt)==NULL; tt++) ;
			if (*tt)  continue;
			l = strlen(t);		/* look for the longest match: */
			if (l>lm && strncmp(s,t,l)==0) {
				if (sscanf(s,"%d%n",&b,&d)>0)  if (d>l)  continue;
				lm = l;  jm = j;  jjm = jj;
			}
		}	/* lm>0: a generator (jjm) substring was found - add the exponent */
		if (lm>0) {
			if (jjm>=pfnumexp) {	/* (an atom string for negative value) */
				jjm -= pfnumexp;  gg = -gg;
			}
			xx.xe[jjm] += (jm<0?1: (jm<PF_NUMSPECPOW?jm+1: PF_NUMSPECPOW-jm-1));
			s += lm;  i++;
			continue;	/* an atom was scanned, continue the cycle again */
		}
			/**
			 * This part looks for integer numbers or special operators on
			 * the input ("minus", "zero", "pow").
			 * Each occurence of "minus" inverts the current number.
			 * Each occurence of "zero" zeroes the current number.
			 * An integer is scanned and tried as a power of the first generator.
			 * A "pow" symbol means that the last scanned generator is
			 * powered to the following integer.
			**/
		if (index(minus,*s)!=NULL) {	/* finds a "minus" sign */
			gg = -gg;  s++;
		} else				/* looks for a numeric value (integer) */
		  if (sscanf(s,"%d%n",&j,&l)>0) {
			jjm = 0;		/* (automatically assigned to the first generator) */
			if (j>=-1 && j<=1) {
				gg *= j;	/* (values 0,+-1 are special) */
			} else if (pfnumexp>jjm) {
						/* tries to find the power to get the scanned number */
				m = pfcnow->pfdef.valmod[jjm];
				bb = pfcnow->pfdef.exbase[jjm];
				for (jj=1,b=bb; jj<1000; b*=bb,jj++) {
					if (m>0)  b = b%m;
					if (b==j)  break;
				}
				if (b==j) {
					xx.xe[jjm] += jj;
				} else {
					ret = -1;  if (ch>0)  {USERERROR("Wrong numeric value %d for generator power %d ^x %% %d.",j,bb,m);}
				}
			} else {		/* this is when we have no generator and the number is not +-1 */
				ret = -1;  if (ch>0)  {USERERROR("Cannot read numeric value with no generator.");}
			}
			if (ret>=0)  i++;
			s += l;
		} else if (index(zero,*s)!=NULL) {	/* finds a "zero" character */
			gg = 0;  i++;
			s++;
		} else if (index(one,*s)!=NULL) {	/* finds a "one" character */
			i++;  s++;
		} else if (index(pow,*s)!=NULL) {	/* finds a "pow" character */
			s++;
			jj = sscanf(s,"%d%n",&j,&l);	/* reads the numeric exponent given after ^ */
			if (jj>0) {
				xx.xe[jjm] *= j;	/* (jjm keeps the index of the last scanned generator) */
				s += l;
				if (abs(j)>PF_MAXEXPONENT/2 || abs(j*xx.xe[jjm])>PF_MAXEXPONENT/2)
					{USERERROR("Too large exponent on input %d.",j);}
			} else {
				ret = -1;		/* invalid numeric exponent */
				if (ch>0)  {USERERROR("Cannot read exponent numeric value ^\"%s\".",s);}
			}
		} else {		/* here all invalid inputs end - error message */
			ret = -1;
			if (ch>0)  {USERERROR("Invalid character on input \"%s\".",s);}
		}
#if CURDLEV<=DEBUGLEV
		if (s-inp>(int)(strlen(inp)+1))  {PROGERROREXIT("Got past input buffer, why??? %s",inp);}
#endif
	}
	if (ret>=0 && n>0 && i<n) {	/* if not all requested atoms were scanned */
		ret = -2;
		if (ch>0)  {USERERROR("Too short input at atom index %d (\"%s\").",i,inp);}
	}
	if (ret>=0)  pfield_tostandform(&xx,&gg);
	if (ret>=0)  DEBUG(CURDLEV+1,"  Scanned  %s\n",pfield_pvalue(28,xx,gg));
	
	if (ret>=0)  ret = s-inp;	/* the length of scanned part of inp */
	if (x)  *x = xx;
	if (g)  *g = gg;	/* returning the result through given pointers */
	return ret;
}



/**
 * This function prints the given pfield value x,g into given buffer tx of max length maxl.
 * If tx==NULL, then the function's internal buffer is used for print.
 * (Be careful - this is a static buffer, and may be overflown easily. Then new prints
 * overwrite the old ones.)
 * The parameter form determines level of formality of print:
 *  2 means very formal (not practical), 1 is formal (all generators are printed as powers),
 *  0 is usual informal print (spec "power" strings may be used, informal strings),
 *  and -1 prints loosely (not guaranteed to be readable again!).
 * 
 * The format of pfield elements is directed by values stored in the structure pf_definition
 * about the current pfield.
 * The general format of printing (/scanning) pfield numbers is implicitely coded-in
 * this function, the pfield_scanvalue_ext() function, and also in the lexical scanner
 * for matrix files emflex.l.
 * This function and pfield_scanvalue_all() should be inverse to each other (for form>=0).
 * Moreover, the lexical scanner provides input for general symbolic expressions
 * (including bracketing) that are defined in the pfield.
 * 
 * The return value is the string to which the number was printed (equal to tx if given).
 * 
**/

char*	pfield_printvalue_ext(char *tx, int maxl, exp_t x, sign_t g, int form) {
	static char	buf[8190];
	static int	bufi = 0;
	int		i,j,z, m,b,u,v,nz,nb, ll,ff, xi,go=g;
	char		*s, *tt;
	
	if (form>=0)  DEBUG(CURDLEV+2,"Going to print given number (form=%d, tx=%p, maxl=%d)...\n",form,tx,maxl);
	if (!tx) if (maxl<4 || 4*maxl>=(int)(sizeof(buf)/sizeof(buf[0])))
		{ PROGERROR("Not enough space for printing  %d",maxl);  return ""; }
	
	if (tx==NULL) {		/* prepares internal buffer for printing if no one provided */
		if (bufi>=(int)(sizeof(buf)/sizeof(buf[0]))-maxl)  bufi = 0;
		tx = buf+bufi;  bufi += maxl;
	}
	tt = tx;  ll = maxl;
	pfield_tostandform(&x,&g);	/* must print in the standard form! */
	
	if (g==0) {		/* printing "zero" value */
		strncpy(tx,form>0?PFPRINTZERO1:PFPRINTZERO,maxl-1);
	}
	if (g<0) {		/* printing "minus" sign separately (but may be overwritten later) */
		strncpy(tx,form>0?PFPRINTMINUS1:PFPRINTMINUS,maxl-1);
		tx[maxl-1] = 0;  maxl -= strlen(tx);  tx += strlen(tx);
	}
	if (g!=0 && pfield_iszeroexp(x)) {	/* printing "one" value, but this can be overwritten */
		strncpy(tx,form>0?PFPRINTONE1:PFPRINTONE,maxl-1);
	}
	for (i=nz=0; i<pfnumexp; i++)	/* the number of nonzero terms in the number x,g */
		nz += (x.xe[i]!=0);
	
			/**
			 * After the above 0 1 or - characters were possibly printed,
			 * it is time to print the generators with their powers.
			 * We prefer to print the positive powers first, then the negative ones.
			 * Possibly, we overwrite the '-' sign if an informal string exists
			 * for the negative value of the generator (only the first term j==0).
			 * 
			 * The multiplication separator is printed between generators (after j>0).
			 * Then we look for a possible informal string (if form<=0)
			 * for the current power of a generator, or for its negative value.
			 * If we do not succeed, then we use the generator string with a ^power.
			**/
	if (g!=0)  for (z=j=0; z<2; z++)
	  for (i=0; i<pfnumexp; i++) {
		xi = x.xe[i];
		ff = (form>=2&&z==0) || (form==1&&z==0&&xi!=0) || (form<=0&&z==0&&xi>0) || (form<=0&&z==1&&xi<0);
		if (!ff)  continue;	/* (determines when this exponent is to be printed) */
		s = NULL;
		if (j>0) {	/* printing the multiplication separator after a previous atom (j>0) */
			strncpy(tx,form>0?PFPRINTTIMES1:PFPRINTTIMES,maxl-1);
			tx[maxl-1] = 0;  maxl -= strlen(tx);  tx += strlen(tx);
		}
				/* if allowed to use short power-strings, then tries to find them */
		if (!s && form<=pfcnow->pfdef.optstrf && form<=1
				&& xi>=-PF_NUMSPECPOW && xi<=PF_NUMSPECPOW) {
			if (g<0 && !j && xi>0)  s = pfcnow->pfdef.exstrp[xi-1][i+pfnumexp];
			if (g<0 && !j && xi<0)  s = pfcnow->pfdef.exstrm[-xi-1][i+pfnumexp];
			if (s) { tx = tt;  g = -g; }	/* (optional string for the first -gen^i found) */
			if (!s && xi>0)  s = pfcnow->pfdef.exstrp[xi-1][i];
			if (!s && xi<0)  s = pfcnow->pfdef.exstrm[-xi-1][i];
			if (s)  xi = 1;	/* (the exponent is already hiden in the opt string above) */
		}
				/* even when the prev is not allowed, try the negative generator */
		if (!s && g<0 && j<=0 && nz==1 && xi==1 && form<=1
				&& (form<=pfcnow->pfdef.optstrf || form<=0)) {
			s = pfcnow->pfdef.exstr[i+pfnumexp];
			if (s) { tx = tt;  g = -g; }	/* (optional string for the only -gen found) */
		}
				/* the default generator strings are taken (may not exist!) */
		if (!s && form<0)  s = pfcnow->pfdef.exstrp[0][i];
		if (!s)  s = pfcnow->pfdef.exstr[i];
			/**
			 * Here the string s obtained above is printed with the exponent xi.
			 * Notice that xi was possibly adjusted above, and that the negative
			 * sign could be erased a g inverted for the first term.
			 * We first decide whether brackets are really necessary (for !nb).
			 * Then we print the term.
			**/
		nb = ((s?*s:0)? s[1]==0 || s[0]=='[' :0);
		nb = nb || (nz==1 && g>0 && xi==1);
		if (s!=NULL) {
		  if (form>=2) {
			snprintf(tx,maxl,"(%s)%s(%d)",s,PFPRINTPOWER1,xi);
		  } else if (form>pfcnow->pfdef.nobracef && !nb) {
			if (xi==1)  snprintf(tx,maxl,"(%s)",s);
			else  snprintf(tx,maxl,"(%s)%s%d",s,PFPRINTPOWER1,xi);
		  } else {
			if (xi==1)  snprintf(tx,maxl,"%s",s);
			else  snprintf(tx,maxl,"%s%s%d",s,form>0?PFPRINTPOWER1:PFPRINTPOWER,xi);
		  }
			/**
			 * Finally, if no generator string is given (only for i=0 !),
			 * then the numeric value of the generator is printed instead.
			 * This depends on the numeric base and modulo given in pfdefs.
			**/
		} else {
			m = pfcnow->pfdef.valmod[i];  b = pfcnow->pfdef.exbase[i];
			for (u=1,v=b; u<xi; u++)
				if (m>0)  v = (v*b)%m;  else  v *= b;
			if (i==0)  snprintf(tx,maxl,"%d",v);
			if (xi<0)  {USERERROR("Cannot print negative power %d in numeric format.",xi);}
			if (i>0)  {USERERROR("Only the first generator may have numeric value.");}
			if (v<0 || v>1000000)  {USERERROR("Possibly too large numeric value for the power^%d:  %d.",xi,v);}
		}
		tx[maxl-1] = 0;  maxl -= strlen(tx);  tx += strlen(tx);
		j = 1;
	}
	tt[ll-1] = 0;		/* (tt,ll saves the original values of tx,maxl) */
	if (form>=0 && ll>25 && (int)strlen(tt)>ll-4)
		{USERERROR("Possibly truncated output when printing: \"%s\"",tt);}
#ifndef	FASTPROG
	if (form>=0 && form<2 && IFRANDDEBUGLESS(222)) if ((int)strlen(tt)<=ll-3) {
		exp_t	xx;	sign_t	gg;
		z = pfield_scanvalue_all(1,tt,&xx,&gg);
		if (z<0 || !pfield_isequal(x,go,xx,gg))  {PROGERROR("Incorrectly printed number -- printed\"%s\": \"%s\" x \"%s\", z=%d.",tt,pfield_pvalue(15,x,go),pfield_pvalue(15,xx,gg),z);}
	}
#endif
	return tt;
}






































