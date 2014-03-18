
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
 * Some theory for start... about "magic values" for (elements of) matrices and matroids.
 * (Read general description of "structural" functions on the top of ../include/struct.h.)
 * 
 * When working with matrices or with the underlying matroids, it is often useful to
 * consider some "magic numbers" associated with these objects.
 * The magic numbers are, by default, invariant on the operations that we consider
 * (that means scaling and permuting lines for matrices, or taking various representations
 * for an abstract matroid in the second case).
 * We may then apply various tricks using symmetric functions to these numbers.
 * 
 * Here we provide several useful and reasonably fast functions for computing such
 * magic values for matrix lines, or for matroid elements and whole matroids.
 * These functions capture the structure of small zero subdeterminants, or of
 * small-rank flats in a matroid.
 * (A similar "matroid hash" function is implemented inside ematrix_printmore()
 * in ../ematrix/ematexp.c, but that hash function is exponentially slow to compute,
 * so it should not be used inside other computations here.)
 * 
**/






#include "macek.h"
#include "str.h"












/******************	Miscelanelous	*******************/
/**********************************************************/


#undef CURDLEV
#define CURDLEV		7		/* (need <=DEBUGLEV+2 to perform paranoic tests) */


/**
 * This function decides whether there is an injection from the array fr[] to the array to[]
 * such that elements of value other than nomag are mapped to elements of the same value.
 * The return value is 0 for no injection, or d+1 where d is the number of distinct values.
 * (We typically use nomag==0, ie. only nonzero values are considered, in both fr,to!)
**/

int	strmag_isinjection_ext(long *fr, int frn, long *to, int ton, long nomag) {
	int	i,j, p,m,r;
	long	*ar,arstack[30];
	
	if (!fr || !to || frn<0 || ton<0)  {PROGERROR("both magic arrays must be given here!");}
	if (frn>ton)  return 0;
	if (frn+ton<28)  ar = arstack;
	else  ar = MMALLOC((frn+ton+2)*2*sizeof(ar[0]));
	
	for (i=0; i<frn; i++)  ar[i] = fr[i];
	for (i=0; i<ton; i++)  ar[i+frn] = to[i];
	r = 1;
	for (i=0; i<frn && r>0; i++)  if (ar[i]!=nomag) {
		m = 0;  r++;  p = ar[i];
		for (j=i; j<frn+ton; j++)
			if (ar[j]==p) {
				(j<frn? m--: m++);
				ar[j] = nomag;
			}
		if (m<0)  r = 0;
	}
	if (ar && ar!=arstack)  FREE(ar);
	return r;
}


/**
 * Here we compute a simple symmetric function of the values in the array ar[] of length arn.
 * The result is always nonzero.
**/

long	strmag_symmfunction_ext(int lev, long *ar, int arn) {
	int	i;
	long	f1,f2,f3,f4, l;
	
	f1 = f2 = f3 = f4 = 0;
	for (i=0; i<arn; i++) {
		f1 += ar[i];  f2 += (l=ar[i]*ar[i]);
		if (lev>2)  f3 += l*ar[i];
		if (lev>3)  f3 += l*l;
	}
	l = f1%1001l + f2*1001l;
	if (lev>2)  l += f3*100001l;
	if (lev>3)  l += f4*10000001l;
	if (!l)  l = 1;
	return l;
}


/**
 * 
 * combining two- to one- hash values (straut.old l. 440)........
 * 
**/














/******************	Matrix values	*******************/
/**********************************************************/


#undef CURDLEV
#define CURDLEV		7		/* (need <=DEBUGLEV+2 to perform paranoic tests) */


/**
 * This function computes the zero-pattern of the lines (rows tr==0, columns tr==1) of
 * the matrix e; storing the zeros of each line into bits of an element of zp[],
 * and the number of zeros +1 in the line to zn[].
 * The total number of zeros +1 is returned from the function.
 * Only lines indexed s1...s2-1 are considered.
 * The arrays zn,zp must be allocated by the caller to length >=s2-s1.
 * If zn or zp is NULL, then these information are not stored.
**/

long	strmag_zerolines_ext(ematrix *e, int tr, int s1, int s2, long *zp, long *zn) {
	int	i,j;
	long	p,k,l;
	
	if (!tr)  ematrix_transpose(e);
	for (l=0,j=s1; j<s2; j++) {
		k = 0;  if (zp)  zp[j-s1] = 0l;
		for (i=0; i<ROWSM(e); i++) {
			if (SIGNM(e,i,j)==0)  k++;
			else if (zp)  zp[j-s1] ^= (1l<<(i%NUM_LONG_BITS));
		}
		l += k;  if (zn)  zn[j-s1] = k+1;
	}
	if (!tr)  ematrix_transpose(e);
	return l+1;
}


/**
 * Here we find the clusters of zero 2x2 subdeterminants of the two given rows r1,r2 of ee.
 * (The zero subdeterminants form cliques if viewed as edges of a graph.)
 * The same is computed for columns if tr==1.
 * The result is stored in ar[]; ar[i] contains the index of the first column j<i for
 * which the subdeterminant of columns i,j is zero, or ar[i]=i if there is no zero.
 * The columns of both-0 entries are treated specially - they form their own cluster.
 * Subdeterminants may be undefined, then they are treated as nonzeros.
 * 
 * The function returns the number of zero 2x2 subdeterminants found +1.
 * If ar==NULL is given, then only this number is returned.
 * (The implementation may look complicated here, but it is actually faster than the trivial method.)
**/

long	strmag_twolinedet_ext(ematrix *e, int tr, int r1, int r2, int *ar) {
	int	i,j,i0,k, p,q, *nr, nrstack[30];
	sign_t	gg;
	
	if (tr)  ematrix_transpose(e);
	if (COLSM(e)<12)  nr = nrstack;
	else  nr = MMALLOC((2*COLSM(e)+4)*sizeof(nr[0]));
	if (!ar)  ar = nr+COLSM(e);
	
	i0 = -1;	/* i0 keeps the first both-0 column, ar[0] starts with 0 automatically in the cycle */
	for (i=0; i<COLSM(e); i++) {
		if (SIGNM(e,r1,i)==0 && SIGNM(e,r2,i)==0) {
			ar[i] = i0 = (i0==-1? i:i0);
			continue;		/* (both-0 columns need special treatment!) */
		}
		q = (SIGNM(e,r1,i)!=0 && SIGNM(e,r2,i)!=0);
		for (j=0; j<i; j++)  if (ar[j]==j && ar[j]!=i0) {
			if (q!=(SIGNM(e,r1,j)!=0 && SIGNM(e,r2,j)!=0))
				continue;	/* (det cannot be 0 with exactly one 0 entry) */
			p = ematrix_determinant_2x2ch(0,e,r1,r2,j,i,NULL,&gg);
			if (p>=0 && gg==0)  break;
		}	/* the subdeterminant may not be defined here, but then treated as nonzero */
		ar[i] = j;
	}
	for (i=0; i<COLSM(e); i++)  nr[i] = 0;
	for (i=k=0; i<COLSM(e); i++)  k += nr[ar[i]]++;
	if (i0>=0)  k += nr[i0]*(COLSM(e)-nr[i0]);
	
#ifndef FASTPROG
	if (IFRANDDEBUGLESS(222)) {
	  for (i=q=0; i<COLSM(e); i++) for (j=0; j<i; j++) {
		p = ematrix_determinant_2x2ch(0,e,r1,r2,j,i,NULL,&gg);
		q += (p>=0 && gg==0);
	  }if (q!=k)  {PROGERROR("something wrong with twolinedet -- number of zeros %d!=%d",k,q);}
	}
	if (CURDLEV+3<=DEBUGLEV) {
		SDEBUG(CURDLEV,"\t\t\ttwolinedet %d (+1) at %d,%d: ar[] = ",k+1,r1,r2);
		for (i=0; i<COLSM(e); i++)  SDEBUG(CURDLEV," %d,",ar[i]); SDEBUG(CURDLEV,"\n"); }
#endif
	if (nr && nr!=nrstack)  FREE(nr);
	if (tr)  ematrix_transpose(e);
	return k+1;
}


/**
 * This is a magic number computed from all entries of the given matrix.
 * It is invariant only on entry-wise equality.
**/

long	strmag_matrixentries_ext(ematrix *e) {
	long	ar[40], l=38;
	int	i,j;
	
	for (i=0; i<l; i++)  ar[i] = 0l;
	for (i=0; i<ROWSM(e); i++)  for (j=0; j<COLSM(e); j++)
		ar[(4*i+j)%l] += (long)pfield_hashvalue(EXPM(e,i,j),SIGNM(e,i,j));
	return strmag_symmfunction(ar,l);
}



/**
 * 
 * 
**/















/******************	Abstract matroid values	*******************/
/******************************************************************/


#undef CURDLEV
#define CURDLEV		6		/* (need <=DEBUGLEV+2 to perform paranoic tests) */


/**
 * This function computes magic numbers for matroid elements based on small flats in
 * the matroid.
 * These numbers are matroid-invariants!
 * The numbers are returned in the array flt[] (indexed first by rows, then by columns).
 * 
 * This function tries (co)closures of all k-tuples of elements (k = minr...maxr-1),
 * and records numbers based on the sizes of nontrivial ones.
**/

long	strmag_flatlines_ext(ematrix *e, int minr, int maxr, int tf, long flt[]) {
	ematrix		*ec, **el,**x;
	long		r,a;
	int		i,j, k,l, p,t,ro,co;
	
	if (maxr<=minr || maxr>minr+6)  maxr = minr+4;
	ro = ROWSM(e);  co = COLSM(e);
	if (flt)  for (i=0; i<ro+co; i++)  flt[i] = 0l;
	DEBUG(CURDLEV+2,"Computing flatlines for the matroid %p [%.12s] (%dx%d)...\n",e,EMNAME(e)?EMNAME(e):"",ro,co);
	r = 0;
	for (k=minr, p=0; k<maxr && k<=ro && k<=co && p<tf; k++) {
		el = ematrix_submatrices_sub(e,k);
			/**
			 * Here we cycle all k-element submatrices of e, look at nontrivial
			 * closures and coclosures, and store hash-values based on their
			 * sizes to all closure elements.
			**/
		for (x=el; x?*x:0; x++)  for (t=0; t<2; t++) {
			ec = ematrix_closure_tr(e,*x,t);
			i = ROWSM(ec)+COLSM(ec);
			a = (i*i+t*i+1l)<<(4*(k-minr));
			if (i>k)  for (p++, j=0; j<i; j++) {
				r += a;
				l = (j<ROWSM(ec)? GETREFMROW(ec,j): GETREFMCOL(ec,j-ROWSM(ec))+ro);
				if (flt)  flt[l] += a;
			}
			dispose_ematrix(ec);
		}
		dispose_alist_mats(el);
	}
#ifndef FASTPROG
	if (CURDLEV+2<=DEBUGLEV && flt) {
		SDEBUG(CURDLEV+2,"\t\t(fll up to k=%d, p=%d) flt[] =\t",k-1,p);
		for (i=0; i<ROWSM(e)+COLSM(e); i++)  SDEBUG(CURDLEV+2,"%ld, ",flt[i]);
		SDEBUG(CURDLEV+2,"\n"); }
#endif
	return r;
}















/******************	Matroid isomorphism	*******************/
/******************************************************************/


#undef CURDLEV
#define CURDLEV		6		/*6 (need <=DEBUGLEV+2 to perform paranoic tests) */


/**
 * This function tests whether the two matroids e1,e2 are isomorphic while fixing
 * the displayed basis of each one of them.
 * If xpf2>=0 is given, then the matrix e2 is considered over the pfield xpf2.
 * If he1,he2[] are given, then the searched isomorphisms must respect the hash codes
 * in heI over eI, and the same applies with lab1,lab2[] codes.
 * (The difference here is that heI[] are supposed to keep supplementary
 * matroid-invariant hash codes, while labI[] are additional "labels" for elements.)
 * 
 * For ch>=1, the isomorphism is printed out when found, all isomorphisms are
 * counted for ch>=2, and all ismorphisms are even printed out for ch>=3.
 * The return value is -2 for different-size matrices, -1 when zero entries do
 * not match, 0 when no isomorphism is found, and >0 for isomorphic.
 * For ch>=2, ret>0 gives the number of all fixed-basis isomorphisms.
**/

static int	p_sdn=0, p_nzn=0, p_ison=0;	/* (for profiling purposes...) */

int	strmag_isomorphic_fixb(int ch, ematrix *e1, ematrix *e2, long he1[], long he2[],
							long lab1[], long lab2[], int xpf2) {
	ematrix	*re2=NULL;
	char	buf[100];
	int	i,j,k,x, ro,co, ret=0, *used=NULL,*stu=NULL, ustack[80];
	long	*nzl1=NULL,*nzl2=NULL,*nlab=NULL, nstack[120];
	
	ro = ROWSM(e1);  co = COLSM(e1);
 	DEBUG(CURDLEV+2,"Testing fixed-bas isom. of %s%p[%.8s] %dx%d and %p[%.8s] (xpf=%d)...\n",
 			lab1?"labeled ":"",e1,EMNAME(e1)?EMNAME(e1):"",ro,co,e2,EMNAME(e2)?EMNAME(e2):"",xpf2);
	EMATDEBUGS(CURDLEV+3,e1,"\t\t\t"); //EMATDEBUG(CURDLEV+3,e2,"\t\t  ~\t");
	if (ISREFMAT(e2))  {PROGERROR("do not call for a refering matrix e2!");}
	if (ro!=ROWSM(e2) || co!=COLSM(e2))  ret = -2;
	
			/**
			 * We first look at zero entries which must match in an isomorphism
			 * with the current basis (regardless of the pfield of e2).
			 * (For that we allocate the arrays nzl1,nzl2[].)
			 * Then other function data are prepared for computation then...
			**/
	if (ret==0) {
		if (ro+co>36)  nzl1 = MMALLOC(3*(ro+co+2)*sizeof(nzl1[0]));
		else  nzl1 = nstack;
		nzl2 = nzl1+(ro+co+1);
		nlab = nzl2+(ro+co+1);	/* (nlab is used for debug-testing only) */
		i = strmag_numzero_rows(e1,nzl1);  j = strmag_numzero_rows(e2,nzl2);
		strmag_numzero_cols(e1,nzl1+ro);  strmag_numzero_cols(e2,nzl2+ro);
#ifndef FASTPROG
		SDEBUG(CURDLEV+3,"\t\t zeros  nzl1 = [ ");
		for (k=0; k<ro+co; k++)  SDEBUG(CURDLEV+3,"%ld,",nzl1[k]);
		SDEBUG(CURDLEV+3," ];  i=%d\n\t\t zeros  nzl2 = [ ",i);
		for (k=0; k<ro+co; k++)  SDEBUG(CURDLEV+3,"%ld,",nzl2[k]);
		SDEBUG(CURDLEV+3," ];  j=%d\n",j);
#endif
		if (i!=j || !strmag_isinjection(nzl1,ro,nzl2,ro) ||
				!strmag_isinjection(nzl1+ro,co,nzl2+ro,co)) {
			if (ch>=0)  p_nzn++;  ret = -1;
			DEBUG(CURDLEV+2,"Rejected fixed-bas isom. of %p[%.8s] %dx%d and %p[%.8s] by zerolines.\n",e1,EMNAME(e1)?EMNAME(e1):"",ro,co,e2,EMNAME(e2)?EMNAME(e2):"");
		}
	}
	if (ret==0) {
		if (ro+co>36)  used = MMALLOC(2*(ro+co+2)*sizeof(used[0]));
		else  used = ustack;
		stu = used+(ro+co+1);
		for (i=0; i<ro+co; i++) { used[i] = 0;  stu[i] = -1; }
		re2 = ematrix_refer_empty(e2);
	}
			/**
			 * Here the main computation starts:
			 * We try all row and column permutations of the matrix e2 which
			 * respect the hash codes in he2[] against he1[], the "labels"
			 * of elements in lab1,lab2[], and which also preserve
			 * the (local) numbers of zeros in lines nzl2[] against nzl1[] (e1).
			 * Each selected permutation is then passed to ematrix_havesamedets()
			 * to find out whether their two matroids are identical.
			 * 
			 * We cycle the search until an isomorphism is found, or we find
			 * all isomorphisms if ch>=2 is given (all printed out when ch>=3).
			**/
	k = 0;
	while (k>=0 && re2 && (ret==0 || ch>=2)) {
		DEBUG(CURDLEV+4,"%*sk=%d, st=%d;  %dx%d\n",2*k,"",k,stu[k],ROWSM(re2),COLSM(re2));
		if (k>ro+co)  {PROGERROR("cannot get here!! k=%d\n",k); break;}
			/**
			 * Here the next choice for line k (in re2) is searched, using
			 * the hash codes in he?[] and in nzl?[].
			 * The array used[] marks matrix lines (of e2) which are already mapped to.
			 * If nothing is found as the next choice, then we go back
			 * and remove the last added line from re2 and from used[].
			**/
		j = -1;
		if (k<ro+co) for (i=stu[k]+1; i<ro+co && j<0; i++)
			if (((he1&&he2)? he2[i]==he1[k]:1) && ((lab1&&lab2)? lab2[i]==lab1[k]:1) &&
					(nzl1? nzl2[i]==nzl1[k]:1) && !used[i])
				if ((k<ro)==(i<ro))  j = i;
		if (j<0) {
			k--;
			if (k>=0 && k<ro) {
				used[GETREFMROW(re2,k)] = 0;
				ematrix_remove_row(re2,k);
			} else  if (k>=ro && k<ro+co) {
				used[GETREFMCOL(re2,k-ro)+ro] = 0;
				ematrix_remove_col(re2,k-ro);
			}
			continue;
		}
			/**
			 * Here we pass the successful next choice for the line k of re2.
			 * We add this line, and we proceed further in backtracking
			 * until we get the whole matrix re2 -> e2.
			**/
		stu[k] = j;  used[j] = 1;
		if (j>=0 && j<ro)  ematrix_refadd_row(re2,j);
		if (j>=ro && j<ro+co)  ematrix_refadd_col(re2,j-ro);
		DEBUG(CURDLEV+3,"%*s -> k=%d, nx=%d;  %dx%d\n",2*k,"",k,stu[k],ROWSM(re2),COLSM(re2));
		stu[++k] = -1;
		if (k<ro+co)  continue;
			/**
			 * Here re2 is a whole permutation of e2, and we test it against e1.
			 * If we find out that their subdeterminants match (ematrix_havesamedets),
			 * then we record this new isomorphism, and possibly print it out.
			**/
		if (ROWSM(re2)!=ro || COLSM(re2)!=co)  {PROGERROR("something wrong with matrix permutation!! %d!=%d || %d!=%d\n",ROWSM(re2),ro,COLSM(re2),co);}
		i = ematrix_havesamedets(e1,re2,xpf2);
		if (ch>=0)  p_sdn++;	/* (for profiling statistics...) */
		if (i)  ret = (ret<=0? 1:ret+1);
		
			/**
			 * Next is for debug and output printing, no more computation...
			 * The first isomorphic matrix is printed out if ch>=1, and all
			 * isomorphic matrices are printed out if ch>=3 (ch=2 only counts all).
			**/
#ifndef FASTPROG
		if (i && CURDLEV+1+(ret+p_ison>1)<=DEBUGLEV) {
			DEBUG(CURDLEV+2,"Found a #%d permutation of [%.11s] with same subdeterminants.\n",ret,EMNAME(e2)?EMNAME(e2):"");
			x = pfield_curindex();  if (xpf2>=0) pfield_switchto_fast(xpf2);
			EMATDEBUGS(CURDLEV+3,re2,"\t\t  =\t");  if (xpf2>=0) pfield_switchto_fast(x);
		}
#endif			/* (p_ison keeps the total number of isomorphisms found, it is used to prevent too much printing out) */
		if (i && ch>=1+2*(ret+p_ison>1)) {
			OUTPUT("  Found an isomorphism #%d of %smatroids [%.12s] to [%.12s] (%dx%d)\n",
					ret,lab1?"labeled ":"",EMNAME(e1)?EMNAME(e1):"",EMNAME(e2)?EMNAME(e2):"",ro,co);
			if (ret==1)  EMATOUTPUTS(e1,printoutpref);
			snprintf(buf,90,"%s=",printoutpref);
			x = pfield_curindex();  if (xpf2>=0) pfield_switchto_fast(xpf2);
			EMATOUTPUTS(re2,buf);
			if (xpf2>=0) pfield_switchto_fast(x);
		}
	}
	if (ch>=2 && ret>0)  OUTPUT("  Found %d isomorphisms (%d so far) for a fixed basis of [%.12s] to [%.12s]%s.\n",
					ret,ret+p_ison,EMNAME(e1)?EMNAME(e1):"",EMNAME(e2)?EMNAME(e2):"",lab1?" labeled":"");
	if (ch>=3 && ret<=0)  OUTPUT("  Found -NO- isomorphism for a fixed basis of [%.12s] to [%.12s]%s.\n",
					EMNAME(e1)?EMNAME(e1):"",EMNAME(e2)?EMNAME(e2):"",lab1?" labeled":"");
#ifndef FASTPROG
 	if (ret==0)  DEBUG(CURDLEV+2,"No fixed-bas isom. of%s %p[%.8s] %dx%d to %p[%.8s] found.\n",
 				lab1?" labeled":"",e1,EMNAME(e1)?EMNAME(e1):"",ro,co,e2,EMNAME(e2)?EMNAME(e2):"");
	if (ch>=0 && IFRANDDEBUGLESS(333) && ret>=-1 && ro<=5 && co<=5) {
		ematrix	*ee = ematrix_copy(e2);
		DEBUG(CURDLEV+2,"== Recursive fixed-basis isomorphism test...\n");
		if (EMNAME(e2))  EMSETNAME(ee,EMNAME(e2));
		x = pfield_curindex();  if (xpf2>=0) pfield_switchto_fast(xpf2);
		ematrix_swaprows(ee,0,ro-1);  ematrix_swapcols(ee,0,co-1);
		if (lab2) {		/* (element labels must be switched with matrix lines...) */
			for (i=0; i<ro+co; i++)  nlab[i] = lab2[i];
			nlab[0] = lab2[ro-1];  nlab[ro-1] = lab2[0];
			nlab[ro] = lab2[ro+co-1];  nlab[ro+co-1] = lab2[ro];
		}
		i = strmag_isomorphic_fixb(-1,ee,e1,NULL,NULL,(lab2?nlab:NULL),lab1,x);
		/* the hash codes he?[] are not used recursively, so a failure may mean they were wrong! */
		if (xpf2>=0) pfield_switchto_fast(x);
		if ((i>0)!=(ret>0))  {PROGERROR("Failed recursive fixed-basis%s isomorphism test!! %d X %d %s",
						lab2?" labeled":"",i,ret,he1?"\n\tMay the element \"hash codes\" be given wrong??":"");}
		dispose_ematrix(ee);
	}
#endif
	if (re2)  dispose_ematrix(re2);
	if (used && used!=ustack)  FREE(used);
	if (nzl1 && nzl1!=nstack)  FREE(nzl1);
	return ret;
}


/**
 * This function tests whether the two matroids e1,e2 are isomorphic (weaker than equivalence).
 * If xpf2>=0 is given, then the matrix e2 is considered over the pfield xpf2.
 * The return value is >0 for isomorphic, and <=0 otherwise.
 * Moreover, one may provide matroid-invariant hash code for matrix lines in hl1[],hl2[]
 * (indexed first by rows, then by columns in one array).
 * Giving these codes (as much distinguishing as possible) speeds up the computation...
 * Similar are line-codes (possibly) given in lab1,lab2[], but these codes must be
 * respected by all isomoprhisms (like giving "labels" to matroid elements).
 * 
 * For ch>=1, the isomorphism is printed out when found, all isomorphisms are
 * counted for ch>=2, and all ismorphisms are even printed out for ch>=3.
 * The return value is -2 for different-size matrices, -1 when zero entries do
 * not match, 0 when no isomorphism is found, and >0 for isomorphic.
 * For ch>=2, ret>0 gives the number of all (computed) isomorphisms.
**/

int	strmag_isomorphic_ext(int ch, ematrix *e1, ematrix *e2, int xpf2,
					long lab1[], long lab2[], long hl1[], long hl2[]) {
	ematrix		*ee1=NULL, **x, **bl=NULL;
	int		i,j, xpf1,ro,co, ret=0;
	long		*he1=NULL,*hee1=NULL,*he2=NULL, *labb1=NULL;
	
	if (!e1 || !e2)  {PROGERROR("The matrices e1,e2 must be given here!"); return -11;}
	xpf1 = pfield_curindex();
	if (xpf2==xpf1)  xpf2 = -1;
	if (!e1 || !e2)  return -4;
	ro = ROWSM(e1);  co = COLSM(e1);
	p_sdn = p_nzn = p_ison = 0;	/* (used for profiling statistics...) */
#ifndef FASTPROG
	DEBUG(CURDLEV+1,"Testing isomorphism of%s %p[%.8s] %dx%d and %p[%.8s] (xpf=%d)...\n",
			lab1?" labeled":"",e1,EMNAME(e1)?EMNAME(e1):"",ro,co,e2,EMNAME(e2)?EMNAME(e2):"",xpf2);
	if (ch>=0) EMATDEBUGS(CURDLEV+2,e1,"\t\t");
	if (lab1) for (i=0; i<ro+co; i++)  SDEBUG(CURDLEV+2,"%s%ld%s",!i?"\t\tlab  (":"",lab1[i],i<ro+co-1?", ":")\n");
	if (xpf2>=0) { pfield_switchto_fast(xpf2); DEBUG(CURDLEV+2,"  ** over another pfield %.8s :\n",pfield_curname()); }
	if (ch>=0) EMATDEBUGS(CURDLEV+2,e2,"\t  ~\t");
	if (lab2) for (i=0; i<ro+co; i++)  SDEBUG(CURDLEV+2,"%s%ld%s",!i?"\t\tlab  (":"",lab2[i],i<ro+co-1?", ":")\n");
	if (xpf2>=0) pfield_switchto_fast(xpf1);
#endif
	if (ro!=ROWSM(e2) || co!=COLSM(e2))  ret = -3;
	if (ch<=0 && !lab1 && !lab2) {
		if (ret==0 && xpf2<0) if (ematrix_isequal(e1,e2))  ret = 3;
		if (ret==0) if (ematrix_havesamedets(e1,e2,xpf2))  ret = 2;
	}
	if (ch>=0 && ret>0)  DEBUG(CURDLEV+1,"  - easy detection of matroid identity %d...\n",ret);
	
			/**
			 * We compute matroid-invariant hash codes for the matrix lines first.
			 * (Or, these codes may be given in hl1[], hl2[] for faster computation...)
			 * The codes here are simpler than those given from strmag_isomorphlist_().
			 * These codes must match for the two matrices to have an isomorphism.
			**/
	if (ret==0) {
		he1 = MMALLOC(4*(ro+co+2)*sizeof(he1[0]));
		he2 = he1+(ro+co+1);  hee1 = he2+(ro+co+1);
		labb1 = (lab1? hee1+(ro+co+1):NULL);
		if (hl1)  for (i=0; i<ro+co; i++)  he1[i] = hl1[i];
		else  strmag_flatlines(e1,2,he1);
		if (lab1&&lab2)  for (i=0; i<ro+co; i++)  he1[i] += lab1[i]*(he1[i]+3);
		if (xpf2>=0)  pfield_switchto_fast(xpf2);
		if (hl2)  for (i=0; i<ro+co; i++)  he2[i] = hl2[i];
		else  strmag_flatlines(e2,2,he2);
		if (lab1&&lab2)  for (i=0; i<ro+co; i++)  he2[i] += lab2[i]*(he2[i]+3);
		if (xpf2>=0)  pfield_switchto_fast(xpf1);
		if (!strmag_isinjection(he1,ro+co,he2,ro+co))  ret = -2;
	}
			/**
			 * If the line codes match over the whole matroids, then we generate all
			 * bases of e1 which match by rows the line codes in e2.
			 * For each one basis we call strmag_isomorphic_fixb() to find out,
			 * whether there is an isomorphism fixing that basis.
			 * To pass the line hash codes for e1 through, we have to rearrange them
			 * to match the lines in the particular basis, using line id's.
			 * The same applies to (optional) element labels in lab1,lab2[].
			 * (However, we do not pass lab1,lab2 separately to ematrix_getbases_hh()
			 *  since the labels are already included in the hash codes he1,he2.)
			 * When ch<=1, we finish the cycle after the first isomorphism.
			**/
	if (ret==0) {
		ee1 = ematrix_copy_notr(e1);  ematrix_resetid(ee1);
		bl = ematrix_getbases_hh(ee1,he1,he2);
		if (ROWSID(ee1,0)!=1 || COLSID(ee1,0)!=-1)  {PROGERROR("We expect standard line labelling of ee1!");}
		DEBUG(CURDLEV+1," - isomorphism trying list of %d (same-hash%s) bases of %p[%.8s] %dx%d...\n",
				alist_getlength(bl),lab1?" labeled":"",e1,EMNAME(e1)?EMNAME(e1):"",ROWSM(e1),COLSM(e1));
		
		for (x=bl; x?*x:0; x++) {
			if (EMNAME(e1))  EMSETNAME(*x,EMNAME(e1));
			for (i=0; i<ro+co; i++) {
				j = (i<ro? ROWSID(*x,i): COLSID(*x,i-ro));
				hee1[i] = he1[j>0? j-1: ro-j-1];
				if (lab1)  labb1[i] = lab1[j>0? j-1: ro-j-1];
			}
			i = strmag_isomorphic_fixb(ch,*x,e2,hee1,he2,labb1,lab2,xpf2);
			if (i>0)  p_ison = (ret += i);
			if (i>0)  DEBUG(CURDLEV+1,"  - found %d isomorphisms for the basis #%d of [%.8s].\n",i,(int)(x-bl)+1,EMNAME(e1)?EMNAME(e1):"");
			if (ch<=1 && ret>0)  break;
		}
		if (ret==0)  ret = -1;
	}
			/**
			 * Output, recursive debug testing, and final cleaning...
			**/
 	if (ch>0)  OUTPUT(" Isomorphism test of%s [%.15s] (%dx%d over %.7s) to [%.15s] (over %.7s) found %s.\n",
 			lab1?" labeled":"",EMNAME(e1)?EMNAME(e1):"",ro,co,pfield_curname(),EMNAME(e2)?EMNAME(e2):"",xpf2>=0?pfield_ixname(xpf2):"cur",(ret>0?"+YES+":"-NO-"));
 	if (ch>=2 && ret>0)  OUTPUT("   (Total %d isomorphisms found...)\n",ret);
 	
#ifndef FASTPROG
 	if (ch>=0)  DEBUG(CURDLEV+0,"Testing isomorphism of%s [%.8s] %dx%d and [%.8s] (xpf=%d) found %s.\n",
			lab1?" labeled":"",EMNAME(e1)?EMNAME(e1):"",ro,co,EMNAME(e2)?EMNAME(e2):"",xpf2,(ret>0?"+YES+":"-NO-"));
 	if (ch>=0 && bl)  DEBUG(CURDLEV+0,"\t(used %d bases, %d rej by zeros, %d samedet calls)\n",alist_getlength(bl),p_nzn,p_sdn);
 	if (ch>=0 && IFRANDDEBUGLESS(555) && ret>=-2) {
 		x = (bl? bl+alist_getlength(bl)/2: NULL);
 		DEBUG(CURDLEV+1,"== Recursive full isomorphism test for %p, %p...\n",x?*x:e1,e2);
 		if (!lab1 && !lab2) {
	 		i = strmag_isomorphic_ext(-1,((x?*x:0)?*x:e1),e2,xpf2,NULL,NULL,NULL,NULL);
	 	} else if (xpf2<0) {
	 		i = strmag_isomorphic_ext(-1,e2,e1,-1,lab2,lab1,NULL,NULL);
	 	} else {
	 		i = strmag_isomorphic_ext(-1,e1,e2,xpf2,lab1,lab2,NULL,NULL);
	 	}
		if ((i>0)!=(ret>0))  {PROGERROR("Failed full recursive isomorphism test!! %d X %d\n",i,ret);}
 	}
#endif
	if (bl)  dispose_alist_mats(bl);
	if (ee1)  dispose_ematrix(ee1);
	if (he1)  FREE(he1);
	return ret;
}


/**
 * This is a variant of the above isomorphism test in which the element mp1 must be
 * mapped to the element mp2 of e2 (mp? considered as indices, first rows then columns).
 * See strmag_isomorphic_ext() for further description...
**/

int	strmag_isomorphic_map(int ch, ematrix *e1, ematrix *e2, int xpf2,
					int mp1, int mp2, long hl1[], long hl2[]) {
	long	*lab1=NULL, *lab2=NULL;
	int	i,ro,co;
	
	ro = ROWSM(e1);  co = COLSM(e1);
	lab1 = MMALLOC(2*(ro+co+2)*sizeof(lab1[0]));  lab2 = lab1+(ro+co+1);
	for (i=0; i<ro+co; i++)  lab1[i] = lab2[i] = 0;
	if (mp1<0 || mp1>=ro+co || mp2<0 || mp2>=ro+co)  {PROGERROR("Isom map indices %d,%d out of bounds!",mp1,mp2);}
	else { lab1[mp1] = 7;  lab2[mp2] = 7; }
	return  strmag_isomorphic_ext(ch,e1,e2,xpf2,lab1,lab2,hl1,hl2);
}



/**
 * This function tests given list(s) for isomoprhic pairs of matroids.
 * If only el1 is given, then all pairs of this list are tested together.
 * If both el1, el2 are given, then each matroid of el1 is tested against all matroids of el2.
 * The arrays xpf1[],xpf2[] give the indices of pfields over which the given
 * matroids are represented (no restrictions), or NULL for the current pfield.
 * 
 * If ch<=0, then nothing is printed out.
 * If ch>=2, then brief info about unique (el2==NULL) or isomorphic to second l. (el2!=NULL)
 * is printed out during computation.
 * If ch>=4, then the isomorphisms (first one for each pair) are printed out.
 * If ch>=5(6), then all isomorphisms for each pair are printed out.
 * 
 * Moreover, the function returns through *out info about isomorphic classes.
 * (Do not forget to free the memory from *out later!)
 * For el2==NULL, index in (*out)[i]=j means that j<i is the first matroid isomorphic to #i.
 * For el2!=NULL, index in (*out)[i]=j means that matroid #i is isomorphic to #j in el2,
 * only the first one of el2 is stored.
 * The return value is the number of isomorphic pairs we have found...
**/

int	strmag_isomorphlist_ext(int ch, ematrix **el1, int xpf1[],
					ematrix **el2, int xpf2[], int **out) {
	ematrix	**x,**y, **el;
	int	i,j,l,r, xpf,xpx1,xpx2, *ra=NULL, ret=0, uq=0;
	long	**hl1=NULL,**hl2=NULL, *hdt=NULL;
	
	if (!el1)  {PROGERROR("The first list el1 must be given!");}
	DEBUG(CURDLEV+0,"Testing isomorphisms list%s (len %d,%d) at ch=%d...\n",el2?"s":"",alist_getlength(el1),alist_getlength(el2),ch);
	DEBUG(CURDLEV+0,"  (mat names  %.8s, %.8s, %.8s,.. X %.8s, %.8s,..)\n",el1[0]?EMNAME(el1[0]):"",(el1[0]&&el1[1])?EMNAME(el1[1]):"",
			(el1[0]&&el1[1]&&el1[2])?EMNAME(el1[2]):"",(el2?el2[0]:0)?EMNAME(el2[0]):"",(el2?el2[0]&&el2[1]:0)?EMNAME(el2[1]):"");
	if (xpf1 && xpf2)  DEBUG(CURDLEV+1,"  (mat pfields  %.8s, %.8s, %.8s,.. X %.8s, %.8s,..;  current %.8s)\n",el1[0]?pfield_ixname(xpf1[0]):"",(el1[0]&&el1[1])?pfield_ixname(xpf1[1]):"",
			(el1[0]&&el1[1]&&el1[2])?pfield_ixname(xpf1[2]):"",(el2?el2[0]:0)?pfield_ixname(xpf2[0]):"",(el2?el2[0]&&el2[1]:0)?pfield_ixname(xpf2[1]):"",pfield_curname());
	if (!el2 && ch>=1)
		OUTPUT(" Testing matroid isomorphism pairs in list {%.8s %.8s %.8s}%d...\n",
			el1[0]?EMNAME(el1[0]):"",(el1[0]&&el1[1])?EMNAME(el1[1]):"",(el1[0]&&el1[1]&&el1[2])?EMNAME(el1[2]):"",alist_getlength(el1));
	if (el2 && ch>=1) if (alist_getlength(el1)*alist_getlength(el2)>1)
		OUTPUT(" Testing matroid isomorphism in lists {%.8s %.8s %.8s}%d {%.8s %.8s}%d...\n",
			el1[0]?EMNAME(el1[0]):"",(el1[0]&&el1[1])?EMNAME(el1[1]):"",(el1[0]&&el1[1]&&el1[2])?EMNAME(el1[2]):"",
			 alist_getlength(el1),(el2?el2[0]:0)?EMNAME(el2[0]):"",(el2?el2[0]&&el2[1]:0)?EMNAME(el2[1]):"",alist_getlength(el2));
	xpf = pfield_curindex();
	
			/**
			 * We first prepare the computation data:
			 *  - hl1,hl2 store the line hash codes for matrices in el1,el2.
			 *  - hdt is the array hl1,hl2 point to.
			 *  - ra records the isomorphic pairs as found here, returned through *out.
			**/
	el = alist_applist(alist_copy(el1),(el2?alist_copy(el2):NULL));
	l = 0;
	for (x=el; x?*x:0; x++)  l += ROWSM(*x)+COLSM(*x);
	hdt = MMALLOC((l+10)*sizeof(hdt[0]));
	hl1 = MMALLOC((alist_getlength(el)+2)*sizeof(hl1[0]));
	hl2 = hl1+(el2? alist_getlength(el1):0);
	ra = MMALLOC((alist_getlength(el1)+2)*sizeof(ra[0]));
			/**
			 * Here we pre-compute the line hash codes ("flatlines") for the
			 * lines of matrices in el1,el2 (lists concatenated in el),
			 * and store them in hl1[] (and in hl2[] implicitly).
			 * We have to switch to their pfields (see xpf1,2) before computation!
			**/
	DEBUG(CURDLEV+1,"Computing thorough flatlines for mats");
	for (i=l=0, x=el; x?*x:0; i++, x++) {
		SDEBUG(CURDLEV+1,".");
		hl1[i] = hdt+l;
		l += ROWSM(*x)+COLSM(*x);
		if (xpf1 && i<alist_getlength(el1))  xpx1 = xpf1[i];
		else if (xpf1 && xpf2)  xpx1 = xpf2[i-alist_getlength(el1)];
		else  xpx1 = -1;
		if (xpx1==xpf)  xpx1 = -1;
		if (xpx1>=0)  pfield_switchto_fast(xpx1);
		strmag_flatlines_more(*x,2,hl1[i]);
		if (xpx1>=0)  pfield_switchto_fast(xpf);
	}
	SDEBUG(CURDLEV+1,"\n");
	for (i=0; i<alist_getlength(el1); i++)  ra[i] = -1;
	
			/**
			 * The main computation starts here:
			 * We consider all matrix pairs in el1, or el1 against el2,
			 * and test these pairs for isomorphism using strmag_isomorphic_().
			 * We skip tests which follow from transitivity, or tests to next
			 * members of el2 when an isomorphic one hes been found for ch<=3.
			**/
	DEBUG(CURDLEV+1,"Testing isomorphic pairs from the lists (%d,%d)...\n",alist_getlength(el1),alist_getlength(el2));
	for (x=el1, i=0; x?*x:0; x++, i++) {
		if (ch>=4) {
			if (ch>=5)  SOUTPUT("\n");
			OUTPUT("  Looking for all isomorphisms of the matroid #%d [%.18s] %dx%d in the list.\n",
						i+1,EMNAME(*x)?EMNAME(*x):"",ROWSM(*x),COLSM(*x));
			if (ch>=5) { EMATOUTPUT(*x,printoutpref);  SOUTPUT("\n"); }
		}
		if (!el2) {
			if (ra[i]<0 && ch>=3)  OUTPUT("Matroid #%d [%.20s] %dx%d is (isom) unique (%d) in the list so far.\n",
							i+1,EMNAME(*x)?EMNAME(*x):"",ROWSM(*x),COLSM(*x),uq+1);
			if (ra[i]<0) { ra[i] = i;  uq++; }
		}
	  for (y=el2?el2:x+1, j=el2?0:i+1; y?*y:0; y++, j++) {
		if (el2? (ra[i]>=0 && ch<=4): (ra[j]>=0 || (ra[i]>=0&&ra[i]<i)))  continue;
			/**
			 * The isomorphism test is called here.
			 * Current pfield is switched to that one of *x, and the second one 
			 * of *y is passed to strmag_isomorphic_().
			 * The above precomputed line hash codes are passed in hl1[i],hl2[j].
			 * (hl2 equals hl1 when el2==NULL.)
			 * Then the test results are recorded, and possibly printed out.
			**/
		if (ROWSM(*x)!=ROWSM(*y) || COLSM(*x)!=COLSM(*y))  continue;
		if (xpf1) {
			xpx1 = (xpf1[i]==xpf? -1:xpf1[i]);
			if (xpx1>=0)  pfield_switchto_fast(xpx1);
			xpx2 = (el2? (xpf2?xpf2[j]:-1): xpf1[j]);
			if (xpx1>=0)  DEBUG(CURDLEV+2,"   Now using this pfield %s, second one %s.\n",pfield_curname(),xpx2>=0?pfield_ixname(xpx2):"same");
		} else  xpx1 = xpx2 = -1;
		r = strmag_isomorphic_ext((ch>3? ch-3:0),*x,*y,xpx2, NULL,NULL,hl1[i],hl2[j]);
		if (xpx1>=0)  pfield_switchto_fast(xpf);
		
		if (r>0)  ret++;
		if (r>0 && ch>=2) if (alist_getlength(el1)*alist_getlength(el2)>1)
			OUTPUT(" Matroid #%d [%.20s] is isomorphic to mat #%d [%.20s] (ct %d/%d).\n",
					i+1,EMNAME(*x)?EMNAME(*x):"",j+1,EMNAME(*y)?EMNAME(*y):"",el2?uq+1:uq,i+1);
		if (ch>=2) if (alist_getlength(el1)*alist_getlength(el2)==1)
			OUTPUT(" Matroid [%.20s] (over %.7s) %s isomorphic to matroid [%.20s] (over %.7s).\n",
					EMNAME(*x)?EMNAME(*x):"",(xpx1>=0?pfield_ixname(xpx1):"cur"),(r>0?" +IS+ ":"is -NOT-"),EMNAME(*y)?EMNAME(*y):"",(xpx2>=0?pfield_ixname(xpx2):"cur"));
		if (r>0 && el2) {
			DEBUG(CURDLEV+0,"Matroid #%d [%.8s] of el1 is isomorphic to matroid #%d [%.8s] of el2.\n",i+1,EMNAME(*x)?EMNAME(*x):"",j+1,EMNAME(*y)?EMNAME(*y):"");
			if (ra[i]<0) { ra[i] = j;  uq++; }
		} else if (r>0) {
			DEBUG(CURDLEV+0,"Matroid #%d [%.8s] is isomorphic to prev matroid #%d [%.8s].\n",j+1,EMNAME(*y)?EMNAME(*y):"",i+1,EMNAME(*x)?EMNAME(*x):"");
			ra[j] = i;
		}
	}}
	if (el2 && ch>=1) if (alist_getlength(el1)*alist_getlength(el2)>1)
		OUTPUT("Total %d matroids in {%.8s %.8s %.8s}%d have isomorphic mats in {%.8s %.8s}%d.\n",
			uq,el1[0]?EMNAME(el1[0]):"",(el1[0]&&el1[1])?EMNAME(el1[1]):"",(el1[0]&&el1[1]&&el1[2])?EMNAME(el1[2]):"",
			 alist_getlength(el1),(el2?el2[0]:0)?EMNAME(el2[0]):"",(el2?el2[0]&&el2[1]:0)?EMNAME(el2[1]):"",alist_getlength(el2));
	if (!el2 && ch>=1)  OUTPUT("Total %d matroids from {%.8s %.8s %.8s}%d are unique to matroid isomorphism.\n",
			uq,el1[0]?EMNAME(el1[0]):"",(el1[0]&&el1[1])?EMNAME(el1[1]):"",(el1[0]&&el1[1]&&el1[2])?EMNAME(el1[2]):"",alist_getlength(el1));
	
	DEBUG(CURDLEV-1,"Found %d isomorphic pairs of matroids in the computation...\n",ret);
	if (xpf!=pfield_curindex())  {PROGERROR("The current pfield has changed!");}
	if (el)  alist_free(el);
	if (hl1)  FREE(hl1);  if (hdt)  FREE(hdt);
	if (out)  *out = ra;
	else if (ra)  FREE(ra);
	return ret;	/* (isn't it better to return uq ???) */
}






































