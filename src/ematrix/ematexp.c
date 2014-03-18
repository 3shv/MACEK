
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
 * This file contains sources for certain exponential and brute-force computations
 * with our matrices.
 * Be careful as the functions here may take really long time to compute, and use
 * huge amounts of memory (that is why they are placed separately).
 * 
 * Read about matrices in ../include/ematrix.h ...
 * 
**/





#include "macek.h"  
#include "emx.h"















/******************	Generating submatrices	*******************/
/******************************************************************/


#undef CURDLEV
#define CURDLEV		7		/*7*/


/**
 * The following function prepares a list of all submatrices of the given matrix e
 * according to several given parameters:
 * If sq>0, then all submatrices must be squares.
 * If nz>0, then the list may (but may not!) skip those singular ones...
 * If tot>=0, then the submatrices have total tot lines (rows+columns).
 * If mrk>=0, then only submatrices of rank <=mrk are generated (no others here).
 * The parameters mins,maxs are the lower and upper bounds for the submatrix rows/cols
 * (not the total size!).
 * The parameters kr,kc >=0 request to include these row/column to the submatrices.
 * If rnd>0, then only random partial list is generated:
 *  - all submatrices of size up to 3 (up to 2) are generated if rnd==1 (rnd==2);
 *  - the list gets sparser with higher values of rnd>2.
 * 
 * The submatrices are returned in a list of refering submatrices to e.
 * Do not forget to dispose them all later.
**/

ematrix**	ematrix_submatrices_ext(ematrix *e, int sq, int nz,
				int rnd, int tot, int mrk, int mins, int maxs, int kr, int kc) {
	ematrix		*rr, **outl;
	int		i,j,k,l, p,s, *rows;
	
	if (ISREFMAT(e))  {PROGERROR("Do not call for a refering matrix e!"); return NULL;}
	DEBUG(CURDLEV+1,"generating all submatrices of %p %dx%d with params\n\t\t\tsq=%d, nz=%d, rnd=%d, tot=%d, mrk=%d, mins=%d, maxs=%d, kr=%d, kc=%d.\n",
			e,ROWSM(e),COLSM(e),sq,nz,rnd,tot,mrk,mins,maxs,kr,kc);
	rows = MMALLOC((ROWSM(e)+5)*sizeof(rows[0]));
	outl = NULL;  k = -1;
	while (1) {			/* the main cycle generating row selections starts here */
		if (k>=0 && rnd>0 && rnd+k>=4) {	/* choosing randomized submatrices */
			rows[k] += RANDOM()%(ROWSM(e)/7+k/2+rnd/2);
		}
		if (k>=0) if (++rows[k]>=ROWSM(e)) {
			if (--k<0)  break;  	/* this breaks the main cycle at the end */
			continue;
		}
		DEBUG(CURDLEV+2,"%*s choice rows[%d]=%d\n",2*(k+1)," ",k,(k>=0?rows[k]:-1));
		if (kr>=0) {			/* request to keep the row kr in the selection */
			for (i=p=0; i<=k; i++)  if (rows[i]==kr)  p = 1;
		} else  p = 1;
		if (mrk==0 && k>=1) {		/* request for rank 0 of the submatrix */
			for (j=s=0; j<COLSM(e) && !s; j++) {
				for (i=l=0; i<=k; i++)  l = (l || SIGNM(e,rows[i],j)!=0);
				s = (s || l==0);
			}
		} else  s = 1;
		if (p && s)			/* passing the row selection for columns below (even k=-1) */
		  if ((k>=mins-1 || mins<0) && (k<=maxs-1 || maxs<0) 
		  			&& (k<=tot-1 || tot<0) && (k+COLSM(e)>=tot-1 || tot<0)) {
		  				/* (we check if it is possible to fill the submatrix wrt mins,maxs,tot) */
			DEBUG(CURDLEV+2,"\t\t+ row choice (%c%c%c%c..) accepted at k=%d +\n",(k>=0?rows[0]+'0':'.'),(k>=1?rows[1]+'0':'.'),(k>=2?rows[2]+'0':'.'),(k>=3?rows[3]+'0':'.'),k+1);
			rr = ematrix_refer_empty(e);
			for (i=0; i<=k; i++)  ematrix_refadd_row(rr,rows[i]);
			//********* use ematrix_matrank(rr)>=k+1 for nz ??? - no, must add all columns!!!
			outl = ematrix_submatrices_rows(e,rr,outl,sq,nz,rnd,tot,mrk,mins,maxs,kr,kc);
			dispose_ematrix(rr);
		}				/* when trying larger subset as a selection */
		if (k<0?1: ((p || rows[k]<kr) && (k<maxs-1 || maxs<0) && s
				&& (ROWSM(e)-rows[k]+k>=mins || mins<0)	&& (!sq || k<COLSM(e)-1))) {
			++k;
			rows[k] = (k>0? rows[k-1]:-1);
		}				/* (must always go up if k<0 !) */
	}
	FREE(rows);
	DEBUG(CURDLEV+1," - generated %d submatrices.\n",alist_getlength(outl));
	return outl;
}

/**
 * Separate selection of columns for the above selection of rows...
**/

ematrix**	ematrix_submatrices_rows(ematrix *e, ematrix *rr, ematrix **outl, int sq, int nz,
				int rnd, int tot, int mrk, int mins, int maxs, int kr, int kc) {
	ematrix		*rc, **oo;
	int		i,k, p,q, cy, *cols, colstack[20];
	
	DEBUG(CURDLEV+3,"generating all submatrices of %p %dx%d col with params\n\t\t\tsq=%d, nz=%d, rnd=%d, tot=%d, mrk=%d, mins=%d, maxs=%d, kr=%d, kc=%d.\n",e,ROWSM(e),COLSM(e),sq,nz,rnd,tot,mrk,mins,maxs,kr,kc);
	if (COLSM(e)<16)  cols = colstack;
	else  cols = MMALLOC((COLSM(e)+5)*sizeof(cols[0]));
	junk = nz;  junk = kr;
	cy = 1;				/* setting (and checking) requested matrix dimensions here */
	if (sq) if ((mins>=0 && mins>ROWSM(rr)) || (maxs>=0 && maxs<ROWSM(rr))
					|| (tot>=0 && tot!=2*ROWSM(rr)))  cy = 0;
	if ((mins>=0 && tot>=0 && mins>tot-ROWSM(rr)) || (maxs>=0 && tot>=0 && maxs<tot-ROWSM(rr))
					|| (mins>=0 && maxs>=0 && mins>maxs))  cy = 0;
	if (sq)  mins = maxs = ROWSM(rr);
	if (tot>=0)  mins = maxs = tot-ROWSM(rr);
	k = -1;
	while (cy) {			/* the main cycle generating column selections starts here */
		if (k>=0 && rnd>0 && rnd+k>=4) {
			cols[k] += RANDOM()%(COLSM(e)/8+k/2+rnd/2);
		}
		if (k>=0) if (++cols[k]>=COLSM(e)) {
			if (--k<0)  break;  	/* this breaks the main cycle at the end */
			continue;
		}
		DEBUG(CURDLEV+3,"\t%*s choice cols[%d]=%d\n",2*(k+1)," ",k,(k>=0?cols[k]:-1));
		if (kc>=0) {			/* request to keep the column kc in the selection */
			for (i=p=0; i<=k; i++)  if (cols[i]==kc)  p = 1;
		} else  p = 1;
		if (mrk>=0) {			/* request of maximal rank of the output submatrix */
			rc = ematrix_refer_all(rr);
			for (i=0; i<=k; i++)  ematrix_drefadd_col(rc,cols[i]);
			q = (ematrix_matrank(rc)<=mrk);
		} else { rc = NULL;  q = 1; }
		if (p && q)			/* passing the column selection for a submatrix */
		  if ((k>=mins-1 || mins<0) && (k<=maxs-1 || maxs<0)) {
			if (!rc) {		/* (rc may already be prepared above in a rank test) */
				rc = ematrix_refer_all(rr);
				for (i=0; i<=k; i++)  ematrix_drefadd_col(rc,cols[i]);
			}
			DEBUG(CURDLEV+2,"\t\t\t++ col choice (%c%c%c%c..) accepted at k=%d ++\n",(k>=0?cols[0]+'0':'.'),(k>=1?cols[1]+'0':'.'),(k>=2?cols[2]+'0':'.'),(k>=3?cols[3]+'0':'.'),k+1);
			outl = alist_append(oo=outl,rc);
			rc = NULL;
			if (oo!=outl) if (alist_getlength(outl)>EM_MAXSUBMAT) {
				PROGERROR("too many output submatrices >%d to handle",alist_getlength(outl));
				break;		/* (to limit huge memory usage) */
			}
		}
		if (rc)  dispose_ematrix(rc);
						/* when trying larger subset as a selection */
		if (k<0?1: (q && (p || cols[k]<kc) && (k<maxs-1 || maxs<0)
				&& (COLSM(e)-cols[k]+k>=mins || mins<0))) {
			++k;
			cols[k] = (k>0? cols[k-1]:-1);
		}				/* (must always go up if k<0 !) */
	}
	if (cols && cols!=colstack)  FREE(cols);
	return outl;
}




/**
 * The following function generates all (displayed) bases for the matroid given by e.
 * This means all non-singular square submatrices of e are tried,
 * and they are moreover pivoted to a new basis if piv==1.
 * The new matrices equivalent to e (piv==1) or refering to e (piv==0)
 * are returned in a list of matrix copies, in the same transposition state as e.
 * (If two distinct bases happen to produce the same matrix, then this one occurs twice.)
 * Remember to free the list and the matrices later.
 * 
 * If ch>0, some more information are printed, and extra paranoic (debug) tests performed.
 * If rnd>0, then only a random sublist of bases is returned - for quick randomized pre-checks.
 * If rnd>=4, then the generated random bases are at least 3 lines from the current one.
 * If hh[],hhb[] are given, then only those bases of e for which the element codes in hh[]
 * give the same row collection (multiset) as the row codes in hhb[] are generated.
**/

ematrix**	ematrix_getbases_ext(int ch, ematrix *e, int piv, int rnd, long hh[], long hhb[]) {
	ematrix		*ee,*rre, **x, **el,**outl;
	int		i,j,k,l, ro,co, *mh,mhstack[60];
	sign_t		gg;
	
	if (ISREFMAT(e))  {PROGERROR("Do not call for a refering matrix e=%p.",e);}
	if (rnd>0 && hh)  {PROGERROR("Do not use line codes hh[] together with randomized.");}
	if (ch>=0)  DEBUG(CURDLEV+0,"going to generate bases (piv=%d, rnd=%d) for the matrix %p (%dx%d)...\n",
			piv,rnd,e,ROWSM(e),COLSM(e));
	ro = ROWSM(e);  co = COLSM(e);
	if (2*ro+co<50)  mh = mhstack;
	else  mh = MMALLOC((2*ro+co+6)*sizeof(mh[0]));
			/**
			 * We generate all square submatrices (possibly only nonzero ones).
			 * Then we check each of them for proper codes hh,hhb if given,
			 * rejecting submatrices making different collection of row-codes than hhb.
			 * The we compute nonzero subdeterminants.
			**/
	el = ematrix_submatrices_bases(e,rnd);
	outl = NULL;
	for (x=el; x?*x:0; x++) {	/* cycling all (nonzero?) square submatrices of e: */
		if (rnd>=4 && ROWSM(*x)<=2 && ROWSM(*x)>0)
			continue;
		if (hh && hhb) {
			for (i=0; i<ro+co; i++)  mh[i] = (i<ro?1:0);
			for (i=0; i<ROWSM(*x); i++)  mh[GETREFMROW(*x,i)] = 0;
			for (i=0; i<COLSM(*x); i++)  mh[GETREFMCOL(*x,i)+ro] = 1;
#ifndef FASTPROG
			for (i=l=0; i<ro+co; i++)  if (mh[i])  l++;
			if (l!=ro)  {PROGERROR("wrong row selection ?!? %d!=%d",l,ro);}
			if (x==el && DEBUGLEV>=CURDLEV+3) {
				SDEBUG(CURDLEV+3,"\t\thh[]  = "); for (i=0; i<ro+co; i++)  SDEBUG(CURDLEV+3,"%s%ld, ",i==ro?"| ":"",hh[i]);
				SDEBUG(CURDLEV+3,"\n\t\thhb[] = "); for (i=0; i<ro+co; i++)  SDEBUG(CURDLEV+3,"%s%ld, ",i==ro?"| ":"",hhb[i]); SDEBUG(CURDLEV+3,"\n");
			}
			DEBUG(CURDLEV+3," submatrix %dx%d:\n",ROWSM(*x),COLSM(*x));
#endif
			for (i=k=0; i<ro+co; i++)  if (mh[i]==1) {	/* checking hh[]-codes: */
				for (j=i,k=0; j<ro+co; j++)
					if (mh[j] && hh[j]==hh[i]) { mh[j] = -1;  k++; }
				l = k;
				for (j=0; j<ro; j++)
					if (hhb[j]==hh[i])  k--;
				SDEBUG(CURDLEV+3,"\t\t\tcode %ld appears %d %d times -> k=%d\n",hh[i],l,l-k,k);
				if (k!=0)  break;
			}
			if (k!=0)  continue;	/* (different hh[] codes than hhb[] found in rows of *x) */
		}
		if (!piv) {
			gg = 1;
			if (COLSM(*x)>0)  ematrix_determinant(*x,NULL,&gg);
			if (gg!=0)  outl = alist_append(outl,ematrix_refer_all(*x));
			continue;
		}
			/**
			 * The refering matrix *x provides us with a square submatrix
			 * of e, that is "inverted" here to display a new basis consisting
			 * of the rows not in *x and columns in *x.
			 * This works iff det(*x)!=0 which is implicitly tested during pivoting.
			**/
		ee = ematrix_copy(e);
		rre = ematrix_refer_all(*x);  ematrix_rerefer(rre,ee);
		while (COLSM(rre)>0) {
			for (i=ROWSM(rre)-1; i>=0; i--)
				if (SIGNM(rre,i,COLSM(rre)-1)!=0)  break;
			if (i<0)  break;
			ematrix_pivot(ee,GETREFMROW(rre,i),GETREFMCOL(rre,COLSM(rre)-1));
			ematrix_remove(rre,i,COLSM(rre)-1);
		}
		j = COLSM(rre);
		dispose_ematrix(rre);
		if (j==0)  outl = alist_append(outl,ee);
		else  dispose_ematrix(ee);
	}
	if (piv && !outl && !hh)	/* (at least the current basis is always returned) */
		outl = alist_append(outl,ematrix_copy(e));
	if (el)  dispose_alist_mats(el);
	if (mh && mh!=mhstack)  FREE(mh);
#ifndef FASTPROG
	i = alist_getlength(outl);
	if (ch>=0)  DEBUG(CURDLEV+0,"total %d bases generated for %p (%dx%d)\n",i,e,ROWSM(e),COLSM(e));
					/* extra check of bases for an equivalent matrix */
	if (IFRANDDEBUGLESS(222) && ch>=0 && rnd==0 && i>2 && i<EM_MAXSUBMAT/4 && !hh) {
		ee = ematrix_copy(piv? outl[RANDOM()%i]: e);
		ematrix_transpose(ee);
		DEBUG(CURDLEV+0,"recursive check of the number of generated bases...\n");
		el = ematrix_getbases_ext(-1,ee,0,0,NULL,NULL);
		if (alist_getlength(el)!=i)  {PROGERROR("incorrect computation of bases %d!=%d",alist_getlength(el),i);}
		dispose_alist_mats(el);  dispose_ematrix(ee);
	}
#endif
	return outl;
}



/**
 * The next function tests whether the given matrix em is in the current pfield, i.e.
 * whether all of its subdeterminants (of any size) are defined in the pfield numbers.
 * The computation is done by brute force - quite slow... (But how to do better?)
 * If klr,klc>0, then only subdeterminants containing the last row or column are considered.
 * 
 * If ch>0, then an error message is printed for a non-pfield matrix,
 * and the minimal bad matrix is printed to the output (not debug).
 * So if ch<=0, the test (if) fails rather quickly, but positive answer is always slow.
 * For fast (and almost accurate) positive test, there is a randomized variant
 * of the test for ch==-7. All 2x2 subdeterminants are tested even in the random variant.
 * 
 * The function returns 0 for a pfield matrix, and -1 otherwise.
 * If one works in a pfield with total sum, then nothing is computed and 0 is returned.
**/

int	ematrix_inpfield_ext(int ch, ematrix *em, int klr, int klc) {
	int		ret,r, i,ii,j,jj, sizmin = 99999;
	ematrix		*e, **sd, **x, *estmin = NULL;
	
	if (!pfispartial) {
		DEBUG(CURDLEV+1,"NOT checking for pfield since the sum is total, not partial\n");
		return 0;
	}
	e = ematrix_copy_to(em,ROWSM(em),COLSM(em));
	DEBUG(CURDLEV+1,"checking for pfield %s the [%dx%d] matrix %p (with last %d,%d) %s...\n",pfield_curname(),ROWSM(e),COLSM(e),e,klr,klc,ch==-7?"randomized":"");
	estmin = NULL;  sd = NULL;
	klr = klr>0 ?1:0;  klc = klc>0? 1:0;
	
	ret = 0;		/* all 2x2 subdeterminants must be always checked first! */
	for (i=0; i<ROWSM(e)-1; i++)  for (j=0; j<COLSM(e)-1; j++)
		for (ii=klr? ROWSM(e)-1:i+1; ii<ROWSM(e); ii++)
		  for (jj=klc? COLSM(e)-1:j+1; jj<COLSM(e); jj++) {
			ret = ematrix_determinant_2x2check(e,i,ii,j,jj)<0? -1:ret;
			if (ret<0)  i=j=ii=jj = 9999;
	}
	if (ret>=0 || ch>0) {	/* continue with checking larger subdeterminants (if ch>0, then we want to print the smallest one) */
		sd = ematrix_submatrices_inpf(e,(ch<=-7?3:0),(ch>0?2:3),klr,klc);
				/* here ch<=-7?3:0 is the randomized submatrix selection, OK? */
		DEBUG(CURDLEV+1," - checking %d subdeterminants...\n",alist_getlength(sd));
		for (x=sd; (x?*x:0) && (ret>=0||ch>0); x++) {
			r = ematrix_determinant_check(*x);
			ret = (r<0? -1:ret);
			if (r<0 && (ch>0?1:(ch>=0&&IFRANDDEBUGLESS(222))) && ROWSM(*x)<sizmin) {
				if (estmin)  dispose_ematrix(estmin);
				estmin = ematrix_copy(*x);	/* (smaller non-pfield matrix is recorded) */
				sizmin = ROWSM(*x);
			}
		}
	}			/* handling the result... */
	DEBUG(CURDLEV+1," ...%s checking for pfield %s found %s\n",ch==-7?"randomized":"",pfield_curname(),ret<0?"-NO-":"+YES+");
	if (ch>0 && ret<0) {
		OUTPUT("Minimal non-%s submatrix of size %d found here:\n",pfield_curname(),sizmin);
		if (estmin)  EMATOUTPUT(estmin,"\tnopf\t");
	}
	
#ifndef FASTPROG		/* extra testing that the non-pfield matrix is really wrong */
	if (ret<0 && estmin && ch>-6) {		/* (estmin was set up if ch>0 or randomly if ch==0) */
		jj = ematrix_inpfield_ext(-6,estmin,0,0);	/* (ch==-6 to prevent inf recursion) */
		if (jj>=0)  {PROGERROR("incorrect computation (not)inpfield, ret=%d",ret);}
	}
	if (ret>=0 && IFRANDDEBUGLESS(222) && ch>-6 && klr==0 && klc==0 && ROWSM(e)>1 && COLSM(e)>1) {
		for (ii=jj=0; ii<4; ii++) {	/* extra testing with pivoted matrix when not random */
			i = RANDOM()%ROWSM(e);  j = RANDOM()%COLSM(e);
			if (SIGNM(e,i,j)!=0)  jj = ematrix_pivot_check(0,e,i,j)<0?-1:jj;
		}
		ematrix_transpose(e);
		jj = ematrix_inpfield_ext(-6,e,0,0)<0?-1:jj;	/* (ch==-6 to prevent inf recursion) */
		if (jj<0)  {PROGERROR("incorrect computation inpfield, ret=%d",ret);}
	}		/* (e is modified here !!!) */
#endif
	if (sd)  dispose_alist_mats(sd);
	if (estmin)  dispose_ematrix(estmin);
	dispose_ematrix(e);
	return ret<0?-1:0;
}



























/******************	Subdeterminant testing	*******************/
/******************************************************************/


#undef CURDLEV
#define CURDLEV		6


/**
 * This function tests whether the two given matrices define the same matroid.
 * (That means whether all subdeterminants have simultaneously zero/nonzero values.)
 * When xpf2>=0 is given, then the second matrix e2 is represented over another pfield xpf2.
 * When klr,klc>=0 are given, then only subdeterminants containing the row klr / column klc
 * are tested in the function.
 * 
 * The return value is 1 for identical dets, 0 for different size matrices, and -d
 * when the smallest (if ch>0) distinct subdeterminat pair is of size d.
 * If an undefined subdeterminant is found, and undf>=0, then the result is 0 (with an error
 * for undf>0); but undefined subdeterminants are considered as nonzero if undf<0.
**/

int	ematrix_samedets_ext(int ch, ematrix *e1, ematrix *e2, int xpf2, int undf, int klr, int klc) {
	sign_t		gg,gg2;
	int		d1,d2,xpf, ret;
	ematrix		*ee1,*ee2, **sd=NULL, **x;
	
	if (ROWSM(e1)!=ROWSM(e2) || COLSM(e1)!=COLSM(e2))  return 0;
	xpf = pfield_curindex();
	ee1 = ematrix_copy_ext(e1,ROWSM(e1),COLSM(e1),0);	/* (copy to 0-transp state!) */
	ee2 = ematrix_copy_ext(e2,ROWSM(e2),COLSM(e2),0);
	sd = ematrix_submatrices_dets(ee1,0,klr,klc);
	//******* randomized lists for negative ch ???
	ret = 1;
			/**
			 * We cycle through all subdeterminants given by submatrices in sd,
			 * and we compute their zero/nonzero values.
			 * (Notice that the rerefer-call is valid here - both ee1,ee2
			 *  are made to have the same max dimensions...)
			**/
	for (x=sd; x?*x:0; x++) {
		d1 = ematrix_determinant_ch(0,*x,NULL,&gg);
		if (xpf2>=0)  pfield_switchto_fast(xpf2);
		ematrix_rerefer(*x,ee2);
		d2 = ematrix_determinant_ch(0,*x,NULL,&gg2);
		if (xpf2>=0)  pfield_switchto_fast(xpf);
		
		if ((d1<0 || d2<0) && undf>=0) {
			ret = 0;
			if (undf>0)  {PROGERROR("undefined subdeterminant in pfield, size %d",ROWSM(*x));}
		}
		if (undf<0 && d1<0)  gg = 2;
		if (undf<0 && d2<0)  gg2 = 2;
		if ((gg==0)!=(gg2==0))
			if (ret>0 || ret<-ROWSM(*x))  ret = -ROWSM(*x);
		if (ch<=0 && ret<=0)  break;
	}
	
	DEBUG(CURDLEV+1,"Checking same subdeterminants in %dx%d matrices has found %d.\n",ROWSM(e1),COLSM(e1),ret);
	if (sd)  dispose_alist_mats(sd);
	if (ee1!=e1)  dispose_ematrix(ee1);
	if (ee2!=e2)  dispose_ematrix(ee2);
	return ret;
}


























/******************	Extended printing	*******************/
/******************************************************************/


#undef CURDLEV
#define CURDLEV		7


/**
 * This is a function for an extended printing of a matroid.
 * Many interesting matroid-invariant numbers are printed here, derived from the
 * list of all bases, and from the closure operator, etc...
 * The output is meant to help you with identifying and comparing matroids,
 * but it has no other meaning in this program.
 * 
 * The amount of printed information is controlled by the value of lev.
 * If lev==-1, then nothing is printed, and the hash values only is computed.
 * If bto is given, then the output (reduced) is printed to the buffer bto,
 * up to given length mx.
 * See also the next function ematrix_printbasecirc_()...
**/

#define	EM_MAXTOPRINT	(NUM_LONG_BITS/2-1)	/* (must be smaller than half long bits!) */
#define	EM_MAXTOPRINTSQ		65
#define	LIID(e,i)	((i)<ROWSM(e)? ROWSID(e,i): COLSID(e,(i)-ROWSM(e)))
int     ematrix_printmore_flats(ematrix *rf, int rk, int wh) ;

#define	SNPRINTF(b,l,m,farg...)	(snprintf(b+l,m-l-8,farg),b[m-1]=0,l+=strlen(b+l))

						/* (these are taken from struct/strmagic.c) */
long    strmag_flatlines_ext(ematrix *e, int minr, int maxr, int tf, long fl[]) ;
#define strmag_flatlines_pr(e,fl)	strmag_flatlines_ext(e,2,6,ROWSM(e)*COLSM(e)+11,fl)


long	ematrix_printmore_ext(ematrix *e, int lev, char *bto, int mx) {
	ematrix		**bas,**flt,**sep, **x, *ee,*xe;
	int		i,ii, j,jj,k,kk, a,ro,co, *li,*basi,**basii,
			 *orb=NULL,*orbx, *ses=NULL, bbm,bbl;
	long		hash = 0, *flts=NULL;
	
	if (bto)  bto[0] = 0;
	bbm = mx;  bbl=0;
	if (!e)  return -1;
	if (lev>0) if (ROWSM(e)*COLSM(e)>EM_MAXTOPRINTSQ || ROWSM(e)>EM_MAXTOPRINT || COLSM(e)>EM_MAXTOPRINT) {
		OUTPUT("Too big matroid for extended printing, sorry.\n");  return -1;
	}
	ee = ematrix_copy(e);
	ro = ROWSM(ee);  co = COLSM(ee);
	basi = MMALLOC(2*(ro+co+2)*sizeof(basi[0]));  li = basi+ro+co;
	basii = malloc_twodim(sizeof(basii[0][0]),ro+co,ro+co);
	for (i=0; i<ro+co; i++)  basi[i] = 0;
	for (i=0; i<ro+co; i++) for (ii=0; ii<ro+co; ii++)  basii[i][ii] = 0;
	flts = MMALLOC((ro+co+2)*sizeof(flts[0]));
	if (lev>0 && !bto)  SOUTPUT("\n");
	hash = ro+10*co;
	hash = 90909l*hash+111l*hash*hash*hash;

		/**
		 * Here we collect and print all bases - their total number, and numbers
		 * per elements (basi[]) and per element pairs (basii[][]).
		 * The bases are given as refering square submatrices of ee (not pivoted!).
		 * We add the numbers to the matroid hash value in a symmetric way.
		 *   ....... more to be added to the hash value - small flatlines???
		 * 
		 * The bases are also used later, and they are freed at the end.
		**/
	bas = ematrix_getbases_sq(ee);
	hash += alist_getlength(bas)*10101l;
	for (x=bas; x?*x:0; x++) {
		xe = ematrix_refextract_xrow(ee,*x);
		if (ROWSM(xe)+COLSM(xe)!=ro)  {PROGERROR("Something wrong with the basis size!");}
		for (i=0; i<ROWSM(xe); i++)
			li[i] = GETREFMROW(xe,i);
		for (i=0; i<COLSM(xe); i++)
			li[i+ROWSM(xe)] = GETREFMCOL(xe,i)+ro;
		for (i=0; i<ro; i++) {
			basi[li[i]]++;
			for (ii=0; ii<ro; ii++)  basii[li[i]][li[ii]]++;
		}
		dispose_ematrix(xe);
	}
	for (i=0; i<ro+co; i++) {
		hash += 505l*basi[i]+7l*basi[i]*basi[i];
		for (ii=0; ii<i; ii++)  hash += 13l*basii[i][ii]+1l*basii[i][ii]*basii[i][ii];
		if (basii[i][ii]!=basii[ii][i])  {PROGERROR("Something wrong - nonsymmetric basis pairs!");}
	}
	//************* adding small flatlines to the hash value? - here, not below, since the element
	//		magic below is slow to compute and not computed always(!)
	
	if (lev>=0) {
		if (!bto)  OUTPUT("Number of matroid [%.25s] bases:  %d\n",EMNAME(e)?EMNAME(e):"",alist_getlength(bas));
		else  SNPRINTF(bto,bbl,bbm,"Matroid  %d x %d [%.25s],  %d bases.\n",ro,co,EMNAME(e)?EMNAME(e):"",alist_getlength(bas));
	}
	if (lev>0 && !bto) {
		OUTPUT("  - per elements ");
		for (i=0; i<ro+co; i++)  SOUTPUT(" [%d: %d]%s",LIID(ee,i),basi[i],i%6==5?"\n\t\t\t":"");
		SOUTPUT("\n");
	}
	if (lev>4 && !bto) {
		SOUTPUT("\n");  OUTPUT("  - per element pairs\n");
		for (i=0,a=1; i<ro+co; i++) {
			if (a)  SOUTPUT("    ~\t ");  a = 0;
			for (ii=0; ii<ro+co; ii++)  if (basii[i][ii]>0) {
				a = 1;
				SOUTPUT("[%2d'%2d: %-3d]",LIID(ee,i),LIID(ee,ii),basii[i][ii]);
			}
			if (a)  SOUTPUT("\n");
		}
	}
		/**
		 * Here we (try to) distinguish matroid elements up to isomorphism.
		 * We either compute "magic values" for the elements - faster,
		 * or we rigorously compute the orbits of the automorphism group.
		 * For the orbits, we use the numbers of bases in basi[] and the
		 * flatline values flts[] for rough distinction, and then we use
		 * strmag_isautmap() to see which pairs of elements are really
		 * mapped to each other by the aut group.
		**/
	if (lev>0) {
		strmag_flatlines_pr(ee,flts);	/* (flts[] is globally allocated) */
		for (i=0; i<ro+co; i++)  flts[i] += 7*basi[i];
	}
	if (lev>0 && lev<=2 && !bto) {
		OUTPUT("  - elem magic ");
		for (i=0; i<ro+co; i++)  SOUTPUT(" [%d: %ld]%s",LIID(ee,i),flts[i],i%6==5?"\n\t\t\t":"");
		SOUTPUT("\n");
	}	
	if (lev>2 && ro>0) {
		if (ro+co>6)  DEBUG(CURDLEV-3,"Warning - aut orbit computation may take very long.\n");
		if (!bto)  OUTPUT("Aut group orbits of [%.25s] are (via first elem id):\n",EMNAME(e)?EMNAME(e):"");
		else  SNPRINTF(bto,bbl,bbm,"Aut group orbits ");
		if (ematrix_checkid(ee)<0 && !bto) {
			ematrix_resetid(ee);  SOUTPUT("\t\tResetting line id's in the matrix!!!\n");
		}
		orb = MMALLOC(2*(ro+co+2)*sizeof(orb[0]));
		orbx = orb+(ro+co+1);
		for (i=0; i<ro+co; i++)  orbx[i] = 1;
		if (!bto) {
			for (i=0; i<ro; i++)  orb[i] = ROWSID(ee,i);
			for (i=ro; i<ro+co; i++)  orb[i] = COLSID(ee,i-ro);
		} else  for (i=0; i<ro+co; i++)  orb[i] = i;
			/* (we print elem ids on output, but their indices to bto...) */
		if (!bto)  SOUTPUT("\t\t(%d",orb[0]);
		else  SNPRINTF(bto,bbl,bbm,"[%d",orb[0]);
		kk = ro+co;
		for (i=1; i<ro+co; i++) {
			for (j=0; j<i; j++) if (orbx[j]) {
				if (basi[i]!=basi[j] || flts[i]!=flts[j])  continue;
				//********** may we efficiently use basii[][] here???
				if (!strmag_isautmap_h(e,i,j,flts))  continue;
				orb[i] = orb[j];  orbx[i] = 0;  kk--;
			}
			if (!bto)  SOUTPUT(",%s %d",i==ro?" ":"",orb[i]);
			else  SNPRINTF(bto,bbl,bbm,",%s %d",i==ro?" ":"",orb[i]);
		}
		if (!bto)  SOUTPUT(") =%d\n",kk);
		else  SNPRINTF(bto,bbl,bbm,"] %d.\n",kk);
		FREE(orb);
	}
		/**
		 * Here we list all nontrivial flats of the matroid up to rank
		 * depending on lev, or up to the first flats found.
		 * See also ematrix_printmore_flats() below.
		 * (We currently do not use information about the printed flats
		 * in the matroid hash-value since the flats are expensive to compute.
		 * Also, the flat comp implementation works only for bounded size!)
		**/
	ematrix_printmore_flats(NULL,0,1);
	if (lev>0)
	  for (k=0,a=1; k<lev+a && k<ro; k++) {
		flt = ematrix_submatrices_sub(ee,k);
		for (x=flt,i=0; x?*x:0; x++) {
			xe = ematrix_closure(ee,*x);
			kk = ematrix_setrank(ee,xe);
			if (kk>k)  {PROGERROR("something wrong with the flat rank %d>%d",kk,k);}
			if (ROWSM(xe)+COLSM(xe)>k && k==kk)
			  if (ematrix_printmore_flats(xe,k,0)>=0) {
				if (i==0 && lev>1 && !bto)  SOUTPUT("\n");
				if (i==0 && !bto)  OUTPUT("Listing all (nontrivial) flats in [%.25s] of rank %d:\n",
							EMNAME(e)?EMNAME(e):"",k);
				if (i==0 && bto)  SNPRINTF(bto,bbl,bbm,"Flats of rank %d:",k);
				i++;
				if (i%5==0 && bto)  SNPRINTF(bto,bbl,bbm,"\n");
				if (!bto) {
				  SOUTPUT("    ~\t - rank-%d flat (%d)\t{",k,i);
				  for (j=0; j<ROWSM(xe); j++)  SOUTPUT(" %d,",ROWSID(xe,j));
				  for (j=0; j<COLSM(xe); j++)  SOUTPUT("%c %d",!j?' ':',',COLSID(xe,j));
				  SOUTPUT(" }\n");
				} else {
				  SNPRINTF(bto,bbl,bbm,"  %df{",i);
				  for (j=0; j<ROWSM(xe); j++)  SNPRINTF(bto,bbl,bbm,"%d,",GETREFMROW(xe,j));
				  for (j=0; j<COLSM(xe); j++)  SNPRINTF(bto,bbl,bbm,"%c%d",!j?' ':',',GETREFMCOL(xe,j)+ro);
				  SNPRINTF(bto,bbl,bbm,"}");
				}
			}
			dispose_ematrix(xe);
		}
		if (i==0 && k==0 && lev>1 && !bto)  SOUTPUT("\n");
		if (i==0 && !bto)  OUTPUT("There are -NO- (nontrivial) flats in [%.25s] of rank %d.\n",EMNAME(e)?EMNAME(e):"",k);
		if (i>0 && bto)  SNPRINTF(bto,bbl,bbm,"\n");
		if (flt)  dispose_alist_mats(flt);
		if (i!=0 && bto && k>=3)  break;
		if (i==0 && k>0)  a++;
	}
	ematrix_printmore_flats(NULL,0,-1);
	
		/**
		 * Here we list all nontrivial separations of the matroid up to lambda
		 * depending on lev, or up to the first separations found.
		 * (We currently do not use these information in the hash-value.)
		**/
	if (lev>1 && !bto) {
		sep = ematrix_submatrices_all(ee);
		ii = sep? alist_getlength(sep):0;
		ses = MMALLOC((ii+2)*sizeof(ses[0]));
		for (x=sep; x?*x:0; x++)
			ses[x-sep] = ematrix_whatsep(ee,*x);
		for (k=0,a=-1; k<lev+a && k<ro; k++) {
			for (jj=i=0; jj<ii; jj++) if (ses[jj]==k) {
				xe = sep[jj];
				if (ROWSM(xe)+COLSM(xe)<=k || ROWSM(xe)+COLSM(xe)>(ro+co)/2)
					continue;
				if (i==0 && lev>1)  SOUTPUT("\n");
				if (i==0)  OUTPUT("Listing all exact separations in [%.25s] of lambda %d:\n",
							EMNAME(e)?EMNAME(e):"",k+1);
				SOUTPUT("    ~\t - %d-separation (%d)\t(",k+1,++i);
				for (j=0; j<ROWSM(xe); j++)  SOUTPUT(" %d,",ROWSID(xe,j));
				SOUTPUT(" ");
				for (j=0; j<COLSM(xe); j++)  SOUTPUT(" %d,",COLSID(xe,j));
				SOUTPUT(" )\n");
			}
			if (i==0 && k==0)  SOUTPUT("\n");
			if (i==0)  OUTPUT("There are -NO- exact separations in [%.25s] of lambda %d.\n",EMNAME(e)?EMNAME(e):"",k+1);
			if (i==0 && k>0)  a++;
		}
		if (sep)  dispose_alist_mats(sep);
		FREE(ses);
	}
	
#ifndef FASTPROG		/* paranoic testing of the hash-value computation: */
	if (lev>=0 && IFRANDDEBUGLESS(111)) {
		for (ii=0; ii<4; ii++) {	/* extra testing hash with pivoted matrix */
			i = RANDOM()%ROWSM(ee);  j = RANDOM()%COLSM(ee);
			if (SIGNM(ee,i,j)!=0)  ematrix_pivot(ee,i,j);
		}
		if (hash!=ematrix_printmore_ext(ee,-1,NULL,0))  {PROGERROR("incorrect computation of matroid hash, ret=%ld",hash);}
	}		/* (ee is modified here !!!) */
#endif
	if (bas)  dispose_alist_mats(bas);
	if (flts)  FREE(flts);
	if (basi)  FREE(basi);	if (basii)  FREE(basii);
	
		/**
		 * Some other final characteristics are printed here.
		 * Representability is surveyed for some basic fields, depending on lev.
		 * Among them the matroid hash-value is printed out and always returned.
		 * So far, the matroid hash-value collects information about rank,
		 * number of bases, numbers of bases per each element and each pair of elements.
		 * You must update the version number if you change the collected information!
		**/
	if (lev>0 && !bto) {
		SOUTPUT("\n");  k = 0;
		OUTPUT("Matroid [%.25s] connectivity is %d",
				EMNAME(e)?EMNAME(e):"",kk=struct_iconnectivity(e,&k));
		if (kk==3 && k)  SOUTPUT(" (internally %d-connected).\n",kk+1);
		else  SOUTPUT(".\n");
		OUTPUT("Matroid [%.25s] girth (shortest cycle) is %d.\n",
				EMNAME(e)?EMNAME(e):"",struct_matgirth(e));
	}
	if (lev>0 && bto) {
		SNPRINTF(bto,bbl,bbm,"Connectivity  %d, girth (shortest cycle)  %d.\n",
				struct_connectivity(e),struct_matgirth(e));
	}
	if (lev>2) {
		char	*xx_fields[] = {"GF(2)","GF(3)","GF(4)","GF(5)","GF(7)","GF(8)","GF(9)"};
		if (!bto)  OUTPUT("Matroid [%.25s] representability:",EMNAME(e)?EMNAME(e):"");
		else  SNPRINTF(bto,bbl,bbm,"Representability over:");
		for (j=0; j<(int)(sizeof(xx_fields)/sizeof(xx_fields[0])); j++) {
			ii = pfield_curindex();
			pfield_switchto(xx_fields[j]);	/* (we expect all the fields defined!) */
			a = grepr_isrepresented(ee,ii)>0;
			if (!bto)  SOUTPUT(" %c%s%c",a?'+':'-',xx_fields[j],a?'+':'-');
			else  SNPRINTF(bto,bbl,bbm," %c%s%c",a?'+':'-',xx_fields[j],a?'+':'-');
			pfield_switchto_fast(ii);
		}
		if (!bto)  SOUTPUT("\n\n");  else  SNPRINTF(bto,bbl,bbm,"\n");
	}
	if (lev>=0 && !bto)  OUTPUT("Overall matroid [%.25s] hash-value (version %s):  %ld\n",
					EMNAME(e)?EMNAME(e):"",EM_HASHVER,hash);
	dispose_ematrix(ee);
	return hash;
}


/**
 * This is a supplementary function used to prevent repetition of flats, and trivial
 * flats obtained just by adding free points to smaller flats.
**/

int	ematrix_printmore_flats(ematrix *rf, int rk, int wh) {
	static long	*save=NULL, sl=0,si;
	ematrix		*rrf;
	long		l;
	int		i,j,k, r = 1;
	
	if (wh!=1) if (!save || !sl)  {PROGERROR("must be initialized first!"); return 0;}
	if (wh==1) {
		sl = 1000;  save = MMALLOC(sl*sizeof(save[0]));
	} else if (wh==-1) {
		FREE(save);  save = NULL;  si = sl = 0;
	} else {
		if (si>=sl-4)  save = MREALLOC(save,(sl*=2)*sizeof(save[0]));
		rrf = REFMAT(rf);
			/**
			 * The previously allowed flats are stored as set bitmaps to the list
			 * save, and the new flat given by rf is compared against all saved
			 * ones bit by bit.
			 * If it equals, or it is just by the rank difference bigger, to saved
			 * one, the new flat is rejected (and not stored).
			**/
		for (i=0,l=0; i<ROWSM(rf)+COLSM(rf); i++)
			l |= 1l<<(i<ROWSM(rf)? GETREFMROW(rf,i):
					GETREFMCOL(rf,i-ROWSM(rf))+ROWSM(rrf));
		for (i=0; i<si/2 && r>=0; i++) {
			for (j=k=0; j<ROWSM(rrf)+COLSM(rrf); j++)
				k += ((l&(1l<<j))!=(save[2*i]&(1l<<j)));
			if ((k==1 && rk==save[2*i+1]) || (k==0 && rk!=save[2*i+1]))
				{PROGERROR("something wrong with the saved flats... (k=%d, %d~%ld)",k,rk,save[2*i+1]);}
			if (k<=rk-save[2*i+1])  r = -1;
		}
		if (r>=0) {
			save[si] = l;  save[si+1] = rk;
			si += 2;
		}
	}
	return r;
}




/**
 * This function is used to print out all matroid bases and all circuits.
 * The matroid is given in e, parameter whp determines the function - <=10 for printing bases,
 * >=10 for printing circuits, and lev is the verbosity printing level.
 * One may possibly specify elements (via their id's) that must be contained in the printed
 * bases/circuits, in the array cix[], where NULL or 0 value mean no element.
 * If bto is given, then the output (reduced) is printed to the buffer bto,
 * up to given length mx.
**/

void	ematrix_printbasecirc_ext(ematrix *e, int whp, int lev, int *cix, char *bto, int mx) {
	
	ematrix		**bas, **be, *ee;
	int		i=0,j,jj,k,n, ro,co, blen,clen;
	long		l, *crc=NULL;
	
	if (bto)  {PROGERROR("Not implemented yet!"); return;}
	if (!e)  return;
	if (lev>0 && !bto)  SOUTPUT("\n");
	if (whp>1) if (ROWSM(e)>EM_MAXTOPRINT || COLSM(e)>EM_MAXTOPRINT) {
		OUTPUT("Too big matroid for extended printing, sorry.\n");  return;
	}
	ee = ematrix_copy(e);
	ro = ROWSM(ee);  co = COLSM(ee);
	if (ematrix_checkid(ee)<0 || ematrix_checkid(ee)>EM_MAXTOPRINT) {
		ematrix_resetid(ee);  SOUTPUT("\t\tResetting line id's in the matrix!!!\n");
	}
	bas = ematrix_getbases(ee);
	blen = alist_getlength(bas);
	if (!bto)  OUTPUT("Number of matroid [%.25s] bases:  %d\n",EMNAME(e)?EMNAME(e):"",blen);
	if (cix && !bto)  for (i=0; i<ro+co && cix[i]!=0; i++) {
		if (i==0)  OUTPUT("Listing all matroid [%.25s] %s containing elements (id): ",EMNAME(e)?EMNAME(e):"",whp<=10?"bases":"circuits");
		//if (cix[i]<-2*EM_MAXTOPRINT || cix[i]>2*EM_MAXTOPRINT)  cix[i] = 0;
		SOUTPUT(" %d,",cix[i]);
	}
	if (cix && !bto && i>0)  SOUTPUT("\n");
	        /**
	         * This part is used to print all bases (by their line ids).
	         * We test that the (possible) required elements from cix[] are there.
	        **/
	if (whp<=10) for (be=bas,k=0; be?*be:0; be++) {
		if (cix)  for (i=0; cix[i]; i++) {
			for (j=0; j<ro; j++)  if (ROWSID(*be,j)==cix[i])  break;
			if (j>=ro)  break;
		}
		if (cix?cix[i]:0)  continue;
		SOUTPUT("    ~\t - base (%d)\t{",++k);
		for (j=0; j<ro; j++)  SOUTPUT(" %d%s",ROWSID(*be,j),(j<ro-1?",":""));
		SOUTPUT(" }\n");
	}
	        /**
	         * This part is used to print all circuits (by their line ids).
	         * We test that the (possible) required elements from cix[] are there.
	         * Then we test that the same circuit (as a set) has not been printed
	         * out before - we store all printed circuits as bitmaps in crc[].
	        **/
	if (whp>=10) {
		if (!cix)  OUTPUT("Listing all matroid [%.25s] circuits:\n",EMNAME(e)?EMNAME(e):"");
		crc = MMALLOC(blen*co*sizeof(crc[0]));
		clen = 0;
	  for (be=bas,k=0; be?*be:0; be++) for (jj=0; jj<co; jj++) {
		if (cix)  for (i=0; cix[i]; i++) {
			for (j=0; j<ro; j++)
				if (SIGNM(*be,j,jj) && ROWSID(*be,j)==cix[i])  break;
			if (j>=ro && COLSID(*be,jj)!=cix[i])  break;
		}
		if (cix?cix[i]:0)  continue;
				/* preventing duplicated circuits (stored as bitmaps here) */
		l = 1l<<(COLSID(*be,jj)+EM_MAXTOPRINT);
		for (i=0; i<ro; i++)  if (SIGNM(*be,i,jj))
			l |= 1l<<(ROWSID(*be,i)+EM_MAXTOPRINT);
		for (j=0; j<clen; j++)  if (l==crc[j])  break;
		if (j<clen)  continue;
		crc[clen++] = l;
		n = 0;		/* printing... */
		SOUTPUT("    ~\t - circuit (%d)\t{",++k);
		for (i=0; i<ro; i++) if (SIGNM(*be,i,jj))
			{ n++;  SOUTPUT(" %d,",ROWSID(*be,i)); }
		SOUTPUT(" %d }  \tlen %d,\n",COLSID(*be,jj),n+1);
		
	}}
	if (lev>1 && !bto)  SOUTPUT("\n");
	if (crc)  FREE(crc);
	if (bas)  dispose_alist_mats(bas);
	dispose_ematrix(ee);
}

































