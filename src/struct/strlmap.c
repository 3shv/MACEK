
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
 * Some theory for start... about matrix line maps:
 * (Read general description of "structural" functions on the top of ../include/struct.h.)
 * 
 * We want to produce all (unscaled or scaled) injective line maps of the given matrix E to E'.
 * (that means injective maps from rows/columns of E to rows/columns of E').
 * To state that VERY clearly, an entry Eij at (i,j) is mapped by a line map to an entry
 * Eij*sr[pr[i]]*sc[pc[j]] at E'(pr[i],pc[j]), where pr,pc are the row and column maps
 * and sr,sc are the row and column scales (scales indexed by the destination matrix indices!).
 * Hence Eij*sr[pr[i]]*sc[pc[j]]==E'(pr[i],pc[j]) for all choices of i,j in E.
 * Lines of E' not covered by the map may be arbitrary, but sr[],sc[] is 0 on uncovered lines.
 * 
 * These maps are useful in describing matrix "automorphisms", displayed minors,
 * and similar structural properties that are tied with a matroid representation
 * against a particular basis.
 * There is a possibility of "recycling" the line maps, but this should be used with caution.
 * 
 * Generaly, the functions here are not supposed for public use, but only for other
 * structural functions.
 * Look at strminor.c as well...
**/




#include "macek.h"
#include "str.h"













/******************	Line map structure	*******************/
/******************************************************************/


#undef CURDLEV
#define CURDLEV		6		/* (need <=DEBUGLEV+2 to perform paranoic tests) */



/**
 * This function allocates and returns one structure for describing a line map of a matrix.
 * The whole structure uses one chunk of memory, so it is enough to "free(a)".
 * Read more in str.h.
**/

emlinemap*	strmap_new_ext(ematrix *es, ematrix *ed, int init) {
	emlinemap	*a;
	int		sr,sc,dr,dc, i,k,t;
	
	sr = ROWSM(es);  sc = COLSM(es);
	dr = ROWSM(ed);  dc = COLSM(ed);
	k = MEMALIGNSIZE(sizeof(*a));
	k += MEMALIGNSIZE(sr*sizeof(a->prc[0][0]))+MEMALIGNSIZE(sc*sizeof(a->prc[0][0]));
	k += MEMALIGNSIZE(dr*sizeof(a->xrc[0][0]))+MEMALIGNSIZE(dc*sizeof(a->xrc[0][0]));
	k += MEMALIGNSIZE(dr*sizeof(a->grc[0][0]))+MEMALIGNSIZE(dc*sizeof(a->grc[0][0]));
	k += sizeof(int);	/* (for a magic number at the end) */
	a = CALLOC(k+4,1);
	a->sz = k;  EMLMMAGIC(a) = EMLMMAGICNUM;
	
	t = ISTRANSPM(es);
	a->prc[t] = (short*)MEMALIGNPOINT(a+1);
	a->prc[!t] = (short*)MEMALIGNPOINT(a->prc[t]+sr);
	a->xrc[t] = (exp_t*)MEMALIGNPOINT(a->prc[!t]+sc);
	a->xrc[!t] = (exp_t*)MEMALIGNPOINT(a->xrc[t]+dr);
	a->grc[t] = (sign_t*)MEMALIGNPOINT(a->xrc[!t]+dc);
	a->grc[!t] = (sign_t*)MEMALIGNPOINT(a->grc[t]+dr);
	if ((long)(a->sz-sizeof(int))!=(long)((byte1*)(a->grc[!t])-(byte1*)a+MEMALIGNSIZE(dc*sizeof(a->grc[0][0]))))
		{PROGERROREXIT("something wrong with allocation of emlinemap %ld %ld",(long)(a->sz-sizeof(int)),(long)((byte1*)(a->grc[1])-(byte1*)a)+MEMALIGNSIZE(dc*sizeof(a->grc[0][0])));}
	
	a->xtr = (ISTRANSPM(es)!=ISTRANSPM(ed));
	EMLMSR(a,es) = sr;  EMLMSC(a,es) = sc;
	EMLMDR(a,es) = dr;  EMLMDC(a,es) = dc;
	if (init>0) {
		for (i=0; i<sr; i++) {
			EMLMPR(a,es,i) = i;
			pfield_setzeroexp(&EMLMXR(a,es,i));  EMLMGR(a,es,i) = 1;
		}
		for (i=0; i<sc; i++) {
			EMLMPC(a,es,i) = i;
			pfield_setzeroexp(&EMLMXC(a,es,i));  EMLMGC(a,es,i) = 1;
		}
	}
	return a;
}

emlinemap*	strmap_copy(emlinemap *a) {
	emlinemap	*an;
	long		df;
	
	an = MMALLOC(a->sz+4);
	MEMCPY(an,a,a->sz);
	df = (byte1*)an-(byte1*)a;
	an->prc[0] = (short*)((byte1*)(a->prc[0])+df);  an->prc[1] = (short*)((byte1*)(a->prc[1])+df);
	an->xrc[0] = (exp_t*)((byte1*)(a->xrc[0])+df);  an->xrc[1] = (exp_t*)((byte1*)(a->xrc[1])+df);
	an->grc[0] = (sign_t*)((byte1*)(a->grc[0])+df);  an->grc[1] = (sign_t*)((byte1*)(a->grc[1])+df);
	return an;
}


/**
 * This function checks the given line map a against the matrices es,ed.
 * Errors are printed if ch>0, -1 is returned for incorrect a.
 * For ch>=3, the line map of es is computed and compared with the submatrix of ed.
**/

int	strmap_check_ext(int ch, emlinemap *a, ematrix *es, ematrix *ed) {
	int	i,j, p;
	exp_t	xx;
	sign_t	gg;
	
	if (a? (EMLMMAGIC(a)!=EMLMMAGICNUM):1)  {PROGERROREXIT("wrong or damaged emlinemap structure %p",a); return -1;}
	if (a->xtr!=(ISTRANSPM(es)!=ISTRANSPM(ed)))  {if (ch>0) PROGERROR("wrong transposition of the destination matrix for emlinemap %p",a); return -1;}
	if (EMLMSR(a,es)!=ROWSM(es) || EMLMSC(a,es)!=COLSM(es) || EMLMDR(a,es)!=ROWSM(ed) || EMLMDC(a,es)!=COLSM(ed))
		{if (ch>0) PROGERROR("wrong matrix dimensions for emlinemap %p",a); return -1;}
	if (ch<0 || (ch<=1 && !IFRANDDEBUG(55)))  return 0;
	p = 0;
	for (i=0; i<ROWSM(es); i++)  if (EMLMPR(a,es,i)<0 || EMLMPR(a,es,i)>=ROWSM(ed))  p = 1;
	for (i=0; i<COLSM(es); i++)  if (EMLMPC(a,es,i)<0 || EMLMPC(a,es,i)>=COLSM(ed))  p = 1;
	if (p)  {if (ch>0) PROGERROR("line mapping out of bounds in emlinemap %p",a); return -1;}
	//******** check an injection as well ???
	p = 0;
	for (i=0; i<ROWSM(es); i++)  if (EMLMGR(a,es,EMLMPR(a,es,i))==0)  p = 1;
	for (i=0; i<COLSM(es); i++)  if (EMLMGC(a,es,EMLMPC(a,es,i))==0)  p = 1;
	if (p)  {if (ch>0)  PROGERROR("line mapping zero-scaled in emlinemap %p",a);
			EMLMAPDEBUGS(0,a,es,"\t!!\t"); return -1;}
	if (ch<3)  return 1;
	
	for (i=0; i<ROWSM(es); i++)  for (j=0; j<COLSM(es); j++) {
		xx = EXPM(es,i,j);  gg = SIGNM(es,i,j);
		pfield_mul(1,xx,gg, 1,EMLMXR(a,es,EMLMPR(a,es,i)),EMLMGR(a,es,EMLMPR(a,es,i)), &xx,&gg);
		pfield_mul(1,xx,gg, 1,EMLMXC(a,es,EMLMPC(a,es,j)),EMLMGC(a,es,EMLMPC(a,es,j)), &xx,&gg);
		if (!pfield_isequal(xx,gg,EXPM(ed,EMLMPR(a,es,i),EMLMPC(a,es,j)),SIGNM(ed,EMLMPR(a,es,i),EMLMPC(a,es,j))))
			{if (ch>0)  PROGERROR("line mapping bad entry %d,%d in emlinemap %p",i,j,a); return -1;}
	}
	return 1;
}


/**
 * This function prints the given line map a to the output fout with a prefix pref.
 * 
**/

void	strmap_fprint(FILE *fout, emlinemap *a, ematrix *e, char *pref) {
	int	i;
	
	fprintf(fout,"%sLine map %p of %dx%d to %dx%d:\n",pref,a,EMLMSR(a,e),EMLMSC(a,e),EMLMDR(a,e),EMLMDC(a,e));
	fprintf(fout,"%s row",pref);
	for (i=0; i<EMLMSR(a,e); i++)  fprintf(fout,"\t%d->%d (*%s) ",i,EMLMPR(a,e,i),
				pfield_pvalue(7,EMLMXR(a,e,EMLMPR(a,e,i)),EMLMGR(a,e,EMLMPR(a,e,i))));
	fprintf(fout,"\n%s col",pref);
	for (i=0; i<EMLMSC(a,e); i++)  fprintf(fout,"\t%d->%d (*%s) ",i,EMLMPC(a,e,i),
				pfield_pvalue(7,EMLMXC(a,e,EMLMPC(a,e,i)),EMLMGC(a,e,EMLMPC(a,e,i))));
	fprintf(fout,"\n");
}

void	struct_lmapfprint(FILE *fout, void *a, ematrix *e, char *pref) {
	if (!pref)  pref = printoutpref;
	strmap_fprint(fout,a,e,pref);
}


/**
 * This function translates the given line map a (es->ed) into a refering submatrix to ed.
 * The map a is stored in the resulting refering matrix, so it must NOT be freed separately.
 * The new matrix is returned.
**/

ematrix*	strmap_tosubmatrix(ematrix *es, ematrix *ed, emlinemap *a) {
	ematrix		*re;
	int		i;
	
	if (ISREFMAT(ed))  {PROGERROR("do not call for a refering matrix ed!");}
	re = ematrix_refer_empty(ed);
	EMDATA(re) = a;
	if (strmap_check_size(a,es,ed)<0)  return re;
	for (i=0; i<ROWSM(es); i++)  ematrix_refadd_row(re,EMLMPR(a,es,i));
	for (i=0; i<COLSM(es); i++)  ematrix_refadd_col(re,EMLMPC(a,es,i));
	return re;
}




















/******************	Restricted matrix line maps	*******************/
/**************************************************************************/


#undef CURDLEV
#define CURDLEV		6		/*6 (need <=DEBUGLEV+2 to perform paranoic tests) */


/**
 * This function is called to find injective line maps of columns of the matrix efr to columns
 * of the matrix eto (unscaled mappings if scal==0) that respect fixed rows of these matrices.
 * The matrices to map are in efr, etox - they need not have the same number of columns,
 * but they must have the same number of rows.
 * Alternatively, etox is a larger matrix from which we select the rows given in the row
 * map ax if given (and this is reflected in the output maps in *al).
 * The arrays hcfr[],hcto[] give column hash codes - if a column j of efr can be mapped to
 * a column jj of eto, then it must be hcfr[j]==hcto[jj].
 * (Code 0 is not tested in the mapping.)
 * 
 * Line maps are consider in the sense described in  str.h  and above.
 * If al is given, then all the maps found here are appended to the list *al.
 * (The list must be initialized prior to calling this function!)
 * If the list afo is given, then its col map is supposed to be a permutation group
 * on the columns of efr, and these permutations are factored out from all col maps here.
 * 
 * The return value is -1 for no mapping found, or 0 for existing mappings but al==NULL,
 * or the number of mappings found and stored in the list *al.
 * If ch>=3, then all mappings are found and printed out.
 * The parameter ll gives the number of previously printed maps for printing results.
**/

int	strmap_fixedrows_ext(int ch, ematrix *efr, ematrix *etox, int scal, emlinemap **afo,
				long hcfr[], long hcto[], emlinemap *ax, emlinemap ***al, int ll) {
	emlinemap	*a=NULL, **x;
	long		*zpfr=NULL,*zpto, zpstack[100];
	int		*ccfr, *mpto, *rosc,*cosc;
	short		*cmap=NULL;
	ematrix		*eto, *esc, **escal=NULL,*escstack[20];
	exp_t		xxc,xxr;
	sign_t		ggc,ggr;
	int		i,j,k, p,r, cfr,cto,ro;
	extern int	sm_callfixrows, sm_callfixrows2;	/* (for debug statistics only) */
	
	if ((efr&&etox)? (ROWSM(efr)!=ROWSM(etox) && !ax):1)
		{PROGERROR("both efr,eto must be given, on the same rows"); return -1;}
	cfr = COLSM(efr);  cto = COLSM(etox);
	ro = ROWSM(efr);
	if (afo?(alist_getlength(afo)<=1):0)  afo = NULL;
#ifndef FASTPROG
	if (ax? (strmap_check_size(ax,efr,etox)<0):0)
		{PROGERROR("given map ax must agree with the dimensions of efr,etox"); return -1;}
	if (afo?(afo[0]&&afo[1]):0) {
		k = RANDOM()%alist_getlength(afo);  i = RANDOM()%ROWSM(efr);
		if (strmap_check_size(afo[k],efr,efr)<0 || EMLMPR(afo[k],efr,i)!=i)
			{EMATDEBUGS(0,efr," !\t"); EMLMAPDEBUGS(0,afo[k],efr,"   !\t"); afo = NULL;
			 PROGERROR("given factor-out maps in afo must be row-fixing automorphisms of efr");}
	}
	DEBUG(CURDLEV+3,"Looking for %s-rows %sline maps of %p (%dx%d) to %p (%dx%d).\n",
			ax?"given":"fixed",scal?"":"unscaled ",efr,ro,cfr,etox,ROWSM(etox),cto);
	if (afo)  DEBUG(CURDLEV+3," (Factoring-out list of %d maps afo.)\n",alist_getlength(afo));
	EMATDEBUGS(CURDLEV+3,efr,"\t\t\t\tfr<\t");  EMATDEBUGS(CURDLEV+3,etox,"\t\t\t\tto>\t");
	if (ax)  EMLMAPDEBUGS(CURDLEV+3,ax,efr,"\t\t\tax-\t");
#endif
	sm_callfixrows++;
	r = 0;
	if (cfr>cto || ((hcfr&&hcto)? !strmag_isinjection(hcfr,cfr,hcto,cto):0))
		r = -1;
			/**
			 * First we look at the destination matrix etox - if it is on 
			 * the same rows as efr (for ax==NULL), then we rename it to eto,
			 * otherwise we extract the rows of etox according to ax row map.
			 * We look at the zero patterns of the columns, whether they
			 * allow an injective fixed-row mapping from efr to eto.
			 * If not, we are finished here.
			**/
	if (ax && r>=0) {
		eto = ematrix_refer(etox,-1,-1,0,cto);
		for (i=0; i<ro; i++)
			ematrix_refadd_row(eto,EMLMPR(ax,efr,i));
		EMATDEBUGS(CURDLEV+3,eto,"\t\t\t\tto>-\t");
	} else  eto = etox;
	if (cfr+cto+cfr+cto+ro+cfr<96)  zpfr = zpstack;
	else  zpfr = MMALLOC((cfr+cto)*sizeof(zpfr[0])+(2*cfr+cto+ro+5)*sizeof(ccfr[0]));
	zpto = zpfr+cfr;
	if (r>=0) {
		strmag_zeropattern_cols(efr,zpfr);
		strmag_zeropattern_cols(eto,zpto);
		DEBUG(CURDLEV+4,"Checking zero-patterns in the matrix columns.\n");
		if (!strmag_isinjection(zpfr,cfr,zpto,cto))  r = -1;
	}
			/**
			 * Then we find the connected order of columns of efr in ccfr[].
			 * Based on this order, we determine when are which rows/cols scaled
			 * when fitting columns later:
			 *  rosc[i] holds the first nonzero index in row i in this order,
			 *  cosc[j] holds the first intersection of j-th column in this order
			 *   with a row that has previous nonzero, or -1 for nothing.
			 * (Notice that both rosc,cosc index columns by the connected order,
			 * not by their positions in efr!)
			 * We also initialize the line map a (against etox!) for later copying,
			 * and the matrix copies escal[] used for scaled fitting columns.
			 * If the maps are not scaled, or we are in GF(2), we do not use escal[].
			**/
	ccfr = (int*)(zpto+cto);  mpto = ccfr+cfr;
	rosc = mpto+cto;  cosc = rosc+ro;
	if (r>=0) {
		k = struct_connorder_col(efr,ccfr);
		for (i=0; i<ro; i++) {
			if (scal==0)  j = cfr;
			else  for (j=0; j<cfr && SIGNM(efr,i,ccfr[j])==0; j++) ;
			rosc[i] = j;
		}
		for (j=0; j<cfr; j++) {
			cosc[j] = -1;
			if (scal!=0)  for (i=0; i<ro; i++)
				if (rosc[i]<j && SIGNM(efr,i,ccfr[j])!=0) {
					cosc[j] = i;  break;
			}
			if (j>0 && k>=0 && scal && cosc[j]==-1)  {PROGERROR("wrong matrix connectivity here");}
		}
		for (j=0; j<cto; j++)  mpto[j] = 0;
		a = strmap_new(efr,etox);
		for (i=0; i<ro; i++)  EMLMPR(a,efr,i) = (ax? EMLMPR(ax,efr,i):i);
		for (j=0; j<cto; j++)  EMLMGC(a,efr,j) = 0;
		for (i=0; i<ROWSM(etox); i++)  EMLMGR(a,efr,i) = 0;
		for (i=0; i<ro; i++) {	/* (uncovered rows are 0, covered ones are 1 by default) */
			EMLMGR(a,efr,EMLMPR(a,efr,i)) = 1;
			pfield_setzeroexp(&EMLMXR(a,efr,EMLMPR(a,efr,i)));
		}
		p = (pfnumexp==0 && !pfield_hasnegative());	/* (identifies GF(2) ) */
		if (scal==0 || p)  escal = NULL;
		else if (cfr<18)  escal = escstack;
		else  escal =  MMALLOC((cfr+2)*sizeof(escal[0]));
		if (escal)  for (j=0; j<cfr; j++)  escal[j] = NULL;
	}
			/**
			 * This is the main cycle that tries to fit the columns of efr to eto.
			 * All possibilities allowed by the given hash hcfr,hcto and by the
			 * zero pattern in zpfr,zpto are tried here in a backtracking "tree".
			 *  cmap[] is the column mapping from a,  mpto[] indicates which
			 *  columns of eto are a;ready covered (at level k),
			 *  escal[] keeps scaled copies of the matrix efr.
			 * Remember that the columns are ordered by the connected order ccfr[]!
			**/
	if (r>=0) {
		sm_callfixrows2++;
		cmap = &EMLMPC(a,efr,0);  cmap[ccfr[0]] = -1;
	}
	DEBUG(CURDLEV+3,"- Starting the main column-fitting cycle.\n");
	k = 0;
	while (r>=0 && k>=0) {
		if (k>=cfr)  {PROGERROREXIT("cannot get here!!! %d",k);}
		if (escal?escal[k]:0) {
			dispose_ematrix(escal[k]);  escal[k] = NULL;
		}
		if (cmap[ccfr[k]]>=0)  mpto[cmap[ccfr[k]]] = 0;
					/* finding next suitable map of col ccfr[k] */
		for (j=cmap[ccfr[k]]+1; j<cto; j++)
			if (mpto[j]==0 && zpfr[ccfr[k]]==zpto[j] \
					&& ((hcfr&&hcto)? hcfr[ccfr[k]]==hcto[j]:1))
				break;
		if (j>=cto || r>SM_MAXLINEMAPS) {	/* no more choices at this level, or finished */
			--k;  continue;
		}
		cmap[ccfr[k]] = j;  mpto[j] = 1;
		DEBUG(CURDLEV+5,"\t\t%*s. col map try (k=%d) cmap[%d]=%d\n",2*k,"",k,ccfr[k],cmap[ccfr[k]]);
		
			/**
			 * Here we deal with one column try at level k (0...k-1 already mapped).
			 * We scale the column according to its nonzero intersection with an
			 * already scaled row (as precomputed in cosc,rosc above) if escal is given.
			 * The parameter scal was already considered when preparing cosc,rosc.
			 * Alternativaly, the column is not scaled and 1 is stored.
			**/
		if (escal)  escal[k] = ematrix_copy(k>0? escal[k-1]:efr);
		esc = (escal? escal[k]:efr);
		if ((escal && cosc[k]>=0)? (SIGNM(eto,cosc[k],cmap[ccfr[k]])!=0) :0) {
			pfield_mul(-1,EXPM(esc,cosc[k],ccfr[k]),SIGNM(esc,cosc[k],ccfr[k]), \
				1,EXPM(eto,cosc[k],cmap[ccfr[k]]),SIGNM(eto,cosc[k],cmap[ccfr[k]]), \
				&xxc,&ggc);
			if (ggc==0)  {PROGERROR("cannot get zero col scale here");}
			if (ggc!=1 || !pfield_iszeroexp(xxc))
				ematrix_multiply_col(esc,ccfr[k], 1,xxc,ggc);
		} else {
			pfield_setzeroexp(&xxc);  ggc = 1;
		}
		EMLMXC(a,efr,cmap[ccfr[k]]) = xxc;  EMLMGC(a,efr,cmap[ccfr[k]]) = ggc;
			/**
			 * We continue the previous task with scaling indicated rows rosc[i]==k
			 * (that have first nonzero appearance in this column).
			 * All entries in the column ccfr[k] are then (after possible row scaling)
			 * compared for exact equality.
			**/
		for (i=0; i<ro; i++) {
			if (escal && rosc[i]==k && SIGNM(eto,i,cmap[ccfr[k]])!=0) {
				pfield_mul(-1,EXPM(esc,i,ccfr[k]),SIGNM(esc,i,ccfr[k]), \
					1,EXPM(eto,i,cmap[ccfr[k]]),SIGNM(eto,i,cmap[ccfr[k]]), \
					&xxr,&ggr);
				if (ggr==0)  {PROGERROR("cannot get zero row scale here");}
				if (ggr!=1 || !pfield_iszeroexp(xxr))
					ematrix_multiply_row(esc,i, 1,xxr,ggr);
				EMLMXR(a,efr,EMLMPR(a,efr,i)) = xxr;
				EMLMGR(a,efr,EMLMPR(a,efr,i)) = ggr;
				p = 1;		/* (row scales are 1 by default, as set above) */
			} else  p = 0;
			if (!ematrix_isequal_entry(esc,i,ccfr[k],eto,i,cmap[ccfr[k]])) {
				if (p)  {PROGERROR("wrong scaling of a row %d at k=%d (%d->%d)",i,k,ccfr[k],cmap[ccfr[k]]);}
				break;
			}
		}
		if (i<ro)  continue;	/* if the entry in some row was not equal (after possible scaling) */
		DEBUG(CURDLEV+4,"\t\t%*s: col map to next (k=%d) cmap[%d]=%d (*%s)\n",2*k,"",
				k,ccfr[k],cmap[ccfr[k]],pfield_pvalue(7,xxc,ggc));
		EMLMAPDEBUGS(CURDLEV+5,a,efr,"\t\t\t\t\t..\t");
		
			/**
			 * We get here if the column map at level k has passed all tests.
			 * We continue deeper if k<cfr-1 with next columns.
			 * When all columns are mapped, we look whether the map is "minimal"
			 * with respect to the col permutations given in afo.
			 * If not, then the this map is factored-out and not recorded.
			 * We record a copy of a in the list al otherwise.
			 * Notice that we have to finish the cycle in gradual steps k--
			 * in order to release the above allocated data (escal[]).
			**/
		if (k<cfr-1) {
			cmap[ccfr[++k]] = -1;	/* passed, go to the next level */
			continue;
		}
		for (x=afo, p=0; (x?*x:0) && p>=0; x++) {
			for (j=0; j<cfr; j++)
				if ((p = cmap[EMLMPC(*x,efr,j)]-cmap[j])!=0)  break;
			if (p<0)  DEBUG(CURDLEV+3,"\tcol map factored-out by the list afo\n");
		}
		if (p<0)  continue;
#ifndef	FASTPROG
		for (i=j=0; j<cto; j++)  i += mpto[j];
		if (i!=cfr)  {PROGERROR("not all columns of efr were mapped, why? %d!=%d",i,cfr); continue;}
		EMLMAPDEBUGS(CURDLEV+3,a,efr,"\t\tcol:\t");
		if (IFRANDDEBUGLESS(111))  strmap_check_comp(a,efr,etox);
		else  strmap_check(a,efr,etox);	/* (errors are reported inside) */
#endif
		if (ch>=4) {
			OUTPUT("Found a #%d %sline map of matrix %p (%dx%d) to matrix %p (%dx%d).\n",
					r+1+ll,scal?"":"unscaled ",efr,ro,cfr,etox,ROWSM(etox),cto);
			EMLMAPOUTPUT(a,efr,NULL);
		}
		if (!al) {
			if (ch<4)  r = SM_MAXLINEMAPS+1;  else  r++;
		} else {
			r++;  *al = alist_append(*al,strmap_copy(a));
			if (r>SM_MAXLINEMAPS)  {PROGERROR("Too many line maps %d>=%d are generated, quitting.\n",r,SM_MAXLINEMAPS);}
		}
	}	/* (end of the main "while (r>=0 && k>=0)" cycle) */
	if (r==0 && k<0)  r = -1;
	
	DEBUG(CURDLEV+2+(r<0),"-> Found %s%d %s-rows %sline maps of %p (%dx%d) to %p (%dx%d).\n",
			r<0?"-NO-":(al?"":">="),r<0?0:(al?r:1),ax?"given":"fixed",scal?"":"unscaled ",efr,ro,cfr,etox,ROWSM(etox),cto);
	if (zpfr && zpfr!=zpstack)  FREE(zpfr);
	if (escal && escal!=escstack)  FREE(escal);
	if (eto!=etox)  dispose_ematrix(eto);
	if (a)  FREE(a);
	return (r<0? -1 : ((al||ch>=4)?r:0));
}




/**
 * This function is called to find injective line maps from the matrix efr to
 * the matrix eto(x) (unscaled mappings if scal==0) on the same number of columns.
 * (We try all possible row and column maps here, the condition of the same number
 * of columns here is for faster computation, there is no pre-fixed col map.)
 * The matrices to map are in efr, etox.
 * Alternatively, etox is a larger matrix from which we select the columns given in the
 * unordered column map ax if given (and this is reflected in the output maps in *al).
 * The arrays hcfr[],hcto[] give column hash codes - if a column j of efr can be mapped to
 * a column jj of etox, then it must be hcfr[j]==hctox[jj].
 * The same applies to row hash codes hrfr[],hrto[].
 * (Code 0 is not tested in the mapping.)
 * 
 * Line maps are consider in the sense described in  str.h  and above.
 * If al is given, then all the maps found here are appended to the list *al.
 * (The list must be initialized prior to calling this function!)
 * If the list afo is given, then its maps are supposed to be a permutation group
 * on the rows and columns of efr, and these permutations are factored out from all maps here.
 * 
 * The return value is -1 for no mapping found, or 0 for existing mappings but al==NULL,
 * or the number of mappings found and stored in the list *al.
 * If ch>=3, then all mappings are found and printed out (need ch>=4 if al==NULL).
**/

int	strmap_samecols_ext(int ch, ematrix *efr, ematrix *etox, int scal, emlinemap **afo,
				long hrfr[], long hrto[], long hcfr[], long hctox[],
				emlinemap *ax, emlinemap ***al) {
	emlinemap	*a=NULL,*aa, **alm, **afoi=NULL, **x;
	ematrix		*eto, *eend;
	long		*hcto=NULL, hcstack[20], *nzfr=NULL,*nzto=NULL, nzstack[30];
	short		*rmap;
	long		**dtfr=NULL,**dtto=NULL, *mpto=NULL;
	int		i,ii,j,k, p,q,r, rfr,rto,co, nend, db,dbf,dbf2;
	extern int	sm_callfixrows, sm_callfixrows2, sm_callsamecols, sm_callsamecols2;	/* (for debug statistics only) */
	
	if ((efr&&etox)? (COLSM(efr)!=COLSM(etox) && !ax):1)
		{PROGERROR("both efr,etox must be given, on the same number of cols"); return -1;}
	rfr = ROWSM(efr);  rto = ROWSM(etox);
	co = COLSM(efr);
	if (afo?(alist_getlength(afo)<=1):0)  afo = NULL;
#ifndef FASTPROG
	if (ax? (strmap_check_size(ax,efr,etox)<0):0)
		{PROGERROR("given map ax must agree with the dimensions of efr,etox"); return -1;}
	if (afo?(afo[0]):0) if (strmap_check_size(afo[0],efr,efr)<0)
		{EMATDEBUGS(0,efr," !\t"); EMLMAPDEBUGS(0,afo[0],efr,"   !\t"); afo = NULL;
		 PROGERROR("given factor-out maps in afo must be automorphisms of efr");}
	DEBUG(CURDLEV+2,"Looking for %s-cols %sline maps of %p (%dx%d) to %p (%dx%d).\n",
			ax?"slct":"same",scal?"":"unscaled ",efr,rfr,co,etox,rto,COLSM(etox));
	if (afo)  DEBUG(CURDLEV+2," (Factoring-out list of %d maps afo.)\n",alist_getlength(afo));
	EMATDEBUGS(CURDLEV+2,efr,"\t\t\tfr<\t");  EMATDEBUGS(CURDLEV+2,etox,"\t\t\tto>\t");
	if (ax)  EMLMAPDEBUGS(CURDLEV+2,ax,efr,"\t\tax-\t");
	if (hrfr && rfr>3)  DEBUG(CURDLEV+2," (Row hash codes hrfr=[%ld,%ld,%ld,%ld] -> hrto=[%ld,%ld,%ld,%ld].\n",hrfr[0],hrfr[1],hrfr[2],hrfr[3],hrto[0],hrto[1],hrto[2],hrto[3]);
#endif
	dbf = sm_callfixrows;  dbf2 = sm_callfixrows2;
	sm_callsamecols++;
	r = 0;
	if (rfr>rto || ((hrfr&&hrto)? !strmag_isinjection(hrfr,rfr,hrto,rto):0))
		r = -1;
			/**
			 * First we look at the destination matrix etox - if no ax is given,
			 * then we rename etox to eto,
			 * otherwise we extract the columns of etox according to ax col map.
			 * Then we possibly update the column hash according to ax.
			 * We look at the numbers of zeros in the rows, whether they
			 * allow an injective mapping from efr to eto.
			 * If not, we are finished here.
			**/
	if (ax && r>=0) {
		eto = ematrix_refer(etox,0,rto,-1,-1);
		for (i=0; i<co; i++)
			ematrix_refadd_col(eto,EMLMPC(ax,efr,i));
	} else  eto = etox;
	if (ax && hctox && r>=0) {
		if (co<18)  hcto = hcstack;
		else  hcto = MMALLOC((co+2)*sizeof(hcto[0]));
		for (i=0; i<co; i++)  hcto[i] = hctox[EMLMPC(ax,efr,i)];
	} else  hcto = hctox;
	if (r>=0 && hcfr && hcto) {
		if (!strmag_isinjection(hcfr,co,hcto,co))  r = -1;
	}
	if (r>=0) {
		if (rfr+rto<27)  nzfr = nzstack;
		else  nzfr = MMALLOC((rfr+rto+2)*sizeof(nzfr[0]));
		nzto = nzfr+rfr;
		strmag_numzero_rows(efr,nzfr);
		strmag_numzero_rows(eto,nzto);
		DEBUG(CURDLEV+3,"Checking numbers of zeros in the matrix rows.\n");
		//******** more can be done here if scal==0 .......
		if (!strmag_isinjection(nzfr,rfr,nzto,rto))  r = -1;
	}
			/**
			 * Then we look at the number of zero 2x2-subdeterminants in all
			 * pairs of rows of efr,eto.
			 * (This number is invariant on scaling and column permutations.)
			 * We store these numbers into allocated 2-dim arrays dtfr,dtto.
			 * The numbers are collected not only from efr,eto themselves,
			 * but also from their images under two endomorphisms of the pfield.
			 * (If a subdeterminant is undefined, then it is taken as nonzero.)
			**/
	if (r>=0) {
		dtfr = malloc_twodim(sizeof(dtfr[0][0]),rfr,rfr);
		dtto = malloc_twodim(sizeof(dtto[0][0]),rto,rto);
		for (i=0; i<rfr; i++)  dtfr[i][i] = 0;
		for (i=0; i<rto; i++)  dtto[i][i] = 0;
		nend = pfield_endomorph_number();  if (nend>2)  nend = 2;
		DEBUG(CURDLEV+3,"Computing twolinedet for pairs of rows (end# %d).\n",nend);
		for (k=0; k<=nend; k++) {
			eend = ematrix_copy(efr);
			if (k>0)  ematrix_pfendomorph(k-1,eend);
			for (i=0; i<rfr; i++)  for (ii=0; ii<i; ii++) {
				p = strmag_twolinedet_numz(eend,0,ii,i);
				dtfr[ii][i] = dtfr[i][ii] = p+ (k==0?0: 1001*dtfr[i][ii]);
			}
			if (eend!=efr)  dispose_ematrix(eend);
			eend = ematrix_copy(eto);
			if (k>0)  ematrix_pfendomorph(k-1,eend);
			for (i=0; i<rto; i++)  for (ii=0; ii<i; ii++) {
				p = strmag_twolinedet_numz(eend,0,ii,i);
				dtto[ii][i] = dtto[i][ii] = p+ (k==0?0: 1001*dtto[i][ii]);
			}
			if (eend!=eto)  dispose_ematrix(eend);
		}
		mpto = MMALLOC((rto+2)*sizeof(mpto[0]));
		for (i=0; i<rto; i++)  mpto[i] = 0;
	}
			/**
			 * This is the main cycle that tries to fit the rows of efr to eto.
			 * All possibilities allowed by the given hash hrfr,hrto, and by the
			 * arrays nzfr,nzto and dtfr,dtto are tried here in a backtracking "tree".
			 *  rmap[] is the row mapping from a,  mpto[] indicates which
			 *  rows of eto are already covered (at level k).
			**/
	a = strmap_new(efr,eto);
	rmap = &EMLMPR(a,efr,0);  rmap[0] = -1;
	if (r>=0)  DEBUG(CURDLEV+3,"== Starting the main row-fitting cycle.\n");
	k = db = 0;
	while (r>=0 && k>=0) {
		if (k>=rfr)  {PROGERROREXIT("cannot get here!!! %d",k);}
		if (rmap[k]>=0)  mpto[rmap[k]] = 0;
					/* finding next suitable map of row k */
		for (i=rmap[k]+1; i<rto; i++)
			if (mpto[i]==0 && nzfr[k]==nzto[i] && ((hrfr&&hrto)? hrfr[k]==hrto[i]:1))
				break;
		if (i>=rto || r>SM_MAXLINEMAPS) {	/* no more choices at this level, or finished */
			--k;  continue;
		}
		rmap[k] = i;  mpto[i] = 1;
		DEBUG(CURDLEV+4,"\t%*s: row map try (k=%d) rmap[%d]=%d\n",2*k,"",k,k,rmap[k]);
		if (k<=rfr/4) {	/* (injection on the second index of two-line determinants) */
			if (!strmag_isinjection(dtfr[k],rfr,dtto[rmap[k]],rto))
				continue;
		}
		for (i=0; i<k; i++) {
			if (dtfr[i][k]!=dtto[rmap[i]][rmap[k]])  break;
		}
		if (i<k)  continue;
		DEBUG(CURDLEV+3,"\t%*s: row map next (k=%d) rmap[%d]=%d\n",2*k,"",k,k,rmap[k]);
		if (k>rfr/5)  db = 1;	/* (for debug statistics only) */
			/**
			 * We proceed to the next level, until all rows of efr are mapped.
			 * Then we check the row map for "minimality" under the row permutations
			 * given in the list afo, and we (first time) record in afoi those
			 * maps of afo that are identical on the rows for next use.
			 * Then we call strmap_fixedrows_() to find column mappings respecting
			 * the selected row mapping (given by rmap[] through a).
			 * If ax was given, then we must also re-map the columns according
			 * to the columns of etox selected by ax.
			**/
		if (k<rfr-1) {
			rmap[++k] = -1;		/* passed, go to the next level */
			continue;
		}
		for (x=afo, p=0, q=!afoi; (x?*x:0) && p>=0; x++) {
			for (i=0; i<rfr; i++)
				if ((p = rmap[EMLMPR(*x,efr,i)]-rmap[i])!=0)  break;
			if (p<0)  DEBUG(CURDLEV+2,"\trow map factored-out by the list afo\n");
			if (q && p==0)  afoi = alist_append(afoi,*x);
		}
		if (p<0)  continue;
		
		alm = NULL;			/* rows are all mapped, now try the columns: */
		p = strmap_fixedrows_ext(ch,efr,eto,scal,afoi,hcfr,hcto,a,(al?&alm:NULL),r);
		if (p<0)  continue;
		if (p>0)  r += p;  else  r++;
		DEBUG(CURDLEV+2,"Found %d new %s-cols %sline maps for the row mapping.\n",
				p>0?p:1,ax?"slct":"same",scal?"":"unscaled ");
		if (ch>=3 && p>0)  OUTPUT("Found %d new %s-cols %sline maps for the selected row map of %p [%s] (%dx%d) to %p [%s] (%dx%d).\n",
				p,ax?"select":"same",scal?"":"unscaled ",efr,EMNAME(efr)?EMNAME(efr):"",ROWSM(efr),COLSM(efr),etox,EMNAME(etox)?EMNAME(etox):"",ROWSM(etox),COLSM(etox));
		if (p<SM_MAXLINEMAPS && r>SM_MAXLINEMAPS) {
			PROGERROR("Too many line maps %d>=%d are generated, quitting.\n",r,SM_MAXLINEMAPS);
		}
		if (!al) {
			if (ch<3)  r = SM_MAXLINEMAPS+1;	/* (already done with one map) */
		} else if (!ax) {
			*al = alist_applist(*al,alm);
			alm = NULL;
		} else for (x=alm; x?*x:0; x++) {	/* (all maps need to be adapted to ax) */
			aa = strmap_new(efr,etox);
			for (i=0; i<rfr; i++)  EMLMPR(aa,efr,i) = EMLMPR(*x,efr,i);
			for (j=0; j<COLSM(etox); j++)  EMLMGC(aa,efr,j) = 0;
			for (i=0; i<rto; i++) {
				EMLMGR(aa,efr,i) = EMLMGR(*x,efr,i);
				EMLMXR(aa,efr,i) = EMLMXR(*x,efr,i);
			}
			for (j=0; j<co; j++) {
				ii = EMLMPC(aa,efr,j) = EMLMPC(ax,efr,EMLMPC(*x,efr,j));
				EMLMGC(aa,efr,ii) = EMLMGC(*x,efr,EMLMPC(*x,efr,j));
				EMLMXC(aa,efr,ii) = EMLMXC(*x,efr,EMLMPC(*x,efr,j));
			}
			*al = alist_append(*al,aa);
#ifndef	FASTPROG
			EMLMAPDEBUGS(CURDLEV+3,aa,efr,"\t\tcolx:\t");
			if (IFRANDDEBUGLESS(111))  strmap_check_comp(aa,efr,etox);
			else  strmap_check(aa,efr,etox);
#endif
			if (ch>=3)  EMLMAPOUTPUT(aa,efr,NULL);
		}
		if (alm)  dispose_alist(alm);
	}	/* (end of the main "while (r>=0 && k>=0)" cycle) */
	if (r==0 && k<0)  r = -1;
	if (ch>=3)  OUTPUT("Found total %s%d %s-cols %sline maps from matrix %p [%s] (%dx%d) to %p [%s] (%dx%d).\n",
			r<0?"-N":"",r<0?0:r,ax?"select":"same",scal?"":"unscaled ",efr,EMNAME(efr)?EMNAME(efr):"",ROWSM(efr),COLSM(efr),etox,EMNAME(etox)?EMNAME(etox):"",ROWSM(etox),COLSM(etox));
	sm_callsamecols2 += db;
#ifndef	FASTPROG
	DEBUG(CURDLEV+1+(r<0)+(ch<0),"=>> Found %s%d %s-cols %sline maps of %p (%dx%d) to %p (%dx%d).\n",
			r<0?"-N":(al?"":">="),r<0?0:(al?r:1),ax?"slct":"same",scal?"":"unscaled ",efr,rfr,co,etox,rto,COLSM(etox));
	DEBUG(CURDLEV+2+(ch<0),"  (Called fixedrows %d times, %d fully.).\n",dbf-sm_callfixrows,dbf2-sm_callfixrows2);
	if (ch>=0 && IFRANDDEBUGLESS(111) && !afo) {
				/* recursive testing with exchanged efr lines, respecting h?fr[] if given */
		eend = ematrix_copy(efr);  EMSETNAME(eend,"rec-test");
		for (ii=0; ii<4; ii++) {
			i = RANDOM()%rfr;  j = RANDOM()%rfr;
			if (hrfr) if (hrfr[i]!=hrfr[j])  continue;
			ematrix_swaprows(eend,i,j);
		}
		for (ii=0; ii<4; ii++) {
			i = RANDOM()%co;  j = RANDOM()%co;
			if (hcfr) if (hcfr[i]!=hcfr[j])  continue;
			ematrix_swapcols(eend,i,j);
		}
		q = strmap_samecols_ext(-1,eend,etox,scal,NULL,hrfr,hrto,hcfr,hctox,ax,NULL);
		if ((q<0)!=(r<0))  {PROGERROR("non-symmetric result for modified (line-exchange efr) matrices");}
		dispose_ematrix(eend);
	}
#endif
	if (hcto && hcto!=hcstack && hcto!=hctox)  FREE(hcto);
	if (nzfr && nzfr!=nzstack)  FREE(nzfr);
	if (eto!=etox)  dispose_ematrix(eto);
	if (a)  FREE(a);
	if (afoi)  alist_free(afoi);	/* (the list elements belong to afo, so they are not disposed!) */
	if (dtfr)  FREE(dtfr);
	if (dtto)  FREE(dtto);
	if (mpto)  FREE(mpto);
	return (r<0? -1 : ((al||ch>=3)?r:0));
}













/******************	General matrix line maps	*******************/
/**************************************************************************/


#undef CURDLEV
#define CURDLEV		6		/* (need <=DEBUGLEV+2 to perform paranoic tests) */


/**
 * This function is called to find all injective line maps of the matrix efr to the matrix
 * eto (unscaled mappings if scal==0).
 * The arrays hcfr[],hcto[] give column hash codes - if a column j of efr can be mapped to
 * a column jj of etox, then must be hcfr[j]==hctox[jj].
 * The same applies to row hash codes hrfr[],hrto[].
 * (Code 0 is not tested in the mapping.)
 * All (unordered) subsets of columns are tried here, so consider this for a transposition.
 * 
 * Line maps are consider in the sense described in  str.h  and above.
 * If al is given, then all the maps found here are appended to the list *al.
 * (The list must be initialized prior to calling this function!)
 * If the list afo is given, then its maps are supposed to be a permutation group
 * on the rows and columns of efr, and these permutations are factored out from all maps here.
 * 
 * The return value is -1 for no mapping found, or 0 for existing mappings but al==NULL,
 * or the number of mappings found and stored in the list *al.
 * If ch>=3(2), then all mappings are found and printed out.
**/

int	strmap_allbysubs_ext(int ch, ematrix *efr, ematrix *eto, int scal, emlinemap **afo,
			long hrfr[], long hrto[], long hcfr[], long hcto[], emlinemap ***al) {
	emlinemap	*a=NULL;
	int		k, p,r, cfr,cto, dbs,dbf,db;
	short		*subs;
	extern int	sm_callfixrows2, sm_callsamecols2;	/* (for debug statistics only) */
	
	if (!efr || !eto)  {PROGERROR("both efr,eto must be given here"); return -1;}
	cfr = COLSM(efr);  cto = COLSM(eto);
	a = strmap_new(efr,eto);
	subs = &EMLMPC(a,efr,0);
#ifndef FASTPROG
	DEBUG(CURDLEV+1+(ch<0),"Looking for all %sline maps of %p (%dx%d) to %p (%dx%d).\n",
			scal?"":"unscaled ",efr,ROWSM(efr),cfr,eto,ROWSM(eto),cto);
	if (afo && ch>=0)  DEBUG(CURDLEV+1," (Factoring-out list of %d maps afo.)\n",alist_getlength(afo));
	if (hrfr && ROWSM(efr)>2 && ch>=0)  DEBUG(CURDLEV+1," (Using row hash codes in hrfr=[%ld,%ld,%ld,..] -> hrto=[%ld,%ld,%ld,..].)\n",
			hrfr[0],hrfr[1],hrfr[2],hrto[0],hrto[1],hrto[2]);
	if (hcfr && COLSM(efr)>2 && ch>=0)  DEBUG(CURDLEV+1," (Using column hash codes in hcfr=[%ld,%ld,%ld,..] -> hcto=[%ld,%ld,%ld,..].)\n",
			hcfr[0],hcfr[1],hcfr[2],hcto[0],hcto[1],hcto[2]);
	EMATDEBUGS(CURDLEV+2,efr,"\t\tfr<<\t");  EMATDEBUGS(CURDLEV+2,eto,"\t\tto>>\t");
#endif
	dbs = sm_callsamecols2;  dbf = sm_callfixrows2;  db = 0;
			/**
			 * Here we simply try all subsets of efr columns in eto.
			 * The column hash codes are actually checked inside strmap_samecols_().
			**/
	k = r = 0;  subs[0] = -1;
	if (COLSM(efr)>COLSM(eto) || ROWSM(efr)>ROWSM(eto))
		r = -1;
	while (r>=0 && k>=0) {
		if (++subs[k]>=cto-cfr+k+1) {
			k--;  continue;
		}
		DEBUG(CURDLEV+4,"\t%*s: col subset next (k=%d) subs[%d]=%d\n",2*k,"",k,k,subs[k]);
		if (k<cfr-1) {
			++k;  subs[k] = subs[k-1];  continue;
		}
		db++;
		p = strmap_samecols_ext(ch,efr,eto,scal,afo,hrfr,hrto,hcfr,hcto,a,al);
		if (p<0)  continue;
		if (p>0)  r += p;  else  r++;
		if (!al && ch<2)  break;
		if (ch==2 && p>0)  OUTPUT("Total %d new %sline maps found from %p [%s] (%dx%d) to selected %d-column subset in %p [%s] (%dx%d).\n",
				p,scal?"":"unscaled ",efr,EMNAME(efr)?EMNAME(efr):"",ROWSM(efr),COLSM(efr),cfr,eto,EMNAME(eto)?EMNAME(eto):"",ROWSM(eto),COLSM(eto));
	}
	if (r==0 && k<0)  r = -1;
	if (ch>=2)  OUTPUT("Altogether %s%d %sline maps found from matrix %p [%s] (%dx%d) to matrix %p [%s] (%dx%d).\n",
			r<0?"-N":"",r<0?0:r,scal?"":"unscaled ",efr,EMNAME(efr)?EMNAME(efr):"",ROWSM(efr),COLSM(efr),eto,EMNAME(eto)?EMNAME(eto):"",ROWSM(eto),COLSM(eto));
	
#ifndef	FASTPROG
	DEBUG(CURDLEV+(r<0)+(ch<0),"==>> Found total %s%d %sline maps of %p (%dx%d) to %p (%dx%d).\n",
			r<0?"-NO-":(al?"":">="),r<0?0:(al?r:1),scal?"":"unscaled ",efr,ROWSM(efr),cfr,eto,ROWSM(eto),cto);
	DEBUG(CURDLEV+1+(ch<0)," (Tested %d column subsets, through %d (deep) same-col tries, and %d fix-row tries.).\n",db,sm_callsamecols2-dbs,sm_callfixrows2-dbf);
	if (ch>=0 && IFRANDDEBUGLESS(111)) {
			/* recursive testing with transposed matrices (afo is OK when transposed!) */
		ematrix	*efx,*etx;
		efx = ematrix_copydual(efr);  etx = ematrix_copydual(eto);
		p = strmap_allbysubs_ext(-1,efx,etx,scal,afo,hcfr,hcto,hrfr,hrto,NULL);
		if ((p<0)!=(r<0))  {PROGERROR("non-symmetric result for transposed matrices");}
		dispose_ematrix(efx);  dispose_ematrix(etx);
	}
#endif
	if (a)  FREE(a);
	return (r<0? -1 : ((al||ch>=2)?r:0));
}


















#if 0

************** computing by the old dispminor algorithm...........

		/* macros for access to quadruple-indexed hash-values for mapping pairs of rows - see ..rowhash() */
#define	STR_MAXABSVALIX	2	/* (the max number of multiplicative endomorphisms used) */
#define IHSETVARACCESS(em,eb)	long  *pti; \
				int   rbi=ROWSM(eb), rmi=ROWSM(em),cmi=COLSM(em)
#define IHGETPOINTER		pti
#define IHGETSIZE		((rmi)*(rmi)*(rbi)*(rbi)*(cmi))
#define IHGET(rm1,rm2,rb1,rb2)	((pti)+((((rm1)*rmi+(rm2))*rbi+(rb1))*rbi+(rb2))*(cmi))


int	strmap_allbyorder_ext() {
	
	***********
}



#undef CURDLEV
#define	CURDLEV		7		/*7 (need <=DEBUGLEV+2 to perform paranoic tests) */


/**
 * This function looks in the matrix e for displayed minors strongly-equivalent to emin.
 * In other words, it looks in all ordered subsets of rows and columns corresponding
 * to the dimensions of emin, whether the respective submatrix is a scale of emin.
 * Possible automorphic row- and column-permutations of emin are factored out,
 * so that one displayed minor comes out only once.
 * The function works only for connected minors(!), otherwise an error is issued.
 * The return list of matrices consists of refering (to e) matrices showing the emin-minor
 * exactly as it corresponds to the entries of emin, up to scaling of course.
 * The list is NULL terminated.
 * (Remember to free the list and the matrices later.)
 * 
 * Call this function only for a standalone (not refering) matrix, since otherwise the
 * returned list is not well-defined.
 * If emin is empty or bigger than e, then nothing is returned in a list.
 * If ch>0, some more information are printed, and extra paranoic tests are performed if ch>=0.
 * If ch==-3, then we only test for a minor, but do not generate any list (retun 0 or 1).
 * If ch==-2, then no minors are collected in the list, but all of them are still tested.
 * 
 * When looking for a displayed minor, we use quite complex hash-values for maps of rows
 * and of pairs of rows (minor -> matrix), which are described in estruct_dispminor_rowhash()
 * and estruct_dispminor_fitldet().
 * These values consider zero entries and zero 2x2 subdeterminants in the rows and the pairs.
 * The subdeterminants are computed also with respect to the "absolute values".
 * (Any mapping satisfying |a*b|=|a|*|b| works here - it preserves zero 2x2 determinants.)
 * In result, the hash values tell us what columns may be used in the displayed minor respecting
 * the selected rows (columns are considered unordered for now).
 * Once having a possible (ordered) row map from the minor to the matrix, we call the function
 * estruct_fitcols() to properly map the columns of the minor.
 * 
 * How may one deal with isomorphic bases when looking for displayed minors ???
 * The best would be to find a way how to "guess" correct basis for the searched minor...
 * 
**/

#define	MXMN	(NUM_LONG_BITS/2-2)		/* (must fit in long bits!) */
int	callhash, rejhash, rejrowf, rejcolf;	/* (for debugging statistics, here and in struct.c) */

FILE	*struct_printdispm = NULL;		/* request to print the displayed minors to the file */
char	*struct_printdispp = NULL;		/*  prefix for the print request */
int	struct_printdispa = 0;			/*  request to print all the displayed minors */
int	struct_printdispix = 0;			/*  minor-list index for the print request */


ematrix**	struct_dispminor_ext(int ch, ematrix *e, ematrix *emin) {
	ematrix		**mout;
	
	if (ch<-3)  {PROGERROR("Do not call with ch==%d <-3",ch);}
	if (ROWSM(e)>=MXMN || COLSM(e)>=MXMN)	/* (not practical for larger matrices anyway) */
		{PROGERROREXIT("not designed to deal with %d,%d > %d rows or columns",ROWSM(e),COLSM(e),MXMN);}
	if (ISREFMAT(e))  {PROGERROR("Do not call for a refering matrix e=%p.",e);}
	DEBUG(CURDLEV-1,"looking for %s displayed %p-minors (%dx%d) in %p (%dx%d)...\n",
			ch>-3?"all":"exi",emin,ROWSM(emin),COLSM(emin),e,ROWSM(e),COLSM(e));
	if (ch>0) { EMATDEBUG(CURDLEV+2,e,"\t");  EMATDEBUG(CURDLEV+2,emin,"\t>\t"); }
	mout = NULL;
	
	if (ROWSM(emin)<=ROWSM(e) && COLSM(emin)<=COLSM(e) && COLSM(emin)>0 && ROWSM(emin)>0) {
		mout = struct_dispminor_core(ch,e,emin,mout);	/* here is the real computation */
	}
	if (ch>-2) {
		DEBUG(CURDLEV-1-(mout!=NULL)," ... done, found %d displayed %p-minors in the matrix %p.\n",alist_getlength(mout),emin,e);
		if (mout && mout==(void*)1)  {PROGERROR("Should have collected the minors! %p",mout);}
	} else {
		DEBUG(CURDLEV-1-(mout!=NULL)," ... done, asked only for an existence of a minor, found %s.\n",mout!=NULL?"YES":"NO");
		if (mout && mout!=(void*)1)  {PROGERROR("Should have only indicated the existence of a minor! %p",mout);}
	}
	return mout;
}

			/* here is the real computation for struct_dispminor() ... */
ematrix**	struct_dispminor_core(int ch, ematrix *e, ematrix *emin, ematrix **mout) {
	int		i,j,k, br, hrr=0, pr[MXMN];
	long		hse[MXMN],hsmin[MXMN], hce[MXMN],hcmin[MXMN];
	long		*subca[MXMN+1],**subc, dsra[MXMN+1],*dsr;
	emxaut		**aut, **auti, **a;
	IHSETVARACCESS(emin,e);	/* supplementary variables for two-row mapping hash, see str.h */
	
	aut = auti = NULL;	/* (obtained only later on demand...) */
	callhash++;
		/**
		 * We first compute hash-values for pairs of rows of emin against e:
		 * Each entry IHGET(rm1,rm2,rb1,rb2) is an array assigned to the (possible)
		 * mapping of rows rm1,rm2 of emin (ordered!) to rows rb1,rb2 of e.
		 * The j-th entry of this array keeps a bitfield showing columns of e
		 * that may become j-th images of some columns of emin (no respect to their
		 * order in emin).
		 * Additionally, single-row hash codes are derived from the above values
		 * and stored in hsmin[],hse[] -- hsmin[i]==hse[j] if we can map row i->j.
		**/
	IHGETPOINTER = MMALLOC(IHGETSIZE*sizeof(IHGETPOINTER[0]));
	br = struct_dispminor_rowhash(ch,emin,e,hsmin,hse,IHGETPOINTER);
	if (br<0) {				/* if there is no row mapping possible */
		FREE(IHGETPOINTER);  if (br<0)  rejhash++;
		DEBUG(CURDLEV,"- row-pair hash test does not allow a minor occurence\n");
		return mout;
	}
	subca[0] = MMALLOC((ROWSM(emin)+1)*COLSM(emin)*sizeof(subca[0][0]));
	for (i=0; i<ROWSM(emin)+1; i++)  subca[i] = subca[0]+i*COLSM(emin);
	subc = subca+1;
	for (j=0; j<COLSM(emin); j++)  subc[-1][j] = (1l<<COLSM(e))-1;
	dsr = dsra+1;  dsr[-1] = 0;
	
		/**
		 * We now try all ordered mappings of rows of emin to rows of e that respect
		 * the above introduced hash-values for rows and their pairs.
		 * We use a standard backtracking technique to cycle through all mappings.
		 * We also check against the permutations of emin rows, see below.
		 * (The array subc[k][] keeps "accumulated" row-pair hash codes by the columns.)
		**/
	pr[0] = -1;  k = 0;
	while (k>=0) {
		if (k>=ROWSM(emin))  {PROGERROR("cannot get here!!! %d",k);}
		if (++pr[k]>=ROWSM(e)) {
			--k;  continue;			/* no more choices on this level */
		}
		if ((dsr[k-1]&(1l<<pr[k]))!=0)  continue;	/* (already mapped to pr[k]) */
		dsr[k] = dsr[k-1]|(1l<<pr[k]);
		DEBUG(CURDLEV+3,"%*smap considered pr[%d]=%d\n",2*k,"",k,pr[k]);
		if (hsmin[k]!=hse[pr[k]])  continue;	/* failed one-row hash */
		
		for (j=0; j<COLSM(emin); j++)  subc[k][j] = subc[k-1][j];
		for (j=0; j<COLSM(emin); j++) {
			for (i=0; i<k; i++)		/* check two row hash for all added pairs */
				subc[k][j] &= IHGET(i,k,pr[i],pr[k])[j];
			if (subc[k][j]==0)  break;
		}
		if (j<COLSM(emin))  continue;		/* when failed two-row hash */
#ifndef FASTPROG
		DEBUG(CURDLEV+2,"%*s  map found pr[%d]=%d%s\n",2*k,"",k,pr[k],k<ROWSM(emin)-1?"":"   \t--");
		if (CURDLEV+4<=DEBUGLEV) {
			SDEBUG(CURDLEV+4,"\t\t\t\t\t\t subc[] =");
			for (j=0; j<COLSM(emin); j++)  SDEBUG(CURDLEV+4," %lX,",subc[k][j]);
			SDEBUG(CURDLEV+4,"\n"); }
#endif
		if (k<ROWSM(emin)-1) {
			pr[++k] = -1;  continue;	/* continuing the search up */
		}
			/**
			 * Since we output the dispminors with their automorphisms factored out,
			 * we need to try one row mapping up to row permutations only.
			 * This is achieved by allowing only the unique "minimal" row selection
			 * from all permutations of the considered row map.
			 * (We later also factor out all emin-automorphisms fixing rows from
			 * selected columns.)
			**/
		if (aut==NULL) {		/* getting automorphisms of the minor emin */
			aut = struct_automorphisms_ch(ch<-1?-1:DECRCH(ch),emin,1);
			auti = NULL;		/* those automorphisms identical on the rows */
			for (a=aut; *a; a++) {
				for (j=br=0; j<ROWSM(emin); j++)  br |= ((*a)->pr[j]!=j);
				if (br==0)  auti = alist_append(auti,*a);
			}
			if (!auti)  {PROGERROR("Where is the identical automorphisms???");}
			DEBUG(CURDLEV+1,"  - using %d total and %d row-preserving automorphisms\n",alist_getlength(aut),alist_getlength(auti));
		}
		br = 0;				/* check whether an automorphism gives a "smaller" mapping */
		for (a=aut; *a && !br; a++)
			for (i=0; i<ROWSM(emin); i++) {
				if (pr[(*a)->pr[i]]>pr[i])  break;
				if (pr[(*a)->pr[i]]<pr[i]) { br = 1;  break; }
			}
		if (br)  DEBUG(CURDLEV+2,"X map not minimal with respect to row automorphisms\n");
		if (br)  continue;
		
			/**
			 * If all these tests pass well, then we call estruct_dispminor_rowfix()
			 * to look for the emin minor (by pr[]) in columns of e by brute force.
			 * The column hash-values hce,hcmin are computed there, not here for now.
			 * The row-fixing column permutations auti[] are used there.
			**/
#ifndef FASTPROG
		DEBUG(CURDLEV+1,"passed row selection, pr[] =");
		for (i=0; i<ROWSM(emin); i++)  SDEBUG(CURDLEV+1," %d,",pr[i]);
		SDEBUG(CURDLEV+1,"\n");
#endif
		/* the column hash values hcmin[],hce[] are not yet used here - set in rowfix() */
		mout = struct_dispminor_rowfix(ch,emin,e,hcmin,hce,pr,subc[k],auti,mout);
		hrr = 1;	/* (for debugging statistics only) */
		if (ch<=-3 && mout)  break;	/* (no output list is requested, just existence) */
	}
	
	if (!hrr)  rejrowf++;  else if (!mout)  rejcolf++;
	FREE(IHGETPOINTER);			/* frees the row hash memory allocated above */
	FREE(subca[0]);
	if (auti)  alist_free(auti);		/* frees all the automorphisms */
	if (aut)  struct_automorphisms_recycle(emin,aut);
	return mout;		/* (returns to struct_dispminor_ext()...) */
}



/**
 * Here we compute hash values that could help to decide whether an (ordered) pair of rows
 * of the matrix em can be mapped to a pair of rows of eb within an occurence of a displayed
 * em-minor in eb; the values are stored in hspt (which is indexed by the macro IHGET()).
 * The entries of output hspt are indexed by (ordered) pairs of rows of em times (ordered) pairs of
 * rows of eb, and each entry format is described at estruct_dispminor_fitldet() below.
 * (Briefly speaking, the entries show to which columns of eb on these rows the columns of em
 * may be mapped - a bitfield for each j-th occurence of em in the columns of eb.)
 * The return value is normally 0, but it is -1 if no mapping of rows is possible.
 * Be careful that the matrices obtained here by applying "endomorphisms" may not be
 * properly represented over the pfield!
 * For use only in estruct_dispminor_ext().
**/

static int	maxabsvalix = 0;	/* the number of "abs values" (endomorphisms) we use here */

int	struct_dispminor_rowhash(int ch, ematrix *em, ematrix *eb, long hsm[], long hsb[], void *hspt) {
	
	ematrix		*ema[1+STR_MAXABSVALIX], *eba[1+STR_MAXABSVALIX];
	int		i,ii,j,jj, k,l, id,jd, ad,ax,ld, ret=0;
	int		rmat[MXMN][MXMN], rmin[MXMN][MXMN], ro[MXMN];
	int		*zprm[MXMN][MXMN][1+STR_MAXABSVALIX], *zpb[1+STR_MAXABSVALIX];
	long		zm[MXMN],zb[MXMN];
	emxaut		**aut,**a;	/* (just for debug testing...) */
	IHSETVARACCESS(em,eb);	/* supplementary variables for accessing two-row mapping hash */
	IHGETPOINTER = hspt;
	
	maxabsvalix = pfield_endomorph_number();
	if (maxabsvalix>STR_MAXABSVALIX)  maxabsvalix = STR_MAXABSVALIX;
	DEBUG(CURDLEV+1,"computing row-hash for %p and minor %p (using %d absv) ...\n",eb,em,maxabsvalix);
		/**
		 * We first identify zero entries and zero 2x2 subdeterminants of rows and
		 * pairs of rows in em,eb and in the endomorphism-valued copies ema,eba.
		 * (There may be more than one way to compute "endomorphism" - see maxabsvalix.)
		 * Be careful that the matrices obtained here by applying "endomorphisms"
		 * may not be properly represented!
		 * However, it is still correct to look at those invalid determinants as nonzeros.
		 * 
		 * We then use these zero-patterns to compute the IHGET (=*hspt) hash codes
		 * for pairs of rows of em times eb (see estruct_dispminor_fitldet).
		 * We also determine rmat[][] values showing which pairs of rows of eb
		 * may be images of pairs of rows of em.
		 * After that, we turn the hash-pairs into simplified one-line hash codes.
		**/
	ema[0] = em;  eba[0] = eb;
	for (l=0; l<maxabsvalix; l++) {
		ema[1+l] = ematrix_copy(em);  eba[1+l] = ematrix_copy(eb);
		for (i=0; i<ROWSM(em); i++)  for (j=0; j<COLSM(em); j++)
			pfield_endomorph(l,&EXPM(ema[1+l],i,j),&SIGNM(ema[1+l],i,j));
		for (i=0; i<ROWSM(eb); i++)  for (j=0; j<COLSM(eb); j++)
			pfield_endomorph(l,&EXPM(eba[1+l],i,j),&SIGNM(eba[1+l],i,j));
	}
			/* first compute zero-pattern for each row of em and of eb to bitfields */
	for (i=0; i<ROWSM(em); i++)
		for (j=0, zm[i]=0; j<COLSM(em); j++)
			if (SIGNM(em,i,j)!=0)  zm[i] |= (1l<<j);
	for (i=0; i<ROWSM(eb); i++)
		for (j=0, zb[i]=0; j<COLSM(eb); j++)
			if (SIGNM(eb,i,j)!=0)  zb[i] |= (1l<<j);
	
			/* then compute the pattern of zero 2x2-determinants for pairs of rows */
	zpb[0] = MMALLOC((1+maxabsvalix)*COLSM(eb)*sizeof(zpb[0][0]));
	for (l=0; l<1+maxabsvalix; l++)
		zpb[l] = zpb[0]+l*COLSM(eb);
	zprm[0][0][0] = MMALLOC(ROWSM(em)*ROWSM(em)*(1+maxabsvalix)*COLSM(em)*sizeof(zprm[0][0][0][0]));
	for (i=0; i<ROWSM(em); i++)  for (ii=i+1; ii<ROWSM(em); ii++) {
		for (l=0; l<1+maxabsvalix; l++)
			zprm[i][ii][l] = zprm[ii][i][l] =	/* (always symmetric on i,ii) */
				zprm[0][0][0]+((i*ROWSM(em)+ii)*(1+maxabsvalix)+l)*COLSM(em);
	}
	DEBUG(CURDLEV+4,"computing linedet for minor:\n");
	for (i=0; i<ROWSM(em); i++)  for (ii=0; ii<i; ii++) {
		for (l=0; l<1+maxabsvalix; l++)
			struct_dispminor_linedet(ema[l],0,i,ii,zprm[i][ii][l]);
	}
			/* finally, record what choices in eb agree by the zero patterns with em */
	DEBUG(CURDLEV+4,"computing linedet for major and pair-maps (min->maj):\n");
	for (i=0; i<ROWSM(em); i++)  for (ii=0; ii<ROWSM(em); ii++)
		rmin[i][ii] = 0;
	for (j=0; j<ROWSM(eb); j++)  for (jj=0; jj<j; jj++) {
		rmat[jj][j] = rmat[j][jj] = 0;
		for (l=0; l<1+maxabsvalix; l++)
			struct_dispminor_linedet(eba[l],0,j,jj,zpb[l]);
		for (i=0; i<ROWSM(em); i++)  for (ii=0; ii<ROWSM(em); ii++)  if (i!=ii) {
					/* (this cycle must use all pairs i<ii and i>ii) */
			DEBUG(CURDLEV+4,"@ at maj %d,%d min %d,%d(+rev) :\n",j,jj,i,ii);
			struct_dispminor_fitldet(COLSM(em),zm[i],zm[ii],zprm[i][ii],
						COLSM(eb),zb[j],zb[jj],zpb, IHGET(i,ii,j,jj));
			for (k=0; k<COLSM(em); k++)	/* (mapping is the same for both pairs reversed) */
				IHGET(ii,i,jj,j)[k] = IHGET(i,ii,j,jj)[k];
					/* marks when the pair j,jj can be mapped to: */
			if (IHGET(i,ii,j,jj)[0]!=0)  rmat[jj][j] = rmat[j][jj] = 1;
					/* marks when the pair i,ii can be mapped from: */
			if (IHGET(i,ii,j,jj)[0]!=0)  rmin[ii][i] = rmin[i][ii] = 1;
		}
	}
#ifndef FASTPROG
	DEBUG(CURDLEV+1,"... done computing row-hash for matrix %p and minor %p\n",eb,em);
	if (ch>0 && CURDLEV+2<=DEBUGLEV) {
		SDEBUG(CURDLEV,"\t\t\ttwo-row hash of em over all pairs mapped to:");
		for (j=0; j<ROWSM(em); j++)  for (jj=0; jj<ROWSM(em); jj++)
			SDEBUG(CURDLEV,"%s%c",jj==0?"\n\t\t\t\t":"  ",j==jj?' ':(rmin[j][jj]?'*':'.'));
		SDEBUG(CURDLEV,"\n\t\t\tend table\n");
		SDEBUG(CURDLEV,"\t\t\ttwo-row hash of eb over all pairs of the minor:");
		for (j=0; j<ROWSM(eb); j++)  for (jj=0; jj<ROWSM(eb); jj++)
			SDEBUG(CURDLEV,"%s%c",jj==0?"\n\t\t\t\t":"  ",j==jj?' ':(rmat[j][jj]?'*':'.'));
		SDEBUG(CURDLEV,"\n\t\t\tend table\n");
	}
	if (ch>=0 && IFRANDDEBUG(1111)) {
		aut = struct_automorphisms(em);		/* (aut-invariant test of IHGET) */
		for (a=aut,ax=0; *a; a++)  ax++;
		if (CURDLEV>DEBUGLEV) { id = 1;  jd = 1+(random()%ROWSM(eb))/3; }
		else { id = jd = 1; }
		for (j=0; j<ROWSM(eb); j+=jd)  for (jj=j+1; jj<ROWSM(eb); jj++)
		  for (i=0; i<ROWSM(em); i+=id)  for (ii=i+1; ii<ROWSM(em); ii++) {
			if (CURDLEV>DEBUGLEV) { ad = 1+random()%(ax+1);  ld = 1+random()%(COLSM(em)+1); }
			else { ad = ld = 1; }
			for (a=aut+1; a-aut<ax; a+=ad)  for (l=0; l<COLSM(em); l+=ld)
				if (IHGET((*a)->pr[i],(*a)->pr[ii],j,jj)[l]!=IHGET(i,ii,j,jj)[l])
					{PROGERROR("two-row hash is not automorphism-invariant:  j=%d, jj=%d,  i=%d->%d, ii=%d->%d,  l=%d;  IH %lX!=%lX",
						j,jj,i,(*a)->pr[i],ii,(*a)->pr[ii],l,IHGET(i,ii,j,jj)[l],IHGET((*a)->pr[i],(*a)->pr[ii],j,jj)[l]);}
		}
		struct_automorphisms_recycle(em,aut);
	}
#endif
		/**
		 * Here we look at all pairs of rows of em that may be mapped to some pairs
		 * of rows of eb - these form a graph described by rmin[][] wich must be complete.
		 * Then we look at all pairs of rows of eb that may be images of some pairs
		 * of rows of em - these form a graph described by rmat[][].
		 * An occurence of em-minor in eb must induce a clique in rmat, so we look
		 * whether rmat contains cliques of size >=ROWSM(em) by brute force.
		**/
	ret = 0;  l = 1;
	for (i=ii=0; i<ROWSM(em) && l; i++)  for (ii=0; ii<i && l; ii++)
		if (rmin[i][ii]==0)  l = 0;
	if (l==0) {
		DEBUG(CURDLEV+1,"no minor possible since the \"from\" two-row hash is not complete (%d,%d)\n",i-1,ii-1);
		ret = -1;
	}
	k = l = 0;  ro[0] = -1;
	while (k>=0 && ret>=0) {	/* (checking all ROWSM(em)-tuples of rows of eb for a clique) */
		if (++ro[k]>=ROWSM(eb)-ROWSM(em)+1+k) { --k;  continue; }
		for (j=0; j<k; j++)  if (rmat[ro[k]][ro[j]]==0)  break;
		if (j<k)  continue;
		if (k<ROWSM(em)-1) { ++k;  ro[k] = ro[k-1]; }
		else { l = 1;  break; }
	}
	if (l==0 && ret>=0) {
		DEBUG(CURDLEV+1,"no minor possible since no large clique in the \"to\" two-row hash\n");
		ret = -1;
	}
			/* we try to compute one-row hash ............ */
	for (i=0; i<ROWSM(em); i++)  hsm[i] = 0;
	for (j=0; j<ROWSM(eb); j++)  hsb[j] = 0;
	//********** what can be computed here ? - so that hsm[i]==hsb[j] if i->j in some mapping
	//	- the hash should be computed into a matrix hsm X hsb, not just arrays...
	
	ld = ax = ad = jd = id = ch = 0;
	FREE(zprm[0][0][0]);	/* disposing local computation data */
	FREE(zpb[0]);
	for (l=0; l<maxabsvalix; l++) {
		dispose_ematrix(eba[1+l]);  dispose_ematrix(ema[1+l]);
	}			/* (do not dipose eba[0]==eb, ema[0]==em !!!) */
	return ret;
}


/**
 * Here we find the clusters of zero 2x2 subdeterminants of the two given rows r1,r2 of ee.
 * (The zero subdeterminants form cliques if viewed as edges of a graph!)
 * The same is computed for columns if tr==1.
 * The result is stored in ar[] -- ar[i] contains the index of the first column j<i for
 * which the subdeterminant of columns i,j is zero, or ar[i]=i if there is no zero.
 * The columns of both-0 entries are treated specially - they form their own cluster.
 * For use only in estruct_dispminor_rowhash().
**/

void	struct_dispminor_linedet(ematrix *ee, int tr, int r1, int r2, int ar[]) {
	int	i,j,p, i0;
	sign_t	gg;
	
	if (tr)  ematrix_transpose(ee);
	i0 = -1;	/* (i0 keeps the first both-0 column, ar[0] starts with 0 automatically in the cycle...) */
	for (i=0; i<COLSM(ee); i++) {
		if (SIGNM(ee,r1,i)==0 && SIGNM(ee,r2,i)==0) {
			ar[i] = i0 = (i0==-1? i:i0);
			continue;	/* (both-0 columns need special treatment!) */
		}
		for (j=0; j<i; j++)  if (ar[j]==j && ar[j]!=i0) {
			p = ematrix_determinant_2x2ch(0,ee,r1,r2,j,i,NULL,&gg);
			if (p>=0 && gg==0)  break;
		}	/* (the subdeterminant may not be defined after endomorphisms, but then like nonzero) */
		ar[i] = j;
	}
#ifndef FASTPROG
	if (IFRANDDEBUG(44))  for (i=0; i<COLSM(ee); i++)
		if (ar[i]<0 || ar[i]>i)  {PROGERROR("something wrong with linedet");}
	if (CURDLEV+4<=DEBUGLEV) {
		SDEBUG(CURDLEV,"\t\t\tlinedet at %d,%d: ar[] = ",r1,r2);
		for (i=0; i<COLSM(ee); i++)  SDEBUG(CURDLEV," %d,",ar[i]);
		SDEBUG(CURDLEV,"\n");
	}
#endif
	if (tr)  ematrix_transpose(ee);
}

/**
 * Here we examine the zero-entry and zero-subdeterminant (see above) patterns zm1,zm2,
 * zb1,zb2 and *zpm[],*zpb[] that were previously taken from two row-pairs of two matrices.
 * (The zero-subdeterminant patterns are given as fields *zpm[],*zpb[] indexed first by
 * the index (<1+maxabsvalix) used to compute "absolute values" of entries for determinants
 *  - this is provided to reveal more differences between the two rows.)
 * The patterns zm1,zm2,zb1,zb2 simply have bits 1 for nonzero entries, the meaning
 * of *zpm[],*zpb[] is the same as described in estruct_dispminor_linedet().
 * 
 * If there is a possibility of a mapping from all cm columns of the minor to a subset cl[]
 * of cm (out of all cb) columns of the big matrix, the we record in bitfield ih[j]
 * the positions of the j-th elements of the subsets cl[].
 * "Possibility" is tested for the numbers of zeros in the single rows (rows are ordered here!),
 * and for the numbers of zero 2x2 subdeterminants in the row pairs.
 * (Both these numbers are permutation- and scale-invariants for proper "absolute values".)
 * For use only in estruct_dispminor_rowhash().
**/

void	struct_dispminor_fitldet(int cm, long zm1, long zm2, int *zpm[],
				int cb, long zb1, long zb2, int *zpb[], long ih[]) {
	int	j,k,l, pb, cl[MXMN], mm[MXMN];
	int	nzm1,nzm2, nzb1[MXMN],nzb2[MXMN], nzpm[1+STR_MAXABSVALIX];
	
	if (cm>cb || cb>MXMN)  {PROGERROR("something wrong with the numbers of columns %d > %d > %d",cm,cb,MXMN);}
	if (zm1==0 || zm2==0 || zb1==0 || zb2==0)  {PROGERROR("something wrong - rows may not be all zeros %lX,%lX,%lX,%lX",zm1,zm2,zb1,zb2);}
	
			/* initialize the output hash, and count the numbers of zeros in the minor rows */
	for (nzm1=nzm2 = j=0; j<cm; j++) {
		ih[j] = 0l;
		if ((zm1&(1l<<j))==0)  nzm1++;
		if ((zm2&(1l<<j))==0)  nzm2++;
	}		/* also count the numbers of zero-subdeterminants in the minor rows by l */
	for (l=0; l<1+maxabsvalix; l++) {
		for (j=0; j<cm; j++)  mm[zpm[l][j]] = 0;
		nzpm[l] = 0;
		for (j=0; j<cm; j++)  nzpm[l] += mm[zpm[l][j]]++;
	}
	k = 0;  cl[0] = -1;
	while (k>=0) {
			/* here we generate all cm-elem subsets of cb columns for the hash computation */
		++cl[k];
		if (cl[k]>=cb-cm+1+k) { --k;  continue; }
			/* we compute the numbers of zeros in the selected columns along the way */
		nzb1[k] = (k>0? nzb1[k-1]:0) + ((zb1&(1l<<cl[k]))==0);
		nzb2[k] = (k>0? nzb2[k-1]:0) + ((zb2&(1l<<cl[k]))==0);
		if (nzb1[k]>nzm1 || nzb2[k]>nzm2)  continue;
			/* while not all columns are selected, we continue further */
		if (k<cm-1) {
			++k;  cl[k] = cl[k-1];  continue;
		}
			/* the final test whether the numbers of zeros agree for the selection */
		if (nzb1[k]!=nzm1 || nzb2[k]!=nzm2)  continue;
		
			/* now we compare numbers of zero subdeterminants (at index l=0..maxabsvalix) */
		for (l=0; l<1+maxabsvalix; l++) {
			for (j=0; j<cm; j++)  mm[ zpb[l][cl[j]] ] = 0;
			pb = 0;
			for (j=0; j<cm; j++)  pb += mm[ zpb[l][cl[j]] ]++;
			if (pb!=nzpm[l])  break;	/* (the numbers must be equal at all l's to pass) */
		}
		if (l<1+maxabsvalix)  continue;
		
		for (j=0; j<cm; j++)	/* if there was no contradiction above, we record the subset */
			ih[j] |= (1l<<cl[j]);
	}
#ifndef FASTPROG
	if (IFRANDDEBUG(140)) {
		SDEBUG(CURDLEV+4,"\t\t\t\tfitldet ih[] = ");
		if (CURDLEV+4<=DEBUGLEV)  for (j=0; j<cm; j++)  SDEBUG(CURDLEV+4," %lX,",ih[j]);
		SDEBUG(CURDLEV+4,"\n");
		for (j=0; j<cm; j++)  if (ih[j]<0 || ih[j]>=(1l<<cb))  {PROGERROR("Wrong output in fitldet");}
	}
#endif
}



#endif








































