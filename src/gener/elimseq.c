
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
 * This part defines the basic handling functions for the "elimination sequence" structure.
 * See in gen.h for the definition of "elimseqm".
 * The elimination sequence is an internal structure for use only in this source directory.
 * See in gener.c how the elim sequence is extracted from an external representation
 * of a matroid (a frame).
 * 
 * 
 * Formally, an elimination sequence consists of the base minor S, of the resulting matroid M,
 * the order of lines of M as they are added to S, and the "signature" of the sequence telling
 * us which lines of M are extended and which are co-extended.
 * (The sequence signature is taken separate from the order since the order naturaly follows
 * from the matrix representation of M, while the signature does not.)
 * The base minor S is displayed in the matrix representation of M.
 * The sequence signature is encoded by a sequence of bits where 1 means that a column
 * is extended, and 0 that a row is co-extended (against the 0-trasposition state).
 * The lines closest to the minor M are taken for the highest bits.
 * 
 * An elimination sequence can have several more attributes controlling what matrices may
 * be generated from the sequence and how (including restrictions on steps of the sequence).
 * Examples include forbidden minors, tight majors, no large fans along, etc.
 * Each generated sequence should be produced "canonically minimal", as described in gener.h.
 * The tests for comparing sequences are provided at the end of this file.
 * 
**/


#include "macek.h"
#include "gen.h"













/******************	Handling elimination sequence	*******************/
/**************************************************************************/


#undef CURDLEV
#define CURDLEV		7



/**
 * Here we create a new empty elimination sequence.
 * (It is enough to "free()" the sequence at the end.)
 * The function returns the new sequence esq.
**/

elimseqm*	new_elimseqm(void) {
	elimseqm	*esq;
	
	esq = MMALLOC(sizeof(*esq));
	esq->sz = sizeof(*esq);
	ESELIMSG(esq) = 0;
	ESMATRIX(esq) = ESMINOR(esq) = NULL;
	ESTIGHT(esq) = ESFORBID(esq) = NULL;
	ESPREFIXM(esq) = NULL;
	ESNOFAN(esq) = ESNOTAIL(esq) = -1;
	ESCONNECT(esq) = ESDEFCONN;
	if (!ESREQ3CONN(esq) || !ESREQSIMPLE(esq) || !ESREQCOSIMPLE(esq))  {PROGERROR("by default, sequences must be 3-connected, simple, and cosimple");}
	ESNAME(esq) = NULL;
	return esq;
}


/**
 * This function sets a new resulting matrix and a new elimination seq signature for
 * the sequence q, and it kepps all other data in q intact.
 * Be careful when calling with a new minor emin!=NULL -- the old minor
 * is not freed now, and so memory leaks may result when calling the function.
 * The function returns the same q as given in.
**/

elimseqm*	elimseqm_setfor_ext(int ch, elimseqm *q, ematrix *e, ematrix *emin, eseqs_t qs) {
	
	if (emin && ESMINOR(q))  {PROGERROR("setting an elim seq %p for another base minor %p",q,emin);}
	ESMATRIX(q) = e;
	if (emin)  ESMINOR(q) = emin;
	ESELIMSG(q) = e? qs:0;
	if (e) if (ESLENGTH(q)>GEN_MAXELIMLENGTH)  {PROGERROR("cannot work with elim seq longer than %d > %d",ESLENGTH(q),GEN_MAXELIMLENGTH);}
#ifndef	FASTPROG
	if (e && ch>=0 && IFRANDDEBUGLESS(111)) {
		int  tri = ISTRANSPM(e);  ematrix  *ec;
		tri = ISTRANSPM(e);  ematrix_transpose_set(e,0);
		ec = ematrix_refer(e,0,ROWSM(ESMINOR(q)),0,COLSM(ESMINOR(q)));
		if (!ematrix_isequal(ec,ESMINOR(q)))  {PROGERROR("setting an elim seq %p for a wrong base minor",q);
				EMATDEBUG(0,e,"\t");EMATDEBUG(0,ec,"\t"); EMATDEBUG(0,ESMINOR(q),"\t\t");}
		ematrix_transpose_set(e,tri);  dispose_ematrix(ec);
	}
	if (e && ch>=0 && (ch>0 || IFRANDDEBUG(33)))  if (elimchm_signbits_check(q)<0)
		{PROGERROR("setting an elim seq %p for a wrong seq signature %lX",q,qs);}
#endif
	return q;
}


/**
 * This function creates a copy of the given elimination sequence structure q,
 * pointing now to the resulting matrix e and with the sequence signature qs.
 * All other stored or refered information from q is identically set in the new structure.
 * It is supposed that e,qs are compatible with the current q structure.
 * An error is reported when sz==1 and the matrix e has different dimensions than in q.
 * The new sequence is returned.
 * (The old one may be then freed, but do not free the supplementary parameters!)
**/

elimseqm*	elimseqm_copyfor_ext(int sz, elimseqm *q, ematrix *e, eseqs_t qs) {
	elimseqm	*esq;
	
	if (q) {
		if (sz>0 && e && ESMATRIX(q))
		  if (ROWSM(ESMATRIX(q))!=ROWSM(e) || COLSM(ESMATRIX(q))!=COLSM(e)
			  		|| ISTRANSPM(ESMATRIX(q))!=ISTRANSPM(e))
			{PROGERROR("cannot copy an elim seq %p for a different size matrix %p",q,e);}
		esq = MMALLOC(sizeof(*esq));
		memcpy(esq,q,sizeof(*esq));
	} else {
		esq = new_elimseqm();
	}		/* (the base minor is not changed - the same one is referred from the new sequence) */
	return elimseqm_setfor_ext(-1,esq,e,NULL,qs);
}

/**
 * The function elimseqm_next(q,e,tr) is defined in the header using elimseqm_copyfor()...
 * It stores the "next step" of the elim sequence q with resulting matrix e and step by tr.
**/




/**
 * This function prints information about the elimination sequence structure q.
 * If ch>0, then the resulting matrix of q is printed; if ch>1, then the base minor is printed
 * as well; and if ch>2, then all tight and forbidden minors are printed, too.
 * The printout is prefixed with the string pref.
**/

void	elimseqm_fprint_ext(FILE *fo, int ch, elimseqm *q, char *pref) {
	int	i,t;
	char	buf[200];
	
	if (!pref)  pref = printoutpref;
	t = ISTRANSPM(ESMATRIX(q));
	ematrix_transpose_set(ESMATRIX(q),0);
	fprintf(fo,"%s=====================================================================\n%s"
		"elimseq %p [%s];  mat %p (tr%d), base %dx%d, sign %lX,\n%s"
		"     nofan %d, notail %d, conn%d%s, forbid# %d, tight# %d\n",
		pref,pref,q,ESNAME(q)?ESNAME(q):EMNAME(ESMATRIX(q)),ESMATRIX(q),t,
		ROWSM(ESMINOR(q)),COLSM(ESMINOR(q)),ESELIMSG(q),
		pref,ESNOFAN(q),ESNOTAIL(q),ESCONNECT(q),ESREQ3CONN(q)?"":" (not 3-conn),",
                alist_getlength(ESFORBID(q)),alist_getlength(ESTIGHT(q)));
	
	strcat(strncpy(buf,pref,100)," ");
	if (ch>0)  ematrix_fprint_pref(fo,ESMATRIX(q),buf);
	strcat(buf,">\t");
	if (ch>1)  ematrix_fprint_nofr(fo,ESMINOR(q),buf);
	if (ch>1)  fprintf(fo,"%s\n",pref);
	buf[strlen(buf)-2] = 0;  strcat(buf,">tigh\t");
	if (ch>2 && ESTIGHT(q))  for (i=0; ESTIGHT(q)[i]; i++)
		ematrix_fprint_nofr(fo,ESTIGHT(q)[i],buf);
	if (ch>2)  fprintf(fo,"%s\n",pref);
	buf[strlen(buf)-5] = 0;  strcat(buf,">forb\t");
	if (ch>2 && ESFORBID(q))  for (i=0; ESFORBID(q)[i]; i++)
		ematrix_fprint_nofr(fo,ESFORBID(q)[i],buf);
	
	elimseqm_fprint_incl(fo,ch,q,pref);
	fprintf(fo,"%s=====================================================================\n",pref);
	ematrix_transpose_set(ESMATRIX(q),t);
}








/******************	Debugging elimination sequence	*******************/
/**************************************************************************/


#undef CURDLEV
#define CURDLEV		7


/**
 * This function debug-tests the given elimination sequence for internal consistency,
 * and it prints errors if found.
 * (There is NO testing for the canonical-minimal form or matroidal-minor properties,
 * only for pfield, 3-connectivity, unit-scale, etc.)
 * If ch<=0, then some tests (pfield, nofan, tight) are not performed - for use when generating
 * extensions before the process is finished.
 * If ch==1, then some tests (tight, forbid) are done only randomized - for faster debug
 * checking.
 * Full long tests are performed for ch>=2.
 * This function is intended for debug purposes only, not for testing generated matrices!
**/

void	test_elimseqm_ext(int ch, elimseqm *q) {
	int		i=0, er=0;
	eseqs_t		k;
	ematrix		*ee;
	elimseqm	*qn;
	
	if (!q || !ESMATRIX(q) || !ESMINOR(q))  {PROGERROREXIT("Null sequence or matrix given to test");}
	if (ESLENGTH(q)>GEN_MAXELIMLENGTH)  {PROGERROREXIT("cannot work with elim seq longer than %d > %d",ESLENGTH(q),GEN_MAXELIMLENGTH);}
	if (ISTRANSPM(ESMINOR(q)))  {PROGERROREXIT("The base minor of a sequence must not be transposed");}
	if (q->sz!=sizeof(*q))  {PROGERROREXIT("The size of elimseqm structure does not agree! %ld!=%ld",q->sz,(long)sizeof(*q));}
	
	printlev -= 2;		/* (to suppress debug printing inside the next code) */
	ee = ematrix_copy(ESMATRIX(q));
	qn = elimseqm_copyfor(q,ee);
	ematrix_transpose_set(ee,0);		/* the sequence works against 0-transpose state! */
	DEBUG(CURDLEV+3,"Going to debug-test the elimination sequence %p (%d x %d)...\n",q,ROWSM(ee),COLSM(ee));
	if (ch>1)  ELIMQDEBUG(CURDLEV+3,q,"\t");
	
	if (ESPREFIXM(q)) if (ROWSM(ESPREFIXM(q))!=1) {
		er = 1;  PROGERROR("The elimination seq extension prefix must have one row! %d x %d",ROWSM(ESPREFIXM(q)),COLSM(ESPREFIXM(q)));
	}
	if (elimchm_signbits_check(q)<0) {
		er = 1;  PROGERROR("The elimination seq signature %lX does not correspond to the matrices!",ESELIMSG(q));
	}
	if (ch>0) {
		if (ch>1? elimchm_represented(q,NULL)<0:elimchm_represented_rand(q,NULL)<0) {
			er = 1;  PROGERROR("The resulting matrix of the sequence %p is not in the pfield!",q);
			ematrix_inpfield_printed(ee); }
	}
	if (ch>0 && ESTIGHT(q)) {
		if ((i= (ch>1? struct_hasdelcontr(ee,ESTIGHT(q)):struct_hasdelcontr_rand(ee,ESTIGHT(q)))))
			{ er = 1;  PROGERROR("The matrix is not a tight major for the list! (elem %d)",i-1); }
	}
	if (ch>0 && ESFORBID(q)) {
		i = (ch>1? struct_hasminorlist(ee,ESFORBID(q)):struct_hasminorlist_rand(ee,ESFORBID(q)));
		if (i>0)  { er = 1;  PROGERROR("The matrix has a forbidden minor #%d from the list!",i); }
	}
	k = ESELIMSG(q);
	for (i=0; i<ESLENGTH(q); i++, k/=2) {
		qn = elimseqm_setfor(qn,ee,k);
		if ((k%2)? ESREQSIMPLE(q): ESREQCOSIMPLE(q))
		  if (elimchm_3connlastline(ee,k%2)<0) {
			er = 1;  PROGERROR("The elimination sequence extension is not x-connected at pos %d !",ESLENGTH(q)-i);
		}
		//********* add elimchm_connlastline() test here...
		if (elimchm_unitlastline(0,ee,k%2,ESREQCONN(q))<0) {
			er = 1;  PROGERROR("The elimination sequence is not unit-line at pos %d !",ESLENGTH(q)-i);
		}
		if (ch>0)  if (genstep_sequencecheck_lev(-3,qn,ee,k%2)<0) {
			er = 1;  PROGERROR("The elimination sequence fails a seq check at pos %d !",ESLENGTH(q)-i);
		}
		ematrix_remove(ee, k%2?-1:ROWSM(ee)-1, k%2?COLSM(ee)-1:-1);
	}
	if (!ematrix_isequal(ESMINOR(q),ee)) {
		er = 1;  PROGERROR("The matrix does not start with the base minor!");
	}
	if ((ESREQ3CONN(q)? !struct_isconnected(ee,3):0) ||
			(ESREQCONN(q)? !struct_isconnected(ee,2):0)) {
		er = 1;  PROGERROR("The base minor of the sequence is not %d-connected!",ESCONNECT(q));
	}
	if (ch>1) {
		qn = elimseqm_setfor(qn,ee,0);
		if (genstep_sequencecheck_lev(-3,qn,NULL,-1)<0 || genstep_structurecheck_lev(-3,qn,NULL,-1)<0)
			{ er = 1;  PROGERROR("The elimination sequence fails a total seq or struct check for the base minor!"); }
	}
	dispose_ematrix(ee);  FREE(qn);
	
	if (ch>1)  if (genstep_sequencecheck_lev(-3,q,NULL,-1)<0) {
		er = 1;  PROGERROR("The elimination sequence fails a total seq check!");
	}
	if (ch>1)  if (genstep_structurecheck_lev(-3,q,NULL,-1)<0) {
		er = 1;  PROGERROR("The elimination sequence fails a total struct check!");
	}
	if (er) {
		elimseqm_fprint_full(errorout,q," !!!\t");
		PROGERROREXIT("## There was an error in the elimination sequence %p (sig %lX).",q,ESELIMSG(q));
	}
	printlev += 2;
}



























/******************	Sequence-related tests	*******************/
/******************************************************************/


#undef CURDLEV
#define	CURDLEV		7



/**
 * This function checks whether the sequence signature qs may correspond to the resulting
 * matrix of q with respect to the base minor - i.e. counts the number of 1's in the signature.
 * (Works with the 0-transp. state of ESMATRIX(q).)
 * The return value is 1 if possible, -1 if not.
**/

int	elimchm_signbits_ext(elimseqm *q, eseqs_t qs) {
	int	i,j,r=0,t;
	
	if (ESMATRIX(q)) {	/* (necessary to take transp state into account - but not of the minor!) */
		if (ISTRANSPM(ESMATRIX(q)))  t = ROWSM(ESMATRIX(q))-COLSM(ESMINOR(q));
		else  t = COLSM(ESMATRIX(q))-COLSM(ESMINOR(q));
		for (i=j=0; i<ESLENGTH(q); i++, qs/=2)
			if (qs%2==1)  j++;
		r = (j==t && qs==0);
	}
	DEBUG(CURDLEV+1+2*r," - signbits test for %p, %lX found %s\n",q,qs,r?"pass":"NO");
	return r? 1:-1;
}



/**
 * This function checks/sets the given line lo (row if tr==0, col if tr==1) in/to a unit-leading
 * form, where the "first 1 (unit)" is searched from the position st on the line lo.
 * If the matrix is disconnected, then the first nonzero entry of each component should be 1.
 * If lo==-1, then the last line is considered.
 * The action is controlled by set - set>0 means to set the unit form, set==0 only checks for unit.
 * If e is a refering matrix and set>1, then the referred matrix(!!!) by e is scaled instead of e.
 *  (That is used mainly in the canonical test - gcanon_dispseqtest_scal().)
 * If acn==1 is given, then the matrix (except lo) is silently assumed connected (one component).
 * 
 * The return value is 1 if the line is or has been set the unit-leading form, and -2 otherwise.
 * (If there is no nonzero entry at all, then the function returns 0.)
**/

int	elimchm_unitline_ext(int set, ematrix *e, int tr, int lo, int st, int acn) {
	int		i,j,r, ro,co, un,*uni=NULL;
	int		cn, *cpi=NULL,*cxpi=NULL,cpistack[64];
	sign_t		gg=1, *gun=NULL;
	ematrix		*ed,*eref;
	exp_t		xx, *xun=NULL;
	pfield_setzeroexp(&xx);
	
	if (ISREFMAT(e) && set>1)  if (ISTRANSPM(e)!=ISTRANSPM(REFMAT(e)))
		{PROGERROR("does not work for different transp state of refering matrix...");}
	if (!ISREFMAT(e) && set>1)  {PROGERROR("cannot ref-scale for no refering matrix...");}
	if (tr<0)  {PROGERROR("No last line is given in tr=%d<0",tr);}
	eref = ISREFMAT(e)? REFMAT(e):NULL;
	if (!tr)  ematrix_transpose(e);	/* (tr is applied here against the current state!) */
	if (!tr && set>1 && eref)  ematrix_transpose(eref);
	ro = ROWSM(e);  co = COLSM(e);
	if (lo==-1)  lo = co-1;
	if (st<0)  st = 0;		/* (starting index for the unit-scale test) */
#ifndef	FASTPROG
	//if (eref && ROWSM(e)+COLSM(e)>5) { EMATDEBUGS(2,eref,"\t>unitline>\t"); EMATDEBUGS(2,e,"\t\t>>>\t"); }
	if (acn && IFRANDDEBUGLESS(111)) {
		ed = ematrix_copy(e);  ematrix_remove_col(ed,lo);
		if (!struct_isconnected(ed,2)) 	{PROGERROR("The matrix, assumed connected, is actually disconnected!"); EMATDEBUGS(2,e,"\t>e>\t"); EMATDEBUGS(2,ed,"\t>disconn>\t");}
		dispose_ematrix(ed);
	}
#endif	
	if (acn)  cn = 1;
	else {		/* connectivity unknown - need the indices of connected matrix-row components */
		if (ro+co<30)  cpi = cpistack;
		else  cpi = MMALLOC((ro+co+ro+4)*sizeof(cpi[0]));
		cxpi = cpi+ro+1;
		ed = ematrix_copy(e);
		for (i=0; i<ro; i++)  SIGNM(ed,i,lo) = 0;
		//******** ematrix_remove_col(ed,lo);	we cannot remove the column since that would harm column indexing...
		cn = struct_connorder_rcxp(ed,cpi,cxpi);
		if (cxpi[lo]!=0)  {PROGERROR("Must have added (now zero) line of component index 0, but %d here!",cxpi[lo]);}
		dispose_ematrix(ed);
	}
		/**
		 * The task of unit scaling is easy on connected matrices,
		 * and so we treat that separately (assumed acn, or connorder OK).
		 * On a connected matrix, it is enough to find the first nonzero
		 * entry, and to invert the added line with it.
		**/
	if (cn>0) {
	  for (i=st; i<ro && !SIGNM(e,i,lo); i++) ;
	  if (i<ro) {
		xx = EXPM(e,i,lo);  gg = SIGNM(e,i,lo);
		un = (pfield_iszeroexp(xx) && gg==1);
		if (set>0 && !un) {	/* if not 1 and should set unit, then divide the column */
			if (set>1 && eref)
				ematrix_multiply_col(eref,GETREFMCOL(e,lo),-1,xx,gg);
			else  ematrix_multiply_col(e,lo,-1,xx,gg);
			un = 1;
		}
		r = (un? 1:-2);
	  } else  r = 0;		/* special - if no nonzero entry was found in lo */
	}
		/**
		 * The following part handles unit-scaling on disconnected matrices.
		 * Above, we have got the component indices cpi[] for rows,
		 * and cxpi[] for columns (where all-zero cols have index 0).
		 * Next, we find the first nonzeros in each row component (uni[]),
		 * and we either test them to 1, or store them (xun[],gun[]) for further
		 * inversion of the added column.
		 * The inversion is rather complicated here because of a possible
		 * refering matrix eref - we actually scale all lines of e except
		 * the added column: columns by the leading values of their components
		 * and rows by their inversions.
		**/
	else {
	  uni = cxpi+co+1;
	  xun = MMALLOC(ro*(sizeof(xun[0])+sizeof(gun[0])+7));
	  gun = (sign_t*)(xun+ro+1);
	  for (i=0; i<=ro; i++)  uni[i] = -1;
	  r = 0;
	  for (i=st; i<ro; i++)	/* finding the first nonzeros in components */
		if (SIGNM(e,i,lo)!=0 && uni[cpi[i]]<0) {
			r = 1;
			xx = EXPM(e,i,lo);  gg = SIGNM(e,i,lo);
			uni[cpi[i]] = (pfield_iszeroexp(xx) && gg==1);
			DEBUG(CURDLEV+2,"%sting unit component #%d at line %d, value %s, -> un=%d...\n",set<=0?"Tes":"Set",cpi[i],i,pfield_pvalue(13,xx,gg),uni[cpi[i]]);
			if (uni[cpi[i]]>0)  continue;
			if (set<=0) { r = -2;  break; }	/* (just checking unit scale - not good) */
			xun[cpi[i]] = xx;  gun[cpi[i]] = gg;
		}
	  if (set>0) {		/* scaling the components to invert the first nonzeros */
		for (i=0; i<ro; i++)  if (uni[cpi[i]]==0) {
			if (set>1 && eref)
				ematrix_multiply_row(eref,GETREFMROW(e,i),-1,xun[cpi[i]],gun[cpi[i]]);
			else  ematrix_multiply_row(e,i,-1,xun[cpi[i]],gun[cpi[i]]);
			DEBUG(CURDLEV+2,"multiply row %d(%d) of comp %d by 1/ %s ...\n",i,eref?GETREFMROW(e,i):i,cpi[i],pfield_pvalue(13,xun[cpi[i]],gun[cpi[i]]));
		}
		for (i=0; i<co; i++)  if (uni[cxpi[i]]==0 && i!=lo) {
			if (set>1 && eref)
				ematrix_multiply_col(eref,GETREFMCOL(e,i),1,xun[cxpi[i]],gun[cxpi[i]]);
			else  ematrix_multiply_col(e,i,1,xun[cxpi[i]],gun[cxpi[i]]);
			DEBUG(CURDLEV+2,"multiply col %d by %s ...\n",i,pfield_pvalue(13,xun[cpi[i]],gun[cpi[i]]));
		}
          }
	}
	//if (eref && ro+co>5) { EMATDEBUGS(2,eref,"\t<unitline<\t"); EMATDEBUGS(2,e,"\t\t<<<\t"); }
	if (!tr)  ematrix_transpose(e);
	if (!tr && set>1 && eref)  ematrix_transpose(eref);
	if (cpi && cpi!=cpistack)  FREE(cpi);  
	if (xun)  FREE(xun);  
	DEBUG(CURDLEV+1+2*(r>0)," - unitline (set=%d) test for %p tr=%d found %s\n",set,e,tr,r>0?"pass":(r>=0?"zero":"NO"));
	return r;
}



/**
 * This function checks the given line lo (row if tr==0, col if tr==1) against the
 * given line prefix in the first row of epr of length len.
 * If lo==-1, then the last line of e is considered.
 * The return value is 1 if the line agrees with epr, and -1 otherwise.
**/

int	elimchm_prefixline_ext(ematrix *e, int tr, int lo, ematrix *epr, int len) {
	int		i,r;
	
	if (tr<0)  {PROGERROR("No last line is given in tr=%d<0",tr);}
	if (!tr)  ematrix_transpose(e);	/* (tr is applied here against the current state!) */
	if (lo==-1)  lo = COLSM(e)-1;
	if (len<0)  len = COLSM(epr);
	
	r = 1;		/* we compare the entries of the ext line with the first row of epr */
	for (i=0; i<len && i<COLSM(epr) && i<ROWSM(e); i++)
		if (!pfield_isequal(EXPM(e,i,lo),SIGNM(e,i,lo),EXPM(epr,0,i),SIGNM(epr,0,i)))
			r = 0;
	
	if (!tr)  ematrix_transpose(e);
	DEBUG(CURDLEV+1+2*r," - prefixline (len %d: %s) test for %p tr=%d lo=%d found %s\n",
			len,gener_extprintline(0,epr,0),e,tr,lo,r?"pass":"NO");
	return r? 1:-1;
}


/**
 * This function checks whether the line lo (row if tr==0, col if tr==1) is a 3-connected
 * (co)extension to the rest of the matrix - i.e. not parallel or serial to anything before, not zero.
 * If lo==-1, then the last line is considered.
 * (This interpretation implicitly assumes that the previous matr. in the sequence is 3-connected.)
 * 
 * Actually, it is nothing tied with 3-connectivity here, the same function (in a diff
 * interpretaion of the situation) is used to check that the line lo extends a simple
 * matroid into a simple matroid...
 * 
 * The return value is 1 if 3-connected extension, and -1 if not.
**/

int	elimchm_3connline_ext(ematrix *e, int tr, int lo) {
	int	i,j,r;
	
	if (tr<0)  {PROGERROR("No last line is given in tr=%d<0",tr);}
	if (!tr)  ematrix_transpose(e);	/* (tr is applied here against the current state) */
	if (lo==-1)  lo = COLSM(e)-1;
	for (i=j=0; i<ROWSM(e); i++)	/* (needs at least 2 nonzeros to be 3-connected) */
		if (SIGNM(e,i,lo)!=0)  j++;
	if (j>=2)  r = ((i=ematrix_parallel_col(e,lo))<0);
	else  r = 0;
	if (!tr)  ematrix_transpose(e);
	DEBUG(CURDLEV+1+2*r," - 3connline (%d tr=%d) test for %p found (nz %d, par %d) %s\n",lo,tr,e,j,i,r?"pass":"NO");
	return r? 1:-1;
}

/**
 * Functions elimchm_connline_...() are defined as macros in the header
 * (using ematrix_linezero).
**/

/**
 * This function checks whether the sequence q is properly "connected" along all steps
 * of elimination (wrt. value cn, see ESCONNECT), including the base minor.
 * (This has no lesser interpretation, unlike elimchm_3connline_ext() above!)
 * The return value is 1 if properly connected, and -1 if not.
 * Used for an input sequence check.
**/

int	elimchm_xconnsequence(elimseqm *q, int cn) {
	eseqs_t		k;
	int		i,r, cnx;
	elimseqm	*qn;
	ematrix		*ee;
	
	if (ESLENGTH(q)>GEN_MAXELIMLENGTH)  {PROGERROR("cannot work with elim seq longer than %d > %d",ESLENGTH(q),GEN_MAXELIMLENGTH); return -1;}
	ee = ematrix_copy(ESMATRIX(q));
	qn = elimseqm_copyfor(q,ee);
	ESCONNECT(qn) = cn;	/* (assumed to be the same, anyway) */
	k = ESELIMSG(q);
	for (i=0, r=1; i<ESLENGTH(q); i++, k/=2) {
		qn = elimseqm_setfor(qn,ee,k);
		if ((k%2)? ESREQSIMPLE(qn): ESREQCOSIMPLE(qn))
			if (elimchm_3connlastline(ee,k%2)<0)  r = -1;
		ematrix_remove(ee, k%2?-1:ROWSM(ee)-1, k%2?COLSM(ee)-1:-1);
	}
	if (ESREQCONN(q)) {
		if (!struct_isconnected(ee,cn))  r = -1;
	}
	//********** missing tests for simple/cosimple base minor....
	dispose_ematrix(ee);  FREE(qn);
	DEBUG(CURDLEV-1+2*(r>0),"%d-conn sequence test for %p [%s] found %s.\n",cn,q,ESNAME(q),r>0?"pass":"-NO-");
	return r;
}


/**
 * Functions elimchm_represented_...() are defined as macros in the header (by inpfield).
 *  - elimchm_represented(q,e) and others...
 * We refer to the sequence q for possible further use, but we actually test the matrix e
 * and do not look at the entries of the resulting matrix in q.
**/









/******************	Sequence ordering	*******************/
/******************************************************************/


#undef CURDLEV
#define	CURDLEV		7



/**
 * This function checks whether the sequence signature in q is "smaller" than the signature qs.
 * (So far "smaller" means smaller as binary numbers, but that may change later...)
 * Returns 1 if ESELIMSG(q)">"qs, 0 if ESELIMSG(q)"=="qs, and -1 if ESELIMSG(q)"<"qs.
 * This defines the "official way" to compare two elimination sequences - the heaviest key.
 * If the signature qs is not admissible for the sequence q, then the result is undefined.
**/

int	elimchm_signbitscomp_ext(elimseqm *q, eseqs_t qs) {
	int	r;
	
	r = 0;
	if (ESELIMSG(q)!=qs)  r = (ESELIMSG(q)<qs)? -1:1;
	return r;
}


/**
 * This function compares two entries/lines/whole of matrices e1,e2, using EM_COMPARE().
 * If r1,c1>=0 and r2,c2>=0, then the entry r1,c1 of e1 is compared with the entry r2,c2 of e2.
 * If r1==-1 (c1==-1), then the whole column c1 (row r1) is used in a lexicographic
 * comparsion from the first to the last entry.
 * 
 * The return value is 0 when equal, -1 when the first is smaller, and 1 when the second is smaller.
 * If the two vectors are not of the same length; and if ch<0, then only the common part
 * of the vectors is compared (so a result of 0 need not mean really equal), otherwise error.
**/

/* EM_COMPARE(xa,ga,xb,gb) (((ga)==0&&(gb)==0)?0: ((ga)!=(gb)? (ga)-(gb): pfield_compareexp(xa,xb))) */

int	elimchm_compare_ext(int ch, ematrix *e1, int r1, int c1, ematrix *e2, int r2, int c2) {
	exp_t	x1,x2;
	sign_t	g1,g2;
	int	i1,ii1, i2,ii2, j1,jj1, j2,jj2, ix1,ix2,jx1,jx2, r;
	
	if (r1>=0 && r2>=0 && c1>=0 && c2>=0) {		/* comparing single entries */
		x1 = EXPM(e1,r1,c1);  g1 = SIGNM(e1,r1,c1);
		x2 = EXPM(e2,r2,c2);  g2 = SIGNM(e2,r2,c2);
		pfield_tostandform(&x1,&g1);  pfield_tostandform(&x2,&g2);
		r = EM_COMPARE(x1,g1,x2,g2);
		if (r!=0)  return (r<0?-1:1);
	} else {					/* comparing whole lines */
		i1 = (r1>=0? r1:0);  ii1 = (r1>=0? r1+1:ROWSM(e1));
		j1 = (c1>=0? c1:0);  jj1 = (c1>=0? c1+1:COLSM(e1));
		i2 = (r2>=0? r2:0);  ii2 = (r2>=0? r2+1:ROWSM(e2));
		j2 = (c2>=0? c2:0);  jj2 = (c2>=0? c2+1:COLSM(e2));
		if ( ii1-i1!=ii2-i2 || jj1-j1!=jj2-j2 ) {
			if (ch>=0)  {PROGERROR("cannot compare vectors of different lengths!");}
			else  DEBUG(CURDLEV-1,"comparing different length vectors %d,%d against %d,%d"
					" - really wanted?\n",ii1-i1,jj1-j1, ii2-i2,jj2-j2);
		}
		for (ix1=i1,ix2=i2; ix1<ii1 && ix2<ii2; ix1++,ix2++)
			for (jx1=j1,jx2=j2; jx1<jj1 && jx2<jj2; jx1++,jx2++) {
				x1 = EXPM(e1,ix1,jx1);  g1 = SIGNM(e1,ix1,jx1);
				x2 = EXPM(e2,ix2,jx2);  g2 = SIGNM(e2,ix2,jx2);
				pfield_tostandform(&x1,&g1);  pfield_tostandform(&x2,&g2);
				r = EM_COMPARE(x1,g1,x2,g2);
				if (r!=0)  return (r<0?-1:1);	/* first non-equal decides */
		}
	}
	return 0;	/* (should not get here unless equal) */
}



/**
 * This function compares the matrices of two given elimination sequences q1 and q2,
 * proceeding line by line as they are added along their sequences (against 0-trasp).
 * It assumes that the base minors and the sequence signatures are the same for both.
 * Moreover, the lines added along the sequence must be scaled to unit-leading form,
 * cf. genstep_linecheck_(), elimchm_unitline_() functions.
 * (There is a debug check for that property, but no regular check, so it must be given such.)
 * 
 * The return value is 0 when equal, -1 when the first is smaller, and 1 when the second
 * is smaller.
 * If len>=0 is given, then only the subsequences of length len are compared.
 * 
 * It defines the "official way" to compare two matrices of two elimination sequences
 * (see gener.h for general description) line by line (lines scaled as described), where each
 * single line is (and must be!) compared using elimchm_compare_...() above.
**/

int	elimchm_seqcompare_ext(int ch, elimseqm *q1, elimseqm *q2, int len) {
	int		i,il,xr,xc, t1,t2, r = 0;
	ematrix		*re1,*re2;
	
	junk = ch;	/* works against the 0-transp state of both matrices */
	t1 = ISTRANSPM(ESMATRIX(q1));  ematrix_transpose_set(ESMATRIX(q1),0);
	t2 = ISTRANSPM(ESMATRIX(q2));  ematrix_transpose_set(ESMATRIX(q2),0);
			/* for adding the lines of along the sequence using refering matrices */
	re1 = ematrix_refer(ESMATRIX(q1),0,ROWSM(ESMINOR(q1)),0,COLSM(ESMINOR(q1)));
	re2 = ematrix_refer(ESMATRIX(q2),0,ROWSM(ESMINOR(q2)),0,COLSM(ESMINOR(q2)));
#ifndef	FASTPROG
	if (len<0)  if (ROWSM(ESMATRIX(q1))!=ROWSM(ESMATRIX(q2)) || COLSM(ESMATRIX(q1))!=COLSM(ESMATRIX(q2)))
		{PROGERROR("cannot compare seq matrices of different dimensions");}
	if (len<0)  if (ESLENGTH(q1)!=ESLENGTH(q2) || ESELIMSG(q1)!=ESELIMSG(q2))
		{PROGERROR("cannot compare seq matrices for different signatures");}
	if (len>=0)  if (ESLENGTH(q1)<len || ESLENGTH(q2)<len ||
			(ESELIMSG(q1)>>(ESLENGTH(q1)-len))!=(ESELIMSG(q2)>>(ESLENGTH(q2)-len)))
		{PROGERROR("cannot partially compare seq matrices for incompatible signatures");}
	if (ch>=0 && IFRANDDEBUGLESS(111))  if (!ematrix_isequal(re1,re2))
		{PROGERROR("cannot compare seq matrices with different base minors");}
	if (ch>=0 && IFRANDDEBUGLESS(222)) {	/* (tests unit scale here as well, but the matrices may not be represented) */
		test_elimseqm_notough(q1);  test_elimseqm_notough(q2); }
#endif
	
	r = 0;			/* here we test the lines as they are added along the sequence */
	xr = ROWSM(ESMINOR(q1));  xc = COLSM(ESMINOR(q1));
	il = (len>=0? (ESLENGTH(q1)-len):0);
	for (i=ESLENGTH(q1)-1; i>=il && r==0; i--) {
		if ((ESELIMSG(q1)&(1l<<i))==0) {
				/* the result is determined by comparing the related lines lexicographically */
			ematrix_refadd_row(re1,xr);  ematrix_refadd_row(re2,xr);
			r = elimchm_compare_row(re1,xr,re2,xr);
			xr++;
		} else {
			ematrix_refadd_col(re1,xc);  ematrix_refadd_col(re2,xc);
			r = elimchm_compare_col(re1,xc,re2,xc);
			xc++;
		}		/* (the first difference r!=0 decides) */
	}
#ifndef	FASTPROG
	DEBUG(CURDLEV+1,"   - seq compare for sequences %p %p has found  %d\n",q1,q2,r);
	if (r==0 && len<0)  if (ROWSM(re1)!=ROWSM(ESMATRIX(q1)) || COLSM(re1)!=COLSM(ESMATRIX(q1)) ||
				ROWSM(re2)!=ROWSM(ESMATRIX(q2)) || COLSM(re2)!=COLSM(ESMATRIX(q2)))
		{PROGERROR("something wrong - incomplete sequence computation ???");}
	if (r==0 && len<0 && IFRANDDEBUGLESS(222)) if (!ematrix_isequal(re1,ESMATRIX(q1)) || !ematrix_isequal(re2,ESMATRIX(q2)))
		{PROGERROR("something wrong - bad (sub)sequence refering ?!?");}
	if (xr>ROWSM(ESMATRIX(q1)) || xc>COLSM(ESMATRIX(q1)))
		{PROGERROR("something wrong - got past the compared matrix ?!?");}
#endif
	dispose_ematrix(re1);  dispose_ematrix(re2);
	ematrix_transpose_set(ESMATRIX(q1),t1);
	ematrix_transpose_set(ESMATRIX(q2),t2);
	return r;
}




































