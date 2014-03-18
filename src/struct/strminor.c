
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
 * Some theory about bases and displayed minors for start...
 * Read general description of "structural" functions on the top of ../include/struct.h.
 * 
 * The word "minor" is used in its matroidal sense, but that has little meaning in
 * matrices alone.
 * Thus we speak about a "displayed minor" which means a submatrix displayed in
 * a suitable representation obtained by a choice of a basis (via pivoting and scaling).
 * (We shortly say a "dispminor".)
 * A "represented minor" is a minor that may be displayed in some basis.
 * Keep in mind that when asking about an abstract minor in the given matrix,
 * one has to give a list of all representations of the minor up to unlabeled strong
 * equivalence.
 * 
 * In general, there are two ways how to search for a represented minor --
 * we may look at dispminors among all matrices for all choices of bases for the matroid,
 * or we may try all deletions/contractions in the matrix to get an equivalent representation.
 * For the first way, we implement a function for finding dispminors, using the
 * above line maps of matrices.
 * For the second one, we implement a (fast) equivalence testing function that
 * uses various matroid-invariant properties of a matrix for faster rejections.
 * (However, in some bad cases we have to search through all bases anyway...)
 * 
 * A "tight major" N of F is defined with respect to a (finite) family F of matroids.
 * N is a tight major of F if N has a minor from F and no element of N may be both contracted
 * or deleted keeping some minor isomorphic to one of F.
 * The function struct_hasdelcontr...() is provided for testing a tight major against
 * the list in tghl (in the opposite meaning of the return value).
 * 
**/




#include "macek.h"
#include "str.h"












/******************	Matrix handling	*******************/
/**********************************************************/


#undef CURDLEV
#define	CURDLEV		6		/*7 (need <=DEBUGLEV+2 to perform paranoic tests) */


/**
 * The following function generates all displayed bases for the matroid given by e.
 * (This means that all possible multiple pivotings are applied to copies e.)
 * The new matrices equivalent to e are returned in a NULL terminated list of matrix copies,
 * in the same transposition state as given e.
 * If two distinct bases happen to produce the same matrix, then it occurs twice.
 * (Remember to free the list and the matrices later.)
 * 
 * If ch>0, then ........ nothing happens now.
 * If rnd>=1, then only a random sublist of bases is returned - for quick randomized pre-checks.
 * If rnd>=2, then the random sublist contains only bases that are at least 3 lines from the
 * current basis.
 * 
 * The list returned from here is then supposed to be used in other structural functions.
 * 
 * To save comp time, a list of all bases may be "recycled", i.e. stored for possible next call
 * for the same list (this must not be used for randomized bases or altered lists).
 * !!! Recycled bases are not guaranteed to have the same id labels for lines !!!  
**/

ematrix**	struct_minorbases_ext(int ch, ematrix *e, int rnd) {
	ematrix		**list;
	
	if (!e)  return NULL;
	//******** try recycled here... (but not for rnd<0 !!!)
	
	if (rnd>0) if (ROWSM(e)<=2 || COLSM(e)<=2)  rnd = -1;
	if (rnd<=0)  list = ematrix_getbases(e);
	else  list = (rnd>1? ematrix_getbases_randl(e): ematrix_getbases_rand(e));
	
	if (rnd==0 && STR_CACHEDBASES>0) {
		//********** allow recycling here..... only !rand !!!!!! and even when recycled!
	}
	DEBUG(CURDLEV+1,"Generated %d %s bases of the matrix %p [%s].\n",alist_getlength(list),(rnd>0?"rand":"all"),e,EMNAME(e)?EMNAME(e):"");
	return list;
}

#if STR_CACHEDBASES<=0
void	struct_minorbases_recycle(ematrix **bl) {
	dispose_alist_mats(bl);
}
void	struct_minorbases_finish(void) {
}
#else
*************	not yet implemented, would be really useful?
#endif



/**
 * Additional stuff necessary for recycling matrix self-maps, see below.
 * Notice that we must test equality for matrices including the transposition state
 * since the line maps depend on it.
**/

static reusedescr	reuseselfmaps = {
		AL_SELFMAPS,  "self-maps",  NULL,
		1, {(void*(*)(void*))&ematrix_copy},
		   {(void(*)(void*))&dispose_ematrix}, {(int(*)(void*,void*))&ematrix_isequaltr}
};
void	struct_matrixselfmaps_recycle(void *fls) {
	alist_reuse_rec(&reuseselfmaps,0,fls);
}
void	struct_matrixselfmaps_finish(void) {
	alist_reuse_null(&reuseselfmaps);
}

/**
 * This function generates all self- line maps (i.e. row and column permutations that,
 * possibly after rescaling, identically preserve all entries) of the matrix e.
 * If scal==0, then rescaling after permuting is not allowed.
 * Read the description of structure emlinemap in str.h for general information about
 * line maps.
 * It is enough to "free(a)" to dispose a line-map structure a.
 * There is a possibility of "recycling" the maps, but this should be used with caution.
 * 
 * The return value is the number of maps found, these maps are stored in the list *fls.
 * (We use a void** type since we do not export the declaration of emlinemap...)
 * This function is quite fast unless the matrix has too many symmetries.
 * If ch>=2, then all maps are printed out when they are found (see in strmap_samecols_() if ch>=3).
 * If ch==1, then the number of all maps is returned even if fls==NULL.
**/

int	struct_matrixselfmaps_ext(int ch, ematrix *e, int scal, void ***fls) {
	void	**x,**z;
	int	p;
	long	hm;
	
	if (!e)  {PROGERROR("The matrix e must be given here!"); return -1;}
	if (fls?*fls:0)  {PROGERROR("The list fls must be given initialized to NULL!");}
	if (ch>2)  OUTPUT("  Printing all (%s) self- line maps (\"automorphisms\") of the matrix %p [%s]:\n",
				scal?"scaled":"unscaled",e,EMNAME(e)?EMNAME(e):"");
	if (ch>2)  EMATOUTPUT(e,printoutpref);
	DEBUG(CURDLEV+1,"Generating (%s) self- line maps of the matrix %p[%s] %dx%d.\n",scal?"scaled":"unscaled",e,EMNAME(e)?EMNAME(e):"",ROWSM(e),COLSM(e));
	EMATDEBUGS(CURDLEV+3,e,"\t\t\t\t");
	z = NULL;
	idsparam.id[0] = e;
	hm = strmag_matrixentries(e)+ROWSM(e);
	if (ch==0 && scal) {	/* trying to get a recycled list for the matrix e */
		z = alist_reuse_get(&reuseselfmaps,0,idsparam,hm);
		if (z) if (!z[0]?1: strmap_check_size(z[0],e,e)<0) {
			EMATDEBUGS(0,e," !\t"); if (z[0]) EMLMAPDEBUGS(0,z[0],e,"   !\t"); PROGERROREXIT("wrong recycled self-map returned!"); }
	}
				/* here we call the function that generates all self-maps: */
	if (!z)  p = strmap_samelines_ch((ch<=2?0:ch),e,e,scal,(emlinemap***)(&z));
	else  p = alist_getlength(z);
	if (p<=0)  {PROGERROR("The identical self-map must be generated here! ch=%d, p=%d",ch,p);}
	if (ch>=2)  OUTPUT("  Found %d (%s) self- line maps of the matrix [%s] (%dx%d).\n",
				p,scal?"scaled":"unscaled",EMNAME(e)?EMNAME(e):"",ROWSM(e),COLSM(e));
	if (ch==2)  for (x=z; x?*x:0; x++)  EMLMAPOUTPUT(*x,e,NULL);
	if (fls)  *fls = alist_applist(*fls,z);
	else if (z)  dispose_alist(z);
	
	if (fls && ch==0 && scal) {	/* the generated list must be registered for recycling */
		alist_reuse_set(&reuseselfmaps,idsparam,hm,*fls);
	}
#ifndef	FASTPROG
	DEBUG(CURDLEV+(p<=1),"Generated %d (%s) self- line maps of the matrix %p [%s] %dx%d.\n",
			p,scal?"scaled":"unscaled",e,EMNAME(e)?EMNAME(e):"",ROWSM(e),COLSM(e));
	if (ch!=1 && ch>=0 && IFRANDDEBUGLESS(222)) {
		ematrix	*ee = ematrix_copydual(e);  EMSETNAME(ee,"rec-test");
		ematrix_swaprows(ee,0,ROWSM(ee)-1);  ematrix_swapcols(ee,0,COLSM(ee)-1);
		if (p!=struct_matrixselfmaps_number(ee))  {PROGERROR("non-symmetric result for modified transposed matrix");}
			/* (struct_matrixselfmaps_number is called with ch==1, so recursive call is blocked) */
		dispose_ematrix(ee); }
#endif
	return p;
}


/**
 * This function compares the two given matrices for ``equality'' up to scaling (scal==1),
 * or up to matrix row maps (fixr==1), or up to matrix col maps (fixc==1).
 * The generated maps may be stored in *fls.
 * There is no use of ch so far...
 * Read the description of structure emlinemap in str.h for general information about
 * line maps.
 * The return value is >0 for ``equality'', and -1 for not.
**/

int	struct_mapmatrix_ext(int ch, ematrix *e1, ematrix *e2, int fixr, int fixc, int scal, void ***fls) {
	long	*hr, hstack[50];
	int	i,p, ro,co;
	
	if (ROWSM(e1)!=ROWSM(e2) || COLSM(e1)!=COLSM(e2)) {
		DEBUG(CURDLEV-3,"Do not call for different-size matrices!\n");
		return -1;
	}
	ro = ROWSM(e1);  co = COLSM(e1);
	if (ro+co<=46)  hr = hstack;
	else  hr = MMALLOC((ro+co+4)*sizeof(hr[0]));
	for (i=0; i<ro+co; i++)  hr[i] = i;
	p = strmap_samelines_hh(e1,e2,scal,(fixr?hr:NULL),(fixr?hr:NULL),
			(fixc?hr:NULL),(fixc?hr:NULL),(emlinemap***)fls);
	if (hr && hr!=hstack)  FREE(hr);
	return (p<0? -1:1);
}


















/******************	Matroid minors	*******************/
/**********************************************************/


#undef CURDLEV
#define	CURDLEV		6		/*6 (need <=DEBUGLEV+2 to perform paranoic tests) */


/**
 * These variables are used to collect statistical data about minor testing.
 * They are zeroed in struct_minorbybas_() or in struct_minorbyrem_(), and printed
 * upon the end of the computation, or after significant comp steps.
 * These variables are used here and in strlmap.c.
**/

int	sm_callfixrows, sm_callfixrows2;
int	sm_callsamecols, sm_callsamecols2;
int	sm_calldispmin, sm_calldispmin2, sm_calldispmin3;
int	sm_callmateqiv, sm_callmateqiv2, sm_callmateqivb;


/**
 * This function looks in the matrix e for displayed minors strongly-equivalent to emin.
 * In other words, it searches through all ordered subsets of rows and columns corresponding
 * to the dimensions of emin, whether the respective submatrix is a perm-scale of emin.
 * Possible self-maps of emin given in afo are factored out
 * (so that one displayed minor comes out only once if all self-maps are given).
 * 
 * The return list *mls consists of refering (to e) matrices showing the emin-minor
 * exactly as it corresponds to the entries of emin, up to scaling of course.
 * The list *mls is NULL terminated, and *mls must be initialized prior to calling this function.
 * (Remember to free the list and the matrices later.)
 * Additionally, the EMDATA(re) field of these matrices stores the line map that generated
 * this minor, and it is automatically freed when disposing the matrix.
 * If emin is empty or bigger than e, then nothing is returned in the list.
 * 
 * If the arrays hrem,hre or hcem,hce or both are given, then only those displayed minors
 * that are consistent with these labels (hash codes) are produced.
 * Namely, a column j of e may represent a column jj of emin in a displayed minor
 * only if hcem[jj]==hce[j].
 * If rnd>0, then only randomized test is performed, but this has no effect here now...
 * 
 * If ch>=2(3,4), then the minors are printed as they are found.
 * The return value is -1 for no dispminor, or 0 or the number of minors found here.
**/

int	struct_dispminor_ext(int ch, ematrix *e, ematrix *emin, void **afox, ematrix ***mls,
				long hre[], long hrem[], long hce[], long hcem[], int rnd) {
	long		*hhx;
	int		r, ro,co,rm,cm, db2,db3, spec = 0;
	emlinemap	**afo, **als,***alsp=NULL, **x;
	ematrix		*ea;
	
	if (!e || !emin)  {PROGERROR("both e, emin must be given here"); return -1;}
	ro = ROWSM(e);  co = COLSM(e);
	rm = ROWSM(emin);  cm = COLSM(emin);
	if (rm>ro || cm>co || rm==0 || cm==0)  return -1;
	if (hre && rm>3 && cm>3)  DEBUG(CURDLEV+2," Dispminor uses hash codes hrem=[%ld,%ld,%ld,%ld], hcem=[%ld,%ld,%ld,%ld] -> hre=[%ld,%ld,%ld,%ld], hce=[%ld,%ld,%ld,%ld].\n",
				hrem[0],hrem[1],hrem[2],hrem[3],hcem[0],hcem[1],hcem[2],hcem[3],hre[0],hre[1],hre[2],hre[3],hce[0],hce[1],hce[2],hce[3]);
	alsp = ((mls||ch>=1)? &als:NULL);
	als = NULL;
	sm_calldispmin++;		/* (used for debug statistics only) */
	db2 = sm_callsamecols2;  db3 = sm_callfixrows2;
		/**
		 * Here we select the best way of looking at the displayed minors,
		 * based on the dimensions of the matrices (to be faster...).
		 * If the rows or the columns of these two matrices are equal, then we
		 * call directly strmap_samecols() (after possible transpose).
		 * Otherwise, we call strmap_allbysubs(), transposed so that the number
		 * of line subsets computed there is smaller.
		 * (Matrix transpositions must also exchange given row/column codes!)
		 * 
		 *  ........ should be extended by strmap_allbyorder_() later ........
		**/
	if (afox==(void**)1) {
		if (hrem || hcem)  {PROGERROR("do not use automatic self-maps with the hash codes");}
		afo = NULL;
		struct_matrixselfmaps(emin,(void***)&afo);
	} else  afo = (emlinemap**)afox;
	if (rm==ro)  spec = 2+1;
	else if (cm==co)  spec = 2;
	else if (ro*(ro-rm)<co*(co-cm))  spec = 0+1;
	else  spec = 0;
	if (spec%2==1) {	/* (transpose to choose subset of rows first) */
		ematrix_transpose(e);
		hhx = hre;  hre = hce;  hce = hhx;
		ematrix_transpose(emin);
		hhx = hrem;  hrem = hcem;  hcem = hhx;
	}
	if (2*(spec/2)==2) {
		r = strmap_samecols(emin,e,afo,hrem,hre,hcem,hce,alsp);
	} else {
		r = strmap_allbysubs(emin,e,afo,hrem,hre,hcem,hce,alsp);
	}
	//******** for a small minor call strmap_allbyorder_()...........
	if (spec%2==1) {
		ematrix_transpose(e);
		hhx = hre;  hre = hce;  hce = hhx;
		ematrix_transpose(emin);
		hhx = hrem;  hrem = hcem;  hcem = hhx;
	}
	sm_calldispmin2 += (sm_callsamecols2>db2);
	sm_calldispmin3 += (sm_callfixrows2>db3);
		/**
		 * The displayed minors are returned as matrix line maps, and they are
		 * converted into refering (to e) submatrices here.
		 * This (als!=NULL) happens if the maps were requested above (for mls or for ch>=1),
		 * and if some maps were really found.
		 * The generated maps from als are stored in the submatrices in *mls!
		**/
	if (als) {
		if (mls)  DEBUG(CURDLEV,"Converting %d found dispminor line-maps to submatrices.\n",alist_getlength(als));
		if (ch>=4)  SOUTPUT("\n");
		if (ch==2 || ch==3)  OUTPUT("   Found a [%s]-minor (%dx%d) displayed in the matrix [%s] (%dx%d).\n",
				EMNAME(emin)?EMNAME(emin):"",ROWSM(emin),COLSM(emin),EMNAME(e)?EMNAME(e):"",ROWSM(e),COLSM(e));
		if (ch>=4)  OUTPUT("  Found %d displayed [%s]-minors (%dx%d) in the matrix [%s] (%dx%d).\n",
				alist_getlength(als),EMNAME(emin)?EMNAME(emin):"",ROWSM(emin),COLSM(emin),EMNAME(e)?EMNAME(e):"",ROWSM(e),COLSM(e));
		if (ch>=3 && ch<=4)  EMATOUTPUTS(e,printoutpref);
		for (x=als; x?*x:0; x++) {
			if ((ch==3 && x==als) || ch==4)  EMLMAPOUTPUT(*x,emin,NULL);
			if (ch<=4 && !mls)  continue;
			ea = strmap_tosubmatrix(emin,e,*x);
			if (ch>=5)  EMATOUTPUTS(ea,printoutpref);
			if (mls)  *mls = alist_append(*mls,ea);
			else  dispose_ematrix(ea);
		}
		alist_free(als);
	}
#ifndef	FASTPROG
	DEBUG(CURDLEV-1+(r<0)+(ch==1),"Found %s%d displayed %p[%s]-minors (%dx%d) in the matrix %p (%dx%d).\n",
			(r<0?"-N":(als?"":">=")),(r<0?0:(!r?1:r)),emin,EMNAME(emin)?EMNAME(emin):"",ROWSM(emin),COLSM(emin),e,ROWSM(e),COLSM(e));
	DEBUG(CURDLEV+(ch==1)," (Dispminor called \"%s\"; used %d same-col tries (%d), and %d fix-row tries (%d).).\n",
			spec>1?"samecols":"allbysubs",sm_callsamecols2-db2,sm_callsamecols2,sm_callfixrows2-db3,sm_callfixrows2);
	if (r>0 && ch!=1 && ch>=0 && rnd<=0 && IFRANDDEBUGLESS(111)
			&& (!afox || afox==(void**)1) && !hre && !hce) {
		ea = ematrix_copydual(emin);  EMSETNAME(ea,"rec-test");
		ematrix_swaprows(ea,0,ROWSM(ea)-1);  ematrix_swapcols(ea,0,COLSM(ea)-1);
		ematrix_transpose(e);	/* trying dispminor in transposed matrices */
		if (r!=struct_dispminor_number(e,ea,afox))  {PROGERROR("non-symmetric result for modified transposed matrix");}
			/* (struct_dispminor_number is called with ch==1, so recursive call is blocked) */
		ematrix_transpose(e);  dispose_ematrix(ea);
	}
#endif
	if (afo && (void**)afo!=afox)  struct_matrixselfmaps_recycle(afo);
	return (r<0? -1 : ((mls||ch>=1)? r:0));
}




/**
 * This function looks in the matrix e for all (pivoted and then displayed) minors strongly
 * equivalent to the matrix emin, by examining all displayed bases of e.
 * (Read the description of struct_dispminor_ext() above...)
 * The function is for internal use only.
 * 
 * Generated minors (as refering submatrices) are appended to the given list *mls.
 * Nothing is stored (only existence test) if mls==NULL.
 * A copy of the displayed basis is stored to *mls after the dispminors refering to it.
 * The return value is -1(<0) for no minor found, or 0 or the number of minors found here.
 * 
 * Possible self-maps of emin given in afox are factored out
 * (so that one displayed minor comes out only once if all self-maps are given).
 * Use afox==(void**)1 to request generating all self-maps automatically.
 * If the arrays hhe,hhem are given, then only those minors that are consistent with
 * these labels (hash codes) are produced; hhe(m) is indexed first by rows, then by cols.
 * If rnd>0, then only randomized test is performed -- against a random sublist of bases.
 * So rnd>0 may return NO even if the minor exists, but a YES answer is always sure.
 * If the list basl is given, then only bases in this list are used.
 * 
 * If ch>=3, then the minors are printed as they are found, all of them if mls.
**/

int	struct_minorbybas_ext(int ch, ematrix *e, ematrix *emin, ematrix **basl, void **afox,
				ematrix ***mls,	long hhe[], long hhem[], int rnd) {
	ematrix		*eb, **bl,**x;
	int		i,j, p,r, ro,co,rm,cm, ppx=0;
	long		*hbe, *hbi, hbstack[50];
	char		*popf=NULL, buf[100];
	void		**afo;
	
	if (mls?*mls:0)  {PROGERROR("the list mls must be given empty (null)"); return -3;}
	ro = ROWSM(e);  co = COLSM(e);
	rm = ROWSM(emin);  cm = COLSM(emin);
	if (rm>ro || cm>co)  return -2;
	if (rm==0 || cm==0)  return 0;
	if ((hhem&&hhe)? !strmag_isinjection(hhem,rm+cm,hhe,ro+co):0) {
		DEBUG(CURDLEV-2+(ch<0),"The given labels (hash codes) allow -NO- %p[%s]-minor in %p (%dx%d).\n",
				emin,EMNAME(emin)?EMNAME(emin):"",e,ROWSM(e),COLSM(e));
		return -2;
	}
	sm_callfixrows = sm_callfixrows2 = sm_callsamecols = sm_callsamecols2 = 0;
	sm_calldispmin = sm_calldispmin2 = sm_calldispmin3 = 0;
		/**
		 * We prepare the list afo of all self- line maps of emin if requested.
		 * We also prepare the list of all displayed bases bl if it is not given.
		 * If rnd>0, we use only a randomized sublist of bases for fast pre-tests.
		 * If the line codes hhe[] are given, we must prepare for reassigning them
		 * in the displayed bases below.
		**/
	if (afox==(void**)1) {
		if (hhe || hhem)  {PROGERROR("do not use automatic self-maps with the hash codes");}
		afo = NULL;  p = struct_matrixselfmaps(emin,&afo);
		if (p<=1) { struct_matrixselfmaps_recycle(afo);  afo = NULL; }
	} else  afo = afox;
	if (basl && rnd>0)  {PROGERROR("do not use a list of bases with rnd=%d>0",rnd);}
	if (!basl) {
		if (hhe || hhem)  {PROGERROR("do not use automatic bases with the hash codes");}
		bl = (rnd>0? struct_minorbases_rn(e,rnd): struct_minorbases(e));
	} else  bl = basl;
	if (bl?!bl[0]:1)  {PROGERROR("at least one basis must be given here");}
	hbe = hbi = hbstack;
	if (hhe) {
		ppx = ematrix_checkid(e);
		if (ppx<ro || ppx<co)  {PROGERROR("cannot handle line hash codes with duplicate line id in the matrix"); ppx=ro+co;}
		if (4*ppx>44)  hbe = MMALLOC((4*ppx+4)*sizeof(hbe[0]));
		hbi = hbe+3*ppx;
		for (i=0; i<ro+co; i++)
			hbi[i<ro? ROWSID(e,i):COLSID(e,i-ro)] = hhe[i];
	}
		/**
		 * Here we look for all emin-minors in e that are displayed in some basis of e
		 * from the list bl.
		 * If the line codes hhe[] are given, then we must shift them according to
		 * the line ids in the bases, and we store the shifted codes in hbe[].
		 * 
		**/
	DEBUG(CURDLEV-1+(ch<0),"Start looking for %p[%s]-minors (%dx%d) in %p (%dx%d) on %d bases.\n",
			emin,EMNAME(emin)?EMNAME(emin):"",ROWSM(emin),COLSM(emin),e,ROWSM(e),COLSM(e),alist_getlength(bl));
	EMATDEBUG(CURDLEV+1,e,"\te:\t"); EMATDEBUG(CURDLEV+1,emin,"\tem:\t");
	r = -1;
	for (x=bl; x?*x:0; x++) {
		if (hhe)  for (i=0; i<ro+co; i++) {
			j = (i<ro? ROWSID(*x,i):COLSID(*x,i-ro));
			hbe[i] = ((j<=ppx && j>=-ppx)? hbi[j]:0l);
		}
#ifndef	FASTPROG
		DEBUG(CURDLEV+1+(ch<0),"Looking for minors (%dx%d) displayed in the basis %p of the matrix %p[%s] (%dx%d).\n",
				ROWSM(emin),COLSM(emin),*x,e,EMNAME(e)?EMNAME(e):"",ROWSM(e),COLSM(e));
		if (hhe && cm>3)  DEBUG(CURDLEV+1+(ch<0),"  (Using line hash codes hhem = [%ld,%ld,%ld,%ld,...,%ld,%ld,%ld,%ld] -> hbe[] = [%ld,%ld,%ld,%ld,...,%ld,%ld,%ld,%ld].)\n",
				hhem[0],hhem[1],hhem[2],hhem[3],hhem[rm+0],hhem[rm+1],hhem[rm+2],hhem[rm+3],hbe[0],hbe[1],hbe[2],hbe[3],hbe[ro+0],hbe[ro+1],hbe[ro+2],hbe[ro+3]);
		EMATDEBUGS(CURDLEV+2,*x,"\t\tbas:\t");
#endif
		if (ch>=3) {
			if (ch>=5)  OUTPUT("  Looking for minors displayed in the basis %p of the matrix [%s] (%dx%d).\n",
					*x,EMNAME(e)?EMNAME(e):"",ROWSM(e),COLSM(e));
			if (!popf) { popf = printoutpref;  printoutpref = buf;  snprintf(buf,96,"%s\t",popf); }
			if (ch>=5)  EMATOUTPUT(*x,buf);
		}
		if ((hhem&&hhe)? !strmag_isinjection(hhem,rm,hbe,ro):0)
			continue;
			/**
			 * Then we look for displayed emin-minors in *x, respecting possible
			 * shifted line codes in hbe against hhem.
			 * If we are going to store the generated dispminors, we must also
			 * store a copy of the current basis which is referred to (in the same list!).
			**/
		eb = (mls? ematrix_copy(*x):*x);
		if (EMNAME(e) && !EMNAME(eb))  EMSETNAME(eb,EMNAME(e));
		p = struct_dispminor_ext(ch,eb,emin,afo,mls,
				(hhe?hbe:NULL),hhem,(hhe?hbe+ro:NULL),(hhem?hhem+rm:NULL),rnd);
		r = (r<0? p: (p<0? r:p+r));
		if (mls) {
			if (p<0)  dispose_ematrix(eb);
			else  *mls = alist_append(*mls,eb);
		}		/* (only those copies of bases that display minors are stored) */
		if (mls && r>STR_MAXMINORS)  {PROGERROR("Too many minors %d>=%d are generated, quitting.\n",r,STR_MAXMINORS); break;}
		if (!mls && r>=0 && ch<=3)  break;
	}
	if (ch>=4 && r>0)  OUTPUT(" Found %d total [%s]-minors (%dx%d) in the matrix [%s] (%dx%d).\n",
			r,EMNAME(emin)?EMNAME(emin):"",ROWSM(emin),COLSM(emin),EMNAME(e)?EMNAME(e):"",ROWSM(e),COLSM(e));
	/* (print also for ch==3?) */
	
#ifndef	FASTPROG
	DEBUG(CURDLEV-2+(r<0)+(ch<0),"Found %s%d total %p[%s]-minors (%dx%d) in the matrix %p (%dx%d).\n",
			(r<0?"-N":(mls?"":">=")),(r<0?0:(!r?1:r)),emin,EMNAME(emin)?EMNAME(emin):"",ROWSM(emin),COLSM(emin),e,ROWSM(e),COLSM(e));
	SDEBUG(CURDLEV-1+(ch<0),"\t\t(Minor test used %d bases, called dispminor %d,%d,%d, mapped same-col %d,%d, fixed-row %d,%d.)\n",
			alist_getlength(bl), sm_calldispmin,sm_calldispmin2,sm_calldispmin3, sm_callsamecols,sm_callsamecols2,sm_callfixrows,sm_callfixrows2);
	if (ch>=0 && rnd<=0 && IFRANDDEBUGLESS(111) && !hhem && !basl) {
		eb = ematrix_copy(e);  EMSETNAME(eb,"rec-test");
		i = RANDOM()%ROWSM(eb);  j = RANDOM()%COLSM(eb);
		ematrix_swaprows(eb,0,ROWSM(eb)-1);  ematrix_swapcols(eb,0,COLSM(eb)-1);
		if (SIGNM(eb,i,j)!=0)  ematrix_pivot(eb,i,j);
		if ((r<0)!=(struct_minorbybas_ext(-1,eb,emin,NULL,afox,NULL,NULL,NULL,0)<0))  {PROGERROR("non-symmetric result for modified pivoted matrix");}
		dispose_ematrix(eb);
	}
	if (ch>=0 && rnd<=0 && IFRANDDEBUGLESS(555) && !hhem && !basl && ro+co-rm-cm<(rm+cm>15? 4:2))
		if ((r<0)!=(struct_minorbyrem(-1,e,emin)<0))  {PROGERROR("different result by other testing algorithm");}
#endif
	if (popf)  printoutpref = popf;
	if (hbe && hbe!=hbstack)  FREE(hbe);
	if (afo && afo!=afox)  struct_matrixselfmaps_recycle(afo);
	if (bl && bl!=basl) {
		if (rnd==0)  struct_minorbases_recycle(bl);
		else  dispose_alist_mats(bl);
	}
	return (r<0? -1 : ((mls||ch>=2)? r:0));
}


/**
 * This function looks in the matrix e for all represented minors strongly
 * equivalent to the matrix emin, by examining all removals of elements from e.
 * The function is for internal use only.
 * 
 * The return value is -1 for no minor found, or 0 for a minor.
 * This function cannot be used for finding all minors or for a randomized test.
 * It is supposed to be faster than struct_minorbybas_ext() if e is very close to emin.
 * (No larger difference than STR_MAXMINREM is allowed anyway...)
 * 
**/

#define	FLN_MINR_EQUIV	2	/* (to use the same parametr for all minor/equiv testing!) */

int	struct_minorbyrem_ext(int ch, ematrix *e, ematrix *emin, int rnd) {
	int		rma[STR_MAXMINREM+2];
	ematrix		*rme[STR_MAXMINREM+2], *ee;
	long		*fln, flstack[50];
	char		buf[70];
	int		i,k,r, ro,co,rm,cm;
	
	if (rnd>0)  {PROGERROR("random test is not implemented here"); rnd = 0;}
	if (ch>3)  {PROGERROR("output printing is not implemented here"); ch = 0;}
	ro = ROWSM(e);  co = COLSM(e);
	rm = ROWSM(emin);  cm = COLSM(emin);
	if (ro+co-rm-cm>STR_MAXMINREM)  {PROGERROR("too big size difference for minor-by-removing %d>%d",ro+co-rm-cm,STR_MAXMINREM); return -1;}
	if (rm>ro || cm>co)  return -1;
	if (rm==0 || cm==0)  return 0;
	if (EMNAME(e))  sprintf(buf,"%6.6s-%d",EMNAME(e),ro+co-rm-cm);
	else  sprintf(buf,"-mrem-%d",ro+co-rm-cm);
	DEBUG(CURDLEV-1+(ch<0),"Start looking for %p[%s]-minors (%dx%d) in %p (%dx%d) by removing -%d.\n",
			emin,EMNAME(emin)?EMNAME(emin):"",ROWSM(emin),COLSM(emin),e,ROWSM(e),COLSM(e),ro+co-rm-cm);
	EMATDEBUG(CURDLEV+1,e,"\te:\t"); EMATDEBUG(CURDLEV+1,emin,"\tem:\t");
	sm_calldispmin = sm_calldispmin2 = sm_calldispmin3 = 0;
	sm_callmateqiv = sm_callmateqiv2 = sm_callmateqivb = 0;
		/**
		 * Here we try to remove all combinations of elements from e to make it
		 * the same size as emin - trying contractions first and then deletions.
		 * The array rme[] keeps the (partially) removed matrices, and rma[]
		 * the indices of lastly removed lines for each step (against rme[]-matrices!).
		 * We consider all contractions unordered, and all deletions as well.
		**/
	for (i=0; i<STR_MAXMINREM; i++)  rme[i] = NULL;
	for (i=0; i<STR_MAXMINREM; i++)  rma[i] = -1;
	if (rm+cm<46 && rm+cm<ro+co)
		strmag_flatlines(emin,FLN_MINR_EQUIV,fln=flstack);
	else  fln = NULL;
	r = k = 0;  rma[0] = -1;
	while (k>=0) {
		if (rme[k]) {
			dispose_ematrix(rme[k]);  rme[k] = NULL;
		}
			/* (rma[] indexes the lines of the smaller rme[] matrices!!!) */
		if (++rma[k]>=ro+co-k || r!=0) {
			--k;  continue;
		}
		ee = (k>0? rme[k-1]:e);
			/**
			 * Here we are going to remove the selected line rma[k] of ee
			 *  (which may be now smaller than original e!),
			 * using contraction in the first (ro-rm) steps, and later deletion.
			 * However, if a loop/coloop is contracted/deleted, we have to
			 * discard the operation due to duplication and wrong resulting size!
			**/
		//DEBUG(0,"to remove %s (k=%d, rma=%d):  %d x %d ",k>=ro-rm?"del":"con",k,rma[k],ROWSM(ee),COLSM(ee));
		rme[k] = ematrix_removemat_rcsum(ee,(k>=ro-rm),rma[k]);
		if (k>=ro-rm? ROWSM(rme[k])<ROWSM(ee): COLSM(rme[k])<COLSM(ee))
			continue;
		//SDEBUG(0," ->  %d x %d\n",ROWSM(rme[k]),COLSM(rme[k]));
		if (k<ro+co-rm-cm-1) {
			++k;
			if (k==ro-rm)  rma[k] = -1;	/* (start new selection for deletion) */
			else  rma[k] = rma[k-1];
			continue;
		}
		DEBUG(CURDLEV+(k>0),"Minor - trying to remove (%dx%d->%dx%d): %d, %d, %d, %d ...\n",ro,co,rm,cm,rma[0],rma[1],rma[2],rma[3]);
		if (ROWSM(rme[k])!=rm || COLSM(rme[k])!=cm)  {PROGERROR("wrong removal submatrix size: (%dx%d->%dx%d): %d, %d, %d, %d  XX %dx%d",ro,co,rm,cm,rma[0],rma[1],rma[2],rma[3],ROWSM(rme[k]),COLSM(rme[k]));}
		EMSETNAME(rme[k],buf);
			/**
			 * Here we test matrix equivalence between emin and the matrix
			 * rme[k] obtained by removing lines from e.
			 * This is supposed to be a very fast test, and so the magic values
			 * for matrix lines of emin are precomputed in fln[].
			**/
		printlev -= 1;
		if (struct_matequivalence_mbr((ch>0?ch:-1),emin,rme[k],fln,NULL)>=0)
			r = 1;	/* (equivalence self-tests are disabled in this call for ch==0) */
		printlev += 1;
	}
	if (r==0 && k<0)  r = -1;
	
#ifndef	FASTPROG
	DEBUG(CURDLEV-2+(r<0)+(ch<0),"Found %s %p[%s]-minors (%dx%d) in the matrix %p (%dx%d).\n",
			(r<0?"-NO-":"+HAS+"),emin,EMNAME(emin)?EMNAME(emin):"",ROWSM(emin),COLSM(emin),e,ROWSM(e),COLSM(e));
	SDEBUG(CURDLEV-1+(ch<0),"\t\t(Minor test called %d,2%d equiv-tests, with %d bases and %d deep dispminors.)\n",
			sm_callmateqiv,sm_callmateqiv2,sm_callmateqivb,sm_calldispmin3);
	if (ch>=0 && rnd<=0 && IFRANDDEBUGLESS(111))
		if ((r<0)!=(struct_minorbybas(-1,e,emin,NULL,0)<0))  {PROGERROR("different result by other testing algorithm");}
#endif
	return (r>=0? 0:r);
}













/******************	Matroid equivalence	*******************/
/******************************************************************/


#undef CURDLEV
#define	CURDLEV		6		/*6 (need <=DEBUGLEV+2 to perform paranoic tests) */


/**
 * This function determines (unlabeled strong) matrix equivalence between the two
 * given matrices.
 * Normally, the return value says whether an equivalence between es,ed exists
 * (>=0 yes, <0 for no).
 * One may give supplementary hash codes for matrix lines in hhs[],hhd[]
 * (indexed first by rows, then by columns); but these codes are not guaranteed(!)
 * to be reflected in the resulting maps, they are used only as hints.
 * If the arrays fls[] or fld[] are given, then they are used as the resulting
 * magic values from strmag_flatlines() for es or ed, respectively.
 * 
 * If ch>=3 or mls are given, then we ask for all equivalence mappings from es to ed
 * (with no self-maps factored out), and these mappings are returned as from
 * struct_minorbybas_bl(), or printed out.
 * If ch>=3, then the minors are printed as they are found, all of them if mls.
 * If rnd>0, then only a randomized test is performed - for a sublist of bases.
 * So rnd>0 may return NO even if equivalent, but a YES answer is always sure.
 * The random return value is -10 if the NO answer is already sure.
 * 
 * The return value is -3 if the matrix sizes or given hash codes do not agree,
 * -2 if magic codes computed for matroid elements (strmagic.c) do not allow equivalence,
 * -1 if no equivalence was found among displayed bases of es to ed,
 * or >=0 if the matrices are found equivalent.
 * This function is for internal use only, call struct_isequivlist_() otherwise.
**/

int	struct_matequivalence_ext(int ch, ematrix *es, ematrix *ed, ematrix ***mls,
					long hhs[], long hhd[], long fls[], long fld[], int rnd) {
	ematrix		*ee, **bl,**bls,**bld, **x,**y;
	long		*hes=NULL,*hed,hestack[50];
	int		i,j,dd, rs,cs,rd,cd, p,q,r, *nza=NULL,*nzd,nzastack[200];
	
	if (!es || !ed)  {PROGERROR("both matrices es, ed must be given here"); return -3;}
	if (mls?*mls:0)  {PROGERROR("the list mls must be given empty (null)"); return -3;}
	rs = ROWSM(es);  cs = COLSM(es);
	rd = ROWSM(ed);  cd = COLSM(ed);
	hed = hes = NULL;  nzd = NULL;
	r = p = dd = 0;
	if (r==0)  if (rs!=rd || cs!=cd)  r = -3;
	if (r==0)  DEBUG(CURDLEV-0,"Testing (unlabeled strong) matrix equivalence of [%s] (%dx%d) to [%s] (%dx%d).\n",
			EMNAME(es)?EMNAME(es):"",ROWSM(es),COLSM(es),EMNAME(ed)?EMNAME(ed):"",ROWSM(ed),COLSM(ed));
	if (r==0)  if ((hhs&&hhd)? !strmag_isinjection(hhs,rs+cs,hhd,rd+cd):0) {
		r = -3;
		DEBUG(CURDLEV-2+(ch<0),"The given labels (hash codes) allow -NO- equivalence %p [%s] to %p [%s] (%dx%d).\n",
				es,EMNAME(es)?EMNAME(es):"",ed,EMNAME(ed)?EMNAME(ed):"",ROWSM(ed),COLSM(ed));
	}
	sm_callmateqiv++;
		/**
		 * We use the given line-codes hhs,hhd only as hints, but it is not
		 * guaranteed that only the same codes are mapped to each other.
		 * These codes are here enhanced by magic values based on small flats.
		 * We check whether the enhanced values allow for an injection.
		**/
	if (r==0) {
		if (2*rs+2*cs<48)  hes = hestack;
		else  hes = MMALLOC((2*rs+2*cs+4)*sizeof(hes[0]));
		hed = hes+rs+cs+2;
		if (fls)  for (i=0; i<rs+cs; i++)  hes[i] = fls[i];
		else  strmag_flatlines(es,FLN_MINR_EQUIV,hes);
		if (fld)  for (i=0; i<rd+cd; i++)  hed[i] = fld[i];
		else  strmag_flatlines(ed,FLN_MINR_EQUIV,hed);
		if (hhs && hhd)  for (i=0; i<rs+cs; i++) {
			hes[i] += 9999*hhs[i];  hed[i] += 9999*hhd[i];
		}
		if (!strmag_isinjection(hes,rs+cs,hed,rd+cd))
			r = -2;
	}
		/**
		 * If all equivalences are required for the output, we use the above
		 * minor-finding functions with their print-outs (but no self-maps!).
		 * The result of struct_minorbybas_() is then returned, plus possible
		 * refering matrices stored in the list *mls if required.
		**/
	if (r>-2)  if (ch>=3 || mls) {
		if (r==0)  DEBUG(CURDLEV-1+(ch<0),"Looking for all (unlabeled strong) equivalences between matrices [%s] (%dx%d) and [%s] (%dx%d).\n",
					EMNAME(es)?EMNAME(es):"",rs,cs,EMNAME(ed)?EMNAME(ed):"",rd,cd);
		if (r>=0)  r = struct_minorbybas_bl(ch,es,ed,NULL,mls,hhs,hhd);
		if (r==0)  {PROGERROR("cannot get return 0 value here");}
	}
		/**
		 * Otherwise, we only ask whether the given matrices are equivalent or not.
		 * So we use various tricks to speed-up the decision process.
		 * We first try to look for an immediate displayed minor.
		 * Then we have to search through all (suitable) bases for dispminors.
		 * If the line codes give a nontrivial distiction between elements,
		 * then we use all bases respecting these codes (ematrix_getbases_hh) --
		 * but these bases must not be recycled because of their special nature.
		 * The list to search is taken from es.
		**/
	if (r==0)  if (struct_dispminor_eqv(es,ed,hes,hed)>=0)
		r = 1;
	if (r!=0)  DEBUG(CURDLEV+(ch<0),"Equivalence pre-decision was made r=%d\n",r);
	else  sm_callmateqiv2++;
	ee = NULL;  bl = bls = bld = NULL;
	if (r==0 && rnd<=0)
		for (i=dd=1; i<rs+cs; i++)	/* (this is not really the number of distinct values among hes[]!) */
			if (hes[i]!=hes[i-1])  dd++;
	if (r==0 && dd>2) {
		ee = ematrix_copy(es);  if (EMNAME(es))  EMSETNAME(ee,EMNAME(es));
		bl = ematrix_getbases_hh(ee,hes,hed);
		DEBUG(CURDLEV+1,"Equivalence testing uses code-selected sublist of %d bases.\n",alist_getlength(bl));
	} else  dd = 0;
		/**
		 * If the codes are not used for bases above, we get all bases for both es,ed,
		 * compare their numbers, and check the numbers of zeros in them.
		 * We then search through those bases having the least-repeating number
		 * of zeros among all.
		 * These standard lists of bases are finally recycled.
		**/
	if (r==0 && !ee) {
		if (rnd>0) {
			bls = struct_minorbases_rn(es,rnd);
			DEBUG(CURDLEV+1,"Equivalence testing uses random list of %d bases.\n",alist_getlength(bls));
		} else {
			bls = struct_minorbases(es);
			if (alist_getlength(bls)>500)  bld = struct_minorbases(ed);
			if (bld)  if (alist_getlength(bls)!=alist_getlength(bld))
				r = -2;
			DEBUG(CURDLEV+1,"Equivalence testing uses full lists of %d and %d bases.\n",alist_getlength(bls),bld?alist_getlength(bld):1);
		}
	}
	if (r==0 && !ee && bld && rnd<=0) {
		if (rs*cs<94)  nza = nzastack;
		else  nza = MMALLOC((rs*cs+8)*sizeof(nza[0]));
		nzd = nza+rs*cs+2;
		for (i=0; i<rs*cs+2; i++)  nza[i] = nzd[i] = 0;
		for (x=bls,y=bld; x?*x:0; x++,y++) {
			p = strmag_numzeros(*x);  q = strmag_numzeros(*y);
			if (p<1 || p>rs*cs+1 || q<1 || q>rs*cs+1)  {PROGERROR("wrong number of zeros in a matrix ???"); p = 1;}
			nza[p]++;  nzd[q]++;
		}
		SDEBUG(CURDLEV+2,"\t\tBases - numbers by zeros:");
		for (i=p=0,j=999999; i<rs*cs+2 && r==0; i++) {
			if (nza[i]!=nzd[i])  r = -2;
			if (nza[i]>0 && nza[i]<=(5*j)/4)  p = i;
			if (nza[i]>0 && nza[i]<j)  j = nza[i];
			if (nza[i]>0)  SDEBUG(CURDLEV+2," %d,",nza[i]);
		}
		SDEBUG(CURDLEV+2,"\n");
	}
	if (r==0 && !ee) {
		for (x=y=bld; x?*x:0; x++)	/* (the value of p is set in the previous block) */
			if (strmag_numzeros(*x)==p) { y = x;  break; }
		ee = ematrix_copy(y?*y:ed);
		p = strmag_numzeros(ee);
		for (x=bls; x?*x:0; x++)
			if (strmag_numzeros(*x)==p)  bl = alist_append(bl,*x);
		if (p>0 && nza) if (alist_getlength(bl)!=nzd[p])  {PROGERROR("wrong length of selected sublist by zeros");}
		DEBUG(CURDLEV+1,"Equivalence testing uses sublist of %d bases by %d zeros.\n",alist_getlength(bl),p-1);
	}
		/**
		 * Finally, we get to comparing dispminors in a (sub)list of bases bl
		 * for the matrix es that was prepared above.
		 * We collect no resulting dispminors, only use the return value r.
		 * After that, we free or recycle the lists prepared above.
		 * Notice that the matrix ee replaces es in the first (coded) case,
		 * while ee replaces ed in the second (full-lists) case.
		 * The second case uses no hash-codes since they provide no distinction.
		**/
	if (r==0)  if (!ee || (bl?!bl[0]:1))  r = -1;
	if (r==0) {
		DEBUG(CURDLEV-0,"Trying (unlabeled strong) equivalence of [%s] to [%s] (%dx%d) for %d displayed bases.\n",
				EMNAME(es)?EMNAME(es):"",EMNAME(ed)?EMNAME(ed):"",ROWSM(ed),COLSM(ed),alist_getlength(bl));
		sm_callmateqivb += alist_getlength(bl);
		printlev--;
		if (dd>0)  r = struct_minorbybas_bl(ch,ee,ed,bl,NULL,hes,hed);
		else  r = struct_minorbybas_bl(ch,es,ee,bl,NULL,NULL,NULL);
		printlev++;
	}
	
#ifndef	FASTPROG
	DEBUG(CURDLEV-2+(r<0)+(ch<0),"Found %s equivalence for matrices [%s] and [%s] (%dx%d).\n",
			r<0?"-NO-":"+YES+",EMNAME(es)?EMNAME(es):"",EMNAME(ed)?EMNAME(ed):"",rd,cd);
	if (ch>=0 && IFRANDDEBUGLESS(222) && rnd<=0) {
		ematrix	*ex;
		if (hhs || fls || alist_getlength(bl)<2)  ex = es;
		else  ex = bl[RANDOM()%alist_getlength(bl)];
		if ((r<0)!=(struct_matequivalence_ext(-1,ed,ex,NULL,hhd,hhs,fld,fls,0)<0))  {PROGERROR("non-symmetric result for exchanged matrices");}
	}
	if (ch>=0 && ch<2 && IFRANDDEBUGLESS(555) && rnd<=0 && !hhs) {
		if ((r<0)!=(struct_minorbybas(-1,ed,es,NULL,0)<0))  {PROGERROR("different result got by (disp)minor testing");}
		DEBUG(CURDLEV+1,"Finished with recursive self-testing of matrix equivalence.\n");
	}
#endif
	if (dd!=0 && bl)  dispose_alist_mats(bl);	/* (disposing basis sublist by codes) */
	else  alist_free(bl);		/* (or, freeing basis sublist from the full list bld) */
	if (bls && rnd==0)  struct_minorbases_recycle(bls);
	else if (bls)  dispose_alist_mats(bls);
	if (bld && bld!=bls && rnd==0)  struct_minorbases_recycle(bld);
	else if (bld)  {PROGERROR("cannot have bld for randomized");}
	if (ee && ee!=es)  dispose_ematrix(ee);
	if (hes && hes!=hestack)  FREE(hes);
	if (nza && nza!=nzastack)  FREE(nza);
	if (rnd>0 && r<-1)  return -10;	/* this randomized NO result is for sure */
	return (r<0? r: ((ch>=3||mls)? r:0));
}


























/******************	Minor testing interface	*******************/
/******************************************************************/


#undef CURDLEV
#define	CURDLEV		6		/*7 (need <=DEBUGLEV+2 to perform paranoic tests) */


/**
 * This function looks in the matrix e for (represented) minors isomorphic to a member
 * of the given list lmin.
 * (In fact, only those unlabeled-equivalent to the given representation of the minor!).
 * Instead of the list lmin, one may give just one matrix emin.
 * If rnd>0, then only a randomized test is performed - for a sublist of bases.
 * So rnd>0 may return NO even if the minor exists, but a YES answer is always sure.
 * 
 * The return value is i>0 if lmin[i-1] is a minor in e (the first one found), or 0 if no minor.
 * No other output is collected or printed.
**/

int	struct_hasminorlist_ext(int ch, ematrix *e, ematrix *emin, ematrix **lmin, int rnd) {
	int	p, ro,co,rm,cm;
	ematrix	**x;
	
	if (!e || !emin==!lmin)  {PROGERROR("the matrix e and exactly one of emin,lmin must be given here"); return 0;}
	ro = ROWSM(e);  co = COLSM(e);
	DEBUG(CURDLEV-1,"Testing a list of %d minors in the matroid %p [%s] (%dx%d).\n",
			lmin?alist_getlength(lmin):1,e,EMNAME(e)?EMNAME(e):"",ROWSM(e),COLSM(e));
	p = -1;
	for (x=lmin?lmin:&emin; (x?*x:0) && (lmin?1:x==&emin) && p<0; x++) {
		rm = ROWSM(*x);  cm = COLSM(*x);
		if (rm>ro || cm>co)
			continue;
		if (rm==ro && cm==co) {
			p = struct_matequivalence_rn(ch,e,*x,rnd);
		} else if (rnd<=0 && ch<3 &&
				ro+co-rm-cm<=STR_MAXMINREM && ro+co<6*(rm+cm)/5) {
			p = struct_minorbyrem(ch,e,*x);
		} else {
			p = struct_minorbybas(ch,e,*x,NULL,rnd);
		}
		if (p>=0)  p = (lmin? x-lmin:0)+1;
	}
	return (p>0? p:0);
}


/**
 * This function looks in the matrix e for all (represented) minors isomorphic to emin.
 * (In fact, only those unlabeled-equivalent to the given representation of the minor!).
 * If af!=0, then the self-maps of emin are factored out from all (disp)minors.
 * 
 * The return value is -1 for no minor found, or 0 or the number of minors found here.
 * The collected minors are stored to the list *mls (which must be given empty or NULL)
 * as refering matrices, together with copies of bases of e that display them (same list!).
 * (See in struct_minorbybas_ext()...)
 * 
 * If ch==3, then the first minor found is printed out, all of them are printed if ch>=4.
 * 
**/

int	struct_allminors_ext(int ch, ematrix *e, ematrix *emin, int af, ematrix ***mls) {
	int	p;
	char	buf[100];
	
	if (!e || !emin)  {PROGERROR("both matrices e, emin must be given here"); return -1;}
	if (mls?*mls:0)  {PROGERROR("the list mls must be given empty (null)"); return -1;}
	if (ch>=4) {
		SOUTPUT("\n");
		OUTPUT("  Looking for all [%s]-minors (%dx%d) in the matroid [%s] (%dx%d).\n",
				EMNAME(emin)?EMNAME(emin):"",ROWSM(emin),COLSM(emin),EMNAME(e)?EMNAME(e):"",ROWSM(e),COLSM(e));
		EMATOUTPUT(e,printoutpref);
		snprintf(buf,90,"%s >  ",printoutpref);  EMATOUTPUT(emin,buf);
	}
	p = struct_minorbybas_all(ch,e,emin,(af?(void**)1:NULL),mls);
	if (ch>=4)  SOUTPUT("\t============\n\n");
	return p;
}



/**
 * This function determines (unlabeled strong) matrix equivalence between the
 * given matrices es and ed or against the whole list ld.
 * Read more at struct_matequivalence_ext() about equivalence testing.
 * 
 * The return value is i>0 if ld[i-1] is equivalent to es (the first one found), or 0 if no one.
 * If ch>=3, then the equivalences are printed as they are found; all of them for ch>=4.
 * There is no randomized version available for now.
**/

int	struct_isequivlist_ext(int ch, ematrix *es, ematrix *ed, ematrix **ld) {
	long		*fln, flstack[50];
	ematrix		**x;
	int		p,r;
	
	if (!es || !ed==!ld)  {PROGERROR("the matrix es and exactly one of ed,ld must be given here"); return 0;}
	DEBUG(CURDLEV-1,"Testing equivalence of the matroid %p [%.9s] (%dx%d) against a list of %d others.\n",
			es,EMNAME(es)?EMNAME(es):"",ROWSM(es),COLSM(es),ld?alist_getlength(ld):1);
	if (ch>=4 && ld) { SOUTPUT("\n"); OUTPUT("  Looking for (unlab strong) equivalences betw matr [%s.14] (%dx%d) and list of %d.\n",
					EMNAME(es)?EMNAME(es):"",ROWSM(es),COLSM(es),alist_getlength(ld)); }
	
	if (ROWSM(es)+COLSM(es)<46 && (ld?alist_getlength(ld):0)>2)
		strmag_flatlines(es,FLN_MINR_EQUIV,fln=flstack);
	else  fln = NULL;
	p = r = -1;
	for (x=ld?ld:&ed; (x?*x:0) && (ld?1:x==&ed); x++) {
		if (ROWSM(es)!=ROWSM(*x) || COLSM(es)!=COLSM(*x))
			continue;
		if (ch>=4) {
			SOUTPUT("\n");  OUTPUT("  Looking for all (unlabeled strong) equivalences betw matrices [%.14s] (%dx%d) and [%.14s].\n",
					EMNAME(es)?EMNAME(es):"",ROWSM(es),COLSM(es),EMNAME(*x)?EMNAME(*x):"");
			EMATOUTPUT(es,printoutpref);  EMATOUTPUT(*x,printoutpref);
		}
		p = struct_matequivalence_mbr(ch,es,*x,fln,NULL);
		if (p>=0 && r<0)  r = (ld? x-ld:0)+1;
		if (ch>=4)  SOUTPUT("\t============\n\n");
		if (ch<=3 && r>=0)  break;
	}
	if (ch>=2)  OUTPUT(" Found equivalence result %s (%d) for matrices [%.14s] (%dx%d) and [%.14s].\n",
			r>0?"+YES+":"-NO-",r,EMNAME(es)?EMNAME(es):"",ROWSM(es),COLSM(es),ed?(EMNAME(ed)?EMNAME(ed):""):"-list-");
	return (r>0? r:0);
}

/**
 * This function is similar to the previous one --
 * it determines (unlabeled strong) matrix equivalence between the pairs of
 * given matrices in the list le.
 * Read more at struct_matequivalence_ext() about equivalence testing.
 * 
 * The return data are stored through the given array pointer *out
 * (data is allocated here, and must be freed later!).
 * The value of (*out)[i]=j means that #j is the first matrix in le equivalent to #i.
 * For higher values of ch, more information is printed out.
 * (See struct_matequivalence_() for ch-1...)
 * There is no randomized version available for now.
**/

int	struct_equivlistself_ext(int ch, ematrix **le, int **out) {
	long		**fln=NULL, *fdt=NULL;
	ematrix		**x,**y;
	int		i,j,r, ln, ret=0,uq=0, *ra;
	
	if (!le)  {PROGERROR("the matrix list le must be given here"); return 0;}
	ln = alist_getlength(le);
	DEBUG(CURDLEV-1,"Testing equivalence in the matr list {%.8s %.8s %.8s }%d.\n",
			le[0]?EMNAME(le[0]):"",(le[0]&&le[1])?EMNAME(le[1]):"",(le[0]&&le[1]&&le[2])?EMNAME(le[2]):"",ln);
	if (ch>=2-(ln>6))  OUTPUT("  Looking for (unlabeled strong) equivalence in matr list {%.8s %.8s %.8s}%d.\n",
			le[0]?EMNAME(le[0]):"",(le[0]&&le[1])?EMNAME(le[1]):"",(le[0]&&le[1]&&le[2])?EMNAME(le[2]):"",ln);
			/**
			 * The flatline codes for the matrix lines in le are pre-computed
			 * here, and later passed to struct_matequivalence_() for faster computation.
			 * We compute higher-level flatline codes for the matrix lines in le
			 * when the list is reasonably long.
			 * The clusters of equivalent matrices are then stored in ra[].
			**/
	for (i=0, x=le; x?*x:0; x++)  i += ROWSM(*x)+COLSM(*x);
	fdt = MMALLOC((i+2)*sizeof(fdt[0]));
	fln = MMALLOC((ln+2)*sizeof(fln[0]));
	for (i=0, x=le; x?*x:0; x++) {
		fln[x-le] = fdt+i;
		i += ROWSM(*x)+COLSM(*x);
		if (ln>33)  strmag_flatlines_more(*x,2,fln[x-le]);
		else  strmag_flatlines(*x,2,fln[x-le]);
	}
	ra = MMALLOC((ln+2)*sizeof(ra[0]));
	for (i=0; i<ln; i++)  ra[i] = -1;
			/**
			 * We test all pairs of matrices from the given list le for strong
			 * matrix equivalence using struct_matequivalence_mbr().
			 * We record the results in the array ra[]; where ra[i] keeps the
			 * smallest index of a matrix equivalent to #i.
			**/
	for (x=le, i=0; x?*x:0; x++, i++) {
		if (ch>=4 && x[1]) {
			SOUTPUT("\n");  OUTPUT(" Looking for all (further) equivalences of the matroid #%d [%.18s] %dx%d in the list.\n",
						i+1,EMNAME(*x)?EMNAME(*x):"",ROWSM(*x),COLSM(*x));
			if (ch>=5) { EMATOUTPUT(*x,printoutpref);  SOUTPUT("\n"); }
		}
		if (ra[i]<0 && ch>=2)  OUTPUT("Matroid #%d [%.18s] %dx%d is (equiv) unique in the list so far.\n",
						i+1,EMNAME(*x)?EMNAME(*x):"",ROWSM(*x),COLSM(*x));
		if (ra[i]<0) { ra[i] = i;  uq++; }
	  for (y=x+1, j=i+1; y?*y:0; y++, j++) {
		if (ra[j]>=0 || (ra[i]>=0 && ra[i]<i))  continue;
		if (ROWSM(*x)!=ROWSM(*y) || COLSM(*x)!=COLSM(*y))  continue;
		
		r = struct_matequivalence_mbr((ch>0?ch-1:0),*x,*y,fln[i],fln[j]);
		if (r>=0) { ra[j] = i;  ret++; }
	}}
	
	if (ch>=1)  OUTPUT("Total %d mats from {%.7s %.7s %.7s}%d are unique to (strong) matrix equivalence.\n",
			uq,le[0]?EMNAME(le[0]):"",(le[0]&&le[1])?EMNAME(le[1]):"",(le[0]&&le[1]&&le[2])?EMNAME(le[2]):"",alist_getlength(le));
	if (fdt)  FREE(fdt);  if (fln)  FREE(fln);
	if (out)  *out = ra;
	else if (ra)  FREE(ra);
	return ret;	/* (isn't it better to return uq ???) */
}




/**
 * This function tests whether the matrix e has an element that can be both contracted and
 * deleted keeping some minors from the list tghl[] (NULL-terminated) - "removable" element.
 * (If not, then the matroid of e is a "tight major" of the list tghl[].)
 * As a "minor" we mean a displayed subrepresentation in an equivalent matrix to em
 *  - see struct_dispminor_..().
 * All bases (unless ch==-7) of e are tested against all minors in tghl[] in this function.
 * 
 * If ch>0+, then some more information are printed, and paranoic (debug) tests performed if ch>=0.
 * If ch==-7, then only a randomized test is performed - for only a sublist of bases.
 * The function returns 1 if a removable (del+contr) element is found, and 0 if not.
 * 
 * It is possible to request printing the removable element by struct_hasdelcontr_printit().
**/

int	struct_hasdelcontrlist_ext(int ch, ematrix *e, ematrix *tgh, ematrix **tghl, int rnd) {
	ematrix		*em,*rr, **remm,**delm,**conm, **basl,**minl, **x,**y,**z,**w;
	emlinemap	***tgmap, **a,***ax;
	int		i,j,k, ld,p,r,ret, minr,minc;
	
	if (!e || !tgh==!tghl)  {PROGERROR("the matrix e and exactly one of tgh,tghl must be given here"); return 0;}
	DEBUG(CURDLEV-1,"Testing %d \"tight-major\" %sminors in the matroid %p [%.18s] (%dx%d).\n",
			tghl?alist_getlength(tghl):1,rnd>0?"randomized ":"",e,EMNAME(e)?EMNAME(e):"",ROWSM(e),COLSM(e));
	sm_calldispmin = sm_calldispmin2 = sm_calldispmin3 = 0;
	if (ch>2)  OUTPUT("  Testing %d \"tight-major\" %sminors in the matroid [%.18s] (%dx%d).\n",
			tghl?alist_getlength(tghl):1,rnd>0?"randomized ":"",EMNAME(e)?EMNAME(e):"",ROWSM(e),COLSM(e));
		/**
		 * We prepare local data for the computation.
		 * We copy the matrix e to em with resetted line labels, and we also copy
		 * the list tghl, possibly appending the tgh matrix.
		 * We also pre-compute the minor self-maps to the list tgmap.
		 * We use the arrays delm[],conm[] (indexed by line id labels!) to store
		 * what elements may be deleted/contracted -- keeping the displayed minors.
		**/
	em = ematrix_copy_notr(e);
	if (EMNAME(e))  EMSETNAME(em,EMNAME(e));
	ld = ematrix_resetid(em);
	remm = MMALLOC((4*ld+12)*sizeof(remm[0]));
	for (i=0; i<4*ld+12; i++)  remm[i] = NULL;	/* (free the matrix pointers below!) */
	delm = remm+ld+2;  conm = remm+3*ld+6;
	tgmap = NULL;
	tghl = (tghl? alist_copy(tghl):NULL);
	if (tgh)  tghl = alist_append(tghl,tgh);
	minr = minc = INT_MAX;
	for (x=tghl; x?*x:0; x++) {
		a = NULL;  struct_matrixselfmaps(*x,(void***)&a);
		tgmap = alist_append(tgmap,a);
		if (ROWSM(*x)<minr)  minr = ROWSM(*x);
		if (COLSM(*x)<minc)  minc = COLSM(*x);
	}
	basl = (rnd>0? struct_minorbases_rn(em,rnd): struct_minorbases(em));
		/**
		 * Here we search through all displayed minors from tghl among all bases
		 * (possibly randomized) of the given matrix em.
		 * For each dispminor as a refering matrix, we note all removable remaining
		 * lines if they are not marked yet, with a copy of the minor.
		 * When we get a del+contr element, we are done.
		**/
	ret = r = 0;  j = 1;
	if (minr<ROWSM(em) && minc<COLSM(em))
	  for (x=basl; (x?*x:0) && !r; x++)  for (z=tghl; (z?*z:0) && !r; z++) {
		minl = NULL;
		p = struct_dispminor(*x,*z,(void**)(tgmap[z-tghl]),&minl);
		if (p>=0)
		  for (y=minl; (y?*y:0) && !r; y++) {
			rr = ematrix_refextract_xall(*x,*y);
			if (rr)
			  for (i=0; i<ROWSM(rr)+COLSM(rr) && !r; i++) {
				j = (i<ROWSM(rr)? ROWSID(rr,i): COLSID(rr,i-ROWSM(rr)));
				w = (i<ROWSM(rr)? &conm[j]: &delm[j]);
				if (*w==NULL) {
					*w = ematrix_refer_all(*y);
					EMDATA(*w) = strmap_copy(EMDATA(*y));
					if (EMNAME(*z))  EMSETNAME(*w,EMNAME(*z));
				}
				if (ch<=2) if (i<ROWSM(rr)? delm[j]:conm[j])  r = 1;
			}
			if (rr)  dispose_ematrix(rr);
		}
		if (minl)  dispose_alist_mats(minl);
	}
		/**
		 * Finally, we get the results out, and return the possible del+contr element
		 * in the return value -- 1+row or 1+col+ROWSM with respect to given e.
		 * Then all internal data are disposed.
		**/
	if (ch>=2)  for (i=-ld; i<=ld; i++) if (delm[i] && conm[i]) {
		j = i;  r = 1;
		OUTPUT("  A contr+del element %s %d(%d) was found in [%.16s],  contr->[%s], del->[%s].\n",
				j>0?"row":"col",j<0?-j-1:j-1,j<0?COLSID(e,-j-1):ROWSID(e,j-1),EMNAME(e)?EMNAME(e):"",EMNAME(conm[i]),EMNAME(delm[i]));
		//********  possibly add printing of the line maps from EMDATA(delm[]/conm[])
	}
	k = (j<0? ROWSM(e)-j-1:j-1);	/* (j was set above to common removal id) */
	ret = r>0? k+1:0;
#ifndef FASTPROG
	DEBUG(CURDLEV-2+(ret<=0)+(ch<0),"Found %p [%s] (%dx%d) %s a \"tight-major\" of the list %d.\n",
			e,EMNAME(e)?EMNAME(e):"",ROWSM(e),COLSM(e),(!ret?"+IS+":"is -NOT-"),alist_getlength(tghl));
	SDEBUG(CURDLEV-1+(ch<0),"\t\t(Tight-major test used %d %sbases, and called %d (%d,%d) dispminors (deep).)\n",
			alist_getlength(basl),rnd>0?"random ":"",sm_calldispmin,sm_calldispmin2,sm_calldispmin3);
	if ((k<ROWSM(em)?ROWSID(em,k):COLSID(em,k-ROWSM(em)))!=j)  {PROGERROR("wrong line picked up ?!? %d %d",j,k);}
	if (ch>=0 && rnd<=0 && IFRANDDEBUGLESS(111)) {
		rr = basl[alist_getlength(basl)/2];
		if ((!ret)!=(!struct_hasdelcontrlist_ext(-1,rr,NULL,tghl,0)))  {PROGERROR("different result for an equivalent matrix");}
	}
#endif
	if (remm)  for (i=0; i<4*ld+12; i++)
		if (remm[i])  dispose_ematrix(remm[i]);
	if (remm)  FREE(remm);
	if (rnd!=0)  dispose_alist_mats(basl);
	else  struct_minorbases_recycle(basl);
	for (ax=tgmap; ax?*ax:0; ax++)
		struct_matrixselfmaps_recycle(*ax);
	if (tgmap)  alist_free(tgmap);
	if (tghl)  alist_free(tghl);	/* (a private copy of the given list is created above!) */
	dispose_ematrix(em);
	return ret>0? ret:0;
}




































