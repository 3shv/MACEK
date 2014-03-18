
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
 * 
 * This piece of code implements generating one-line extensions of the given matrix.
 * (A job that may seem simple, but it is not so simple in partial fields with infinitely
 * many elements -- then the choice for extension elements gets quite tricky...)
 * All possible extensions are generated; in the partial case even many
 * non-properly represented ones.
 * One may generate only a random sublist, or require representability testing for all
 * matrices.
 * 
 * See gener.c for additional (canonical) filtering applied to the generated extensions.
 * 
**/


#include "macek.h"
#include "gen.h"













/******************	Generating matrix extensions	*******************/
/**************************************************************************/


#define	CURDLEV		7		/*7 (use <=DEBUGLEV+2 to perform extra checks) */



/**
 * This function generates column/row (by tr=1/0) extensions to the matrix ex.
 * All distinct extension lines are generated here, normalized to unit-leading form
 * (the first unit entry for each component of a disconnected matrix).
 * When extending over a partial field, then also some pfield-checks of the matrix
 * are used to restrict the number of generated choices.
 * (But the extension is not guarranteed to be properly pfield-represented unless tst>=1!)
 * 
 * The generated extensions are collected in the list *els (if given).
 * (Remember to free the list and the matrices later.)
 * The function returns the number of generated extensions, or -1 on an error.
 * If ch>1 (>2), then the generated extensions are all checked and printed out.
 * If ez is given, then the extended line of ex must have the same zero pattern as
 * that line in ez (this is for eeq<0).
 * Moreover, if ez and eeq>0 are given, then the first eeq entries of the extension
 * line must match the first eeq entries (in row 0) of ez exactly.
 * 
 * If rnd>0, then only a random sublist of extensions is generated, about rnd
 * choices per each entry.
 * If tst==0, then all generated extensions are tested in a partial field randomly,
 * and if tst==1, then they are fully tested (all aplies only to partial fields).
**/

int	gener_matextens_ext(int ch, ematrix *ex, int tr, int rnd, int tst, 
						ematrix *ez, int eeq, ematrix ***els) {
	long		**ppat,**mpat;
	char		pprf[50];
	int		i,j,k,l, r,lc, ro,co,cco, rnf;
	int		*rcon,*rcpi,**radj, *selsg,*selex;
	ematrix		*ee, *edet, **rre, *rr;
	
	if (els?*els:0)  {PROGERROR("the list els must be given empty (null)"); return -1;}
	if (ch>1 && tst<1)  tst = 1;	/* (must test all extensions when printing output) */
	if (ch>2)  sprintf(pprf,"%.30s \t",printoutpref);
	if (ch>1)  OUTPUT("Going to generate %s extensions of the matrix %p[%s] in %s.\n",
				tr?"col":"row",ex,EMNAME(ex)?EMNAME(ex):"",pfield_curname());
	if (ch>1 && ez && eeq>0)  SOUTPUT("\t (Starting with extension entries as in \'%s\'...)\n",gener_extprintline(0,ez,0));
	if (ch>1)  EMATOUTPUT(ex,printoutpref);
#ifndef FASTPROG
	DEBUG(CURDLEV-1+(ch<0),"Going to generate %s extensions of the matrix %p[%s] %dx%d in %s...\n",
			tr?"col":"row",ex,EMNAME(ex)?EMNAME(ex):"",ROWSM(ex),COLSM(ex),pfield_curname());
	if (rnd>0)  DEBUG(CURDLEV-1+(ch<0),"  - generating only random (%d) subset of extensions\n",rnd);
	if (tst>=0 && pfispartial)  DEBUG(CURDLEV-1+(ch<0),"  - %s testing extensions for proper representation\n",tst>0?"full":"");
	EMATDEBUGS(CURDLEV+1,ex,"\t  ex\t");
	if (ez && eeq<0) {
		DEBUG(CURDLEV-1+(ch<0),"  - using zero pattern as in (not the right pfield!) %s\n",gener_extprintline(0,ez,tr));
		if (ROWSM(ez)<ROWSM(ex)+(!tr) || COLSM(ez)<COLSM(ex)+(!!tr))  {PROGERROR("The matrix ez is too small %dx%d",ROWSM(ez),COLSM(ez));}
	}
	if (ez && eeq>0) {
		DEBUG(CURDLEV-1+(ch<0),"  - using extension prefix as in %s\n",gener_extprintline(0,ez,0));
		if (ROWSM(ez)<1 || COLSM(ez)<eeq)  {PROGERROR("The matrix ez is too small %dx%d",ROWSM(ez),COLSM(ez)); eeq = 0;}
	}
	if (ez)  EMATDEBUGS(CURDLEV+1,ez,"\t  ez\t");
#endif
	
	if (!tr)  ematrix_transpose(ex);
	if (!tr && ez && eeq<0)  ematrix_transpose(ez);
	r = 0;
	j = (rnd>0? rnd: pfield_numfundamentals())+1;
	j = (pfispartial? j/(pfnumexp+1): j);
	for (i=0,k=1; i<ROWSM(ex)-(eeq>0?eeq:1); i++)
		if (k<GEN_MAXMATEXTENS)  k *= (j<2? 2:j);
	if (k>GEN_MAXMATEXTENS)  {PROGERROR("there are too many expected extensions to generate"); r = -11;}
	rnf = (rnd>0? (pfield_numfundamentals()+1)/rnd: 0);
	
		/**
		 * We prepare here working matrices ee (for the extended column) and
		 * edet (for choosing pfield value according to previous nonzero subdeterminants).
		 * Then we get the "connected" order rcon[] of rows of ex, component
		 * indexing in rcpi[], and the adjacency matrix of rows in radj[][].
		 * Here rcpi[] is indexed same as rcon[] -- rcpi[i] is the component of row rcon[i],
		 * but radj[][] is indexed in absolute indices.
		 * (We then fill in the entries of ee in order rcon[0],rcon[1],...)
		**/
	ro = ROWSM(ex);  cco = ro*ro;  co = COLSM(ex);
	ee = ematrix_copy_to(ex,ro,co+1);
	ee = ematrix_append_col(ee);
	edet = new_ematrix(ro,cco+1,ro,cco+1);
	rre = NULL;
	radj = pfispartial? malloc_twodim(sizeof(radj[0][0]),ro,ro):NULL;
	rcon = MMALLOC((2*ro+6)*sizeof(rcon[0]));
	rcpi = rcon+ro+2;
	j = struct_connorder_radj(ex,rcon,rcpi,radj);
	if (j<0)  DEBUG(CURDLEV-0+(ch<0),"  - generating extensions for a disconnected matrix(!!!)\n");
	DEBUG(CURDLEV-1+(ch<0),"  - using \"connecting\" order %d,%d,%d,%d,... (%d)\n",rcon[0],rcon[1],rcon[2],rcon[3],ro);
	ppat = mpat = NULL;
		/**
		 * For a partial field, we prepare induced paths for all pairs of rows -
		 * the array ppat[][] contains the internal vertices of these paths,
		 * and mpat[][] the connecting columns for these paths.
		 * The rows/columns i are encoded by the bits (1l<<i).
		 * Additionally, we prepare refering matrices rre[] to ee which will
		 * be used for representability pre-checks in a partial field.
		**/
	if (pfispartial && r>=0) {
		ppat = malloc_twodim(sizeof(ppat[0][0]),ro,ro);
		mpat = malloc_twodim(sizeof(mpat[0][0]),ro,ro);
		j = struct_connpaths(ro,radj,rcon,ppat,mpat);
		if (j<0)  r = -2;
		else  j = gener_matextens_cdet(ex,rcon,ppat,mpat,edet);
		if (j<0)  r = -2;
		rre = MMALLOC(ro*sizeof(rre[0]));
		rr = ematrix_refer(ee,-1,-1,0,COLSM(ee));
		for (i=0; i<ro; i++) {
			ematrix_refadd_row(rr,rcon[i]);
			if (i<=0 || i>=ro-2 || GEN_USEEXTPREPF<=0)  rre[i] = NULL;
			else  rre[i] = ematrix_refer_all(rr);
		}
		dispose_ematrix(rr);
	}
	
		/**
		 * We construct the k-th entries of the added line for k=0,1,2,... here:
		 * We use functions pfield_get[pa]value() to guess possibilities
		 * for the entry inside gener_matextens_cycl();
		 * the arrays selsg[],selex[] keep our choices (-1 means the initial choice).
		 * Moreover, we check that the first nonzero entry in each conn. component
		 * of rows of ee is always 1, other nonzeros are then arbitrary.
		 * (This is in real index order of ee, not "through" rcon[].)
		 * 
		 * The matrix edet is used here over and over as the working storage
		 *  - its rows are indexed straight 012.., not via rcon[]!
		 * The new entries are stored into ee later, in order given by rcon[].
		 * Possible representability pre-check in a partial field is carried out
		 * later, and then we continue in the backtracking search.
		**/
	selsg = MMALLOC((2*ro+4)*sizeof(selsg[0]));
	selex = selsg+ro+2;
	k = 0;
	selsg[0] = selex[0] = -1;
	while (k>=0 && r>=0) {
				/* whether this is the first nonzero in a comp. - to be set 1 */
				/* (we set the first one in order rcon[], but then we must rescale the line) */
		for (j=k-1; j>=0 && SIGNM(edet,j,cco)==0; j--) ;
		lc = (j<0?1: rcpi[k]>rcpi[j]);
				/* the actual selection of the next entry happens in this call: */
		l = gener_matextens_cycl(k,ppat,edet,selsg,selex,lc,rnf);
		if (l<=0) {	/* here cycling (l=0) or going back in backtracking (l<0) */
			k += l;  continue;
		}
			/**
			 * These are some supplementary or required pre-checks applied to
			 * the generated extensions along the generating process.
			 * The zero-entry and representability checks do not depend on scaling,
			 * and so they are applied directly on edet.
			 * For the prefix check we first need to scale the extended column to ee,
			 * and to apply this test only if the scaling was successfull.
			**/
		if (ez && eeq<0)	/* zero-value check if ez is given */
			if ((SIGNM(edet,k,cco)==0)!=(SIGNM(ez,rcon[k],co)==0)) {
				DEBUG(CURDLEV+1+(ch<0),"  zero-value (ez) check rejected at k=%d(%d), co=%d\n",k,rcon[k],co);
				continue;
		}
		if (rre?rre[k]:0) {	/* representability pre-check for pfields (see rre[] above) */
			for (j=0; j<=k; j++)  COPYEXSIGM(ee,rcon[j],co,edet,j,cco);
			if (!pfield_isequal(EXPM(edet,k,cco),SIGNM(edet,k,cco),EXPM(rre[k],k,co),SIGNM(rre[k],k,co)))
				{PROGERROR("Wrong rre[%d] matrix for repr pre-check!",k);}
			if (ematrix_inpfield_randlast(rre[k],1,1)<0) {
				DEBUG(CURDLEV+1+(ch<0),"  X repr pre-check rejected at k=%d.\n",k);
				continue;
			}
		}
		if (ez && eeq>0) {	/* exact-prefix check if ez is given (must copy to ee and scale first) */
			j = gener_matextens_cpscal(k+1,edet,cco,rcon,rcpi,ee,ro,co);
			if (j>rcpi[k])  {PROGERROR("cannot have larger component index here %d>%d",j,rcpi[k]);}
			if (j>=rcpi[k])  for (l=0; l<=k; l++) if (rcon[l]<eeq)
				if (!pfield_isequal(EXPM(ee,rcon[l],co),SIGNM(ee,rcon[l],co),
						EXPM(ez,0,rcon[l]),SIGNM(ez,0,rcon[l])))  break;
			if (j>=rcpi[k] && l<=k) {
				DEBUG(CURDLEV+0+(ch<0),"  same-prefix (ez) check rejected at k=%d,l=%d(%d); %s!=ez[%d]=%s.\n",
						k,l,rcon[l],pfield_printvalue(16,EXPM(ee,rcon[l],co),SIGNM(ee,rcon[l],co)),rcon[l],pfield_printvalue(16,EXPM(ez,0,rcon[l]),SIGNM(ez,0,rcon[l])));
				continue;
			}
		}
		if (k<ro-1) {		/* going forward in backtracking - new cycle at level k+1 */
			++k;  selsg[k] = selex[k] = -1;
			continue;
		}
		DEBUG(CURDLEV+3+(ch<0)," \t extension #%d in edet:  %s\n",r+1,gener_extprintline(-1,edet,1));
			/**
			 * After we choose all entries of the extended column (stored in edet so far),
			 * we copy the ext column to ee and scale it so that the first nonzeroes
			 * are equal 1 for each component of rows.
			 * The unit scaling is rather complicated here since we want the
			 * really first unit, not the first one in order of rcon[].
			 * (That is why we could not set the first nonzero in advance!!!)
			 * Then we possibly apply the final representability test,
			 * and output or store the extension.
			 * The return value r counts the number of extensions found here.
			**/
		j = gener_matextens_cpscal(ro,edet,cco,rcon,rcpi,ee,ro,co);
		if (j<rcpi[ro-1])  {PROGERROR("Why the extension line is not correctly scaled???");}
		
		if (tst>=0 && pfispartial) {
			j = (tst>0? ematrix_inpfield_last(ee,0,1):ematrix_inpfield_randlast(ee,0,1));
			if (j<0)  continue;
		}
		DEBUG(CURDLEV+1+(ch<0)," > extension #%d found:  %s\n",r+1,gener_extprintline(-1,ee,1));
		r++;
		if (ch>=2)  OUTPUT(" > Extension #%d to [%s] found:  %s\n",r,EMNAME(ex)?EMNAME(ex):"",gener_extprintline(1,ee,1));
		rr = ematrix_copy(ee);
		if (!tr)  ematrix_transpose(rr);
		if (ch>2)  EMATOUTPUTS(rr,pprf);
		if (els)  *els = alist_append(*els,rr);
		else  dispose_ematrix(rr);
	}
	
		/**
		 * Finally, we perform some late tests for correctness, and clean
		 * all data and matrices used in the generating process above.
		 * However, we must first dispose all rre[] matrices since they refer to
		 * ee which is going to be modified here.
		**/
	if (rre) {
		for (i=0; i<ro; i++)  if (rre[i])  dispose_ematrix(rre[i]);
		FREE(rre);
	}
	ematrix_remove_col(ee,co);
	if (!ematrix_isequal(ee,ex))  {PROGERROR("The matrix ee was damaged here!");}
	if (!tr && ez && eeq<0)  ematrix_transpose(ez);
	if (!tr)  ematrix_transpose(ex);
	if (ch>0)  OUTPUT("Generated %s %d %s-extensions of the matrix %p [%s] (%dx%d).\n",
			rnd>0?"random":"total",r,pfield_curname(),ex,EMNAME(ex)?EMNAME(ex):"",ROWSM(ex),COLSM(ex));
	if (ch>2)  OUTPUT("\n");
#ifndef FASTPROG
	DEBUG(CURDLEV-1+(ch<0),"... All done, generated %d %s-extensions of %p [%s] (%dx%d).\n",
			r,pfield_curname(),ex,EMNAME(ex)?EMNAME(ex):"",ROWSM(ex),COLSM(ex));
	if (ch>=0 && r>=0 && els) {
		for (i=0; i<r; i++) for (j=0; j<r; j++) {	/* (random equal pairs among extensions?) */
			for (l=0; l<ro; l++)  if (!(tr?
					pfield_isequal(EXPM((*els)[i],l,co),SIGNM((*els)[i],l,co),EXPM((*els)[j],l,co),SIGNM((*els)[j],l,co)):
					pfield_isequal(EXPM((*els)[i],co,l),SIGNM((*els)[i],co,l),EXPM((*els)[j],co,l),SIGNM((*els)[j],co,l)) ))
				break;
			if (i!=j && l>=ro)  {PROGERROR("Cannot produce equal (after scaling) extensions! #%d==#%d",i,j); break;}
			if (r>25)  j += RANDOM()%(r/10);
			if (r>25 && j>=r-1)  i += RANDOM()%(r/10);
		}
	}
	if (ch>=0 && r>=0 && els && ez && eeq>0) {
		for (i=0; i<r; i++) {			/* (random prefix test) */
			for (l=0; l<eeq; l++)  if (!pfield_isequal(EXPM((*els)[i],tr?l:co,tr?co:l),SIGNM((*els)[i],tr?l:co,tr?co:l),EXPM(ez,0,l),SIGNM(ez,0,l)))
				break;
			if (l<eeq)  {PROGERROR("The generated extension #%d has wrong prefix!",i); break;}
			if (r>35)  i += RANDOM()%(r/20);
		}
	}
	if (ch>=0 && r>=0 && pfispartial && rnd<=0 && tst>=1 && !ez) {
		rr = ematrix_copy(ex);  ematrix_swaprows(rr,0,ROWSM(rr)-1);
		if (ROWSM(rr)>3)  ematrix_swaprows(rr,1+RANDOM()%3,ROWSM(rr)-2);
		j = gener_matextens_ext(-1,rr,tr,0,1,NULL,0,NULL);
		dispose_ematrix(rr);
		if (j!=r)  {PROGERROR("Wrong number of extensions in a recursive test! %d!=%d",r,j);}
	}
			//************ check recursively extensions with prefix... ???
	
	if (r>=0 && !pfispartial && rnd<=0 && rcpi[ro-1]==rcpi[0]) {
		j = pfield_numfundamentals()+1;
		k = l = 0;	/* (number of nonzeros in the supplementary matrix ez, if given) */
		if (ez && eeq<0)  for (i=0; i<ro; i++)
			if ((tr?SIGNM(ez,i,co):SIGNM(ez,co,i))!=0)  l++;
		if (ez && eeq>0)  for (i=0; i<eeq; i++)
			if (SIGNM(ez,0,i)!=0)  l++;
		if (ez && eeq<0) {	/* extensions with given zero pattern */
			for (i=0,k=1; i<l-1; i++)  k *= j-1;
		} else if (l==0) {		/* normal extensions, possibly with given few entries in ez */
			for (i=0,k=1; i<ro-(eeq>0?eeq:0); i++)  k *= j;
			k = (l>0? k/(j-1): (k-1)/(j-1)+1);
			//************ check also the number of extensions with non-zero prefix...
		}
		if (k>0 && r!=k)  {PROGERROR("Wrong total numer of extensions ?!? %d!=%d\n",r,k);}
	}
#endif
	if (r<0)  {PROGERROR("Generating matrix extensions has failed! %d",r);}
	if (ee)  dispose_ematrix(ee);		/* (used in generating extensions) */
	if (edet)  dispose_ematrix(edet);
	if (selsg)  FREE(selsg);
	if (ppat)  FREE(ppat);		/* (used to store "connected order" information) */
	if (mpat)  FREE(mpat);
	if (rcon)  FREE(rcon);
	if (radj)  FREE(radj);
	return r;
}



/**
 * Here we compute a "connecting subdeterminants" for generating extensions of the matrix ex.
 * This function is usable only for a real partial field.
 * We store the connecting subdeterminants in the matrix edet which was allocated before.
 * The rows of edet correspond to the rows of ex, but in the connected order given in rcon[].
 * Each column of edet corresponds to a pair of rows of ex.
 * 
 * The practical meaning of edet and ppat,mpat is the following:
 * If r1,r2 are two rows of edet, both entries l1,l2 of rows rcon[r1],rcon[r2] in the last column
 * in ex are nonzero, and other entries in the last column of ex corresponding to 1's in the
 * bitfield ppat[r1][r2] are zero, then the entries d1, d2 in
 * the column (r1,r2) of edet at rows r1, r2 are nonzero and d1*l2-l1*d2 is equal to
 * a (+-) subdeterminat of the extended matrix ee.
 * (This allows to guess possible values of l2 knowing l1 and d1,d2 in a partial field.)
 * 
 * The return value is -1 for an error here, or 1 after normal computation.
**/

int	gener_matextens_cdet(ematrix *ex, int rcon[], long **ppat, long **mpat, ematrix *edet) {
	int		i,j, z,m,n, ro,co;
	long		zp,zm;
	ematrix		*re;
	exp_t		xx;
	sign_t		gg;
	
	ro = ROWSM(ex);  co = COLSM(ex);
	if (ro>=NUM_LONG_BITS-2 || co>=NUM_LONG_BITS-2)  {PROGERROR("The matrix is too large to fit into long bits!"); return -1;}
	DEBUG(CURDLEV+0," Computing \"connecting\" subdeterminants for guessing pfield values...\n");
		/**
		 * For each pair of rows i,j, we compute the "connecting" subdeterminants
		 * of ex which have columns given by the bits of mpat[i][j], and rows by the
		 * bits of ppat[i][j] plus the row j or i, resp.
		 * These subdeterminants are then used later to guess pfield entries in ex.
		 * The tricky part here is with signs of these subdeterminants since the rows
		 * in rcon[] are swapped -- so we compute the subdeterminants like if the
		 * end-rows of the path rcon[i],rcon[j] were after all internal rows.
		 * (Only relative sign is significat for us, so this is a correct approach.)
		**/
	for (i=0; i<ro; i++)  for (j=i+1; j<ro; j++) {
		zp = ppat[i][j];  zm = mpat[i][j];
		if (zp<0)  continue;
		for (z=0,m=i; z<2; z++,m=j) {
			re = ematrix_refer_empty(ex);
			for (n=0; n<ro; n++)  if ((zp&(1l<<n))!=0)
				ematrix_refadd_row(re,rcon[n]);
			ematrix_refadd_row(re,rcon[m]);	
			for (n=0; n<co; n++)  if ((zm&(1l<<n))!=0)
				ematrix_refadd_col(re,n);
			if (ROWSM(re)!=COLSM(re))  {PROGERROR("The subdeterminant must be square for an induced path!"); return -1;}
			ematrix_determinant(re,&xx,&gg);
			if (gg==0)  {PROGERROR("The subdeterminant cannot be zero for an induced path!"); return -1;}
			SETEXSIGM(edet,m,i+j*ro,xx,gg);
			SETEXSIGM(edet,m,j+i*ro,xx,gg);
			dispose_ematrix(re);
		}
	}
#ifndef FASTPROG
	DEBUG(CURDLEV+1,"  - relative 2x2determinant connections :\n");
	if (CURDLEV+1<=DEBUGLEV)  for (i=0; i<ro; i++)  for (j=i+1; j<ro; j++)
		SDEBUG(CURDLEV+1,"\t\t [#%d (%d), #%d (%d)]:\t ppat=%lX mpat=%lX,  det1=%s  det2=%s\n",
			i,ROWSID(ex,rcon[i]), j,ROWSID(ex,rcon[j]), ppat[i][j],mpat[i][j],
			pfield_printvalue(16,EXPM(edet,i,i+j*ro),SIGNM(edet,i,i+j*ro)),
			pfield_printvalue(16,EXPM(edet,j,i+j*ro),SIGNM(edet,j,i+j*ro)));
	DEBUG(CURDLEV+1,"Done with \"connecting\" subdeterminants.\n");
#endif
	return 1;
}


/**
 * This function implements one choice of an entry in gener_matextens_ext() above.
 * One call decides an entry at the row k of the extended column of edet.
 * The values selsg[k],selex[k] tell the last/next choice of a sign and exponent for the entry.
 * (Values selsg[k],selex[k]==-1 start choosing.)
 * The "connecting subdeterminants" (see gener_matextens_cdet()) edet and path patterns ppat[]
 * are used to guess allowable entries over a partial field.
 * Moreover, the function pre-checks representability over a partial field using other
 * suitable connecting paths.
 * If sel1>0 is given, then the only nonzero entry choice is 0 -- used for the first
 * nonzero in each component, but that is handled above.
 * If rnf>0 is given, then only random subselection is considered, see below.
 * 
 * The return value -1 means to go back in backtracking (no more choices here),
 * 0 means to try a new iteration,
 * and 1 means a successful choice -- to go forward or record the whole extension.
**/

int	gener_matextens_cycl(int k, long **ppat, ematrix *edet,
					int selsg[], int selex[], int sel1, int rnf) {
	int	i,j,jf=0,r, ro,cco;
	long	zcp;
	exp_t	xx,xo;
	sign_t	gg,go;
	
	if (pfispartial && (!edet || !ppat))  {PROGERROREXIT("Must get edet when in partial field!");}
	ro = ROWSM(edet);  cco = COLSM(edet)-1;
	for (zcp=0,i=0; i<ro; i++)		/* (get zero pattern of the extended column) */
		if (i>=k || SIGNM(edet,i,cco)!=0)  zcp |= (1l<<i);
	
		/**
		 * The current values of selsg[k] and selex[k] tell what is our next choice:
		 * selsg[k]==1 chooses a nonzero entry indexed by selex[k] (1...max) in
		 * pfield_get[pa]value(), and selsg[k]==0 chooses a zero entry.
		 * The select values of -1 mean the initial choice for selsg[k] and selex[k].
		 * Here we (possibly) move to the next selection of the sign, but we set the
		 * number of exponent choices below.
		 * We also try a random jump on the exponent choices for rnf>0.
		**/
	if (selex[k]<=0)
		if (++selsg[k]>1) {	/* the next value of the sign (trying in order 0,1, but -1 is no longer separately used!) */
			DEBUG(CURDLEV+3,"%*s calling \"extcycle\" at k=%d:  ^back^\n",2*k,"",k);
			return -1;	/* (choices are exhausted after 0,1...) */
		}
	DEBUG(CURDLEV+3,"%*s calling \"extcycle\" at k=%d: selsg=%d, selex=%d, sel1=%d, zcp=%lX\n",2*k,"",k,selsg[k],selex[k],sel1,zcp);
		/**
		 * The meaning of selex[],selsg[] is explained above, and
		 * the meaning of edet and ppat,mpat for a partial field is explained in
		 * gener_matextens_cdet().
		 * Here we either choose a zero entry at k (selsg[k]==0), or we choose the
		 * fixed 1 as a nonzero entry (sel1>0), or we choose the
		 * selex[k]-th nonzero entry at k, possibly relative to the "connecting"
		 * subdeterminants stored in edet when in a partial field.
		 * In the latter case, we use the first induced path joining k with a previous
		 * nonzero across all internal zeroes.
		**/
	if (selsg[k]<=0) {
		selex[k] = 1;  r = 0;
		go = 0;  pfield_setzeroexp(&xo);
	} else if (sel1>0) {
		selex[k] = 1;  r = 0;
		go = 1;  pfield_setzeroexp(&xo);
	} else if (!pfispartial) {
		if (selex[k]<=0)  selex[k] = pfield_getavalue_num();
		if (rnf>0)  selex[k] -= RANDOM()%(rnf+1);
		r = (selex[k]<=0?-1: pfield_getavalue(selex[k],&xo,&go));
	} else {
		for (jf=0; jf<k; jf++)
			if ((zcp&(1l<<jf))!=0 && (zcp&ppat[jf][k])==0)  break;
		if (jf>=k)  {PROGERROR("No induced connecting path found from the next entry to a previously set nonzero!,"
					" k=%d, zcp=%ld\nHave you already set some element in this component?",k,zcp); return -1;}
		pfield_mul(1,EXPM(edet,jf,cco),SIGNM(edet,jf,cco),
				1,EXPM(edet,k,jf+k*ro),SIGNM(edet,k,jf+k*ro),&xx,&gg);
		pfield_mul(-1,EXPM(edet,jf,jf+k*ro),-SIGNM(edet,jf,jf+k*ro),
				1,xx,gg,&xx,&gg);
		if (selex[k]<=0)  selex[k] = pfield_getpvalue_num(xx,gg);
		if (rnf>0) if (RANDOM()%(pfnumexp+1)==0)
			selex[k] -= RANDOM()%(rnf+1);
		r = (selex[k]<=0?-1: pfield_getpvalue(selex[k],xx,gg,&xo,&go));
	}
	--selex[k];			/* (the next choice for the exponent on the next call) */
	SETEXSIGM(edet,k,cco,xo,go);
	
		/**
		 * After choosing the extended entry, we try few more pfiel-oriented tests.
		 * We use the remaining (suitable) connecting determinants, for pfield-def checks.
		 * If an undefined subdeterminant is found here, then this extension is dead.
		**/
	if (go!=0 && !sel1 && pfispartial && edet && ppat) {
		for (j=0; j<k && r>=0; j++) {
			if ((zcp&(1l<<j))!=0 && (zcp&ppat[j][k])==0)
				r = ematrix_determinant_2x2check(edet,j,k,j+k*ro,cco);
			if (j<=jf && r<0)  {PROGERROR("Wrong setting of the new entry ?!?");}
		}		/* the first choice (at row jf found above) is supposed to be set OK */
	}
	if (r<0)  return 0;		/* (if the current choice at row k has failed, then try again) */
	DEBUG(CURDLEV+2,"%*s    -> chosen at #%d : %s (selex %d)\n",2*k,"",k,pfield_pvalue(16,xo,go),selex[k]+1);
	return 1;		/* successfull extension is returned */
}


/**
 * This function copies the last column (#cco) from edet to the last column (#co) of ee.
 * Then the copied column is scaled so that each connected component of rows of ee
 * starts with 1 as the first nonzero (each component scaled separately).
 * Here rcpi[] gives the indices of components for the rows, where rcpi[i] is the index
 * for the row rcon[i].
 * However, the words "the first nonzero" refer to the basic order of rows 0,1,...
 * 
 * The total number of rows of the matrices is ro.
 * But if k<ro is given, then only the first k rows are considered in order rcon[0],...
 * rcon[k-1], and the scaling is done only for those components for which the first
 * nonzero element can be determined (recall that rcon[] jumps through the rows,
 * and some rows missed by ..rcon[k-1] may apper before those determined entries).
 * The return value is the number of components that were successfully scaled.
 * (This applies mainly to k<ro...)
**/

int	gener_matextens_cpscal(int k, ematrix *edet, int cco,
				int rcon[], int rcpi[], ematrix *ee, int ro, int co) {
	int	i,j,p,ic;
	exp_t	xx;
	sign_t	gg;
	
	DEBUG(CURDLEV+1+(k<ro),"  calling unit-scale at level %d (ro=%d)\n",k,ro);
	if (k>ro)  k = ro;
	for (i=0; i<ro; i++)  SIGNM(ee,i,co) = 0;
	for (j=0; j<k; j++)  COPYEXSIGM(ee,rcon[j],co,edet,j,cco);
			/**
			 * Here we go through all defined entries of ee in order rcon[]
			 * (corresponds to the first k rows of edet), and we look for
			 * the first nozero defined entry in each component.
			 * We also consider a component finished if we pass all zero entries in it.
			 * The numnber of finished components is stored in ic.
			**/
	for (i=ic=0; i<k; i++) {
		if (i==0? rcpi[0]!=1: rcpi[i-1]>rcpi[i])  {PROGERROR("The components of rows of ee must be indexed in order from 1!, i=%d",i); break;}
		if (rcpi[i]<=ic)  continue;
		if (rcpi[i]>ic+1)  break;
			/* this cycle indentifies whether rcon[i] is the first nonzero in its component: */
		for (j=p=0; j<ro; j++) if (rcon[j]<rcon[i])
			if (rcpi[j]!=rcpi[i] || (j<k && SIGNM(ee,rcon[j],co)==0))
				p++;
		if (p<rcon[i] || (SIGNM(ee,rcon[i],co)==0 && i<ro-1 && rcpi[i+1]==rcpi[i]))
			continue;	/* (we finish also a full component of zeroes) */
		ic = rcpi[i];	/* (actually, always makes ic++, look above) */
			/* the next part scales all entries in the same component, up to rcon[k-1]: */
		xx = EXPM(ee,rcon[i],co);  gg = SIGNM(ee,rcon[i],co);
		if (gg==0 || (gg==1 && pfield_iszeroexp(xx)))  continue;
		DEBUG(CURDLEV+3,"  scaling component %d at first nonzero #%d(%d)\n",rcpi[i],i,rcon[i]);
		for (j=0; j<k; j++) if (rcpi[j]==rcpi[i])
			pfield_mul(1,EXPM(ee,rcon[j],co),SIGNM(ee,rcon[j],co),
					-1,xx,gg,&EXPM(ee,rcon[j],co),&SIGNM(ee,rcon[j],co));
		if (SIGNM(ee,rcon[i],co)!=1 || !pfield_iszeroexp(EXPM(ee,rcon[i],co)))  {PROGERROR("Wrong component scaling ?!?");}
	}
	return ic;
}



/**
 * This is a supplementary function for an easy printing of the extended line in the matrix e.
 * The last row is printed (as comma-separated numbers) if trq==0, the last column if trq==1.
 * The resulting string is returned in a static buffer - do not free(!!!) it.
**/

char*	gener_extprintline(int ch, ematrix *e, int trq) {
	int		i;
	static char	buf[200];
	
	if (!trq)  ematrix_transpose(e);
	buf[0] = 0;
	for (i=0; i<ROWSM(e) && strlen(buf)<160; i++) {
		if (ch<0)  pfield_pvalue_to(buf+strlen(buf),7,EXPM(e,i,COLSM(e)-1),SIGNM(e,i,COLSM(e)-1));
		else  pfield_printvalue_to(buf+strlen(buf),(ch>0?14:9),EXPM(e,i,COLSM(e)-1),SIGNM(e,i,COLSM(e)-1));
		if (i<ROWSM(e)-1)  strcat(buf,",");
	}
	if (i<ROWSM(e))  strcat(buf,"...");
	if (!trq)  ematrix_transpose(e);
	return buf;
}

































