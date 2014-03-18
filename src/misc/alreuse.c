
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
 * 
 * various general functions for the package.....
 * 
 * 
**/



#include "macek.h"














/******************	Recycling arbitrary lists	*******************/
/**************************************************************************/


#undef CURDLEV
#define CURDLEV		6		/* (need <=DEBUGLEV+2 to perform paranoic tests) */


/**
 * We provide a general way of "recycling" arbitrary lists generated in the program.
 * Instead of disposing a list, we cache it here together with copies of its
 * identification objects (like matrices from which the list was generated).
 * Next time, when comparable objects are to generate the list again, the cached
 * one may be returned.
 * 
 * There is a number of "classes" of lists available for recycling.
 * Each class is described by the structure reusedescr which carries an index
 * in 0...AL_MAXCLASSES-1 for the class, and information about the kind and number
 * of identification objects used for each list in the class.
 * The list cache is then stored in static variables in  alreuse.c .
 * 
**/


#if	AL_ONETABLE>1

/**
 * These are the static variables used to store the list caches.
 * The alr_*save variables are used to store information about (just) generated lists
 * to allow their recycling later:
 *  This is the suggested hash for storage position of the list once it is returned,
 *  the list pointer and length, and the identification structure together with
 *  the same structure carrying copies of the identification objects (originals may change).
 * The alr_*cache variables are used to store the cached lists, their identification
 *  structures and copies, and an additional arbitrary description number for the list
 *  (which affects the storing position as well).
 * 
 * Read below in alist_reuse() about the use of these variables.
**/

static int	alr_firstcallc[AL_MAXCLASSES], alr_firstcall = 1;
static int	alr_usedn[AL_MAXCLASSES], alr_usedlen[AL_MAXCLASSES];
static int	alr_calls[AL_MAXCLASSES], alr_callrec[AL_MAXCLASSES];
static int	alr_callreg[AL_MAXCLASSES], alr_callbad[AL_MAXCLASSES], alr_callbk[AL_MAXCLASSES];

static long		alr_hsave[AL_MAXCLASSES][AL_ONETABLE+2];
static void**		alr_gsave[AL_MAXCLASSES][AL_ONETABLE+2];
static int		alr_lsave[AL_MAXCLASSES][AL_ONETABLE+2];
static reuseid		alr_idsave[AL_MAXCLASSES][AL_ONETABLE+2];
static reuseid		alr_icsave[AL_MAXCLASSES][AL_ONETABLE+2];
static void**		alr_gcache[AL_MAXCLASSES][AL_ONETABLE+2];
static int		alr_ncache[AL_MAXCLASSES][AL_ONETABLE+2];
static reuseid		alr_idcache[AL_MAXCLASSES][AL_ONETABLE+2];
static reuseid		alr_iccache[AL_MAXCLASSES][AL_ONETABLE+2];

reuseid		idsparam, idsnull;	/* (used in macros) */


void	alist_clearinit(reusedescr *descr, int lev) {
	int		i;
	void**		*gsave = alr_gsave[descr->ix];
	void**		*gcache = alr_gcache[descr->ix];
	
	if (lev>1) {
		for (i=0; i<AL_ONETABLE; i++)  gsave[i] = NULL;
		i = descr->ix;
		alr_callrec[i] = alr_calls[i] = alr_callbad[i] = alr_callreg[i] = alr_callbk[i] = 0;
	}
	for (i=0; i<AL_ONETABLE; i++)  gcache[i] = NULL;
	//alr_firstcallc[descr->ix] = 0;
	alr_usedn[descr->ix] = alr_usedlen[descr->ix] = 0;
}


/**
 * This function allows to "re-use" once computed arbitrary list.
 * For this to work, one has first to call this function (what==1)
 * from the generator function of the list
 * with the generated list in  al  and with the idn,ids and hash information
 * (hash is not repeated in the recycling call, so it must be rememered here).
 * 
 * Then the user should call this function (what==-1) with the list to recycle in  al,
 * and with the same idn number (ids or hash is not taken here).
 * 
 * After that, the next call (what==2) with the same idn,hash and with comparable ids
 * should recover the cached list and return it back, or NULL if not found.
 * The call (what==1) must be repeated after that if next recycling should take place.
 * See the structure reusedescr for ways how the identification objects are compared.
 * 
 * Finally, the call of (what==0) erases and frees all stored data in the cache.
 * (Must be called for each class separately.)
 * 
**/

void*	alist_reuse(reusedescr *descr, int what, int idn, reuseid ids, long hash, void *al) {
	int		i,j,k,l,q, pos,clidi;
	void		*clls, **list = al;
	reuseid		ics, clid;
	
		/**
		 * These are local access pointers to the above variables of the list cache.
		 * They are determined by the list class given in descr->ix.
		 * Next, we initialize (clear) the cache for the first time.
		**/
	long		*hsave = alr_hsave[descr->ix];
	int		*lsave = alr_lsave[descr->ix];
	void**		*gsave = alr_gsave[descr->ix];
	reuseid		*idsave = alr_idsave[descr->ix], *icsave = alr_icsave[descr->ix];
	void**		*gcache = alr_gcache[descr->ix];
	int		*ncache = alr_ncache[descr->ix];
	reuseid		*idcache = alr_idcache[descr->ix], *iccache = alr_iccache[descr->ix];
	
	if (alr_firstcall) {			/* first-call initialization of the recycling system */
		for (i=0; i<AL_MAXCLASSES; i++)  alr_firstcallc[i] = 2;
		alr_firstcall = 0;
	}
	if (alr_firstcallc[descr->ix]>0) {	/* first-call initialization of this class */
		alist_clearinit(descr,alr_firstcallc[descr->ix]);
		alr_firstcallc[descr->ix] = 0;
	}
	clidi = 0;  clls = NULL;
					/* too long list for recycling is simply disposed */
	if (what==-1) if (alist_getlength(list)>AL_MAXLTOTAL/4) {
		clls = list;  what = -111;
	}
					/* when the cache grows too much, then it is cleared first */
	if (what==-1) if (alr_usedn[descr->ix]>AL_ONETABLE/2
			|| alr_usedlen[descr->ix]+alist_getlength(list)>AL_MAXLTOTAL) {
		alist_reuse_cnull(descr);	/* (only gcache is cleared here, not gsave!) */
	}
	DEBUG(CURDLEV+1+(what<-1),"Calling %s [%d/%d] list reuse at what==%d, al = %p.\n",descr->name,alr_usedn[descr->ix],alr_usedlen[descr->ix],what,al);
	
		/**
		 * This part is called by the generator function to register a list for recycling.
		 * The pointer and length of the list, and the identification and suggested
		 * hash for (later) storing of the list are stored here, for later use.
		 * The identification objects from ids must be copied locally because they may change.
		**/
	if (what==1) {
		DEBUG(CURDLEV+1,"Registering a %d %s list %p for recycling, hash=%ld, id=(%p,%p)\n",alr_callreg[descr->ix]+1,descr->name,list,hash,ids.id[0],descr->idlen>1?ids.id[1]:NULL);
		alr_callreg[descr->ix]++;
		pos = HASHPOINTER(list,AL_ONETABLE);
		for (i=0; i<AL_ONETABLE/6 && gsave[(pos+i)%AL_ONETABLE]; i++) ;
		if (gsave[(pos+i)%AL_ONETABLE])  i = alr_callreg[descr->ix]%AL_ONETABLE/6;
		pos = (pos+i)%AL_ONETABLE;
		gsave[pos] = list;
		lsave[pos] = alist_getlength(list);
		hsave[pos] = hash;
		idsave[pos] = icsave[pos] = ids;
		for (j=0; j<descr->idlen; j++)  if (icsave[pos].id[j])
			icsave[pos].id[j] = (*descr->copyid[j])(icsave[pos].id[j]);
	}
		/**
		 * Here the recycled list is stored to the cache, to a position determined
		 * by the previous suggested hash, and by the current value of idn.
		 * Possible previous list and identification at this position are freed later.
		**/
#	define	AL_POSHASH(hh,dn)	((((unsigned)(hh)%AL_ONETABLE)*8+(unsigned)(dn))%AL_ONETABLE)
	if (what==-1) {
		alr_callbk[descr->ix]++;
		pos = HASHPOINTER(list,AL_ONETABLE);
		for (i=0; i<AL_ONETABLE/6 && gsave[pos]!=list; i++)
			pos = (pos+1)%AL_ONETABLE;
		if (gsave[pos]!=list || lsave[pos]!=alist_getlength(list)) {
			alr_callbad[descr->ix]++;
			DEBUG(CURDLEV,"Trying to recycle (%d) an unregistered or altered %s list %p, len %d.\n",alr_callbad[descr->ix],descr->name,list,lsave[pos]);
			clls = list;
		} else {
			hash = hsave[pos];
			ids = idsave[pos];  ics = icsave[pos];
			gsave[pos] = NULL;
			pos = AL_POSHASH(hash,idn);
			for (i=0; i<AL_ONETABLE/2+22 && gcache[(pos+i)%AL_ONETABLE]; i++) ;
			l = (pos+i)%AL_ONETABLE;
			if (gcache[l])  {PROGERROR("must find an empty spot among half of the table!");}
			if (gcache[pos]) {
				gcache[l] = gcache[pos];  ncache[l] = ncache[pos];
				idcache[l] = idcache[pos];  iccache[l] = iccache[pos];
				DEBUG(CURDLEV+1,"Old %s list %p moved for hash=%ld, idn=%d, pos=%d.\n",descr->name,gcache[pos],hash,idn,pos);
			}
			gcache[pos] = list;
			ncache[pos] = idn;
			idcache[pos] = ids;  iccache[pos] = ics;
			alr_usedn[descr->ix]++;
			alr_usedlen[descr->ix] += alist_getlength(list);
			DEBUG(CURDLEV,"Recycling a %s list %p at hash=%ld, pos=%d, id=(%p,%p), idn=%d.\n",descr->name,list,hash,pos,ids.id[0],descr->idlen>1?ids.id[1]:NULL,idn);
		}
	}
		/**
		 * This part is called to (possibly) retrieve a saved list from the position
		 * hash, identified by idn, ids.
		 * The list to be retrieved must have the same idn and the same each ident
		 * object according to comparing functions from descr->cmpid.
		 * The list is searched from the position given by hash and idn,
		 * where the first AL_HASHOFFS idents are compared all, and then next up
		 * to AL_HASHOFFS idents are compared if they have the same (ident) pointers.
		 * If the list is found, then it is returned and erased from the cache.
		**/
	if (what==2) {
		if (list)  {PROGERROR("no list is supposed to be given now");}
		pos = AL_POSHASH(hash,idn);
		list = NULL;
		for (i=k=l=0; i<AL_ONETABLE && l<2*AL_HASHOFFS; i++, pos=(pos+1)%AL_ONETABLE) {
			if (!gcache[pos] || ncache[pos]!=idn)
				continue;
			for (j=q=0; j<descr->idlen && !q; j++)
				if (idcache[pos].id[j]!=ids.id[j])  q = 1;
			if (k++>=AL_HASHOFFS && q)
				continue;
			if (!q)  l++;
			for (j=q=0; j<descr->idlen && !q; j++)
				if (!(*descr->cmpid[j])(iccache[pos].id[j],ids.id[j]))
					q = 1;
			if (!q) {
				list = gcache[pos];  gcache[pos] = NULL;
				clidi = 1;  clid = iccache[pos];
				alr_usedn[descr->ix]--;
				alr_usedlen[descr->ix] -= alist_getlength(list);
				DEBUG(CURDLEV,"Recycled %s list %p returned from pos=%d for id=(%p,%p), idn=%d.\n",descr->name,list,pos,ids.id[0],descr->idlen>1?ids.id[1]:NULL,idn);
				break;
			}
		}
		alr_calls[descr->ix]++;
		if (list)  alr_callrec[descr->ix]++;
		else  DEBUG(CURDLEV+1,"Recycled %s list for id=(%p,%p), idn=%d was NOT found.\n",descr->name,ids.id[0],descr->idlen>1?ids.id[1]:NULL,idn);
		SDEBUG(CURDLEV+1,"\t\t(Success rate for %s recycling: %d out of %d calls.)\n",descr->name,alr_callrec[descr->ix],alr_calls[descr->ix]);
	}
		/**
		 * Here we dispose all data (lists and identification copies) held
		 * in the whole table.
		 * If idn==1 is given, then even all saved list identifications are cleared.
		 * 
		 * The next part contains code disposing no longer used lists and idents,
		 * as given above in the variables clls, clid.
		**/
	if (what==0) {
		if (idn>(k=0))  for (i=0; i<AL_ONETABLE; i++)  if (gsave[i]) {
			k++;
			for (j=0; j<descr->idlen; j++)  if (icsave[i].id[j]) {
				if (descr->dispid[j])
					(*descr->dispid[j])(icsave[i].id[j]);
				else  FREE(icsave[i].id[j]);
			}
		}
		if (idn>0)  DEBUG(CURDLEV+1-2*(k>0),"Freeing %d saved %s list identifications...\n",k,descr->name);
		
		for (i=k=l=0; i<AL_ONETABLE; i++)  if (gcache[i]) {
			k++;  l += alist_getlength(gcache[i]);
			for (j=0; j<descr->idlen; j++)  if (iccache[i].id[j]) {
				if (descr->dispid[j])
					(*descr->dispid[j])(iccache[i].id[j]);
				else  FREE(iccache[i].id[j]);
			}
			dispose_alist_ext(gcache[i],descr->displs);
		}
		DEBUG(CURDLEV+1-2*(k>0),"Freeing %d cached %s lists of combined length %d...\n",alr_usedn[descr->ix],descr->name,alr_usedlen[descr->ix]);
		if (idn)  SDEBUG(CURDLEV+1-2*(alr_callrec[descr->ix]>0),"\t\t(Success rate for %s recycling: %d out of %d calls, %d reg calls, and %d bad rec calls.)\n",
				descr->name,alr_callrec[descr->ix],alr_calls[descr->ix],alr_callreg[descr->ix],alr_callbad[descr->ix]);
		if (alr_usedn[descr->ix]!=k || alr_usedlen[descr->ix]!=l)  {PROGERROR("wrong accounting in recycling cache! %d!=%d or %d!=%d",alr_usedn[descr->ix],k,alr_usedlen[descr->ix],l);}
		if (alr_callbad[descr->ix]>100)  DEBUG(CURDLEV-3,"Too many %d bad recycling calls occured ?!?\n",alr_callbad[descr->ix]);
		if (alr_callreg[descr->ix]-alr_callbk[descr->ix]>30)  DEBUG(CURDLEV-3,"Too many %d registered %s list were not recyced ?!?\n",alr_callreg[descr->ix]-alr_callbk[descr->ix],descr->name);
		alist_clearinit(descr,1);		/* (this does not clear gsave[]!) */
		
	} else {		/* cleaning the list clls and ident clid (that were set above): */
		if (clls!=NULL)
			dispose_alist_ext(clls,descr->displs);
		if (clidi>0)
			for (j=0; j<descr->idlen; j++)  if (clid.id[j]) {
				if (descr->dispid[j])  (*descr->dispid[j])(clid.id[j]);
				else  FREE(clid.id[j]);
			}
		if (clls==list)  list = NULL;
	}
	return list;
}


#else	/* of #if AL_ONETABLE>1 */

void*   alist_reuse(reusedescr *descr, int what, int idn, reuseid ids, long hash, void *al) {
	if (al && what==-1)  dispose_alist_ext(al,descr->displs);
	return NULL;
}

#endif































































