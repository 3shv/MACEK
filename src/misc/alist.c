
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
 * alist implementation....
 * 
 * 
**/



#include "macek.h"













/******************	Handling arbitrary lists	*******************/
/**************************************************************************/


#undef CURDLEV
#define CURDLEV	7


/**
 * A "list" is a structure of general type  **list_type, actually terminated by NULL.
 * (However, to avoid pointer compatability problems, we transfer lists also as *void...)
 * To allow lists to grow, we store additional information right after the terminating NULL:
 *  If ls[i]==NULL, then ls[i+1] is a pointer to the start of the list, and ls[i+2] (as long)
 *  is the length remaining after the end of the list (minus the terminating NULL and these two).
 * A special magic number is stored into one element just before the list starts,
 * and a NULL value and a pointer to the terminating NULL are stored at positions -3,-2.
 * Moreover, a list may be given pointing to an arbitrary internal element, and the code
 * here can recover the list's true start and end.
 * 
 * No internal element of a list may be set to NULL !!!
 * 
 * The internal function alist_getend() returns the end pointer (terminating NULL) of the
 * given list.
 * To speed-up this operation, we use a simple hash-table storing, for given list pointer
 * in alisttb_in[], the pointer to the true start of the list (if the list pointer is not
 * the true start, alisttb_mark[]==1), or the pointer to the terminating NULL of the list
 * (if the list pointer is the true start, alisttb_mark[]==2).
 * 
 * The internal function alist_setend() sets in the hash-table the end of the given list.
 * (If not the true start of the list is given, then it stores the true start first,
 * then the list end.)
 * If no etn (==NULL) is given, then it clears the apropriate entry from the table.
 * 
 * The functions alist_alloc() and alist_free() are provided for (re)allocating
 * new lists and freeing them; alist_free() is for public use.
 * These functions keep care of proper extra structure of lists, and of setting
 * the list end references in the hashtable via alist_setend().
 * 
 * The function alist_finishall() is for calling at the finish of the whole program.
**/

#ifndef	ALIST_MAXLIST
#define	ALIST_MAXLIST	5000000l	/* will not create longer list to save memory */
#endif
#ifndef	ALIST_HASHSZ
#define	ALIST_HASHSZ	3507		/* the table of list ends, and its size */
#endif
#define	ALISTHASHP(l)	HASHPOINTER((l),ALIST_HASHSZ)

static short	alisttb_mark[ALIST_HASHSZ],
		alisttb_first = 1;
static void	**alisttb_in[ALIST_HASHSZ],
		**alisttb_out[ALIST_HASHSZ];
static int	alistalloc = 0, alisthashend = 0,  alistslowend = 0;


void**	alist_getend(void **ls) {
	void	**t, **x,**y;
	int	i;
	
	if (alisttb_first) {	/* (first run of the function) */
		alisttb_first = 0;
		for (i=0; i<ALIST_HASHSZ; i++)  alisttb_mark[i] = 0;
	}
	
	if (!ls)  return NULL;
	if (!*ls)  return ls;	/* possible (lucky) quick access to the list end */
	if (!ls[-3] && ls[-1]==(void*)12345678l) {
		x = ls[-2];  t = x[1];
		if (x[0]!=NULL || (t?t[-1]:t)!=(void*)12345678l)  {PROGERROREXIT("Fatal error - probably broken list %p (%p %p) !",ls,x,t);}
		return x;
	}
#if ALIST_HASHSZ>1
		/**
		 * Here we look into the hash-table, whether  ls  is at position i there.
		 * If alisttb_mark[i]==1 (stores the true start of the list), then we look
		 * into the table again to find the end of the list.
		 * If alisttb_mark[i]==2, we have got the end, which is then returned.
		**/
	alisthashend++;
	i = ALISTHASHP(ls);
	if (alisttb_mark[i]>0 && alisttb_in[i]==ls) {
		if (alisttb_mark[i]==1) {
			t = alisttb_out[i];  i = ALISTHASHP(t);
			if (alisttb_in[i]!=t)  i = -1;
		}
		if (i>=0)  if (alisttb_mark[i]==2) {
			x = alisttb_out[i];  t = x[1];
			if (x[0]!=NULL || t[-1]!=(void*)12345678l)  {PROGERROREXIT("Fatal error - probably broken list %p !",ls);}
			return x;
		}
	}
#endif
		/**
		 * This is the default simple (but slow) catch-up - find the terminating
		 * NULL of the list, and return the pointer (plus some checks).
		**/
	alistslowend++;
	for (x=ls; *x; x++) ;
	t = x[1];
	if (t[-1]!=(void*)12345678l)  {PROGERROREXIT("Fatal error - probably broken list %p: %p !",ls,t[-1]);}
	return x;
}

void	alist_setend(void **ls, void **et, void **etn) {
	long    x,*y;
	void    **st,**z;
	int	i;
	
	if (alisttb_first)  alist_getend(NULL);	/* (first run of the function) */
	if (!et)  et = etn;
	st = et[1];  x = *((long*)(et+2));
	if (ls<st || ls>et || et[0]!=NULL || et[-1]==NULL
			|| st[-2]!=et || st[-1]!=(void*)12345678l || x<0)
		{PROGERROREXIT("Fatal error - probably broken list %p (st %p, et %p, etn %p)!",ls,st,et,etn);}
	if (x>ALIST_MAXLIST)
		{PROGERROREXIT("Fatal error - too long %ld list %p (st %p, et %p, etn %p)!",x,ls,st,et,etn);}
		/**
		 * Here we set the necessary extra data after the list end and at the start.
		 * The new list end is given in etn, and we find out the new free (remaining)
		 * size from the difference between et and etn.
		 * The extra elements of the list are filled with non-NULL garbage.
		**/
	if (etn) {
		etn[0] = NULL;  etn[1] = st;
		y = (long*)(etn+2);  *y = x+(et-etn);
		st[-2] = etn;
		for (z=et; z<etn; z++)  *z = et[-1];
		if (*y<0 || etn[-1]==NULL)  {PROGERROREXIT("Fatal error - probably broken list %p %ld",ls,*y);}
	}
		/**
		 * If the given ls is not the true start of the list (as found at the end et),
		 * then we must first store a pointer to the true start, and mark 1.
		**/
	if (etn && ls!=st) {
		i = ALISTHASHP(ls);
		alisttb_mark[i] = 1;
		alisttb_in[i] = ls;  alisttb_out[i] = st;
	}
		/**
		 * If not etn is given, then we clear the table entry  st  if it exists.
		 * Otherwise, we store the pointer to the end of the list, and mark 2.
		**/
	i = ALISTHASHP(st);
	if (etn==NULL) {
		if (alisttb_in[i]==st)  alisttb_mark[i] = 0;
	} else {
		alisttb_mark[i] = 2;
		alisttb_in[i] = st;  alisttb_out[i] = etn;
	}
}


void**	alist_alloc(void **ls, long sz) {
	void	**pt, **et,**st;
	
		/**
		 * If ls==NULL is given, then we simply allocate a new list.
		 * Otherwise, we reallocate the list ls if there is less than sz space left.
		 * When reallocating, we take proper care of the list structure and
		 * of the hashtable entries for the new list.
		**/
	if (!ls) {
		alistalloc++;
		pt = MMALLOC((sz+10)*sizeof(pt[0]));  pt += 3;
		pt[-1] = (void*)12345678l;  pt[-3] = pt[0] = NULL;
		ls = st = et = pt;
	} else {
		et = alist_getend(ls);  st = et[1];
		if (sz<=*((long*)(et+2)))  return ls;
		sz += (et-st);
		alist_setend(ls,et,NULL);
		pt = MREALLOC(st-3,(sz+(et-st)+10)*sizeof(pt[0]));  pt += 3;
		ls = pt+(ls-st);
		et = pt+(et-st);  st = pt;
	}
	st[-2] = et;
	et[1] = st;  *((long*)(et+2)) = sz;
	alist_setend(ls,NULL,et);
	return ls;
}

void	alist_free(void *ls) {
	void	**et,**st;
	
	if (!ls)  return;
	et = alist_getend(ls);  st = et[1];
	if (st[-1]!=(void*)12345678l || st[-3])  {PROGERROREXIT("Fatal error - probably broken list %p: %p %p !",ls,st[-1],st[-3]);}
	alist_setend(ls,et,NULL);	/* (to clear the hash-table entry of ls, also consistency check) */
	st[-1] = NULL;
	alistalloc--;
	FREE(st-3);
}


void	alist_finishall(void) {
	alisttb_first = 1;
	if (alisthashend>100)  DEBUG(CURDLEV-1,"Alist_ has used %d hash-calls to look for a list end.\n",alisthashend);
	if (alistslowend>10)  DEBUG(CURDLEV-2,"Alist_ has used %d SLOW calls to look for a list end.\n",alistslowend);
	if (alistalloc>5)  DEBUG(CURDLEV-3,"Alist_ finished with %d lists not freed!\n",alistalloc);
}




/**
 * Here we create and dispose lists - these functions are for general use.
 * - new_alist() creates a new empty list of the given initial size.
 * - alist_free() above frees the given list (possible to give any element as a pointer),
 *  but NOT the list elements.
 * - dispose_alist_ext() frees all the given list elements (using an optional disposing function
 *  given in *ff), and then frees the list itself.
**/

void*	new_alist(long sz) {
	void	**pt;
	pt = alist_alloc(NULL,sz);
	return pt;
}

void	dispose_alist_ext(void *ls, void (*ff)(void*)) {
	void	**pt,**et,**x;
	
	if (!ls)  {PROGERROREXIT("Cannot dispose NULL list!"); return;}
	et = alist_getend(ls);  pt = et[1];
	for (x=pt; x?*x:0; x++) {	/* frees all elements of the list, possibly using supplied (*ff)() */
		if (ff)  (*ff)(*x);
		else  FREE(*x);
	}
	if (x!=et)  {PROGERROREXIT("Fatal error - probably broken list %p: %p %p !",ls,x,et);}
	alist_free(et);
}


/**
 * Here we get the start position of a list (recall than a list may be referred to
 * by any internal position), and the total length.                               
**/

void*	alist_getstart(void *ls) {
	void	**et;
	
	if (!ls)  return NULL;
	et = alist_getend(ls);
	return et[1];
}

int	alist_getlength(void *ls) {
	void	**et;
	
	if (!ls)  return 0;
	et = alist_getend(ls);
	return et-(void**)et[1];
}

int	alist_getlength_part(void *ls) {
	void	**et;
	
	if (!ls)  return 0;
	et = alist_getend(ls);
	return et-(void**)ls;
}


/**
 * This function appends a new element ap to the end of the given list.
 * (Again, it is possible to give any element of the list as a pointer ls.)
 * If the list is not sufficiently large, it is reallocated to twice its size.
 * The list is returned back (possibly different address, so it must be stored!).
**/

void*	alist_append(void *ls, void *ap) {
	void	**et;
	
	if (!ap)  {PROGERROR("Cannot append NULL element to a list!"); return ls;}
	ls = alist_alloc(ls,1);
	et = alist_getend(ls);
	alist_setend(ls,et,et+1);
	et[0] = ap;
	return ls;
}

void*	alist_applist(void *ls, void *apl) {
	void	**x;
	
	if (!ls)  return apl;
	for (x=apl; x?*x:0; x++)  ls = alist_append(ls,*x);
	if (apl)  alist_free(apl);
	return ls;
}


/**
 * The function alist_delete() deletes one element stored at address ls from the list
 * (the element itself is not freed).
 * The remaining elements in the list are simply shifted by one to the front.
 * The list is returned back (address is unchanged).
 * 
 * The function alist_delete_val() deletes all elements of the list starting from ls
 * that point to the given value (pointer) val.
 * The list is returned back (address is unchanged).
**/

void*	alist_delete(void *ls) {
	void	**et;
	
	if (!ls)  {PROGERROR("Cannot delete NULL position from a list!"); return ls;}
	for (et=ls; *et; et++)  if (et[1])  et[0] = et[1];
	alist_setend(ls,et,et-1);
	return ls;
}

void*	alist_delete_val(void *ls, void *val) {
	void	**et;
	int	i;
	
	if (!ls || !val)  {PROGERROR("Cannot delete NULL value from a list!"); return ls;}
	for (et=ls,i=0; et[i]; ) {
		if (et[i]!=val) {
			if (i>0 && et[i])  et[0] = et[i];
			et++;
		} else  i++;
	}
	if (i>0)  alist_setend(ls,et+i,et);
	return ls;
}


/**
 * This function inserts one element to the position pointed by ins in the list.
 * The list is possibly enlarged like in alist_append(), and the remaining elements
 * are shifted by one to the end.
 * The list is returned back (possibly different address, so it must be stored!).
**/

void*	alist_insert(void *ls, void *ins, void *ai) {
	void	**et;
	long	di;
	
	if (!ins || !ai)  {PROGERROR("Cannot insert NULL element or NULL position to a list!"); return ls;}
	if (!ls)  ls = alist_getstart(ins);
	di = (void**)ins-(void**)ls;	/* (must keep the difference, since ls may change next) */
	ls = alist_append(ls,ai);
	et = alist_getend(ls);  et--;
	if (di<0 || di>et-(void**)ls)  {PROGERROR("The insert position %p (%ld) is not in the list %p-%p?!",ins,di,ls,et);}
	for ( ; et-(void**)ls>di; et--)  et[0] = et[-1];
	et[0] = ai;	/* the inserting position is now prepared */
	return ls;
}


/**
 * This function returns a copy of the given list.
 * (The list ls may be given by an arbitrary internal position.)
**/

void*	alist_copy(void *ls) {
	void	**x, **st,**ln;
	
	if (!ls)  {PROGERROR("Cannot copy NULL list!"); return ls;}
	ln = new_alist(20+alist_getlength(ls));
	st = alist_getstart(ls);
	for (x=st; *x; x++)  ln = alist_append(ln,*x);
	return  ln+ ((void**)ls-st);
}

































