
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
 * This source holds some general-use variables and functions for the Macek project.
 * They are not related directly to the project, and they may be useful elsewhere...
 * 
**/












/******************	Handling arbitrary lists	*******************/
/**************************************************************************/


/**
 * A "list" is a structure of general type  **list_type, terminated by NULL.
 * (However, to avoid pointer compatability problems, we transfer lists also as *void...)
 * 
 * To allow lists to grow, we store additional information right after the terminating NULL:
 *  If ls[i]==NULL, then ls[i+1] is a pointer to the start of the list, and ls[i+2] (as long)
 *  is the length remaining for the list growth (minus the terminating NULL and these two).
 * A special magic number is stored into the element -1 just before the list starts,
 * and a NULL value and a pointer to the terminating NULL are stored at positions -3,-2.
 * However, user need not worry about this internal list structure at all.
 * 
 * Moreover, a list may be given pointing to an arbitrary internal element, and the code
 * here can recover the list's true start and end.
 * No internal element of a list may be set to NULL !!!
**/


/**
 * Here we create and dispose lists - these functions are for general use.
 * - new_alist() creates a new empty list of the given initial size.
 * - alist_free() frees the given list (possible to give any element as a pointer), but
 *   NOT(!!!) the list elements.
 * - dispose_alist_ext() frees all the given list elements (using an optional disposing function
 *   given in *ff), and then frees the list itself.
 *   (All list elements from the start of the list are freed here!!!)
**/

void*   new_alist(long sz) ;
void    alist_free(void *ls) ;
void    alist_finishall(void) ;

void    dispose_alist_ext(void *ls, void (*ff)(void*)) ;

#define	dispose_alist_copts(ls)	dispose_alist_ext(ls,(void (*)(void*))&dispose_comoption)
#define	dispose_alist_mats(ls)	dispose_alist_ext(ls,(void (*)(void*))&dispose_ematrix)
#define	dispose_alist_frams(ls)	dispose_alist_ext(ls,(void (*)(void*))&dispose_frame_ls)
	/* (if NULL disposing function for the elements is given, then free() is used) */
#define	dispose_alist(ls)	dispose_alist_ext(ls,NULL)


/**
 * Here we get the start position of a list (recall than a list may be referred to
 * by any internal position), and the total length.
 * 
 * Then a function appends a new element ap to the end of the given list.
 * (Again, it is possible to give any element of the list as a pointer ls.)
 * If the list is not sufficiently large, it is reallocated to twice its size.
 * The list is returned back (possibly different address, so it must be stored!).
 * 
 * Next similar function appends a whole list apl to the given list ls, and frees
 * the list apl (not its elements).
**/

void*   alist_getstart(void *ls) ;
int     alist_getlength(void *ls) ;
int     alist_getlength_part(void *ls) ;

void*   alist_append(void *ls, void *ap) ;
void*   alist_applist(void *ls, void *apl) ;


/**
 * This function deletes one element pointed by ls from the list (the element itself is not freed).
 * The remaining elements in the list are simply shifted by one to the front.
 * The list is returned back (address is unchanged).
 * The function alist_delete_val() deletes all elements of the list starting from ls
 * that point to the given value (pointer) val.
 **/

void*   alist_delete(void *ls) ;
#define	alist_delete_last(ls)	alist_delete(((void**)(ls))+alist_getlength(ls)-1)

void*   alist_delete_val(void *ls, void *val) ;


/**
 * This function inserts one element to the position pointed by ins in the list.
 * The list is possibly enlarged like in alist_append(), and the remaining elements
 * are shifted by one to the end.
 * The list is returned back (possibly different address, so it must be stored!).
 * 
 * The function alist_copy() returns a copy of the given list.
 * (The list ls may be given by an arbitrary internal position.)
**/

void*   alist_insert(void *ls, void *ins, void *ai) ;

void*   alist_copy(void *ls) ;









/******************	Recycling lists 	*******************/
/******************************************************************/


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


/**
 * The following defines restrict the number of recycling classes, the size of one table
 * for recycled lists, the offset tried in the table if the hash position is occupied,
 * and the maximal number of objects in a list identification.
 * 
 * The next structure describes one recycling class - its index, name, the disposing
 * function for list elements, number of ident objects, and copying, disposing, and
 * comparing functions for the identification objects (separately for each one).
 * The structure reuseid is then given as the identification of a list.
 * 
**/

enum	AL_CLASSES { AL_SELFMAPS=0, AL_BASES, 	AL_LASTX };
#define	AL_MAXCLASSES	AL_LASTX

#define	AL_ONETABLE	1901	/* (use strange odd numbers here; or 0 to disable caching) */
#define	AL_HASHOFFS	5
#define	AL_MAXREUSEID	3
#define	AL_MAXLTOTAL	300000	/* the maximal sum of list length that may be stored in one cache */


typedef struct {
	
	int	ix;
	char	*name;
	void	(*displs)(void*);	/* (the disposing function for an element of the list, free() if NULL) */
	int	idlen;
	void	*(*copyid[AL_MAXREUSEID])(void*);
	void	(*dispid[AL_MAXREUSEID])(void*);
	int	(*cmpid[AL_MAXREUSEID])(void*,void*);
	
}	reusedescr;

typedef struct {
	void	*id[AL_MAXREUSEID];
}	reuseid;



/**
 * These are the different ways how to call the recycling function alist_reuse().
 * This function is not really for a public use, you should prepare separate recycling
 * functions for each recycled entity, which then call these functions.
 * The type of recycled entity is given by the structure descr for each type.
 * 
 * Read more in  ../misc/alreuse.c .
**/

void*   alist_reuse(reusedescr *descr, int what, int idn, reuseid ids, long hash, void *al) ;
#define	alist_reuse_set(descr,ids,hash,al)	alist_reuse(descr,1,0,ids,hash,al)
#define	alist_reuse_rec(descr,idn,al)		alist_reuse(descr,-1,idn,idsnull,0l,al)
#define	alist_reuse_get(descr,idn,ids,hash)	alist_reuse(descr,2,idn,ids,hash,NULL)
#define	alist_reuse_null(descr)			alist_reuse(descr,0,1,idsnull,0l,NULL)
#define	alist_reuse_cnull(descr)		alist_reuse(descr,0,0,idsnull,0l,NULL)

extern reuseid	idsparam, idsnull;











/******************	Other...................	*******************/
/**************************************************************************/


/**
 * 
 * 
 * some global variables, used mainly in debug and error messages....
 * 
**/


extern	int	junk;
extern	char    *ssss,  debugoutbuf[10000];

extern	int	uerror_occured, perror_occured;
extern	int	crash_all_errors;

extern	FILE	*errorout, *debugout, *printout;

extern  char	*printoutpref;


/**
 * 
 * 
 * some simple file functions........
 * 
**/


char*   fname_extract(char *fn) ;

char*   fdir_extract(char *fn) ;


FILE*   fopen_path(char *fn, char *mod, char **path, char *dir, char *ext) ;

FILE*   fopen_dir(char *fn, char *mod, int wr, int pre, int pix) ;
int     fopen_safename(char *fn, int wr) ;


/**
 * 
 * 
 * some other functions........
 * 
**/

void*   malloc_twodim(int sz, int d1, int d2) ;
void*   malloc_threedim(int sz, int d1, int d2, int d3) ;


































