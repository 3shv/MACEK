
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
 * This header file "duplicates" some hard-limit and performance defines
 * from other source files.
 * It is enough to set a new value here to overwrite the default.
 * 
**/




#ifndef	PROFILE_H
#define	PROFILE_H






/******************	Hard-limit settings	*******************/
/******************************************************************/


/**
 * The hard limits listed here are (except PF_NUMEXPS, GEN_MAXELIMLENGTH) not for
 * fixed data allocations, but just for protection from memory exhaustion.
 * You may increase or decrease these limits if you need...
 * 
**/

/* the maximal number of generators for a pfield (pfdefs.c) */
/*#define	PF_NUMEXPS	2*/

/* maximal number of generated submatrices in ematrix_submatrices_() (ematexp.c) */
/*#define	EM_MAXSUBMAT	800000*/

/* the number of allocated matrices to keep for recycling (ematrix.c) */
/*#define	EM_ALCMAX	70000*/


/* the maximal number of line-maps generated from one matrix (strlmap.c) */
/*#define	SM_MAXLINEMAPS  400000*/

/* the maximal number of refering minors generated for one matrix (strminor.c) */
/*#define	STR_MAXMINORS	400000*/

/* the maximal size difference for which we may compute minors by removing (strminor.c) */
/*#define	STR_MAXMINREM   6*/


/* the maximal number of (pre-)generated matrix extension vectors (gmatext.c) */
/*#define	GEN_MAXMATEXTENS	1000000*/

/* the length of an elimination sequence may NOT be bigger than NUM_LONG_BITS !!! */
/*#define	GEN_MAXELIMLENGTH	(NUM_LONG_BITS-2)*/




/**
 * These defines are for handling lists in the program (void** alist), as implemented
 * in misc/alist.c.
 * The defines include control of a general-purpose "recycling" cache for previously
 * generated lists, as implemented in alreuse.c, and in particular source files.
 * 
**/

/* the size of the table of arbitrary-list ends (alist.c), and the maximal total list length */
/*#define	ALIST_HASHSZ	3507*/
/*#define	ALIST_MAXLIST	5000000l*/

/* see more "recycling" defines in misc.h for alreuse.c ... */
/*#define	AL_MAXLTOTAL	300000*/

/* "recycling" of self-map lists for matrices are enabled always (strminor.c) */

/* whether to "recycle" generated bases lists for matrices - not implemented now... (strminor.c) */
/*#define	CACHEDBASES	0*/



/**
 * These are some defines for control of the extension generating process, mainly
 * for the canonical testing.
 * They speed-up the rejection process, and they are desired...
 * 
**/

/* whether to use pfield pre-checking when generating extensions (gen.h, ematext.c) */
/*#define GEN_USEEXTPREPF	1*/

/* a speedup choice - use randomized structure checks at level 1 for at least this (9)
    number of elements (gen.h, genstep.c) */
/*#define GEN_RANDSTRUCTCH	9*/

/* a speedup choice - use randomized canonical checks at level 1 for at least this (9)
    number of elements (gen.h, genstep.c) */
/*#define GEN_RANDCANON		9*/

/* a speedup choice - use some pre-computed data for tests when filtering extensions (gener.c, genstep.c) */
/*#define GEN_USEPRECHECK	1*/

/* a speedup choice - compare (canon) incomplete subsequences when looking for a smaller one,
    when at least this less (2) from the full length (genstep.c) */
/*#define GEN_SUBSEQCOMP	2*/

/* a speedup choice - added lines of the same type may be considered unordered (genstep.c) */
/*#define GEN_SUBSEQUNORD	1*/
                







#endif





























