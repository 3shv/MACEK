
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
 * This header contains the system-dependent defines and checks...
 * 
**/


/******************	Checking data types	*******************/


	/* we need 1-byte type "byte1" for correct pointer arithmetics... */
/*#if CHAR_BIT != 8	- is it true that sizeof(char) is independent of the number of bits? */
#if UCHAR_MAX != 255
#error "The program expects 8-bit char type !"
#endif
typedef	char	byte1;		/* (see also MEMALIGNMOD below...) */
/*typedef	int8_t	byte1;*/
#ifndef CHAR_BIT
#define	CHAR_BIT	8
#endif

	/* integer and long types - we need them of certain size to store all bits of our computation */
#if INT_MAX < 2147483647
#error "The program needs >=32-bit integer type!"
#endif
#define NUM_INT_BITS	((int)sizeof(int)*CHAR_BIT)
#define NUM_LONG_BITS	((int)sizeof(long)*CHAR_BIT)		/* (some program limits are derived from this) */


	/* memory alignment for our pointer arithmetics */
#define	MEMALIGNMOD		(sizeof(long))
#define	MEMALIGNSIZE(l)		(((l)/MEMALIGNMOD+1)*MEMALIGNMOD)
#define	MEMALIGNPOINT(p)	((byte1*)(p)+(MEMALIGNMOD-(long)(p)%MEMALIGNMOD))


	/* to close input files after flex is finished with them (why is this not automatic?) */
#define CLOSE_AFTER_FLEX	1




/******************	System-dependent settings	*******************/


	/* not all systems have snprintf(), so we have to emulate it... */
#if (defined __osf__)
extern char	snprintf_buf[];
#define snprintf(st,sz,ar...)	\
		(junk=sprintf(snprintf_buf,ar),strncpy(st,snprintf_buf,sz),(st)[(sz)-1]=0,junk)
#endif

	/* solaris does not have fdatasync()? - not really used anyway */
#if (defined sun)
#define	fdatasync(d)	0
#endif




/******************	Own memory allocations	*******************/


#define	MALLOC(sz)	malloc(sz)
#define	CALLOC(n,sz)	calloc(n,sz)
#define	FREE(pt)	free(pt)

	/* this macro normally allocates memory, but stops the program with an error if no mem is available */
extern void	*mmmm;
#define	MMALLOC(sz)	(((mmmm=malloc((sz)))!=NULL)?mmmm: (void*)(PROGERROREXIT("cannot allocate size %d",(int)(sz)),NULL))
#define	MREALLOC(pt,sz)	(((mmmm=realloc(pt,(sz)))!=NULL)?mmmm: (void*)(PROGERROREXIT("cannot reallocate size %d for %p",(int)(sz),pt),NULL))


#define	MEMCPY(d,s,n)	memcpy(d,s,n)
#define	MEMSET(d,c,n)	memset(d,c,n)

#define MSTRDUP(s)	strdup(s)




/******************	Random numbers	*******************/


#define	RANDOM()	random()
#define	SRANDOM(in)	srandom(in)


#define	TIME()		(int)(time(NULL))
#define	TIME000		(int)(time(NULL)%999)


#define HASHPOINTER(p,mod)   ((unsigned)(((long)(p)+((long)(p)*(long)(p))/333)&0x7FFFFF)%(mod))
















