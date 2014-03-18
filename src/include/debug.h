
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
 * This header contains the debugging defines and messages.....
 * 
 * This file contains some debug-related functions and variables for the MACEK program.
 * These functions are not used for MACEK computation, but for its debug-testing.
 * Among the defines, there are ?DEBUG?() macros for printing debug messages,
 * which can be controlled via CURDLEV, and printlev.
 * 
 * In general, the Macek project contains lots of debugging code right inside its
 * functions (this can be switched off with the define FASTPROG).
 * Such inline debug-code is usually executed only randomly, not to slow down
 * computation too much.
 * 
**/




#ifndef	DEBUG_H
#define	DEBUG_H








/******************	"Debug" definitions	*******************/



#ifdef ERROREXITCRASH
#define	CRASH			(fsync(fileno(debugout)),fsync(fileno(errorout)),mmmm=NULL,*(int*)mmmm=0)
#else
#define	CRASH			0
#endif

#define	DUPLERROR(out)		(errorout!=(out) && (!isatty(fileno(out)) || !isatty(fileno(errorout))))
#define	DERRORPRINT(args...)	(fprintf(errorout,args),DUPLERROR(debugout)?fprintf(debugout,args):0,\
				 (printout!=debugout && DUPLERROR(printout))?fprintf(printout,args):0)

#define PROGERROR(args...)	(sprintf(debugoutbuf,"\n****** Program ERROR [bad:-(]: %.25s, in %.25s() line %d ******\n  ",__FILE__,__FUNCTION__,__LINE__),\
				 snprintf(debugoutbuf+strlen(debugoutbuf),5000,args),\
				 sprintf(debugoutbuf+strlen(debugoutbuf),"\n================================================================\n"),\
				 DERRORPRINT(debugoutbuf),perror_occured=(crash_all_errors?CRASH:-1))
#define PROGERROREXIT(args...)	(PROGERROR(args),junk=CRASH,exit(-1),-1)

#define USERERROR(args...)	(sprintf(debugoutbuf,"*** ERROR: (reported by %.25s, in %.25s() l.%d) ***\n  ",__FILE__,__FUNCTION__,__LINE__),\
				 snprintf(debugoutbuf+strlen(debugoutbuf),5000,args),sprintf(debugoutbuf+strlen(debugoutbuf),"\n"),\
				 DERRORPRINT(debugoutbuf),uerror_occured=-1)
#define USERERROREXIT(args...)	(USERERROR(args),exit(-1),-1)

#define USERERRORWH(wh,wl,args...)	(sprintf(debugoutbuf,"*** ERROR (by %.12s, in %.15s() %d) in \"%.20s\" l.%d: ***\n  ",\
					 __FILE__,__FUNCTION__,__LINE__,(wh),(wl)),\
					 snprintf(debugoutbuf+strlen(debugoutbuf),5000,args),sprintf(debugoutbuf+strlen(debugoutbuf),"\n"),\
					 DERRORPRINT(debugoutbuf),uerror_occured=-1)


extern int	timeprecout, timeprecdiv;
#define	OUTPUT(args...)		(timeprecout>0?fprintf(printout,"~%*.*d~\t",timeprecout,timeprecout,\
					(TIME()/timeprecdiv)%(timeprecout<=3?1000:10000)):\
					fprintf(printout,"   m~\t"),SOUTPUT(args))
#define	SOUTPUT(args...)	(fprintf(printout,args))

/* (other output defines are contained in their specific headers...) */






/******************	Debug messages	*******************/


extern int	printlev;	/* the level for printing messages (adjusted by '-gN') */
#define	TESTDLEV(l)		(((l)<=DEBUGLEV && (l)<=printlev) || (l)<=0)

#define	DEBUGLEV	6

#define	IFRANDDEBUG(f)		(DEBUGLEV>=CURDLEV || (DEBUGLEV>=CURDLEV-3 && RANDOM()%(f)==1))
#define	IFRANDDEBUGLESS(f)	(DEBUGLEV>=CURDLEV+2 || (DEBUGLEV>=CURDLEV-3 && RANDOM()%(f)==1))

#define	DEBUG(l,args...)	(junk = TESTDLEV(l)? \
					(fprintf(debugout,"%*s[%-8.8s:%17.17s()%-4d~%3.3d] ",\
					((l)>2?(l)-2:(l)),((l)>2?"":"#"),\
					__FILE__,__FUNCTION__,__LINE__,TIME000),\
					fprintf(debugout,args),0):0)
#define	SDEBUG(l,args...)	(junk = TESTDLEV(l)? (fprintf(debugout,args),0):0)

/* (other debug messages are contained in their specific headers...) */







/******************	Fast program re-definitions	*******************/


/* uncomment the following to skip (almost) all debug tests and messages in the program */
/*#define	FASTPROG	*/


#ifdef	FASTPROG
#undef	DEBUGLEV
#define DEBUGLEV	-1
#undef	CRASH
#define CRASH		0

#undef	DEBUG
#define DEBUG(l,args...)	junk=0
#undef	SDEBUG
#define SDEBUG(l,args...)	junk=0

#undef	IFRANDDEBUG
#define	IFRANDDEBUG(f)		0
#undef	IFRANDDEBUGLESS
#define	IFRANDDEBUGLESS(f)	0
#endif


#endif	/* (of #ifndef DEBUG_H) */


#ifdef	DEBUG_H2


#ifdef	FASTPROG
#undef	EMATDEBUG
#define EMATDEBUG(l,e,pf)	junk=0
#undef	EMATDEBUGS
#define EMATDEBUGS(l,e,pf)	junk=0
#undef	FRDEBUG
#define FRDEBUG(l,f,pf)		junk=0
#undef	FRDEBUGS
#define FRDEBUGS(l,f,pf)	junk=0
#undef	ELIMQDEBUG
#define ELIMQDEBUG(l,q,pf)	junk=0
#undef	ELIMQDEBUGS
#define ELIMQDEBUGS(l,q,pf)	junk=0
#endif


#endif	/* (of #ifdef DEBUG_H2) */














