
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
 * This header declares supplementary functions and values used when handling frames
 * in the program.
 * These are internal declarations only.
 * Read more in ../include/frame.h.
 * 
**/


#ifndef FR_H
#define FR_H







/******************	Internal functions of the analyzer	*******************/
/**********************************************************************************/


/**
 * These functions handle includes and subframes in the input lexical analyzer.
 * See emflex.l and emflexsup.c, also emflexb.l...
 * 
**/

char*   yy_specrepl(char *st, char *repl[]) ;
char*   yy_ovaluesubst(char *st) ;
char*   yyb_processvalue(framestruc *cfr, char *st) ;
void    yyb_appendout(char *s) ;
void    yyb_appendsubst(char *sb) ;

int     yyincl(char *fn, int fil) ;
int     yy_entrysubst(char *st, int wh, char **rep) ;
void    yy_startframe(char *name, framestruc *frst) ;
framestruc*     yy_endframe(void) ;


/**
 * These functions and constants are used to evaluate pfield expressions when reading
 * a matrix (frame) in the input lexical analyzer.
 * See emflex.l and emflexsup.c...
 * 
**/

enum YY_ENTRYOPS {
	YYOP_NOP, YYOP_EVAL, YYOP_RBRACE,
		YYOP_PLUS, YYOP_MINUS, YYOP_TIMES, YYOP_OVER,
	YYOP_POWER, YYOP_LBRACE
};
int     yy_entryop(char *nst, int op, int pow) ;
int     yy_entryevalcheck() ;
int     yy_mxentryeval(int r, int c) ;
int     yy_mxentryline(int r, int c, int ln) ;
int     yy_valentryeval(int cd) ;
void    yy_valentrylast(exp_t *x, sign_t *g) ;


/**
 * These functions read and store options and commands from the input frame.
 * See emflex.l and emflexsup.c...
 * 
**/

void    yy_optionstart(char *nm, int cd, int ln, int ct) ;
void    yy_optionval(int ix, char *val, long dv) ;
void    yy_optionend(int ix, int ln) ;
optionstruc*    yy_optionproc(void) ;
void    yy_optionprocpar(int ix, char *val) ;
void    yy_optionpfval(int ix, exp_t x, sign_t g) ;
		/* (immediate application of some options during input) */
void    yy_applyoption(optionstruc *op, int ln) ;






/******************	Option descriptions	*******************/
/******************************************************************/


/**
 * There is a list of all options used in the program, for help and input-checking
 * purposes only.
 * The list is defined in fropts.c, and it contains simple descriptions as well.
 * The options listed in optnowrite[] by names are not written to files.
**/

typedef struct option_descriptions {
	
	char		*name;		/* the name of this option (case-insensitive) */
	char		*description;
	int		numpar;		/* the number of the option values (-1 for unlimited) */
	
}	optdescript;

extern optdescript	optdescdefs[];

extern char		*optnowrite[];






/******************	Command-handling declarations	*******************/
/**************************************************************************/


/**
 * The following structure describes a command handling function and its parameters.
 * (See frame.h and frame_applycommand() in frameop.c for more details on its use.)
 * The declarations provided here are private for the command-handling code.
 * You find the master description and other declarations in ../include/frame.h.
 * Command handles themselves are defined in frcoms.c.
 * 	......
**/

#define	COMFDECL(ff)	int (*ff)(char**,void**[],int,long[],int,int,void**[],long[])
					/* (see master COMFUNCTION(ff) in frame.h) */

#define	FCMAXPARAMS	(MAXOPTPARAMS-1)	/* the maximal allowed number of parameters (<=MAXOPTPARAMS !) */

	/* this structure is used to define handling of an arbitrary command, see frcoms.c */
typedef struct command_handles {
	
	char		*name;		/* the name of this command (case-insensitive) */
	char		*fnformat;	/* the format for getting file name for the result */
	char		*defpar[FCMAXPARAMS];	/* possible default parameter values, or NULL */
	char		*defret;	/* possible default return description, or NULL */
	char		*description;
	
	COMFDECL(fhandle);		/* the function to handle this command */
	
	int		numpar;		/* the number and types of the command parameters */
	int		par[FCMAXPARAMS];
	
	int		ret;		/* the type of the command result */
	
}	commhandles;

extern commhandles	comhdefs[];	/* (defines all command handling functions) */

	/* the default frame-name format after applying a command - see frcoms.c */
#define	FCDEFNFORMAT	"%s-%d"



/**
 * These functions collect the input parameters for a command based on a text description.
 * Read about the description syntax in frame.h.
**/

ematrix**       frame_getparammx_list(framestruc *fr, int ip, int ip0, char *pr, int one, int gmark, framestruc ***remfrl) ;
framestruc**    frame_getparamfr_list(framestruc *fr, int ip, int ip0, char *pr, int one, int gmark, int needmx) ;

framestruc**    frame_recoverlres(framestruc **fl, int needmx) ;

int     frame_getparamfr_recur(framestruc *fr, framestruc ***list,
                                 int ip, char *pr, int one, int gmark, int needmx, int depth) ;
                                        
int     frame_paramcloseup(char *pr) ;


/**
 * These functions are used to erase used frames, and to store command output
 * according to a text description.
 * Read about the description syntax in frame.h.
**/

void    frame_deletemark_recur(framestruc *fr, int del) ;

framestruc**    frame_storeres_list(framestruc *fr, framestruc **oul, char *pr, int one, char *nn) ;
int	frame_storeres_recur(framestruc *fr, framestruc ***list, char *pr, int one, char *nn, int depth) ;


/**
 * One may control output printing verbosity within commands by a command and by an option...
 * 
**/

extern int	frcom_printbrief;
int     frcom_verbose_out(void) ;







/******************	Input/output declarations	*******************/
/**************************************************************************/


/**
 * Supplementary functions for frame writing...
 * 	......
**/

int     frame_write_recur(FILE *file, framestruc *fr, int rec, int opt, int inh, int onm, int depth) ;

int     frame_writecomment(FILE *fp, char *comchar, char *comment) ;


void    frame_gettree_recur(framestruc **fin, int ord, framestruc ***fout) ;


















#endif	/* (of #ifndef FR_H) */
























