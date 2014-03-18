
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
 * A "frame" is the basic general data entity in the program.
 * Read more about a frame and its data in ../include/frame.h.
 * 
 * This file contains a simple description of the options recognized by the program.
 * (If an option is not described here, it is rejected on input!)
 * 
**/





#include "macek.h"  
#include "fr.h"








/******************	Option descriptions	*******************/
/******************************************************************/


//#undef CURDLEV
//#define CURDLEV         4



/**
 * See the declaration of optdescript in fr.h.
 * The first value is the name of the option (case-insensitive), 
 * the second one is a short description for line-help,
 * and the third is the number of values accepted by this option, -1 for arbitrary.
 * 
 * Each option must have an entry here to be accepted by the program!
 * (At least, an exact prefix of the option must be here.)
**/

optdescript	optdescdefs[] = {
	
	{
		"replace",	
		"Use: \"@replace XX repl_expr\"; \tall next XX in expressions are text-replaced",
		2
	},{
		"repl-",	
		"Use: \"@repl-PF XX repl_expr\"; \tlike @replace for specific pfield PF only",
		2
	},{
		"require",	
		"Use: \"@require expression [.01]\"; \trequires the expression defined (and zero/nonzero)",
		2
	},{
		"transpose",	
		"Use: \"@transpose\"; \t(immediately) transposes the frame matrix",
		0
	},{
		"name",	
		"Use: \"@name string\"; \tgives the name to the current frame",
		1
	},{
		"comment",	
		"Use: \"@comment string\"; \tstores the comment with the current frame",
		1
	},{
		"nopfcheck",	
		"Use: \"@nopfcheck\"; \tdisables the quick pfield check for input matrices",
		0
	},{
		"inputpf",	
		"Use: \"@inputpf PF\"; \tinput matrices of all next subframes in the pfield PF",
		1
	},{
		"finherit",	
		"Use: \"@finherit optname optname ...\"; \tinherits these options when file-writing",
		-1
	},{
		"finheritall",	
		"Use: \"@finheritall optname ...\"; \tinherits all values of the options when file-writing",
		-1
	},{
		"extinherit",	
		"Use: \"@extinherit optname ...\"; \tinherits these opts when extending (bsize, signature always)",
		-1
	},{
		"extinheritall",
		"Use: \"@extinheritall optname ...\"; \tinherits all these options when extending a matrix",
		-1
	},{
		FR_ERASEOPT,
		"Use: \"@"FR_ERASEOPT" optname ...\"; \terases the last previous option(s) optname value",
		-1
	},{
		FR_ERASEALLOPT,
		"Use: \"@"FR_ERASEALLOPT" optname ...\"; \terases all previous option(s) optname values",
		-1
	},{
		"prbrief",
		"-Use: \"@prbrief [+-n]\"; \tbrief printing of command output (-1 default)",
		-1
	},{
		FR_SUBSTPREFIX,
		"Use: \"@"FR_SUBSTPREFIX"XX xx\"; \tfurther $XX (in com, opt, incl) will be substituted with xx",
		1
	},{
		FR_SUBSTDEFPREFIX,
		"Use: \"@"FR_SUBSTDEFPREFIX"XX xx\"; \t\"default\" value for substitutions",
		1
	},{
		
		GEN_OPTPREFIX "bsize",
		"Use: \"@"GEN_OPTPREFIX"bsize r c\"; \tthe base-minor size when generating extensions",
		2
	},{
		GEN_OPTPREFIX "signature",
		"Use: \"@"GEN_OPTPREFIX"signature sig\"; \tthe elim-sequence signature when gener. extensions",
		1
	},{
		GEN_OPTPREFIX "forbid",
		"Use: \"@"GEN_OPTPREFIX"forbid fname ...\"; \tforbidden minors when generating extensions",
		-1
	},{
		GEN_OPTPREFIX "tight",
		"Use: \"@"GEN_OPTPREFIX"tight fname ...\"; \t\"tight\" minors when generating extensions",
		-1
	},{
		GEN_OPTPREFIX "nofan",
		"Use: \"@"GEN_OPTPREFIX"nofan fsize\"; \tthe elim-sequence must not contain a fan >=f",
		1
	},{
		GEN_OPTPREFIX "connect",
		"-Use: \"@"GEN_OPTPREFIX"connect c\"; \telim-sequence connectivity 0-3 (use with caution!)",
		1
	},{
		GEN_OPTPREFIX "simple",
		"Use: \"@"GEN_OPTPREFIX"simple\"; \tthe elim-sequence has all simple steps (also disconn.)",
		0
	},{
		GEN_OPTPREFIX "cosimple",
		"Use: \"@"GEN_OPTPREFIX"cosimple\"; \tthe elim-sequence has all cosimple steps (disconn.)",
		0
	},{
		GEN_OPTPREFIX "connected",
		"Use: \"@"GEN_OPTPREFIX"connected\"; \tthe elim-sequence has all connected steps",
		0
	},{
		GEN_OPTPREFIX "3connected",
		"Use: \"@"GEN_OPTPREFIX"3connected\"; \tthe elim-sequence has all 3-connected steps (default)",
		0
	},{
		GEN_OPTPREFIX,
		"Use: \"@"GEN_OPTPREFIX"option value(s)\"; \tgeneric option for the elim-sequence",
		-1
	},{
		
/*		"oname",	{"defpars", ...}
		"Use: \"@\"; \tsample...",
		-1
	},{*/

	/**
	 * Possible third-party additions to the option list may be added here...
	**/
#include "fropts-more.inc"
	//},{
		NULL,NULL,0
	}
};


/**
 * 
 * These are option names (prefixes...) that should not be written to files with their frames:
 * 
**/

char	*optnowrite[] = {
	"sub-par",
	"transpose",	/* transpose must never be written since it has already affected the matrix */
	"inputpf",	/* input pfield is written before the matrix when different... */
	"repl",
	"require",	/* the replace and require options already affected the matrix when reading it */
	
	NULL	/* (must end with NULL !) */
};



































