
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
 * This is a supplementary lexical analyzer for the main one in emflex.l.
 * It purpose is to implement text-replacements in the option/command values and includes.
 * Read about frames in ../include/frame.h and frame.c.
 * 
 * The text macros are replaced in the following way:
 * If an input like $mac or ${mac} appears elsewhere than on a matrix line,
 * then it is replaced with the first word (use quotas for longer text with spaces)
 * following the latest option @samp{@@sub-macro}.                            
 * If no such option is found, then @samp{@@subd-macro} is taken.
 * If even this is not found, then the replacement text is empty.
 * 
**/







/******************	Lex declarations	*******************/
/******************************************************************/


/**
 * The following declarations define what inputs are passed unchanged to the output string,
 * what "escaped" active characters are passed as well,
 * and what are the macros to be substituted.
 * 
**/

OVALESCAPED	("\\"[$]|[$][$])
OVALPASS	([^$\\]+|"$"[^[:alpha:]_$]|"\\"[^$])
OVALUESUB	"$"[[:alpha:]_][[:alnum:]_]*
OVALUESUBB	"${"[^{}]*"}"


/**
 * These are some options for the generated scanner...
 * (We need a different scanner prefix from emflex.l!)
**/

/*%option stack*/
/*%option debug*/
%option noyywrap
%array
%option prefix="yyb"


/**
 * Here we name the scanner's main function.
 * ...
**/

%{

#include "macek.h"
#include "fr.h"

#define YY_DECL		int	frame_flexoval(void)

%}

%%






 /******************	Lexical sub-scanner	*******************/
 /*****************************************************************/


%{
#define CURDLEV		6
%}

 /**
  * Here the lexical analysis starts.
  * There is not much to do so far, most text is simply copied to the output string.
  * Macros "$mac" are substituted with the value of "@sub-mac".
  * 
  * **********
  * more should be allowed for macros, like expansion of multiple-values.........
 **/


{OVALPASS}	{
		yyb_appendout(yytext);
	}
{OVALESCAPED}	{
		yyb_appendout(yytext+1);
	}
{OVALUESUB}	{
		DEBUG(CURDLEV+1,"getting sub-stitution for \'%s\'.\n",yytext+1);
		yyb_appendsubst(yytext+1);
	}
{OVALUESUBB}	{
		yytext[yyleng-1] = 0;
		DEBUG(CURDLEV+1,"getting sub-stitution for \'%s\'.\n",yytext+2);
		yyb_appendsubst(yytext+2);
	}

<*>.		{
		if (yyleng>1)  {PROGERROR("Invalid character \'%s\' in option-value, should not happen!",yytext);}
	}		/* (there can be trailing \, silently ignored...) */



%%




/**
 * This supplementary scanner is called from yy_ovaluesubst() in emflexsup.c to process
 * macro substitutions in the option/command values or include names.
 * The call to yyb_processvalue() gives the current frame (saved here in cframe), and the
 * string st to process (converted here into an input buffer).
 * The output string is collected in the variable outs, and then returned back.
 * The function yyb_appendsubst() makes the substitution to happen, by reading
 * appropriate "sub-" option values from the frame cframe.
 * 
**/

static framestruc	*cframe;
static char		*outs;
static int		outl;


char*	yyb_processvalue(framestruc *cfr, char *st) {
	char	*s;
	
	if (!st)  return NULL;
	cframe = cfr;
	outs = MMALLOC((outl = 5000+50*strlen(st)));
	outs[0] = 0;
	yyb_scan_string(st);
	frame_flexoval();	/* this calls the above lexical analyzer preprocessed by flex */
	yyb_delete_buffer(YY_CURRENT_BUFFER);
	s = MSTRDUP(outs);  FREE(outs);
	return s;
	//******* extract multiple-values from one substitution??? - return list of strings...
}


void	yyb_appendout(char *s) {
	
	if (s? ((int)(strlen(s)+strlen(outs))>outl):1)  return;
	strcat(outs,s);
}

void	yyb_appendsubst(char *sb) {
	char	buf[110], **repl,**repd;
	
	strncpy(buf,FR_SUBSTPREFIX,30);  buf[30] = 0;
	strncat(buf,sb,60);  buf[99] = 0;
	repl = frame_getoptionval_last(cframe,buf);
				/* deafult substitution value for no or empty string */
	if (repl?(repl[0]?!repl[0][0]:1):1) {
		strncpy(buf,FR_SUBSTDEFPREFIX,30);  buf[30] = 0;
		strncat(buf,sb,60);  buf[99] = 0;
		repd = frame_getoptionval_last(cframe,buf);
		if (repd && repl)  dispose_alist(repl);
		if (repd)  repl = repd;
	}
	if (!repl)  DEBUG(CURDLEV-2,"Empty substitution value for \'$%s\'.\n",sb);
	if (repl)  yyb_appendout(repl[0]);
	//******** handle multiple-values from one substitution......... ???
	if (repl)  dispose_alist(repl);
}




































