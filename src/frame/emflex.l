
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
 * This is the lexical analyzer for the general input entity "frame" to the program.
 * Read about frames in ../include/frame.h and frame.c.
 * 
 * 
 * The scanner starts scanning of a given string (by a call to frame_doinput_ext()),
 * and processes the content plus file-includes and subframes.
 * It is internally called as frame_flex() once for the whole input (i.e. there is no return
 * call during scanning, all scanned values and properties are stored into static variables).
 * Of course, the scanner is not re-entrant.
 * In general, the given first string forms the resulting root frame (returned back)
 * together with the regular includes; and subframes or possible consequent strings form
 * frames that are sons of the root, or sons of the sons, etc.
 * So the result is a tree-structure of frames with the root frame returned back.
 * (Subframes may be given in separate files, as well as all in one file.)
 * 
 * The following information are read for a new frame:
 * Matrix entries are given (row by row) on specific matrix lines.
 * Option lines give (arbitrary) named options and their values.
 * Command lines give (arbitrary) commands with their parameters (also procedures).
 * Moreover, comment and blank lines are ignored from the input.
 * For the option/command values and include filenames, additional special replacements
 * apply, read their description in emflexb.l (a supplementary scanner implementing these).
 * 
 * It is possible to include any files to the scanned input.
 * One may also start and end subframes anywhere through the input, and even sub-input
 * subframes on one frame line, see below.
 * When reading the input, an end-of-line is appended to each input file.
 * If an error happens, then it returns NULL with an error message.
 * 
 * Read the scanner definitions below for details about how different line-types
 * are recognized, and how their content is scanned.
 * Kepp in mind that the scanner definitions here are tightly related to the
 * pfield_scanvalue..() and pfield_printvalue..() functions in ../pfield/pfield.c,
 * and to the frame-output code in emfout.c.
 * 
**/







/******************	Lex declarations	*******************/
/******************************************************************/


/**
 * The following declarations define what kind of a line is scanned.
 * For each kind of a line, a separate exclusive scanner is defined below.
 * The declaration must "eat" the whole determining prefix, but nothing more.
 * (Subframe lines are somehow special - one line may start or end a regular subframe,
 * but also several strings may be given on one line as subframes "{ x1 x2 .. }".)
 * 
**/

COMMENTLINE	^[[:blank:]]*[#%]+
MATRIXLINE	^([[:blank:]]+[^@!&<}]|[01])
MATRIXLINES	[=]
OPTIONLINE	^[[:blank:]]*([@]|[@]?[oO]ption)
COMMANDLINE	^[[:blank:]]*([!]|[!]?[cC]ommand)
COMPROCLINE	^[[:blank:]]*([&]|[&]?[Pp]rocedure)
INCLUDELINE	^([[:blank:]]*[<]|[<]?[iI]nclude)

	/* (add the begin-lines '^' separately when switching to the subframe condition) */
SUBFRBEGIN	[/>]?(SUBFR(AME)?|[sS]ubfr(ame)?|"{")
SUBFREND	[/>]?(EOFR(AME)?|[eE]ofr(ame)?|[[:blank:]]*"}")

/**
 * Special additional rules applies to the first input string (and to line-subframes):
 * Any line starting with a number is a matrix line, any line starting with
 * a letter is an include line (controlled by the value of inpspecial).
 * Moreover, all ',' are replaced by space, all ';' are replaced by eol,
 * which happens already before the scanning starts using SPECREPLACE.
**/

SPECMATRIXLINE	^[[:digit:]]
SPECINCLUDELINE	^[[:alpha:]./~]
%{
#define	SPECREPLACE	{ ",;", " \n" }
%}

/**
 * If an unrecognized character appear on the input, we may try to do something
 * else in special cases identified here...
 * (The expression should catch a whole line.)
 * The first expression works only in the INITIAL state.
 * 
**/

ERRORCATCH_A	([^[:blank:]$]*[$][^$[:blank:]][^\n]*)

/**
 * These declarations provide separators for matrix entries on a matrix line,
 * and for option values / command parameters on an option / command line.
 * (The NO-version of each is supplied for matching non-separated pieces of the input,
 * notice that "\n" must be included there.)
**/

MATRIXSEP	[[:blank:],"]
MATRIXNOSEP	[^[:blank:]\n,"]
/*MATRIXNOSEP	[^[:blank:]\n,()\[\]]*/

OPTIONNAME	[[:alpha:]_][[:alnum:]_.()+-]*
 /* (see also the declaration of "substitutes" for option/command values in emflexb.l) */
OPTIONSEP	[[:blank:],]
OPTIONNOSEP	[^[:blank:]\n,]
%{
char	*STROPTIONSEP = " \t,\"";	/* what characters require the option value to be quoted */
%}

/**
 * These declarations provide separators for include lines and subframes.
 * (A subframe-flags may be alone on separate lines, but it is also possible to start
 * a full subframe line that contains subframe strings and the subframe end, even
 * nested subframes...)
**/

INCLUDESEP	[[:blank:],]
INCLUDENOSEP	[^[:blank:]\n,]

SUBFRSEP	[[:blank:],]
SUBFRNOSEP	[^[:blank:]\n,}{]

/**
 * This expression matches option names which are supposed to have a pfield-valued first
 * parameter (which is handled separately in Ioptionpfval).
**/

OPTIONVALPF	([rR]equire)


/**
 * These declarations control how one pfield entry is scanned.
 * In general, an entry may be an arbitrary algebraic expression, including bracketing,
 * which is defined in the pfield (using + - * / ^).
 * Usual precedence and bracketing rules apply.
 * Moreover, the values 0 or 1 may be written using special chracters SPECZERO or SPECONE,
 * but only when there is nothing else in the entry.
 * The multiplication operator * may be skipped if it is clear where the atoms end.
 * 
 * Keep in mind that symbolic atoms are tried first for the longest match (but at most
 * 10 characters!), only then the scanner tries to decode an expression.
 * (I.e. "1/a-1" may be declared as a symbolic atom.)
 * NO floating numbers are allowed at all - see that '.' is optionally used for multiplication.
 * 
**/

OPENBRACE	[(\[]
CLOSEBRACE	[)\]]
STARTNUMBER	[0-9]
STARTSYMBOL	[[:alpha:]1-9]
SPECZERO	[0Oo]
SPECONE		[1Il]
OPPLUS		[+]
OPMINUS		[-]
OPMUL		[*.]
OPDIV		[/:]
OPPOWER		[&^]


/**
 * These are exclusive conditions for the line-type subscanners.
 * A special condition Ierror is provided for "breaking" a subscanner when an error occurs.
 * 
**/

%x	Icomment
%x	Imxnumbers0 Ivalnumbers0 Imxnumbers Ivalnumbers
%x	Ioption Ioptionval Ioptionpf Ioptionpfval
%x	Icommand Icomproc Ioptionvalnosub
%x	Iinclude Isubframe Isubframeend
%x	Ierror Ierrorcatch

%option stack
/*%option debug*/
%array

/**
 * We use local variables inpfile, inpline, and errhint when reporting errors.
 * The variables bracelev, exprix, mxline, mxcolumn are used when scanning matrix
 * entries, they keep the bracketing level, number of atom in an entry, row and
 * column of the matrix.
 * The variable fullineno counts the logical lines (except comments) across all included
 * files, and it is used later to define the precedence of scanned commands and options.
 * The variable framelev keeps track of the nested level of subframes (it is used just
 * to catch unbalanced subframe marks).
 * The variable optvalix counts the values for an option or a command (excluding the name).
 * The variable inpspecial controls whether special input-lines apply to the input.
 * The variable include_stack_ptr is taken from yywrap(), it handles file inclusion.
**/

%{
static char	*inpfile = NULL;
static int	inpline = 1, errhint = -1;
static int	bracelev = 0, exprix = 0;
static int	mxline = 0, mxcolumn = 0;
static int	fullineno = 1;
static int	framelev = 0, lineframe = 0;
static int	optvalix = 0;
static int	inpspecial = 0;
static int	include_stack_ptr = -1;             


#include "macek.h"
#include "fr.h"

#define YY_DECL		int	frame_flex(void)

/**
 * A special error reporting scheme is used here - it prints where the bad input happened...
 * After any input error, the scanner returns whenever possible, and the result of
 * scanning is null instead of a frame.
 * A modification BEGINLINE of the macro BEGIN is provided to return from an error at
 * the beginning of the next line.
**/

static int	lerror_occured = 0;
#define LERROR(args...) (USERERRORWH(inpfile,inpline,args),lerror_occured=-1)

#define	BEGINLINE(cd)	if (lerror_occured<0) return lerror_occured; else BEGIN(cd)
%}

%%











 /******************	Lexical scanner		*******************/
 /*****************************************************************/


%{
#define CURDLEV		6
%}

 /**
  * Here the lexical analysis starts.
  * In the initial condition (or inclusive), we decide what type of line we are going to read.
  * (Blank lines and '\'-escaped spaces are silently ignored.)
  * We also count logical lines of the input.
  * 
  * Special additional rules applies to the first input string (inpspecial):
  * Any line starting with a number is a matrix line, any line starting with
  * a letter is an include line.
 **/


<*>"\\"([[:blank:]])+			/* '\'-escaped blanks are ALWAYS silently discarded */
<*>"\\"([[:blank:]])*"\n"	++inpline;

[[:blank:]]*"\n"	{		/* empty lines are not scanned (unless in an excl condition) */
		++inpline;  fullineno++;
	}
{COMMENTLINE}	{			/* start of a comment line */
		BEGINLINE(Icomment);
	}
{MATRIXLINE}	{			/* start of a line of matrix entries: ("pfnumbers") */
		yyless(0);		/* (starting number must be scanned again) */
		BEGINLINE(Imxnumbers0);
	}
^{MATRIXLINES}	{
		BEGINLINE(Imxnumbers0);
	}
{OPTIONLINE}	{			/* start of an option line */
		BEGINLINE(Ioption);
	}
{COMMANDLINE}	{			/* start of a command line */
		BEGINLINE(Icommand);
	}
{COMPROCLINE}	{			/* start of a procedure line (works like an include of additional commands) */
		BEGINLINE(Icomproc);
	}
{INCLUDELINE}	{			/* start of an include line */	
		BEGINLINE(Iinclude);
	}
 /* (the subframe lines are handled differently in the include section below) */
^{SUBFRBEGIN}		{
		yyless(0);  BEGINLINE(Isubframe);
	}
^{SUBFREND}		{
		yyless(0);  BEGINLINE(Isubframeend);
	}

{SPECMATRIXLINE}	{		/* special start of a matrix line in the input string */
		if (!inpspecial)  REJECT;
		yyless(0);
		BEGINLINE(Imxnumbers0);
	}
{SPECINCLUDELINE}	{		/* special start of an include line in the input string */
		if (!inpspecial)  REJECT;
		yyless(0);
		BEGINLINE(Iinclude);
	}


	/**
	 * We simply discard all comments up to (including) the newline.
	 * (To make logical long lines, use '\n' before eol, not comment...)
	**/
<Icomment>{
[^\n\\]*	DEBUG(CURDLEV+2,"Comment: \"%s\".\n",yytext);
"\n"		++inpline;  BEGIN(INITIAL);
}









 /******************	Options and Commands	*******************/
 /*****************************************************************/


	/**
	 * Here we read options - first the name of the option, then the list of values.
	 * The option line ends with an eol, unless the eol is escaped with '\'.
	 * Option values are separated by {OPTIONSEP} (it is possible to enclose the value in "").
	 * There may be arbitrarily many values present.
	 * If the value looks numeric, then it is also converted to a number.
	 * If the option name matches something in {OPTIONVALPF}, then the first value(s)
	 * is supposed to be a pfield expression, and so it is given also to the
	 * number-subscanner for reading, see below at <Ioptionpfval>.
	 * 
	 * The same code is used also to scan commands and their parameters,
	 * since there is no distinction between these two from the lex's point of view.
	**/
<Ioption,Ioptionpf,Icommand,Icomproc>{
[[:blank:]]+				/* ignored blanks before the option name */

{OPTIONNAME}	{			/* the option name is read here: */
		int	stc = YY_START;
		if (stc==Ioptionpf) {
			/* (if requested a pfield argument, then the calling condition is on the stack) */
			yy_pop_state();  stc = YY_START;
			BEGIN(Ioptionpfval);
		} else {		/* switching to the values-reading condition (common for all) */
			/*** BEGIN(stc!=Icomproc?Ioptionval:Ioptionvalnosub);	not used now ***/
			BEGIN(Ioptionval);
		}
		DEBUG(CURDLEV+2,"Option/command %s \"%s\":\n",YY_START==Ioptionpfval?"(pf)":"",yytext);
					/* creating the new option/command with the given name */
		yy_optionstart(yytext,stc,fullineno,0);
		optvalix = 0;		/* to count the number of following values/parameters */
	}
}
<Ioptionval,Ioptionpfval,Ioptionvalnosub>{
{OPTIONSEP}+				/* (ignored value-separators) */

{OPTIONNOSEP}+			|	/* reading values:   one-word strings */
["]([^"]*)("\\\""[^"]*)*["]	{	/*	""-enclosed strings, possible \" inside */
		long	d,k;
		char	*sv0,*sv=NULL;
		if (YY_START!=Ioptionvalnosub)	/* substituting possible macros in the value string */
			sv0 = yy_ovaluesubst(yytext);
		else  sv0 = MSTRDUP(yytext);
		if (sv0)  sv = yy_specrepl(sv0,NULL);
		DEBUG(CURDLEV+3,"  \'%s\' -> \'%s\'\n",yytext,sv);
		if (!sv)  break;
		k = sscanf(sv,"%ld",&d);	/* (trying to read as an integer) */
		if (k>0)  SDEBUG(CURDLEV+2," (%ld)",d);
		else  d = OPTNOPARAM;
					/* recording one value for the current option/command */
		yy_optionval(optvalix,sv,d);
		++optvalix;
		if (YY_START==Ioptionpfval) {
			yyless(0);	/* (to read as a pfield entry as well...) */
			yy_push_state(Ivalnumbers0);
		}
		FREE(sv);  FREE(sv0);
	}
}
<Ioptionval,Ioptionvalnosub>"\n" {	/* ending the option line */
		optionstruc	*proc;
		int		i;
		++inpline;  fullineno++;
		BEGIN(INITIAL);		/* (begin init first, so that possible procedure is included from newline) */
		proc = yy_optionproc();
		if (proc!=NULL) {	/* if this is a "&procedure" call, then it translates to an include: */
			for (i=0; i<MAXOPTPARAMS; i++)
				yy_optionprocpar(i,OPTPARAM(proc,i));
			yyincl(OPTNAME(proc),1);
		}
		yy_optionend(optvalix,fullineno-1);
	}

	/**
	 * This part of the option subscanner controls reading of the first option value
	 * as a pfield expression.
	 * (Since the patterns here are longer than those above, they are matched first.)
	 * At the beginning, we identify the special option name(s) with condition Ioptionpf,
	 * and return the name back for regular scanning.
	 * Later, after a value is read above in the Ioptionpf condition, it is pushed back
	 * for further scanning there, and caught in Ivalnumbers0.
	 * Then the option separator(s) (after the value scanning ends) is caught here, and the
	 * pfield value is stored to the option.
	 * 
	**/
<Ioption,Icommand>{
{OPTIONVALPF}{OPTIONSEP}	{	/* initial recognition of a pfield-valued option */
		yyless(0);
		yy_push_state(Ioptionpf);
	}
}
<Ioptionpfval>{
{OPTIONSEP}*({OPTIONNOSEP}|"\n")	{	/* this is called after the pfield scan ends */
		exp_t	xx;	sign_t	gg;
		if (optvalix>0) {	/* getting the last scanned pfield value */
			yy_valentrylast(&xx,&gg);
			SDEBUG(CURDLEV+2," =%s",pfield_pvalue(14,xx,gg));
			yy_optionpfval(optvalix-1,xx,gg);	/* (optvalix was already incremented above) */
		}
		if (optvalix>=1 || yytext[yyleng-1]=='\n') {
			yyless(0);
			BEGIN(Ioptionval);	/* when no more pfield expressions are expected */
		} else {
			REJECT; }	/* rejects the next option for a regular scan */
	}
}








 /******************	Includes and Subframes	*******************/
 /*****************************************************************/

%{
#undef CURDLEV
#define CURDLEV		6
%}

	/**
	 * This part handles subframe and end-of-subframe lines.
	 * One may start and end each suframe on separate lines, but it is also possible
	 * to use the following shortcut:
	 *  All space-separated strings (use "for longer") on the subframe-line are
	 *  included as separate subframes, even deeply nested.
	 * ...
	**/
<Isubframe,Isubframeend>{
"\n"		{		/* nothing special when ending a line - subframes may continue on next lines */
		++inpline;  fullineno++;
		BEGIN(INITIAL); 
	}
{SUBFRSEP}+				/* ignored separators */
}

<Isubframe,Isubframeend>{SUBFRBEGIN}	{
		++framelev;
		lineframe = 0;
		DEBUG(CURDLEV+1,"Input subframe (lev %d) ...\n",framelev);
		yy_startframe(NULL,NULL);
	}
<Isubframe,Isubframeend>{SUBFREND}	{
		if (framelev<=0)  {LERROR("Too many ending subframes!"); break;}
		DEBUG(CURDLEV+1,"Input subframe end (lev %d).\n",framelev);
		yy_endframe();
		--framelev;
	}

<Isubframe>{
["]([^"]*)("\\\""[^"]*)*["]	|	/* (""-enclosed strings, possible \" inside) */
{SUBFRNOSEP}+			{	/* non-space line-subframe strings are read here */
		char	*sv, *repl[2] = SPECREPLACE;
		sv = yy_specrepl(yytext,repl);
		yyincl(sv,0);
		FREE(sv);
		if (lineframe++>0) {	/* (separating subframes on the input frame-line) */
			yy_endframe();  yy_startframe(NULL,NULL);
		}
		DEBUG(CURDLEV+1,"Input subframe middle-separation %d (lev %d).\n",lineframe,framelev);
	}
}


	/**
	 * This part handles include lines (possibly with multiple files).
	 * Each file sequentially scanned from the line is included into the input stream.
	 * (Not creating a subframe!)
	**/
<Iinclude>{
{INCLUDESEP}*				/* ignored separators */
"\n"		{			/* ending the include line */
		++inpline;  fullineno++;  BEGIN(INITIAL);
	}
["]([^"]*)("\\\""[^"]*)*["]	|	/* ""-enclosed strings, possible \" inside */
{INCLUDENOSEP}+			{		/* non-space filenames */
		char	*sv0,*sv=NULL;
		sv0 = yy_ovaluesubst(yytext);	/* substituting possible macros in the include-name */
		if (sv0)  sv = yy_specrepl(sv0,NULL);
		DEBUG(CURDLEV-2,"Including from \'%s\' (->\'%s\')\n",yytext,sv);
		if (!sv)  break;
		yyincl(sv,1);
		FREE(sv);  FREE(sv0);
	}
}









 /******************	Input matrix lines	*******************/
 /*****************************************************************/

%{
#undef CURDLEV
#define CURDLEV		5	/*5*/
%}

	/**
	 * Here we read matrix number lines - {MATRIXSEP}-separated (space) matrix entries.
	 * The same sub-scanner (Ivalnumbers) is also used to read an option-values that are
	 * in the form of pfield entries.
	 * (First of all, blanks are skipped, and some internal variables are set.)
	 * 
	 * In general, an entry may be an arbitrary algebraic expression, including bracketing,
	 * which is defined in the pfield (using + - * / ^).
	 * Usual precedence and bracketing rules apply.
	 * Moreover, the values 0 or 1 may be written using special chracters SPECZERO or SPECONE,
	 * but only when there is nothing else in the entry.
	 * If no operation separates two atoms, then * is substituted.
	 * If there is no number preceeding an operation, then 0 (+ -) or 1 (* /) is added.
	 * 
	 * Keep in mind that symbolic atoms are tried first for the longest match (but at most
	 * 10 characters!), only then the scanner tries to decode an expression.
	 * (I.e. "1/a-1" may be declared as a symbolic atom.)
	 * NO floating numbers are allowed at all - see that '.' is used for multiplication.
	 * 
	 * ........
	**/
<Imxnumbers0,Ivalnumbers0>{		/* starting - skip blanks, and initialize */
{MATRIXSEP}*{MATRIXNOSEP}	{
		yyless(yyleng-1);
		mxcolumn = exprix = bracelev = 0;
		errhint = -1;
		BEGIN(YY_START==Imxnumbers0? Imxnumbers:Ivalnumbers);
	}
{MATRIXSEP}*	{		/* no entries to scan on the line - return */
		if (YY_START==Imxnumbers0)  BEGIN(INITIAL);
		else { yyless(0);  yy_pop_state(); }
	}
}
<Imxnumbers,Ivalnumbers>{
		/**
		 * We first deal with separators - spaces ending matrix entries,
		 * and eol's ending whole matrix lines.
		**/
{MATRIXSEP}*"\n"	|		/* here are blanks that end the previous entry */
{MATRIXSEP}+		{
		if (exprix<=0) { REJECT; }	/* (if not ending a meaningful expression) */
		if (bracelev>0)  {LERROR("Unbalanced bracketing in the matrix-entry expression (at row %d col %d)!",mxline+1,mxcolumn+1);}
		else if (yy_entryevalcheck(YY_START)<0 && YY_START==Imxnumbers)
			{LERROR("Cannot evaluate the matrix-entry expression over %s (at row %d col %d)!",pfield_curname(),mxline+1,mxcolumn+1);}
		exprix = bracelev = 0;
		if (lerror_occured<0)  return lerror_occured;
		
		if (YY_START==Imxnumbers) {	/* using the evaluated entry in the matrix */
			yy_mxentryeval(mxline,mxcolumn);
			++mxcolumn;			/* (newline must be processed separately later) */
			if (yytext[yyleng-1]=='\n') { REJECT; }
		} else {		/* (different evaluation for Ivalnumbers - used elsewhere) */
			yy_valentryeval(YY_START);
			yyless(0);		/* (give chance to process the separators again) */
			yy_pop_state();		/* returns back to the calling-state */
		}
	}
{MATRIXSEP}*"\n"	{		/* here the eol ends one matrix line (and begins 0 cond) */
		++inpline;
		if (YY_START==Imxnumbers) {
			yy_mxentryline(mxline,mxcolumn,fullineno);
			++mxline;
		} else  REJECT;		/* (no other condition than Imxnumbers should get here) */
		fullineno++;
		BEGIN(INITIAL);
	}
{MATRIXSEP}+		{		/* this happens when an entry ends without a value */
		 LERROR("Something wrong; meaningless entry? (at row %d col %d)",mxline+1,mxcolumn+1);
	}

		/**
		 * We allow bracketing of a symbolic expression.
		 * This has precedence over all other operations.
		 * (Notice that bracketing here does not apply to exponents...)
		**/
{OPENBRACE}		{
		++bracelev;
		yy_entryop(NULL,YYOP_LBRACE,0);
	}
{CLOSEBRACE}		{
		if (--bracelev<0)  {LERROR("Unbalanced bracketing in the matrix-entry expression (at row %d col %d)!",mxline+1,mxcolumn+1); break;}
		yy_entryop(NULL,YYOP_RBRACE,0);
	}

		/**
		 * This part scans numeric entries, inluding special representations of 0,1.
		 * (The special 0 or 1 expression must be the whole entry, possibly with +-!)
		**/
{SPECONE}({MATRIXSEP}|"\n")	|	/* special symbol "one" */
{SPECZERO}({MATRIXSEP}|"\n")	{	/* special symbol "zero" */
		int	k;
		if (exprix>0) { REJECT; }
		++exprix;
		k = yy_entryop(yytext,YYOP_NOP,0);
		if (k<yyleng-1) { REJECT; }
		else  yyless(yyleng-1);
		DEBUG(CURDLEV+3,"\tFound ZERO or ONE \"%s\"\n",yytext);
	}
{STARTNUMBER}[[:digit:]]*	{	/* scanning a numeric string - see symbols below */
		int	k;
		++exprix;
		k = yy_entryop(yytext,YYOP_NOP,0);
		if (k<=0)  {LERROR("An invalid number \"%s\" was scanned.",yytext); break;}
		if (k<yyleng)  {LERROR("Not the whole number \"%s\" was scanned, why?",yytext); break;}
		DEBUG(CURDLEV+3,"\tFound number (%d) \"%s\"\n",exprix,yytext);
	}

		/**
		 * This is the most tricky part - scanning a symbolic generator for the pfield.
		 * A whole substring (of at most 10 characters!) starting with {STARTSYMBOL}
		 * is tried for a longest prefix match in yy_entryop().
		 * If a match is found, then the rest of the string is returned back,
		 * otherwise the whole string is rejected (and the rejection is remembered
		 * in yyst, not to send the same string for scanning again).
		 * 
		 * It is OK that this rule is after scanning a numeric string, since a symbol
		 * should (despite starting with a digit) be longer than a number.
		 * Moreover, this rule allows for replacement of substrings by given
		 * substitutes, which happens before scanning continues.
		 * (Replacement works only in included files, not in the top string.)
		 * One may prevent such replacements by prefixing the symbol with '_'.
		**/
[_]		
[_]?{STARTSYMBOL}{MATRIXNOSEP}{0,9}	{
		int		i=0,k;
		char		*rep=NULL;
		static char	yyst[52] = "";	/* (used to prevent repeated rejections of the same input) */
		static int	lmx=-1, lmn=0;	/* (used to prevent cycled entry replacements) */
		
		if (include_stack_ptr==0) { k = 0;  errhint = 5; }
		else  k = yy_entrysubst(yytext,0,&rep);
				/* this code prevents infinite loops when replacing entries */
		if (lmx!=mxline*mxcolumn) { lmx = mxline*mxcolumn;  lmn = 0; }
		if (k>0) if (++lmn>=50) {
			if (lerror_occured>=0)
				{LERROR("Probably cycled entry replacement \"%s\" in %s\n",rep,yytext);}
			k = 0;
		}
		if (k>0) {	/* a replacement for a string prefix was requested */
			DEBUG(CURDLEV+2,"\t\tEntry replacement by \"%s\" %d\n",rep,k);
			yyless(k);
			/*****i = YY_START;  yyincl(rep,0);  BEGIN(i);	not good - adds newline */
			for (i=strlen(rep)-1; i>=0; i--)  unput(rep[i]);
			break;
		}
		++exprix;	/* scan only if not the same input was scanned previously */
		if (strncmp(yytext,yyst,yyleng<50?yyleng:50)==0)  k = 0;
		else  if (lerror_occured>=0) {
			i = (yytext[0]=='_'? 1:0);
			k = yy_entryop(yytext+i,YYOP_NOP,0);
		}
		if (k<=0) {	/* if no match, rejects, for other rules to apply */
			strncpy(yyst,yytext,50);
			REJECT;
		} else {	/* if matched a prefix, returns the rest back for scanning */
			DEBUG(CURDLEV+3,"\tFound symbol (%d) \"%.*s\"\n",exprix,k,yytext+i);
			yyless(k+i);
		}
	}

		/**
		 * This part reads operation symbols from the entry expression.
		**/
{OPPLUS}	{
		DEBUG(CURDLEV+3,"\t\tFound PLUS\n");
		yy_entryop(NULL,YYOP_PLUS,0);
	}
{OPMINUS}	{
		DEBUG(CURDLEV+3,"\t\tFound MINUS\n");
		yy_entryop(NULL,YYOP_MINUS,0);
	}
{OPMUL}	{
		if (exprix==0) { REJECT; }
		DEBUG(CURDLEV+3,"\t\tFound MUL\n");
		yy_entryop(NULL,YYOP_TIMES,0);
	}
{OPDIV}	{
		if (exprix==0) { REJECT; }
		DEBUG(CURDLEV+3,"\t\tFound DIV\n");
		yy_entryop(NULL,YYOP_OVER,0);
	}
		/**
		 * The power ^ operation is somehow special - it may accept only integer
		 * exponents, and so the exponent is scanned already here as an int.
		 * (Optional brackets may appear around the exponent.)
		**/
{OPPOWER}[(\[][+-]?[[:digit:]]+[)\]]	|
{OPPOWER}[+-]?[[:digit:]]+		{
		int	a,d;
		if (exprix==0) { REJECT; }
		a = sscanf(yytext,"%*[ ([^]%d",&d);
		if (a<=0) { REJECT; }
		if (abs(d)>PF_MAXEXPONENT/2)  {LERROR("Too large exponent on input %d",d); d=0;}
		DEBUG(CURDLEV+3,"\t\tFound POWER ^ %d\n",d);
		yy_entryop(NULL,YYOP_POWER,d);
	}
}







 /******************	Syntax error handling	*******************/
 /*****************************************************************/


	/**
	 * This is error-catching for all unmatched input, with optional hints.
	 * (In certain context, we may guess what type of error will soon happen,
	 * but the code is still rejected for normal scanning above.)
	**/

^({STARTSYMBOL}|[+-])				errhint = 0;  REJECT;
<Ioptionval>["][^"]+"\n"			errhint = 6;  BEGIN(Ierror);  yyless(0);
<Imxnumbers,Ivalnumbers>{
[-+]?[[:digit:]]+"."[[:digit:]]+		errhint = 1;  REJECT;
[[:space:]]+({OPMUL}|{OPDIV}|{OPPOWER})		errhint = 2;  REJECT;
{OPPOWER}[-+]?[^[:digit:]]+			errhint = 3;  REJECT;
{OPPOWER}({OPMUL}|{OPDIV}|{OPPOWER})					|
({OPPLUS}|{OPMINUS})({OPPLUS}|{OPMINUS}|{OPMUL}|{OPDIV}|{OPPOWER})	|
({OPMUL}|{OPDIV})({OPPLUS}|{OPMINUS}|{OPMUL}|{OPDIV}|{OPPOWER})	{
		errhint = 4;  BEGIN(Ierror);  yyless(0);	}
}
<Ierrorcatch>{ERRORCATCH_A}	{
		char	*sv,*sv1, *repl[2] = SPECREPLACE;
		
		sv = yy_ovaluesubst(yytext);   /* trying macro substitution to correct the error */
		if (sv?!strcmp(sv,yytext):1)  { FREE(sv);  REJECT; }
		BEGIN(INITIAL);
		DEBUG(CURDLEV-3,"Syntax error catching - try \'%s\' -> \'%s\'\n",yytext,sv);
		sv1 = yy_specrepl(sv,repl);
		yyincl(sv1,0);		/* possibly changed line is tried again recursively */
		FREE(sv1);  FREE(sv);
	}
<*>.		{
		static int	lhint = -1, lline = -1;
		static char	*errhints[] = {
					"(Start a matrix line with blanks...?)",
					"(Floating numbers are not allowed...?)",
					"(Starting entry by an operator...?)",
					"(Non-integral ^exponent...?)",
					"(Invalid combination of operators...?)",
					"(Replacement is not working on the input line...?)",
					"(Unclosed double-quotas...?)"
				};
		if (YY_START==INITIAL) {	/* try some more error-catching before reporting */
			DEBUG(CURDLEV-2,"Try syntax error catching after unknown char \'%s\'\n",yytext);
			yyless(0);  BEGIN(Ierrorcatch);  break;
		}
		/****if (lerror_occured<0)  return lerror_occured; ???	*****/
		if (errhint!=lhint || fullineno!=lline)
			{LERROR("Invalid input at \"%s\".\n\t%s",yytext,errhint>=0?errhints[errhint]:"");}
		lhint = errhint;  lline = fullineno;
	}


%%














/* See emflexsup.c for supporting scanner functions... */
/* (must be included here to have flex declarations available) */

#include "emflexsup.c"
























