
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
 * This is just a continuation of the file frcoms.c...
 * (Included from there, NOT compiled separately!)
**/




/******************	Definition table for command handles	*******************/
/**********************************************************************************/


#undef CURDLEV
#define CURDLEV         4


/**
 * Extra third-party handles may be included here.......
 * 
**/

#define	FR_COMSEXTRA	0
#include "frcoms-more.inc"
#undef	FR_COMSEXTRA



/**
 * Here is the definition of the command handling array...
 * Read more about input/output parameter descriptions in ../include/frame.h,
 * including the definitions of the FCPAR_* and FCRET_* macros.
 * 
 * The definition of the whole structure commhandles is in fr.h .
 * In particular, the second field is the pattern for getting new frame name.
 * It gets the output name from the first input frame, and
 * values of index and command parameters (use %i$ conversion to change the order):
 	snprintf(buf,490,cmf.fnformat, nmret,i,
                  OPTPARAM(com,0),OPTPARAM(com,1),OPTPARAM(com,2),OPTPARAM(com,3));
 * (No second 'i' after nmret is given when printing name for FCRET_REPLMX[1L] !)
 * The default pattern FCDEFNFORMAT is in fr.h.
 * 
**/


commhandles	comhdefs[] = {
	
	{
	
		/*********	Frame i/o and handling	*********/
		/************************************************/
		
		"pr",	NULL,	{NULL},	NULL,
		"Use: \"!pr frame-list\"; \tprints the given frames.",
		&frcom_prframe,
		1, { FCPAR_LISTFR0 },
		FCRET_PRINTED
	},{
		"print",	NULL,	{NULL},	NULL,
		"Use: \"!print frame-list\"; \tprints the matrices of the given frames.",
		&frcom_print,
		1, { FCPAR_LISTFRMP0 },
		FCRET_PRINTED
	},{
		"printmore",	NULL,	{NULL},	NULL,
		"-Use: \"!printmore mat-list\"; \tprints info about the matroids in frames.",
		&frcom_printmore,
		1, { FCPAR_LISTMX0 },
		FCRET_PRINTED
	},{
		"prmore",	NULL,	{NULL},	NULL,
		"Use: \"!prmore mat-list\"; \tprints info about the matroids in frames.",
		&frcom_printmore,
		1, { FCPAR_LISTMX0 },
		FCRET_PRINTED
	},{
		"prbases",	NULL,	{ NULL, "0","0","0","bas", NULL },	NULL,
		"Use: \"!prbases mat-list [id ..]\"; \tprints all bases of the matroids in frames.",
		&frcom_printbasecirc,
		5, { FCPAR_LISTMX0, FCPAR_INTEGER,FCPAR_INTEGER,FCPAR_INTEGER, FCPAR_STRING },
		FCRET_PRINTED
	},{
		"prcircuits",	NULL,	{ NULL, "0","0","0","circ", NULL },	NULL,
		"Use: \"!prcircuits mat-list [id ..]\"; \tprints all circuits of the matroids in frames.",
		&frcom_printbasecirc,
		5, { FCPAR_LISTMX0, FCPAR_INTEGER,FCPAR_INTEGER,FCPAR_INTEGER, FCPAR_STRING },
		FCRET_PRINTED
	},{
		"prtree",	NULL,	{ "(T)", NULL},	NULL,
		"Use: \"!prtree frame-list\"; \tprints the frame subtree (not matrices).",
		&frcom_printtree,
		1, { FCPAR_LISTFR0 },
		FCRET_PRINTED
	},{
		"prtext",	NULL,	{NULL},	NULL,
		"Use: \"!prtext text\"; \tprints the given text.",
		&frcom_printtext,
		1, { FCPAR_STRING },
		FCRET_PRINTED
	},{
		"verbose",	NULL,	{ "1", NULL},	NULL,
		"Use: \"!verbose [+-n]\"; \tadjusts the current output verbosity (default up).",
		&frcom_setverbose,
		1, { FCPAR_INTEGER },
		FCRET_NOTHING
	},{
		"quiet",	NULL,	{ "-1", NULL},	NULL,
		"Use: \"!quiet [-+n]\"; \tadjusts the current output verbosity (default down).",
		&frcom_setverbose,
		1, { FCPAR_INTEGER },
		FCRET_NOTHING
	},{
		"setname",	NULL,	{NULL},	NULL,
		"Use: \"!setname name frame(s)\"; \tnames the given frames (%%d to number them).",
		&frcom_setname,
		2, { FCPAR_STRING, FCPAR_LISTFR },
		FCRET_NOTHING
	},{
		"write",	NULL,	{NULL},	NULL,
		"Use: \"!write frame(s)\"; \twrites the frame to a file (with matrices).",
		&frcom_write,
		1, { FCPAR_LISTFR },
		FCRET_PRINTED
	},{
		"writeto",	NULL,	{NULL},	NULL,
		"Use: \"!writeto fn frame\"; \twrites the frame to the file fn.",
		&frcom_write,
		2, { FCPAR_STRING, FCPAR_LISTFR },
		FCRET_PRINTED
	},{
		"writecom",	NULL,	{NULL, "1", NULL},	NULL,
		"Use: \"!writecom frame\"; \twrites the frame to a file with matroid struct comments.",
		&frcom_write,
		2, { FCPAR_LISTFR, FCPAR_INTEGER },
		FCRET_PRINTED
	},{
		"writetree",	NULL,	{ "(T)", NULL},	NULL,
		"Use: \"!writetree frame(s)\"; \twrites the frame subtree to a file (with matrices).",
		&frcom_writetree,
		1, { FCPAR_LISTFR },
		FCRET_PRINTED
	},{
		"writetreeto",	NULL,	{ NULL, "(T)", NULL},	NULL,
		"Use: \"!writetreeto fn frame\"; \twrites the frame subtree to the file fn.",
		&frcom_writetree,
		2, { FCPAR_STRING, FCPAR_LISTFR },
		FCRET_PRINTED
	},{
		"read",	NULL,	{ "NOname", NULL},	NULL,
		"Use: \"!read fname >dest\"; \treads the given file fname to a new subtree.",
		&frcom_read,
		1, { FCPAR_STRING },
		FCRET_ONE
	},{
		"mread",	"%s",	{ "NOname", NULL},	NULL,
		"Use: \"!mread fname >dest\"; \treads all matrices from the given file fname.",
		&frcom_mread,
		1, { FCPAR_STRING },
		FCRET_LIST
	},{
		"append",	NULL,	{NULL},	NULL,
		"Use: \"!append frame input\"; \tappends the given input to the frame (experimental!!!).",
		&frcom_append,
		2, { FCPAR_LISTFR, FCPAR_STRING },
		FCRET_PRINTED
	},{
	



		/*********	Simple matrix manipulation	*********/
		/********************************************************/
	
	
		"mark",		NULL,	{NULL},	NULL,
		"Use: \"!mark fr-list\"; \tmarks frames within the tree (for use with ~N,^N).",
		&frcom_marklist,
		1, { FCPAR_LISTFR0 },
		FCRET_LREMFR
	},{
		"move",		NULL,	{NULL},	"()",
		"Use: \"!move src-list [>dest-list]\"; \tmoves (copies, deletes) frames within the tree.",
		&frcom_move,
		1, { FCPAR_LISTFR0 },
		FCRET_PARTLIST
	},{
		"mmove",	"%s",	{NULL},	"()",
		"Use: \"!mmove src-list [>dest-list]\"; \tmoves (copy, del) matrices-only within the tree.",
		&frcom_mmove,
		1, { FCPAR_LISTMX0 },
		FCRET_PARTLIST
	},{
		"flatten",	NULL,	{NULL},	"()",
		"Use: \"!flatten src-list >dest-list\"; \tflattens frame-subtrees within the tree.",
		&frcom_flatten,
		1, { FCPAR_LISTFR0 },
		FCRET_LIST
	},{
		"pivot",	"%s~p%s%s",	{NULL},	NULL,
		"Use: \"!pivot row col matrix\"; \tthe pivot is applied to one given matrix.",
		&frcom_pivot,
		3, { FCPAR_INTEGER, FCPAR_INTEGER, FCPAR_ONEMX },
		FCRET_REPLMX1
	},{
		"pivot-pr",	NULL,	{NULL},	NULL,
		"-Use: \"!pivot-pr row col matrix\"; \tprinting the pivot applied to the matrix.",
		&frcom_pivot,
		3, { FCPAR_INTEGER, FCPAR_INTEGER, FCPAR_ONEMX },
		FCRET_PRINTED
	},{
		"dual",		"%s#",	{NULL},	NULL,
		"Use: \"!dual matrix(ces)\"; \ttranspose is applied to given matrix(ces).",
		&frcom_dual,
		1, { FCPAR_LISTMX },
		FCRET_REPLMXL
	},{
		"delete",	"%s~d%s",	{NULL},	NULL,
		"Use: \"!delete id matrix\"; \tthe line of given (id) is deleted from the matrix.",
		&frcom_delete,
		2, { FCPAR_INTEGER, FCPAR_ONEMX },
		FCRET_REPLMX1
	},{
		"contract",	"%s~c%s",	{NULL},	NULL,
		"Use: \"!contract id matrix\"; \tthe line of given (id) is contracted from the matrix.",
		&frcom_contract,
		2, { FCPAR_INTEGER, FCPAR_ONEMX },
		FCRET_REPLMX1
		
	},{
		"deleach",	"%s~de%d",	{"((T))", "1", NULL},	"(((0t)))",
		"Use: \"!deleach mat(s) >out-list\"; \tmakes list by deleting each element from the matrix.",
		&frcom_removeeach,
		2, { FCPAR_LISTMX, FCPAR_INTEGER },
		FCRET_LIST
	},{
		"coneach",	"%s~ce%d",	{"((T))", "0", NULL},	"(((0t)))",
		"Use: \"!coneach mat(s) >out-list\"; \tmakes list by contracting each element from the matrix.",
		&frcom_removeeach,
		2, { FCPAR_LISTMX, FCPAR_INTEGER },
		FCRET_LIST
	},{
		"remeach",	"%s~re%d",	{"((T))", "2", NULL},	"(((0t)))",
		"Use: \"!coneach mat(s) >out-list\"; \tmakes list by removing each element from the matrix.",
		&frcom_removeeach,
		2, { FCPAR_LISTMX, FCPAR_INTEGER },
		FCRET_LIST
	},{
		"pfield",	NULL,	{"", NULL},	NULL,
		"Use: \"!pfield newpf\"; \tswitches (temporarily) to the given new pfield.",
		&frcom_pfield,
		1, { FCPAR_STRING },
		FCRET_PRINTED
	},{
		"import",	"%s<",	{"Id0", NULL},	NULL,
		"Use: \"!import transl_nm frame(s)\"; \timports the frames by translation transl_nm.",
		&frcom_import,
		2, { FCPAR_STRING, FCPAR_LISTFR },
		FCRET_REPLMXL
	},{
	
	


		/*********	Structural matrix tests	*********/
		/************************************************/
	
	
		"fan",		NULL,	{NULL},	NULL,
		"Use: \"!fan mat-list\"; \tprints the max fans in the matroids.",
		&frcom_hasfan,
		1, { FCPAR_LISTFRM },
		FCRET_PRINTED
	},{
		"hasfan",	NULL,	{"-1", NULL},	NULL,
		"Use: \"!hasfan f mat-list\"; \twhether there is a fan of len f in the matroids.",
		&frcom_hasfan,
		2, { FCPAR_INTEGER, FCPAR_LISTFRM },
		FCRET_PRINTED_F
	},{
		"msize",	NULL,	{NULL, "0", "0", "==.", "0", NULL},	NULL,
		"Use: \"!msize mat-list R C cmp [S]\"; \ttest matrix sizes RxC, or elems S, cmp is \"[<=>.]+\".",
		&frcom_hasmsize,
		5, { FCPAR_LISTFR0, FCPAR_INTEGER, FCPAR_INTEGER, FCPAR_STRING, FCPAR_INTEGER },
		FCRET_PRINTED_F
	},{
		"girth",	NULL,	{NULL},	NULL,
		"Use: \"!girth mat-list\"; \tprints the girth (min cycle) of the matroids.",
		&frcom_hasgirth,
		1, { FCPAR_LISTFRM },
		FCRET_PRINTED
	},{
		"hasgirth",	NULL,	{"-1", NULL},	NULL,
		"Use: \"!hasgirth g mat-list\"; \twhether girth is >=g in the matroids.",
		&frcom_hasgirth,
		2, { FCPAR_INTEGER, FCPAR_LISTFRM },
		FCRET_PRINTED_F
	},{
		"ispaving",	NULL,	{NULL, "11111"},	NULL,
		"Use: \"!ispaving mat-list\"; \twhether the matroid(s) is paving (girth>=rank).",
		&frcom_hasgirth,
		2, { FCPAR_LISTFRM, FCPAR_INTEGER },
		FCRET_PRINTED_F
	},{
		"connectivity",	NULL,	{NULL},	NULL,
		"Use: \"!connectivity mat-list\"; \tprints connectivity of the given matroids.",
		&frcom_isconnected,
		1, { FCPAR_LISTFRM },
		FCRET_PRINTED
	},{
		"isconn",	NULL,	{NULL, "2"},	NULL,
		"Use: \"!isconn mat-list [cn]\"; \twhether the matroids are connected (connectivity >=cn).",
		&frcom_isconnected,
		2, { FCPAR_LISTFRM, FCPAR_INTEGER },
		FCRET_PRINTED_F
	},{
		"isconn3",	NULL,	{NULL, "3"},	NULL,
		"Use: \"!isconn3 mat-list\"; \twhether the given matroids are 3-connected.",
		&frcom_isconnected,
		2, { FCPAR_LISTFRM, FCPAR_INTEGER },
		FCRET_PRINTED_F
	},{
		"isconni4",	NULL,	{NULL, "3", "i"},	NULL,
		"Use: \"!isconni4 mat-list\"; \twhether the given matroids are internally 4-connected.",
		&frcom_isconnected,
		3, { FCPAR_LISTFRM, FCPAR_INTEGER, FCPAR_STRING },
		FCRET_PRINTED_F
	},{
		"issimple",	NULL,	{NULL, "1", "s"},	NULL,
		"Use: \"!issimple mat-list\"; \twhether the given matroids are simple.",
		&frcom_isconnected,
		3, { FCPAR_LISTFRM, FCPAR_INTEGER, FCPAR_STRING },
		FCRET_PRINTED_F
	},{
		"iscosimple",	NULL,	{NULL, "1", "c"},	NULL,
		"Use: \"!iscosimple mat-list\"; \twhether the given matroids are cosimple.",
		&frcom_isconnected,
		3, { FCPAR_LISTFRM, FCPAR_INTEGER, FCPAR_STRING },
		FCRET_PRINTED_F
	},{
		"bwidth3",	NULL,	{NULL},	NULL,
		"Use: \"!bwidth3 mat-list\"; \twhich of (3-conn!) matroids have branch-width 3.",
		&frcom_bwidth3test,
		1, { FCPAR_LISTFRM },
		FCRET_PRINTED_F
	},{
		"minor",	NULL,	{NULL},	NULL,
		"Use: \"!minor mat-list min-list\"; \twhich matroids have minor in min-list.",
		&frcom_minorlist,
		2, { FCPAR_LISTFRM, FCPAR_LISTMX },
		FCRET_PRINTED_F
	},{
		"minorusi",	NULL,	{NULL},	NULL,
		"Use: \"!minorusi mat-list minor id\"; \twhich mats have minors using elem id.",
		&frcom_minorusing,
		3, { FCPAR_LISTFRM, FCPAR_ONEMX, FCPAR_INTEGER },
		FCRET_PRINTED_F
	},{
		"tmajor",	NULL,	{NULL},	NULL,
		"Use: \"!tmajor mat-list min-list\"; \twhich matroids are tight majors of min-list.",
		&frcom_tmajorlist,
		2, { FCPAR_LISTFRM, FCPAR_LISTMX },
		FCRET_PRINTED_F
	},{
		"inpfield",	NULL,	{NULL},	NULL,
		"Use: \"!inpfield mat-list\"; \twhich of the matrices are properly represented.",
		&frcom_pfreprlist,
		1, { FCPAR_LISTFRM },
		FCRET_PRINTED_F
	},{
		"equiv",	NULL,	{NULL},	NULL,
		"Use: \"!equiv mat-list1 mat-list2\"; \twhich of list1 are (unlab) equiv to some in list2.",
		&frcom_equivalent,
		2, { FCPAR_LISTFRM, FCPAR_LISTMX },
		FCRET_PRINTED_F
	},{
		"eqpair",	NULL,	{NULL},	NULL,
		"Use: \"!eqpair mat-list\"; \twhich mats in the list are (unlab) pairwise equiv.",
		&frcom_equivalent,
		1, { FCPAR_LISTFRM },
		FCRET_PRINTED_F
	},{
		"unique",	NULL,	{NULL},	NULL,
		"-Use: \"!unique mat-list alt-list\"; \tmake list unique using altern. repres. from alt-list.",
		&frcom_makeunique,
		2, { FCPAR_LISTFRM, FCPAR_LISTMX },
		FCRET_PRINTED_F
	},{
		"mhash",	NULL,	{"-1", NULL},	NULL,
		"Use: \"!mhash h-value mat-list\"; \twhich of list have matroid hash-value = h-value.",
		&frcom_matroidhash,
		2, { FCPAR_INTEGER, FCPAR_LISTFR },
		FCRET_PRINTED_F
	},{
		"isomorph",	NULL,	{NULL, NULL},	NULL,
		"Use: \"!isomorph mat-list iso-list\"; \ttest matroids in mat-list for isomorphic in iso-list.",
		&frcom_isomorph,
		2, { FCPAR_LISTFRMP0, FCPAR_LISTFRMP0 },
		FCRET_PRINTED_F
	},{
		"isompair",	NULL,	{NULL },	NULL,
		"Use: \"!isompair mat-list\"; \ttest matroids in mat-list for pairwise isomorphisms.",
		&frcom_isomorph,
		1, { FCPAR_LISTFRMP0 },
		FCRET_PRINTED_F
	},{
	
	


		/*********	Matrix constructions	*********/
		/************************************************/
	
	
		"autom",	NULL,	{NULL},	NULL,	/* replaced by !selfmap */
		"-Use: \"!autom matrix\"; \tprints all automorphisms of the matrix.",
		&frcom_automorph,
		1, { FCPAR_ONEMX },
		FCRET_PRINTED
	},{
		"selfmap",	NULL,	{NULL},	NULL,
		"Use: \"!selfmap matrix\"; \tprints all line self-maps of the matrix.",
		&frcom_automorph,
		1, { FCPAR_ONEMX },
		FCRET_PRINTED
	},{
		"extend",	NULL,	{"b", "((T))", NULL},	"(((0t)))",
				/* name-pattern GEN_?EXTNAME for extensions is in gener.h */
		"Use: \"!extend [rcb]* matrix\"; \tgenerates (co)extensions of the matrix by \"rcb..\".",
		&frcom_extend,
		2, { FCPAR_STRING, FCPAR_LISTFRM },
		FCRET_LIST
	},{
		"extendsize",	NULL,	{NULL,NULL, "((T))", NULL},	"(((0t)))",
		"-Use: \"!extendsize r c matrix\"; \tgenerates (co)extensions of the matrix to size r x c.",
		&frcom_extend,
		3, { FCPAR_INTEGER, FCPAR_INTEGER, FCPAR_ONEFRM },
		FCRET_LIST
	},{
		"extendto",	NULL,	{NULL,NULL, "((T))", "==", NULL},	"(((0t)))",
		"Use: \"!extendto r c matrix\"; \tgenerates (co)extensions of the matrix to size r x c.",
		&frcom_extend,
		4, { FCPAR_INTEGER, FCPAR_INTEGER, FCPAR_ONEFRM, FCPAR_STRING },
		FCRET_LIST
	},{
		"extendtor",	NULL,	{NULL,NULL, "((T))", "=<", NULL},	"(((0t)))",
		"Use: \"!extendtor r c matrix\"; \tgenerates (co)extensions to r rows and <=c cols.",
		&frcom_extend,
		4, { FCPAR_INTEGER, FCPAR_INTEGER, FCPAR_ONEFRM, FCPAR_STRING },
		FCRET_LIST
	},{
		"extendtoc",	NULL,	{NULL,NULL, "((T))", "<=", NULL},	"(((0t)))",
		"Use: \"!extendtoc r c matrix\"; \tgenerates (co)extensions to c cols and <=r rows.",
		&frcom_extend,
		4, { FCPAR_INTEGER, FCPAR_INTEGER, FCPAR_ONEFRM, FCPAR_STRING },
		FCRET_LIST
	},{
		"extendupto",	NULL,	{NULL,NULL, "((T))", "<<", NULL},	"(((0t)))",
		"Use: \"!extendupto r c matrix\"; \tgenerates (co)extensions to <=r rows and <=c cols.",
		&frcom_extend,
		4, { FCPAR_INTEGER, FCPAR_INTEGER, FCPAR_ONEFRM, FCPAR_STRING },
		FCRET_LIST
	},{
		"repres",	NULL,	{"", NULL, NULL},	NULL,
		"Use: \"!repres pfield mat-list\"; \ttest pfield representability of the given mats.",
		&frcom_reprtest,
		2, { FCPAR_STRING, FCPAR_LISTFRMP0 },
		FCRET_PRINTED_F
	},{
		"represgen",	"%sR%d",	{"((T))", "one", NULL},	"(((0t)))",
		"Use: \"!represgen mat-list [all[q]]\"; \tgenerate representations of the given mats.",
		&frcom_reprgener,
		2, { FCPAR_LISTFRMP0, FCPAR_STRING },
		FCRET_LIST
	},{




		/*********	Command-flow Control	*********/
		/************************************************/
	
	
		"restart",	NULL,	{NULL},	NULL,
		"Use: \"!restart\"; \trestart command processing in the whole tree.",
		&frcom_restart,
		0, { 0 },
		FCRET_PRINTED	/* FCRET_NOTHING ?? */
	},{
		"exit",	NULL,	{"0", NULL},	NULL,
		"Use: \"!exit num\"; \tstop processing of all commands immediately, exit with num.",
		&frcom_comexit,
		1, { FCPAR_INTEGER },
		FCRET_PRINTED	/* FCRET_NOTHING ?? */
	},{
		"skip",	NULL,	{"1", NULL},	NULL,
		"Use: \"!skip num\"; \tskip the next num commands in the current frame.",
		&frcom_comskip,
		1, { FCPAR_INTEGER },
		FCRET_PRINTED	/* FCRET_NOTHING ?? */
	},{
		"iflist",	NULL,	{NULL},	NULL,
		"Use: \"!iflist len [<=>!] fr-list\"; \twhether the frame-list has len?? frames.",
		&frcom_iflist,
		3, { FCPAR_INTEGER, FCPAR_STRING, FCPAR_LISTFR0 },
		FCRET_PRINTED	/* FCRET_NOTHING ?? */
	},{
		"ifshell",	NULL,	{ "echo", "0", NULL },	NULL,
		"Use: \"!ifshell command [stat]\"; \ttesting return status of given shell command.",
		&frcom_ifshell,
		2, { FCPAR_STRING, FCPAR_INTEGER },
		FCRET_PRINTED	/* FCRET_NOTHING ?? */
	},{
		"iffile",	NULL,	{ "", "r", NULL },	NULL,
		"Use: \"!iffile filename\"; \ttesting readability of the given file.",
		&frcom_iffile,
		1, { FCPAR_STRING, FCPAR_STRING },
		FCRET_PRINTED	/* FCRET_NOTHING ?? */
	},{
	
	
	
	
	
/*		
		"t",		NULL,	{NULL},	NULL,
		"Use: \"!\"; \tsample...",
		&frcom_test,
		1, { FCPAR_ONEFR },
		FCRET_LIST
	},{
*/
			/*****		*****/
	

#define	FR_COMSEXTRA	1	/* (extra command handles included) */
#include "frcoms-more.inc"
#undef	FR_COMSEXTRA

		NULL,NULL,{NULL},NULL,NULL,NULL,0,{0},0	/* must end with NULL's! */
	}
};





















