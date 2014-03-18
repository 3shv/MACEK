
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













#include "macek.h"
#include "pf.h"


/******************	Partial field definitions	*******************/
/**************************************************************************/


/**
 * First, read the partial field theory in ../include/pfield.h, and the declaration of 
 * pf_definition in pf.h.
 * This file contains the initialized array pf_definition that defines partial fields.
 * 
 * A quick overview of the definition format is here:
 * 
	{ "partial_field_name 1", "...name 2", ... }	(all names are equivalent, but put the formal one first)
	"short_description",
	characteristic,  number_of_exponents,  is_partial_field,
	{ modulo_of_exponent, ... },		(use 0 for no limit on exponent - Z subgroup)
	
	{ minus_one_value },			(exponential representation of -1, if not the same)
	{
		{ fundamental_element },	{ fundamental_element-1 },
		...
	},
	{ numeric_of_generators, ... },		(must give some numbers at least for hash and input)
	{ modulos_of_generators, ... },		(modulo by which the numeric generator is used, optional)
	
	{ "strings_for_generators", ... },	(if not given spec string at [0], then numeric input is used)
	{
		{ "short_strings_for_gen^i", ... },	(optional, like "2","4", etc, followed by negatives)
		...
	},{
		{ "short_strings_for_gen^-i", ... },	(optional, like "1/2","1/4", etc)
		...
	},
	format_to_use_opt_strings,  format_not_use_brackets,  flag_no_atomic_input
 * 
 * Format levels are 1 for file saving, 0 for screen printing, and -1 (loose) for debug messages.
 * 
 * Several important notes about the strings representing generators and their informal powers:
 *  The generator string is always scanned from the input, and so we suggest to enclose
 *  it in square brackets if it looks like an expression -- for example "[a-1]".
 *  If an expression-atom is used anywhere, it should be correctly parsed as an expression.
 *  The optional strings may be strange symbols (like "ai" for a-1), or also expressions.
 *  However, optional informal strings are not scanned from the input if they contain
 *  characters from PFNOATOMINPUT (like "+-*^/") and flag_no_atomic_input>0 is given.
 * 
 * There are some hard allocation limits on the above subarrays, see in pf.h.
 * (Like PF_NUMEXPS which should be always checked, or PF_NUMNAMES, PF_NUMFUNDAMENTALS, etc.)
 * 
 * The initialization of pfield definitions follows:
**/


pf_definition	pfdefins[] = {
	{
	
	/**
	 * This is the GF(2) field - i.e. computing modulo 2.
	**/
#	if PF_NUMEXPS>=0
		{ "GF(2)", "GF2", "binary" },
		"GF(2): computing modulo 2.",
		2,	0,	0,
		{ 1 },
		{ 1, {{}} },
		{
			{ 1, {{}} },	{ 0, {{}} },
		},
		{ 1 },	{ 2 },
		{ NULL },
		{
		},{
		},
		-1,  1,  0
	},{
#	endif
	
	/**
	 * This is the REGular pfield - i.e. computing with only -1,0,1.
	 * The only fundamental element is 1 (0 is not counted as fundamental).
	**/
#	if PF_NUMEXPS>=0 && !defined(BINARYONLY)
		{ "regular", "reg" },
		"Regular: partial computing with -1,0,1.",
		-1,	0,	1,
		{ 1 },
		{ -1, {{}} },
		{
			{ 1, {{}} },	{ 0, {{}} },
		},
		{ 1 },	{ 2 },
		{ NULL },
		{
		},{
		},
		-1,  1,  0
	},{
#	endif
	
	/**
	 * This is the GF(3) field - i.e. computing 0,1,2 mod3.
	**/
#	if PF_NUMEXPS>=1
		{ "GF(3)", "GF3", "ternary" },
		"GF(3): computing 0,1,2 mod3.",
		3,	1,	0,
		{ 2 },
		{ 1, {{1}} },
		{
			{ 1, {{0}} },	{ 0, {{0}} },
			{ 1, {{1}} },	{ 1, {{0}} },
		},
		{ 2 },	{ 3 },
		{ NULL },
		{
		},{
		},
		-1,  1,  0
	},{
#	endif
	
	/**
	 * This is the GF(5) field - i.e. computing 0,1,2,3,4 mod5.
	**/
#	if PF_NUMEXPS>=1
		{ "GF(5)", "GF5" },
		"GF(5): computing 0,1,2,3,4 mod5.",
		5,	1,	0,
		{ 4 },
		{ 1, {{2}} },
		{
			{ 1, {{0}} },	{ 0, {{0}} },
			{ 1, {{1}} },	{ 1, {{0}} },
			{ 1, {{2}} },	{ 1, {{3}} },
			{ 1, {{3}} },	{ 1, {{1}} },
		},
		{ 2 },	{ 5 },
		{ NULL },
		{
		},{
		},
		-1,  1,  0
	},{
#	endif
	
	/**
	 * This is the GF(7) field - i.e. computing 0,1,2,3,4,5,6 mod7.
	**/
#	if PF_NUMEXPS>=1
		{ "GF(7)", "GF7" },
		"GF(7): computing 0,1,..,6 mod7.",
		7,	1,	0,
		{ 6 },
		{ 1, {{3}} },
		{
			{ 1, {{0}} },	{ 0, {{0}} },
			{ 1, {{1}} },	{ 1, {{2}} },
			{ 1, {{2}} },	{ 1, {{0}} },
			{ 1, {{3}} },	{ 1, {{5}} },
			{ 1, {{4}} },	{ 1, {{1}} },
			{ 1, {{5}} },	{ 1, {{4}} },
		},
		{ 3 },	{ 7 },
		{ NULL },
		{
		},{
		},
		-1,  1,  0
	},{
#	endif
	
	/**
	 * This is the GF(4) field - i.e. computing 0,1,w,w+1.
	**/
#	if PF_NUMEXPS>=1
		{ "GF(4)", "GF4", "quaternary" },
		"GF(4): computing 0,1,w, where 1+1=0 & w^2=w+1.",
		2,	1,	0,
		{ 3 },
		{ 1, {{0}} },
		{
			{ 1, {{0}} },	{ 0, {{0}} },
			{ 1, {{1}} },	{ 1, {{2}} },
			{ 1, {{2}} },	{ 1, {{1}} },
		},
		{ 2 },	{ 0 },
		{ "w" },
		{
			{ "w" }, { "w+1" }
		},{
			{ "w+1" }, { "w" }
		},
		1,  1,  1	/* (optional strings are printed, but not scanned) */
	},{
#	endif
	
	/**
	 * This is the GF(9) field - i.e. computing 0,1,w,w+1.
	**/
#	if PF_NUMEXPS>=1
		{ "GF(9)", "GF9" },
		"GF(9): computing 0,1,w, where 3=0 & w^2=w+1.",
		3,	1,	0,
		{ 8 },
		{ 1, {{4}} },
		{
			{ 1, {{0}} },	{ 0, {{0}} },
			{ 1, {{1}} },	{ 1, {{7}} },
			{ 1, {{2}} },	{ 1, {{1}} },
			{ 1, {{3}} },	{ 1, {{5}} },
			{ 1, {{4}} },	{ 1, {{0}} },
			{ 1, {{5}} },	{ 1, {{6}} },
			{ 1, {{6}} },	{ 1, {{3}} },
			{ 1, {{7}} },	{ 1, {{2}} },
		},
		{ 3 },	{ 0 },
		{ "w" },
		{
			{ "w" }, { "w+1" }, { "1-w" }, { "2" },
			{ "-w" }, { "2-w" }, { "w+2" }
		},{
			//{ "w+1" }, { "w" }
		},
		1,  1,  1	/* (optional strings are printed, but not scanned) */
	},{
#	endif
	
	/**
	 * This is the GF(8) field - i.e. computing 0,1,w,w+1.
	**/
#	if PF_NUMEXPS>=1
		{ "GF(8)", "GF8" },
		"GF(8): computing 0,1,w, where 1+1=0 & w^3=w+1.",
		2,	1,	0,
		{ 7 },
		{ 1, {{0}} },
		{
			{ 1, {{0}} },	{ 0, {{0}} },
			{ 1, {{1}} },	{ 1, {{3}} },
			{ 1, {{2}} },	{ 1, {{6}} },
			{ 1, {{3}} },	{ 1, {{1}} },
			{ 1, {{4}} },	{ 1, {{5}} },
			{ 1, {{5}} },	{ 1, {{4}} },
			{ 1, {{6}} },	{ 1, {{2}} },
		},
		{ 2 },	{ 0 },
		{ "w" },
		{
			{ "w" }, { NULL }, { "w+1" }
		},{
			//{ "w+1" }, { "w" }
		},
		1,  1,  1	/* (optional strings are printed, but not scanned) */
	},{
#	endif
	
	/**
	 * The "dyadic" partial field consists of +- all integral powers of 2 and 0.
	 * This pfield has four fundamental elements  1, 2, 1/2, -1.
	**/
#	if PF_NUMEXPS>=1
		{ "dyadic" },
		"Dyadic: partial computing with 0 and all integer powers of 2.",
		-1,	1,	1,
		{ 0 },
		{ -1, {{0}} },
		{
			{ 1,  {{0}} },	{ 0,  {{0}} },
			{ 1,  {{1}} },	{ 1,  {{0}} },
			{ 1,  {{-1}} },	{ -1, {{-1}} },
			{ -1, {{0}} },	{ -1, {{1}} },
		},
		{ 2 },	{ 0 },
		{ "2" },
		{
			{ "2" }, { "4" }, { "8" }
		},{
			{ "1/2" }, { "1/4" }, { "1/8" }
		},
		1,  1,  1	/* (optional strings are printed, but not scanned) */
	},{
#	endif
	
	/**
	 * The "near-regular" partial field consists of all expressions of the form
	 *  +- (a)^m * (a-1)^n and 0  (distinguished positive and negative values).
	 * This pfield has seven fundamental elements  1, a, 1-a, (a-1)/a, 1/a,
	 * 1/(1-a), a/(a-1).
	**/
#	if PF_NUMEXPS>=2
		{ "near-reg", "near-regular", "nreg" },
		"Near-regular: partial computing with 0 and symbolic a^*.(a-1)^*.",
		-1,	2,	1,
		{ 0, 0 },
		{ -1, {{0,0}} },
		{
			{ 1,  {{0,0}} },	{ 0,  {{0,0}} },
			{ 1,  {{1,0}} },	{ 1,  {{0,1}} },
			{ -1, {{0,1}} },	{ -1, {{1,0}} },
			{ 1,  {{-1,1}} },	{ -1, {{-1,0}} },
			{ 1,  {{-1,0}} },	{ -1, {{-1,1}} },
			{ -1, {{0,-1}} },	{ -1, {{1,-1}} },
			{ 1,  {{1,-1}} },	{ 1,  {{0,-1}} },
		},
		{ 5, 4 },	{ 0, 0 },
		{ "a", "[a-1]",		NULL, "1-a" },
		{
			{ "a", "a-1",		NULL }
		},{
			{ "1/a", "1/a-1",	NULL }
		},
		-1,  -1,  1	/* (no optional strings plus bracketing always, except -1 form) */
	},{
#	endif
	
	/**
	 * The "golden mean" partial field consists of all expressions of the form
	 *  +- (a)^m and 0  (distinguished positive and negative values) where a^2=a+1.
	 * This pfield has seven fundamental elements  1, a, -a, a^2, 1/a, -1/a, a^-2.
	**/
#	if PF_NUMEXPS>=1
		{ "gmean", "golden-mean", "goldm" },
		"Golden mean: partial computing with 0 and +-a^* where a^2=a+1.",
		-1,	1,	1,
		{ 0 },
		{ -1, {{0}} },
		{
			{ 1,  {{0}} },	{ 0,  {{0}} },
			{ 1,  {{1}} },	{ 1,  {{-1}} },
			{ -1, {{1}} },	{ -1, {{2}} },
			{ 1,  {{2}} },	{ 1,  {{1}} },
			{ 1,  {{-1}} },	{ -1, {{-2}} },
			{ -1, {{-1}} },	{ -1, {{1}} },
			{ 1,  {{-2}} },	{ -1, {{-1}} },
		},
		{ 3 },	{ 0, 0 },
		{ "a" },
		{
			{ "a" },			{ "a+1" }
		},{
			{ "a-1",	"1-a" }		/*, { "2-a" }*/
		},
		0,  -1,  1	/* (optional strings are allowed, but they are bracketed) */
	},{
#	endif
	
	/**
	 * The "6-root-of-unity" partial field consists of all expressions of the form
	 * 0 and a^* where a^6=1.
	 * This pfield has three fundamental elements  1, a, a^5=1/a.
	**/
#	if PF_NUMEXPS>=1
		{ "sqrt-6", "root6", "6root", "6sqrt" },
		"6th-root-of-unity: partial computing with 0 and a^*, where a^6=1 (a^2=a-1).",
		-1,	1,	1,
		{ 6 },
		{ 1, {{3}} },
		{
			{ 1,  {{0}} },	{ 0,  {{0}} },
			{ 1,  {{1}} },	{ 1,  {{2}} },
			{ 1, {{5}} },	{ 1, {{4}} },
		},
		{ 2 },	{ 9 },
		{ "a" },
		{
			{ "a" }, { "a-1" }, { "-1" },{ "-a" }, { "1-a" }
		},{
		},
		1,  0,  1	/* (optional strings are printed, formal brackets, but not scanned) */
	},{
#	endif
	
	
	/**
	 * Sample entry... (copy and fill with your values)
	**/
#	if PF_NUMEXPS>=9999
		{ "name1", "name2" },
		"description",
		ch,	x,	0,
		{ 0 },
		{ -1, {{}} },
		{
			{ 1, {{0}} },	{ 0, {{0}} },
			{ f, {{}} },	{ f-1, {{}} },
		},
		{ g },	{ 0 },
		{ "gstr" },
		{
			{ "g" }
		},{
			{ "1/g" }
		},
		0, 0, 0
	},{
#	endif

	/**
	 * User-side include...
	 * 
	 * 
	**/
#include "pfdef-more.inc"
	
		{ 0 },0, 0,0,0, {}, {0,{{}}}, {{0,{{}}}}, {},{}, {},{{}},{{}}, 0,0,0
	}
};	/* (of pf_definition pfdefs[]) */













/******************	Partial field translations	*******************/
/**************************************************************************/


/**
 * First, read the partial field theory in ../include/pfield.h, and the declaration of 
 * pf_translation in pf.h.
 * This file contains the array pf_translation that defines partial field translations.
 * 
 * A quick overview of the definition format is here:
 * 
	{ "partial_field_name 1", "...name 2", ... }	(only prefix is compared, so use "" for any)
	{ "translation_name 1", "...name 2", ... }	(all names are equivalent, but put the formal one first)
	"short_description",
	number_from_exponents,		(how many exponents are translated from the numbers)
	{				(images of the generators, indexed by exponents)
		{ image_of_generator1 },
		  ...
	}
 * 
**/

pf_translation	pftrans[] = {
	{

	/**
	 * This is the default identical translation preserving the first 0 exponents.
	 * (Can be used in any current pfield.)
	**/
#	if PF_NUMEXPS>=0
		{ "" },
		{ "Id0", "nothing0" },
		"Default translation 0: first 0 exponents identical.",
		0,
		{
		},
	},{
#	endif

	/**
	 * This is the default identical translation preserving the first 1 exponent.
	 * (Can be used in any current pfield.)
	**/
#	if PF_NUMEXPS>=1
		{ "" },
		{ "Id1", "nothing1" },
		"Default translation 1: first 1 exponent identical.",
		1,
		{
			{ 1, {{1}} },
		},
	},{
#	endif

	/**
	 * This is the default identical translation preserving the first 2 exponents.
	 * (Can be used in any current pfield.)
	**/
#	if PF_NUMEXPS>=2
		{ "" },
		{ "Id2", "nothing2" },
		"Default translation 2: first 2 exponents identical.",
		2,
		{
			{ 1, {{1}} },
			{ 1, {{0,1}} },
		},
	},{
#	endif

	/**
	 * These are translations from near-regular to other pfields.
	 * The first one translates to odd fields and to all others for which 1+1
	 * is defined nonzero.
	 * The second one translates to GF(4) and other pfields with x^2=x-1.
	 * 
	**/
#	if PF_NUMEXPS>=2
		{ "GF(3)", "GF(5)", "GF(7)", "dyadic" },
		{ "Nreg-tr" },
		"General translation from near-regular to odd fields.",
		2,
		{
			{ 2, {{0}} },
			{ 1, {{0}} },
		},
	},{
		{ "GF(4)", "GF(8)", "sqrt-6" },
		{ "Nreg-tr4" },
		"Translation from near-regular to GF(4) or pfields with x^2=x-1.",
		2,
		{
			{ 1, {{1}} },
			{ 1, {{2}} },
		},
	},{
#	endif

	/**
	 * Sample entry... (copy and fill with your values)
	**/
#	if PF_NUMEXPS>=9999
		{ "pf_name1", "pf_name2" },
		{ "tr_name1", "tr_name2" },
		"description",
		0,
		{
			{ 1, {{0}} },
			{ 1, {{0}} },
		},
	},{
#	endif

	/**
	 * User-side include...
	 * 
	 * 
	**/
#include "pftran-more.inc"

		{NULL}, {NULL}, NULL, 0, {{0,{{}}}}
	}
};	/* (of pf_translation pftrans[]) */

































