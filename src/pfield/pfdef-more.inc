
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
 * This file is a continuation of the definition of the pfdefins[] array from
 * pfdefs.c, for third-party additions...
 * Read ../include/pfield.h.
 * 
 * When adding a new partial field here, you may need to increase the value of
 * PF_NUMEXPS if your pfield has more generators.
 * Write your new definition of PF_NUMEXPS to include/macek-more.inc .
 * However, use values of PF_NUMEXPS>2 only when it is really neccessary,
 * and compile the extended executable separately, to keep access to the default version.
**/










/**
 * First, read the partial field theory in ../include/pfield.h, and the declaration of 
 * pf_definition in pf.h.
 * This file contains the initialized array pf_definition that defines partial fields.
 * 
 * A quick overview of the definition format is here:
 * 
#if PF_NUMEXPS>=0
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
		{ "short_strings_for_gen^i", ... },	(optional, like "2","4", etc)
		...
	},{
		{ "short_strings_for_gen^-i", ... },	(optional, like "1/2","1/4", etc)
		...
	},
	format_to_use_opt_strings,  format_not_use_brackets
},{
#endif
 * 
 * There are some hard allocation limits on the above subarrays, see in pf.h.
 * (Like PF_NUMEXPS which should be always checked, or PF_NUMNAMES, PF_NUMFUNDAMENTALS, etc.)
**/




























