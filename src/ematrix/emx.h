
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
 * 
 * generating submatrices.....
 * 
**/

#ifndef EM_MAXSUBMAT		/* to avoid memory exhaustion, we stop generating after more submatrices */
#define	EM_MAXSUBMAT	800000
#endif

ematrix**	ematrix_submatrices_rows(ematrix *e, ematrix *rr, ematrix **outl, int sq, int nz,
				int rnd, int tot, int mrk, int mins, int maxs, int kr, int kc) ;















































