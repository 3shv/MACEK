
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
 * This is the main header file for the Macek project.
 * There is nothing interesting here, just a collection of other specific
 * header files.
 * Temporary third-party includes should be put into macek-more.inc ...
**/




#ifndef	MACEK_H
#define	MACEK_H


/* the system includes */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <limits.h>
#include <unistd.h>
#include <time.h>
#include <sys/stat.h>


/* the local includes */

#include "version.h"

#include "profile.h"
#include "debug.h"
#include "sysdep.h"
#include "misc.h"

#include "pfield.h"
#include "ematrix.h"
#include "frame.h"
#include "struct.h"
#include "gener.h"

#define	DEBUG_H2
#include "debug.h"


/* temporary third-party includes for additional functions... */

#include "macek-more.inc"


#endif






















