
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
 * This file is provided for third-party additions to the extension-generating process.
 * See master description of generating in ../include/gener.h.
 * 
 * This file (possibly) adds more fields to the elimination-sequence structure elimseqm,
 * adds processing for more frame-options that control the generating process,
 * and adds more tests for the sequence and structure checks (see genstep.c).
 * 
 * This file is included from eseq.h, and also from gener.c
 * (controlled by the value of GENER_ELIMQDEF).
 * 
**/


#include "macek.h"
#include "gen.h"










/******************	Additions to elimseqm	*******************/
/******************************************************************/


#if	GENER_ELIMQDEF>0

                        /* additions to the sequence structure type: */
/* typedef struct {
	..............
*/
	
	/* additions come here... */
	
	
/*
	..............
}       elimseqm;
*/


#else









/******************	Additional generating functions	*******************/
/**************************************************************************/


#undef CURDLEV
#define	CURDLEV		6


/**
 * This is for debug-printing of the additional fields of elimseqm to the file fo.
 * See elimseqm_fprint_ext() in elimseq.c.
 * 
**/

void	elimseqm_fprint_incl(FILE *fo __attribute__ ((unused)), int ch __attribute__ ((unused)),
		elimseqm *q __attribute__ ((unused)), char *pref __attribute__ ((unused))) {
	//
}



/**
 * These functions add more tests to the sequence and structure checks at level lev.
 * Read carefully the theory description in ../include/gener.h and in genstep.c
 * before adding anything to these functions!!!
 * 
**/

int	genstep_sequencecheck_incl(int lev __attribute__ ((unused)), elimseqm *q __attribute__ ((unused)),
		ematrix *e __attribute__ ((unused)), int tr __attribute__ ((unused))) {
	//
	return 1;
}

int	genstep_structurecheck_incl(int lev __attribute__ ((unused)), elimseqm *q __attribute__ ((unused)),
		ematrix *e __attribute__ ((unused)), int tr __attribute__ ((unused))) {
	//TIMEELAPSED();
	return 1;
}


/**
 * This function may be used to pre-compute some data used later in the above checks.
 * Read in genstep.c...
 * 
**/

int	genstep_pretest_incl(elimseqm *q __attribute__ ((unused)),
		int tr __attribute__ ((unused)), int wh __attribute__ ((unused))) {
	//
	return 1;
}


/**
 * Some special tests may be added here, be careful.
 * Read in genstep.c...
**/

int	genstep_specialcheck_incl(int lev __attribute__ ((unused)), elimseqm *q __attribute__ ((unused)),
		ematrix *e __attribute__ ((unused)), int tr __attribute__ ((unused))) {
	//
	return 1;
}



/**
 * This function is called from gener_extframe() in gener.c to convert more frame-options
 * in fr to local information held in the structure elimseqm during generating,
 * and to free the data later.
 * 
 * The firs call with wh==1 should set the additional information to the sequence q
 * (returned back) acoording to fr. (No qel,frc,flist used here!)
 * The second call with wh==0 gets the generated list qel of sequences, and a new frame frc.
 * It should set possible options in frc to values copied to all output frames.
 * The last call with wh==-1 gets the whole output list of frames flist for final
 * polishing, and this call should also free all local data refered in the sequence q.
 * 
**/

elimseqm*	gener_extframe_incl(framestruc *fr __attribute__ ((unused)), int trq __attribute__ ((unused)),
				int ctrl __attribute__ ((unused)), elimseqm *q,
				elimseqm **qel __attribute__ ((unused)), framestruc **frc __attribute__ ((unused)),
				framestruc **flist __attribute__ ((unused)), int wh) {
	//
	junk = wh;
	return q;
}




#endif


























