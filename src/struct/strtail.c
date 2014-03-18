
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
 * Some theory for start... about 3-separations, fans, branch-width, and others:
 * (Read general description of "structural" functions on the top of ../include/struct.h.)
 * 
 * 
 * ....... not used now ..........
 * 
 * 
 * 
**/






#include "macek.h"
#include "str.h"






#if 0


/******************	Checking for tails	*******************/
/******************************************************************/


#undef CURDLEV
#define CURDLEV		4


/**
 * This function computes the full closure of the given submatrix re in the matrix e.
 * A full closure is defined by a sequence of closures and coclosures applied to the submatrix,
 * until this process stops.
 * The full-closure of re is returned in a new refering matrix to e.
 * 
 * ........ maybe store also the sequence of matrices as the closures are generated ???
**/

ematrix*	struct_fclosure_ext(int ch, ematrix *e, ematrix *re) {
	ematrix		*ec1,*ec2;
	int		i,om;
	
	if (ISREFMAT(e) || !ISREFMAT(re) || REFMAT(re)!=e)  {PROGERROREXIT("call only for a refering matrix re -> e");}
	if (ISTRANSPM(e)!=ISTRANSPM(re))  {PROGERROR("call only for matrices e,re of the same transpose");}
	DEBUG(CURDLEV+1,"looking for a full closure of re: %p %dx%d  in e: %p %dx%d\n",re,ROWSM(re),COLSM(re),e,ROWSM(e),COLSM(e));
	
	i = om = 0;
	ec1 = re;		/* a simple loop switching between closures and co-closures */
	while (i++<2 || om<ROWSM(ec1)+COLSM(ec1)) {
		ec2 = ematrix_closure_ext(e,ec1,i%2);
		om = ROWSM(ec1)+COLSM(ec1);
		if (ec1!=re)  dispose_ematrix(ec1);
		ec1 = ec2;
	}
#ifndef FASTPROG
	if (ch>0) { EMATDEBUG(CURDLEV+2,e,"\t\t\t"); EMATDEBUGS(CURDLEV+2,re,"\t\t\t fc\t"); }
	else if (CURDLEV+1<=DEBUGLEV) {
		DEBUG(CURDLEV+1," given re:\t");  for (i=0; i<ROWSM(re); i++)  SDEBUG(CURDLEV+1,"(%d) ",ROWSID(re,i));
		SDEBUG(CURDLEV+1," x  ");  for (i=0; i<COLSM(re); i++)  SDEBUG(CURDLEV+1,"(%d) ",COLSID(re,i));
		SDEBUG(CURDLEV+1,"\n");
	}
	if (ch>0)  EMATDEBUGS(CURDLEV+2,ec2,"\t\t\t->\t");
	else if (CURDLEV+1<=DEBUGLEV) {	
		DEBUG(CURDLEV+1,"-> full cl:\t");  for (i=0; i<ROWSM(ec2); i++)  SDEBUG(CURDLEV+1,"(%d) ",ROWSID(ec2,i));
		SDEBUG(CURDLEV+1," x  ");  for (i=0; i<COLSM(ec2); i++)  SDEBUG(CURDLEV+1,"(%d) ",COLSID(ec2,i));
		SDEBUG(CURDLEV+1,"\n");
	}
#endif
	return ec2;
}



/**
 * This function looks whether the matrix e has a "width-3 tail" of the given length.
 * Moreover, if tr>=0 is given, then it assumes that the searched tail did not exist
 * before the last line lo (row for tr==0, column for tr==1) was added to the matrix.
 * A "tail of width 3" is a piece of the matroid having path-width 3, and it is actually
 * obtained by a full closure of a triangle or a triad.
 * See struct_fclosure_ext() above.
 * The matrix must be properly represented and 3-connected when calling this function!
 * 
 * If tr==-1, then the the whole matrix is searched without respect to the (last) line lo.
 * If tail==-1, then the length of the longest tail is returned.
 * If lo==-1, then the last line is considered (determined row/col by tr) as the added one.
 * If txl!=NULL, then only the list of starting submatrices from txl is considered when
 * making full closures.
 * 
 * The value of tmod allows to request modifications to the tail size for particular tail
 * types (size up or down), as defined by the macros TMOD_xxx (see in the header).
 * This tail type must be self-dual.
 * Precisely speaking, when a sub-tail of the particular type is found inside the whole tail
 * (or if the whole tail is of that particular type for TMOD_xxxONLY flag),
 * then the value of TMOD_xxx modifier is ADDed to the length of the particular-type tail, 
 * and the resulting tail length is adjusted accordingly up.
 * For negative modifier value, if the whole tail length plus mod falls below the particular
 * type length, then the resulting tail length is adjusted accordingly down.
 * See more inside struct_hastail3_core() ...
 * 
 * The function returns 0 if there is NO such tail, and 1 if there is some;
 * or the length of the longest tail if tail==-1 is given (subject to modifications by tmod).
**/

int	struct_hastail3_ext(int ch, ematrix *e, int tr, int lo, ematrix **txl, int tail, int tmod) {
	ematrix		*ee, *re;
	int		i,k, x,y,z, tl, ret=0;
	
	if (tail>=0 && tail<3)  {PROGERROR("What is a tail of size %d ???",tail);}
	if (tr>=0 || txl!=NULL)  {PROGERROR("Not implemented yet");  return lo=0;}
	if (lo<0)  lo = tr? COLSM(e)-1:ROWSM(e)-1;
#ifndef FASTPROG
	if (ch>=0 && IFRANDDEBUGLESS(222))
		if (ematrix_inpfield_rand(e)<0)  {PROGERROR("The matrix must be properly represented here!");}
	if (!inrecur)  DEBUG(CURDLEV,"looking for tail3 of len %d (tmod %X) in %p %dx%d (tr=%d, lo=%d%s)\n",tail,tmod,e,ROWSM(e),COLSM(e),tr,lo,txl?", txl":"");
	if (ch>0)  EMATDEBUG(CURDLEV+1,e,"\t\t-..\t");
#endif
	ee = ematrix_copy(e);	/* (works on a copy of the matrix, but does not change it) */
	
		/**
		 * A tail of width 3 must end with a triangle or a triad, so we try to look
		 * for all such triples in three loops (x,y,z).
		 * We use the above function struct_hasfan_triax() to test for them.
		 * If a triangle or a triad is found, we compute its full closure in the matrix,
		 * which gives us the whole tail.
		 * Then we call struct_hastail3_core() to compare the size of the tail
		 * with the given parameter (using tmod as well), or to find the longest tail.
		**/
	ret = tl = 0;
	for (x=0; x<ROWSM(ee)+COLSM(ee) && !ret; x++)
	  for (y=x+1; y<ROWSM(ee)+COLSM(ee) && !ret; y++)
		for (z=y+1; z<ROWSM(ee)+COLSM(ee) && !ret; z++) {
			//********* how this computation should be when not 3-connected???
			if (!struct_hasfan_triax(0,ee,0,x,y,z,0) && !struct_hasfan_triax(0,ee,1,x,y,z,0))
				continue;	/* (only triples that are triangles or triads are interesting) */
			re = ematrix_refer_empty(ee);
			ematrix_refadd_rc(re,x<ROWSM(ee)?x:-1,x<ROWSM(ee)?-1:x-ROWSM(ee));
			ematrix_refadd_rc(re,y<ROWSM(ee)?y:-1,y<ROWSM(ee)?-1:y-ROWSM(ee));
			ematrix_refadd_rc(re,z<ROWSM(ee)?z:-1,z<ROWSM(ee)?-1:z-ROWSM(ee));
			
			k = struct_hastail3_core(ch,ee,re,tail,tmod);
			if (k>tl)  tl = k;	/* (the longest tail found so far) */
			dispose_ematrix(re);
			if (tail>0 && tl>=tail)  ret = 1;
		}
	//******* for txl!=NULL, use the list for triangles and triads, also consider the last line...
	
	if (tail<0)  ret = tl>=3? tl:0;
	else  ret = ret? 1:0;	/* return value indicates existence of a tail, or the longest tail */
	i = 0;
#ifndef FASTPROG
	DEBUG(CURDLEV-1+(!ret)+(!!inrecur),"- hastail(%d) test for %p (l %d, tr=%d, tmod %X) found max %d - %s\n",tail,e,lo,tr,tmod,tl,ret?"HAS":"no tail");
	if (tr<0 && txl==NULL && ch>=0) if (IFRANDDEBUG(111) && !inrecur) {
		inrecur = 1;		/* (debug-testing the full-test for diff basis) */
		DEBUG(CURDLEV-0," ... hastail3 test checking another basis recursively now ...\n");
		for (x=0; x<5; x++) {
			i = random()%ROWSM(ee);  y = random()%COLSM(ee);
			if (SIGNM(ee,i,y)!=0)  ematrix_pivot(ee,i,y);
		}
		ematrix_transpose(ee);
		if (ret!=struct_hastail3_ext(ch,ee,-1,-1,NULL,tail,tmod))  {PROGERROR("The result must be the same! ret=%d",ret);}
		ematrix_transpose(ee);  inrecur = 0;
	}		/* (ee is changed here !!!) */
#endif
	dispose_ematrix(ee);	/* (ee is a copy of e used in the computation) */
	return ret;
}


/**
 * The computation of one tail is performed here, including tmod modifications for special tail types.
 * The tail length is returned, after tmod modifications (see TMOD... in frame.h).
 * 
 * For use only in struct_hastail3_ext()...
**/

int	struct_hastail3_core(int ch, ematrix *ee, ematrix *re, int tail, int tmod) {
	ematrix		*rf;
	int		i, lc,ln,lnn,ll, tm;
	
	rf = struct_fclosure_ch(DECRCH(ch),ee,re);
	lc = ln = lnn = ROWSM(rf)+COLSM(rf);
	ll = 0;
	tm = ISTMOD_FAN(tmod);		/* modify the length if the tail is a fan (or close to) */
	if ((tm!=0 && tail<0) || (ln<tail && ln+tm>=tail) || (ln>=tail && ln+tm<tail)) {
		//printlev -= 1;
		ll = struct_hasfan_part(ch>1?ch-1:0,ee,rf,-1);	/* (ch<0 has spec. meaning in hasfan!) */
		//printlev += 1;
		if (ll==ln || !ISTMOD_FANONLY(tmod)) {
			if (tm>0 && ll+tm>ln)  lnn = ll+tm;
			if (tm<0 && ll>ln+tm)  lnn = ln+tm;
		}
	}
	//******** other tail-length modifications similarly... - but against the original lc !!!
	ln = lnn;
	i = 0;
#ifndef FASTPROG
	if (lc!=ln)  DEBUG(CURDLEV+1,"  - tail3 tmod-changed (%X) length %d -> %d (ll=%d)\n",tmod,lc,ln,ll);
	if (tail>0 && ln>=tail && !inrecur) {	/* if a sufficiently long tail is found... */
		DEBUG(CURDLEV,"  - tail3 of length %d>=%d (tmod %X) found here\n",ROWSM(rf)+COLSM(rf),tail,tmod);
		if (ch>0 && CURDLEV<=DEBUGLEV)  EMATDEBUGS(CURDLEV,rf,"\t\t-->\t");
		else if (CURDLEV-1<=DEBUGLEV && ch>=0) {
			DEBUG(CURDLEV-1,"  --> tail3 %d~%d: ",lc,ln);  for (i=0; i<ROWSM(rf); i++)  SDEBUG(CURDLEV-1,"(%d) ",ROWSID(rf,i));
			SDEBUG(CURDLEV-1," x  ");  for (i=0; i<COLSM(rf); i++)  SDEBUG(CURDLEV-1,"(%d) ",COLSID(rf,i));
			SDEBUG(CURDLEV-1,"\n");	}
	}
#endif
	dispose_ematrix(rf);
	return ln;
}




#endif
































