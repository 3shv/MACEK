
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
 * Some theory about generating matroids for start...
 * The master description of the generating process is in ../include/gener.h.
 * 
 * This file contains the main generating functions.
 * The function gener_extframe_ext() is the general interface for this process.
 * (It converts relevant information from the given frame to the elim-sequence,
 * and calls the internal extension-functions.
 * The resulting elim-sequences are then stored back into a list of frames.)
 * Read the program manual about extension-handling options in frames.
 * 
 * The remaining functions here implement the generating process, (see also genstep.c).
 * The function gener_matextens_() from gmatext.c is used to generate all possible
 * extension lines to a matrix (!up to scaling! which cannot be filtered out later).
 * Potential elim sequences for these extensions are created in gener_extensions().
 * Then these extensions are filtered through all required tests (at various levels),
 * as implemented in gener_extensions_pass().
 * One may control the test classes applied in this filtering with the ctrl parameter,
 * as described in include/gener.h (GCTRL_* macros).
 * 
**/


#include "macek.h"
#include "gen.h"















/******************	Extensions of a frame	*******************/
/******************************************************************/


#undef CURDLEV
#define	CURDLEV		6		/*6 (use <=DEBUGLEV+2 to perform extra checks) */



/**
 * This function extends the matrix(oid) of the given frame fr by a column if trq==1,
 * or by a row if trq==0, respecting other extension parameters given in the frame.
 * Keep in mind that only non-equivalent "properly" connected (see ESCONNECT)
 * extensions satisfying the additional parameters are generated here.
 * The optional parameter ctrl may be used to control the generating process
 * (default 0, set also by GEN_OPTPREFIX"ctrl", see gener_extensions() for more).
 * 
 * The return list of new frames holding the generated matrices is stored in *frout,
 * these frames have no parents or sons;
 * and the return value is the number of generated frames, or -1 for an error.
 * If ch>1(2,3) is given, then the generated extensions are printed out.
 * 
 * The parameters of the generating process are given by options in the frame fr.
 * The two main options (that must be present right in the frame fr, no inheritance)
 * are GEN_OPTPREFIX"bsize" and GEN_OPTPREFIX"signature".
 * They determine the base minor dimensions (top-left submatrix of fr) and the sequence
 * signature, as defined in ../include/gener.h.
 * These options are automatically stored to the new frames (signature updated),
 * plus additional options that are requested by "extinherit[all]".
 * The new frames also get names determined by the pattern GEN_??EXTNAME, and comments.
 * The sequence signature should never be set by hand(!), but taken from previous computations.
 * 
 * More attributes of the elimination sequence are obtained next from other options.
 * (These options are already inherited, and they are not automatically written to the
 * new frames.)
 * For example, one may give list of forbidden minors, length of forbidden fan, etc.
 * See the description of options in ../frame/fropts.c.
 * It is possible to implement more options for an elimination sequence in the function
 * gener_extframe_incl() in gener-more.inc.
 * To store requested options in the newly generated frames, one must list their option
 * names within the options "extinherit" or "extinheritall".
 * 
 * Find more information in  include/gener.h ...
**/

int	gener_extframe_ext(int ch, framestruc *fr, int trq, int ctrl, framestruc ***frout) {
	ematrix		*e,*emin, *frem, *ee,*eq, **eel, **eml;
	elimseqm	*q, **qel,**qq;
	framestruc	*frc,*frn, **flist;
	optionstruc	*op, **oih;
	int		i,j, r;
	long		b, bsz[3], cnn[3];
	char		**ml,**s,*sa, buf[400],buf2[30];
	
	if (fr?FRMATRIX(fr)==NULL:1)  {PROGERROR("No frame/matrix given to gener_extframe()!");  return -1;}
	if (frout?*frout:0)  {PROGERROR("The given list frout is supposed to be NULL-initialized!");  *frout = NULL;}
	frem = FRMATRIX(fr);
	DEBUG(CURDLEV-3,"Calling to get %sextensions of the seq %p[%.10s] (%dx%d) in %s...\n",
			trq?"column ":"row co-",fr,FRNAME(fr),ROWSM(frem),COLSM(frem),pfield_curname());
	EMATDEBUG(CURDLEV-1,frem,"\t\t");
	r = 0;
		/**
		 * Here we set up the basic attributes of a new elimination sequence
		 * according to the current frame fr.
		 * These attributes are kept together by the elim structure q,
		 * which is based on the matrix e copied from the frame.
		 * Notice that e is set in the 0-transp state as required by elim sequences.
		 * No transposition is applied to e here, but trq is passed through.
		 * We find out the base-minor size in the frame, and the sequence signature.
		 * We also set the name for the sequence by fr, but this name is not copied(!).
		 * (The matrix e and the base minor emin are freed from q at the end.)
		**/
	q = new_elimseqm();
	e = ematrix_copy_notr(frem);	/* (FRMATRIX(fr)) */
	i = frame_getoptionnum_noanc(fr,GEN_OPTPREFIX"bsize",bsz,2,0);
	if (!i) {
		DEBUG(CURDLEV-2,"## no base minor is given, assuming equal to the whole matrix...\n");
		emin = e;
	} else {
		if (bsz[0]<=0)  bsz[0] += ROWSM(e);  if (bsz[1]<=0)  bsz[1] += COLSM(e);
		DEBUG(CURDLEV-2,"- base minor of size %ld x %ld is given...\n",bsz[0],bsz[1]);
		emin = ematrix_refer(e,0,bsz[0],0,bsz[1]);
	}
	i = frame_getoptionnum_noanc(fr,GEN_OPTPREFIX"signature",bsz,1,0);
	if (!i) {
		bsz[0] = (1<<(COLSM(e)-COLSM(emin)))-1;
		DEBUG(CURDLEV-2,"## no elim-sequence signature is given, assuming %ld...\n",bsz[0]);
	}
	q = elimseqm_setfor_noch(q,e,ematrix_copy_notr(emin),bsz[0]);	/* (no check is performed on q here) */
	if (emin && emin!=e)  dispose_ematrix(emin);
	ESNAME(q) = FRNAME(fr);		/* (do not ever free the name from q !!!) */
	if (ESLENGTH(q)>GEN_MAXELIMLENGTH) {
		USERERROR("This elimination sequence is too long for handling.");
		return -20;
	}
		/**
		 * Optionally, different connectivity requirement (than 3-conn) may be given
		 * here, see the description at ESCONNECT().
		**/
	if (frame_getoptionnum(fr,GEN_OPTPREFIX"connect",cnn,1,0)) {
		ESCONNECT(q) = cnn[0];
	}
	{static char	*conoptx[] = {"cosimple","simple","connected","3connected",NULL};
	for (j=0; conoptx[j]; j++) {
		strcpy(buf,GEN_OPTPREFIX);  strcat(buf,conoptx[j]);
		if (frame_getoptionnum(fr,buf,NULL,0,0))
			ESCONNECT(q) = j;
	}
	if (ESCONNECT(q)>3 || ESCONNECT(q)<0)  {PROGERROR("No other connect values than 0,1,2,3 are implemented so far!, %d",ESCONNECT(q));}
	if (ESCONNECT(q)!=3)  DEBUG(CURDLEV-3,"- the sequence is not generated as 3-connected, only %s (%d)...\n\t\t\t\tDo you know what you are doing?\n",conoptx[ESCONNECT(q)],ESCONNECT(q));
	}
		/**
		 * Additional attributes are taken from fr here as matrix-lists.
		 * That is a list of forbidden minors, and a list of minors for a "tight major".
		 * These lists are not copied or modified in the generating process.
		 * (The lists are disposed at the end.)
		**/
	ml = frame_getoptionval_all(fr,GEN_OPTPREFIX"forbid");
	for (eml=NULL, s=ml; s?*s:0; s++) {
		eel = frame_inputmatrices(*s);
		if (eel)  eml = alist_applist(eml,eel);
	}
	if (ml)  dispose_alist(ml);
	if (eml)  DEBUG(CURDLEV-1,"- a list of forbidden minors of length %d is given... [%s %s ]\n",
			alist_getlength(eml),eml[0]?EMNAME(eml[0]):"",eml[0]&&eml[1]?EMNAME(eml[1]):"");
	if (eml)  ESFORBID(q) = eml;
	ml = frame_getoptionval_all(fr,GEN_OPTPREFIX"tight");
	for (eml=NULL, s=ml; s?*s:0; s++) {
		eel = frame_inputmatrices(*s);
		if (eel)  eml = alist_applist(eml,eel);
	}
	if (ml)  dispose_alist(ml);
	if (eml)  DEBUG(CURDLEV-1,"- a list of \"tight-major\" minors of length %d is given... [%s %s ]\n",
			alist_getlength(eml),eml[0]?EMNAME(eml[0]):"",eml[0]&&eml[1]?EMNAME(eml[1]):"");
	if (eml)  ESTIGHT(q) = eml;
		/**
		 * One may give a prefix entries for the generated extension lines
		 * as one matrix row, as in '@ext-entries "e n t r i e s"'.
		 * These are stored in ESPREFIXM(q) as a matrix, and later used in
		 * gener_matextens_ext() and in elimchm_prefixline_ext().
		**/
	ml = frame_getoptionval_last(fr,GEN_OPTPREFIX"entries");
	buf[0] = 0;  j = 0;
	for (s=ml; (s?*s:0); s++) {
		if (*ml[0] && (j>0 || *ml[0]!='<'))  buf[j++] = ' ';
		for (sa=*s; *sa && j<380; sa++)
			if (*sa==',' || *sa==';' || *sa=='\n')  buf[j++] = ' ';
			else  buf[j++] = *sa;
	}
	buf[j] = 0;  if (ml)  dispose_alist(ml);
	if (j>1) {
		eel = frame_inputmatrices(buf);
		if (eel) if ((eel[0]?ROWSM(eel[0])*COLSM(eel[0]):0)>0)
			ESPREFIXM(q) = ematrix_copy(eel[0]);
		if (eel)  dispose_alist_mats(eel);
	}
	if (ESPREFIXM(q))  DEBUG(CURDLEV-1,"- a \"prefix\" for generating ext lines is given: %s\n",gener_extprintline(0,ESPREFIXM(q),0));
	
		/**
		 * More attributes (number-valued) are obtained here (like forbidden fans).
		 * 
		 * Other parameters will be implemented later......
		**/
	i = frame_getoptionnum(fr,GEN_OPTPREFIX"nofan",bsz,1,0);
	if (i>0 && bsz[0]<11111 && bsz[0]>=3)  ESNOFAN(q) = bsz[0];
	if (ESNOFAN(q)>0)  DEBUG(CURDLEV-1,"- no fan of length >=%d along the elim sequence, given...\n",ESNOFAN(q));
	//********** no-tails, and more parameters.........
		/**
		 * The given value of ctrl may be further modified (in OR mode) by options
		 * of the frame fr.
		 * (See gener_extensions() for the meaning of the ctrl parameter.)
		 * This is an undocumented feature of the program....
		**/
	i = frame_getoptionnum(fr,GEN_OPTPREFIX"ctrl",bsz,1,0);
	if (i>0 && bsz[0]>0)  ctrl |= bsz[0];
	//********* more specific ctrl bits come in here (if needed)......
	if (!ESPREFIXM(q))  ctrl |= GCTRL_NOSPECTEST;
	if (ctrl!=0 && ctrl!=GCTRL_NOSPECTEST) {
		DEBUG(CURDLEV-3,"- using special control parameter ctrl = %x=%d...\n",ctrl,ctrl);
		SDEBUG(CURDLEV-2,"\t\t\t\t(%s%s%c%s%s)\n",(ctrl&GCTRL_NOPRETEST)?"nopretest, ":"",
				GCTRL_CANONLEV(ctrl)?"canonlev ":"",GCTRL_CANONLEV(ctrl)?'0'+GCTRL_CANONLEV(ctrl):' ',
				(ctrl&GCTRL_NOSEQTEST)?", noseqtest":"",(ctrl&GCTRL_NOSTRUCTEST)?", nostructest":"");
	}
	
		/**
		 * Now we check that the elim sequence is given correctly, which includes
		 * proper signature, excluded minors, connectivity, and canonical minimality.
		 * Passing these properties is required for correct computation of extensions!
		**/
	if (r>=0) if (elimchm_signbits_check(q)<0) {
		USERERROR("The elimination sequence signature of the frame %p [%s] is wrong!",fr,FRNAME(fr));
		r = -11; }
	if (r>=0 && ESPREFIXM(q)) if (ROWSM(ESPREFIXM(q))>1) {
		USERERROR("Extension entries (\"prefix\") may be given as one matrix row only!");
		ROWSM(ESPREFIXM(q)) = 1; }
#ifndef FASTPROG
	printlev -= 1;
	if (r>=0 && ESFORBID(q)) if (struct_hasminorlist(e,ESFORBID(q))) {
		USERERROR("The matrix to extend in [%s] already has a forbidden minor in the given list!",FRNAME(fr));
		r = -20; }
	if (r>=0 && ESTIGHT(q)) if (struct_hasdelcontr(e,ESTIGHT(q))) {
		USERERROR("The matrix to extend in [%s] is not a tight major of the given list!",FRNAME(fr));
		r = -20; }
	if (r>=0) if (elimchm_xconnsequence(q,ESCONNECT(q))<0) {
		USERERROR("The elimination sequence of the frame [%s] is not properly(%d) connected!",FRNAME(fr),ESCONNECT(q));
		r = -11; }
	if (r>=0 && ESNOFAN(q)>0) if (struct_hasfan_max(e)>=ESNOFAN(q)) {
		USERERROR("The matrix to extend in [%s] has a forbidden fan %d>=%d!",FRNAME(fr),struct_hasfan_max(e),ESNOFAN(q));
		r = -11; }
	if (r>=0 && ESLENGTH(q)>0) if (gcanon_allseqtest_rand(q,e)<0) {
		USERERROR("The elimination sequence of the frame [%s] is not (canonically) minimal!",FRNAME(fr));
		r = -11; }
	if (ESLENGTH(q)<=0 && r==-11)  {USERERROR("Remark: The extension generating process requires one to start"
			" from a sufficiently\n  \"%d-connected\" matroid, and it cannot work for an arbitrary matroid by its nature.",ESCONNECT(q));}
	if (ESLENGTH(q)>0 && r==-11)  {USERERROR("Remark: The extension generating attributes (especially the sequence signature)"
			" are very fragile.\n  So do not set them by hand, but use a previously generated frame instead.");}
	printlev += 1;
#endif
		/**
		 * This part calls a third-party modification function gener_extframe_incl(),
		 * for more attributes to the sequence q.
		 * Then the generating function gener_extensions() is called, storing the
		 * generated list of (next) sequences in qel.
		 * For easier further use of the generated matrices, they are all transposed
		 * so that the extended line is the last column.
		 * Possible output (ch>1) is printed out here, and also later if ch>3.
		 * The next call to a third-party modification function gener_extframe_incl()
		 * may then prepare a "master copy" of the output frame frc.
		**/
	q = gener_extframe_incl(fr,trq,ctrl,q,NULL,NULL,NULL,1);
	qel = NULL;
	if (r>=0)  r = gener_extensions((ch>4?ch-4:0),ctrl,q,trq,(void***)&qel);
	if (ch>0 && r>=0) {
		OUTPUT("Generated %d non-equiv %s %sextens of the sequence [%.15s] (%dx%d|%dx%d).\n",
				r,ESCONNECT_S(q),trq?"column ":"row co-",FRNAME(fr),ROWSM(frem),COLSM(frem),ROWSM(ESMINOR(q)),COLSM(ESMINOR(q)));
		if (ctrl!=0 && ctrl!=GCTRL_NOSPECTEST)  OUTPUT("  (Used special control parameter ctrl = %x=%d.)\n",ctrl,ctrl);
		if (ESPREFIXM(q))  OUTPUT("  (Started with extension entries as in \'[%s..]\'.)\n",
					gener_extprintline(0,ESPREFIXM(q),0));
		if (ch>2)  ELIMQOUTPUTS(q,printoutpref);
		else if (ch>1)  EMATOUTPUT(frem,printoutpref);
	}
	if (r>=0 && ch>1 && ch<=3)  for (qq=qel,i=1; qq?*qq:0; qq++,i++)
		OUTPUT("  #%d extension  %s\n",i,gener_extprintline(1,ESMATRIX(*qq),trq));
	frc = new_frame(NULL);
	q = gener_extframe_incl(fr,trq,ctrl,q,qel,&frc,NULL,0);
	FRPFINDEX(frc) = FRPFINDEX(fr);
	
		/**
		 * This part inherits the options requested by "@extinherit[all]"
		 * for the master frame frc.
		**/
	ml = frame_getoptionval_all(fr,"extinherit");
	oih = frame_getoptionlist_last(fr,ml);
	if (ml)  dispose_alist(ml);
	if (oih)  FROPTIONS(frc) = alist_applist(FROPTIONS(frc),oih);
	ml = frame_getoptionval_all(fr,"extinheritall");
	oih = frame_getoptionlist_all(fr,ml);
	if (ml)  dispose_alist(ml);
	if (oih)  FROPTIONS(frc) = alist_applist(FROPTIONS(frc),oih);
		/**
		 * This part converts the generated list of elim sequences to output frames.
		 * The master output frame prepared above is cloned for each extension,
		 * and the matrices of the generated extensions are each copied to a new
		 * matrix obtained from FRMATRIX(fr) (since we made other changes to
		 * the matrices in elim sequences).
		 * We set the "signature" options accordingly for the new frames, also
		 * give the new names by GEN_??EXTNAME and set simple comments.
		 * (The options previously set to the master frame frc are also copied.)
		**/
	flist = NULL;
	if (r>=0)  for (qq=qel,i=1; qq?*qq:0; qq++,i++) {
		frn = frame_copy_full(frc);
		ee = ematrix_copy(FRMATRIX(fr));  ee = ematrix_append_rc(ee,trq);
		eq = ESMATRIX(*qq);
		if (trq)  for (j=0; j<ROWSM(ee); j++)
			COPYEXSIGM(ee,j,COLSM(ee)-1,eq,j,COLSM(eq)-1);
		else  for (j=0; j<COLSM(ee); j++)
			COPYEXSIGM(ee,ROWSM(ee)-1,j,eq,ROWSM(eq)-1,j);
		frame_setmatrix(frn,ee);
		op = new_option_onenum(GEN_OPTPREFIX"signature",ESELIMSG(*qq));
		frame_addoption(frn,op,2);
		op = new_option_twonum(GEN_OPTPREFIX"bsize",ROWSM(ESMINOR(*qq)),COLSM(ESMINOR(*qq)));
		frame_addoption(frn,op,3);
		snprintf(buf,380,trq?GEN_EXTNAME:GEN_COEXTNAME,FRNAME(fr),i);
		buf[380] = 0;  frame_setname(frn,buf);
		if (ch>3) {
			OUTPUT(" #%d extension %s:\n",i,FRNAME(frn));  EMATOUTPUTS(ee,"\t\t");
		}
		for (j=ESLENGTH(*qq)-1, b=ESELIMSG(*qq); j>=0; j--, b/=2)
			buf2[j] = (b%2)?'c':'r';
		buf2[ESLENGTH(*qq)] = 0;
		snprintf(buf,380,"m#%d %sext \'%.13s\', seq %dx%d, %s(%ld)",
			i,trq?"c-":"r-",FRNAME(fr),ROWSM(ESMINOR(*qq)),COLSM(ESMINOR(*qq)),buf2,(long)ESELIMSG(*qq));
		buf[380] = 0;  FRCOMMENT(frn) = MSTRDUP(buf);
		if (frout)  flist = alist_append(flist,frn);
		else  dispose_frame(frn);
	}
	if (r>=0 && r!=i-1)  {PROGERROR("wrong number of output frames converted ?!?");}
	if (FRPARENT(frc) || FRSONS(frc))  {PROGERROR("the master frame frc should not have parent or sons %p %p %p",frc,FRPARENT_XD(frc),FRSONS_XD(frc));}
	dispose_frame(frc);
	if (ch>1 && r>=0)  OUTPUT("Total %d non-equiv %s %sextensions of the sequence [%.10s] generated.\n",r,ESCONNECT_S(q),trq?"column ":"row co-",FRNAME(fr));
	if (ch>0 && r<0)  OUTPUT("Generating %s %sextensions of the frame [%.10s] has failed %d.\n",ESCONNECT_S(q),trq?"column ":"row co-",FRNAME(fr),r);
	
		/**
		 * Here the third-party modification is called again to (possibly) add
		 * more options to the output frames in flist, and to destroy its local data.
		 * Finally, the data associated with the elimination sequence are freed.
		 * (Recall that the extra data is not copied, just refered in the ext sequences.)
		 * The matrices from qel were already copied to the output frames.
		**/
	q = gener_extframe_incl(fr,trq,ctrl,q,qel,NULL,flist,-1);
	for (qq=qel; qq?*qq:0; qq++)
		if (ESMATRIX(*qq))  dispose_ematrix(ESMATRIX(*qq));
	if (qel)  dispose_alist(qel);	/* (the elim-sequences are freed with their list here) */
	if (ESFORBID(q))  dispose_alist_mats(ESFORBID(q));
	if (ESTIGHT(q))  dispose_alist_mats(ESTIGHT(q));
	if (ESPREFIXM(q))  dispose_ematrix(ESPREFIXM(q));
	dispose_ematrix(ESMINOR(q));  dispose_ematrix(ESMATRIX(q));
		/* (ESMINOR(q) is a private copy, and ESMATRIX(q)==e here, both to be freed) */
	FREE(q);
	if (frout)  *frout = flist;
	return r;
}













/******************	Generating matrix extensions	*******************/
/**************************************************************************/


#undef CURDLEV
#define	CURDLEV		6	/* (>=5 is OK for randomized checks, use <=DEBUGLEV+2 for extra checks) */


		/* these variables are used for profiling output of the generating process: */
int	stseq3, ststruc3, stcanon3, ststruc2, stcanon2, ststruc1, stcanon1, stseq0, ststruc0, stcanon0;


/**
 * This function generates (co)extensions to the matrix in the given elim sequence q.
 * All distinct extension lines are generated by calling gener_matextens_..().
 * The generated extensions are then tested by the genstep_????check() functions from
 * genstep.c (see gener_extensions_pass()), and only proper extensions pass through.
 * Internal parameters of the sequence structure q are in control of these tests.
 * 
 * If trq==1, then we compute column extensions, while for trq==0 row co-extensions
 * (both against the 0-trasnposition state of the matrix!).
 * We use "column implementation", so the matrices are first transposed for row coextension.
 * The parameter ctrl controls what tests are performed:
 * ctrl==0 means all tests, ctrl>0 only some (or no) canonical tests or no pre-tests.
 * (Use macros like GCTRL_NOPRETEST, GCTRL_CANONLEV(l) in "or" expr.)
 * Read more about performing canonical tests in gener.h ...
 * 
 * If ch>1, then the sequences are printed along.
 * The function returns a list *qout of pointers to new sequence structures holding
 * all generated and passed matrices, with internal parameters copied from the sequence q.
 * (Remember to free the list and the structures later; but not the refered data.)
 * The number of extensions (or <0 for an error) is returned.
**/

int	gener_extensions(int ch, int ctrl, void *vq, int trq, void ***qout) {
	int		i,k,r, npas, tri, pre;
	ematrix		*ei,*ee, **exl=NULL,**x;
	elimseqm	*q=vq, *qn, *qo;
	
	ei = ESMATRIX(q);  tri = ISTRANSPM(ei);
#ifndef FASTPROG
	if (!ei || (ESMINOR(q)? (ROWSM(ESMINOR(q))<=0||COLSM(ESMINOR(q))<=0):1))
		{PROGERROR("Call only for nonempty matrix %p and nonempty base minor %p !",ei,ESMINOR(q)); return -1;}
	if (qout?*qout:0)  {PROGERROR("Expect an empty list *qout to be given.");}
	DEBUG(CURDLEV-2,"Generating %sextensions of a (%dx%d|%dx%d) sequence [%s] over %s...\n",
			trq?"column ":"row co-",ROWSM(ei),COLSM(ei),ROWSM(ESMINOR(q)),COLSM(ESMINOR(q)),ESNAME(q),pfield_curname());
	if (ctrl!=0)  DEBUG(CURDLEV-2,"   - canonlev=%d, nopretest %c\n",GCTRL_CANONLEV(ctrl),(GCTRL_NOPRETEST&ctrl)?'Y':'N');
	ELIMQDEBUGS(CURDLEV-0,q,"  =\t");
	stseq3 = ststruc3 = stcanon3 = ststruc1 = ststruc2 = stcanon2 = stcanon1= stseq0 = ststruc0 = stcanon0 = 0;
#endif		/* (the static variables stseq?... are used for debug statistics of the computation) */
	if (ch>1)  OUTPUT("Generating all %s %sextens (ctrl=%d) to sequence [%s]:\n",ESCONNECT_S(q),trq?"column ":"row co-",ctrl,ESNAME(q));
	
		/**
		 * We prepare the variables for generating extensions:
		 * ei is the given matrix in q, qn is the new elim sequence which is
		 * later copied for the output, and exl gets the list of all possible
		 * matrix extensions generated in gener_matextens().
		 * The matrix ei is transposed so that we are going to extend a column,
		 * and the line id's of ei are reset to default labels.
		**/
	ematrix_transpose_set(ei,!trq);  ematrix_resetid(ei);
	ee = ematrix_append_col(ematrix_copy(ei));
	qn = elimseqm_next(q,ee,trq);
	r = 0;  exl = NULL;
	if (genstep_seqsigncheck(qn,ESELIMSG(qn))<0) {
		DEBUG(CURDLEV-3,"!!! NOT accepting the next sequence signature %lX !!!\n",ESELIMSG(qn));
		r = -11;
	}
		/**
		 * We cycle through all generated matrix extensions, and look whether
		 * they pass the admissibility tests in gener_extensions_pass().
		 * If so, then a new copy of the elim sequence qn is created for the
		 * matrix, and stored in the output list (in the original transp-state).
		 * Finally, some quick tests of the sequence are run.
		 * 
		 * An optional prefix for the generated extensions may be specified
		 * in ESPREFIXM(qn) (as the entries in the first-one row of ESPREFIXM(qn)).
		 * Possibly, genstep_pretest_start() is called first to prepare some data
		 * used in faster admissibility tests inside gener_extensions_pass().
		**/
	if (r>=0) {
		if (!ESPREFIXM(qn))  r = gener_matextens_col(ei,&exl);
		else  r = gener_matextens_prefix(ei,ESPREFIXM(qn),&exl);
	}
	if (r>=0 && GEN_USEPRECHECK && (ctrl&GCTRL_NOPRETEST)==0) {
		genstep_precheck_start(qn,1);  pre = 1;
	} else  pre = 0;
	DEBUG(CURDLEV-2,"== cycling through %d potential extension vectors...\n",alist_getlength(exl));
	npas = 0;
	for (x=exl,i=1; (x?*x:0) && r>=0; x++,i++) {
		DEBUG(CURDLEV-2,"  - extension %s:\n",gener_extprintline(0,*x,1));
		k = gener_extensions_pass(ch,ctrl,qn,1,*x,i);
		if (k<0)  continue;
		++npas;
		qo = elimseqm_copyfor(qn,ematrix_copy(*x));
		ematrix_transpose_set(ESMATRIX(qo),tri);
#ifndef FASTPROG
		DEBUG(CURDLEV-2,"==>>  extension #%d (%d) of [%s] accepted:  %s.\n",npas,i,ESNAME(q),gener_extprintline(0,*x,1));
		EMATDEBUGS(CURDLEV-0,ESMATRIX(qo),"\t++\t");
				/* final debug testing (some even in the fastprog mode!) */
		if (IFRANDDEBUGLESS(111))  test_elimseqm(qo);
		else if (IFRANDDEBUG(22))  test_elimseqm_rand(qo);  else
#endif
		test_elimseqm_notough(qo);	/* (always at least the fast test...) */
		if (ch>1)  OUTPUT(" =>> accepted extension #%d to [%s]:  %s.\n",npas,ESNAME(q),gener_extprintline(1,*x,1));
		if (qout)  *qout = alist_append(*qout,qo);
		else { dispose_ematrix(ESMATRIX(qo));  FREE(qo); }
	}
	if (exl)  dispose_alist_mats(exl);
	if (ch>0)  OUTPUT("Total %d %s %sextens to elimseq [%s] %dx%d generated.\n",npas,ESCONNECT_S(q),trq?"column ":"row co-",ESNAME(q),ROWSM(ei),COLSM(ei));
	if (ch>0 && ctrl>0)  OUTPUT("(Used special ctrl=%d:  %s%s%c%s%s)\n",ctrl,(ctrl&GCTRL_NOPRETEST)?"nopretest, ":"",
				GCTRL_CANONLEV(ctrl)?"canonlev ":"",GCTRL_CANONLEV(ctrl)?'0'+GCTRL_CANONLEV(ctrl):' ',
				(ctrl&GCTRL_NOSEQTEST)?", noseqtest":"",(ctrl&GCTRL_NOSTRUCTEST)?", nostructest":"");
	
		/**
		 * At last, we free the above precomputed information for admissibility tests
		 * (indicated by pre==1), and the matrices used in generating extensions.
		**/
	if (pre)  genstep_precheck_stop(qn,1);
	FREE(qn);  dispose_ematrix(ee);
	ematrix_transpose_set(ei,tri);
	DEBUG(CURDLEV-3,"Gener - passed %d out of %d (co)extensions of %dx%d matrix %p[%s].\n",npas,(int)(x-exl),ROWSM(ei),COLSM(ei),ei,ESNAME(q));
	SDEBUG(CURDLEV-3,"\t(seq3=%d, struc3=%d, canon3=%d, struc2=%d, canon2=%d, struc1=%d, canon1=%d, seq0=%d, struc0=%d, canon0=%d)\n",stseq3,ststruc3,stcanon3,ststruc2,stcanon2,ststruc1,stcanon1,stseq0,ststruc0,stcanon0);
	SDEBUG(CURDLEV-2,"\t====================================================================================================\n\n");
	return npas;
}



/**
 * Here we get a generated extension ee of the sequence qn, and we test whether this
 * extension passes all requirements of the sequence qn.
 * Recall that the last added line is a column of ee, and so tr==1.
 * The parameter ctrl controls what tests are performed:
 * ctrl==0 mean all tests, ctrl>0 only some (or no) canonical tests, or no pre-tests.
 * Read more about performing canonical tests in include/gener.h ...
 * 
 * The function returns -1 for a rejected matrix, and 1 for an accepted one.
 * For internal use only!
**/

int	gener_extensions_ix = 0;

int	gener_extensions_pass(int ch, int ctrl, elimseqm *qn, int tr, ematrix *ee, int ix) {
	int		i, r=0;
	ematrix		*eqo;
	
	junk = ch;		/* (not really useful here) */
	if (!ee || tr!=1)  {PROGERROR("Wrong parameters %p %d",ee,tr); return -1;}
#ifndef FASTPROG
	if (elimchm_unitlastline(0,ee,1,0)<0) {
		if (!ESPREFIXM(qn))  DEBUG(CURDLEV-2,"Extension (%d) %s of [%s] is not given in unit-leading form!\n",ix,gener_extprintline(-1,ee,1),ESNAME(qn));
		else  {PROGERROR("Extension (%d) %s of [%s] must be given in unit-leading form when req. prefix!",ix,gener_extprintline(-1,ee,1),ESNAME(qn));}
	}	/* when filtering by prefix entries, we must not unit-scale the line, as the filter was already pre-applied above! */
#endif
	gener_extensions_ix = ix;
	eqo = ESMATRIX(qn);
	qn = elimseqm_setfor_mat(qn,ee);
		/* genstep_pretest_start()/stop() is possibly called in gener_extensions() above... */
	
		/**
		 * First we set the extended line into canonical scale -- unit-leading form.
		 * (The line-check has been already performed together with unit-scaling.)
		 * Actually, for the right function of '@ext-entries ...', the extension
		 * lines already should be in unit-leading form, which is now checked above
		 * at debug level / error level when using '@ext-entries ...'.
		 * 
		 * Then we continue with partial-tests for faster rejection (at levels lev>0):
		 * If these tests reject a matrix, then it is for sure, but some bad matrices
		 * may still pass through.
		 * The full tests are carried out in the next section.
		**/
	if (r>=0)  r = genstep_linecheck_set(qn,ee,1);
	if (r>=0 && (ctrl&GCTRL_NOSPECTEST)==0)  r = genstep_specialcheck_lev(3,qn,ee,1);
	
	if (r>=0 && (ctrl&GCTRL_NOSEQTEST)==0)  r = genstep_sequencecheck_lev(3,qn,ee,1);
	else if (r>=0)  r = genstep_sequencecheck_lev(4,qn,ee,1);
		/* (we must always test connectivity of the last line at level 4) */
	if (r>=0)  stseq3++;
	if (r>=0 && (ctrl&GCTRL_NOSTRUCTEST)==0)  r = genstep_structurecheck_lev(3,qn,ee,1);
	else if (r>=0)  r = genstep_structurecheck_lev(4,qn,ee,1);
		/* (we must always test representability of the matrix at level 4) */
	if (r>=0)  ststruc3++;
	if (r>=0 && GCTRL_CANONLEV(ctrl)<=3)  r = genstep_canonicalcheck_lev(3,qn,ee,1);
	if (r>=0)  stcanon3++;
	
	if (r>=0 && (ctrl&GCTRL_NOSEQTEST)==0)  r = genstep_sequencecheck_lev(2,qn,ee,1);
	if (r>=0 && (ctrl&GCTRL_NOSTRUCTEST)==0)  r = genstep_structurecheck_lev(2,qn,ee,1);
	if (r>=0)  ststruc2++;
	if (r>=0 && GCTRL_CANONLEV(ctrl)<=2)  r = genstep_canonicalcheck_lev(2,qn,ee,1);
	if (r>=0)  stcanon2++;
	
	if (r>=0)  DEBUG(CURDLEV-1,">>  extension (%d) of [%s] trying more: %s\n",ix,ESNAME(qn),gener_extprintline(-1,ee,1));
	if (r>=0 && IFRANDDEBUGLESS(555))  if (elimchm_represented_last(qn,ee,1)<0)
		{PROGERROR("The matrix must already be proper in the pfield at levels <=1 !");}
	if (r>=0 && (ctrl&GCTRL_NOSEQTEST)==0)  r = genstep_sequencecheck_lev(1,qn,ee,1);
	if (r>=0 && (ctrl&GCTRL_NOSTRUCTEST)==0)  r = genstep_structurecheck_lev(1,qn,ee,1);
	if (r>=0)  ststruc1++;
	if (r>=0 && GCTRL_CANONLEV(ctrl)<=1)	/* (only if randomized may help - some other full tests are present) */
		if (ESFORBID(qn) || ESTIGHT(qn) || GCTRL_CANONLEV(ctrl)==1)
			r = genstep_canonicalcheck_lev(1,qn,ee,1);
	if (r>=0)  stcanon1++;
	
		/**
		 * These are the full tests at levels 0 - must be carried out(!):
		 * (See that genstep_seqsigncheck() was called in gener_extensions() above,
		 * and the line-check has been already performed at the start of this function...)
		 * Although the full tests should not be based on the last added line in general,
		 * we need to examine the last line for connectivity in the sequence-check.
		**/
	if (r>=0 && (ctrl&GCTRL_NOSPECTEST)==0)  r = genstep_specialcheck(qn,ee,1);
	if (r>=0 && (ctrl&GCTRL_NOSEQTEST)==0)  r = genstep_sequencecheck(qn,ee,1);
	if (r>=0)  stseq0++;
	if (r>=0)  r = genstep_basischeck(qn,ee);
	if (r>=0 && (ctrl&GCTRL_NOSTRUCTEST)==0)  r = genstep_structurecheck(qn);
	if (r>=0)  ststruc0++;
	if (r>=0 && GCTRL_CANONLEV(ctrl)==0)  r = genstep_canonicalcheck(qn);
	if (r>=0)  stcanon0++;
	
	ESMATRIX(qn) = eqo;
	return r>=0? 1:-1;
}

















/******************	Additions to generating	*******************/
/******************************************************************/


/**
 * The file gener-more.inc is provided for third-party additions to the generating process.
 * (See master description of generating in ../include/gener.h.)
 * The file (possibly) adds more fields to the elimination-sequence structure elimseqm
 * when included from eseq.h with GENER_ELIMQDEF>0.
 * When included from here with GENER_ELIMQDEF==0, the file
 * adds processing for more frame-options that control the generating process,
 * and adds more tests for the sequence and structure checks (see genstep.c).
 * 
**/

#define	GENER_ELIMQDEF 0
#include "gener-more.inc"
#undef	GENER_ELIMQDEF





































