
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
 * A "frame" is the basic general data entity in the program.
 * Read more about a frame and its data in ../include/frame.h.
 * 
 * This file defines general functions for frame command-handling.
 * (Commands are given to frames on the input, and they manipulate the frames and their matrices.)
 * A command is identified by its (string) name, and the attributes (parameters, return
 * value, etc) are declared in the all-command list in frcoms.c.
 * The name is case-insensitive.
 * Here you find just the general-level functions that identify the command, prepare
 * its parameters, and then colect the command results and store them back.
 * 
 * Read the section on handling frame commands in ../include/frame.h to learn more.
 * You find out how command parameters and results are determined, and other...
 * 
**/





#include "macek.h"  
#include "fr.h"














/******************	Applying frame commands	*******************/
/******************************************************************/


#undef CURDLEV
#define CURDLEV         5


/**
 * The strings here are used as default parameter descriptions when nothing else is given.
 * (There are different defaults for one- or for list- frames...)
 * For debug purposes, the replaced default parameters are refered in prstrings[].
**/

static char	*fcpar_oneframed[FCMAXPARAMS] = { "((t))", "(()(T))", "((2)(T))", "((3)(T))",
						"((4)(T))", "((5)(T))", NULL };
static char	*fcpar_listframed[FCMAXPARAMS] = { "((s))", "(()(S))", "((2)(S))", "((3)(S))",
						"((4)(S))", "((5)(S))", NULL };
static char	*fcpar_oneframens[FCMAXPARAMS] = { "(t)", NULL };
static char	*fcpar_listframens[FCMAXPARAMS] = { "((t))", "(()(T))", "((2)(T))", "((3)(T))", NULL };

static char	*fcpar_oneret = "((t))",	*fcpar_listret = "((s))";

static char	*prstrings[FCMAXPARAMS];	/* (no string copies here, just informal references) */


/**
 * These variables define allowable return-type alterations.
 * (It is possible to alter the return type of some commands by adding a certain prefix to it.)
 * So far, the only replacements here are for the type FCRET_PRINTED_F that can be changed
 * to a filter, or inverse filter, or just marking to remember the resulting list (without
 * storing anything to the frame tree).
 * See in frame_applycommand() below...
**/

static char	*comprefixes[] = {FRES_FILTPREFIX, FRES_INVFILTPREFIX, FRES_REMPREFIX, FRES_INVREMPREFIX, FRES_INVREMPREFIX2};
static int	comretorig[] = {FCRET_PRINTED_F, FCRET_PRINTED_F, FCRET_PRINTED_F, FCRET_PRINTED_F, FCRET_PRINTED_F},
		comretnew[] =  {FCRET_YESNOFR, FCRET_NOYESFR, FCRET_LREMFR, FCRET_LREXFR, FCRET_LREXFR};


/**
 * The first frame parameter found for the command is stored here for possible later use.
 * (Like for getting new names, or for output matrix replacement.)
 * Also, the (several) last resulting lists of frames are stored for later quick access
 * by variables ~1,~2 or ^1,^2,... (going back, ~1 is the most recent one, see FRGET_PREVRES[KEEP]).
 * 
 * However, to allow re-entrance to the function frame_processcommands(), we have
 * to implement the static variables with a local "stack" shifted by calls to frame_processcommands().
 * This is a transparent solution obtained through macros and frame_procstack_??() functions...
**/

#define	PROCSTACKMAX	101
static int	procstackix = 0;

static framestruc	*firstfrparam_a[PROCSTACKMAX], **firstfrparamlist_a[PROCSTACKMAX];
#define	firstfrparam		firstfrparam_a[procstackix<0?(PROGERROR("uninitialized stack!"),0):procstackix]
#define	firstfrparamlist	firstfrparamlist_a[procstackix]

static framestruc	**lastresults_a[PROCSTACKMAX][10];
static int		lastresultix_a[PROCSTACKMAX];
#define	lastresults		lastresults_a[procstackix]
#define	lastresultix		lastresultix_a[procstackix<0?(PROGERROR("uninitialized stack!"),0):procstackix]
#define	lastresultix_l		lastresultix_a[procstackix]

static int	stverbose[PROCSTACKMAX];

void	frame_procstack_enter(void) {
	int		i;
	extern int	frcom_verbose;
	
	stverbose[procstackix] = frcom_verbose;
	procstackix = (procstackix+1)%PROCSTACKMAX;
	for (i=0; i<10; i++)  lastresults[i] = NULL;
	lastresultix = 0;
	firstfrparam = NULL;  firstfrparamlist = NULL;
	frcom_verbose = 1;
}

void	frame_procstack_exit(void) {
	extern int	frcom_verbose;
	procstackix = (procstackix+PROCSTACKMAX-1)%PROCSTACKMAX;
	frcom_verbose = stverbose[procstackix];
}


/**
 * This function processes all commands stored in the frame tree of fr.
 * The commands are executed in the reverse depth-first order, by lists.
 * (The root frame of a command is never deleted, so fr stays valid after this function...)
 * If a frame-handle returns negative value (error), then the processing of commands stops.
 * (The return value is the last return value of a command handle.)
 * Additionally, the whole frame tree is written to the disk ("/tmp..") after an error.
 * Commands are deleted after execution from the frame.
 * Notice that the above static variables are protected here with a local "stack".
 * 
**/

int	frame_processcommands(framestruc *fr) {
	framestruc	**y, **frl;
	optionstruc	*cm,**c;
	int		tim, pi, r=0, skip;
	char		buf[300];
	
	if (!fr)  {PROGERROR("The frame must be given here!"); return -1;}
	DEBUG(CURDLEV-1,"Calling command-processing (rev depth-first) in the frame %p [%s]...\n",fr,FRNAME(fr)?FRNAME(fr):"");
	tim = TIME();
	frame_procstack_enter();
	
	frl = frame_gettree_rdep(fr);		/* process all frames in the subtree */
	for (y=frl; (y?*y:0) && r>=0; (y?y++:0)) {
		pi = pfield_curindex();  skip = 0;
		while (r>=0 && (y?*y:0)) {
			c = FRCOMMANDS(*y);	/* process all commands in each frame */
			if (c?!*c:1)  break;
			cm = *c;		/* (some commands may be skipped if previously requested) */
			if (--skip<0)  r = frame_applycommand(*y,cm);
			else  r = 0;
			for (c=FRCOMMANDS(*y); *c && *c!=cm; c++) ;
			if (*c==cm)  alist_delete(c);	/* (the commands are deleted after execution) */
			DEBUG(CURDLEV-0," proc - frame %p [%s]%s, return %x, %d com left\n",*y,FRNAME(*y)?FRNAME(*y):"",skip>=0?" skip":"",r,alist_getlength(FRCOMMANDS(*y)));

			if (r>0 && (r&FRES_SKIPCOM)) {
				if (frame_skipcommands==-1)  y = NULL;
				if (frame_skipcommands>0)  skip = frame_skipcommands;
			}
			if (r>0 && (r&FRES_RESTART)) {	/* restart the whole command tree processing */
				alist_free(frl);
				y = frl = frame_gettree_rdep(fr);
			}
		}			/* the starting pfield is returned for other frames in the tree */
		if (pi!=pfield_curindex())  pfield_switchto_fast(pi);
	}
	if (frl)  alist_free(frl);
	frame_procstack_exit();
	
	if (r<0 && TIME()-tim>30 && frame_autosave_pref) {
		snprintf(buf,290,"%smacek-err-%d",frame_autosave_pref,(int)(RANDOM()%1000));
		buf[290] = 0;  frame_setname(fr,buf);
		frame_write_tree(fr,"Automatically saved the frame tree after an error...\n");
		OUTPUT("Automatically saved the frame tree (after an error) to \"%s\"\n",FRNAME(fr));
	}
	return FRES_PROCSTRIP2(r);
}

int	frame_skipcommands = 0;	/* used to transfer the number of commands to skip from processing */









/******************	One command handling	*******************/
/******************************************************************/


#undef CURDLEV
#define CURDLEV         5


/**
 * This function is called to apply the given command com to the frame fr.
 * The result is the value returned from the command handling function, and stripped from high bits.
 * (Negative return value means an error happened.)
 * 
 * The command handling functions are defined in frcoms.c, associated with the command name.
 * This function prepares the command parameters according to the parameter values in com,
 * and to the requested parameter types for the command in frcoms.c.
 * Then it also stores the command results back to the frame subtree (or prints).
 * The major parameter-description defines are in ../include/frame.h, and their use in
 * commented in the functions frame_getparamfr_list() and frame_storeres_list() below.
 * Depending on the parameter description, the input frames may be deleted after use.
 * (However, the given frame fr is never deleted!)
 * 
 * It is possible to modify command names with special prefixes (like FRES_FILTPREFIX) that
 * change the return type (like from printing to a yes/no filter).
 * See in ../include/frame.h.
 * It is also possible to pick one of the previous resulting lists with special variables
 * ~1-~9, ^1-^9, which is done in frame_getparamfr_list() together with normal parameters.
 * 
 * ......
**/

int	frame_applycommand(framestruc *fr, optionstruc *com) {
	commhandles		cmh;
	extern int		ret_yes,ret_no,ret_def;
	int		i,j,z, jn,jf,k,m,n,retv, clpr=0, eoc=0, fmm[FCMAXPARAMS];
	long		lk,lj, npar[FCMAXPARAMS], outn[FCMAXPARAMS];
	char		*s, **stpar, *sret, *nmret=NULL, **nmretlist=NULL, buf[550];
	void		**x,**fmp, *xx, **fmpar[FCMAXPARAMS], **outf[FCMAXPARAMS];
	framestruc	**y, **fro, *frn;
	
	if (!fr)  {PROGERROR("The frame must be given here!"); return -1;}
	if (com?!OPTISCOMMAND(com):1)  {PROGERROR("Call only for a command %p \"%s\" !",com,com?OPTNAME(com):"null"); return -1;}
	if (OPTCONTINUED(com))  {PROGERROR("Cannot handle \"continued\" commands yet %p \"%s\"...",com,OPTNAME(com)); return -1;}
	DEBUG(CURDLEV-1,"# Calling handle function for \"%s\" in frame %p [%s]...\n",OPTNAME(com),fr,FRNAME(fr)?FRNAME(fr):"");
	eoc = (uerror_occured<0 || perror_occured<0);
	for (i=0; i<FCMAXPARAMS; i++) {
		fmpar[i] = outf[i] = NULL;
		npar[i] = fmm[i] = outn[i] = 0;
	}
	
		/**
		 * Looking for a handle for the given command com by its name.
		 * The command is matched by its given name; but, optionally, the name
		 * may be modified by a prefix (see comprefixes[] above) to generate
		 * a different output type (see comretnew[z])...
		**/
	for (j=z=-1; j<0 && z<(int)(sizeof(comprefixes)/sizeof(comprefixes[0])); z++)
		for (i=0; j<0 && comhdefs[i].name!=NULL; i++) {
			buf[0] = 0;
			if (z>=0)  strncpy(buf,comprefixes[z],100);
			strncat(buf,comhdefs[i].name,200);
			if (OPTISNAME(com,buf))  j = i;
		}
	cmh = comhdefs[j<0?0:j];  z--;
	if (j<0 || !cmh.fhandle) {
		USERERROR("The requested handling function for command \"%s\" is not defined in the program!",OPTNAME(com));
		return -1;
	}
	if (z>=0 && cmh.ret==comretorig[z])  cmh.ret = comretnew[z];
	else if (z>=0)  {USERERROR("Wrong prefix \"%s\" for the handling function of \"%s\".",comprefixes[z],cmh.name);}
	DEBUG(CURDLEV," - handle for the command #%d \"%s\" applies here (%s%d)...\n",j,cmh.name,z>=0?"modified by #":"",z+1);
	
		/**
		 * After the command handling description (cmh) has been found, we must
		 * prepare all requested parameters for the execution.
		 * The types of parameters are given in cmh.par[], and their values are taken
		 * from the command parameters in com.
		 * If the command parameter string starts with FRGET_RESULT, then it describes
		 * the return value type, and it is not processed here but stored in sret.
		 * 
		 * A string-type parameter is passed unchanged.
		 * An integer-type parameter is passed from the numeric value in com.
		 * A matrix-type parameter may be either one matrix, or a matrix list (possibly 0).
		 * The matrix/ces is/are extracted from the current subtree of frames according
		 * to the parameter description string.
		 * The return list is checked to conform with the requested type in frame_getparammx_list().
		**/
	stpar = NULL;  sret = NULL;
	firstfrparam = NULL;	/* (this later points to the first frame of the first parameter) */
	firstfrparamlist = NULL;	/* (this later points to the frame-list of the first parameter) */
	jn = jf = 0;
	for (i=0; i<FCMAXPARAMS; i++) {
		s = OPTPARAM(com,i);
		if (s) if (s[0]==FRGET_RESULT && index(s+1,FRGET_DOWN)!=NULL) {
			sret = s+1;	/* (the return value type - not used here, but later) */
			s = NULL;
		}
		if (i>=cmh.numpar)  continue;	/* (no input parameters after cmh.numpar) */
		if (!s)  s = cmh.defpar[i];
		prstrings[i] = s;		/* (remembers the (possibly replaced) parameter for debugging) */
		
				/* a string parameter - passed unchanged */
		if (cmh.par[i]==FCPAR_STRING) {
			if (!s) {
				USERERROR("The requested string parameter \"!%s #%d\" is missing.",OPTNAME(com),i+1);
				s = "";
			}
			stpar = alist_append(stpar,MSTRDUP(s));
				/* an integer parameter - passed numeric */
		} else if (cmh.par[i]==FCPAR_INTEGER) {
			lk = OPTNPARAM(com,i);
			if (lk==OPTNOPARAM && (s?sscanf(s,"%ld",&lj):0)>0)
				lk = lj;	/* (possible default value) */
			if (lk==OPTNOPARAM) {
				USERERROR("The requested numeric parameter \"!%s #%d\" is missing.",OPTNAME(com),i+1);
				lk = 0;
			}
			npar[jn++] = lk;
			
				/* a matrix (list) parameter - getting the matrix list by OPTPARAM */
		} else if (FCPAR_MXGROUP(cmh.par[i])) {
			k = FCPAR_ONEREQ(cmh.par[i]);	/* (how many matrices are requested) */
			m = FCRET_SMARKREQ(cmh.ret);	/* (marking parent of son-list requested) */
			xx = ((cmh.ret==FCRET_REPLMXL&&!firstfrparamlist)?
				&firstfrparamlist:NULL);
			fmp = (void**)frame_getparammx_list(fr,i,jf,s,k,m,(framestruc***)xx);
			fmm[jf] = 1;	/* (remember that this parameter has matrices (not frames) */
			fmpar[jf++] = fmp;
				/* a frame (list) parameter - getting the matrix list by OPTPARAM */
		} else if (FCPAR_FRGROUP(cmh.par[i])) {
			k = FCPAR_ONEREQ(cmh.par[i]);
			n = FCPAR_MATREQ(cmh.par[i]);
			m = FCRET_SMARKREQ(cmh.ret);
			fmp = (void**)frame_getparamfr_list(fr,i,jf,s,k,m,n);
			if (cmh.ret==FCRET_REPLMXL && !firstfrparamlist)
				firstfrparamlist = alist_copy(fmp);
			fmm[jf] = 0;	/* (remember that this parameter holds frames */
			fmpar[jf++] = fmp;
		
		} else  {USERERROR("Unknown parameter %d type %d!",i,cmh.par[i]); return -1;}
	}
	if (cmh.numpar+1<FCMAXPARAMS) if (OPTPARAM(com,cmh.numpar+1))
		{USERERROR("Extra parameter to the command \"!%s ... %s\".",OPTNAME(com),OPTPARAM(com,cmh.numpar+1));}
	
	
		/**
		 * A header is printed first for commands with printed output, marked by clpr.
		 * 
		 * Now the command parameters are prepared:
		 * - stpar is the list of (generic) string parameters,
		 * - fmpar,jf is an array of the matrix/frame parameter lists
		 *   (matrices are given as copies, while frames are not copied),
		 * - npar,jn is an array of the integral parameters.
		 * Moreover, an array of list pointers outf[], and an array of integers outn[]
		 * are prepared for the command output, both cleared to zeros.
		 * 
		 * The return value retv of the command handle indicates what happened
		 * (in general) and what should be done with the results, see below.
		 * The matrices/frames in outf[] must be already returned as private
		 * copies - ready to be arranged within the (sub)frame tree.
		 * The return frames (if any) are saved for possible later recall
		 * (as parameter) in the static variable lastresults[lastresultix].
		**/
	DEBUG(CURDLEV,"  - call with parameters  %d ints, %c chars, %d lists\n",jn,stpar?'?':'0',jf);
	if (FCRET_PRINTREQ(cmh.ret) || frcom_verbose_out()>2) {
		frame_getoptionnum(fr,"prbrief",&lj,1,0);
		frcom_printbrief = (lj<100?lj:-1);	/* (possible request for a brief output) */
		clpr = (frcom_verbose_out()>0);
		if (clpr)  SOUTPUT(" vv===========================================================================================vv\n");
		if (clpr)  OUTPUT("\tOutput of the command \"!%s %.25s %.20s %.8s %.8s%s[%d]\":\n",
				OPTNAME(com),cmh.numpar>0?prstrings[0]:"",cmh.numpar>1?prstrings[1]:"",cmh.numpar>2?prstrings[2]:"",
				cmh.numpar>3?prstrings[3]:"",cmh.numpar>4?" ...":"",cmh.numpar);
		if (frcom_verbose_out()>4) {
			OUTPUT("\t\t%d list parameters: ",jf);
			for (i=0; i<jf; i++)  SOUTPUT(" %c %d,",fmm[i]?'m':'f',alist_getlength(fmpar[i]));
			SOUTPUT("\n");  OUTPUT("\t\t%d int parameters: ",jn);
			for (i=0; i<jn; i++)  SOUTPUT("  %ld,",npar[i]);
			SOUTPUT("\n");  OUTPUT("\t\t%d string parameters: ",stpar?alist_getlength(stpar):0);
			for (x=(void**)stpar; x?*x:0; x++)  SOUTPUT("  %s,",(char*)*x);
			SOUTPUT(".\n");
		}
		if (clpr)  SOUTPUT("\n");
	}
	for (i=0; i<FCMAXPARAMS; i++) { outf[i] = NULL;  outn[i] = -1; }	/* (empty output lists) */
					/* (saves the name(s) for the output frame list...) */
	nmret = firstfrparam? FRNAME(firstfrparam): FRNAME(fr);
	nmret = MSTRDUP(nmret?nmret:FR_DEFNONAME);
	for (y=firstfrparamlist; y?*y:0; y++) {
		s = FRNAME(*y)? FRNAME(*y):nmret;
		nmretlist = alist_append(nmretlist,MSTRDUP(s));
	}
	lastresultix_l = (lastresultix+1)%10;
	if (lastresults[lastresultix])  alist_free(lastresults[lastresultix]);
	lastresults[lastresultix] = NULL;	/* (clears the previous output frame list for storing the new one) */
	
	retv = (*cmh.fhandle)(stpar,fmpar,jf,npar,jn,cmh.ret,outf,outn);	/* (the command handle) */
	if (retv<0)  {USERERROR("An error %d occured while calling a handling function for \"!%s\".",retv,OPTNAME(com));}
	
	
		/**
		 * The return value retv of the command handle indicates what happened
		 * (in general) and what should be done with the results,
		 * in addition to the cmh.ret return type:
		 * - retv<0 indicates an error,
		 * - retv==0 indicates OK, but nothing is returned as a result,
		 * - retv>0 indicates OK, and a result is returned,
		 * - retv&FRES_STOREMX indicates that matrices should be stored,
		 * - retv&FRES_STOREFR indicates that frames should be stored,
		 * - retv&FRES_DELETEMARK requests deletion of (previously!) marked frames,
		 * - 
		 * (If the actions are not requested by the return value retv, they do not happen.)
		 * 
		 * The first result handling is a simple one - there is one output matrix
		 * that replaces the matrix in the first input frame firstfrparam
		 * (which is remembered in frame_getparamfr_list()).
		 * The second one stores the return number in outnum[0] to the first input frame,
		 * unless the return type is FCRET_NUMBER which is stored separately.
		 * The third one replaces matrices in a whole list of input frames.
		 * 
		**/
	if (cmh.ret==FCRET_REPLMX1 && retv>0 && (retv&FRES_STOREMX)) {
		if ((outf[0]?(!outf[0][0]||outf[0][1]):1) || outf[1] || !firstfrparam) {
			PROGERROR("Exactly one resulting matrix is expected for \"!%s\".",OPTNAME(com));
		} else {
			frame_setmatrix(firstfrparam,outf[0][0]);
			FRDELMARK(firstfrparam) = 0;
			if (cmh.fnformat) {	/* new frame name if the command requests that */
				snprintf(buf,490,cmh.fnformat, nmret,OPTPARAM(com,0),OPTPARAM(com,1),OPTPARAM(com,2),OPTPARAM(com,3));
				buf[490] = 0;  frame_setname(firstfrparam,buf);
			}
			if (retv&FRES_STOREREM)	/* (remember the output frame for possible later use) */
				lastresults[lastresultix_l] = alist_append(NULL,firstfrparam);
		}
	}
	if (cmh.ret!=FCRET_NUMBER && retv>0 && (retv&FRES_STORENUM)) {
		FRNUMBER(firstfrparam) = outn[0];
	}
	if (cmh.ret==FCRET_REPLMXL && retv>0 && (retv&FRES_STOREMX)) {
		for (i=0,y=firstfrparamlist; y?*y:0; i++,y++) {
			if (outf[0]?!outf[0][i]:1)  {PROGERROR("must give output matrices for all input frames! outf[0]=%p, i=%d",outf[0],i); break;}
			//******** should we allow skipping some output matrices?? (carefully with the length of the list!)
			frame_setmatrix(*y,outf[0][i]);	/* (no need to copy the matrix, outf[] elems are not disposed later) */
			FRDELMARK(*y) = 0;
			if (cmh.fnformat) {	/* new frame name if the command requests that */
				snprintf(buf,490,cmh.fnformat, nmretlist[i],OPTPARAM(com,0),OPTPARAM(com,1),OPTPARAM(com,2),OPTPARAM(com,3));
				buf[490] = 0;  frame_setname(*y,buf);
			}
		}
		if ((retv&FRES_STOREREM) && firstfrparamlist)	/* (remember the output frame for possible later use) */
			lastresults[lastresultix] = alist_copy(firstfrparamlist);
		if (outf[0]?outf[0][i]:0)  {PROGERROR("there are more output matrices %d than we may replace, why?",i+1);}
	}
		/**
		 * This part handles the result giving a yes/no list for printing (matrices or frames).
		 * The yes/no lists in outf[] must exactly match the parameter lists in fmpar[].
		 * (If an outf[] list is NULL, then the default value from outn[] is taken.)
		 * ...
		**/
	if (cmh.ret==FCRET_YESNOPRINT && retv>0) {
		for (z=0; z<2; z++) {
			if (z==0)  OUTPUT("Passed \"!%s\" :\n",OPTNAME(com));
			else  OUTPUT("NOT passed \"!%s\" :\n",OPTNAME(com));
			for (i=0; i<jf; i++)  if (fmpar[i] && outn[i]>=0)
			  for (x=fmpar[i],j=0; *x; x++,j++) {
				k = outf[i]? *(int*)(outf[i][j]): outn[i];
				k = (k!=ret_def? (k==ret_yes): (1));
				if ((k!=0)==z)  continue;
				s = FRMAGICTEST((framestruc*)(*x))? FRNAME((framestruc*)(*x)):"mx?";
				SOUTPUT("\t%p [%s] ;",*x,s);
			}
			SOUTPUT("\n");
		}
	}
		/**
		 * This part handles the result giving a yes/no list for filtering input frames.
		 * The yes/no lists in outf[] must exactly match the parameter lists in fmpar[].
		 * (If an outf[] list is NULL, then the default value from outn[] is taken.)
		 * A special value ret_def means that the current del-mark of the frame should be kept.
		 * Here, delete-marks are set according to the yes/no, and the frames are
		 * deleted later in FRES_DELETEMARK...
		**/
	if (FCRET_ISYESNO(cmh.ret) && retv>0 && (retv&FRES_STOREFR)) {
		for (i=0; i<jf; i++)
		  if (fmpar[i] && fmm[i]==0 && outn[i]>=0) {	/* (only frame-parameter lists that were marked with the default value in outn[]) */
			for (x=fmpar[i],j=z=0; *x; x++,j++) {
				if (!FRMAGICTEST((framestruc*)(*x)))  {PROGERROR("The parameter %d in the list %d is not a frame!",j,i); break;}
				m = (FCRET_ISYESNOY(cmh.ret)? ret_yes:ret_no);
				k = outf[i]? *(int*)(outf[i][j]): outn[i];
				if (k!=ret_def)  FRDELMARK((framestruc*)(*x)) = (k!=m);
				if (FRDELMARK((framestruc*)(*x)))  z++;
				else if (retv&FRES_STOREREM) {
					lastresults[lastresultix_l] = alist_append(lastresults[lastresultix],*x);
				}		/* (remember the resulting frame for possibe later use) */
			}
			DEBUG(CURDLEV,"Result: %d frames (of %d, param list %d) %smarked\n",
					z,alist_getlength(fmpar[i]),i,(retv&FRES_DELETEMARK)?"del-":"");
		}
	}
		/**
		 * Here the del-marked frames are removed from the frame tree.
		 * (That means input frames that were not requested for keeping, and possibly
		 * those frames that were filtered out by the FCRET_YESNOFR handle.)
		 * However, if FRES_DELETEMARK is not returned from the command, then
		 * nothing is deleted and all marks are cleared.
		**/
	FRDELMARK(fr) = 0;
	frame_deletemark_recur(fr,(retv&FRES_DELETEMARK)?1:0);
	
	
		/**
		 * Finally, we process command handles that return "real results" - one or a
		 * list of frames or matrices (i.e. not just printing, filtering or modification).
		 * The description of where the output should be stored is already in sret.
		 * Moreover, the list in outf[0] contains the (new copies of!) matrices or frames.
		 * We decide between matrices or frames based on (retv&FRES_STOREMX) flag.
		 * If a matrix list is on the output, then we must also create frames for them.
		 * 
		 * We give names to the output frames (if not given yet) using the pattern
		 * in cmh.fnformat, and comment them with a simple command description.
		 * Moreover, we prepare in buf[] a name for (possible) new frames created
		 * when storing the output - this name is the same for all.
		 * A then we call frame_storeres_list() to actually store the resulting frames.
		 * (The frames not stored there are returned, so they must be disposed after.)
		 * 
		**/
	if ((FCRET_LISTOUT(cmh.ret) || cmh.ret==FCRET_NUMBER) && retv>0 && 
			(retv&(FRES_STOREMX|FRES_STOREFR|FRES_STORENUM))) {
		if (FCRET_LISTOUT(cmh.ret)) if ((!outf[0]&&0) || outf[1])
			{USERERROR("Exactly one resulting list is expected for \"!%s\".",OPTNAME(com));}
		
		if (cmh.ret==FCRET_NUMBER) {	/* creating frame for the output number */
			fro = NULL;
			frn = new_frame(NULL);	FRLASTOPTI(frn) = 10000;
			FRNUMBER(frn) = outn[0];
			fro = alist_append(fro,frn);
		} else if (retv&FRES_STOREMX) {
			fro = NULL;		/* creating frames for the output matrices */
			for (x=outf[0]; x?*x:0; x++) {
				frn = new_frame(NULL);	FRLASTOPTI(frn) = 10000;
				frame_setmatrix(frn,*x);
				fro = alist_append(fro,frn);
			}
		} else {
			fro = (framestruc**)outf[0];
			outf[0] = NULL;		/* should not be freed twice (already as fro) */
		}
		if (!sret) {			/* the default output parameter description */
			sret = cmh.defret;
			if (!sret)  sret = (cmh.ret==FCRET_ONE||cmh.ret==FCRET_NUMBER? fcpar_oneret:fcpar_listret);
		}
		for (y=fro,i=1; y?*y:0; y++,i++) {	/* filling frame name (if required or not given) */
		  if (!FRNAME(*y) || cmh.fnformat) {
			snprintf(buf,490, cmh.fnformat?cmh.fnformat:FCDEFNFORMAT,
				(retv&FRES_STOREMXNM)&&FRMATRIX(*y)?EMNAME(FRMATRIX(*y)):nmret, i,
				OPTPARAM(com,0),OPTPARAM(com,1),OPTPARAM(com,2),OPTPARAM(com,3));
			buf[490] = 0;  frame_setname(*y,buf);
		  } if (!FRCOMMENT(*y)) {		/* filling frame comment if not existing */
			snprintf(buf,490,"fr #%d got by \'!%s %.8s %.6s %.4s%s\', to \'%.8s\'",
				i,OPTNAME(com),prstrings[0],cmh.numpar>1?prstrings[1]:"",cmh.numpar>2?prstrings[2]:"",cmh.numpar>3?"..":"",sret);
			buf[490] = 0;  FRCOMMENT(*y) = MSTRDUP(buf);
		} }
		for (y=fro; y?*y:0; y++)	/* (remembers the output frame for possible later use) */
			if (retv&FRES_STOREREM)
				lastresults[lastresultix_l] = alist_append(lastresults[lastresultix],*y);
		
		snprintf(buf,490,cmh.fnformat?cmh.fnformat:FCDEFNFORMAT, nmret,0,OPTPARAM(com,0),OPTPARAM(com,1),OPTPARAM(com,2),OPTPARAM(com,3));
		buf[490] = 0;
		DEBUG(CURDLEV-1,"Result: %d frames from \"!%s %s %s %s..\", stored in \"%s\".\n",i-1,OPTNAME(com),prstrings[0],cmh.numpar>1?prstrings[1]:"",cmh.numpar>2?prstrings[2]:"",sret);
		
		k = FCRET_ONEREQ(cmh.ret);
		fro = frame_storeres_list(fr,fro,sret,k,buf);
		if (fro) {	/* disposing the remaining frames, the list fro is freed after: */
			for (y=fro; *y; y++)  dispose_frame_recur(*y);
		}	/* (must not use dispose_alist_frams() since only part of the list is disposed!) */
		if (fro)  alist_free(fro);
		frame_deletemark_recur(fr,0);	/* (del-marks were locally used to mark the new frames) */
	}
	
		/**
		 * Here we complete the command handling routine.
		 * Mainly, we free the parameter lists that are no longer needed.
		**/
	if (clpr)  SOUTPUT(" ^^===========================================================================================^^\n");
	
			/* disposing the above collected parameters (they are locally copied) */
	if (stpar)  dispose_alist(stpar);
	if (firstfrparamlist)  alist_free(firstfrparamlist);
	for (i=0; i<jf; i++)  if (fmpar[i]) {
		if (fmm[i]==1)  dispose_alist_mats(fmpar[i]);
		else  alist_free(fmpar[i]);	/* (must not dispose frames, only matrices!) */
	}
	if (nmretlist)  dispose_alist(nmretlist);
	if (nmret)  FREE(nmret);
	/* (sret must not be freed!) */
	for (i=0; i<FCMAXPARAMS; i++)  if (outf[i]) {
		alist_free(outf[i]);		/* (unused output elements were already disposed) */
	}
	if (!eoc && retv>0) if (uerror_occured<0 || perror_occured<0) {
		USERERROR("An error (u%d, p%d) occured while calling a handling function for \"!%s\".",uerror_occured,perror_occured,OPTNAME(com));
		retv = (uerror_occured<0? uerror_occured:perror_occured);
	}
	return FRES_PROCSTRIP(retv);		/* the return value of the command handle (lower bits) */
}


















/******************	Getting command parameters	*******************/
/**************************************************************************/


#undef CURDLEV
#define CURDLEV         5


/**
 * These functions obtain the list of frames or the list of matrices requested by the
 * given parameter description pr from the subtree of the frame fr.
 * The value of ip is the index of this parameter (0...).
 * Defines FRGET_... are used to decrypt the parameter description string.
 * The frames collected as parameters may be marked for future deletion ("del" vs "keep").
 * 
 * If one==1 is given, then precisely one frame is expected from the parameter.
 * If one==0, then an empty list of frames is permitted, while this is not so if one==-1.
 * If needmx==1, then only frames containing matrices are collected.
 * If gmark==1, then also the parent of a son list is (possibly) marked for deletion.
 * See frame_getparamfr_recur() below for more comments.
 * If the parameter description is empty, then one of fcpar_oneframed[] or
 * fcpar_listframed[] is used.
 * 
 * The function frame_getparamfr_list() returns the list of requested frames, or NULL.
 * These frames are not copied, and so they MUST NOT be modified.
 * Possibly, the previous result-lists ^1,^2,... are picked here in frame_recoverlres().
 * 
 * It is better to call the function frame_getparammx_list() which returns a list
 * of matrix copies (however, the frame options are lost here).
 * The second function frame_getparammx_list() always returns a new list, even when
 * it is empty (so one may be sure there is no NULL returned here).
 * 
**/

framestruc**	frame_getparamfr_list(framestruc *fr, int ip, int ip0, char *pri, int one,\
								int gmark, int needmx) {
	framestruc	**gfr, **x;
	char		*pr = pri;
	int		i,j=0;
	
	if (pr) while (*pr==' ')  pr++;
	if (pr?!pr[0]:1) {	/* (default frame parameter description strings) */
		if (one==1)  pr = fcpar_oneframed[ip0];
		else  pr = fcpar_listframed[ip0];
		if (FRNUMSONS(fr)<=0 && fcpar_oneframens[ip0])  pr = fcpar_oneframens[ip0];
		if (one!=1 && fcpar_listframens[ip0])  if (alist_getlength(FRSONS(fr))>ip0)
			if (FRNUMSONS(FRSONS(fr)[ip0])<=0)  pr = fcpar_listframens[ip0];
		prstrings[ip] = pr;
		DEBUG(CURDLEV+1,"   (default parameter %d description \'%s\')\n",ip,pr);
	}
	DEBUG(CURDLEV+2,"Getting frame list for \"%s\".\n",pr);
	gfr = NULL;
	while (pr?*pr:0) {
				/* here we possibly pick the previous resulting list by ~^N: */
	  if (pr[0]==FRGET_PREVRES || pr[0]==FRGET_PREVRESOLD || pr[0]==FRGET_PREVRESKEEP) {
		i = pr[1]-'1';
		if (i<0 || i>9)  {USERERROR("Wrong \"previous result\" reference \"%s\".",pr); i=0;}
		i = (lastresultix+10-i)%10;
		gfr = frame_recoverlres(lastresults[i],needmx);
		if (pr[0]==FRGET_PREVRES)	/* ^N should be deleted after command */
			for (x=gfr; x?*x:0; x++)  FRDELMARK(*x) = 1;
		j = 1;
	  } else if (*pr==FRGET_DOWN) {
				/* this is the function that actually collects the parameters: */
		j = frame_getparamfr_recur(fr,&gfr,ip,pr+1,one,gmark,needmx,0);
	  } else {
	  	j = 0;
	  	if (*pr!=FRGET_CONCAT)  USERERROR("Invalid character in the parameter description \'%s\' at %d.",pr,(int)(pr-pri));
	  }
	  pr += 1+j;
	}
	if (gfr && !firstfrparam)	/* (we store the first frame parameter for some output storing) */
		if (gfr[0])  firstfrparam = gfr[0];
	
	if (!gfr && one==1)  {USERERROR("The requested parameter #%d \'%s\' cannot be found in the subtree.",ip,pri);}
	if (!gfr && one<0)  DEBUG(0,"Empty input list was found for parameter #%d \'%s\'.\n",ip+1,pri);
	if ((gfr?gfr[1]:0) && one==1)  {USERERROR("The requested one-parameter #%d \'%s\' has got more matrices.",ip,pri);}
	DEBUG(CURDLEV+2,"Got frame list %p of length %d, after scanning %d chars.\n",gfr,alist_getlength(gfr),(int)(pr-pri));
	
	return gfr;		/* the list of extracted frames is returned (do not modify!!!) */
}

int	frame_check_needmx(framestruc *fr, int ndmx) {
	
	if (!fr)  return 0;
	if (!FRMATRIX(fr))  {USERERROR("The requested frame %p has no matrix.",fr);  return 0;}
	if (ndmx>0)
	  if (FRPFINDEX(fr)!=pfield_curindex())  {USERERROR("The requested frame %p has matrix over a different pfield %d!=%d.",fr,FRPFINDEX(fr),pfield_curindex());  return 0;}
	return 1;
}

framestruc**	frame_recoverlres(framestruc **fl, int needmx) {
	framestruc	**x, **flo = NULL;
	
	for (x=fl; x?*x:0; x++) if (FRMAGICTEST(*x)) {
		if (needmx) if (frame_check_needmx(*x,needmx)==0)  continue;
		flo = alist_append(flo,*x);
	}
	DEBUG(CURDLEV-0,"Got frame list %p of length %d from the previous result.\n",flo,alist_getlength(flo));
	return flo;
}


/**
 * It is better to call the function frame_getparammx_list() which returns a list
 * of matrix copies (however, the frame options are lost here).
 * The function frame_getparammx_list() always returns a new list, even when
 * it is empty (so one may be sure there is no NULL returned here).
 * If remfrl is given, then the original list of frames is stored there.
**/

ematrix**	frame_getparammx_list(framestruc *fr, int ip, int ip0, char *pr, int one,
							int gmark, framestruc ***remfrl) {
	framestruc	**gfr, **x;
	ematrix		**eml, *em;
	
	gfr = frame_getparamfr_list(fr,ip,ip0,pr,one,gmark,1);
	eml = new_alist(6);
	if (gfr)  for (x=gfr; *x; x++) {	/* extracting matrix copies from the list of frames */
		em = frame_extractmatrix(*x);
		if (em)  eml = alist_append(eml,em);
	}
	if (remfrl)  *remfrl = gfr;
	else if (gfr)  alist_free(gfr);
	return eml;
}




/**
 * This is the recursive routine for frame_getparamfr_list().
 * The subframes requested by the (part of) parameter description pr are appended to the
 * list in *list, using the same additional inputs as in frame_getparamfr_list().
 * 
 * When this function finds a new opening '(', it dives into a recursion (possibly cycled).
 * When this function finds the closing ')' (the first opening '(' is already taken out),
 * it returns back, with additions to the *list.
 * See inside comments for a description of the algorithm.
**/

int	frame_getparamfr_recur(framestruc *fr, framestruc ***list,
				int ip, char *pr, int one, int gmark, int needmx, int depth) {
	framestruc	**gfr, **xs,**y;
	int		i,j,jj,k,l,p, cycle;
	char		buf[90];
	
	DEBUG(CURDLEV+2,"%*sCalling recur (depth %d) fr=%p, pr=\'%s\'\n",depth,"",depth,fr,pr);
	gfr = *list;
	xs = FRSONS(fr);
	cycle = 0;
	for (i=0; pr[i]; ) {
		/**
		 * The characters FRGET_UP, FRGET_DOWN are used to move up/down the frame tree.
		 * After FRGET_DOWN, an optional integer c meaning repetition of the next part
		 * below may be present.
		 * Such repeated part is sequentially applied to the sons starting from the
		 * current one, or to all remaining sons up to -c last ones if c<=0.
		 * (Applies to all remaining sons if c==0.)
		 * This is controlled by the variable cycle.
		 * Moreover, it is possible to move in the frame list directly to a frame of
		 * the given "name" which is enclosed in FRGET_BYNAME.
		**/
		if (cycle>0) {			/* cycling through sons: */
			if (xs?!xs[0]:1) {
				USERERROR("Missing #%d requested subframe of %p[%s] in \'%s\', depth %d.",ip+1,fr,FRNAME(fr),pr+i,depth);
				i += frame_paramcloseup(pr+i);  break;
			}
			j = frame_getparamfr_recur(xs[0],&gfr,ip,pr+i,one,gmark,needmx,depth+1);
			xs++;
			if (--cycle<=0)  i += j;
			continue;
		}
		if (pr[i]==FRGET_UP)		/* returning up in the tree */
			break;
		if (pr[i]==FRGET_CLOSE) {	/* returning all the way up in the tree */
			if (depth>0) i--;  break;
		}
		if (pr[i]==FRGET_DOWN) {	/* to go down - sets the cycle (above) */
			if (pr[i+1]==FRGET_BYNAME) {
				for (j=2; j<80 && pr[i+j]!=FRGET_BYNAME; j++)
					buf[j-2] = pr[i+j];	/* find a frame by (possible) given name */
				buf[j-2] = 0;
				for (jj=0; (xs?xs[jj]:0) && !FRISNAME(xs[jj],buf); jj++) ;
				if (xs?xs[jj]:0) {
					i += j;  xs += jj;
				} else {
					for (; pr[i+1] && pr[i+1]!=FRGET_UP && pr[i+1]!=FRGET_DOWN; i++) ;
					USERERROR("Cannot find the requested frame named \"%s\".",buf);
				}
			}
			i++;  j = 0;
			if (!pr[i])  break;
			if (sscanf(pr+i,"%d%n",&k,&j)>0) {	/* possible repetition of the frame */
				i += j;
				if (k<=0) {
					for (y=xs,l=0; y?*y:0; y++,l++) ;
					cycle = k+l;
				} else  cycle = k;
			} else  cycle = 1;
			if (cycle<0)  {USERERROR("Missing subframes of %p[%s] for \"negative\" repetition in \'%s\', depth %d.",fr,FRNAME(fr),pr+i-j,depth);}
			if (cycle<=0)		/* skip the descendants since 0-cycling is requested */
				i += frame_paramcloseup(pr+i);
			continue;
		}
		/**
		 * The characters FRGET_THIS.., FRGET_SONS.., etc request including the current
		 * frame or its (some?) subframes to the parameter list.
		 * Moreover, if needmx==1 is given, then the requested frames must contain matrices.
		 * The ...KEEP versions of the characters indicate that the current frame fr
		 * should be kept in the tree after execution.
		 * Otherwise, the frame fr is then (after execution) deleted;
		 * for gmark>0 even if only sons(!) of fr are included in parameters.
		**/
		if (pr[i]==FRGET_THIS || pr[i]==FRGET_THISKEEP) {
			p = 0;
			if (needmx) if (frame_check_needmx(fr,needmx)==0)  p = 1;
			if (!p)  gfr = alist_append(gfr,fr);
			if (pr[i]==FRGET_THIS)  FRDELMARK(fr) = 1;
			DEBUG(CURDLEV+2,"%*s     recur gets this frame %p.\n",depth,"",fr);
		
		} else  if (pr[i]==FRGET_SONS || pr[i]==FRGET_SONSKEEP) {
			y = FRSONS(fr);
			if (y?!y[0]:1) {
				if (one>=1)  {USERERROR("Missing #%d requested subframes of %p[%s] in \'%s\', depth %d.",ip+1,fr,FRNAME(fr),pr+i,depth);}
				if (one<0)  DEBUG(0,"Missing par#%d requested subframes of %p[%s] in \'%s\', depth %d.\n",ip+1,fr,FRNAME(fr),pr+i,depth);
			} else  for ( ; *y; y++) {
				if (needmx) if (frame_check_needmx(*y,needmx)==0)  continue;
				gfr = alist_append(gfr,*y);
				if (pr[i]==FRGET_SONS)  FRDELMARK(*y) = 1;
			}
			if (pr[i]==FRGET_SONS && gmark>0)  FRDELMARK(fr) = 1;
			if (gmark>0 && !firstfrparam)  firstfrparam = fr;
			DEBUG(CURDLEV+2,"%*s     recur gets sons of the frame %p.\n",depth,"",fr);
		}
		else  {USERERROR("Invalid character in #%d parameter description \'%s\'.",ip+1,pr+i); break;}
		i++;
	}
	if (!pr[i])  {--i; USERERROR("Unfinished #%d parameter description in \'%s\', depth %d.",ip+1,pr,depth);}
	
	*list = gfr;
	DEBUG(CURDLEV+2,"%*s   - recur (depth %d) finds %p, returns %d \'%s\'.\n",depth,"",depth,gfr,i+1,pr+i+1);
	return i+1;
}


/**
 * This function returns the number of characters that must be skipped in the given pr
 * to get one level up in a bracketed expression.
 * For use in frame_getparamfr_recur().
**/

int	frame_paramcloseup(char *pr) {
	int	i,j;
	for (i=0,j=1; pr[i] && j>0; i++)
		j += (pr[i]==FRGET_DOWN?1: (pr[i]==FRGET_UP?-1: 0));
	return i;
}











/******************	Command return values	*******************/
/******************************************************************/


#undef CURDLEV
#define CURDLEV         5


/**
 * Here the delete-marked frames are removed from the frame tree of fr.
 * However, if del==0 is given, then nothing is deleted, and all marks are cleared.
 * 
**/

void	frame_deletemark_recur(framestruc *fr, int del) {
	framestruc	**sc,**x;
	
	DEBUG(CURDLEV+3,"Call fr=%p, del=%d,  mark %d\n",fr,del,FRDELMARK(fr));
	if (del==0)  FRDELMARK(fr) = 0;
	if (FRDELMARK(fr)) {
		dispose_frame_recur(fr);	/* deleting this subtree */
	} else if (FRSONS(fr)) {
		sc = alist_copy(FRSONS(fr));	/* (must use a copy since deletion changes the son-list!) */
		for (x=sc; *x; x++)  frame_deletemark_recur(*x,del);
		alist_free(sc);
	}
}


/**
 * This function is called to store the resulting list of frames oul back to the frame
 * subtree of fr.
 * All frames in the list oul must be unique copies and not having parents.
 * (However, they may have their own frame subtrees.)
 * The description of where to store the results is in the string pr.
 * The description is similar to that in frame_getparamfr_list(), but there is an important
 * detail - all nonexisting sons requested by pr are automatically created (with a name nn).
 * These created frames may be later replaced by the resulting frames, or they remain in the tree.
 * If the frames are to be stored as ((s)), then they are arranged as sons of one added
 * (new!) parent frame.
 * The distinction between "keep" and "del" versions of parameters is ignored.
 * 
 * If the output parameter description is empty, then one of fcpar_oneret or fcpar_listret
 * is used.
 * If one==1 is given, then precisely one frame is stored.
 * If one==0, then only those resulting frames that have storage room are really stored,
 * and possible more storage room than frames is ignored.
 * While if one==-1, then there must be exactly the same number of frames as the storage room.
 * The return value is the sublist of frames that were not stored (to be deleted).
 * See in frame_storeres_recur() below for more comments.
 * 
**/

framestruc**	frame_storeres_list(framestruc *fr, framestruc **oul, char *pr, int one, char *nn) {
	int	j;
	
	if (!oul)  return NULL;
	if (!fr || !pr)  {PROGERROR("fr and pr must be given here!"); return NULL;}
	while (*pr==' ')  pr++;
	/*DEBUG(0,"pr=\'%s\"\n",pr);*/
	if (!strncmp(pr,"(s",2) || !strncmp(pr,"(S",2)) {
		USERERROR("Cannot store frames to sons of the current node with >(S), use >((0T)) instead.");
		pr = "((0T))";
	}
	if (pr[0]) {	/* this is the call that actually stores the parameters: */
		do {
			if (pr[0]!=FRGET_DOWN)  {USERERROR("Invalid character starting the result description \'%s\'.",pr);}
			j = frame_storeres_recur(fr,&oul,pr+1,one,nn,0);
			pr += 2+j;
		} while (pr[-1]==FRGET_CONCAT);
		if (pr[-1]!=0)  {USERERROR("Invalid character \'%s\' continues the result description.",pr-1);}
	}
	if ((oul?oul[0]:0) && one!=0)  {USERERROR("There is not enough room to store the resulting frame list. %p %p",oul,oul[0]);}
	return oul;
}


/**
 * This is the recursive part of frame_storeres_list().
 * It works similarly as frame_getparamfr_recur()....
 * 
**/

int	frame_storeres_recur(framestruc *fr, framestruc ***list, char *pr, int one, char *nn, int depth) {
	framestruc	**oul, **xs,**y,*nf;
	int		i,j,jj,k,l, xsi,cycle;
	char		buf[90], *nnn=NULL;
	
	DEBUG(CURDLEV+2,"%*sCalling recur (depth %d) fr=%p, pr=\'%s\'\n",depth,"",depth,fr,pr);
	if (!nnn)  nnn = nn;
	oul = *list;
	cycle = xsi = 0;
	for (i=0; pr[i]; ) {
		/**
		 * The characters FRGET_UP, FRGET_DOWN are used to move up/down the frame tree.
		 * Missing tree nodes are created on the fly, as requested (without matrices).
		 * After FRGET_DOWN, an optional integer c meaning repetition of the next part
		 * below may be present (variable cycle).
		 * Negative cycle means all but the last sons, or all but the last results.
		 * Notice that we must carefully access the son list of fr since it is changed
		 * as new result frames are stored!
		**/
		xs = FRSONS(fr);
		if (cycle>0) {			/* cycling through sons: */
			if (xs?!xs[xsi]:1) {
				nf = new_frame(fr);  FRDELMARK(nf) = 1;
				frame_setname(nf,nnn);  nnn = nn;
					/* (setting name with nnn works only for newly created frames!) */
				xs = FRSONS(fr);
				if (xs[xsi]!=nf)  {PROGERROR("Wrong son %p!=%p of the parent %p created.",xs[xsi],nf,fr);}
			}
			j = frame_storeres_recur(xs[xsi],&oul,pr+i,one,nn,depth+1);
			xsi++;
			if (--cycle<=0)  i += j;
			continue;
		}
		if (pr[i]==FRGET_UP)		/* returning up in the tree */
			break;
		if (pr[i]==FRGET_CLOSE) {	/* returning all the way up in the tree */
			if (depth>0) i--;  break;
		}
		if (pr[i]==FRGET_DOWN) {	/* to go down - sets the cycle (above) */
			if (pr[i+1]==FRGET_BYNAME) {
				for (j=2; j<80 && pr[i+j]!=FRGET_BYNAME; j++)
					buf[j-2] = pr[i+j];	/* find a frame by (possible) given name */
				buf[j-2] = 0;
				if (j<80)  i += j;
				else  for (; pr[i+1] && pr[i+1]!=FRGET_UP && pr[i+1]!=FRGET_DOWN; i++) ;
				
				for (jj=0; (xs?xs[jj]:0) && !FRISNAME(xs[jj],buf); jj++) ;
				if (xs?!xs[jj]:1)  nnn = buf;	/* if name not existent, then create a new node */
				else  xsi += jj;
			}
			i++;  j = 0;
			if (!pr[i])  break;
			if (sscanf(pr+i,"%d%n",&k,&j)>0) {	/* number-repetition scanned here */
				i += j;
				if (k<=0) {
					for (y=FRGET_ISGET(pr[i])?oul:xs, l=0; *y; y++,l++) ;
					cycle = k+l;	/* (counts against sons or output list) */
				} else  cycle = k;
			} else  cycle = 1;
			if (cycle<0)  {USERERROR("Missing subframes of %p[%s] for \"negative\" repetition in \'%s\', depth %d.",fr,FRNAME(fr),pr+i-j,depth);}
			if (cycle<=0)  i += frame_paramcloseup(pr+i);
			continue;
		}
		/**
		 * The characters FRGET_THIS, FRGET_SONS request storing the results.
		 * Resulting frames are stored sequentially, as they are available.
		 * If FRGET_THIS appears, then the next frame is stored.
		 * If FRGET_SONS appears, then a new frame is created that has all remaining
		 * result frames as its sons, and this one frame is then stored.
		 * If the current node fr is one of those created above, then it is replaced
		 * by the new frame nf (transparently to the rest).
		**/
		nf = NULL;
		if (pr[i]==FRGET_THIS || pr[i]==FRGET_THISKEEP) {
			if (one!=0 && (oul?!oul[0]:1))  {USERERROR("There is more storage room \'%s\' than available resulting frames.",pr);}
			else  if ((nf = oul[0])!=NULL)  oul++;
		} else  if (pr[i]==FRGET_SONS || pr[i]==FRGET_SONSKEEP) {
			if (one>=1 && (oul?!oul[0]:1))  {USERERROR("There is more storage room \'%s\' than available resulting frames.",pr);}
			nf = new_frame(NULL);  frame_setname(nf,nn);
			for (y=oul; y?*y:0; y++)  frame_setson(*y,nf);
			oul = y;
		}
		else  {USERERROR("Invalid character in parameter description \'%s\'.",pr+i); break;}
		if (nf) {
			if (FRPARENT(nf))  {PROGERROR("Cannot have %p parent %p here!",nf,FRPARENT(nf)); break;}
			DEBUG(CURDLEV+2,"Storing output frame %p to the parent %p, position %p.\n",nf,FRPARENT(fr),fr);
			frame_setson_wh(nf,FRPARENT(fr),fr);
			if (FRDELMARK(fr)>0 && FRNUMSONS(fr)<=0) {
				dispose_frame_recur(fr);	/* (replaced previously added frame) */
				fr = nf;
			}	/* now the new frame nf takes over the role of fr - transparently */
		}
		i++;
	}
	if (!pr[i])  {--i; USERERROR("Unfinished output parameter description in \'%s\', depth %d.",pr,depth);}
	*list = oul;
	return i+1;
}




































