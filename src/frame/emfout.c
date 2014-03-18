
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
 * A "frame" is the basic general data (and i/o) entity in the program.
 * Read more about a frame and its data in ../include/frame.h.
 * 
 * This file contains file-related function for the frame (mainly frame writing).
 * The input reading of a frame is done in emflex.l and emflexsup.c, not here.
 * There are also some short-help functions here...
 * 
**/



#include "macek.h"  
#include "fr.h"








/******************	Frame filenames		*******************/
/******************************************************************/


#undef CURDLEV
#define CURDLEV         6


/**
 * The following variables define paths for reading and for writing files,
 * and the default file extension.
 * Note that all the paths must be NULL-terminated!
 * Notice that frame_autosave_pref is also used in program consistency check!
 * See in frame.h...
 * 
**/

char	*frame_path_read[FR_PATH_LENGTH+2] = { FR_PATH_READ, NULL,NULL };

char	*frame_path_write[FR_PATH_LENGTH+2] = { FR_PATH_WRITE, NULL,NULL };

char	*frame_extension = FR_DEFEXTENSION;	/* (should contain '.') */

char	*frame_autosave_pref = FR_AUTOWRITEPREF;

int	file_overwrite = 9999;	/* whether files should be overwritten without warning
				   (allowing only at the first X write-path entries) */

int	safexec = 0;		/* flag for a restriction to safe execution only */












/******************	Frame output	*******************/
/**********************************************************/


#undef CURDLEV
#define CURDLEV         6



/**
 * This function writes the given frame fr to the file whose name is obtained from this
 * frame name (or from ancestors???, or noname).
 * If dir is given, then it is used as the directory for writing (a prefix to fname).
 * 
 * The options are written only if opt==1 (then @finherit and @finheritall apply).
 * The options listed in optnowrite[] (fropts.c) are never written to the file.
 * If onm==1, then the @name, @comment options are written first, before the matrix,
 * and they are not repeated among other options (onm>1 even for the root).
 * 
 * If rec==1, then the descendants of the frame fr are written to the same file as subframes.
 * (In the rec==1 case, only the root frame inherits options from its ancestors!)
 * 
 * The return value is 0 fo OK, and negative values for an error.
**/

int	frame_write_ext(framestruc *fr, int rec, int opt, int onm, char *dir, char *comment) {
	FILE		*file;
	char		*fn, fname[300], *fext, buf[100];
	int		ret = 0;
	
	fn = FRNAME(fr);
	fext = FREXTENS(fr)? FREXTENS(fr):frame_extension;
	if (!fn) {
		//******** read @name options from ancestors......
	}
	if (!fn)  fn = FR_DEFNONAME;
	if (dir)  sprintf(fname,"%.150s/",dir);
	else  fname[0] = 0;
	strncat(fname,fn,140);
	DEBUG(CURDLEV,"Going to write frame %p to \"%s\" (%s,%s,%s).\n",fr,fname,rec?"rec":"",opt?"opt":"",onm?"oname":"");
	
	file = fopen_path(fname,"w",frame_path_write,NULL,fext);
	ret = file?0:-1;
	if (ret>=0 && comment)  ret = frame_writecomment(file,FR_COMMENT,comment);
	snprintf(buf,48,"(Written in the (p)field %s.)\n",pfield_curname());  buf[49] = 0;
	if (ret>=0)  ret = frame_writecomment(file,FR_COMMENT,buf);
	if (ret>=0)  fprintf(file,"\n");
	if (ret>=0)  ret = frame_write_recur(file,fr,rec,opt,1,onm,0);
	snprintf(buf,88,"\n(Successfully finished writing frame %s in %s.)\n",fname,pfield_curname());  buf[89] = 0;
	if (ret>=0)  ret = frame_writecomment(file,FR_COMMENT,buf);
	if (file)  ret = fclose(file)<0? -1:ret;
	if (ret<0) {
		USERERROR("Error writing frame %p to \"%s\": %s.",fr,fname,strerror(errno));
	} else {
		DEBUG(CURDLEV,"- Written %p OK.\n",fr);
	}
	return ret;
}


/**
 * This is the recursive part of frame_write().
 * It writes the matrix and the options of the frame fr, and then (if rec==1) it writes
 * the sons of fr as subframes (to the same file).
 * 
 * The options are printed only if opt==1.
 * If inh==1 is given, then it also writes the inherited options from ancestors.
 * (In the rec==1 case, only the root frame inherits options from its ancestors!)
 * If onm==1, then the @name, @comment options are written first, before the matrix,
 * and they are not repeated among other options.
 * The options listed in optnowrite[] (fropts.c) are never written to the file.
 * 
**/

int	frame_write_recur(FILE *file, framestruc *fr, int rec, int opt, int inh, int onm, int depth) {
	ematrix		*e;
	framestruc	**x;
	optionstruc	**ol;
	int		i,j,k,l, ret,xpf, onmw=0;
	char		*buf,*sb, **inhn,**inhna, **inhl, **y,**yy;
	
	DEBUG(CURDLEV+2,"%*sCalling recur,  fr=%p, rec=%d, opt=%d, inh=%d, onm=%d\n",depth*3,"",fr,rec,opt,inh,onm);
	buf = MMALLOC(5100*sizeof(buf[0]));
	
	if (onm && depth+onm>1) {	/* writing options @name and @comment: */
		onmw = 1;	/* (not to repeat these later in local options) */
		if (FRNAME(fr))  fprintf(file,"%sname \"%s\"\n",FR_OPTION,FRNAME(fr));
		if (FRCOMMENT(fr))  fprintf(file,"%scomment \"%s\"\n\n",FR_OPTION,FRCOMMENT(fr));
	}
	
	e = FRMATRIX(fr);		/* writing the frame matrix: */
	if (e && FRPFINDEX(fr)!=pfield_curindex()) {
		DEBUG(CURDLEV-0,"The matrix of %s is written in different pfield %d!=%d.\n",FRNAME(fr)?FRNAME(fr):"",FRPFINDEX(fr),pfield_curindex());
		xpf = pfield_curindex();
		pfield_switchto_fast(FRPFINDEX(fr));
		fprintf(file,"%sinputpf \"%s\"\n",FR_OPTION,pfield_curname());
	} else  xpf = -1;		/* (the matrix should be written in the right pfield!) */
	if (e)  for (i=0; i<ROWSM(e); i++) {
		for (j=0; j<COLSM(e); j++) {
			pfield_printvalue_formal(buf,5000,EXPM(e,i,j),SIGNM(e,i,j));
			fprintf(file,"\t%s",buf);
		}
		fprintf(file,"\n");
		if (ferror(file))  break;
	} else {
		fprintf(file,"%s\t(no matrix)\n",FR_COMMENT);
	}
	if (xpf>=0)  pfield_switchto_fast(xpf);
	fprintf(file,"\n");
	
	if (opt) {			/* writing this (and possibly inherited) frame options: */
		if (inh && FRPARENT(fr)) {
			inhn = frame_getoptionval_all(fr,"finherit");
			inhna = frame_getoptionval_all(fr,"finheritall");
		} else  inhn = inhna = NULL;
		if (inhn)  for (y=inhn; *y; y++) {
				/* the options local to fr are not inherited (for single) */
			for (i=0,ol=FROPTIONS(fr); ol?*ol:0; ol++)
				if (OPTISNAME(*ol,*y))  i = 1;
			if (i)  continue;
			DEBUG(CURDLEV+2,"%*sInherited option \"%s\"...\n",depth*3,"",*y);
			inhl = frame_getoptionprint(FRPARENT(fr),*y);
			for (yy=inhl; yy?*yy:0; yy++)
				fprintf(file,"%s\n",*yy);
			if (inhl)  dispose_alist(inhl);
		}
		if (inhna)  for (y=inhna; *y; y++) {
			DEBUG(CURDLEV+2,"%*sInherited-all option \"%s\"...\n",depth*3,"",*y);
			inhl = frame_getoptionprint_all(FRPARENT(fr),*y);
			for (yy=inhl; yy?*yy:0; yy++)
				fprintf(file,"%s\n",*yy);
			if (inhl)  dispose_alist(inhl);
		}
		for (ol=FROPTIONS(fr); ol?*ol:0; ) {
			DEBUG(CURDLEV+2,"%*sLocal option \"%s\"...\n",depth*3,"",OPTNAME(*ol));
			ol = option_print_to(ol,buf,5000);
			sb = buf+strlen(FR_OPTION);
				/* this code checks not to write name again, and not to write those options in optnowrite[]: */
			if (onmw && strncmp(sb,"name",4)==0)  continue;
			if (onmw && strncmp(sb,"comment",7)==0)  continue;
			for (j=0; optnowrite[j]; j++)
				if (strncmp(sb,optnowrite[j],strlen(optnowrite[j]))==0)  break;
			if (optnowrite[j]) {
				for (k=l=0, y=inhn; k<2; k++, y=inhna)
				  for (; (y?*y:0) && !l; y++)
					if (strncmp(sb,*y,strlen(*y))==0)  l = 1;
				if (!l)  continue;	/* (skipping no-write options unless inherited) */
			}
			fprintf(file,"%s\n",buf);
			if (ferror(file))  break;
		}
		if (inhn)  dispose_alist(inhn);
		if (inhna)  dispose_alist(inhna);
	}
	
	if (rec && FRSONS(fr))		/* writing the frame descendants recursively */
	  for (x=FRSONS(fr); *x; x++) {
		fprintf(file,"\n%s%*s## subframe lev %d.%d [%s] ##\n",FR_COMMENT,3*depth+1,"",
				depth+1,(int)(x-FRSONS(fr))+1,FRNAME(*x)?FRNAME(*x):"");
		fprintf(file,"%s\n",FR_SUBFRAME);
		ret = frame_write_recur(file,*x,rec,opt,0,onm,depth+1);
		fprintf(file,"%s\n",FR_ENDFRAME);
		fprintf(file,"%s%*s#  endframe lev %d.%d  #\n",FR_COMMENT,3*depth+1,"",
				depth+1,(int)(x-FRSONS(fr))+1);
		if (ret<0 || ferror(file))  break;
	}
	FREE(buf);
	return (ferror(file)?-1:0);
}


/**
 * This function writes the given one frame similarly as frame_write_ext() above,
 * but the output file will additionally contain comments on the written matroid
 * structure - like bases, small flats, connectivity, aut group orbits.
 * Do not be surprised when the function takes long time, that is because of
 * computation of the matroid aut group.
**/

int	frame_write_matcom(framestruc *fr, int opt, int onm, char *dir, char *comment) {
	ematrix		*ee;
	char		*buf=NULL;
#	define		XBUFMAX 30000
	int		k;
	
	ee = frame_extractmatrix(fr);
	if (!ee) {
		USERERROR("Missing matrix in the frame when writing with matroid comments.");
		return frame_write_ext(fr,0,opt,onm,dir,comment);
	}
	buf = MALLOC(XBUFMAX*sizeof(buf[0]));
	strncpy(buf,comment,XBUFMAX/3);
	strcat(buf,"\n\n");
	ematrix_printmore_to(ee,3,buf+strlen(buf),XBUFMAX/2);
			/* (need level 3 printing to compute the automorhpism group) */
	buf[XBUFMAX-5] = 0;  strcat(buf,"\n");
	
	k = frame_write_ext(fr,0,opt,onm,dir,buf);
	if (ee)  dispose_ematrix(ee);
	if (buf)  FREE(buf);
	return k;
}



/**
 * This function prints the string comment into the file fp,
 * prefixing each line with the string comchar.
**/

int	frame_writecomment(FILE *fp, char *comchar, char *comment) {
	char	*sc,*sc0,*ss;
	
	if (!comment || !fp )  return 0;
	sc = sc0 = MSTRDUP(comment);
	while (sc[0]) {
		if ((ss = strchr(sc,'\n'))!=NULL)  *ss++ = 0;
		else  ss = sc+strlen(sc);
		fprintf(fp,"%s%s\n",comchar,sc);
		sc = ss;
	}
	FREE(sc0);
	return (ferror(fp)?-1:0);
}
















/******************	Short help messages	*******************/
/******************************************************************/


#undef CURDLEV
#define CURDLEV         6


/**
 * These functions print short help messages on the frame options and commands.
 * The names and descriptions are taken from arrays optdescdefs, comhdefs which are
 * defined in fropts.c, frcoms.c.
 * The message is slightly longer if ch>0.
 * 
**/

void	frame_foptionhelp(int ch, FILE *fo) {
	int	i,j,jj,k;
	
	fprintf(fo,"Short help on %s %s frame-options:\n",PROGNAME,PROGVER);
	
	for (i=0; optdescdefs[i].name; i++) {
		for (j=0; optdescdefs[j].name; j++) {
			for (jj=k=0; optdescdefs[jj].name; jj++)
				if (strcasecmp(optdescdefs[j].name,optdescdefs[jj].name)>0)  k++;
			if (k==i)  break;
		}
		if (optdescdefs[j].description[0]=='-')  continue;
			/* (starting '-' means a deprecated option) */
		if (ch>0) {
			fprintf(fo,"  %5s \t[par#%d]\n",optdescdefs[j].name,optdescdefs[j].numpar);
			fprintf(fo,"\t(%s)\n",optdescdefs[j].description);
		} else {
			k = index(optdescdefs[j].description,';')-index(optdescdefs[j].description,':');
			fprintf(fo,"  %5s \t[par#%d] \t%.*s\n",optdescdefs[j].name,optdescdefs[j].numpar,
				k,index(optdescdefs[j].description,':'));
		}
	}
}

void	frame_fcommandhelp(int ch, FILE *fo) {
	int		i,z, j,jj,k;
	static int	rettp[] = { FCRET_NOTHING, FCRET_PRINTED, FCRET_PRINTED_F, FCRET_YESNOFR, FCRET_REPLMX1, FCRET_ONE, 111111 };
	static char	*retst[] = { "nothing", "printed", "print/filter", "filter", "modify matrix", "new list", NULL };
	
	fprintf(fo,"Short help on %s %s frame-commands:\n",PROGNAME,PROGVER);
	fprintf(fo," (Give connected matrices with all determinants defined on input.)\n");
	
	for (z=0; z<(int)(sizeof(rettp)/sizeof(rettp[0]))-1; z++) {
	  fprintf(fo,"%s- commands of return value type  \'%s\':\n",ch>0?"\n":"",retst[z]);
	  if (rettp[z]==FCRET_PRINTED_F) {
	  	fprintf(fo,"  (use prefixes \"%s\" or \"%s\" change printing to filtering)\n",
	  			FRES_FILTPREFIX,FRES_INVFILTPREFIX);
	  	fprintf(fo,"  (use prefixes \"%s\" or \"%s\" to change printing to remember-list)\n",
	  			FRES_REMPREFIX,FRES_INVREMPREFIX);
	  }
	  for (i=0; comhdefs[i].name; i++) {
		for (j=0; comhdefs[j].name; j++) {
			for (jj=k=0; comhdefs[jj].name; jj++)
				if (strcasecmp(comhdefs[j].name,comhdefs[jj].name)>0)  k++;
			if (k==i)  break;
		}
		if (comhdefs[j].ret<rettp[z] || comhdefs[j].ret>=rettp[z+1])
			continue;
		if (comhdefs[j].description[0]=='-')  continue;
			/* (starting '-' means a deprecated option) */
		if (ch>0) {
			fprintf(fo,"    %4s \t[par#%d] %s\n",comhdefs[j].name,comhdefs[j].numpar,
				FCRET_LISTOUT(comhdefs[j].ret)?">[ret]":"");
			fprintf(fo,"\t(%s)\n",comhdefs[j].description);
		} else {
			k = index(comhdefs[j].description,';')-index(comhdefs[j].description,':');
			fprintf(fo,"    %4s \t[par#%d] %6s\t%.*s\n",comhdefs[j].name,comhdefs[j].numpar,
				FCRET_LISTOUT(comhdefs[j].ret)?">[ret]":"",
				k,index(comhdefs[j].description,':'));
		}
	  }
	}
	if (ch>0)  fprintf(fo,"\n");
	fprintf(fo,"Write complex commands to a file \"procx\" (a procedure), and call them \"&procx par1 ...\".\n");
	fprintf(fo,"  The procedure parameters are then text-substituted for \"$param1\" etc.\n");
}

void	frame_framehelp(int ch, FILE *fo) {
	int	i;
	
	fprintf(fo,"Short help on %s %s input frames:\n",PROGNAME,PROGVER);
	if (ch>0)  fprintf(fo,"\n");
	fprintf(fo,"  A \"frame\" is the basic general data entity in the program.\n");
	fprintf(fo,"  One frame keeps one matrix (of a matroid), and a list of options and commands.\n");
	fprintf(fo,"  Frames may have lists of descendant (sub)frames forming a tree-like structure.\n");
	fprintf(fo,"Frame search paths (read / write \'*\' means not create dirs):\n  <");
	for (i=0; frame_path_read[i]; i++)  fprintf(fo," %s",frame_path_read[i]);
	fprintf(fo,"\n  >");
	for (i=0; frame_path_write[i]; i++)  fprintf(fo," %s",frame_path_write[i]);
	fprintf(fo,"\n  Default extension:  %s\n",frame_extension);
	fprintf(fo,"The frame input consists of lines of the following types:\n");
	fprintf(fo,"  Start a matrix line with a space, then write space-separated matrix entries.\n");
	fprintf(fo,"  Start an option / command line with \"%s\" / \"%s\".\n",FR_OPTION,FR_COMMAND);
	fprintf(fo,"  Write your comments after starting \"%s\"...\n",FR_COMMENT);
	fprintf(fo,"  Start a subframe with line \"%s\" or \"{\", end it with \"%s\" or \"}\".\n",FR_SUBFRAME,FR_ENDFRAME);
	fprintf(fo,"  Include files \"<fn\", include files as subframes in a line \"{ fn1 fn2 ... }\".\n");
	fprintf(fo,"  On command-line or in subframe-line, each \';\' is replaced by a newline.\n");
	fprintf(fo,"How values and parameters are given to commands (or options):\n");
	fprintf(fo,"  The input parameters for commands are addressed in the frame tree using\n");
	fprintf(fo,"   bracketing - \"(()((%c)))\" means the son of the second son of the root.\n",FRGET_THIS);
	fprintf(fo,"   Here \'%c\' or \'%c\' picks one frame, \'%c\' or \'%c\' picks all sons.\n",FRGET_THIS,FRGET_THISKEEP,FRGET_SONS,FRGET_SONSKEEP);
	fprintf(fo,"   (Use lower-case if the picked frames should be deleted - only for certain commands.)\n");
	fprintf(fo,"   You may close all open brackets in the parameter with single \"%c\".\n",FRGET_CLOSE);
	fprintf(fo,"  Moreover, \"%c1\" alone picks the previous command output, \"%c2\" one before, etc.\n",FRGET_PREVRESKEEP,FRGET_PREVRESKEEP);
	fprintf(fo,"   Alternative \"%c1\", \"%c2\", etc indicate that the frames should be deleted after.\n",FRGET_PREVRES,FRGET_PREVRES);
	fprintf(fo,"  The place to store output frame is addressed similarly, like \"%c((%c))\".\n",FRGET_RESULT,FRGET_THIS);
	fprintf(fo,"   When storing, all non-existent tree nodes are automatically created.\n");
	fprintf(fo,"  A shurtcut like \"(...(5xxx)...)\" repeats 5-times the part \"xxx\" which is as above.\n");
	fprintf(fo,"  An expression \"(...(%cname%c)...)\" is for a direct access to the son \"name\".\n",FRGET_BYNAME,FRGET_BYNAME);
	fprintf(fo,"   You may combine these like \"(...(%cname%c3xxx)...)\".\n",FRGET_BYNAME,FRGET_BYNAME);
	fprintf(fo,"  Macros like \"$macr\" are expanded to the value of option \"@%smacr xxx..\".\n",FR_SUBSTPREFIX);
	fprintf(fo,"Play with the command \"!move in out\" to learn more about parameter addressing...\n");
	fprintf(fo,"(Do not forget that many characters have special meaning in your shell, so protect them.)\n");
	//fprintf(fo,"\n");
}






































