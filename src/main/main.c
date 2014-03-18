
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
 * This file contains the main starting function of the package,
 * and code handling the command-line options given to the program.
 * 
 * You find nothing much interesting here...
 * 
 * 
**/




#include "macek.h"








/******************	Main program start	*******************/
/******************************************************************/



/**
 * This is the main function of the program.
 * It task is to process the command-line options and inputs,
 * scan the given frames, and process commands in the frame-tree.
 * The resulting tree is automatically saved after a long computation.
 * The function macek_begin() is called before anything else, and macek_final()
 * is called at the end.
 * 
**/

char**  macek_readopts(int argc, char *argv[]) ;
int     macek_begin(void) ;
int     macek_final(void) ;
int	begin_check(void) ;
int	final_check(void) ;

int  main(int argc, char *argv[]) {
	framestruc      *fr=NULL;
	char		**inparg, buf[300];
	long		tim;
	int		ret = 0;
	
	macek_begin();
	if (argc>2)  SDEBUG(2,"  %s %s (%s) starting...\n",PROGNAME,PROGVER,PROGDATE);
	inparg = macek_readopts(argc,argv);
	if (argc<=2)  SDEBUG(2,"  %s %s (%s) starting...\n",PROGNAME,PROGVER,PROGDATE);
	if (!pfield_curname())  pfield_switchto("gf2");
#ifndef FASTPROG
	if (begin_check()<0) {	/* possible internal consistency check, randomized inside */
		PROGERROR("Internal error occured during program consistency check, exitting!");
		exit (-1);
	}
#endif
		/**
		 * All non-option arguments are scanned as input frames, the first one
		 * as the root frame.
		 * Then the frame commands are processed, in reverse depth-first.
		**/
	if (inparg)  fr = frame_doinput_list(inparg);
	SDEBUG(4,"\n====\n");
	tim = TIME();
	if (fr)  ret = frame_processcommands(fr);
	
	if (fr && TIME()-tim>50 && frame_autosave_pref) {
		snprintf(buf,290,"%smacek-out-%d",frame_autosave_pref,(int)(RANDOM()%100));
		frame_setname(fr,buf);
		frame_write_tree(fr,"Automatically saved the frame-tree after long computation...\n");
		OUTPUT("Automatically saved the frame-tree to \"%s\"\n",FRNAME(fr));
	}
	if (fr)  dispose_frame_recur(fr);
	if (inparg)  dispose_alist(inparg);
	
	macek_final();		/* (calls final check inside) */
	return ret;
}




/**
 * Reading command-line options for the program.
 * If a help message or version info is requested, or if an unknown option is given,
 * then the program immediately ends here.
 * 
 * All given non-option arguments are returned in a string list from macek_readopts().
**/

char*	macek_readoptval(int *ix, int argc, char *argv[], int next) {
	char	*val = NULL;
	int	j,i = *ix;
	
	if (argv[i][1]!='-') {	/* simple -X options, read immediately next chars: */
		if (argv[i][1]?argv[i][2]:0)
			val = argv[i]+2;
	} else {		/* long --Xxx options, read after '=' or next: */
		for (j=2; argv[i][j] && argv[i][j]!='='; j++) ;
		if (argv[i][j]?argv[i][j+1]:0)
			val = argv[i]+j+1;
	}
	if (!val && i+1<argc && next) {	/* taking the next input as an option value */
		val = argv[i+1];
		(*ix)++;
	}
	return val;
}


char**	macek_readopts(int argc, char *argv[]) {
	int	i,j,jj;
	char	**plist=NULL, *oval, *x;
	FILE	*lpout = stdout;
	
	if (argc<2) {
		fprintf(lpout,"***   %s %s (%s)   ***\n",PROGNAME,PROGVER,PROGDATE);
		fprintf(stderr," Try \"%s -h\" or \"%s -H\" for help.\n",argv[0],argv[0]);
		exit(0);
	}
	
	for (i=1; i<argc; i++) {	/* cycling all command-line options... */
		
#define	PRINTDESC()	fprintf(lpout,"Program for practical computations with represented matroids.\n");\
		fprintf(lpout,PHCOPYRIGHT ".\n");\
		fprintf(lpout,"  This program is free software; you can redistribute it and/or modify\n");\
		fprintf(lpout,"  it under the terms of the GNU General Public License, version 2.\n");\
		fprintf(lpout,"  (See the file COPYING for details.)\n")

	if (strncmp(argv[i],"-h",2)==0 || strncmp(argv[i],"--h",3)==0) {
		PRINTDESC();
		fprintf(lpout,"\nUse:\t    %s [-ghHprRsSwtTvx] frame [subframe ...]\n",argv[0]);
		fprintf(lpout,"  -g N\tincrease/decrease verbosity level of debug messages (--debug=),\n");
		fprintf(lpout,"  -h\tthis quick help,\n");
		fprintf(lpout,"  -H [H][pfco]\tmore help on `pfields, `frames, `commands, and `options,\n");
		fprintf(lpout,"  -p pfield\tswitch to this partial field (--pfield=),\n");
		fprintf(lpout,"  -r dir/\tappend this directory to the read path (--readpathapp=),\n");
		fprintf(lpout,"  -R dir/\tinsert first this directory to the read path (--readpathins=),\n");
		fprintf(lpout,"  -s \tonly safe execution is allowed - no shell access (--safe),\n");
		fprintf(lpout,"  -S \tonly very safe exect. - no shell access, no absol. paths (--safest),\n");
		fprintf(lpout,"     \t\tno file overwrite except for the first path entry,\n");
		fprintf(lpout,"  -w dir/\tuse this directory for writing files (--writedir=),\n");
		fprintf(lpout,"  -t dir/\tuse this temp directory for auto-save (--temps=),\n");
		fprintf(lpout,"  -T [qsmh]\t4-digit time in output, opt. minutes or hours (--time),\n");
		fprintf(lpout,"  -v\tprint the version info (--version),\n");
		fprintf(lpout,"  -x ext\tuse ext as the extensions (.mck) of frame files (--extension=).\n");
		//fprintf(lpout,"  - and more to come...\n");
		fprintf(lpout,"\nSee  " MACEKWEB ".\n");
		//fprintf(lpout,"\n");
		exit(0);
		
	} else if (strncmp(argv[i],"-H",2)==0 || strncmp(argv[i],"--H",3)==0) {
		oval = macek_readoptval(&i,argc,argv,0);
		if (!oval)  oval = "pfco";
		if (strcmp(oval,"H")==0)  oval = "Hpfco";
		j = (index(oval,'H')!=NULL);
		
		if (index(oval,'p')!=NULL) {
			pfield_pfieldhelp(j,lpout);
			fprintf(lpout,"\n");
		}
		if (index(oval,'f')!=NULL) {
			frame_framehelp(j,lpout);
			fprintf(lpout,"\n");
		}
		if (index(oval,'o')!=NULL) {
			frame_foptionhelp(j,lpout);
			fprintf(lpout,"\n");
		}
		if (index(oval,'c')!=NULL) {
			frame_fcommandhelp(j,lpout);
			fprintf(lpout,"\n");
		}
		if (!j)  fprintf(lpout,"Try \"%s -HH [pcfo]\" for even more description.\n",argv[0]);
		exit(0);
		
	} else if (strncmp(argv[i],"-v",2)==0 || strncmp(argv[i],"--ver",5)==0) {
		fprintf(lpout,"***   %s %s (%s)   ***\n",PROGNAME,PROGVER,PROGDATE);
		PRINTDESC();
		exit(0);
		
	} else if (strncmp(argv[i],"-p",2)==0 || strncmp(argv[i],"--pf",4)==0) {
		j = pfield_switchto(macek_readoptval(&i,argc,argv,1));
		if (j<0)  {USERERROR("Wrong value for '-p'."); exit(-1);}
		
	} else if (strncmp(argv[i],"-g",2)==0 || strncmp(argv[i],"--deb",5)==0) {
		j = jj = 0;  oval = macek_readoptval(&i,argc,argv,1);
		if (oval)  sscanf(oval,"%d%n",&j,&jj);
		if (oval?(jj>=(int)strlen(oval)-1):0)
			printlev += j;
		else  {USERERROR("Missing parameter for '-g'."); exit(-1);}
	
	} else if (strncmp(argv[i],"-T",2)==0 || strncmp(argv[i],"--tim",5)==0) {
		oval = macek_readoptval(&i,argc,argv,1);
		timeprecout = 4;
		if (oval?oval[0]=='m':0)  timeprecdiv = 60;
		if (oval?oval[0]=='h':0)  timeprecdiv = 3600;
		if (oval?oval[0]=='q':0)  timeprecout = 0;
	
	} else if (strncmp(argv[i],"-r",2)==0 || strncmp(argv[i],"-R",2)==0 ||
					strncmp(argv[i],"--readpath",10)==0) {
		oval = macek_readoptval(&i,argc,argv,1);
		if (safexec<10) {
		if (strncmp(argv[i],"-r",2)==0 || strncmp(argv[i],"--readpathapp",13)==0) {
			for (j=0; j<FR_PATH_LENGTH-2 && frame_path_read[j]; j++) ;
			frame_path_read[j+1] = NULL;
		} else {
			for (j=FR_PATH_LENGTH-1; j>0; j--)
				frame_path_read[j] = frame_path_read[j-1];
			frame_path_read[FR_PATH_LENGTH-1] = NULL;
		}
		frame_path_read[j] = MSTRDUP(oval);
		} else {USERERROR("Not allowed to change the path in a safe mode!"); exit(-1);}
		
	} else if (strncmp(argv[i],"-w",2)==0 || strncmp(argv[i],"--wri",5)==0) {
		oval = macek_readoptval(&i,argc,argv,1);
		if (safexec<10) {
		for (j=FR_PATH_LENGTH-2; j>0; j--)
			frame_path_write[j] = frame_path_write[j-1];
		frame_path_write[0] = MSTRDUP(oval);
		frame_path_write[FR_PATH_LENGTH-1] = NULL;
		} else {USERERROR("Not allowed to change the path in a safe mode!"); exit(-1);}
	
	} else if (strncmp(argv[i],"-S",2)==0 || strncmp(argv[i],"--safest",8)==0) {
		safexec = 11;	/* here we check all file names (incl. the path),
				   and allow to overwrite only at the first path entry */
		file_overwrite = 1;
		frame_autosave_pref = NULL;  frame_path_write[1] = NULL;
		fprintf(lpout,"\t***  Only very safe Macek execution is allowed.\n");
	} else if (strncmp(argv[i],"-s",2)==0 || strncmp(argv[i],"--safe",6)==0) {
		safexec = 1;
		fprintf(lpout,"\t***  Only safe Macek execution is allowed.\n");
		
	} else if (strncmp(argv[i],"-t",2)==0 || strncmp(argv[i],"--temp",6)==0) {
		oval = macek_readoptval(&i,argc,argv,1);
		if (safexec<10) {
		x = MMALLOC(5+strlen(oval)*sizeof(oval[0]));
		strcpy(x,oval);  jj = strlen(oval);
		if (jj==1 && oval[0]=='-')  x = NULL;
		else if (jj>0 && x[jj-1]!='/')  strcat(x,"/");
		frame_autosave_pref = x;
		} else {USERERROR("Not allowed to change the path in a safe mode!"); exit(-1);}
	
	} else if (strncmp(argv[i],"-x",2)==0 || strncmp(argv[i],"--ext",5)==0) {
		oval = macek_readoptval(&i,argc,argv,1);
		x = MMALLOC(5+strlen(oval)*sizeof(oval[0]));  x[0] = 0;
		if (oval[0] && oval[0]!='.')  strcpy(x,".");
		frame_extension = strcat(x,oval);
	
	//} else if (strncmp(argv[i],"-",2)==0 || strncmp(argv[i],"--",3)==0) {
	
	} else if (argv[i][0]!='-') {
		plist = alist_append(plist,MSTRDUP(argv[i]));
	} else {
		USERERROR("Wrong option given on the command line, try '-h'...");
		exit(-1);
	}
	
	}	/* (end of the whole option cycle) */
	return plist;
}












/******************	Starting and final tasks	*******************/
/**************************************************************************/


/**
 * These are the starting and final functions of the program, called from main().
 * 
 * macek_begin() is called before debug check, there is no pfield selected yet!
 * 
**/

int	macek_begin(void) {
	
	SRANDOM(time(NULL));
	errorout = stderr;
	debugout = printout = stdout;
	return 1;
}


void    ematrix_finishall(void) ;
void    frame_finishall(void) ;

int	macek_final(void) {
	
	final_check();
	
	if (errorout!=stderr && errorout!=stdout)  fclose(errorout);
	if (debugout!=stderr && debugout!=stdout)  fclose(debugout);
	if (printout!=stderr && printout!=stdout)  fclose(printout);
	
	struct_minorbases_finish();
	struct_matrixselfmaps_finish();
	frame_finishall();
	alist_finishall();
	ematrix_finishall();
	return 0;
}






























