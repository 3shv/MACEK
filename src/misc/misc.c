
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
 * various general variables and functions for the Macek package.....
 * 
 * 
**/



#include "macek.h"









/******************	Miscelanelous	*******************/
/**********************************************************/


/**
 * These are various general variables used in the project.
 *  perror_occured,uerror_occured detect errors during program run,
 *  crash_all_errors causes the program to crash with a core image on an error.
 *  errorout,debugout,printout define output streams,
 *  printlev manages output verbosity of debug messages.
 * 
**/


int	perror_occured = 0;
int	uerror_occured = 0;
int	crash_all_errors = 0;

int	junk = 0;

int	timeprecout = 3,
	timeprecdiv = 1;
int	printlev = DEBUGLEV-3;

void	*mmmm;	/* used in macros... */
char	*ssss,  debugoutbuf[10000];


FILE	*errorout,
	*debugout,
	*printout;

char	*printoutpref = "    ~\t";


#ifdef snprintf		/* replacement of snprintf() on systems without it */
char	snprintf_buf[20000];
#endif












/******************	File utilities	*******************/
/**********************************************************/



/**
 * This function tries to open a file named  fn  in mode mod ("w" or not),
 * using given path (and optional extension ext and subdirectory dir).
 * Leading and trailing whitespaces are removed from names.
 * 
 * The file name fn may normally contain a directory prefix, which is
 * used with it -- fn starting with "/", "./", or "../" gives an absolute
 * or fixed relative directory, while others are prefixed with the search path.
 * When ext is given, that is appended as a file extension to fn -- for reading
 * we try first no extension and then with extension, for write extension first.
 * (Only when no extension exists in fn!)
 * The parameter dir if not NULL is tried inserted between the path and fn.
 * 
 * For writing files (mod starts with "w"), the missing subdirectories from fn
 * are always created, and the subdirectories tried from dir and the path are
 * created only when the path entry does not start with special '*'.
 * If file_overwrite==0, then overwriting an existing file causes a warning
 * message; but this does not apply to "frame_autosave_pref/fn".
 * 
 * The return value is the new file pointer.
 * Additionally, the entry of the given path which was used is stored in
 * fopen_pathdir_used, and the actual file name (incl extension) is copied
 * in fopen_fname_used -- do not change these strings!!
 * 
**/

char	*fopen_pathdir_used = NULL;	/* records the directory found in the path during an open call */
char	*fopen_fname_used = NULL;	/* records the complete file name found in an open call */

FILE*	fopen_path(char *fn, char *mod, char **path, char *dir, char *ext) {
	char	*dfn, *buf, *xx,*xbuf;
	FILE	*fp = NULL;
	int	i,j,l,wri, bl1,bl2;
	
	if (fn==NULL || path==NULL) {
		PROGERROR("No filename or path given for reading!");  return NULL;
	}
				/* removing leading and trailing whitespaces, preparing fnames */
	while (*fn==' ' || *fn=='\t' || *fn=='\n')  fn++;
	dfn = MSTRDUP(fn);
	for (i=strlen(dfn)-1; i>=0 && (dfn[i]==' ' || dfn[i]=='\t' || dfn[i]=='\n'); i--) ;
	dfn[i+1] = 0;
	bl1 = strlen(dfn)+(ext?strlen(ext):0)+(dir?strlen(dir):0)+32;
	bl2 = (path[0]?strlen(path[0]):0)+(path[1]?strlen(path[1]):0)+320;
	buf = MMALLOC((bl1+bl2)*sizeof(buf[0]));
	if (fopen_fname_used)  FREE(fopen_fname_used);
	wri = (mod[0]=='w');	/* indicates the write(-create) mode for file openning */
	
				/* given ext is tried when no extension exists in fn: */
	for (xx=dfn+strlen(dfn)-1; xx>dfn && xx[0]!='/' && xx[0]!='.'; xx--) ;
	if (wri && xx[0]=='.')  ext = NULL;
	
	i = -1;	 xbuf = NULL;	/* looking for the file to open, first a complete fname: */
	if (safexec<10 && (dfn[0]=='/' ||\
			(dfn[0]=='.'&&dfn[1]=='/') || (dfn[0]=='.'&&dfn[1]=='.'&&dfn[2]=='/')) )
		for (j=0; j<2 && !fp; j++) {
			strcpy(buf,dfn);		/* (openning absolute file name) */
			if (j==1 && !ext)  continue;
			if (j==(!wri) && ext) { strcat(buf,ext); }
			fp = fopen_dir(xbuf=buf,mod,wri,1,0);
	
	} else			/* else searching directory path, trying also ext and dir: */
		for (i=0; !fp && path[i]; i++) for (j=0; j<4 && !fp; j++) {
			strncpy(buf,path[i],bl2-10);  buf[bl2-9] = 0;
			xbuf = (buf[0]=='*'? buf+1:buf);/* (path starting with * is not created for writing!) */
			l = (buf[0]=='*'? strlen(buf):0);
			if (strlen(xbuf)>0)  strcat(buf,"/");
			if (j%2==1 && !dir)  continue;	/* (inserting optional dir) */
			if (j%2==1) { strcat(buf,dir);  strcat(buf,"/"); }
			strcat(buf,dfn);
			if (j/2==1 && !ext)  continue;	/*  (ext is preferred when writing and no ext exists in fn) */
			if (j/2==(!wri) && ext) { strcat(buf,ext); }
			fp = fopen_dir(xbuf,mod,wri,l,i+1);
	}
	fopen_pathdir_used = (i<0?"": path[i-1]);
	fopen_fname_used = MSTRDUP((xbuf&&fp)? xbuf:"");
	FREE(dfn);  FREE(buf);
	return (fp==(FILE*)-1)?NULL: fp;
}

/**
 * This function opens a full file named  fn  in mode mod (fn includes the directory!).
 * Parameter wr says whether it is a writing mode.
 * All directories from fn starting at position pre are created automatically.
 * The optional index pix determines the index on a file path.
 * 
**/

FILE*   fopen_dir(char *fn, char *mod, int wr, int pre, int pix) {
	FILE    	*fp;
	char		*path,*pp;
	struct stat	stf;
	extern int	file_overwrite;
	
	if (safexec>10) {	/* safe - not allowing absolute or ../ file names, modifies dfn: */
		if (fopen_safename(fn,wr)>0 && pix<=1)
			{USERERROR("Absolute and ../ file names are not allowed in the safest mode, changed to \"%s\".",fn);}
	}
	if (!wr) {	/* normal open for reading, only regular files, nothing is created */
		if (stat(fn,&stf)!=0)  return NULL;
		if (!S_ISREG(stf.st_mode))  return NULL;
		fp = fopen(fn,mod);
		return fp;
	}
			/* different open routine for writing - check the file and create dirs */
	if (pix>file_overwrite? !stat(fn,&stf):0) {
				/* check if the file exists, for path position >file_overwrite */
		if (safexec>10) {
			USERERROR("Not allowed to overwrite \"%s\" in the safest mode!\n",fn);
			return (FILE*)-1;
		}
		if (frame_autosave_pref?strncmp(fn,frame_autosave_pref,strlen(frame_autosave_pref))!=0:1)
			{USERERROR("Warning - going to overwrite an existing file \"%s\".",fn);}
	}
	path = MSTRDUP(fn);	/* special for writing - creating subdirectories, after pre */
	for (pp = strchr(path+pre,'/'); pp; pp = strchr(pp+1,'/')) {
		*pp = 0;	/* (creates all required subdirectories after pre) */
		if (pp!=path && stat(path,&stf)!=0)
			mkdir(path,0755);
		*pp = '/';
	}
	FREE(path);
	fp = fopen(fn,mod);	/* the real file-write openning (finally) happens here */
	return fp; 
}


/**
 * This function makes a safe file name fn (directly changing the string fn).
 * That means the leading / and all | \ and ../ are replaced with _ .
 * 
**/

int	fopen_safename(char *fn, int wr) {
	int	i,l, r=0;
	
	if (fn[0]=='/') { fn[0] = '_';  r = 1; }
	l = strlen(fn);
	if (fn[l-1]=='/' || (fn[l-1]=='.' && (l<2?1: fn[l-2]=='/'))) {
		fn[l-1] = '_';  r = 1111; }
	for (i=0; i<l; i++)
		if (fn[i]=='\\' || fn[i]=='|')  fn[i] = '_';
	for (i=1; i<l; i++)
		if (fn[i-1]=='.' && fn[i]=='.' && fn[i+1]=='/') {
			fn[i] = fn[i-1] = fn[i+1] = '_';  
			r = i+1; 
		}
	return r;
}



/**
 * 
 * extracting the filename (not the dir path and not extension)
 * - return a string copy...
 * 
 * analogous for directory path...
 * 
**/

char*	fname_extract(char *fn) {
	char	*x,*y,*z;
	
	for (x=y=z=fn; *x; x++) {
		if (*x=='/')  y = x+1;
		if (*x=='.')  z = x;
	}
	x = MSTRDUP(y);
	if (z>y)  x[(z-y)] = 0;
	return x;
}

char*	fdir_extract(char *fn) {
	char	*x,*y;
	
	x = MSTRDUP(fn);
	for (y=x; *x; x++) {
		if (*x=='/')  y = x;
	}
	y[0] = 0;
	return x;
}











/******************	Pointer arrays	*******************/
/**********************************************************/


/**
 * These functions implement allocation of "pointer arrays" in 2D and 3D.
 * These are accessed as pt[][] and pt[][][], and are freed simply with FREE(pt),
 * 
**/

void*	malloc_twodim(int sz, int d1, int d2) {
	void	**pt,*ptd;
	int	i;
	
	pt = MMALLOC((d1+d1%2)*sizeof(void*)+d1*d2*sz);
	ptd = pt+(d1+d1%2);
	for (i=0; i<d1; i++)  pt[i] = (byte1*)ptd+i*sz*d2;
	return pt;
}

void*	malloc_threedim(int sz, int d1, int d2, int d3) {
	void	**pt,**ptm,*ptd;
	int	i,j;
	
	d1 += d1%2;
	pt = MMALLOC((d1+d1*d2)*sizeof(void*)+d1*d2*d3*sz);
	ptm = pt+d1;  ptd = ptm+d1*d2;
	for (i=0; i<d1; i++)  pt[i] = ptm+i*d2;
	for (i=0; i<d1; i++)  for (j=0; j<d2; j++)
		((void***)pt)[i][j] = (byte1*)ptd+(i*d2+j)*sz*d3;
	return pt;
}




































