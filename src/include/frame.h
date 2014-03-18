
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
 * One frame keeps one matrix (of a matroid), and a list of options and commands.
 * Additionally, a frame may have a list of descendant (sub)frames -- they form
 * a tree-like structure together.
 * (This is to allow more complex operations with sets of matrices.)
 * "Options" affect how the matrix and commands are used in the program, and
 * "commands" control what operations are to be applied to the frames (matrices).
 * 
 * The input/output code for frames is written in emflex.l, emflexsup.c, and emfout.c.
 * The general frame functions are in frame.c and frameop.c.
 * The file frcoms.c contains the skeleton for executing frame commands,
 * and the (unified) command handling functions.
 * 
 * A detailed descriptions of the frame structure, and the option/command structure, follow next.
 * 
**/


#ifndef FRAME_H
#define FRAME_H






/******************	Definitions for the frame	*******************/
/**************************************************************************/



/**
 * This structure is a common definition for one frame option or frame command.
 * It holds the name of the option/command and the list of values/parameters.
 * Each value/parameter for an option/command is a text string, but it may also have
 * a numeric value which is stored separately (in addition to the string).
 * Specially, the first value may be a pfield element, which is then stored in
 * a separate part of the structure, but such options must be extra declared in
 * the lexical analyzer source(!).
 * 
 * The number of values/parameters is limited by MAXOPTPARAMS, but it is possible
 * to store more values in so called "continued" options/commands of the same name
 * right after the first structure in the frame.
 * All values and the name are stored as private copies of the strings, so that they
 * may be freed when disposing the structure.
 * 
 * Below are the macros for accessing option/command fields...
 * Always use OPTISNAME(op,nm) to compare an option/command name!
 * 
**/

#define	MAXOPTPARAMS	9	/* (the max number of option values is fixed) */

typedef struct struct_option {
	
	char	*oname;		/* this option/command name (a word, case insensitive) */
	int	optcom;		/* distinguishes between an option (0) and a command (1) */
	
	char	*params[MAXOPTPARAMS];	/* the option values/parameters as strings */
	long	nparams[MAXOPTPARAMS];	/* the option values/parameters as integers */
	int	numpar;		/* the number of scanned parameters (set in option_setparam()) */
	
	int	optix;		/* supplementary "index" of when the option was created */
	int	optcont;	/* a mark of a "continued" option (when MAXOPTPARAMS parameters is not enough) */
	
	exp_t	xparam;		/* an optional pfield-value of the first parameter */
	sign_t	gparam;
	
}	optionstruc;


#define	OPTNAME(op)		((op)->oname)
#define	OPTISNAME(op,nm)	(strcasecmp((nm),OPTNAME(op))==0)
	/* (there is also another name-matching by a prefix used in frame_addoption() check) */

#define	OPTCOM(op)		((op)->optcom)
#define	OPTISOPTION(op)		((op)->optcom==0)
#define	OPTISCOMMAND(op)	((op)->optcom==1)

#define	OPTINDEX(op)		((op)->optix)
#define	OPTCONTINUED(op)	((op)->optcont)

#ifdef	FASTPROG
#define OPTAX(op,i)		(i)
#else
#define OPTAX(op,i)		(((i)>=0&&(i)<MAXOPTPARAMS)?i:(PROGERROR("opt index out of bound!"),0))
#endif
#define OPTPARAM(op,i)		((op)->params[OPTAX(op,i)])
#define OPTNPARAM(op,i)		((op)->nparams[OPTAX(op,i)])
#define OPTNOPARAM		LONG_MAX

#define OPTPARNUM(op)		((op)->numpar)

#define OPTPAREX(op)		((op)->xparam)
#define OPTPARSG(op)		((op)->gparam)


/**
 * These functions are provided for creating and disposing option/command structures
 * (oc distinguishes option 0 or command 1, cont is for the "continued" flag).
 * The function option_setparam() sets the ix-th value/parameter of op to the string
 * (a copy of) val, and to the numeric value dv (OPTNOPARAM when number undefined).
 * The number of values in op is updated as well.
 * The function option_print_to() prints the option with its values (as if read from the input)
 * to the given text buffer of length max.
**/

optionstruc*    new_comoption_ext(char *nm, int oc, int cont) ;
#define	new_comoption(nm,oc)		new_comoption_ext(nm,oc,0)
#define	new_comoption_cont(nm,oc)	new_comoption_ext(nm,oc,1)

#define	new_option_onenum(nm,d)		option_setparam(new_comoption(nm,0),0,NULL,d)
#define	new_option_twonum(nm,d1,d2)	option_setparam(option_setparam(new_comoption(nm,0),0,NULL,d1),1,NULL,d2)
#define	new_option_onestr(nm,st,d)	option_setparam(new_comoption(nm,0),0,st,d)

void    dispose_comoption(optionstruc *op) ;

optionstruc*    option_copy(optionstruc *op) ;

optionstruc*    option_setparam(optionstruc *op, int ix, char *val, long dv) ;

optionstruc**   option_print_to(optionstruc **ol, char *buf, int max) ;



/**
 * The declaration of the frame structure follows.
 * The main information kept in a frame is the frame matrix, and lists of frame options and
 * commands.
 * (The distinction between an option and a command is that an option takes an effect
 * immediately; in order as the options are read - later overwrites former.
 * A command is, on the other hand, just stored into the list; and later it is executed
 * on request. Options may modify behavior of commands.)
 * 
 * In addition to these, a frame refers to its parent and to its list of sons,
 * thus creating a tree-like structure of frames.
 * It is user's responsibility to ensure that the frame structure is a simple forest
 * at any time - otherwise strange things would happen!
 * A name (and possibly an extension) may be given to the frame.
 * A "magic" value 7654321 is used to check for a correct access to frame structures.
 * 
 * All data refered in a frame structure must be local copies, not refered anywhere else!
 * 
**/

typedef struct struct_frame {
	
	char		*name, *ext;	/* the name (and an optional extension) of this frame */
	char		*comment;
	
	struct struct_frame	*parent;	/* the parent of this frame, or NULL */
	struct struct_frame	**sons;		/* the list of sons (NULL-terminated), or NULL */
	int		numsons;	/* the number of sons (only for informal use) */
	
	exp_t		xx;		/* the value of the last scanned expression for this frame */
	sign_t		gg;
	int		number;		/* an optional number-value of this frame */
	
	ematrix		*matrix;	/* the matrix held in the frame, or NULL */
	int		mxrows, mxcols;	/* the current matrix dimensions (only for use in scanning!) */
	int		pfindex;	/* the index of the pfield in which the matrix was read */
	int		pfindexsv;	/*  and the saved pfield index for handling option '@pfindex XX' */
	
	optionstruc	**opts;		/* the list of all options for this frame */
	optionstruc	**coms;		/* the list of all commands for this frame */
	int		lastopti;	/* the "index" of the last added option (see frame_addoption()) */
	
	int		delmark;	/* a mark that this frame should be later deleted */
	int		treemark;
	int		magic;		/* a magic value for the frame - to check correct access */
	
}	framestruc;


#define	FRMAGICSET(fr)	((fr)->magic=7654321)
#define	FRMAGICCLR(fr)	((fr)->magic=-1)
#define	FRMAGICTEST(fr)	((fr)->magic==7654321)
#define	FRMAGIC(fr)	(!FRMAGICTEST(fr)?(framestruc*)(PROGERROREXIT("accessing invalid frame")+0l):(fr))

#define	FRNAME(fr)	((fr)->name)
#define FRISNAME(fr,nm)	(strcmp((nm),FRNAME(fr))==0)
#define	FREXTENS(fr)	((fr)->ext)
#define	FRCOMMENT(fr)	((fr)->comment)

//#define	FRPARENT(fr)	(FRMAGIC(fr),(fr)->parent)
#define	FRPARENT(fr)	(FRMAGIC(fr)->parent)
#define	FRPARENT_XD(fr)	((fr)->parent)
//#define	FRSONS(fr)	(FRMAGIC(fr),(fr)->sons)
#define	FRSONS(fr)	(FRMAGIC(fr)->sons)
#define	FRSONS_XD(fr)	((fr)->sons)
#define	FRNUMSONS(fr)	((fr)->numsons)


//#define	FRMATRIX(fr)		(FRMAGIC(fr),(fr)->matrix)
#ifdef	FASTPROG
#define	FRMATRIX(fr)		((fr)->matrix)
#else
#define	FRMATRIX(fr)		(FRMAGIC(fr)->matrix)
#endif
#define	FRMATRIXROWS(fr)	((fr)->mxrows)
#define	FRMATRIXCOLS(fr)	((fr)->mxcols)
#define	FRPFINDEX(fr)		((fr)->pfindex)
#define	FRPFINDEX_SAVE(fr)	((fr)->pfindexsv)

#define	FRVALUEEXP(fr)	((fr)->xx)
#define	FRVALUESIG(fr)	((fr)->gg)
#define	FRNUMBER(fr)	((fr)->number)

#define	FROPTIONS(fr)	((fr)->opts)
#define	FRCOMMANDS(fr)	((fr)->coms)
#define	FRLASTOPTI(fr)	((fr)->lastopti)

#define	FRDELMARK(fr)	((fr)->delmark)
#define	FRTREEMARK(fr)	((fr)->treemark)









/******************	Handling frames and their options/coms	*******************/
/**********************************************************************************/



/**
 * These functions create a new frame (optionally with the given parent frame),
 * and dispose a frame (the whole subtree of the frame if rec==1).
 * The function frame_copy..() creates a copy of the given frame (not arranged in the tree
 * of subframes), possibly copying also options or commands.
 * 
 * The next function sets the given name and possibly an extension for the frame.
 * Similar function is provided for setting a new matrix in a frame (the old one is disposed).
 * (However, unlike for the names, the new matrix must be already provided as a private copy!)
 * 
 * The function frame_setson_ext() sets a new son for the parent frame, either (NULL)
 * at the end of the list, or inserted in the position of the son where.
 * (If there is a current parent, then the son is removed from it.)
 * 
**/

framestruc*     new_frame(framestruc *parent) ;

void    dispose_frame_ext(framestruc *fr, int rec, int par) ;
#define	dispose_frame_recur(fr)	dispose_frame_ext(fr,1,1)
#define	dispose_frame(fr)	dispose_frame_ext(fr,0,0)
void    dispose_frame_ls(framestruc *fr) ;

framestruc*     frame_copy_ext(framestruc *fr, char *nnm, int copt, int ccom) ;
#define	frame_copy(fr)		frame_copy_ext(fr,NULL,0,0)
#define	frame_copy_all(fr)	frame_copy_ext(fr,NULL,1,0)
#define	frame_copy_full(fr)	frame_copy_ext(fr,NULL,1,1)

void    frame_setname_ext(framestruc *fr, char *nm1, char *nm2, char *ext) ;
#define	frame_setname(fr,nm)		frame_setname_ext(fr,nm,NULL,NULL)
#define	frame_setname2(fr,nm1,nm2)	frame_setname_ext(fr,nm1,nm2,NULL)

void    frame_setmatrix(framestruc *fr, ematrix *en) ;
ematrix*        frame_extractmatrix(framestruc *fr) ;

void    frame_setson_ext(framestruc *fr, framestruc *parent, framestruc *where) ;
#define	frame_setson(fr,par)		frame_setson_ext(fr,par,NULL)
#define	frame_setson_wh(fr,par,wh)	frame_setson_ext(fr,par,wh)

framestruc**    frame_gettree_ext(framestruc *fr, int ord) ;
#define	frame_gettree_br(fr)	frame_gettree_ext(fr,0)
#define	frame_gettree_dep(fr)	frame_gettree_ext(fr,1)
#define	frame_gettree_rdep(fr)	frame_gettree_ext(fr,2)


void    frame_print_ext(FILE *fout, framestruc *fr, int verb, int mat, char *pref) ;
#define	frame_fprint(fout,fr,pref)		frame_print_ext(fout,fr,2,0,pref)
#define	frame_fprint_all(fout,fr,pref)		frame_print_ext(fout,fr,2,1,pref)
#define	frame_fprint_inf(fout,fr,pref)		frame_print_ext(fout,fr,1,0,pref)
#define	frame_fprint_brief(fout,fr,pref)	frame_print_ext(fout,fr,0,0,pref)

void    frame_printtree_ext(FILE *fout, char *buf, framestruc *fr, int verb, char *pref,
                                                 int maxbuf, int maxson, int maxdepth) ;
#define	frame_printtree(fo,fr,vb,pr)		frame_printtree_ext(fo,NULL,fr,vb,pr,0,-1,-1)
#define	frame_printtree_norep(fo,fr,vb,nr,pr)	frame_printtree_ext(fo,NULL,fr,vb,pr,0,nr,-1)
#define	frame_printtree_to(buf,max,fr,nr)	frame_printtree_ext(NULL,buf,fr,0,"",max,nr,-1)

#define FRMOUTPUT(f,pf)		frame_fprint_all(printout,f,pf)
#define FROUTPUT(f,pf)		frame_fprint(printout,f,pf)
#define FROUTPUTS(f,pf)		frame_fprint_inf(printout,f,pf)
#define	FRDEBUG(l,f,pf)		(junk = TESTDLEV(l)? (frame_fprint(debugout,f,pf),0):0)
#define	FRDEBUGS(l,f,pf)	(junk = TESTDLEV(l)? (frame_fprint_brief(debugout,f,pf),0):0)


/**
 * The first two functions add an option/command (already must be created before)
 * to the given frame, at the end of list.
 * 
 * The function frame_getoptionval_..() extracts the option values/parameters from
 * (all or the last) options of name "onm" from the frame fr and its ancestors.
 * If nv>0, then exactly nv values/parameters are taken from each option instance (substituting
 * "" for nonexistent parameters), otherwise, all parameters are taken from each instance.
 * The parameters are returned as a list of string copies.
 * (The list and the strings should be freed after use.)
 * 
 * The function frame_getoptionnum_ext() extracts the option numeric parameters from the last
 * option of name onm of the frame fr (and its ancestors only if anc==1).
 * The parameters are returned in the given array nout[] (of size nsz), by default value def.
 * 
 * The function frame_getoptionlist_ext() copies all options of the names given in
 * the string list onml, and returns a new option list.
 * 
**/

void    frame_addoption(framestruc *fr, optionstruc *op, int li) ;
void    frame_addcommand(framestruc *fr, optionstruc *op, int li) ;

char**  frame_getoptionval_ext(framestruc *fr, char *onm, int all, int nv) ;
#define	frame_getoptionval_exval(fr,on,nv)	frame_getoptionval_ext(fr,on,1,nv)
#define	frame_getoptionval_all(fr,on)		frame_getoptionval_ext(fr,on,1,0)
#define	frame_getoptionval_last(fr,on)		frame_getoptionval_ext(fr,on,0,0)

int     frame_getoptionnum_ext(framestruc *fr, char *onm, int anc, long nout[], int nsz, long def) ;
#define	frame_getoptionnum(fr,on,nout,nsz,def)		frame_getoptionnum_ext(fr,on,1,nout,nsz,def)
#define	frame_getoptionnum_noanc(fr,on,nout,nsz,def)	frame_getoptionnum_ext(fr,on,0,nout,nsz,def)

optionstruc**   frame_getoptionlist_ext(framestruc *fr, char **onml, int anc, int all) ;
#define	frame_getoptionlist_all(fr,onl)		frame_getoptionlist_ext(fr,onl,1,1)
#define	frame_getoptionlist_last(fr,onl)	frame_getoptionlist_ext(fr,onl,1,0)


/**
 * The function frame_getoptionprint..() returns a list of strings
 * that are print-outs of the selected options (like as they were read from the input).
**/

char**  frame_getoptionprint_ext(framestruc *fr, char *onm, int anc, int all) ;
#define	frame_getoptionprint(fr,on)		frame_getoptionprint_ext(fr,on,1,0)
#define	frame_getoptionprint_all(fr,on)		frame_getoptionprint_ext(fr,on,1,1)









/******************	Frame input/output	*******************/
/******************************************************************/


/**
 * These are default names and prefixes used for program input/output...
 * 
 * Moreover, there are flags file_overwrite for managing overwriting existing
 * files, and safexec for restrictions to safe execution only.
 * If file_overwrite==0, then every overwriting of an existing file is reported.
 * (In connection with safexec>10, overwriting is disallowed for other than
 * the first file_overwrite path entries!)
 * If safexec>0, then no shell execution is allowed (see !ifshell).
 * If safexec>10, then, in addition, no absolute paths are allowed
 * in file names, and no ../ (implic. changing the default path for writing files).
 * Moreover, the safexec>10 file restriction are applied also to the both paths.
 * 
**/

extern char     *frame_path_read[], *frame_path_write[],
		*frame_extension, *frame_autosave_pref;
extern int      file_overwrite, safexec;

#define	FR_PATH_LENGTH	100
		/* (see the real declarations in frame/emfout.c ...) */
#define	FR_PATH_READ	\
		".", "Matrices", "Matrices/class", "Matrices/generic", \
		"Matrices/alternative", "Matrices/graph", \
		"Procedures", "Procedures/shortcut", "Procedures/examples", \
		"mat", "out", "temp", "/tmp"
#define	FR_PATH_WRITE	\
		"*out", "*temp", ".", "*/tmp"
		/* (the '*' prefix means not to create the directory if not existing yet) */

#define	FR_DEFEXTENSION		".mck"		/* extension - should contain '.' */
#define	FR_AUTOWRITEPREF	"/tmp/.macek/"	/* where to store "auto-saved" frames - prefix */



/**
 * These are the input functions for frames and matrices...
 * The whole frame-reading is implemented by a quite complex lexical scanner
 * defined in the file emflex.l.
 * 
 * The scanner starts scanning of a given string-list (by a call to frame_doinput_ext()),
 * and processes the content plus file-includes and subframes.
 * In general, the first given string forms the resulting root frame (returned back)
 * together with the regular includes; and subframes or possible consequent strings form
 * frames that are sons of the root, or sons of the sons, etc.
 * So the result is a tree-structure of frames with the root frame returned back.
 * The following information are read for each frame:
 * Matrix entries are given (row by row) on specific matrix lines " 1 0 w ..".
 * Option lines give (arbitrary) named options and their values "@optname optval1 ..".
 * Command lines give (arbitrary) commands with their parameters "!com param1 ..".
 * If an error happens, then it returns NULL with an error message.
 * 
 * A special useful option for matrices is "@replace XX repl-text" that text-replaces
 * all occurence of the symbol XX in the matrix entries by the given replacement text.
 * (You may check correctness of the replacement value with "@require expr [01]"...)
 * 
 * It is possible to include other files in the scanner with "<filename", or to
 * define subframes "{","}", or even to include a list as subframes "{ fn1 fn2 .. }".
 * For detailed description of the scanner syntax read the definitions on the
 * top of emflex.l...
 * 
**/

framestruc*     frame_doinput_ext(char *ins, char **insl, framestruc *infr, char *name, int rec) ;
#define	frame_doinput(in)		frame_doinput_ext(in,NULL,NULL,NULL,0)
#define	frame_doinput_list(inl)		frame_doinput_ext(NULL,inl,NULL,NULL,0)
#define	frame_doinput_ap(in,ifr)	frame_doinput_ext(in,NULL,ifr,NULL,0)

			/* default name formats used in yy_endframe(): */
#define	FR_DEFSUBNAME		"-%d"		/* values: son number */
#define	FR_DEFSUBNONAME		"noname%d-%d"	/* values: subfr level, son number */
#define	FR_DEFNONAME		"noname"	/* values:  */


/**
 * The function frame_inputmatrices() uses the frame-reading function, but then
 * it extracts all matrices from the frame tree to a separate list which is returned.
 * 
**/

ematrix**       frame_inputmatrix_ext(char *ins, int one) ;
#define	frame_inputmatrices(ins)	frame_inputmatrix_ext(ins,0)


/**
 * In addition to "@replace .." for matrices, one may give options like "@sub-XX xxx"
 * for macro-substitution in the other than matrix lines (applies to other option values,
 * command parameters and include names).
 * The meaning is that each subsequent occurence of $XX is replaced with xxxx.
 * (Do not mistake these substitutions with $1,$2 parameter descriptions as below...)
 * If "@sub-XX xxx" is not found in the options, then "@subd-XX xxx" is tried as the default.
 * 
 * It is even possible to call command-procedures like "&procedure p1 p2 ..".
 * Such a call simply defines "@sub-param1 p1", etc, and then includes the file "procedure".
 * This is intended for easy calling of user-defined complex procedures from external files.
 * 
**/

#define	FR_SUBSTPREFIX		"sub-"
#define	FR_SUBSTDEFPREFIX	"subd-"
/* the use of '$' or '${}' for macro substitution is implicitly defined in emflexb.l ... */

#define	FR_SUBPARAMNAME		FR_SUBSTPREFIX"param0"
#define	FR_SUBPARAMRES		FR_SUBSTPREFIX"paramres"

/**
 * The following two options are used for erasing other options that appear in the frame.
 * "erase" clears the last previous option occurence, while "eraseall" erases all previous.
**/

#define	FR_ERASEOPT		"erase"
#define	FR_ERASEALLOPT		"eraseall"


/**
 * This function writes the given frame fr to the file whose name is obtained from this
 * frame name (or from ancestors???, or noname).
 * The options are written only if opt==1 (then @inherit and @inheritall apply).
 * If onm==1, then the @name option is written first, before the matrix (onm>1 even for the root).
 * If rec==1, then the descendants of the frame fr are written to the same file as subframes.
 * The return value is 0 fo OK, and negative values for an error.
 * 
 * Moreover, the writing function uses several constant strings presented below
 * when writing frames.
 * Remember that these strings must be compatible with the input scanner in emflex.l.
 * 
**/

int     frame_write_ext(framestruc *fr, int rec, int opt, int onm, char *dir, char *comment) ;
#define	frame_write_nopt(fr)		frame_write_ext(fr,0,0,0,NULL)
#define	frame_write_one(fr,com)		frame_write_ext(fr,0,1,0,NULL,com)
#define	frame_write_tree(fr,com)	frame_write_ext(fr,1,1,1,NULL,com)
#define	frame_writed_one(fr,dir,com)	frame_write_ext(fr,0,1,0,dir,com)
#define	frame_writed_tree(fr,dir,com)	frame_write_ext(fr,1,1,1,dir,com)

		/* symbols used when writing the frame (compare with emflex.l declarations) */
#define	FR_SUBFRAME		"SUBFRAME"
#define	FR_ENDFRAME		"EOFRAME"
#define	FR_OPTION		"@"
#define	FR_COMMAND		"!"
#define	FR_COMMENT		"# "

int     frame_write_matcom(framestruc *fr, int opt, int onm, char *dir, char *comment) ;
#define	frame_writed_onemc(fr,dir,com)	frame_write_matcom(fr,1,0,dir,com)


/**
 * Printing simple help on options and commands...
 * 
**/

void    frame_foptionhelp(int ch, FILE *fo) ;
void    frame_fcommandhelp(int ch, FILE *fo) ;
void    frame_framehelp(int ch, FILE *fo) ;













/******************	Handling frame commands		*******************/
/**************************************************************************/



/**
 * Commands are collected when reading the frame, and they are later sequentially
 * executed by a call to frame_processcommands().
 * The file frameop.c contains a general interface for calling command handling functions.
 * This interface prepares the command parameters according to the parameter values in the
 * command, and to the requested parameter types for the command in frcoms.c.
 * Then it also stores the command results back to the frame subtree (or prints out).
 * 
 * The command handling functions are defined in frcoms.c, associated with the command name.
 * It is possible to modify command names with special prefixes (like FRES_FILTPREFIX) that
 * change the return type (like from printing to a yes/no filter).
 * The major parameter-description defines follow here.
 * 
 * The first define here is a standard declaration of a command-handling function.
 * - instr  gets the list of string-parameters,
 * - inls[], inlsn  get the input lists (of matrices or frames),
 * - inum[], inumn  get the list of input integer numbers,
 * - rettype  is the required return type,
 * - and outls[], outnum[]  are the output list (arbitrary) and output integer array.
**/


	/* to be used as a unified declaration for a command-handling function (return int): */
#define	COMFUNCTION(ff)	ff(char **instr __attribute__ ((unused)), void **inls[] __attribute__ ((unused)),\
				int inlsn __attribute__ ((unused)), long inum[] __attribute__ ((unused)),\
				int inumn __attribute__ ((unused)), int rettype __attribute__ ((unused)),\
				void **outls[] __attribute__ ((unused)), long outnum[] __attribute__ ((unused)))

/**
 * The parameter description addresses nodes in the frame tree.
 * We use the well-known correspondence between balanced bracketings like (()(())())
 * and rooted trees to move through the tree.
 * The first level () points to the tree root.
 * Moreover, to make the expressions shorter, we allow repetitions like (5(..)) which
 * means to repeat the (balanced) part after the number that many times.
 * If we use a non-positive number n, we mean all sons up to the last -n.
 * One may also write (/name/...) to skip preceeding sons of a frame up to (not including)
 * the son named "name".
 * The closing character | closes all open braces in the description.
 * 
 * We pick a one frame by a letter 't' or 'T' from the tree, and all sons of a frame
 * by the letter 's' or 'S' (capital letters mean that the parameters should not be deleted
 * after executing the command).
 * For example, ((s)(-2t)) picks all sons of the first son of the root, plus all remaining
 * sons of the root but the last two.
 * 
 * All picked frames are stored in one input list per each parameter.
 * One may also use '+' to add more selection to the same list (even the same frame again).
 * It is possible to pick one of the previous resulting lists with special variables $1-$9.
 * (This is a different mechanism than the text substitutions applied in the lexical scanner.)
 * 
 * The character '>' starts the description of the return parameter, see below...
**/

#define	FRGET_DOWN	'('
#define	FRGET_UP	')'
#define	FRGET_CLOSE	'|'
#define	FRGET_BYNAME	'/'

#define	FRGET_THIS	't'
#define	FRGET_THISKEEP	'T'
#define	FRGET_SONS	's'
#define	FRGET_SONSKEEP	'S'

#define	FRGET_CONCAT	'+'
#define	FRGET_RESULT	'>'

#define	FRGET_PREVRESOLD	'$'
#define	FRGET_PREVRES		'^'
#define	FRGET_PREVRESKEEP	'~'

		/* all characters that are picking parameter frames from the tree: */
#define	FRGET_ISGET(ch)	((ch)==FRGET_THIS || (ch)==FRGET_THISKEEP || (ch)==FRGET_SONS || (ch)==FRGET_SONSKEEP)


	/* the major types of command parameters (see frame_applycommand() and frcoms.c): */
enum FCPARAMTYPES	{
	
	FCPAR_STRING=100,	/* string and integer parameters are given separately to the handle */
	FCPAR_INTEGER,
				/* one matrix or a list of them (list may be empty if FCPAR_LISTMX0) */
	FCPAR_ONEMX, FCPAR_LISTMX, FCPAR_LISTMX0,
				/* similarly one or list of frames (possibly requiring matrices)... */
	FCPAR_ONEFR, FCPAR_LISTFR, FCPAR_LISTFR0,
	FCPAR_ONEFRM, FCPAR_LISTFRM, FCPAR_LISTFRM0,
	FCPAR_LISTFRMP0,	/* frame list with matrices, but not testing in the same pfield */
	
	FCPAR_
};

		/* when getting matrix / frame input list for a parameter: */
#define	FCPAR_MXGROUP(pr)	((pr)==FCPAR_ONEMX || (pr)==FCPAR_LISTMX || (pr)==FCPAR_LISTMX0)
#define	FCPAR_FRGROUP(pr)	((pr)==FCPAR_ONEFR || (pr)==FCPAR_LISTFR || (pr)==FCPAR_LISTFR0 ||\
				 (pr)==FCPAR_ONEFRM || (pr)==FCPAR_LISTFRM ||\
				  (pr)==FCPAR_LISTFRM0 || (pr)==FCPAR_LISTFRMP0)

		/* when only one matrix/frame is requested, and when an empty list is permitted: */
#define	FCPAR_ONEREQ(pr)	(((pr)==FCPAR_ONEMX||(pr)==FCPAR_ONEFR||(pr)==FCPAR_ONEFRM)? 1:\
				 (((pr)==FCPAR_LISTMX0||(pr)==FCPAR_LISTFR0||\
				  (pr)==FCPAR_LISTFRM0||(pr)==FCPAR_LISTFRMP0)? 0:-1))

		/* when a matrix in the frame is required (default for MXGROUP...): */
#define	FCPAR_MATREQ(pr)	(((pr)==FCPAR_ONEFRM || (pr)==FCPAR_LISTFRM || (pr)==FCPAR_LISTFRM0)? \
				 1: ((pr)==FCPAR_LISTFRMP0? -1: 0))



/**
 * A command may give results in several ways:
 * Printed commands just print (or write to file) their output.
 * Filter-commands tell which of the input frames should be preserved and which deleted.
 * Filter-rem-commands tell which of the input frames should be remembered as the last result.
 * Matrix-replace commands work over one matrix, and they replace it with a new matrix.
 * List-commands read their input (not modify), and return a list of output matrices/frames.
 * 
 * The return value position (applies for list-commands) is determined similarly as
 * the input parameter position in the tree by a "bracketed expression" - see above.
 * However, one important difference is that for output all missing tree nodes are
 * automatically created, possibly with the name given by (/name/...) - only new ones.
 * Also, non-negative repetition in the expression may refer to the number of output
 * frames if it is followed by 't' or 's'.
 * The character 's' for storage means to create a new parent node getting all
 * the results as its sons (rather than using an existing parent!).
 * 
 * Try to play with the command '!move .. ..' to understand parameter descriptions better.
**/

	/* the major types of command results (see frame_applycommand() and frcoms.c): */
enum FCRETURNTYPES	{
	
	FCRET_NOTHING=100,
	FCRET_PRINTED,		/* the result is already printed, nothing to store */
	FCRET_YESNOPRINT,	/* outmx has a list of yes/no values to be printed out */
	FCRET_PRINTED_F,	/* - the same, but may be converted to any yes/no filter */
	
	FCRET_YESNOFR,		/* outmx has a list of yes/no values matching input lists of frames */
	FCRET_NOYESFR,		/* outmx has a list of no/yes values - inverted meaning from the previous */
	FCRET_LREMFR,		/* outmx has a list of yes/no values to remember as the last result, but nothing to delete */
	FCRET_LREXFR,		/* outmx has a list of no/yes values to rememebr as the last result, but nothing to delete */
	
	FCRET_REPLMX1,		/* one resulting matrix stored to the place of the first parameter */
	FCRET_REPLMXL,		/* list of resulting matrices stored to the places of the first list */
	FCRET_NUMBER,		/* returns an integer number (in outnum[0]) */
	
	FCRET_ONE,		/* these are real one- or list-returning functions: */
	FCRET_LIST,		/*  ... */
	FCRET_PARTLIST,		/* part of the return list (only those having out-storage) */
	
	FCRET_
};

		/* when marking also parent of a son-list: (see frame_getparamfr_list()) */
#define	FCRET_SMARKREQ(pr)	((pr)==FCRET_ONE || (pr)==FCRET_LIST || (pr)==FCRET_PARTLIST)

		/* when there is (only?) printed output of the command: */
#define	FCRET_PRINTREQ(pr)	((pr)==FCRET_PRINTED || (pr)==FCRET_PRINTED_F || (pr)==FCRET_YESNOPRINT)
#define	FCRET_PRINTOWN(pr)	((pr)==FCRET_PRINTED || (pr)==FCRET_PRINTED_F)

		/* when the return is of yes/no type output (special) list: */
#define	FCRET_ISYESNO(pr)	((pr)==FCRET_YESNOFR || (pr)==FCRET_NOYESFR || (pr)==FCRET_LREMFR || (pr)==FCRET_LREXFR)
#define	FCRET_ISYESNOPR(pr)	(FCRET_ISYESNO(pr) || (pr)==FCRET_YESNOPRINT)
#define	FCRET_ISYESNOY(pr)	((pr)==FCRET_YESNOFR || (pr)==FCRET_LREMFR)

		/* when there is a real list of return matrices/frames: */
#define	FCRET_LISTOUT(pr)	((pr)==FCRET_ONE || (pr)==FCRET_LIST || (pr)==FCRET_PARTLIST)

		/* when only one frame is stored, and when a partial storing is permitted */
#define	FCRET_ONEREQ(pr)	(((pr)==FCRET_ONE)? 1: (((pr)==FCRET_PARTLIST)? 0:-1))


/**
 * Each handling function should return a value (flags in OR) indicating what should be
 * done with the output (and also to distinguish between frames and matrices).
 * (If the value is not returned in OR, then the results are not stored, and the input
 * parameters of filters are not deleted...)
 * The special flag FRES_STOREREM tells to remember the return list for a possible later
 * use as an input parameter ($1..).
 * 
 * Moreover, the lower bits of the command return value are further used in parental
 * functions (returned via FRES_PROCSTRIP(r)), like for command-flow control, see below.
 * 
**/

#define	FRES_DELETEMARK		0x0100000

#define	FRES_STOREMX		0x0200000
#define	FRES_STOREMXNM		0x2000000
#define	FRES_STOREFR		0x0400000
#define	FRES_STORENUM		0x0800000
#define	FRES_STOREREM		0x1000000

#define	FRES_PROCSTRIP(r)	((r)<0?(r):((r)&0xfffff))

		/* when the command processing should be restarted, or commands skipped */
#define	FRES_RESTART		0x0010000
#define	FRES_SKIPCOM		0x0020000
extern int		frame_skipcommands;

#define	FRES_PROCSTRIP2(r)	(((r)&0xff))


		/* command-name prefixes modifying FCRET_PRINTED_F to FCRET_YESNOFR or FCRET_NOYESFR */
#define	FRES_FILTPREFIX		"filt-"
#define	FRES_INVFILTPREFIX	"filx-"
		/* command-name prefixes modifying FCRET_PRINTED_F to FCRET_YESNOFR or FCRET_LREMFR */
#define	FRES_REMPREFIX		"rem-"
#define	FRES_INVREMPREFIX	"rex-"
#define	FRES_INVREMPREFIX2	"remx-"


/**
 * These functions apply the given command, and execute all commands in the given tree...
 * After executing each command inside frame_processcommands(), the command is immediately
 * deleted from the list at the frame.
 * One may append more commands during execution...
 * The function frame_processcommands() is re-entrant (up to limited depth of recurrence),
 * but be careful with frame_applycommand() alone.
**/

int     frame_applycommand(framestruc *fr, optionstruc *com) ;

int     frame_processcommands(framestruc *fr) ;













#endif	/* (of #ifndef FRAME_H) */




































