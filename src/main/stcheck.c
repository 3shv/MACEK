
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
 * This file contains a self-check code which is called on program startup,
 * and it tries various fast tests before anything else from the program is run.
 * (The pfield of the program is already selected from the command-line parameters.)
 * Those longer tests here are run only at random, so that the program is not so slow.
 * 
**/


#include "macek.h"










/******************	Starting debug function		*******************/
/**************************************************************************/


#undef CURDLEV
#define CURDLEV	1	/* (this level is low, but printlev is also decreased by 4 here) */


int     check_pfieldar(void) ;
int     check_framerw(void) ;
int     check_inpfield(void) ;
int     check_structural(void) ;
int     check_structuralb(void) ;
int     check_extgener(void) ;

/**
 * 
 * the starting thorough debug check, and the final check........
 * 
**/


int	begin_check(void) {
	int	i,j, r=0, er=0;
	
	printlev -= 4;
	SDEBUG(CURDLEV,"%s %s (%s) in \'%s\': checking program consistency...\n",PROGNAME,PROGVER,PROGDATE,pfield_curname());
	
	for (i=j=0; i<50 && r>=0; i++) {
		r = check_pfieldar();  j += r;
	}
	if (r>=0)  SDEBUG(CURDLEV,"\tpf arithm%d",j);
	else  {PROGERROR("Internal pfield arithmetic error."); er=1;}
	perror_occured = uerror_occured = 0;
	
	r = check_framerw();
	if (r>=0)  SDEBUG(CURDLEV,"\tframe r/w%d",r);
	else  {PROGERROR("Internal frame read/write error - see %smacek-test-??.",frame_autosave_pref); er=1;}
	perror_occured = uerror_occured = 0;
	
	for (i=j=0; i<50 && r>=0 && j<3; i++) {
		r = check_inpfield();  j += r;
	}
	if (r>=0)  SDEBUG(CURDLEV,"\tinpf test%d",j);
	else  {PROGERROR("Internal pfield matrix (det,rank) testing error."); er=1;}
	
	r = check_structural();
	if (r>=0)  SDEBUG(CURDLEV,"\tstruct(reg)%d",r);
	else  {PROGERROR("Internal matroidal-structure (minors) error."); er=1;}
	perror_occured = uerror_occured = 0;
	
	r = check_structuralb();
	if (r>=0)  SDEBUG(CURDLEV,"\tstruct(iso)%d",r);
	else  {PROGERROR("Internal matroidal-structure (isom or repres) error."); er=1;}
	perror_occured = uerror_occured = 0;
	
	r = check_extgener();
	if (r>=0)  SDEBUG(CURDLEV,"\text/repr%d",r);
	else  {PROGERROR("Internal matroidal-structure error when generating extensions or representations."); er=1;}
	perror_occured = uerror_occured = 0;
	
	if (!er)  SDEBUG(CURDLEV,"\n... consistency check passed OK.\n");
	printlev += 4;
	return er? -1:1;
}


int	final_check(void) {
	
	if (perror_occured || uerror_occured) {
		printf("\n\n!!!!!!!!!!!!!!! ERROR occured during computation !!!!!!!!!!!!!!\n\n");
		if (printout!=stdout)  fprintf(printout,"\n\n!!!!!!!!!!!!!!! ERROR occured during computation !!!!!!!!!!!!!!\n\n");
		if (errorout!=stdout && errorout!=stderr)  fprintf(errorout,"\n\n!!!!!!!!!!!!!!! ERROR occured during computation !!!!!!!!!!!!!!\n\n");
	} else {
		SDEBUG(1,"\n%s %s finished OK\n",PROGNAME,PROGVER);
	}
	return 0;
}












/******************	Supplementary testing functions	*******************/
/**************************************************************************/


#undef CURDLEV
#define CURDLEV	1


/**
 * A supplementary function generating a random matrix for the tests.
**/

ematrix*	check_randommatrix(int r, int c, int w) {
	ematrix		*e;
	exp_t		xx;
	int		i,j,k;
	
	e = new_ematrix(r,c,r,c);
	for (i=0; i<r; i++)  for (j=0; j<c;j++) {
		SIGNM(e,i,j) = RANDOM()%3-1;
		if (SIGNM(e,i,j)!=0)  for (k=0; k<pfnumexp; k++)
			xx.xe[k] = RANDOM()%(2*w+1)-w;
		EXPM(e,i,j) = xx;
		pfield_tostandform(&EXPM(e,i,j),&SIGNM(e,i,j));
	}
	return e;
}


/**
 * This function checks pfield arithmetic and input of expressions in the frame.
 * Random numbers are printed as an expression into a string.
 * This string is then fed to frame_doinput(), and the resulting number is compared
 * with the number obtained directly from the starting values in a slightly different way.
 * It is necessary to ensure that all parts of computation are defined here!
 * The return value is -1 for an error, 1 if a valid computation happened, and 0 if nothing.
**/

int	check_pfieldar(void) {
	ematrix		*e1,*er;
	framestruc	*fr;
	int		j,r;
	char		buf[500];
	exp_t		xx,xx2;
	sign_t		gg,gg2;
	
	e1 = check_randommatrix(1,10,5);
	sprintf(buf," (%s)*((%s)^-3+(%s)/(%s)^2)",
			pfield_printvalue_formal(NULL,90,EXPM(e1,0,1),SIGNM(e1,0,1)),
			pfield_printvalue_to(NULL,90,EXPM(e1,0,2),SIGNM(e1,0,2)),
			pfield_printvalue_to(NULL,90,EXPM(e1,0,3),SIGNM(e1,0,3)),
			pfield_printvalue_formal(NULL,90,EXPM(e1,0,4),SIGNM(e1,0,4)));
			/* (the mix of formal and regular printing is intended here) */
	r = 0;
	xx = EXPM(e1,0,3);  gg = SIGNM(e1,0,3);
	xx2 = EXPM(e1,0,2);  gg2 = SIGNM(e1,0,2);
	j = (SIGNM(e1,0,1) && SIGNM(e1,0,4))? 0:-1;
	if (j>=0)  j = pfield_mul_ch(0,1,xx,gg,-2,EXPM(e1,0,4),SIGNM(e1,0,4),&xx,&gg);
	if (j>=0)  j = pfield_mul_ch(0,1,xx,gg,1,EXPM(e1,0,1),SIGNM(e1,0,1),&xx,&gg);
	if (j>=0)  j = pfield_mul_ch(0,-3,xx2,gg2,1,EXPM(e1,0,1),SIGNM(e1,0,1),&xx2,&gg2);
	if (j>=0)  j = pfield_sum_ch(0,1,xx,gg,1,xx2,gg2,&xx,&gg);
	if (j>=0) {
		fr = frame_doinput(buf);
		er = fr?FRMATRIX(fr):NULL;
		if (!er?1: !pfield_isequal(xx,gg,EXPM(er,0,0),SIGNM(er,0,0))) {
			r = -1;
			SDEBUG(0,"BAD:  %s  =  %s : %s\n",buf,!er?"XX":pfield_printvalue(90,EXPM(er,0,0),SIGNM(er,0,0)),pfield_printvalue(90,xx,gg));
		} else  r = 1;
		if (fr)  dispose_frame_recur(fr);
	}
	dispose_ematrix(e1);
	return (perror_occured || uerror_occured)? -1:r;
}	


/**
 * This function writes a random matrix into a file, and tests the matrix as read back.
 * It also checks option writing and inheritance.
 * These tests are performed only when frame_autosave_pref is allowed (disable with '-t-').
**/

int	check_framerw(void) {
	ematrix		*e1,*e2;
	framestruc	*fr1,*fr2, *frx=NULL, *frxx=NULL;
	int		r = 0;
	char		buf[300], **ov;
	
	if (!frame_autosave_pref)  return 0;
	e1 = check_randommatrix(10,10,20);
	fr1 = new_frame(NULL);
	frame_setmatrix(fr1,e1);
	fr2 = new_frame(fr1);
	e2 = check_randommatrix(10,10,20);
	frame_setmatrix(fr2,e2);
	frame_addoption(fr1,new_option_onestr("finherit","abc",0),1);
	frame_addoption(fr1,new_option_onestr("finherit","def ................:.............. xxxx",0),2);
	frame_addoption(fr1,new_option_onestr("finherit","finherit",0),3);
	
	snprintf(buf,250,"<%smacek-test-x",frame_autosave_pref);
	buf[250] = 0;  frame_setname(fr1,buf+1);
	r = frame_write_tree(fr1,"This file is created only for debug-testing of frame disk operations in Macek...");
	if (r>=0)  frx = frame_doinput(buf);
	if (r>=0)  r++;
	if (frx) if (FRSONS(frx))  frxx = FRSONS(frx)[0];
	if (frx && frxx) {
		if (FRNUMSONS(frx)!=1)  r = -1;
		if (!ematrix_isequal(e1,FRMATRIX(frx)))  r = -2;
		if (!ematrix_isequal(e2,FRMATRIX(frxx)))  r = -3;
		ov = frame_getoptionval_all(frxx,"finherit");
		if (alist_getlength(ov)!=3)  r = -4;
		if (ov)  dispose_alist(ov);
	} else  r = -10;
	
	strcat(buf,"x");	/* (the file name from above +x) */
	if (r>=0) {
		frame_setname(fr2,buf+1);
		r = frame_write_tree(fr2,"This file is created only for debug-testing of frame disk operations in Macek...");
	}
	if (r>=0)  frxx = frame_doinput(buf);
	else  frxx = NULL;
	if (frxx) {
		ov = frame_getoptionval_all(frxx,"finherit");
		if (alist_getlength(ov)!=1)  r = -6;
		if (ov)  dispose_alist(ov);
	}
	if (r>=0)  r++;
	dispose_frame_recur(fr1);
	if (frx)  dispose_frame_recur(frx);
	if (frxx)  dispose_frame_recur(frxx);
	return (perror_occured || uerror_occured)? -1:r;
}


/**
 * This function checks the inpfield-testing function, and determinant and rank computations.
 * The checks are cairried out on a random matrix, and only if the matrix is properly
 * represented over the pfield.
 * Determinant and rank are tested against dual and line-swaps, inpfield also against pivoting.
 * The return value is -1 for an error, 1 if a valid computation happened, and 0 if nothing.
**/

int	check_inpfield(void) {
	ematrix		*e, *e2;
	int		i,j,r,c, n1=0,sz, ret=0;
	exp_t		xx,xx2;
	sign_t		gg,gg2;
	
	sz = (pfispartial? 4:8);
	e = check_randommatrix(sz,sz,pfispartial?2:10);
	e2 = ematrix_copydual(e);
	j = ematrix_inpfield(e);
	if (j>=0) {
		ret = 1;
		ematrix_determinant(e,&xx,&gg);
		n1 = ematrix_matrank(e);
		if (n1>sz || (gg!=0 && n1<sz))  ret = -1;
	}
	if (j>=0 && ret>=0) {
		ematrix_swaprows(e2,1,3);  ematrix_swapcols(e2,0,1);
		if (ematrix_matrank(e2)!=n1)  ret = -2;
		ematrix_determinant(e2,&xx2,&gg2);
		if (!pfield_isequal(xx,gg,xx2,gg2))  ret = -3;
		for (i=0; i<5; i++) {
			r = RANDOM()%ROWSM(e2);  c = RANDOM()%COLSM(e2);
			if (SIGNM(e2,r,c)!=0)  ematrix_pivot(e2,r,c);
		}
		if (ematrix_inpfield(e2)<0)  ret = -4;
	}
	if (ret<0)  {SDEBUG(0,"\nerr%d:  det=%s %s  rank=%d %d\n",ret,
			pfield_pvalue(14,xx,gg),pfield_pvalue(14,xx2,gg2),n1,ematrix_matrank(e2));
			EMATDEBUG(0,e," !1!\t"); EMATDEBUG(0,e2," !2!\t");}
	dispose_ematrix(e);  dispose_ematrix(e2);
	return ret;
}


/**
 * This function checks some structural functions of the program.
 * It uses well-known regular matroids R10, R12, K5, and K33, and hence it works over all pfields.
 * The tests first compute the number of automorphisms and bases for R12,R10.
 * Then we verify that all 1-elem deletions of R10 are K33, and that no deletion
 * of R12 has R10-minor.
 * We also verify that R12 has branch-width 3, while R10 does not, and that R10 is self-dual.
 * The last test computes girth of R10, K5, and K5*.
**/

int	check_structural(void) {
	ematrix		**el,*e,*ee,*eee, *r10,*r12, *k5,*k33;
	int		i,j, ret = 0;
	
	el = frame_inputmatrices(" -1 1 0 0 1; 1 -1 1 0 0; 0 1 -1 1 0; 0 0 1 -1 1; 1 0 0 1 -1");
	r10 = ematrix_copy(el[0]);
	dispose_alist_mats(el);
	el = frame_inputmatrices(" 1 1 1 0 0 0; 1 1 0 1 0 0; 1 0 0 0 1 0; 0 1 0 0 0 1; 0 0 1 0 -1 -1; 0 0 0 1 -1 -1");
	r12 = ematrix_copy(el[0]);
	dispose_alist_mats(el);
	el = frame_inputmatrices(" 1 1 0 0 0 1; 1 1 0 1 1 1; 1 0 1 1 1 1; 1 0 1 0 1 0");
	k5 = ematrix_copy(el[0]);
	dispose_alist_mats(el);
	el = frame_inputmatrices(" 1 1 0 0; 1 1 1 0; 1 1 1 1; 0 1 1 1; 0 1 0 1");
	k33 = ematrix_copy(el[0]);
	dispose_alist_mats(el);
	
	if (RANDOM()%5==1) {
		if (struct_matrixselfmaps_number(r10)!=10)  ret = -2;
		if (struct_matrixselfmaps_number(r12)!=8)  ret = -3;
		if (ret>=0)  ret++;
	}
	if (RANDOM()%6==1) {
		if (ematrix_getmhash(r10)!=EM_HASHR10)  ret = -5;
		/* (definition of EM_HASHR10 is in ../include/ematrix.h for update) */
		el = ematrix_getbases_sq(r10);
		if (ret>=0)  ret++;
		if (alist_getlength(el)!=162)  ret = -6;
		dispose_alist_mats(el);
	}
	if (RANDOM()%30==1) {
		el = ematrix_getbases_sq(r12);
		if (ret>=0)  ret++;
		if (alist_getlength(el)!=441)  ret = -7;
		dispose_alist_mats(el);
	}
	if (RANDOM()%15==1) {
		e = ematrix_copy(r10);
		i = RANDOM()%10-5;  if (i>=0)  i++;
		ee = ematrix_removematid_del(e,i);
		if (ret>=0)  ret++;
		if (!struct_isequivalent(ee,k33))  ret = -10;
		dispose_ematrix(e);  dispose_ematrix(ee);
	}
	if (RANDOM()%100==1) {
		e = ematrix_copy(r12);
		i = RANDOM()%8-4;  if (i>=0)  i++;
		ee = ematrix_removematid_contr(e,i);
		if (ret>=0)  ret += 5;
		if (struct_hasminor(ee,k5))  ret = -11;
		eee = ematrix_removematid_contr(k5,i);
		if (!struct_hasminor(ee,eee))  ret = -12;
		if (struct_isconnected(eee,3))  ret = -15;
		dispose_ematrix(e);  dispose_ematrix(ee);  dispose_ematrix(eee);
	}
	if (RANDOM()%5==1)
		if (struct_hasbwidth3(r10))  ret = -20;
	if (RANDOM()%25==1) {
		if (ret>=0)  ret += 2;
		if (struct_connectivity(r10)!=4)  ret = -21;
		if (struct_connectivity(k5)!=3)  ret = -22;
	}
	if (RANDOM()%40==1) {
		if (ret>=0)  ret += 2;
		if (!struct_hasbwidth3(r12))  ret = -26;
		if (!struct_hasdelcontr_one(r12,k33))  ret = -28;
		if (struct_hasdelcontr_one(r10,k33))  ret = -29;
	}
	if (RANDOM()%10==1) {
		e = ematrix_copydual(r10);
		i = RANDOM()%5;  j = RANDOM()%5;
		if (SIGNM(e,i,j))  ematrix_pivot(e,i,j);
		if (ret>=0)  ret++;
		if (!struct_isequivalent(e,r10))  ret = -30;
		dispose_ematrix(e);
	}
	if (RANDOM()%15==1) {
		if (ret>=0)  ret += 3;
		if (struct_matgirth(r10)!=4)  ret = -35;
		if (struct_matgirth(k5)!=3)  ret = -36;
		ematrix_transpose(k5);
		if (struct_matgirth(k5)!=4)  ret = -35;
		ematrix_transpose(k5);
	}
	
	dispose_ematrix(r10);  dispose_ematrix(r12);
	dispose_ematrix(k5);  dispose_ematrix(k33);
	if (ret<0)  {SDEBUG(0,"\n\tBAD  ret=%d\n",ret);}
	return (perror_occured || uerror_occured)? -1:ret;
}


/**
 * This function checks some more structural functions of the program.
 * We test abstract isomorphism between matroids R10, and its modification
 * (over different pfields), and automorphism maps between elements.
 * Then we generate all representations of matroids.
**/

int	check_structuralb(void) {
	ematrix		**el, *r10,*r10b,*r10n;
	int		pi,pi5, ret = 0;
	
	el = frame_inputmatrices(" -1 1 0 0 1; 1 -1 1 0 0; 0 1 -1 1 0; 0 0 1 -1 1; 1 0 0 1 -1");
	r10 = r10b = ematrix_copy(el[0]);
	dispose_alist_mats(el);
	pi = pi5 = pfield_curindex();
#ifndef	BINARYONLY
	pfield_switchto("GF(5)");
	pi5 = pfield_curindex();
	el = frame_inputmatrices(" 4 1 0 0 1; 0 1 1 0 1; 1 4 4 1 4; 1 0 0 1 4; 0 0 1 4 1");
	r10b = ematrix_copy(el[0]);
	dispose_alist_mats(el);
	el = frame_inputmatrices(" 4 1 0 0 1; 0 1 1 0 1; 1 4 4 1 4; 1 2 0 1 4; 0 0 1 4 1");
	r10n = ematrix_copy(el[0]);
	dispose_alist_mats(el);
	pfield_switchto_fast(pi);
	
	if (RANDOM()%12==1) {
		if (!strmag_isisomorphic(r10,r10b,pi5))  ret = -12;
		if (strmag_isisomorphic(r10,r10n,pi5))  ret = -13;
		if (ret>=0)  ret += 2;
	}
	pfield_switchto_fast(pi5);
	if (RANDOM()%12==1) {
		if (strmag_isautmap(r10n,1,2))  ret = -16;
		if (!strmag_isautmap(r10n,1,5))  ret = -17;
		if (ret>=0)  ret += 1;
	}
	if (RANDOM()%30==1) {
		if (grepr_generate_all(r10n,pi5,NULL)!=2)  ret = -26;
		if (ret>=0)  ret += 2;
	}
	pfield_switchto_fast(pi);
	dispose_ematrix(r10);  dispose_ematrix(r10n);
#endif
	if (RANDOM()%30==1) {
		if (grepr_generate_all(r10b,pi5,NULL)!=1)  ret = -25;
		if (ret>=0)  ret += 1;
	}
	dispose_ematrix(r10b);
	if (ret<0)  {SDEBUG(0,"\n\tBAD  ret=%d\n",ret);}
	return (perror_occured || uerror_occured)? -1:ret;
}


/**
 * This function checks the extension-generating process.
 * It first tries all extensions of W_4 with forbidden K_3,3^* (only one co-graphic exists),
 * and then all co-extensions of W_4 (two graphic ones exist, one of them not co-graphic).
 * 
 * After those tests, we check how representations of F7 are generated -- those exist
 * iff the characteristic is 2.
**/

int	check_extgener(void) {
	framestruc	*fro, *fr;
	int		x,pi, ret = 0;
	
	pi = pfield_curindex();
#if	!defined(BINARYONLY) && !defined(PRIMEFIELDONLY)
	pfield_switchto("regular");
	if (RANDOM()%50==1) {
		fro = frame_doinput("!quiet;!extend c;@ext-forbid "\
			"\" 1 1 1 0 0; 1 1 1 1 1; 0 1 1 1 0; 0 0 1 1 1\";"\
			"{ \" 1 0 0 -1; -1 1 0 0; 0 -1 1 0; 0 0 -1 1\" }");
		if (fro) {
			frame_processcommands(fro);
			fr = FRSONS(fro)? FRSONS(fro)[0]:NULL;
			if (!fr?1: FRNUMSONS(fr)!=1)  ret = -51;
			dispose_frame_recur(fro);
		} else  ret = -111;
		if (ret>=0)  ret++;
	}
	if (RANDOM()%50==1) {
		fro = frame_doinput("!quiet;!extend r;{ \" 1 0 0 -1; -1 1 0 0; 0 -1 1 0; 0 0 -1 1\" }");
		frame_processcommands(fro);
		if (fro) {
			fr = FRSONS(fro)? FRSONS(fro)[0]:NULL;
			if (!fr?1: FRNUMSONS(fr)!=2)  ret = -52;
			dispose_frame_recur(fro);
		} else  ret = -112;
		if (ret>=0)  ret++;
	}
	pfield_switchto_fast(pi);
#endif
	if (RANDOM()%15==1) {
		fro = frame_doinput("!quiet;!represgen ((T)) all >((0t));{ \"@inputpf binary; 0 1 1 1; 1 0 1 1; 1 1 0 1\" }");
			/* (we assume the GF2 filed defined in every instance of the program!) */
		if (fro)  frame_processcommands(fro);
		x = (!fro?0: FRNUMSONS(fro));
		if (x!=1+(pfield_characteristic()==2))  ret = -53;
		if (fro)  dispose_frame_recur(fro);
		if (ret>=0)  ret++;
	}
	if (!pfispartial) if (RANDOM()%50==1) {
		fro = frame_doinput("!quiet;!represgen ((T)) all >((0t));!filx-isompair (s);{ \" 1 1 1; 1 1 -1; 1 -1 0\" }");
			/* (we assume the GF2 filed defined in every instance of the program!) */
		if (fro)  frame_processcommands(fro);
		x = (!fro?0: FRNUMSONS(fro));
		if (x!=1)  ret = -55;
		if (fro)  dispose_frame_recur(fro);
		if (ret>=0)  ret++;
	}
	return ret;
}




































