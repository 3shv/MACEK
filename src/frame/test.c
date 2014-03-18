






#include "macek.h"
#include "fr.h"



int  main(int argc, char *argv[]) {
	framestruc      *fr=NULL;
	char		**osl, **x;
	void		**y;
	
	macek_begin();
	
	++argv, --argc;  /* skip over program name */
	pfield_switchto("dyadic");
	//pfield_switchto("nreg");
	//pfield_switchto("GF4");
	if (argc>0) {
		fr = frame_doinput_list(argv);
		//ematrix_print_nofr(FRMATRIX(fr),"\t:");
	}
	
	//osl = frame_getoptionval_ext(FRSONS(fr)[0],"op",1,0);
	//if (osl) for (x=osl; *x; x++)  printf("  \'%s\'",*x);
	//printf(";\n");
	//alist_free(osl);
	
	printf("\n====\n");
	frame_processcommands(fr);
	/*
	if (fr && FRCOMMANDS(fr))
		for (y=FRCOMMANDS(fr); *y; y++) {
			//printf("\n====\n");
			frame_applycommand(fr,*y);
		}
	*/
	if (fr)  dispose_frame_recur(fr);
	
	//if (!fopen_path("b/cc","w",frame_path_read,NULL,".exx"))  {PROGERROR("not open!");}
	
	macek_final();
	return 0;
}




/*
int	xxx(COMFDECL(f)) {
//int	xxx(int (*f)(void**[],int,void**[],int*)) {
	void	**a[2], **b[3];
	int	x,y;
	
	(*f)(a,x,b,&y);
}
*/















