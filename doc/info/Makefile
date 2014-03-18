#!/usr/bin/make -f
# 
# 

TOPDIR = ../../
EDITOR = xjoe -geometry 100x40

all:	cinfo cpdf

bak:	backup
backup:
	@$(MAKE) -C $(TOPDIR) backup

edit:
	$(EDITOR)+60+20 *.texi

test:
	pinfo ./macek.info


cinfo:	
	makeinfo macek.texi

cpdf:	
	texi2dvi macek.texi
	dvipdf macek.dvi


clean:	
	rm -f *.aux *.pdf *.cp *.cps *.fn *.fns *.ky *.log *.pg *.toc *.tp *.vr


