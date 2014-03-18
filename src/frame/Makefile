#!/usr/bin/make -f 

TOPDIR = ..

SOURCE_IGNORE = emflexsup.c frcomst.c
OBJS_IGNORE = 

first:	test

include $(TOPDIR)/Make.common

CFLAGS = $(CFLAGS_COMMON)
LFLAGS = $(LFLAGS_COMMON)


edit:
	$(EDITOR)+60+20 $(wildcard *.h) ../include/frame.h
	$(EDITOR)+110+40 $(wildcard frame*.c)
	$(EDITOR)+160+60 frcoms.c frcomst.c fropts.c
	$(EDITOR)+210+80 $(wildcard *.l) emflexsup.c emfout.c

# additional neccessary dependencies (emflex.c should be kept even as intermediate)
emflex.c:	emflexsup.c fr.h
emflexb.c:	fr.h
frcoms.o:	frcomst.c frcoms-more.inc


../misc/misc.a:		misc
../pfield/pfield.a:	pfield
../ematrix/ematrix.a:	ematrix

test:	test.o $(OBJS_COMMON)  ../pfield/pfield.a ../misc/misc.a ../ematrix/ematrix.a  
	$(CC) -o $@ $^



