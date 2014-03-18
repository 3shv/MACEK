#!/usr/bin/make -f 

TOPDIR = ..

SOURCE_IGNORE = pfdefs.c
OBJS_IGNORE =

include $(TOPDIR)/Make.common

CFLAGS = $(CFLAGS_COMMON)

pfmore.o:	pfdefs.c  pfdef-more.inc pftran-more.inc


edit:
	$(EDITOR)+60+20 *.h ../include/pfield.h
	$(EDITOR)+110+40 *.c






