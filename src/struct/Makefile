#!/usr/bin/make -f 

TOPDIR = ..

SOURCE_IGNORE = 
OBJS_IGNORE = 

first:	test

include $(TOPDIR)/Make.common

CFLAGS = $(CFLAGS_COMMON)


edit:
	$(EDITOR)+60+20 $(wildcard *.h) ../include/struct.h
	$(EDITOR)+110+40 strmagic.c strconn.c
	$(EDITOR)+160+60 strlmap.c strminor.c
	$(EDITOR)+210+80 $(filter-out strmagic.c strconn.c strlmap.c strminor.c,$(wildcard *.c))


test:	$(OBJS_COMMON) struct.a
	$(MAKE) -C ../main



