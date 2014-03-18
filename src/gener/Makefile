#!/usr/bin/make -f 

TOPDIR = ..

SOURCE_IGNORE = 
OBJS_IGNORE = 

first:	test

include $(TOPDIR)/Make.common

CFLAGS = $(CFLAGS_COMMON)

$(OBJS_COMMON):		gener-more.inc  ../struct/str.h


edit:
	$(EDITOR)+60+20 $(wildcard *.h) ../include/gener.h
	$(EDITOR)+110+40 elimseq.c gmatext.c
	$(EDITOR)+160+60 $(filter-out elimseq.c gmatext.c,$(wildcard *.c))


test:	$(OBJS_COMMON) gener.a
	$(MAKE) -C ../main




