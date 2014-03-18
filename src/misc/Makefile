#!/usr/bin/make -f 

TOPDIR = ..

first:	test

include $(TOPDIR)/Make.common

CFLAGS = $(CFLAGS_COMMON)


edit:
	$(EDITOR)+60+20 *.h ../include/misc.h
	$(EDITOR)+110+40 *.c



test:	test.o  debug.o misc.o

