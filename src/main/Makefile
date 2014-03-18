#!/usr/bin/make -f 

TOPDIR = ..

SOURCE_IGNORE = 
OBJS_IGNORE = 

first:	mt0

include $(TOPDIR)/Make.common

CFLAGS = $(CFLAGS_COMMON)
LFLAGS = $(LFLAGS_COMMON)


edit:
	$(EDITOR)-100+180 $(wildcard *.h) ../include/debug.h ../include/misc.h
	$(EDITOR)-60+200 $(wildcard *.c)



mt0:	misc frame struct pfield ematrix gener	test

test:	$(OBJS_COMMON) ../misc/misc.a ../pfield/pfield.a ../ematrix/ematrix.a \
	 ../frame/frame.a ../struct/struct.a ../gener/gener.a \
	 ../ematrix/ematrixop.o ../pfield/pfmore.o
	$(CC) -o $@ $^

