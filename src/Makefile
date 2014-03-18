#!/usr/bin/make -f 
# 
# This is the main compilation makefile for the Macek project.
# (Local modifications should be entered in Make.local.)
# 
# Use 'make all' and 'make xall' for the complete compilations.


TOPDIR = .

include $(TOPDIR)/Make.common

CFLAGS = $(CFLAGS_COMMON)

AOUTS = $(shell for i in $(SRC_SUBDIRS); do echo $$i/$$i.a; done)
XAOUTS = $(shell for i in $(SRC_SUBDIRS); do echo $$i/$$i.xa; done)
# notice that $(AOUTS) must be used twice in linking, otherwise future references
# are not found, the other solution with '--whole-archive' is not portable...
#	$(CC) -o $(EXEOUT) -Xlinker --whole-archive $(AOUTS) -Xlinker --no-whole-archive


all:	allobj	exedir
	$(CC) -o $(EXEOUT) $(AOUTS) $(AOUTS)
	ln -s ../src/$(EXEOUT) ../exe/$(EXEOUT) 2>/dev/null; true

test:	allobj	
	$(CC) -o test */*.o

allobj:
	@for i in $(SRC_SUBDIRS); do \
		$(MAKE) -C $$i obj; \
	done

xall:	allxobj  exedir
	$(CC) -o $(EXEOUTND) $(XAOUTS) $(XAOUTS)
	ln -s ../src/$(EXEOUTND) ../exe/$(EXEOUTND) 2>/dev/null; true

allxobj:
	@for i in $(SRC_SUBDIRS); do \
		$(MAKE) -C $$i xobj; \
	done

exedir:
	@mkdir -p ../exe/temp
	@#mkdir -p ../exe/out

edit:
	$(EDITOR)+60+20 $(wildcard Make*[a-z])


flags:	
	@echo CFLAGS_COMMON: $(CFLAGS_COMMON)
	@echo
	@echo LFLAGS_COMMON: $(LFLAGS_COMMON)
	@echo A: $(AOUTS)

xxx:	
	
