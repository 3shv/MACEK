#!/usr/bin/make -f
# 
#   M A C E K  ("MAtroids Computed Efficiently" Kit)
#   Copyright (C) 2001--2011  Petr Hlineny.
# 
# Go into src/ for building the programm,
# and then into exe/ for running it.
# Documentation in doc/.
# 

first:	about
about:
	@les=`which less`; \
	if [ ! -x "$$les" ]; then  les=more;  fi; \
	(cat ./ABOUT; cat doc/version | awk '/PROGVER/{print "Version  " $$3 "."}'; \
		echo; echo "Try 'make qhelp' or 'make help' ..."; echo) | $$les

help:	
	@les=`which less`; \
	if [ ! -x "$$les" ]; then  les=more;  fi; \
	(cat ./doc/README; cat doc/version | awk '/PROGVER/{print "Version  " $$3 "."}') | $$les

qhelp:	
	@les=`which less`; \
	if [ ! -x "$$les" ]; then  les=more;  fi; \
	(cat ./doc/QUICKSTART; cat doc/version | awk '/PROGVER/{print "Version  " $$3 "."}') | $$les

info:
	@inf=`which pinfo`; \
	if [ ! -x "$$inf" ]; then  inf=info;  fi; \
	cd doc; $$inf ./macek
	

comp:	compile
compile:
	@$(MAKE) -C src all
	@$(MAKE) -C src xall

clean:	
	@$(MAKE) -C src clean
	$(MAKE) -C doc/info clean

docall:	
	$(MAKE) -C doc/info all


bak:	backup
backup:
	@$(MAKE) -C pack bak


