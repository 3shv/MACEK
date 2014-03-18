#!/usr/bin/make -f

TOPDIR = ..

bak:	backup
BACKVER = $(shell cat $(TOPDIR)/src/include/version.h| awk '/PROGVER/ {print $$3}' | tr -dc '[0-9.-]')
BACKDATE = $(shell date +'"%y%m%d"')
BACKNAME = Macek_$(BACKVER)_$(BACKDATE).tar.gz
BACKFULLNAME = $(BACKNAME)
#BACKFULLNAME = pack/$(BACKNAME)
BACKCOMMAND = tar -c -X pack/Nobackup
#BACKCOMMAND = gtar -c -X pack/Nobackup -T pack/Plusbackup
backup:
	@touch $(TOPDIR)/_M_A_C_E_K_
	@rm -f $(TOPDIR)/Version-*[0-9]
	touch $(TOPDIR)/Version-$(BACKVER)
	(cd $(TOPDIR); $(BACKCOMMAND) -f - *) | gzip - >$(BACKFULLNAME)
	@echo; echo "The latest backup is in (top)/$(BACKFULLNAME)"; echo
	







