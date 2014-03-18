
Make targets:

bak	- backup of all project files (to TOP/pack/...)

test	- compiling all for a resulting file "test"

all	- compiling all for a resulting file "macek"
		- go to TOP/exe for running the compiled program...
xall	- compiling all without debug-code for a resulting file "macek.nodebug"
		(uses *.xo and *.xa object files and archives for distinction)

edit	- editting source files (see $(EDITOR)...)

clean	- removing compilation results


Read Make.common which defines some general variables and rules
included in every makefile in the project source.

Write your local makefile adjustments to Make.local.

For documenation - see in TOP/doc and TOP/doc/info.

