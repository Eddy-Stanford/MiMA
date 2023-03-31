# Makefile created by mkmf $Id: mkmf,v 13.0 2006/04/10 21:20:01 fms Exp $ 

include mkmf.template.apple


.DEFAULT:
	-touch $@
all: mima.x
SRC =
OBJ =
clean: neat
	-rm -f .cppdefs $(OBJ) mima.x
neat:
	-rm -f $(TMPFILES)
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
mima.x: $(OBJ)
	$(LD) $(OBJ) -o mima.x $(LDFLAGS)
