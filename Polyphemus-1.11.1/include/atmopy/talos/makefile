CC	=	g++
CFLAGS	=	-O2 -Wall -ansi -pedantic -Dlinux
INCPATH	=	-I../../Talos
LINK	=	g++

CCTARGETS =	extract_configuration
CCOBJS	=	$(CCTARGETS:%=%.o)

all: $(CCTARGETS)

$(CCOBJS): %.o : %.cpp
	@$(CC) $(CFLAGS) -c -o $@ $< $(INCPATH)

$(CCTARGETS): % : %.o
	@$(LINK) $(LFLAGS) -o $@ $<

clean:
	@rm -f $(CCOBJS) core

cleanall: clean
	@rm -f $(CCTARGETS)
