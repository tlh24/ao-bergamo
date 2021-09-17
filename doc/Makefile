# makefile. 


OBJS = optical_schematic.pdf layout.pdf
#order matters here! 
SRC = $(OBJS:.pdf=.svg)

.SUFFIXES: .svg .pdf 

.svg.pdf:
	inkscape $< --export-pdf=$@
	
all: $(OBJS)

clean:
	rm -rf $(OBJS)
