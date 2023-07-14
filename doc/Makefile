# makefile. 


OBJS = optical_schematic.pdf layout.pdf CAD_layout.pdf
#order matters here! 
SRC = $(OBJS:.pdf=.svg)

.SUFFIXES: .svg .pdf 

.svg.pdf:
	inkscape $< --export-pdf=$@

.tex.pdf:
	pdflatex $<
	
all: $(OBJS) ao_bergamo.pdf

clean:
	rm -rf $(OBJS)
	rm -rf *.aux progress_report.log
