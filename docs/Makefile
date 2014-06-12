SUBDIRS = source 

html:
	@for i in $(SUBDIRS); do \
        echo "make html in $$i..."; \
        (cd $$i; $(MAKE) $(MFLAGS) $(MYMAKEFLAGS) html); done

pdf:
	@for i in $(SUBDIRS); do\
	echo "make latex in $$i.."; \
	(cd $$i; $(MAKE) $(MFLAGS) $(MYMAKEFLAGS) latex); done
	echo "\def\sphinxdocclass{MO_report}" > build/latex/docs.tex
	echo "\documentclass[letterpaper,10pt,english]{MO_report}" >> build/latex/docs.tex
	tail -n +4 build/latex/Iris.tex >> build/latex/docs.tex
	sed 's/\\tableofcontents/\\tableofcontents\n\\pagenumbering\{arabic\}/' build/latex/docs.tex > build/latex/docs2.tex
	sed 's/subsection{/section{/' build/latex/docs2.tex > build/latex/documentation.tex
	(cd build/latex; pdflatex -interaction=scrollmode documentation.tex)
	# call latex again to get page numbers right...
	(cd build/latex; pdflatex -interaction=scrollmode documentation.tex);

all:
	@for i in $(SUBDIRS); do \
	echo "make all in $$i..."; \
	(cd $$i; $(MAKE) $(MFLAGS) $(MYMAKEFLAGS) all); done

install:
	@for i in $(SUBDIRS); do \
	echo "Installing in $$i..."; \
	(cd $$i; $(MAKE) $(MFLAGS) $(MYMAKEFLAGS) install); done

build:
	@for i in $(SUBDIRS); do \
	echo "Clearing in $$i..."; \
	(cd $$i; $(MAKE) $(MFLAGS) $(MYMAKEFLAGS) build); done

clean:
	@for i in $(SUBDIRS); do \
	echo "Clearing in $$i..."; \
	(cd $$i; $(MAKE) $(MFLAGS) $(MYMAKEFLAGS) clean); done

doctest:
	@for i in $(SUBDIRS); do \
        echo "Running doctest in $$i..."; \
        (cd $$i; $(MAKE) $(MFLAGS) $(MYMAKEFLAGS) doctest); done
