# Makefile thesis.tex
# Programa de Pos-Graduacao LNCC
# Modelagem Computacional
# segunda-feira, junho 4, 2007 at 14:11 - Leandro Gazoni
# Suporte aos seguintes formatos:
#
# all                Cria documento PDF
# pdfview      Cria seu arquivo PDF e VISUALIZAR PDF
# FILENAME.pdf Cria seu arquivo PDF
# clean  Limpa os arquivos desnecessarios !

# ATTENTION!
# File-extensions to delete recursive from here
EXTENSION=aux toc idx ind ilg log out lof lot lol bbl blg bak gls glo


target = thesis
WITH_INDEX = Yes
WITH_GLOSSARY = No
latex_SOURCES = $(target).tex
html = html

all: pdf

latex: $(latex_SOURCES)
		latex $(latex_SOURCES)

bib: $(latex_SOURCES)
		bibtex $(target).aux

build: $(latex_SOURCES)
		pdfview: $(target).pdf
		acroread $(target).pdf

ps: $(latex_SOURCES)
		ps2pdf thesis.ps

dvips: $(latex_SOURCES)
		dvips -t a4 thesis.dvi -o thesis.ps

pdf: $(latex_SOURCES)
		@echo "*"
		@echo "* Executando o o PDFLaTeX*"
		@echo "*"
		pdflatex -shell-escape $(target)
		@echo "*"
		@echo "* Executando o BiBiTex *"
		@echo "*"
		bibtex $(target)
		rm -f $*.log $*.aux $*.blg
		@echo "* Re-Executando o PDFLaTeX *"
		@echo "*"
		pdflatex -shell-escape $(target)
		

%.pdf: $(target).ps $(latex_SOURCES)
		@echo "*"
		@echo "* Running pdflatex *"
		@echo "*"
		pdflatex -shell-escape $(target
		@echo "*"
		@echo "* Rerunning pdflatex *"
		@echo "*"
		pdflatex -shell-escape $(target)
		@echo "*"
		@echo "* Thumbnails *"
		@echo "*"
		thumbpdf $(target).pdf
ifeq ($(WITH_INDEX),Yes)
		@echo "*"
		@echo "* Running makeindex for index *"
		@echo "*"
		makeindex -g -s include/bbbind.idx $(target)
endif
ifeq ($(WITH_GLOSSARY),Yes)
		@echo "*"
		@echo "* Running makeindex for glossary *"
		@echo "*"
		makeindex -s nomencl.ist -o $(target).gls $(target).glo
endif
		@echo "*"
		@echo "* Rerunning pdflatex *"
		@echo "*"
		pdflatex -shell-escape $(target)
		
%.ps: %.dvi
		dvips -o $(*F).ps $(*F).dvi

%.dvi: $(latex_SOURCES)
		@echo "*"
		@echo "* Running latex *"
		@echo "*"
		latex -shell-escape $(target)
		@echo "*"
		@echo "* Rerunning latex *"
		@echo "*"
		latex -shell-escape $(target)
ifeq ($(WITH_INDEX),Yes)
		@echo "*"
		@echo "* Running makeindex for index *"
		@echo "*"
		makeindex -g -s include/bbbind.idx $(target)
endif
ifeq ($(WITH_GLOSSARY),Yes)
		@echo "*"
		@echo "* Running makeindex for glossary *"
		@echo "*"
		makeindex -s nomencl.ist -o $(target).gls $(target).glo
endif
		@echo "*"
		@echo "* Rerunning latex *"
		@echo "*"
		latex -shell-escape $(target)
		
# Clean
clean:
		@echo " "
		@echo "###############################################################"
		@echo "################## LIMPANDO OS ARQUIVOS  ######################"
		@echo "###############################################################"
		@echo " "
		for EXT in ${EXTENSION}; \
		do \
		find `pwd` -name \*\.$${EXT} -exec rm -v \{\} \; ;\
		done
		rm -f ${target}.dvi
		rm -f ${target}.pdf
		rm -f ${target}.ps
		rm -f ${target}.*~
