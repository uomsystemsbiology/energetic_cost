#Makefile for MTT documentation

###############################################################
## Version control history
###############################################################
## $Id: Makefile,v 1.3 2002/04/27 11:07:17 geraint Exp $
## $Log: Makefile,v $
## Revision 1.3  2002/04/27 11:07:17  geraint
## Fixed [ 549265 ] make install-doc fails (cp mtt_toc.html)
## Replaced "mtt.html mtt_toc.html" with "mtt*html".
## This works regardless of local texi2html options.
##
## Revision 1.2  2002/04/26 17:17:20  geraint
## Allow variables to be set outside Makefile.
##
## Revision 1.1  2001/06/04 08:18:52  gawthrop
## Putting documentation under CVS
##
## Revision 1.2  1998/07/02 18:40:07  peterg
## More on install option
##
## Revision 1.1  1998/05/14 09:19:52  peterg
## Initial revision
##
###############################################################


#MTTPATH = /home/peterg/mtt_new/mtt
#DOCPATH = /home/peterg/web-docs/software/MTT/doc
#INFOPATH = /usr/info

MTTPATH ?= /home/peterg/mtt_new/mtt
MTT_COMPONENTS ?= /home/eng4/peterg/mtt_new/mtt/lib/comp
MTT_EXAMPLES ?= /home/eng4/peterg/mtt_new/mtt/lib/examples
DOCPATH ?= /home/peterg/web-docs/software/MTT/doc
#INFOPATH = $(DOCPATH)
INFOPATH ?= /usr/info

all: mtt.info  mtt.html mtt.pdf mtt.ps.gz 

mtt.info: mtt.texi
	echo "Making info manual. Please wait ..."; makeinfo mtt.texi

mtt.dvi: mtt.texi
	echo "Making dvi manual. Please wait ..."; tex mtt.texi; tex mtt.texi

mtt.html: mtt.texi
	echo "Making html manual. Please wait ..."; texi2html -glossary  mtt.texi

mtt.ps.gz: mtt.dvi
	echo "Making ps manual. Please wait ..."; dvips -o mtt.ps mtt; gzip -f mtt.ps

mtt.pdf: mtt.texi
	echo "Making pdf manual. Please wait ..."; texi2pdf mtt.texi

Compound-Components_rep.ps: 
	(cd $(MTT_COMPONENTS); mtt Compound-Components rep ps)
	mv $(MTT_COMPONENTS)/Compound-Components_rep.ps .

Compound-Components_rep: 
	(cd $(MTT_COMPONENTS); mtt Compound-Components rep html)
	mv $(MTT_COMPONENTS)/Compound-Components_rep .

Examples_rep.ps: 
	(cd $(MTT_EXAMPLES); mtt Examples rep ps)
	mv $(MTT_EXAMPLES)/Examples_rep.ps .

Examples_rep: 
	(cd $(MTT_EXAMPLES); mtt Examples rep html)
	mv $(MTT_EXAMPLES)/Examples_rep .

clean:
	rm -f mtt.aux	    mtt.fns	  mtt.pg	mtt.tp
	rm -f mtt.cp	    mtt.pgs	  mtt.vr	mtt.ps
	rm -f mtt.cps	    mtt.vrs	  mtt.dvi       mtt.ps.gz
	rm -f mtt.ky	    mtt.fn	  mtt.log	mtt.toc	
	rm -f 		    mtt.ps        mtt.info-*    mtt.pdf
	rm -f mtt.html	    mtt_toc.html  mtt.info
	rm -f mtt.log       mtt.ky        mtt.toc

tidy:
	rm -f mtt.aux	    mtt.fns	  mtt.pg	mtt.tp
	rm -f mtt.cp	    mtt.pgs	  mtt.vr	
	rm -f mtt.cps	    mtt.vrs	  mtt.dvi       
	rm -f mtt.log       mtt.ky        mtt.toc


install-doc: mtt.html mtt.ps.gz  mtt.info
	cp mtt.info* $(INFOPATH)
	cp mtt*.html mtt.ps.gz  $(DOCPATH)
	chmod -R a+r $(DOCPATH)
	chmod a+x  $(DOCPATH)

install-components:  Compound-Components_rep.ps Compound-Components_rep
	mv Compound-Components_rep.ps $(DOCPATH)/../components
	mv Compound-Components_rep    $(DOCPATH)/../components
	chmod -R a+r $(DOCPATH)/../components
	chmod -R a+x $(DOCPATH)/../components/Compound-Components_re




