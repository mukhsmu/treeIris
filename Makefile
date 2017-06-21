# Makefile
BASEDIR	      = $(shell pwd)/..
INSTALLDIR    = $(BASEDIR)/treeIris
INCLUDEDIR    = $(INSTALLDIR)/include
LIBDIR        = $(INSTALLDIR)/lib
BINARYDIR     = $(INSTALLDIR)/bin
SOURCEDIR     = $(INSTALLDIR)/src
OBJECTDIR     = $(INSTALLDIR)/obj

HEADER = -I$(INCLUDEDIR) 

CXX = g++
LD = g++

# ROOT library
ROOTLIBS  = $(shell $(ROOTSYS)/bin/root-config --libs)  -lXMLParser -lThread -Wl,-rpath,$(ROOTSYS)/lib
ROOTGLIBS = $(shell $(ROOTSYS)/bin/root-config --glibs) -lXMLParser -lThread -Wl,-rpath,$(ROOTSYS)/lib

# Path of ROOT analyzer
ROOTANA = $(HOME)/packages/rootana

CXXFLAGS += -g -O -Wall -Wuninitialized -I./  -I$(INCLUDEDIR) -I$(ROOTSYS)/include  -I$(ROOTANA) -I$(ROOTANA)/include

ROOTCFLAGS    = $(shell root-config --cflags)
CXXFLAGS      += -g -Wall -ansi -Df2cFortran -fPIC $(ROOTCFLAGS) 

ANAOBJECTS  =  $(OBJECTDIR)/HandleMesytec.o $(OBJECTDIR)/HandleV1190.o $(OBJECTDIR)/HandleSTAT.o $(LIBDIR)/libTEvent.so $(OBJECTDIR)/TEventDict.o $(OBJECTDIR)/HandleScaler.o  $(OBJECTDIR)/CalibMesytec.o $(OBJECTDIR)/geometry.o 

SOFLAGS       = -g -shared
LDFLAGS	      = -O2	

all: $(BINARYDIR)/treeIris $(LIBDIR)/libTEvent.so

$(BINARYDIR)/treeIris: $(ANAOBJECTS) $(OBJECTDIR)/treeIris.o $(ROOTANA)/lib/librootana.a 
	$(CXX) -o $@ $(CXXFLAGS) $^ $(ROOTGLIBS) -lm -lz -lutil -lnsl -lpthread -lrt

$(LIBDIR)/libTEvent.so: $(OBJECTDIR)/IDet.o $(OBJECTDIR)/ITdc.o $(OBJECTDIR)/IScaler.o $(OBJECTDIR)/TEventDict.o
	$(LD) $(SOFLAGS) $(LDFLAGS) $^ -o $@
	@echo "$@ done"

$(OBJECTDIR)/treeIris.o: $(SOURCEDIR)/main.cxx 
	$(CXX) $(CXXFLAGS) -c $< -o $@ 

$(OBJECTDIR)/HandleMesytec.o: $(SOURCEDIR)/HandleMesytec.cxx 
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJECTDIR)/HandleV1190.o: $(SOURCEDIR)/HandleV1190.cxx 
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJECTDIR)/HandleSTAT.o: $(SOURCEDIR)/HandleSTAT.cxx 
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJECTDIR)/HandleScaler.o: $(SOURCEDIR)/HandleScaler.cxx 
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJECTDIR)/CalibMesytec.o: $(SOURCEDIR)/CalibMesytec.cxx 
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJECTDIR)/geometry.o: $(SOURCEDIR)/geometry.cxx 
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJECTDIR)/IDet.o: $(SOURCEDIR)/IDet.cxx 
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJECTDIR)/ITdc.o: $(SOURCEDIR)/ITdc.cxx 
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJECTDIR)/IScaler.o: $(SOURCEDIR)/IScaler.cxx 
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJECTDIR)/TEventDict.o: $(LIBDIR)/TEventDict.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(LIBDIR)/TEventDict.cxx:  $(INCLUDEDIR)/IDet.h $(INCLUDEDIR)/ITdc.h $(INCLUDEDIR)/IScaler.h $(INCLUDEDIR)/TEventLinkDef.h
	@echo "Generating dictionary $@..."
	@rootcint -f $@ -c $(HEADER) $^

clean::
	rm -f $(OBJECTDIR)/*.o 
	rm -f $(LIBDIR)/*Dict.cxx 
	rm -f $(LIBDIR)/*Dict.cpp 
	rm -f $(LIBDIR)/*Dict.h 
	rm -f $(LIBDIR)/*.pcm 
	rm -f core 
	rm -f $(SOURCEDIR)/*~ $(INCLUDEDIR)/*~ 
	rm -f $(LIBDIR)/libTEvent.so
	rm -f $(INSTALLDIR)/bin/treeIris
# end
