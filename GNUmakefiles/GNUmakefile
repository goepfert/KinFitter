###############################################################
#
# $Id: GNUmakefile,v 1.5 2008/06/02 09:18:40 petzold Exp $
#
# Makefile for shared libs with dictionary to be loaded in root
# 
###############################################################

# file extensions
# extensions can be defined for each package in the package GNUmakefile

ifndef CXXHDREXT
	CXXHDREXT = hh
endif

ifndef CXXSRCEXT
	CXXSRCEXT = cc
endif

ifndef CXXHDRDIR
	CXXHDRDIR = ./
endif

ifndef CXXSRCDIR
	CXXSRCDIR = ./
endif

ifndef DEPS
	DEPS = KinFitter
endif


# Directories
OBJDIRBASE    = $(SRC)/lib/$(BFARCH)/$(ROOTRELEASE)
LIBDIR        = $(SRC)/lib/$(BFARCH)/$(ROOTRELEASE)
TEMPLATEDIR   = $(OBJDIRBASE)/$(PACKAGE)

# ROOT - Stuff
ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
ROOTLDFLAGS   = $(shell root-config --ldflags)

# Fortran Compiler Flags
FC = g77
FCFLAGS = -g -fno-automatic -fvxt -fdollar-ok -fno-backslash \
          -ffixed-line-length-132 -fno-second-underscore

# C++ - Compiler Flags
CXX = g++
CXXFLAGS = -g  -w -O0  -fno-exceptions -fPIC -I.. -I. 

# BABAR specific settings for TUD. AFAIK nobody uses this but 
# it's better to keep it for backward compatibility
ifeq "tud" "$(BFSITE)"
   CXX = g++
   CXXFLAGS = -g  -w -O0 -DBABAR_COMP_INST -fno-exceptions -fPIC -I.. -I.
endif

LD            =  $(CXX)
LDFLAGS      += $(ROOTLDFLAGS) -g
SOFLAGS       = -shared
CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) -lEG -lProof -lTreePlayer
GLIBS         = $(ROOTGLIBS) -lEG -lProof -lTreePlayer -lMinuit

###################################################################################

# cint dicttionary settings
DICT    = $(PACKAGE)_Cint
DICTO   = $(DICT).o
DICTCXX = $(DICT).$(CXXSRCEXT)

# List of lib and object directories:
LIBA	= libT$(PACKAGE).a
LIBSO	= libT$(PACKAGE).so
OBJDIR  = $(OBJDIRBASE)/$(PACKAGE)

# Lists of objects NOT to include in library:
NOTTOCOMPILE += $(DICT)

# c++ files
ALLBASE   = $(notdir $(wildcard $(CXXSRCDIR)/*.$(CXXSRCEXT)) )
ALL       = $(basename $(ALLBASE) )
NAMELIST  = $(filter-out  $(NOTTOCOMPILE),$(ALL))
HDRS      = $(addsuffix .$(CXXHDREXT), $(NAMELIST) )

# filter out hdrs that should not get into the rootcint dict
# NOTINDICT is usually defined in the PACKAGE GNUmakefile
DICTHDRS  = $(filter-out  $(NOTINDICT), $(HDRS))
OO 	  = $(addprefix $(OBJDIR)/, $(NAMELIST))
OBJS      = $(addsuffix .o, $(OO) )

# Fortran files
ALLF      = $(basename $(wildcard *.F) )
NAMELISTF = $(filter-out  $(NOTTOCOMPILE),$(ALLF))
OOF 	  = $(addprefix $(OBJDIR)/, $(NAMELISTF))
OBJSF     = $(addsuffix .o, $(OOF) )


# Default target build shared lib
all,shlib: $(LIBSO)

statlib: $(LIBA)

info:
	@echo "All"
	@echo "$(ALL)"
	@echo "NAMELIST"
	@echo "$(NAMELIST)"
	@echo "HDRS"
	@echo "$(HDRS)"
	@echo "OBJS"
	@echo "$(OBJS)"

# Static library:
$(LIBA) : $(OBJSF) $(OBJS) $(OBJDIR)/$(DICTO)
	@if [ ! -d $(LIBDIR)/$(PACKAGE) ] ; then echo -e "\n   !!!! ERROR: $(LIBDIR)/$(PACKAGE) doesn't exist!\nMaybe you should run 'make installdirs'!\n"; exit 1; fi
	rm -f $@
	$(AR) $(ARFLAGS)  $(LIBDIR)/$@ $^

# Shared library:
$(LIBSO) : $(OBJSF) $(OBJS) $(OBJDIR)/$(DICTO)
	@if [ ! -d $(LIBDIR)/$(PACKAGE) ] ; then echo -e "\n   !!!! ERROR: $(LIBDIR)/$(PACKAGE)  doesn't exist!\nMaybe you should run 'make installdirs'!\n"; fi
	@echo "   >>>> Linking shared library $@"
	$(LD) -g $(EXTRALIBS) $(SOFLAGS) $(GLIBS) -o $(LIBDIR)/$@ $^

# Rules for Dictionary:
$(OBJDIR)/$(DICTO) : $(DICTCXX) 
	@echo "   >>>> Compiling dictionary $(DICTCXX)"
	@if [ ! -d $(OBJDIR) ] ; then echo -e "\n   !!!! ERROR: $(OBJDIR) doesn't exist!\nMaybe you should run 'make installdirs'!\n"; fi
	$(CXX) $(CXXFLAGS) -I$(PACKAGE) -c $< -o $(OBJDIR)/$(DICTO)


$(DICTCXX) : $(DICTHDRS) $(PACKAGE)_LinkDef.$(CXXHDREXT)
	@echo "   >>>> Creating dictionary $(DICTCXX)"
	rootcint  -f $@ -c  -I.. -I$(PACKAGE) $^


# Static pattern rule for object file dependency on sources:
$(OBJS) : $(OBJDIR)/%.o : %.$(CXXSRCEXT) %.$(CXXHDREXT)
	@if [ ! -d $(OBJDIR) ] ; then echo -e "\n   !!!! ERROR: $(OBJDIR) doesn't exist!\nMaybe you should run make installdirs!\n"; exit 1 ; fi
	@echo "   >>>> Compiling $<";
	$(CXX) $(CXXFLAGS) -I$(DEPS) -c $< -o $@

# Static pattern rule for object file dependency on Fortran sources
$(OBJSF) : $(OBJDIR)/%.o : %.F
	@echo "   >>>> Compiling $<";
	$(FC) $(CXXFLAGS)  $< -c -o $@


# create required directories for libs etc.
installdirs : 
	@if [ ! -d $(SRC)/lib/$(BFARCH) ] ; then mkdir -p $(SRC)/lib/$(BFARCH); fi
	@if [ ! -d $(OBJDIRBASE) ]; then mkdir -p $(OBJDIRBASE); echo "   >>>> Create $(OBJDIRBASE)"; fi
	@if [ ! -d $(OBJDIR) ]; then mkdir -p $(OBJDIR); echo "   >>>> Create $(OBJDIR)"; \
        else echo " $(OBJDIR) already existent"; fi	

install :
	@echo "   >>>> The 'install' target is deprecated! Use 'make installdirs' instead!";

clean : 
	rm -f $(OBJDIRBASE)/$(LIBA) \
        $(OBJDIRBASE)/$(LIBSO) \
	$(OBJDIRBASE)/$(PACKAGE)/*.o \
	$(SRC)/$(PACKAGE)/$(PACKAGE)_Cint.$(CXXSRCEXT) \
	$(SRC)/$(PACKAGE)/$(PACKAGE)_Cint.$(CXXHDREXT)

