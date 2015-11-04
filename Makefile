#
# 'make'	build 'functions.so' file
# 'make clean'	removes all .so files
#

# compiler
CC = gcc
GCCVERSIONGTEQ4 := $(shell expr `gcc -dumpversion | cut -f1` \> 4.2.1)

# compiler flags
CFLAGS = -Wall -O3 -fPIC -g -pipe

# if GCC version > 4.2.1
ifeq "$(GCCVERSIONGTEQ4)" "1"
    CFLAGS +=-march=native
endif

# libraries
#LIBS = -lm
LIBS = -lgsl -lgslcblas

# telluric spectra
TELLDIR = MODELS/TELLURIC
TELLSRC = $(TELLDIR)/LBLRTM_CH4_+0.0.fits.gz \
	$(TELLDIR)/LBLRTM_CO2_+0.0.fits.gz \
	$(TELLDIR)/LBLRTM_H2O_+0.0.fits.gz \
	$(TELLDIR)/LBLRTM_N2O_+0.0.fits.gz \
	$(TELLDIR)/LBLRTM_O2_+0.0.fits.gz \
	$(TELLDIR)/LBLRTM_O3_+0.0.fits.gz
TELLOBJS = $(TELLDIR)/LBLRTM_CH4_+0.0.npy \
	$(TELLDIR)/LBLRTM_CO2_+0.0.npy \
	$(TELLDIR)/LBLRTM_H2O_+0.0.npy \
	$(TELLDIR)/LBLRTM_N2O_+0.0.npy \
	$(TELLDIR)/LBLRTM_O2_+0.0.npy \
	$(TELLDIR)/LBLRTM_O3_+0.0.npy
#MODELS/TELLURIC/LBLRTM*.npy

# Uranium Neon spectra
UNEDIR = MODELS/UNE
UNESRC = $(UNEDIR)/une.fits.gz
UNEOBJ = $(UNEDIR)/une.npy

# Laser Comb spectra
COMBDIR = MODELS/COMB
COMBSRC = $(COMBDIR)/comb.fits.gz
COMBOBJ = $(COMBDIR)/comb.npy

# c source files
SRC1 := functions.c
SRC1HEADER := armparams.h
OBJ1 := functions.so
SRC2 := sampling_cum_v3.c
OBJ2 := sampling_cum_sim.so

# object file
# OBJ := functions.so
OBJS = $(OBJ1) $(OBJ2) $(TELLOBJS) $(UNEOBJ) $(COMBOBJ)

.PHONY: all
all: $(OBJS)
	@echo "Done."

$(OBJ1): $(SRC1) $(SRC1HEADER)
	$(CC) -shared -o $(OBJ1) $(CFLAGS) $(LIBS) $(SRC1)
	@echo  "Compilation of '$(SRC1)' to '$(OBJ1)' successful."

$(OBJ2): $(SRC2)
	$(CC) -shared -o $(OBJ2) $(CFLAGS) $(LIBS) $(SRC2)
	@echo  "Compilation of '$(SRC2)' to '$(OBJ2)' successful."
	
$(TELLOBJS): $(TELLSRC)
	python plant.py tell

$(UNEOBJ): $(UNESRC)
	python plant.py une

$(COMBOBJ): $(COMBSRC)
	python plant.py comb

clean:
#	$(RM) -rf *~ src/*.so src/*.so.dSYM src/*.pyc src/*~
	$(RM) -rf *~ *.so *.so.dSYM *.pyc *~ MODELS/TELLURIC/*.npy MODELS/UNE/*.npy MODELS/COMB/*.npy
	@find . -name ".*_*" -ls -exec rm -rf {} \;
	@echo  "Folder cleaned."
