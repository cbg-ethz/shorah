SHELL = /bin/bash

# common variables
#HOSTNAME = echo hostname
hostname := $(shell hostname)
whoami := $(shell whoami)
host-type := $(shell uname)
CCFLAGS = -Wall -Wextra# -O3# -pedantic
OLIBS = 
CPP = g++
SUBDIRS = ./samtools

# DOCUMENTATION
DOXYGEN = /usr/bin/doxygen
DOXYFILE_1 = dpm_src/doxyfile

###################################
# EDIT LINES BELOW TO INCLUDE GSL #
###################################
CXLIBS = -lgsl -lgslcblas -lpthread
CFLAGS = $(CCFLAGS) # -I/opt/local/include
#XLIBS  = -L/opt/local/lib $(CXLIBS)
XLIBS = $(CXLIBS)
CPPFLAGS = $(CFLAGS)

#SRC_1 = dpm_src/dmm_sampler.cpp
OBJS_1 = dpm_src/dpm_sampler.o
EXE_1 = diri_sampler

SRC_2 = contain_src/contain.cc
EXE_2 = contain
FLAGS_2 = -g -O2 -Wall

#OBJS_3 = freqEst_src/freqEst.o
SRC_3 = freqEst_src/freqEst.cc
EXE_3 = freqEst
FLAGS_3 = -g -O2 -Wall -ffast-math

SRC_4 = b2w_src/b2w.c
EXE_4 = b2w
FLAGS_4 = -Isamtools -Lsamtools -lbam -lm -lz -lpthread

SRC_5 = filter_src/fil.c
EXE_5 = fil
FLAGS_5 = -Isamtools -Lsamtools -lbam -lm -lz $(XLIBS) 
LIB_BAM = samtools/libbam.a

%.o : %.cpp %.h data_structures.h Makefile # only used for diri_sampler
	@echo ''
	@echo '*********************************'
	@echo '   making object: $@ '
	@echo '*********************************'
	$(CPP) $(CFLAGS) $(XLIBS) -c $< -o $@

all: $(EXE_1) $(EXE_2) $(EXE_3) $(LIB_SAMTOOLS) $(EXE_4) $(EXE_5)

$(EXE_1): $(OBJS_1) Makefile #diri_sampler
	@echo ''
	@echo '*********************************'
	@echo ' making executable: $@ '
	@echo '*********************************'
	$(CPP) $(OBJS_1) -g -O2 -o $(EXE_1) $(XLIBS)

	@echo '*******************'
	@echo 'compiled for $(host-type)'

$(EXE_2): $(SRC_2) Makefile #contain
	@echo ''
	@echo '*********************************'
	@echo ' making executable: $@ '
	@echo '*********************************'
	$(CPP) $(SRC_2) -o $(EXE_2) $(FLAGS_2)

	@echo '*******************'
	@echo 'compiled for $(host-type)'

$(EXE_3): $(SRC_3) Makefile #freqEst
	@echo ''
	@echo '*********************************'
	@echo ' making executable: $@ '
	@echo '*********************************'
	$(CPP) $(SRC_3) -o $(EXE_3) $(FLAGS_3)

	@echo '*******************'
	@echo 'compiled for $(host-type)'

$(EXE_4): $(SRC_4) $(LIB_BAM) Makefile #b2w
	@echo ''
	@echo '*********************************'
	@echo ' making executable: $@ '
	@echo '*********************************'
	$(CPP) $(SRC_4) -o $(EXE_4) $(CFLAGS) $(FLAGS_4)
	
$(EXE_5): $(SRC_5) Makefile #fil
	@echo ''
	@echo '*********************************'
	@echo ' making executable: $@ '
	@echo '*********************************'
	$(CPP) $(SRC_5) -o $(EXE_5) $(FLAGS_5)

$(LIB_BAM): Makefile
	@echo ''
	@echo '*********************************'
	@echo 'Building samtools'
	@echo '*********************************'
	cd samtools && $(MAKE) libbam.a


.PHONY : clean doc

clean:
	rm -rf $(OBJS_1) $(EXE_1) $(EXE_2) $(EXE_3) $(EXE_4) $(EXE_5) $(EXE_2).dSYM $(EXE_3).dSYM *pyc ./*pyc
	for i in $(SUBDIRS); do \
	( cd $$i ; make clean ; cd ../ ) ;\
	done
doc:
	$(DOXYGEN) $(DOXYFILE_1)
