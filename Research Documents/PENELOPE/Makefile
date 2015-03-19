######################################################
# Makefile for PENELOPE 2003
#
# Written by Kyung-Wook Shin
# July 28th, 2011
######################################################

#
# Change FCC according to your Fortran compiler.
# 
FCC = gfortran
FFLAGS = -O3 -Wall
PROGRAMS = material tables penslab pencyl penmain
PREFIX = ./
SOURCE = fsource/
PENCYL = mains/pencyl/
PENSLAB = mains/penslab/
PENDBASE = pendbase/
PENMAIN = mains/penmain/
WORK = work/
#OTHER = other/
TABLES = other/tables/

all:$(PROGRAMS)

material: ./$(SOURCE)material.f ./$(SOURCE)penelope.f
	$(FCC) $(FFLAGS) -o ./material ./$(SOURCE)material.f
tables: ./$(TABLES)tables.f ./$(SOURCE)penelope.f
	cp -vf ./$(SOURCE)penelope.f ./$(TABLES)penelope.f
	$(FCC) $(FFLAGS) -o ./tables ./$(TABLES)tables.f
penslab: ./$(SOURCE)penslab.f ./$(SOURCE)penelope.f ./$(SOURCE)penvared.f ./$(SOURCE)timer.f
	$(FCC) $(FFLAGS) -o ./penslab ./$(SOURCE)penslab.f
pencyl: ./$(PENCYL)pencyl.f ./$(SOURCE)penelope.f \
		./$(SOURCE)timer.f
	cp -vf ./$(SOURCE)penelope.f ./$(PENCYL)penelope.f
	cp -vf ./$(SOURCE)penvared.f ./$(PENCYL)penvared.f
	cp -vf ./$(SOURCE)timer.f ./$(PENCYL)timer.f
	$(FCC) $(FFLAGS) -o pencyl ./$(PENCYL)pencyl.f
penmain: ./$(PENMAIN)penmain.f ./$(SOURCE)pengeom.f \
		 ./$(SOURCE)penvared.f ./$(SOURCE)timer.f
	cp -vf ./$(SOURCE)penelope.f ./$(PENMAIN)penelope.f
	cp -vf ./$(SOURCE)pengeom.f ./$(PENMAIN)pengeom.f
	cp -vf ./$(SOURCE)penvared.f ./$(PENMAIN)penvared.f
	cp -vf ./$(SOURCE)timer.f ./$(PENMAIN)timer.f
	$(FCC) $(FFLAGS) -o penmain ./$(PENMAIN)penmain.f

install:$(PROGRAMS)
	mv -vrf ./material $(PREFIX)$(PENDBASE)
	mv -vrf ./tables $(PREFIX)$(TABLES)
	mv -vrf ./penslab $(PREFIX)$(PENSLAB)
	mv -vrf ./pencyl $(PREFIX)$(PENCYL)
	mv -vrf ./penmain $(PREFIX)$(PENMAIN)

clean:
	rm -vrf ./material
	rm -vrf ./tables ./$(TABLES)penelope.f
	rm -vrf ./penslab 
	rm -vrf ./pencyl ./$(PENCYL)penelope.f ./$(PENCYL)timer.f ./$(PENCYL)penvared.f
	rm -vrf ./penmain ./$(PENMAIN)penelope.f ./$(PENMAIN)timer.f ./$(PENMAIN)pengeom.f ./$(PENMAIN)penvared.f

uninstall:
	rm -vrf $(PREFIX)$(PENDBASE)material
	rm -vrf $(PREFIX)$(TABLES)tables ./$(TABLES)penelope.f
	rm -vrf $(PREFIX)$(PENSLAB)penslab 
	rm -vrf $(PREFIX)$(PENCYL)pencyl ./$(PENCYL)penelope.f ./$(PENCYL)timer.f ./$(PENCYL)penvared.f
	rm -vrf $(PREFIX)$(PENMAIN)penmain ./$(PENMAIN)penelope.f ./$(PENMAIN)timer.f ./$(PENMAIN)pengeom.f ./$(PENMAIN)penvared.f

#test:$(PROGRAMS)
#		@echo "Testing: material..."
#		cd $(PREFIX)$(PENDBASE)
#		./material $(PREFIX)$(PENDBASE)al.mat
#		./material $(PREFIX)$(PENDBASE)nai.mat
#		./material $(PREFIX)$(PENDBASE)al2o3.mat
#		./material $(PREFIX)$(PENDBASE)glass.mat
#		./material $(PREFIX)$(PENDBASE)polyst.mat
#		@echo "Testing: material completed."
