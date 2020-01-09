SRCDIR = ./src
OBJDIR = ./bin

all : pre rm_mof_solvents find_space_groups in_cell icsd_classify csd_classify format 

pre : 
	$(shell if [ -d ${OBJDIR} ]; then echo ""; else mkdir ${OBJDIR}; fi;)

rm_mof_solvents : pre
	g++ ${SRCDIR}/rm_mof_solvents.cpp -o rm_mof_solvents -std=c++11
	@mv ./rm_mof_solvents ${OBJDIR}

find_space_groups : pre
	g++ ${SRCDIR}/find_space_groups.cpp ./include/spglib/_build/libsymspg.so -o find_space_groups -I./include/spglib/src -std=c++11
	@mv ./find_space_groups ${OBJDIR}

in_cell : pre
	g++ ${SRCDIR}/in_cell.cpp -o in_cell -std=c++11
	@mv ./in_cell ${OBJDIR}

occupancy_filter : pre
	g++ ${SRCDIR}/cif_occupancy_filter.cpp -o occupancy_filter -std=c++11
	@mv ./occupancy_filter ${OBJDIR}

icsd_classify : pre
	g++ ${SRCDIR}/icsd_classify.cpp -o icsd_classify -std=c++11
	@mv ./icsd_classify ${OBJDIR}

icsd_classify_fp : pre
	g++ ${SRCDIR}/icsd_classify_fp.cpp -o icsd_classify_fp -std=c++11 -L ./lib -static -llapack -lblas -ltmglib -lf2c
	@mv ./icsd_classify_fp ${OBJDIR}

csd_classify : pre
	g++ ${SRCDIR}/csd_classify.cpp -o csd_classify -std=c++11
	@mv ./csd_classify ${OBJDIR}

cod_classify : pre
	g++ ${SRCDIR}/cod_classify.cpp -o cod_classify -std=c++11
	@mv ./cod_classify ${OBJDIR}

format : pre
	g++ ${SRCDIR}/format.cpp -o format -std=c++11
	@mv ./format ${OBJDIR}

splice_molecule : pre
	g++ ${SRCDIR}/splice_molecule.cpp -o splice_molecule -std=c++11
	@mv ./splice_molecule ${OBJDIR}

clean : pre
	@rm ${OBJDIR}/*