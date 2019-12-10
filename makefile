SRCDIR = ./src
OBJDIR = ./bin

all : pre rm_mof_solvents find_space_groups in_cell ICSD_classify CSD_classify format 

pre : 
	$(shell if [ -d ${OBJDIR} ]; then echo ""; else mkdir ${OBJDIR}; fi;)

rm_mof_solvents : pre
	g++ ${SRCDIR}/rm_mof_solvents.cpp -o rm_mof_solvents -std=c++11
	@mv ./rm_mof_solvents ${OBJDIR}

find_space_groups : pre
	g++ ${SRCDIR}/find_space_groups.cpp ./include/spglib/_build/libsymspg.so -o find_space_groups -I./include/spglib/src
	@mv ./find_space_groups ${OBJDIR}

in_cell : pre
	g++ ${SRCDIR}/in_cell.cpp -o in_cell -std=c++11
	@mv ./in_cell ${OBJDIR}

ICSD_classify : pre
	g++ ${SRCDIR}/ICSD_classify.cpp -o ICSD_classify -std=c++11
	@mv ./ICSD_classify ${OBJDIR}

CSD_classify : pre
	g++ ${SRCDIR}/CSD_classify.cpp -o CSD_classify -std=c++11
	@mv ./CSD_classify ${OBJDIR}

format : pre
	g++ ${SRCDIR}/format.cpp -o format -std=c++11
	@mv ./format ${OBJDIR}

clean : pre
	@rm ${OBJDIR}/*