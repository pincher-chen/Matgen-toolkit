SRCDIR = ./src
OBJDIR = ./bin

all : pre rm_mof_solvents find_space_groups in_cell ICSD_classify CSD_classify format 

pre : 
	$(shell if [ -d ${OBJDIR} ]; then echo ""; else mkdir ${OBJDIR}; fi;)

rm_mof_solvents : pre
	g++ ./src/rm_mof_solvents.cpp -o rm_mof_solvents.out -std=c++11
	@mv ./rm_mof_solvents.out ./bin

find_space_groups : pre
	g++ ./src/find_space_groups.cpp ./include/spglib/_build/libsymspg.so -o find_space_groups.out -I./include/spglib/src
	@mv ./find_space_groups.out ./bin

in_cell : pre
	g++ ./src/in_cell.cpp -o in_cell.out -std=c++11
	@mv ./in_cell.out ./bin

ICSD_classify : pre
	g++ ./src/ICSD_classify.cpp -o ICSD_classify.out -std=c++11
	@mv ./ICSD_classify.out ./bin

CSD_classify : pre
	g++ ./src/CSD_classify.cpp -o CSD_classify.out -std=c++11
	@mv ./CSD_classify.out ./bin

format : pre
	g++ ./src/format.cpp -o format.out -std=c++11
	@mv ./format.out ./bin

clean : pre
	@rm ./bin/*.out