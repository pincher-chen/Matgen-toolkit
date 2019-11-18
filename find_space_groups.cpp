#include <iostream>
#include "include/cmdline.h"

using namespace std;

int main(int argc, char *argv[]) {
    // cmd
    cmdline::parser parser;
    parser.add<string>("filename", 'f', "cif file name", true, "");
    parser.add("version", 'v', "return the version of spglib");
    parser.add("why", 'y', "this method is used to see roughly why  spglib failed");
    parser.add("spacegroup", 's', "internatioanl space group short symbol and number are obtained as a string");
    parser.add("symmetry", 'm', "symmetry operations are obtained as a dictionary");
    parser.add("refine", 'r', "standardized crystal structure is obtained as a tuple of lattice (a 3x3 numpy array), atomic scaled positions (a numpy array of [number_of_atoms,3]), and atomic numbers (a 1D numpy array) that are symmetrized following space group type.");

    parser.parse_check(argc, argv);

    if (argc == 1 || parser.exist("help")) {
        cout << parser.usage();
        return 0;
    }




    // try {
        
    // } 
    // catch(Exception err) {
    //     cout << err.msg << endl;
    //     return -1;
    // }

    return 0;
}