#include <iostream>
#include "include/cmdline.h"
#include "include/cif.h"
#include "spglib.h"

using namespace std;

vector<int> get_spglib_info() {
    // spg_get_major_version, spg_get_minor_version, spg_get_micro_version
    int major_version = spg_get_major_version();
    int minor_version = spg_get_minor_version();
    int micro_version = spg_get_micro_version();
    
    SpglibError error;
    error = spg_get_error_code();

    return vector<int>{major_version, minor_version, micro_version, error};
}

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

    CIF cif = CIF(parser.get<string>("filename"));
    cif.parse_file();

    try {
        string filename = parser.get<string>("filename");

        if(parser.exist("version") || parser.exist("why")) {
            vector<int> info = get_spglib_info();
            
            // version
            cout << "The spglib version is: " << info[0] << "." << info[1] << "." << info[2] << endl;
            cout << "Error Message: " << spg_get_error_message((SpglibError)info[3]) << endl;
        }
        else if(parser.exist("spacegroup")) {

        }
    } 
    catch(Exception err) {
        cout << err.msg << endl;
        return -1;
    }

    return 0;
}