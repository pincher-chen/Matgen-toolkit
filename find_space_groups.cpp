#include <iostream>
#include "include/cmdline.h"
#include "include/cif.h"
#include "spglib.h"

using namespace std;

struct Cell {
    // n - array 1*3
    vector<vector<double>> trans_arr;

    // n - array 3*3
    vector<vector<vector<double>>> symm_arr;

    // 3 * 3
    vector<vector<double>> lattice;

    // atom coordinate, n * 3
    vector<vector<doouble>> position;

    // atom type
    vector<int> type;
};


void get_res() noexcept(false) {
    cout << "Getting some known resources..." << endl;
    get_atom_radius();

    get_known_solvent();

    get_elements();

    get_HM2Hall();

    get_Hall2Number();

    get_Number2Hall();

    get_Hall2HM();

    get_Rhomb2HexHall();

    get_AP2Number();

    get_Number2AP();

    get_SymOpsHall();
}

void set_symm_info(CIF &cif, Cell *cell) {
    vector<vector<string>> symm = cif.get_symm();
    string name_HM = cif.get_name_HM();
    string name_Hall = cif.get_name_Hall();

    if(symm.empty()) {
        if(name_Hall.empty() && name_HM.empty()) {
            cout << "P1 symmetry is assumed!" << endl;
            symm.push_back(vector<string>{"x", "y", "z"});
        }
        else if(!name_Hall.empty()) {
            symm = SymOpsHall[name_Hall];
        }
        else if(!name_HM.empty()) {
            symm = SymOpsHall[HM2Hall[name_HM]];
        }
    }

    for(auto item : symm) {
        vector<vector<double>> s(3, vector<double>(3, 0.0));
        vector<double> t(3, 0.0);
        vector<string> ex{"-x", "-y", "-z"};

        // x - y - z
        for(int i = 0; i < 3; i++) {
            if(item[i].find(ex[i]) != string::npos) {
                s[i][i] = -1.0;
            }
            else {
                s[i][i] = 1.0;
            }

            t[i] = frac2double(item[i]);
        }

        cell->symm_arr.push_back(s);
        cell->trans_arr.push_back(t);
    }

    // for(int i = 0; i < symm.size(); i++) {
    //     printVec(symm[i]);
    //     cerr << " ";

    //     cerr << "trans: ";
    //     printVec(trans_arr[i]);

    //     cerr << "Symm: ";
    //     for(auto arr : symm_arr[i]) {
    //         printVec(arr);
    //         cerr << "-";
    //     }
    //     cerr << endl;
    // }
}

void get_cell(CIF &cif, Cell *cell) {

    // get rotation and position
    set_symm_info(cif, cell);
    
    cell->lattice = cif.get_lattice();


}


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

    try {
        string filename = parser.get<string>("filename");
        CIF cif = CIF(filename);
        cif.parse_file();

        get_res();

        Cell cell;
        get_cell(cif, &cell);


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