#include <iostream>
#include <cctype>
#include "./include/cif.h"
#include "./include/cmdline.h"

using namespace std;

const double Threshold = 0.1;

void get_res() noexcept(false);

double MSD(CIF &a, CIF &b);

int main(int argc, char *argv[]) {
    // cmd
    cmdline::parser parser;
    parser.add<string>("input_dir", 'i', "icsd folder location", true, "");
    parser.add<string>("output_dir", 'o', "classification result export location", true);

    parser.parse_check(argc, argv);
    
    if (argc == 1 || parser.exist("help")) {
        cout << parser.usage();
        return 0;
    }

    string input = parser.get<string>("input_dir");
    string output = parser.get<string>("output_dir");

    if(!is_folder_exist(input)) {
        cout << "CSD data not found!" << endl;
    }

    get_res();
    vector<string> files = get_all_files(input, "cif");
    for(auto item : files) {
        try {
            cout << "-----------------------FILE " + item + " -----------------------" << endl;
            CIF cif = CIF(item);
            cif.parse_file();

            vector<string> cmp_files = get_all_files(input, "cif");

            bool flag = false;
            for(auto cmp : cmp_files) {
                if(cmp == item) {
                    flag = true;
                    continue;
                }

                if(!flag) {
                    continue;
                }
                CIF other = CIF(cmp);
                other.parse_file();

                double msd = MSD(cif, other);

                cerr << item + " -- " + cmp << " : \t\t" << msd << endl;
            }
        }
        catch(Exception err) {
            cout << err.msg << endl;
        }
    }
    return 0;
}


void get_res() noexcept(false) {
    cout << "Getting some known resources..." << endl;
    get_atom_radius();

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

double MSD(CIF &a, CIF &b) {
    double msd = 0;
    map<string, vector<double>> a_atom_cd = a.get_atom_cd();
    map<string, vector<double>> b_atom_cd = b.get_atom_cd();

    if(a_atom_cd.size() != b_atom_cd.size()) {
        return -1;
    }

    for(auto a_iter = a_atom_cd.begin(), b_iter = b_atom_cd.begin(); a_iter != a_atom_cd.end() && b_iter != b_atom_cd.end(); a_iter++, b_iter++) {
        if(del_split(a_iter->first, '-')[1] != del_split(b_iter->first, '-')[1]) {
            return -1;
        }

        vector<double> a_cart_cd = a.frac2cart(a_iter->second);
        vector<double> b_cart_cd = b.frac2cart(b_iter->second);

        msd += pow(a_cart_cd[0] - b_cart_cd[0], 2) + 
                pow(a_cart_cd[1] - b_cart_cd[1], 2) + 
                pow(a_cart_cd[2] - b_cart_cd[2], 2);
    }

    msd = msd / a_atom_cd.size();

    double rmsd = sqrt(msd);

    return msd;
}