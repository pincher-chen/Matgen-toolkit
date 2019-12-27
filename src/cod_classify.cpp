#include <iostream>
#include <cctype>
#include <iomanip>
#include "../include/cif.h"
#include "../include/cif_func.h"
#include "../include/cif_sim.h"
#include "../include/cmdline.h"
#include "../include/conf.h"
#include "../include/fp/fplib.h"

using namespace std;

const double Threshold = 0.1;

void get_res() noexcept(false);

double get_cif_similarity(CIF &a, CIF &b);

int main(int argc, char *argv[]) {
    // cmd
    cmdline::parser parser;
    parser.add<string>("input_dir", 'i', "icsd folder location", true, "");
    parser.add<string>("output_dir", 'o', "classification result export location", true);
    parser.add("log", 'l', "print the detail log, no log by default");

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
    
    string log_file = output + "/" + "file_similarity_log_" + to_string(get_cutoff()).substr(0, 4) + ".txt";
    ofstream lout(log_file);

    map<string, set<string>> sim;

    if(!is_folder_exist(output)) {
        make_dir(output);
    }


    bool log = false;
    if(parser.exist("log")) {
        log = true;
    }

    get_res(log);
    vector<string> files = get_all_files(input, "cif");
    for(auto item : files) {
        try {

            if(log) {
                cout << "-----------------------FILE " + item + " -----------------------" << endl;
            }

            CIF cif = CIF(item, log);
            cif.parse_file();

            // component/element type/space group/
            string f1 = "", f2 = "", f3 = "";

            set<string> atoms = cif.get_atoms();
            f1 = to_string(atoms.size());
            for(auto iter = atoms.begin(); iter != atoms.end(); iter++) {
                f2 += (*iter) + "-";
            }
            f2 = f2.substr(0, f2.size()-1);

            string name_Hall = cif.get_name_Hall();
            string name_HM = cif.get_name_HM();
            if(!name_Hall.empty()) {
                f3 = to_string(Hall2Number[name_Hall]);
            }
            else if(!name_HM.empty()) {
                f3 = to_string(Hall2Number[HM2Hall[name_HM]]);
            }
            else {
                cerr << item << " space group not found!" << endl;
                continue;
            }
            
            string base_path = output + "/" + f1 + "/" + f2 + "/" + f3 + "/";

            if(!is_folder_exist(base_path)) {
                make_dir(base_path);
                cp_file(item, base_path);
            }
            else {
                vector<string> files = get_all_files(base_path, "cif");
                
                for(auto cmp_item : files) {
                    CIF other = CIF(cmp_item, log);
                    other.parse_file();
                    
                    double sim = get_cif_similarity(cif, other);

                    if(sim < 0) {
                        continue;
                    }
                    else {
                        string a = get_filename(item);
                        string b = get_filename(cmp_item);

                        lout.flags(ios::fixed);
                        lout.precision(8);
                        lout << setw(50) << left << a << setw(50) << left << b << setw(20) << left << sim << endl;
                    }
                }

                cp_file(item, base_path);
            }
        }
        catch(Exception err) {
            if(log) {
                cerr << err.msg << endl;
            }

            if(err.msg.find("Error:") != string::npos) {
                cp_file(item, output + "/error_file/");
            }
        }
    }
    lout.close();
    return 0;
}

double get_cif_similarity(CIF &a, CIF &b) {
    map<string, vector<double>> a_atom_cd = a.get_atom_cd();
    map<string, vector<double>> b_atom_cd = b.get_atom_cd();

    if(a_atom_cd.size() != b_atom_cd.size()) {
        return -1;
    }

    for(auto a_iter = a_atom_cd.begin(), b_iter = b_atom_cd.begin(); a_iter != a_atom_cd.end() && b_iter != b_atom_cd.end(); a_iter++, b_iter++) {
        if(del_split(a_iter->first, '-')[1] != del_split(b_iter->first, '-')[1]) {
            return -1;
        }
    }

    return get_fp_similarity(a, b);
}