#include <iostream>
#include <cctype>
#include <iomanip>
#include "../include/cif.h"
#include "../include/cif_func.h"
#include "../include/cmdline.h"
#include "../include/conf.h"

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
        cout << "COD data not found!" << endl;
    }

    map<string, vector<string>> fList;

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
            
            string key = f1 + "|" + f2 + "|" + f3;
            fList[key].push_back(get_filename(item));
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

    string log_file = output + "/" + "log_classification.txt";
    ofstream lout(log_file);

    for(auto i = fList.begin(); i != fList.end(); i++) {
        lout.flags(ios::fixed);
        lout.precision(8);
        lout << setw(30) << left << i->first << ":";
        for(auto item : i->second) {
            lout << setw(30) << left << item;
        }
        lout << endl;
    }

    lout.close();


    return 0;
}