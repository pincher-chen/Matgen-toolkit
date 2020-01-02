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

const double Threshold = 0.01;

void get_res() noexcept(false);

double get_cif_similarity(CIF &a, CIF &b);

int main(int argc, char *argv[]) {
    // cmd
    cmdline::parser parser;
    parser.add<string>("input_dir", 'i', "cod folder location", true, "");
    parser.add<string>("output_dir", 'o', "classification result export location", true);
    parser.add<string>("log_file", 'f', "cod classification log file", true, "");
    parser.add<int>("index", '\0', "index of log file", true);
    parser.add("log", 'l', "print the detail log, no log by default");

    parser.parse_check(argc, argv);
    
    if (argc == 1 || parser.exist("help")) {
        cout << parser.usage();
        return 0;
    }

    string input = parser.get<string>("input_dir");
    string output = parser.get<string>("output_dir") + +"/cutoff-" + to_string(get_cutoff()).substr(0, 4);
    string log_file = parser.get<string>("log_file");
    int index = parser.get<int>("index");
    
    bool log = false;
    if(parser.exist("log")) {
        log = true;
    }

    if(!is_folder_exist(input)) {
        cout << "COD data not found!" << endl;
    }

    map<string, vector<string>> files;
    get_res(log);
    
    // parse log
    ifstream in(log_file, ios::in);
    if(!in.is_open()) {
        throw Exception(log_file + " does not exist or fails to open!");
    }
    else {
        string str((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());

        vector<string> blocks = re_split(str, "\n");

        // cal
        auto t = del_split(blocks[index], ':');
        string key = t[0];
        vector<string> files = re_split(t[1], "\\s+");

        string out_log = output + "/log_" + to_string(index) +".txt";
        if(!is_folder_exist(output)) {
            make_dir(output);
        }
        ofstream lout(out_log);
        lout << key << endl;

        if(files.size() <= 1) {
            lout << setw(50) << left << files[0] << setw(50) << left << files[0] << setw(20) << left << 0.0 << endl;
        }
        else {
            for(int i = 0; i < files.size(); i++) {
                CIF a_cif = CIF(input + "/" + files[i], log);
                a_cif.parse_file();
                for(int j = i + 1; j < files.size(); j++) {
                    CIF b_cif = CIF(input + "/" + files[j], log);
                    b_cif.parse_file();
                    double sim = get_cif_similarity(a_cif, b_cif);

                    lout.flags(ios::fixed);
                    lout.precision(8);
                    lout << setw(50) << left << files[i] << setw(50) << left << files[j] << setw(20) << left << sim << endl;
                }
            }
        }

        lout.close();
    }
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