#include <iostream>
#include <cctype>
#include "../include/cif.h"
#include "../include/cif_func.h"
#include "../include/cmdline.h"
#include "../include/conf.h"

using namespace std;

void get_res() noexcept(false);

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
    string output = parser.get<string>("output_dir") + "/";
    
    bool log = false;
    if(parser.exist("log")) {
        log = true;
    }

    if(!is_folder_exist(input)) {
        cout << "CSD data not found!" << endl;
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

            cp_file(item, output);
        }
        catch(Exception err) {
            if(log) {
                cerr << err.msg << endl;
            }
        }
    }
    return 0;
}