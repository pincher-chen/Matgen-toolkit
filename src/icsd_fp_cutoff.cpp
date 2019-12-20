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

void get_res() noexcept(false);

int main(int argc, char *argv[]) {
    // cmd
    cmdline::parser parser;
    parser.add<string>("input_dir", 'i', "icsd folder location", true, "");

    parser.parse_check(argc, argv);
    
    if (argc == 1 || parser.exist("help")) {
        cout << parser.usage();
        return 0;
    }


    string input = parser.get<string>("input_dir");

    bool log = false;
    if(!is_folder_exist(input)) {
        cout << "CSD data not found!" << endl;
    }

    get_res(log);
    vector<string> files = get_all_files(input, "cif");

    set<int> fd;
    for(int i = 0; i < files.size() / 5; i++) {
        int tmp = rand() % files.size();
        while(fd.count(tmp)) {
            tmp = rand() % files.size();        
        }
        fd.insert(tmp);
    }

    for(auto it = fd.begin(); it != fd.end(); it++) {
        string item = files[*it];
        try {
            CIF cif = CIF(item, log);
            cif.parse_file();

            cout << setw(40) << left << get_filename(item);
            get_fp_periodic(cif);
        }
        catch(Exception err) {
            cerr << err.msg << endl;
        }
    }
    return 0;
}
