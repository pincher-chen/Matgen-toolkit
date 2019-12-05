#include <iostream>
#include <iomanip>
#include <fstream>
#include "include/cif.h"
#include "include/cmdline.h"

using namespace std;

const string conf_dir = "./conf/";

void export_gaussion_format_result(string input, string output, CIF &cif);

int main(int argc, char *argv[]) {
    // cmd
    cmdline::parser parser;
    parser.add<string>("input", 'i', "input csd file", true, "");
    parser.add<string>("output", 'o', "output path of the conversion result ", true);
    parser.add<string>("mode", 'm', "the mode of the format conversion, default is asymmetric", false, "asymmetric");
    parser.add<string>("type", 't', "the type of the result format", false, "opt-freq");

    parser.parse_check(argc, argv);
    
    if (argc == 1 || parser.exist("help")) {
        cout << parser.usage();
        return 0;
    }

    string mode = "asymmetric";
    if(parser.exist("mode")) {
        mode = parser.get<string>("mode");
    }

    string type = "opt-freq";
    if(parser.exist("mode")) {
        type = parser.get<string>("type");
    }

    string input = parser.get<string>("input");
    string output = parser.get<string>("output");

    if(!is_folder_exist(input)) {
        cout << "CSD data not found!" << endl;
    }

    try {
        CIF cif = CIF(input);
        cif.parse_file();

        

        if(mode == "asymmetric" && type == "opt-freq") {
            export_gaussion_format_result(input, output, cif);
        }
    }
    catch(Exception err) {
        cout << err.msg << endl;
        return -1;
    }
    return 0;
}

void export_gaussion_format_result(string input, string output, CIF &cif) {
    // export
    cout << "Exporting format result..." << endl;
    if(!is_folder_exist(output)) {
        make_dir(output);
    }

    map<string, vector<double>> atom_cd = cif.get_atom_cd();
    
    string basic = conf_dir + "opt-freq.conf";
    ifstream in(basic, ios::in);
    string str((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    in.close();

    string name = output + "/" + del_split(get_filename(input), '.')[0] + ".gjf";
    ofstream out(name);
    out << str << endl;
    for(auto iter = atom_cd.begin(); iter != atom_cd.end(); iter++) {
        string item = del_split(iter->first, '-')[1];
        vector<double> cart_cd = cif.frac2cart(iter->second);

        // out.setf(ios::showpoint);
        out.flags(ios::fixed);
        out.precision(8);
        out << setw(3) << left << " " + item << setw(20) << right << cart_cd[0] << setw(15) << right << cart_cd[1] << setw(15) << right << cart_cd[2] << endl;
    }
    out << endl;
    out.close();
}
