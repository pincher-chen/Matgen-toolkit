#include <iostream>
#include <iomanip>
#include "../include/cmdline.h"
#include "../include/cif.h"
#include "../include/cif_func.h"
#include "spglib.h"

#define SYMPREC 1e-5

using namespace std;

void set_symm_info(CIF &cif, Cell *cell);

void get_cell(CIF &cif, Cell *cell);

vector<int> get_spglib_info();

void get_spacegroup(CIF &cif, Cell *cell);

void get_symmetry(Cell *cell);

void refine_cell(Cell *cell);

void find_primitive(Cell *cell);

void get_symmetry_dataset(Cell *cell);

void get_symm_from_database(Cell *cell);

void get_spacegroup_type(Cell *cell);

void niggli_reduce(Cell *cell);

void delaunay_reduce(Cell *cell);

void get_ir_reciprocal_mesh(Cell *cell);

void show_cell(double lattice[3][3], double position[][3], const int types[], const int atom_num);

void show_spg_dataset(double lattice[3][3], double position[][3], const int atom_num, const int types[]);


int main(int argc, char *argv[]) {
    // cmd
    cmdline::parser parser;
    parser.add<string>("input", 'i', "input cif file name", true, "");
    parser.add("version", 'v', "return the version of spglib");
    parser.add("why", 'w', "this method is used to see roughly why  spglib failed");
    parser.add("spacegroup", 's', "internatioanl space group short symbol and number are obtained as a string");
    parser.add("symmetry", 'm', "symmetry operations are obtained as a dictionary");
    parser.add("refine", 'r', "standardized crystal structure is obtained as a tuple of lattice (a 3x3 numpy array), atomic scaled positions (a numpy array of [number_of_atoms,3]), and atomic numbers (a 1D numpy array) that are symmetrized following space group type.");
    parser.add("primitive", 'p', "is found, lattice parameters (a 3x3 numpy array), scaled positions (a numpy array of [number_of_atoms,3]), and atomic numbers (a 1D numpy array) is returned.");
    parser.add("dataset", 'd', "dataset,cell and symprec;angle_tolerance;hall_number;number;choice;transformation_matrix;origin shift;wyckoffs;site_symmetry_symbols;equivalent_atoms;mapping_to_primitive;rotations and translations;pointgroup;std_lattice;std_positions;std_types;std_rotation_matrix;std_mapping_to_primitive");
    parser.add("symmfdset", 'c', "A set of crystallographic symmetry operations corresponding to hall_number is returned by a dictionary where rotation parts and translation parts are accessed by the keys rotations and translations, respectively.");
    parser.add("spgfdset", 'f', "This function allows to directly access to the space-group-type database in spglib (spg_database.c). A dictionary is returned. To specify the space group type with a specific choice, hall_number is used.");
    parser.add("niggli", 'n', "Niggli reduction is achieved using this method.");
    parser.add("delaunay", 'l', "Delaunay reduction is achieved using this method.");
    parser.add("irrkpoints", 'k', "Irreducible k-points are obtained from a sampling mesh of k-points");

    parser.parse_check(argc, argv);

    if (argc == 1 || parser.exist("help")) {
        cout << parser.usage();
        return 0;
    }

    try {
        string filename = parser.get<string>("input");
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
            get_spacegroup(cif, &cell);
        }
        else if(parser.exist("symmetry")) {
            get_symmetry(&cell);
        }
        else if(parser.exist("refine")) {
            refine_cell(&cell);
        }
        else if(parser.exist("primitive")) {
            find_primitive(&cell);
        }
        else if(parser.exist("dataset")) {
            get_symmetry_dataset(&cell);
        }
        else if(parser.exist("symmfdset")) {
            get_symm_from_database(&cell);
        }
        else if(parser.exist("spgfdset")) {
            get_spacegroup_type(&cell);
        }
        else if(parser.exist("niggli")) {
            niggli_reduce(&cell);
        }
        else if(parser.exist("delaunay")) {
            delaunay_reduce(&cell);
        }
        else if(parser.exist("irrkpoints")) {
            get_ir_reciprocal_mesh(&cell);
        }
        else {
            cout << "not support" << endl;
        }
    } 
    catch(Exception err) {
        cout << err.msg << endl;
        return -1;
    }

    return 0;
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
    
    vector<vector<double>> tmp = cif.get_lattice();
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            cell->lattice[i][j] = tmp[i][j];
        }
    }

    map<string, vector<double>> atom_cd = cif.get_atom_cd();
    for(auto iter = atom_cd.begin(); iter != atom_cd.end(); iter++) {
        string specie = cif.get_species(iter->first);
        cell->position.push_back(iter->second);
        auto atom = radius_dict.find(specie);
        if(atom == radius_dict.end()) {
            atom = radius_dict.find(Universal);
            // throw Exception(specie + " 's information not found");
        }
        cell->types.push_back(get_num(atom->second.O_metal));
    }

    cell->atom_num = atom_cd.size();
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

void get_spacegroup(CIF &cif, Cell *cell) {
    string name_Hall = cif.get_name_Hall();
    string name_HM = cif.get_name_HM();
    if(!name_Hall.empty() && name_Hall != "P 1") {
        cout << "The space group is: ";
        cout << Hall2HM[name_Hall] << " " << Hall2Number[name_Hall] << endl;
        return;
    }

    if(!name_HM.empty() && name_HM != "P1") {
        cout << "The space group is: ";
        cout << name_HM << " " << Hall2Number[HM2Hall[name_HM]] << endl;
        return;
    }
    
    double position[cell->atom_num][3];
    for(int i = 0; i < cell->position.size(); i++) {
        for(int j = 0; j < 3; j++) {
            position[i][j] = cell->position[i][j];
        }
    }

    int types[cell->atom_num];
    for(int i = 0; i < cell->types.size(); i++) {
        types[i] = cell->types[i];
    }

    char symbol[21];

    int num_spg = spg_get_international(symbol, cell->lattice, position, types, cell->atom_num, 1e-5);
    if(num_spg > 0) {
        cout << "The space group is: ";
        int i = 0;
        while(symbol[i] != '\0') {
            cout << symbol[i++];
        }
        cout << " " << num_spg << endl;
    } else {
        cout << "Space group could not be found." << endl;
    }
}

void get_symmetry(Cell *cell) {
    int max_size = 200;
    int rotation[max_size][3][3];
    double translation[max_size][3];

    double position[cell->atom_num][3];
    for(int i = 0; i < cell->position.size(); i++) {
        for(int j = 0; j < 3; j++) {
            position[i][j] = cell->position[i][j];
        }
    }

    int types[cell->atom_num];
    for(int i = 0; i < cell->types.size(); i++) {
        types[i] = cell->types[i];
    }

    int size = spg_get_symmetry(rotation,
                                translation,
                                max_size,
                                cell->lattice,
                                position,
                                types,
                                cell->atom_num,
                                1e-5);
    
    
    if(size >  0) {
        cout << "The symmetry is: " << endl;
    }
    else {
        cout << "The symetry is not found" << endl;
    }
    for(int i = 0; i < size; i++) {
        cout << "---- " << i+1 << " ----" << endl;
        for(int j = 0; j < 3; j++) {
            cout << rotation[i][j][0] << " "
                 << rotation[i][j][1] << " "
                 << rotation[i][j][2] << " " << endl;
        }
        cout << translation[i][0] << " " << translation[i][1] << " " << translation[i][2] << endl;
    }
}

void refine_cell(Cell *cell) {
    double position[cell->atom_num][3];
    for(int i = 0; i < cell->position.size(); i++) {
        for(int j = 0; j < 3; j++) {
            position[i][j] = cell->position[i][j];
        }
    }

    int types[cell->atom_num];
    for(int i = 0; i < cell->types.size(); i++) {
        types[i] = cell->types[i];
    }

    int num_atom_bravais = spg_refine_cell(cell->lattice,
                                           position,
                                           types,
                                           cell->atom_num,
                                           1e-5);
    show_cell(cell->lattice, position, types, num_atom_bravais);
}

void find_primitive(Cell *cell) {
    double position[cell->atom_num][3];
    for(int i = 0; i < cell->position.size(); i++) {
        for(int j = 0; j < 3; j++) {
            position[i][j] = cell->position[i][j];
        }
    }

    int types[cell->atom_num];
    for(int i = 0; i < cell->types.size(); i++) {
        types[i] = cell->types[i];
    }

    int num_primitive_atom = spg_find_primitive(cell->lattice, position, types, cell->atom_num, SYMPREC);
    if(num_primitive_atom == 0) {
        cout << "Primitive cell was not found." << endl;
    }
    else {
        show_cell(cell->lattice, position, types, num_primitive_atom);
    }
}

void get_symmetry_dataset(Cell *cell) {
    double position[cell->atom_num][3];
    for(int i = 0; i < cell->position.size(); i++) {
        for(int j = 0; j < 3; j++) {
            position[i][j] = cell->position[i][j];
        }
    }

    int types[cell->atom_num];
    for(int i = 0; i < cell->types.size(); i++) {
        types[i] = cell->types[i];
    }

    show_spg_dataset(cell->lattice, position, cell->atom_num, types);
}

void get_symm_from_database(Cell *cell) {
    double position[cell->atom_num][3];
    for(int i = 0; i < cell->position.size(); i++) {
        for(int j = 0; j < 3; j++) {
            position[i][j] = cell->position[i][j];
        }
    }

    int types[cell->atom_num];
    for(int i = 0; i < cell->types.size(); i++) {
        types[i] = cell->types[i];
    }

    SpglibDataset *dataset = spg_get_dataset(cell->lattice,
                                             position,
                                             types,
                                             cell->atom_num,
                                             SYMPREC);

    int hall_num = dataset->hall_number;
    spg_free_dataset(dataset);

    int rotations[192][3][3];
    double translations[192][3];
    int size = spg_get_symmetry_from_database(rotations,
                                              translations,
                                              hall_num);
    for(int i = 0; i < size; i++) {
        cout << "---- " << i+1 << " ----" << endl;
        for (int j = 0; j < 3; j++) {
            cout << rotations[i][j][0] << " " << rotations[i][j][1] << " " << rotations[i][j][2] << endl;
        }
        cout << translations[i][0] << " " << translations[i][1] << " " << translations[i][2] << endl;
    }   
}

void get_spacegroup_type(Cell *cell) {
    double position[cell->atom_num][3];
    for(int i = 0; i < cell->position.size(); i++) {
        for(int j = 0; j < 3; j++) {
            position[i][j] = cell->position[i][j];
        }
    }

    int types[cell->atom_num];
    for(int i = 0; i < cell->types.size(); i++) {
        types[i] = cell->types[i];
    }

    SpglibDataset *dataset = spg_get_dataset(cell->lattice,
                                             position,
                                             types,
                                             cell->atom_num,
                                             SYMPREC);

    int hall_num = dataset->hall_number;
    spg_free_dataset(dataset);

    SpglibSpacegroupType spgtype = spg_get_spacegroup_type(hall_num);

    cout << "Number:\t\t" << spgtype.number << endl;
    cout << "Schoenflies:\t\t" << spgtype.schoenflies << endl;
    cout << "International:\t\t" << spgtype.international << endl;
    cout << "International:\t\t" << spgtype.international_full << endl;
    cout << "International:\t\t" << spgtype.international_short << endl;
    cout << "Hall symbol:\t\t" << spgtype.hall_symbol << endl;
}

void niggli_reduce(Cell *cell) {
    // for(int i = 0; i < 3; i++) {
    //     cout << cell->lattice[i][0] << " " << cell->lattice[i][1] << " " << cell->lattice[i][2] << endl;
    // } 
    int niggli_lattice = spg_niggli_reduce(cell->lattice, SYMPREC);
    if(niggli_lattice == 0) {
        cout << "Niggli reduction failed" << endl;
    }
    cout << "The niggli lattice is: " << endl;
    for(int i = 0; i < 3; i++) {
        cout << cell->lattice[i][0] << " " << cell->lattice[i][1] << " " << cell->lattice[i][2] << endl;
    }
}

void delaunay_reduce(Cell *cell) {
    // for(int i = 0; i < 3; i++) {
    //     cout << cell->lattice[i][0] << " " << cell->lattice[i][1] << " " << cell->lattice[i][2] << endl;
    // } 
    int delaunay_lattice = spg_delaunay_reduce(cell->lattice, SYMPREC);
    if(delaunay_lattice == 0) {
        cout << "Delaunay reduction failed" << endl;
    }
    cout << "The delaunay lattice is: " << endl;
    for(int i = 0; i < 3; i++) {
        cout << cell->lattice[i][0] << " " << cell->lattice[i][1] << " " << cell->lattice[i][2] << endl;
    } 
}

void get_ir_reciprocal_mesh(Cell *cell) {
    // get mesh
    double position[cell->atom_num][3];
    for(int i = 0; i < cell->position.size(); i++) {
        for(int j = 0; j < 3; j++) {
            position[i][j] = cell->position[i][j];
        }
    }

    int types[cell->atom_num];
    for(int i = 0; i < cell->types.size(); i++) {
        types[i] = cell->types[i];
    }

    int num_primitive_atom = spg_find_primitive(cell->lattice, position, types, cell->atom_num, SYMPREC);
    if(num_primitive_atom == 0) {
        throw Exception("No primitive atom");
    }

    double outer_result[3];
    cross(cell->lattice[1], cell->lattice[2], outer_result);

    double omega = cell->lattice[0][0] * outer_result[0] +
                   cell->lattice[0][1] * outer_result[1] +
                   cell->lattice[0][2] * outer_result[2];

    double b1[3], b2[3], b3[3];

    cross(cell->lattice[1], cell->lattice[2], b1);
    cross(cell->lattice[2], cell->lattice[0], b2);
    cross(cell->lattice[0], cell->lattice[1], b3);

    double coff = 2.0 * PI / omega;
    for(int i = 0; i < 3; i++) {
        b1[i] = b1[i] * coff;
        b2[i] = b2[i] * coff;
        b3[i] = b3[i] * coff;
    }
  
    double calc1 = 0.0, calc2 = 0.0, calc3 = 0.0;
    for(int i = 0; i < 3; i++) {
        calc1 += b1[i] * b1[i];
        calc2 += b2[i] * b2[i];
        calc3 += b3[i] * b3[i];
    }

    double rlc_1 = round(sqrt(calc1), 6);
    double rlc_2 = round(sqrt(calc2), 6);
    double rlc_3 = round(sqrt(calc3), 6);
    
    double kpresolv = 0.06;

    int mesh[3];
    mesh[0] = int(round(rlc_1 / (2 * PI * kpresolv)));
    mesh[1] = int(round(rlc_2 / (2 * PI * kpresolv)));
    mesh[2] = int(round(rlc_3 / (2 * PI * kpresolv)));

    int is_shift[] = {0, 0, 0};
    int cnt = mesh[0] * mesh[1] * mesh[2];

    int grid_address[cnt][3];
    int grid_mapping_table[cnt];
    

    int num_ir = spg_get_ir_reciprocal_mesh(grid_address,
                                            grid_mapping_table,
                                            mesh,
                                            is_shift,
                                            1,
                                            cell->lattice,
                                            position,
                                            types,
                                            cell->atom_num,
                                            SYMPREC);

    cout << "Number of ir-kpoints: " << num_ir << endl;

    cout << "[ " << endl;
    for(int i = 0; i < cnt; i++) {
        cout << "  [";
        for(int j = 0; j < 3; j++) {
            cout << setw(4) << left << grid_address[grid_mapping_table[i]][j] * 1.0 / mesh[j];
        }
        cout << "]" << endl;
    }
    cout << "]" << endl;
}

void show_cell(double lattice[3][3], double position[][3], const int types[], const int atom_num) {
    cout << "Lattice parameter:" << endl;
    for (int i = 0; i < 3; i++) {
        cout << lattice[i][0] << " " << lattice[i][1] << " " << lattice[i][2] << endl;
    }

    cout << "Atomic positions:" << endl;
    for (int i = 0; i < atom_num; i++) {
        cout << types[i] << ": " << position[i][0] << " " << position[i][1] << " " << position[i][2] << endl;
    }
}

void show_spg_dataset(double lattice[3][3], double position[][3], const int atom_num, const int types[]) {
    SpglibDataset *dataset;
    char ptsymbol[6];
    
    int pt_trans_mat[3][3];
    int i, j, size;
    const char *wl = "abcdefghijklmnopqrstuvwxyz";

    dataset = spg_get_dataset(lattice,
                              position,
                              types,
                              atom_num,
                              SYMPREC);

    cout << "International: " << dataset->international_symbol << " (" << dataset->spacegroup_number << ")" << endl;
    cout << "Hall symbol:\t" << dataset->hall_number << endl;

    spg_get_pointgroup(ptsymbol,
                       pt_trans_mat,
                       dataset->rotations,
                       dataset->n_operations);

    cout << "Point group:\t" << ptsymbol << endl;
    cout << "Transformation matrix:" << endl;
    for(i = 0; i < 3; i++) {
        cout << dataset->transformation_matrix[i][0] << " " << dataset->transformation_matrix[i][1] << " " << dataset->transformation_matrix[i][2] << endl;
    }

    cout << "Wyckoff letters:" << endl;

    for(i = 0; i < dataset->n_atoms; i++ ) {
        cout << wl[dataset->wyckoffs[i]] << " ";
    }
    cout << endl;

    cout << "Equivalent atoms:" << endl;
    for(i = 0; i < dataset->n_atoms; i++) {
        cout << dataset->equivalent_atoms[i] << " ";
    }
    cout << endl;

    for(i = 0; i < dataset->n_operations; i++) {
        cout << "---- " << i+1 << " ----" << endl;
        for (j = 0; j < 3; j++) {
            cout << dataset->rotations[i][j][0] << " " << dataset->rotations[i][j][1] << " " << dataset->rotations[i][j][2] << endl;
        }
        cout << dataset->translations[i][0] << " " << dataset->translations[i][1] << " " << dataset->translations[i][2] << endl;
    }

    spg_free_dataset(dataset);
}