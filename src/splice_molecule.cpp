#include <iostream>
#include <cmath>
#include <iomanip>
#include "../include/cmdline.h"
#include "../include/func.h"
#include "../include/molecule.h"
#include "../include/cif_func.h"

using namespace std;

#define PI acos(-1)

const string conf_dir = "./conf/";

// 计算向量的模长的平方
double calc_vector_length(vector<double> &vec);

// 计算向量变换的矩阵
vector<vector<double>> calc_matrix(vector<double> &src, vector<double> &des);

// 计算向量叉乘
vector<double> cross_product(vector<double> &a, vector<double> &b);

// 判断两个向量是否平行
bool isParallel(vector<double> &a, vector<double> &b);

// 判断两个向量的是否共向
bool isCollinear(vector<double> &a, vector<double> &b);

// 判断两个浮点型是否相等
bool equal(double a, double b);

vector<vector<double>> calc_rotation_matrix(vector<double> vec, double theta);

vector<double> route_x(vector<double> &point, double theta);

vector<double> route_y(vector<double> &point, double theta);

vector<double> route_z(vector<double> &point, double theta);

int main(int argc, char *argv[]) {

    // cmd
    cmdline::parser parser;
    parser.add<string>("molecule_a", 'a', "path of the molecule A", true, "");
    parser.add<string>("molecule_b", 'b', "path of the molecule B", true, "");
    parser.add<string>("output", 'o', "the output path", true, "");
    parser.add<string>("type", 't', "the type of the result format(gjf/xyz), convert format to gaussion format or xyz format", true, "");

    parser.add<int>("connect_a", 'i', "the serial number of connect site in molecule A", true);
    parser.add<int>("connect_b", 'j', "the serial number of connect site in molecule B", true);

    parser.parse_check(argc, argv);
    
    if (argc == 1 || parser.exist("help")) {
        cout << parser.usage();
        return 0;
    }

    string inputA = parser.get<string>("molecule_a");
    string inputB = parser.get<string>("molecule_b");
    string output = parser.get<string>("output");
    string type = parser.get<string>("type");

    int indexA = parser.get<int>("connect_a") - 1;
    int indexB = parser.get<int>("connect_b") - 1;
    
    try {
        get_res();

        Molecule molA(inputA);
        molA.parse_file();

        Molecule molB(inputB);
        molB.parse_file();

        vector<double> vecA = molA.get_direct(indexA);
        
        /*
            计算A分子与connect_a的延长线上B分子的connect_b对应的坐标
        */
        vector<double> b_point;
        auto atom_a = radius_dict.find(molA.ele[indexA]);
        auto atom_b = radius_dict.find(molB.ele[indexB]);
        if(atom_a == radius_dict.end()) {
            atom_a = radius_dict.find(Universal);
            // throw Exception(atom_a + " 's radius not found");
        }

        if(atom_b == radius_dict.end()) {
            atom_b = radius_dict.find(Universal);
            // throw Exception(atom_b + " 's radius not found");
        }

        double radius_a = get_num(atom_a->second.radius);
        double radius_b = get_num(atom_b->second.radius);

        double length = radius_a + radius_b;
        
        /*
            for length = |vecA * k| 
            => k = length / |vecA|
        */
        double coef = length / calc_vector_length(vecA);
        for(int i = 0; i < vecA.size(); i++) {
            b_point.push_back(coef * vecA[i] + molA.atoms_cd[indexA][i]);
        }

        /*
            B 分子坐标旋转
            B分子中心点与connect_b的向量方向 平行 A分子中心与connect_a的向量方向的逆方向
        */
        // 得到vecA的反向量
        vector<double> oppVecA;
        for(auto &v : vecA) {
            if(!equal(v, 0)) {
                oppVecA.push_back(-1 * v);
            }
            else {
                oppVecA.push_back(v);
            }
        }

        vector<double> vecB = molB.get_direct(indexB);


        // 计算 A X B
        vector<double> product = cross_product(oppVecA, vecB);
        if(equal(0, calc_vector_length(product))) {
            auto a = route_z(vecB, PI / 2);
            auto b = route_y(vecB, PI / 2);

            if(!equal(0, calc_vector_length(a))) {
                product = {0, 0, 1};
            }
            else if(!equal(0, calc_vector_length(b))) {
                product = {0, 1, 0};
            }
            else {
                cerr << "[ERROR] A and B 'normal vector calculates failed" << endl;
            }
        }

        // 计算A B向量夹角
        double theta = 0, tmp = 0;
        for(int i = 0; i < vecB.size(); i++) {
            tmp += oppVecA[i] * vecB[i];
        }
        tmp = tmp / (calc_vector_length(oppVecA) * calc_vector_length(vecB));
        theta = acos(tmp);

        vector<vector<double>> T;
        {
            auto T1 = calc_rotation_matrix(product, theta);
            auto T2 = calc_rotation_matrix(product, 2 * PI - theta);
            auto T3 = calc_rotation_matrix(product, PI);

            vector<double> c1(3, 0.0);
            vector<double> c2(3, 0.0);
            vector<double> c3(3, 0.0);

            for(int i = 0; i < 3; i++) {
                for(int j = 0; j < 3; j++) {
                    c1[i] += T1[i][j] * vecB[j];
                    c2[i] += T2[i][j] * vecB[j];
                    c3[i] += T3[i][j] * vecB[j];
                }
            }

            if(isCollinear(c1, oppVecA)) {
                T = calc_rotation_matrix(product, theta);
                theta = theta;
            }
            else if(isCollinear(c2, oppVecA)) {
                T = calc_rotation_matrix(product, 2 * PI - theta);
                theta = 2 * PI - theta;
            }
            else {
                T = calc_rotation_matrix(product, PI);
                theta = PI;
            }
        }

        // cout << "Theta: " << theta << endl;
        // for(auto i : T) {
        //     printVec(i);
        //     cout << endl;
        // }

        vector<vector<double>> new_atoms_cd;
        for(auto item : molB.atoms_cd) {
            vector<double> p;
            for(int i = 0; i < T.size(); i++) {
                double sum = 0;
                for(int j = 0; j < T[i].size(); j++) {
                    sum += T[i][j] * item[j];
                }
                p.push_back(sum);
            }
            new_atoms_cd.push_back(p);
        }

        // 得到偏移矩阵
        vector<double> offset;

        for(int i = 0; i < b_point.size(); i++) {
            offset.push_back(b_point[i] - new_atoms_cd[indexB][i]);
        }

        // printVec(offset);
        // cout << endl;


        // 坐标加上偏移
        for(auto &item : new_atoms_cd) {
            for(int i = 0; i < offset.size(); i++) {
                item[i] += offset[i];
            }
        }

        double minX, minY, minZ, maxX, maxY, maxZ;
        maxX = maxY = maxZ = -1 * DBL_MAX;
        minX = minY = minZ = DBL_MAX;

        for(auto &item : new_atoms_cd) {
            double x = item[0];
            double y = item[1];
            double z = item[2];

            minX = min(x, minX);
            minY = min(y, minY);
            minZ = min(z, minZ);

            maxX = max(x, maxX);
            maxY = max(y, maxY);
            maxZ = max(z, maxZ);
        }

        vector<double> newVec;
        newVec.push_back(new_atoms_cd[indexB][0] - (minX + maxX) / 2);
        newVec.push_back(new_atoms_cd[indexB][1] - (minY + maxY) / 2);
        newVec.push_back(new_atoms_cd[indexB][2] - (minZ + maxZ) / 2);

        // printVec(newVec);
        // cout << endl;
        
        // printVec(oppVecA);
        // cout << endl;

        // double _tmp = (newVec[0] * oppVecA[0] + newVec[1] * oppVecA[1] + newVec[2] * oppVecA[2]) / (calc_vector_length(oppVecA) * calc_vector_length(newVec));
        // cout << acos(_tmp) << endl;

        // print
        string name;
        if(type == "xyz") {
            name = output + "/molecule" + del_split(get_filename(inputA), '-')[1] + "&" + del_split(get_filename(inputB), '-')[1] + ".xyz";
            ofstream out(name);

            out << molA.atoms_cd.size() + molB.atoms_cd.size() << endl;
            out << "(xyz format) generated by MAT-TOOL." << endl;
            // A atoms
            out.flags(ios::fixed);
            out.precision(8);
            for(int i = 0; i < molA.atoms_cd.size(); i++) {
                out << setw(8) << left << molA.ele[i] << setw(20) << right << molA.atoms_cd[i][0] << setw(15) << right << molA.atoms_cd[i][1] << setw(15) << right << molA.atoms_cd[i][2] << endl;
            }

            for(int i = 0; i < new_atoms_cd.size(); i++) {
                out << setw(8) << left << molB.ele[i] << setw(20) << right << new_atoms_cd[i][0] << setw(15) << right << new_atoms_cd[i][1] << setw(15) << right << new_atoms_cd[i][2] << endl;
            }

            out.close();
        }
        else if(type == "gjf") {
            string basic = conf_dir + "opt-freq.conf";
            ifstream in(basic, ios::in);

            if(!in) {
                cerr << "gaussion opt-freq.conf not found!" << endl;
            }
            string str((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
            in.close();

            name = output + "/molecule" + del_split(get_filename(inputA), '-')[1] + "&" + del_split(get_filename(inputB), '-')[1] + ".gjf";
            ofstream out(name);
            out << str << endl;

            out.flags(ios::fixed);
            out.precision(8);
            for(int i = 0; i < molA.atoms_cd.size(); i++) {
                out << setw(8) << left << " " + molA.ele[i] << setw(20) << right << molA.atoms_cd[i][0] << setw(15) << right << molA.atoms_cd[i][1] << setw(15) << right << molA.atoms_cd[i][2] << endl;
            }

            for(int i = 0; i < new_atoms_cd.size(); i++) {
                out << setw(8) << left << " " + molB.ele[i] << setw(20) << right << new_atoms_cd[i][0] << setw(15) << right << new_atoms_cd[i][1] << setw(15) << right << new_atoms_cd[i][2] << endl;
            }
        }

        cout << "Export file " << name << " successfully!" << endl;
    }
    catch(Exception &e) {
        cerr << e.msg << endl;
    }
    return 0;
}

// 计算向量的模长
double calc_vector_length(vector<double> &vec) {
    double sum = 0;
    for(auto num : vec) {
        sum += num * num;
    }
    return sqrt(sum);
}

// 计算向量变换的矩阵, 3 * 3
vector<vector<double>> calc_matrix(vector<double> &src, vector<double> &des) {
    /*
        for vector A and vector B
        B = T * A

        => T = (B * A^T) / (||A||^2)
    */

    vector<vector<double>> T;
   
    double length = calc_vector_length(src);

    for(int i = 0; i < 3; i++) {
        vector<double> t;
        for(int j = 0; j < 3; j++) {
            t.push_back(des[i] * src[j] / (length * length));
        }
        T.push_back(t);
    }

    // // test
    // cout << "des: ";
    // printVec(des);
    // cout << endl;
    // cout << "calc: ";
    // for(int i = 0; i < 3; i++) {
    //     double cur = 0;
    //     for(int j = 0; j < 3; j++) {
    //         cur += T[i][j] * src[j];
    //     }
    //     cout << cur << " ";
    // }
    // cout << endl;
}

vector<vector<double>> calc_rotation_matrix(vector<double> vec, double theta) {
    /*
        将向量的叉乘得到向量作为旋转轴

        旋转分解：
        1. 将整个坐标轴旋转，使得旋转轴 p 和 z 轴重合
        2. 再将点 P 绕 z 轴旋转 θ 角
        3. 再将整个坐标轴旋转回原位
    */

    // 将向量变为单位向量
    double mold = calc_vector_length(vec);
    for(auto &v : vec) {
        v /= mold;
    }

    double psi, phi;

    /* calculate psi */
    // for xoy
    {
        vector<double> t1 = {vec[0], vec[1], 0};
        vector<double> t2 = {1, 0, 0};
        double cosXY = 0;

        for (int i = 0; i < t1.size(); i++) {
            cosXY += t1[i] * t2[i];
        }

        if (equal(0, calc_vector_length(t1)) || equal(0, calc_vector_length(t2))) {
            psi = 0;
        }
        else {
            cosXY = cosXY / (calc_vector_length(t1) * calc_vector_length(t2));
            psi = acos(cosXY);

            auto a = route_z(t2, psi);
            auto b = route_z(t2, 2 * PI - psi);
            auto c = route_z(t2, 0);
            auto d = route_z(t2, PI);

            if (isCollinear(a, t1)) {
                psi = psi;
            }
            else if (isCollinear(b, t1)) {
                psi = 2 * PI - psi;
            }
            else if (isCollinear(c, t1)) {
                psi = 0;
            }
            else if (isCollinear(d, t1)) {
                psi = PI;
            }
            else {
                cout << "other" << endl;
            }
        }
    }

    /* calculate phi */
    {
        vector<double> t1 = {vec[0], vec[1], vec[2]};
        vector<double> t2 = {0, 0, 1};

        if (equal(0, calc_vector_length(t1)) || equal(0, calc_vector_length(t2))) {
            phi = 0;
        }
        else {
            double tmp = t1[2] / (calc_vector_length(t1) * calc_vector_length(t2));
            phi = acos(tmp);
        }
    }

    // cout << psi << " " << phi << " " << theta << endl;

    double x, y, z;
    x = sin(phi) * cos(psi);
    y = sin(phi) * sin(psi);
    z = cos(phi);

    vector<vector<double>> T;
    /*
        T = [
                [cos + x^2 (1 - cos)        xy(1-cos)-zsin          xz(1 - cos) + ysin]
                [yx(1 - cos) + zsin         cos + y^2 (1 - cos)     yz(1 - cos) - xsin]
                [zx(1 - cos) - ysin         zy(1 - cos) + xsin      cos + z^2 (1 - cos)]
            ]
    */

    T.push_back({cos(theta) + pow(x, 2) * (1 - cos(theta)),
                 x * y * (1 - cos(theta)) - z * sin(theta),
                 x * z * (1 - cos(theta)) + y * sin(theta)});

    T.push_back({y * x * (1 - cos(theta)) + z * sin(theta),
                 cos(theta) + pow(y, 2) * (1 - cos(theta)),
                 y * z * (1 - cos(theta)) - x * sin(theta)});
    
    T.push_back({z * x * (1 - cos(theta)) - y * sin(theta),
                 z * y * (1 - cos(theta)) + x * sin(theta),
                 cos(theta) + pow(z, 2) * (1 - cos(theta))});
    return T;
}

vector<double> cross_product(vector<double> &a, vector<double> &b) {
    // a = (a1, a2, a3), b = (b1, b2, b3)
    // a x b = (a2b3-a3b2, a3b1-a1b3, a1b2-a2b1)

    vector<double> res(3, 0);
    res[0] = a[1] * b[2] - a[2] * b[1];
    res[1] = a[2] * b[0] - a[0] * b[2];
    res[2] = a[0] * b[1] - a[1] * b[0];

    return res;
}

vector<double> route_x(vector<double> &point, double theta) {
    /*
        [
            [1  0   0 ]
            [0  cos -sin]
            [0  sin cos]
        ]
    */
    vector<vector<double>> T;
    T.push_back({1, 0, 0});
    T.push_back({0, cos(theta), -1 * sin(theta)});
    T.push_back({0, sin(theta), cos(theta)});


    vector<double> res(3, 0);
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            res[i] += T[i][j] * point[j];
        }
    }

    return res;
}

vector<double> route_y(vector<double> &point, double theta) {
    /*
        [
            [cos    0   sin]
            [0      1   0]
            [-sin   0   cos]
        ]
    */
    vector<vector<double>> T;
    T.push_back({cos(theta), 0, sin(theta)});
    T.push_back({0, 1, 0});
    T.push_back({-1 * sin(theta), 0, cos(theta)});


    vector<double> res(3, 0);
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            res[i] += T[i][j] * point[j];
        }
    }

    return res;
}

vector<double> route_z(vector<double> &point, double theta) {
    /*
        [
            [cos    -sin    0]
            [sin    cos     0]
            [0      0       1]
        ]
    */
    vector<vector<double>> T;
    T.push_back({cos(theta), -1 * sin(theta), 0});
    T.push_back({sin(theta), cos(theta), 0});
    T.push_back({0, 0, 1});


    vector<double> res(3, 0);
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            res[i] += T[i][j] * point[j];
        }
    }

    return res;
}

bool isParallel(vector<double> &a, vector<double> &b) {
    double k;
    bool flag = false;

    // printVec(a);
    // cout << endl;
    // printVec(b);
    // cout << endl;

    for(int i = 0; i < a.size(); i++) {
        if(equal(a[i], 0.0) || equal(b[i], 0.0)) {
            if(equal(a[i], 0) && equal(b[i], 0)) {
                continue;
            }
            return false;
        }
        else {
            if(!flag) {
                // cout << i << " " << a[i] << " " << b[i] << endl;
                k = a[i] / b[i];
                flag = true;
            }
            else if(!equal(k, a[i] / b[i])) {
                return false;
            }
        }
    }
    return true;
}

bool equal(double a, double b) {
    if((a - b) > -1e-6 && (a - b) < 1e-6) {
        return true;
    }

    return false;
}

bool isCollinear(vector<double> &a, vector<double> &b) {
    double k;
    bool flag = false;

    // printVec(a);
    // cout << endl;
    // printVec(b);
    // cout << endl;

    for(int i = 0; i < a.size(); i++) {
        if(equal(a[i], 0.0) || equal(b[i], 0.0)) {
            if(equal(a[i], 0) && equal(b[i], 0)) {
                continue;
            }
            return false;
        }
        else {
            if(!flag) {
                // cout << i << " " << a[i] << " " << b[i] << endl;
                k = a[i] / b[i];
                flag = true;
            }
            else if(!equal(k, a[i] / b[i])) {
                return false;
            }
        }
    }
    return k > 0;
}