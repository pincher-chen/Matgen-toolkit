#ifndef MOLECULE_H
#define MOLECULE_H

#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include <cfloat>
#include <regex>
#include <algorithm>
#include <vector>
#include <map>
#include <set>

#include "exception.h"
#include "func.h"
#include "conf.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::ios;
using std::map;
using std::ofstream;
using std::pair;
using std::set;
using std::string;
using std::to_string;
using std::vector;

class Molecule
{

public:

    Molecule(string _filename) {
        this->filename = _filename;

        maxX = maxY = maxZ = -1 * DBL_MAX;
        minX = minY = minZ = DBL_MAX;
    }

    inline double min(double a, double b) {
        if(a < b) {
            return a;
        }
        return b;
    }

    inline double max(double a, double b) {
        if(a < b) {
            return b;
        }
        return a;
    }

    void parse_file() noexcept(false) {
        ifstream in(this->filename, ios::in);
        if(!in.is_open()) {
            throw Exception(this->filename + " does not exist or fails to open!");
        }
        else {
            string str((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
            this->buffer = str;

            vector<string> lines = re_split(str, "\\n");            


            for(auto &line : lines) {
                vector<string> items = re_split(line, "\\s+");

                if(items.size() >= 4 && std::find(elements.begin(), elements.end(), items[3]) != elements.end()) {
                    this->ele.push_back(items[3]);

                    double x = get_num(items[0]);
                    double y = get_num(items[1]);
                    double z = get_num(items[2]);

                    minX = min(x, minX);
                    minY = min(y, minY);
                    minZ = min(z, minZ);

                    maxX = max(x, maxX);
                    maxY = max(y, maxY);
                    maxZ = max(z, maxZ);

                    this->atoms_cd.push_back({x, y, z});
                }
            }

            // // show
            // for(int i = 0; i < ele.size(); i++) {
            //     cout << this->ele[i] << "\t";
            //     printVec(this->atoms_cd[i]);
            //     cout << endl;
            // }
        }
    }

    vector<double> get_direct2vertex(int aims) {
        double x = atoms_cd[aims][0];
        double y = atoms_cd[aims][1];
        double z = atoms_cd[aims][2];

        double _x = x, _y = y, _z = z;
        
        // A(maxX, minY, minZ)
        double distance = pow(this->maxX - x, 2) + pow(this->minY - y, 2) + pow(this->minZ - z, 2);
        _x = x - this->maxX;;
        _y = y - this->minY;
        _z = z - this->minZ;

        // B(maxX, maxY, minZ)
        if(distance > pow(this->maxX - x, 2) + pow(this->maxY - y, 2) + pow(this->minZ - z, 2)) {
            distance = pow(this->maxX - x, 2) + pow(this->maxY - y, 2) + pow(this->minZ - z, 2);
            _x = x - this->maxX;;
            _y = y - this->maxY;
            _z = z - this->minZ;
        }

        // C(minX, minY, minZ)
        if(distance > pow(this->minX - x, 2) + pow(this->minY - y, 2) + pow(this->minZ - z, 2)) {
            distance = pow(this->minX - x, 2) + pow(this->minY - y, 2) + pow(this->minZ - z, 2);
            _x = x - this->minX;;
            _y = y - this->minY;
            _z = z - this->minZ;
        }

        // D(minX, maxY, minZ)
        if(distance > pow(this->minX - x, 2) + pow(this->maxY - y, 2) + pow(this->minZ - z, 2)) {
            distance = pow(this->minX - x, 2) + pow(this->maxY - y, 2) + pow(this->minZ - z, 2);
            _x = x - this->minX;;
            _y = y - this->maxY;
            _z = z - this->minZ;
        }

        // E(maxX, minY, maxZ)
        if(distance > pow(this->maxX - x, 2) + pow(this->minY - y, 2) + pow(this->maxZ - z, 2)) {
            distance = pow(this->maxX - x, 2) + pow(this->minY - y, 2) + pow(this->maxZ - z, 2);
            _x = x - this->maxX;;
            _y = y - this->minY;
            _z = z - this->maxZ;
        }

        // F(maxX, maxY, maxZ)
        if(distance > pow(this->maxX - x, 2) + pow(this->maxY - y, 2) + pow(this->maxZ - z, 2)) {
            distance = pow(this->maxX - x, 2) + pow(this->maxY - y, 2) + pow(this->maxZ - z, 2);
            _x = x - this->maxX;;
            _y = y - this->maxY;
            _z = z - this->maxZ;
        }

        // G(minX, minY, maxZ)
        if(distance > pow(this->minX - x, 2) + pow(this->minY - y, 2) + pow(this->maxZ - z, 2)) {
            distance = pow(this->minX - x, 2) + pow(this->minY - y, 2) + pow(this->maxZ - z, 2);
            _x = x - this->minX;
            _y = y - this->minY;
            _z = z - this->maxZ;
        }

        // H(minX, maxY, maxZ)  
        if(distance > pow(this->minX - x, 2) + pow(this->maxY - y, 2) + pow(this->maxZ - z, 2)) {
            distance = pow(this->minX - x, 2) + pow(this->maxY - y, 2) + pow(this->maxZ - z, 2);
            _x = x - this->minX;
            _y = y - this->maxY;
            _z = z - this->maxZ;
        }

        return {_x, _y, _z};
    }

    vector<double> get_direct(int aims) {
        double centerX = (this->minX + this->maxX) / 2;
        double centerY = (this->minY + this->maxY) / 2;
        double centerZ = (this->minZ + this->maxZ) / 2;

        double x = this->atoms_cd[aims][0] - centerX;
        double y = this->atoms_cd[aims][1] - centerY;
        double z = this->atoms_cd[aims][2] - centerZ;

        return {x, y, z};
    }

    vector<double> get_direct2flat(int aims) {
        double x = atoms_cd[aims][0];
        double y = atoms_cd[aims][1];
        double z = atoms_cd[aims][2];
        
        double _x = x, _y = y, _z = z;

        // before
        double distance = pow((this->maxX - x), 2);
        _x = 1;
        _y = 0;
        _z = 0;

        // after
        if(distance > pow(x - this->minX, 2)) {
            _x = -1;
            _y = 0;
            _z = 0;
            distance = pow(x - this->minX, 2);
        }

        // left
        if(distance > pow(y - this->minY, 2)) {
            _x = 0;
            _y = -1;
            _z = 0;
            distance = pow(y - this->maxY, 2);
        }

        // right
        if(distance > pow(this->maxY - y, 2)) {
            _x = 0;
            _y = 1;
            _z = 0;
            distance = pow(this->maxY - y, 2);
        }

        // up
        if(distance > pow(this->maxZ - z, 2)) {
            _x = 0;
            _y = 0;
            _z = 1;
            distance = pow(this->maxZ - z, 2);
        }

        // down
        if(distance > pow(z - this->minZ, 2)) {
            _x = 0;
            _y = 0;
            _z = -1;
            distance = pow(z - this->minZ, 2);
        }
        return {_x, _y, _z};
    }

public:
    vector<vector<double>> atoms_cd;
    vector<string> ele;

private:
    /* 
    box 
          G _____________ H
          /|            /|
         / |           / |
       E/_____________/F |
        |  |C________|___|D
        | /          |  / 
        |/___________|/   
        A             B
        
        A(maxX, minY, minZ)
        B(maxX, maxY, minZ)
        C(minX, minY, minZ) 
        D(minX, maxY, minZ)       
        E(maxX, minY, maxZ)
        F(maxX, maxY, maxZ)
        G(minX, minY, maxZ) 
        H(minX, maxY, maxZ)  
    */

    double minX, minY, minZ, maxX, maxY, maxZ;
    string filename;
    string buffer;
};

#endif