#ifndef FUNC_H
#define FUNC_H

#include <string>
#include <iostream>
#include <regex>
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <io.h>
#include <direct.h>

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::map;
using std::vector;
using std::set;

bool is_folder_exist(string folder) {
    return _access(folder.c_str(), 0) == 0;
}

bool make_dir(string path) {
    string cur = "";
    for(auto c : path) {
        cur += c;
        if(c == '\\' || c == '/') {
            if(!is_folder_exist(cur)) {
                int ret = _mkdir(cur.c_str());
                if(ret == -1) {
                    return false;
                }
            }
        }
    }

    if(!is_folder_exist(cur)) {
        int ret = _mkdir(cur.c_str());
        return ret != -1;
    }

    return true;
}

string get_filename(string path) {
    string::size_type pos = path.find_last_of('\\');
    
    if(pos == string::npos) {
        pos = path.find_last_not_of('/');
        if(pos == string::npos) {
            return path;
        }
    }
    return path.substr(pos+1);
}

// trim string
string trim(const string &str) {
    string::size_type pos = str.find_first_not_of(" \r\n\t");
    if (pos == string::npos) {
        return str;
    }
    string::size_type pos2 = str.find_last_not_of(" \r\n\t");
    if (pos2 != string::npos) {
        return str.substr(pos, pos2 - pos + 1);
    }
    return str.substr(pos);
}

// divide strings according to delimiters
vector<string> del_split(const string &str, char delimiter) {
    vector<string> result;
    string cur = str;
    while (!cur.empty()) {
        int index = cur.find_first_of(delimiter);
        if (index == -1) {
            result.push_back(cur);
            cur.clear();
        }
        else {
            if(trim(cur.substr(0, index)) != "") {
                result.push_back(cur.substr(0, index));
            }
            cur = cur.substr(index + 1, cur.size() - index - 1);
        }
    }
    return result;
}

// divide strings according to regular matches
vector<string> re_split(const string &str, string pattern) {
    string cur = str;
    std::regex rgx(pattern);
    vector<string> result;
    for (std::sregex_token_iterator iter(cur.begin(), cur.end(), rgx, -1), end; iter != end; ++iter) {
        // cerr << iter->str() << endl;
        string s = iter->str();
        if(!s.empty()) {
            result.push_back(iter->str());
        }
    }
    return result;
}

// get real num from string
double get_num(string value) {
    try {
        return atof(value.c_str());
    }
    catch(...) {
        vector<string> nums = del_split(value, '(');

        if(!nums.empty()) {
            return atof(nums[0].c_str());
        }

        return -1;
    }
}

// print vector for debugging
template <typename T> 
void printVec(vector<set<T>> &vec) {
    for(auto item : vec) {
        for(auto it = item.begin(); it != item.end(); it++) {
            cerr << *it << " ";
        }
        cerr << endl << endl;
    }
}

// print vector for debugging
template <typename T>
void printVec(vector<T> &vec) {
    for(auto &item : vec) {
        cerr << item << " ";
    }
    cerr << endl;
}


// print set for debugging
template <typename T>
void printSet(set<T> &set) {
    for(auto iter = set.begin(); iter != set.end(); iter++) {
        cerr << *iter << " ";
    }
    cerr << endl;
}

// print map for debugging
template <typename key, typename member>
void printMap(map<key, vector<member>> &m) {
    for(auto iter = m.begin(); iter != m.end(); iter++) {
        cerr << iter->first << " : ";
        for(auto it = iter->second.begin(); it != iter->second.end(); it++) {
            cerr << *it << " ";
        }
        cerr << endl;
    }
}

// print map for debugging
template <typename key, typename value>
void printMap(map<key, value> &m) {
    for(auto iter = m.begin(); iter != m.end(); iter++) {
        cerr << iter->first << " : " << iter->second << endl;
    }
}

#endif