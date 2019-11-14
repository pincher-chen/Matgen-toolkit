#ifndef FUNC_H
#define FUNC_H

#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <regex>
#include <algorithm>

using std::cerr;
using std::cout;
using std::endl;
using std::map;
using std::string;
using std::vector;


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
void printVec(vector<T> &vec) {
    for(auto &item : vec) {
        cerr << item << " ";
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