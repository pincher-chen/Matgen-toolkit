#ifndef IO_HELPER_H
#define IO_HELPER_H

#include <string>

using namespace std;

class CIF {

public:
    CIF(string filePath) {
        this->filePath = filePath;
    }

    void parseCIF() {
        
    }

private:
    string filePath;
};


#endif