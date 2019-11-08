#ifndef EXCEPTION_H
#define EXCEPTION_H

#include <string>
using std::string;

class exception
{
public:
    exception(string err) {
        this->msg = err;
    }
    
    
private:
    string msg;
};


#endif