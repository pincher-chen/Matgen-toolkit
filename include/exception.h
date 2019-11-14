#ifndef EXCEPTION_H
#define EXCEPTION_H

#include <string>
using std::string;

class Exception
{
public:
    Exception(string err) : msg(err) {}
    string msg;
};

#endif