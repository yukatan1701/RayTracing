#ifndef __PARSE_EXCEPTION__
#define __PARSE_EXCEPTION__

#include <iostream>

class ParseException : public std::exception {
public:
    enum ExceptionType { NOT_ENOUGH, SAME_ARGS };
    ParseException(ExceptionType type) : type(type) {}
    void printMessage() const;
private:
    ExceptionType type;
};

#endif