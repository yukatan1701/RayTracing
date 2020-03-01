#include "ParseException.h"

void ParseException::printMessage() const {
    if (type == ExceptionType::NOT_ENOUGH) {
        std::cerr << "Not enough arguments." << std::endl;
    } else if (type == ExceptionType::SAME_ARGS) {
        std::cerr << "Too many values for the same argument." << std::endl;
    } else {
        std::cerr << "Parse error." << std::endl;
    }
}