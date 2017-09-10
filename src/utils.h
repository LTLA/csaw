#include "Rcpp.h"
#ifndef UTILS_H
#define UTILS_H

template<typename T>
Rcpp::StringVector makeStringVector(T start, T end) {
    Rcpp::StringVector output(end-start);
    Rcpp::StringVector::iterator oIt=output.begin();
    while (start!=end) {
        (*oIt)=*start;
        ++start;
        ++oIt;
    }
    return output;
}   

#endif
