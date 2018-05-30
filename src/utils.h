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

template<typename T>
size_t check_By_vector(T start, T end) {
    if (start==end) {
        return 0;
    }

    size_t total=1;
    T next=start; 
    ++next;

    while (next!=end) { 
        if (*next < *start) {
            throw std::runtime_error("vector of cluster ids should be sorted");
        } else if (*next != *start) {
            ++total;
        }
        ++next;
        ++start;
    }
    return total;
}

bool check_logical_scalar(Rcpp::RObject x, const char* thing);

int check_integer_scalar(Rcpp::RObject x, const char* thing);

double check_numeric_scalar(Rcpp::RObject x, const char* thing);

Rcpp::String check_string(Rcpp::RObject x, const char* thing);

#endif
