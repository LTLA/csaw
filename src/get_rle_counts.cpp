#include "csaw.h"
#include "utils.h"

SEXP get_rle_counts(SEXP start, SEXP end, SEXP nr, SEXP space, SEXP first) {
    BEGIN_RCPP

    const int nrows=check_integer_scalar(nr, "number of rows");
    const int spacing=check_integer_scalar(space, "spacing");
    const int usefirst=check_logical_scalar(first, "first point specification"); // yes, "int" type with "logical" check is deliberate!
	
    const Rcpp::IntegerVector _start(start), _end(end);
    const int n=_start.size();
    if (n!=_end.size()) { 
        throw std::runtime_error("start/end vectors must have equal length"); 
    }

    // Running through output.
    Rcpp::IntegerVector output(nrows);
	for (int i=0; i<n; ++i) {
        const int& curstart=_start[i];
        const int& curend=_end[i];

		// Get the zero-index corresponding to the smallest spacing point larger than the current inclusive start/end.
        if (curend < curstart) { 
            throw std::runtime_error("invalid coordinates for read start/ends"); 
        }
        const int left=(curstart < 2 ? 0 : int((curstart-2)/spacing)+usefirst);
        const int right=(curend < 1 ? 0 : int((curend-1)/spacing)+usefirst);

        // Adding the steps for addition (at the start) and deletion (after the end) of each read.
        if (left<right) { 
            if (left<nrows) { ++output[left]; }
            if (right<nrows) { --output[right]; }
        }
    }

    // Running and computing the RLE, given the steps at each position.
    int cum=0;
    for (auto& o : output) { 
        cum+=o;
        o=cum;
    }
	
    return output;
    END_RCPP
}
