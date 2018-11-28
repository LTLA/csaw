#ifndef INTERSECTOR_H
#define INTERSECTOR_H
#include "utils.h"

/* Designed to remove sorted reads that fall wholly within any entry of a given set of regions.  
 * We use the real alignment width, just in case we have very long reads in the alignment that 
 * are heavily soft-clipped (i.e., they should be reported as within but the read length will 
 * put them out).
 */

class intersector {
public:
    intersector(SEXP, SEXP);
    void advance_to_start(int);
    bool end_is_within(int) const;

private:
    Rcpp::IntegerVector positions, elements;
    int index;

    std::vector<int> open;
    int num_open;

    int lastpos;
};

#endif
