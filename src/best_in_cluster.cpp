#include "csaw.h"

SEXP best_in_cluster(SEXP pval, SEXP by, SEXP weight) {
    BEGIN_RCPP
    Rcpp::NumericVector _pval(pval);
    Rcpp::IntegerVector _by(by);
    Rcpp::NumericVector _weight(weight);

	const int n=LENGTH(pval);
	if (n!=_by.size() || n!=_weight.size()) {
        throw std::runtime_error("input vector lengths are not equal"); 
    }

	// Checking that the 'by' is sorted, counting the number of elements.
	int total=0;
	if (n > 0) {
		total=1;
        Rcpp::IntegerVector::iterator bIt=_by.begin(), bIt_next=bIt+1;
	    while (bIt_next!=_by.end()) {
			if (*bIt_next < *bIt) {
                throw std::runtime_error("vector of cluster ids should be sorted"); 
            } else if (*bIt_next != *bIt) {
                ++total; 
            }
            ++bIt_next;
            ++bIt;
		}
	}

	// Pulling out results.
    Rcpp::NumericVector out_pval(total);
    Rcpp::IntegerVector out_best(total);
    Rcpp::NumericVector::iterator opIt=out_pval.begin();
    Rcpp::IntegerVector::iterator obIt=out_best.begin();

    int i=0;
    while (i<n) {
        int j=i+1;
        double subweight=_weight[i];
        while (j < n && _by[i]==_by[j]) { 
            subweight+=_weight[j];
            ++j; 
        }

		/* Computing the Holm p-value for the best window (basically Bonferroni, if we're taking the minimum).
		 * Weights are defined according to the weighted Bonferroni (see http://arxiv.org/abs/math.ST/0604172,
		 * though some mental arithmetic is needed). These can also be treated as relative frequency weights,
		 * i.e. the total number of tests is rescaled relative to the weight of the current test (so, [10,1] 
		 * weights would consider there to be 1.1 tests for the first one and 11 tests for the second one).
		 */
		int& outi=(*obIt=i);
		double& outp=(*opIt=_pval[i]/_weight[i]);
		double tempp=0;
		for (int x=i+1; x<j; ++x) {
			tempp=_pval[x]/_weight[x];
			if (tempp < outp) { 
				outi=x;
				outp=tempp;
			}
		}
		outp*=subweight;
	    if (outp > 1) { outp=1; }	
		++outi; // For 1-based indexing in R.

		// Setting it up for the next round.
        ++obIt;
        ++opIt;
		i=j;
	}
	
    return Rcpp::List::create(out_pval, out_best);
    END_RCPP
}
