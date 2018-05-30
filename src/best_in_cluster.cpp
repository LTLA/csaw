#include "csaw.h"
#include "utils.h"

SEXP best_in_cluster(SEXP pval, SEXP by, SEXP weight) {
    BEGIN_RCPP
    const Rcpp::NumericVector winp(pval);
    const Rcpp::IntegerVector clustids(by);
    const Rcpp::NumericVector winweight(weight);

    const size_t nwin=winp.size();
	if (nwin!=clustids.size() || nwin!=winweight.size()) {
        throw std::runtime_error("input vector lengths are not equal"); 
    }
    size_t nclust=check_By_vector(clustids.begin(), clustids.end());

	// Pulling out results.
    Rcpp::NumericVector out_pval(nclust);
    Rcpp::IntegerVector out_best(nclust);
    auto opIt=out_pval.begin();
    auto obIt=out_best.begin();

    size_t run_start=0;
    while (run_start < nwin) {
        size_t run_end=run_start+1;
        double subweight=winweight[run_start];
        while (run_end < nwin && clustids[run_start]==clustids[run_end]) { 
            subweight+=winweight[run_end];
            ++run_end; 
        }

		/* Computing the Holm p-value for the best window (basically Bonferroni, if we're taking the minimum).
		 * Weights are defined according to the weighted Bonferroni (see http://arxiv.org/abs/math.ST/0604172,
		 * though some mental arithmetic is needed). These can also be treated as relative frequency weights,
		 * i.e. the total number of tests is rescaled relative to the weight of the current test (so, [10,1] 
		 * weights would consider there to be 1.1 tests for the first one and 11 tests for the second one).
		 */
        size_t best=run_start;
		double& outp=(*opIt=winp[run_start]/winweight[run_start]);
		for (size_t curwin=run_start+1; curwin<run_end; ++curwin) {
			const double tempp=winp[curwin]/winweight[curwin];
			if (tempp < outp) { 
				best=curwin;
				outp=tempp;
			}
		}

        outp*=subweight;
	    if (outp > 1) { outp=1; }	
	    *obIt=int(best+1); // For 1-based indexing in R.

		// Setting it up for the next round.
        ++obIt;
        ++opIt;
		run_start=run_end;
	}
	
    return Rcpp::List::create(out_pval, out_best);
    END_RCPP
}
