#include "csaw.h"
#include "utils.h"

/* This function scans through the track and pulls out local maxima. */

SEXP get_profile(SEXP starts, SEXP ends, SEXP regstarts, SEXP range, SEXP average, SEXP normalize) {
    BEGIN_RCPP

    const Rcpp::IntegerVector _starts(starts), _ends(ends), _regstarts(regstarts);
	const int nfrags=_starts.size();
    if (nfrags!=_ends.size()) {
        throw std::runtime_error("fragment start/end vectors should have same length"); 
    }
    const int nregs=_regstarts.size();
    if (nregs==0) { 
        throw std::runtime_error("no regions supplied"); 
    }
	
    // Setting up scalars.
    const int norm_type=check_integer_scalar(normalize, "normalize specification");
    const int maxrange=check_integer_scalar(range, "range distance");
    const bool use_average=check_logical_scalar(average, "average specification");

	/* Setting up a separate profile for each region. This is necessary
	 * to ensure that the calculations are integer (despite weighting),
	 * in order to preserve numerical stability.
	 */
	const int totallen=2*maxrange+1;
    Rcpp::IntegerMatrix profiles(totallen, nregs);
    std::vector<Rcpp::IntegerMatrix::iterator> all_profiles(nregs);
    if (nregs) { 
        all_profiles[0]=profiles.begin() + maxrange; // so index of 0 = distance of 0.
        for (int r=1; r<nregs; ++r) {
            all_profiles[r]=all_profiles[r-1]+totallen;
        }
    }

	/* Running through the reads. We use a strategy of identifying the
	 * regions for each read and adding that to the profile, rather than
	 * setting up some queue (which would involve lots of insertions/deletions).
	 */
	for (int frag=0; frag<nfrags; ++frag) {
		const int& curstart=_starts[frag];
		const int& curend=_ends[frag];

		/* Getting all regions starting after the fragment. Don't bother
		 * trying to optimize with reducing the binary search by sorting
		 * the fragments beforehand; the earlier sort would take more time
		 * anyway, because it's nfrag*log(nfrag), not nfrag*log(nregs).
		 */
        auto ptr=std::upper_bound(_regstarts.begin(), _regstarts.end(), curend), ptr_copy=ptr;
        int loc=ptr-_regstarts.begin();
		while (ptr!=_regstarts.end()) {
			int dist = *ptr - curend;
			if (dist > maxrange) { 
                break; 
            }
            auto curprof=all_profiles[loc];

			--(*(curprof-dist+1));
			int dist2 = *ptr - curstart;
			++(*(curprof-std::min(dist2, maxrange)));

			++ptr;
            ++loc;
		}

		// Getting all regions starting before the fragment end.
		ptr=ptr_copy;
        loc=ptr-_regstarts.begin();
		while (ptr!=_regstarts.begin()) { 
			--ptr;
            --loc;

			int dist = curstart - *ptr;
			if (dist > maxrange) { 
                break; 
            }
            auto curprof=all_profiles[loc];

			++(*(curprof+std::max(dist, -maxrange)));
			int dist2 = curend - *ptr;
			if (dist2 < maxrange) { 
                --(*(curprof+dist2+1)); 
            }
		}
	}

	// Compiling profile based on addition/subtraction instructions.
	for (int r=0; r<nregs; ++r) {
        auto curprof=all_profiles[r];
		for (int i=-maxrange+1; i<=maxrange; ++i) {
		    *(curprof+i)+=*(curprof+i-1);
		}
	}

	if (!use_average){
        return profiles;
    }

    std::deque<double> weightings(nregs);
    if (norm_type==1) { 
        // No normalization.
        std::fill(weightings.begin(), weightings.end(), 1);
    } else if (norm_type==2) { 
        // Normalizing by the total coverage.
        for (int r=0; r<nregs; ++r) {
            auto curprof=all_profiles[r];
            weightings[r]=std::accumulate(curprof-maxrange, curprof+maxrange+1, 0.0);
        }
    } else if (norm_type==3) { 
        // Normalizing by the maximum height.
        for (int r=0; r<nregs; ++r) {
            auto curprof=all_profiles[r];
            weightings[r]=*std::max_element(curprof-maxrange, curprof+maxrange+1);
        }
    }

    Rcpp::NumericVector ave_out(totallen);
    auto aoIt=ave_out.begin()+maxrange; // 0 is now distance of zero.
	for (int r=0; r<nregs; ++r) {
        auto curprof=all_profiles[r];
		const double& cur_weighting=weightings[r];
		for (int i=-maxrange; i<=maxrange; ++i) { 
            *(aoIt+i)+=double(*(curprof+i))/cur_weighting; 
        }
	}
    return ave_out;
    END_RCPP
}

