#include "csaw.h"
#include "utils.h"

struct region_data {
    region_data(int i, int e, double m) : index(i), endpt(e), metric(m) {}
    int index, endpt;
    double metric;

    friend bool operator<(const region_data& left, const region_data& right) { 
        if (left.metric==right.metric) { 
            if (left.endpt==right.endpt) {
                return left.index < right.index; 
            } 
            return left.endpt < right.endpt;
        }
        return left.metric < right.metric;  
    }
};

typedef std::multiset<region_data> order_set;
struct compare_iterators {
	bool operator() (const order_set::iterator& left, const order_set::iterator& right) const {
		if (left->endpt==right->endpt)  { return (left->index > right->index); } // Opposite, as we want the smallest value to be treated as the largest in the queue.
		return (left->endpt > right->endpt); 
	}
};
typedef std::priority_queue<order_set::iterator, std::deque<order_set::iterator>, compare_iterators> pqueue;

/* This function scans through the track and pulls out local maxima. */

SEXP find_maxima(SEXP chrs, SEXP starts, SEXP ends, SEXP metric, SEXP range) {
    BEGIN_RCPP

    const Rcpp::IntegerVector _chrs(chrs), _starts(starts), _ends(ends);
    const Rcpp::NumericVector _metric(metric);
	const int nlen=_chrs.size();
    if (nlen!=_starts.size() || nlen!=_ends.size() || nlen!=_metric.size()) {
	     throw std::runtime_error("vectors must be of equal length"); 
    }

    const int maxrange=check_integer_scalar(range, "range");
	if (maxrange <= 0) { 
        throw std::runtime_error("range should be a positive integer"); 
    }
    
    // Setting up structures to compute maxima on the fly.
	order_set incoming;
	pqueue first_to_leave;
    if (nlen) { 
	    first_to_leave.push(incoming.insert(region_data(0, _ends[0], _metric[0])));
    }
    Rcpp::LogicalVector output(nlen);

	// Assuming we're sorted by sptr.
	int right_edge=1;
	double cur_max=_metric[0];
	for (int i=0; i<nlen; ++i) {
		if (i) {
			if (_chrs[i] < _chrs[i-1]) { throw std::runtime_error("regions must be sorted by chromosome"); }
			else if (_chrs[i]==_chrs[i-1]) { 
				if (_starts[i] < _starts[i-1]) { throw std::runtime_error("regions on the same chromosome must be sorted by start position"); }
				if (_ends[i] < _ends[i-1]) { throw std::runtime_error("nested regions are not supported"); }
			}
		}

		// Getting rid of regions running off the end.
		while (!first_to_leave.empty() && (_chrs[i] != _chrs[first_to_leave.top()->index] || _starts[i] - first_to_leave.top()->endpt > maxrange)) {
			incoming.erase(first_to_leave.top());
			first_to_leave.pop();
		}

		// Finding the maximum (records first tie, if there are ties).
		double max_right=R_NaReal;
        int is_max=0;
		while (right_edge < nlen && _chrs[i]==_chrs[right_edge] && _starts[right_edge] - _ends[i] <= maxrange) {
			if (ISNA(max_right) || max_right < _metric[right_edge]) {
				is_max=right_edge;
				max_right=_metric[right_edge];
			}
			++right_edge;
		}

		if (!ISNA(max_right)) {
			/* Checking if the new maximum is greater than the current maximum, in which 
			 * case we clear out everything because the new maximum supercedes anything 
			 * already in the set.
			 */	
			if (cur_max < max_right) {
				incoming.clear();
				first_to_leave = pqueue(compare_iterators());
			}
		
			/* Only adding those after the maximum, and which are also reverse
	 		 * cumulative maxima themselves (ties allowed). There's no point having a 
			 * window which is sandwiched by two maxima, as it'll never be its own maxima.
			 */
			int right_copy=right_edge-1;
			max_right=_metric[right_copy];
			while (right_copy >= is_max) { 
                const double& cur_metric=_metric[right_copy];
				if (cur_metric >= max_right) { 
					first_to_leave.push(incoming.insert(region_data(right_copy, _ends[right_copy], cur_metric)));
					max_right=cur_metric;
				}
				--right_copy;
			}
		}

//		Rprintf("Current: %i (%i) %.3f\n", i, right_edge, mptr[i]);
//		for (order_set::const_iterator itx=incoming.begin(); itx!=incoming.end(); ++itx) { 
//			Rprintf("\t%i %.3f\n", *itx, mptr[*itx]);
//		}

		// Checking if we're currently the max (some allowance for exactly tied values).
		if (incoming.size()) { 
            auto it=incoming.end();
			--it;
			cur_max=it->metric;
			do {
				if (it->index==i) { 
					output[i]=1;
					break;
				}
				if (it==incoming.begin()) { break; } 
				--it;
			} while (it->metric==cur_max);
		} else {
			throw std::runtime_error("empty set during maxima detection");
		}
	}
	
    return output;
    END_RCPP
}

