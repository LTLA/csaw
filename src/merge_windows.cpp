#include "csaw.h"
#include "utils.h"

int split_cluster(const Rcpp::IntegerVector&, const Rcpp::IntegerVector&, Rcpp::IntegerVector&, const int&, const int&, const int&, const int&);

/* We assume that incoming elements are sorted by chr -> start -> end. We then proceed 
 * to aggregate elements by chaining together elements that are less than 'tolerance'
 * apart and, if required, have the same 'sign'. We also split them if the difference
 * between the first start and the last end is greater than 'max_size'.
 */

SEXP merge_windows(SEXP chrs, SEXP start, SEXP end, SEXP sign, SEXP tolerance, SEXP max_size) {
    BEGIN_RCPP

    Rcpp::IntegerVector _chrs(chrs), _start(start), _end(end);
    Rcpp::LogicalVector _sign(sign);
	const int n = _chrs.size();
	if (n!=_start.size() || n!=_end.size() || n!=_sign.size()) { 
        throw std::runtime_error("lengths of vectors are not equal"); 
    }

    const int tol=check_integer_scalar(tolerance, "tolerance");
    Rcpp::IntegerVector _max_size(max_size);
	if (_max_size.size() > 1) { 
        throw std::runtime_error("maximum size should be an integer scalar or NULL"); 
    }
	const bool limit_size=(_max_size.size()==1);
	const int maxs=(limit_size ? _max_size[0] : 0);
	
	// Proceeding with the merge operation.
    Rcpp::IntegerVector out_index(n);
	int start_index=0;
    int ngroups=0;
    int last_end, last_sign;
    if (n) { 
        ngroups=1;
        last_end=_end[0];
        last_sign=_sign[0];
        out_index[0]=1;
    }
    bool warned=false;

	for (int i=1; i<n; ++i) {
        bool diffchr=(_chrs[i]!=_chrs[i-1]);
        bool diffsign=(_sign[i]!=last_sign);
        if (!diffchr && _start[i] < _start[i-1]) { 
            throw std::runtime_error("regions should be sorted by start position");
        }

		/* Fully nested regions don't have a properly defined interpretation when it comes
		 * to splitting things by sign. We only support consideration of nested regions where
		 * either of the boundaries are the same. That can be considered to break the 
		 * previous stretch if it had opposite sign. Otherwise, the next window would have to
		 * make the decision to match the parent window or its nested child.
		 *
		 * If the start is the same, the window with the earlier end point should be ordered
		 * first, so that won't pass the first 'if'. This means that it'll only enter with
		 * a fully nested window. Start and end-point equality might be possible at the ends 
		 * of chromosomes where trimming enforces sameness, but full nesting should not be observed.
		 *
		 * If the nested region has the same sign as the parent, then everything proceeds
		 * normally i.e. same cluster. We make sure to keep 'last_end' as the parent end in
		 * such cases. This ensures that the following windows get a change to compute
		 * distances to the parent end (which should be closer).
		 */
        bool nested=(!diffchr && _end[i] < last_end);
	    if (nested) { 
		   	if (diffsign && !warned) { 
                Rcpp::warning("fully nested windows of opposite sign are present and ignored"); 
                warned=true;
                diffsign=false;
            }
		} 

		if (diffchr 											// Obviously, changing if we're on a different chromosome.
			|| _start[i]-last_end-1 > tol						// Space between windows, start anew if this is greater than the tolerance.
			|| diffsign 										// Checking if the sign is consistent.
	   	) {
		    if (limit_size) { 
                // Splitting the cluster, if desired.
                ngroups=split_cluster(_start, _end, out_index, last_end, start_index, i, maxs); 
            } 
			++ngroups;
			out_index[i]=ngroups; 
			start_index=i;
		} else {
			out_index[i]=out_index[i-1];
		}

        // Using the parent window as the endpoint if it's nested, but otherwise bumping up the last stats.
        if (!nested) { 
            last_end=_end[i]; 
            last_sign=_sign[i];
        }
	}

	// Cleaning up the last cluster, if necessary.
  	if (limit_size) { 
        ngroups=split_cluster(_start, _end, out_index, last_end, start_index, n, maxs); 
    }

	// Now, identifying the chromosome, start and end of each region.
    Rcpp::IntegerVector out_chr(ngroups, -1), out_start(ngroups, -1), out_end(ngroups, -1);

	for (int i=0; i<n; ++i) { 
		int curgroup=out_index[i]-1;
		if (out_chr[curgroup]<0) { 
			out_chr[curgroup]=_chrs[i];
			out_start[curgroup]=_start[i]; // Sorted by start, remember; only need this once.
			out_end[curgroup]=_end[i];
		} else if (out_end[curgroup] < _end[i]) { 
		    out_end[curgroup]=_end[i]; 
		}
	}
	
    return Rcpp::List::create(out_index, out_chr, out_start, out_end);
    END_RCPP
}

int split_cluster(const Rcpp::IntegerVector& starts, const Rcpp::IntegerVector& ends, Rcpp::IntegerVector& output, 
        const int& actual_end, const int& xs, const int& xe, const int& width) {

	double full_width=actual_end-starts[xs]+1;
	if (full_width <= width) { return output[xs]; }
	int mult=int(std::ceil(full_width/width));

	/* There can only be `mult` subclusters. At the worst, `cur_diff`
	   will be equal to `actual_end - starts[xs]`. Division by `subwidth` will
	   give an expression of `(actual_end - starts[xs]) / [ (actual_end - starts[xs] + 1) / mult ]`.
	   This will always be less than `mult`, so flooring it will give `mult-1`,
	   i.e., the last index of `instantiated`.
	 */
    std::vector<int> instantiated(mult, 0);
	int output_index=output[xs];

	// Allocating windows into subclusters, based on their midpoints.
	double subwidth=full_width/mult, cur_diff;
	for (int i=xs; i<xe; ++i) {
		cur_diff = double(starts[i]+ends[i])*0.5 - starts[xs];
		output[i] = int(cur_diff/subwidth);	
		if (!instantiated[output[i]]) { instantiated[output[i]] = 1; }
	}

	/* Allocating output indices to the subclusters. This approach avoids
	   situations where you get nonconsecutive cluster indices, e.g., when
	   `tol` is greater than the maximum width.	 
	 */
	for (int i=0; i<mult; ++i) { 
		if (!instantiated[i]) { continue; }
		instantiated[i]=output_index;
		++output_index;
	}

	// Assigning indices back to the output vector.
	for (int i=xs; i<xe; ++i) { output[i]=instantiated[output[i]]; }

	// Returning the last group index that was used.
	return output_index-1;	
}

