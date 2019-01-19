#include "csaw.h"
#include "utils.h"

int split_cluster(const Rcpp::IntegerVector& starts, const Rcpp::IntegerVector& ends, Rcpp::IntegerVector& output, 
    const int& first, const int& last, const int& width);

/* We assume that incoming elements are sorted by chr -> start -> end. We then proceed 
 * to aggregate elements by chaining together elements that are less than 'tolerance'
 * apart and, if required, have the same 'sign'. 
 */

SEXP merge_windows(SEXP win_chrs, SEXP win_start, SEXP win_end, SEXP win_sign, SEXP tolerance, SEXP max_size) {
    BEGIN_RCPP

    const Rcpp::IntegerVector chrs(win_chrs), start(win_start), end(win_end);
    const Rcpp::LogicalVector sign(win_sign);
    const int nwin = chrs.size();
    if (nwin!=start.size() || nwin!=end.size() || nwin!=sign.size()) { 
        throw std::runtime_error("lengths of vectors are not equal"); 
    }

    const int tol=check_integer_scalar(tolerance, "tolerance");
    const Rcpp::IntegerVector maxregion(max_size);
    if (maxregion.size() > 1) { 
        throw std::runtime_error("maximum size should be an integer scalar or NULL"); 
    }
    const bool limit_size=(maxregion.size()==1);
    const int maxs=(limit_size ? maxregion[0] : 0);

    // Proceeding with the merge operation.
    Rcpp::IntegerVector out_index(nwin);
    int first_of_run=0;
    int ngroups=0;
    if (nwin) { 
        ngroups=1;
        out_index[0]=1;
    }

    for (int i=1; i<nwin; ++i) {
        if (chrs[i]!=chrs[i-1]                                  // Obviously, changing if we're on a different chromosome.
            || start[i] - end[i-1] - 1 > tol                    // Space between windows, start anew if this is greater than the tolerance.
            || sign[i]!=sign[i-1]                               // Checking if the sign is consistent.
           ) {
            if (limit_size) {
                ngroups = split_cluster(start, end, out_index, first_of_run, i, maxs);
            }
            ++ngroups;
            first_of_run=i;
        }

        out_index[i]=ngroups; 
    }

    // Cleaning up the last cluster, if necessary.
    if (limit_size && nwin) {
        ngroups = split_cluster(start, end, out_index, first_of_run, nwin, maxs);
    }

    // Now, identifying the chromosome, start and end of each region.
    Rcpp::IntegerVector out_chr(ngroups, -1), out_start(ngroups), out_end(ngroups);

    for (int i=0; i<nwin; ++i) { 
        int curgroup=out_index[i]-1;
        if (out_chr[curgroup]<0) { 
            out_chr[curgroup]=chrs[i];
            out_start[curgroup]=start[i]; // Sorted by start, remember; only need this once.
            out_end[curgroup]=end[i];
        } else if (out_end[curgroup] < end[i]) { 
            out_end[curgroup]=end[i]; 
        }
    }

    return Rcpp::List::create(out_index, out_chr, out_start, out_end);
    END_RCPP
}

int split_cluster(const Rcpp::IntegerVector& starts, const Rcpp::IntegerVector& ends, Rcpp::IntegerVector& output, 
        const int& first, const int& last, const int& width) {
    
    if (last - first == 1) {
        return output[first];
    }
    const int true_end=*std::max_element(ends.begin() + first, ends.begin() + last);
    const int full_width=true_end - starts[first] + 1;
    if (full_width <= width) { 
        return output[first];
    }

    /* We split up the cluster into evenly sized partitions where each partition is less than 'width'.
     * There are 'mult' partitions in total, each of width 'subwidth'.
     * Windows are assigned to partitions based on their midpoints. 
     * 
     * Note that there can only be `mult` subclusters. 
     * In the most extreme case, `cur_diff` will be equal to `true_end - starts[first]`. 
     * Division by `subwidth` will give an expression of `(true_end - starts[first]) / [ (true_end - starts[first] + 1) / mult ]`.
     * This will always be less than `mult`, so flooring it will give `mult-1`, i.e., the last index of `instantiated`.
     */
       const int mult=int(std::ceil(double(full_width)/width));
    std::vector<int> instantiated(mult, 0);
    int cur_group=output[first];

    const double subwidth=double(full_width)/mult;
    for (int i=first; i<last; ++i) {
        const double cur_diff = double(starts[i]+ends[i])*0.5 - starts[first];
        output[i] = int(cur_diff/subwidth);    
        if (!instantiated[output[i]]) { 
            instantiated[output[i]] = 1; 
        }
    }

    // Allocating output indices to the subclusters. This avoids nonconsecutive cluster indices.
    for (int i=0; i<mult; ++i) { 
        if (!instantiated[i]) { 
            continue; 
        }
        instantiated[i]=cur_group;
        ++cur_group;
    }

    // Assigning indices back to the output vector.
    for (int i=first; i<last; ++i) { 
        output[i]=instantiated[output[i]]; 
    }

    // Returning the new number of groups.
    return cur_group - 1;
}

