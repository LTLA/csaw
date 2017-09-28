#include "csaw.h"
#include "utils.h"

enum posttype { START, MIDSTART, MIDEND, END };

struct signpost {
	signpost(int p, posttype t, int l, int i) : position(p), type(t), library(l), index(i) {}
	int position;
	posttype type;
	int library, index;
	bool operator> (const signpost& right) const {
		return (position > right.position);
	}
};

SEXP check_bimodality (SEXP all, SEXP regstart, SEXP regend, SEXP priorcount, SEXP invert) {
    BEGIN_RCPP

	// Setting structures for the data.
    const Rcpp::List _all(all);
    const int nlibs=_all.size();
    std::vector<Rcpp::IntegerVector> left1(nlibs), right1(nlibs), left2(nlibs), right2(nlibs), strand(nlibs);
    std::vector<int> nums(nlibs), indices(nlibs);
    std::priority_queue<signpost, std::deque<signpost>, std::greater<signpost> > next;
	
	for (int i=0; i<nlibs; ++i) {
        const Rcpp::List current=_all[i];
        if (current.size()!=5) { 
			throw std::runtime_error("fragment data must be supplied as a list of length 5"); 
        }
		
		for (int j=0; j<5; ++j) {
            const Rcpp::IntegerVector current_col=current[j];
			switch (j) {
				case 0: 
					left1[i]=current_col;
					nums[i]=current_col.size();
                    break;
				case 1:
					right1[i]=current_col;
                    break;
				case 2: 
					left2[i]=current_col; 
                    break;
				case 3:
					right2[i]=current_col;
                    break;
				case 4:
					strand[i]=current_col;
                    break;
			}
            if (current_col.size()!=nums[i]) { 
                throw std::runtime_error("length of vectors must be equal"); 
            }
		}
		
		// Populating the priority queue.
		if (nums[i]) { 
            next.push(signpost(left1[i][0], START, i, 0)); 
            next.push(signpost(left2[i][0], MIDSTART, i, 0)); 
        }
	}

	// Setting up structures for the regions.
    Rcpp::IntegerVector _regstart(regstart), _regend(regend);
	const int nregs=_regstart.size();
	if (nregs!=_regend.size()) { 
        throw std::runtime_error("length of region vectors should be equal"); 
    }
	int reg_index=0, next_regstart=-1;
	if (nregs) { 
		next.push(signpost(_regstart[0], START, -1, 0));
		next_regstart=_regstart[0];
	}

	// Setting up the prior count and inversion flag.
    const double pc=check_numeric_scalar(priorcount, "prior count");
	const bool inv=check_logical_scalar(invert, "invert specification");

    // Reporting bimodality output.
    Rcpp::NumericVector output(nregs);
	std::set<int> current_regs;
	int left_forward=0, left_reverse=0, right_forward=0, right_reverse=0;
	double current_score=1; // default when no reads added; prior counts cancel out.
	std::deque<int> new_regs;

	// Running through the set; stopping when there are no regions, or when everything is processed.
	while (!next.empty()) { 
		int current_position=next.top().position;
        bool modified_stats=false;
//		Rprintf("At position %i\n", current_position);

		// Pulling out all features at this current position.
		do {
			int current_library=next.top().library;
			posttype current_type=next.top().type;
			int current_index=next.top().index;
//			Rprintf("\t\t processing %i: %i -> %i\n", current_type, current_library, current_index);
			next.pop();

			if (current_library<0) {
				// Using negative library values to encode the region index.
				if (current_type==START) {
					current_regs.insert(current_index);
					new_regs.push_back(current_index);
					if ((++reg_index) < nregs) { 
						next.push(signpost(_regstart[reg_index], START, -1, reg_index));
						next_regstart=_regstart[reg_index];
					} else {
						next_regstart=-1;
					}
					next.push(signpost(_regend[current_index]+1, END, -1, current_index));
				} else {
					current_regs.erase(current_index);
				}

			} else {
				// Recording the number of forward/reverse reads to the left or right of the current position.
				modified_stats=true;
				const int& isforward=strand[current_library][current_index];
				switch (current_type) { 
					case START:
						if (isforward) { ++right_forward; }
						else { ++right_reverse; }
						next.push(signpost(right1[current_library][current_index] + 1, MIDEND, current_library, current_index));
						break;
					case MIDSTART:
						if (isforward) { ++left_forward; } 
						else { ++left_reverse; }
						next.push(signpost(right2[current_library][current_index] + 1, END, current_library, current_index));
						break;
					case MIDEND:
						if (isforward) { --right_forward; }
						else { --right_reverse; }
						break;
					case END:
						if (isforward) { --left_forward; }
						else { --left_reverse; }
						break;
					default:
						break;
				}

				// Adding the next thing (skipping intervening reads that don't affect any regions).
				if (current_type==START) {
					int& next_index=indices[current_library];
					while ((++next_index) < nums[current_library]) { 
						if (current_regs.empty() && next_regstart >=0 && next_regstart >= right2[current_library][next_index] + 1) {
 							continue; 
						}
						next.push(signpost(left1[current_library][next_index], START, current_library, next_index));
						next.push(signpost(left2[current_library][next_index], MIDSTART, current_library, next_index));
						break;
					}
				}
			}
		} while (!next.empty() && next.top().position==current_position);

		// Running through and asking if the updated bimodality score is higher or lower.
		if (modified_stats) {
			if (inv) {
				current_score = std::min((double(left_reverse)+pc)/(double(left_forward)+pc), 				
					(double(right_forward)+pc)/(double(right_reverse)+pc));
			} else {
				current_score = std::min((double(left_forward)+pc)/(double(left_reverse)+pc), 				
					(double(right_reverse)+pc)/(double(right_forward)+pc));
			}
		}
//		Rprintf("\t values are Left forward/reverse: %i/%i, right reverse/forward %i/%i\n", left_forward, left_reverse, right_reverse, right_forward);
//		Rprintf("\t score is %.3f\n", current_score);
		if (!new_regs.empty()) { 
            for (const auto& nr : new_regs) { 
                output[nr]=current_score; 
            }
			new_regs.clear();
		}
		if (modified_stats) { 
			for (const auto& r : current_regs) { 
                double& curout=output[r];
 			    if (curout < current_score) { curout=current_score; }
			}
		}
 	   	
		// Quitting if there's no more regions to process.
		if (current_regs.empty() && next_regstart < 0) { break; }
	}

	return output;
    END_RCPP
}
