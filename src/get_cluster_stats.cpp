#include "csaw.h"
#include "utils.h"

SEXP get_cluster_stats (SEXP fcdex, SEXP pvaldex, SEXP tab, SEXP by, SEXP weight, SEXP fcthreshold) {
    BEGIN_RCPP

	// Checking indices.
    Rcpp::IntegerVector _fcdex(fcdex), _pvaldex(pvaldex); 
	if (_pvaldex.size()!=1) { throw std::runtime_error("only one index should be supplied for the p-value column"); }
	const int pdex=_pvaldex[0];
	const int fcn=_fcdex.size();

	// Setting up the p-value column.
    Rcpp::List _tab(tab);
    if (pdex < 0 || pdex >= _tab.size()) { 
        throw std::runtime_error("p-value index out of range"); 
    }
    Rcpp::NumericVector pval(_tab[pdex]);
	const int n=pval.size();

	// Setting up the log-FC columns.
    std::vector<Rcpp::NumericVector> fcs(fcn);
    for (int i=0; i<fcn; ++i) { 
        const int& curfcdex=_fcdex[i];
        if (curfcdex < 0 || curfcdex >= _tab.size()) { 
            throw std::runtime_error("log-FC index out of range"); 
        }

		fcs[i]=Rcpp::NumericVector(_tab[curfcdex]);
		if (n!=fcs[i].size()) { 
            throw std::runtime_error("vector lengths are not equal"); 
        }
	}

    Rcpp::NumericVector _fcthreshold(fcthreshold);
	if (_fcthreshold.size()!=1) { 
        throw std::runtime_error("log-fold change threshold should be a numeric scalar"); 
    }
	const double fcthresh=_fcthreshold[0];

	// Setting up the remaining inputs. 
    Rcpp::IntegerVector _by(by);
    Rcpp::NumericVector _weight(weight);
	if (n!=_by.size() || n!=_weight.size()) { 
        throw std::runtime_error("vector lengths are not equal"); 
    }
    int total=checkByVector(_by.begin(), _by.end());

	// Pulling out results.
    Rcpp::IntegerVector out_nwin(total);
    Rcpp::IntegerMatrix out_nfc(total, fcn*2);
    Rcpp::NumericVector out_p(total);
    Rcpp::IntegerVector out_dir(fcn==1 ? total : 0);
    int element=0;

    // Various temporary values.
    int i=0;
    std::deque<std::pair<double, int> > psorter;
    while (i<n) {
        psorter.push_back(std::make_pair(pval[i], i));
        int j=i+1;
        double subweight=_weight[i];
		while (j < n && _by[i]==_by[j]) { 
			subweight+=_weight[j];
            psorter.push_back(std::make_pair(pval[j], j));
			++j; 
		}

		// Computing the total number of windows, and that up/down for each fold change.
        out_nwin[element]=j-i;
        Rcpp::IntegerMatrix::Row cur_nfc=out_nfc.row(element);
		for (int k=0; k<fcn; ++k) { 
			int& allup=cur_nfc[k*2];
			int& alldown=cur_nfc[k*2+1];

			for (int x=i; x<j; ++x) {  
                const double& curfc=fcs[k][x];
				if (curfc > fcthresh) { 
                    ++allup; 
                } else if (curfc < -fcthresh) { 
                    ++alldown; 
                }
			}
		}

		/* Computing the weighted Simes value. The weights are implemented as frequency 
		 * weights, e.g., if you had 2 tests with a weight of 10 to 1, you'd consider the
		 * one with the higher weight 10 more times to try to reject the global null (i.e.,
		 * expanding it in-place in the sorted vector of p-values).
		 */
		std::sort(psorter.begin(), psorter.end());
        double remaining=_weight[psorter.front().second];
		double& outp=(out_p[element]=psorter.front().first/remaining); 
        std::deque<std::pair<double, int> >::iterator itps=psorter.begin()+1;
 
        int minx=i;
		for (int x=i+1; x<j; ++x, ++itps) {
	    	const int& current=itps->second;
			remaining+=_weight[current];
			const double more_temp=(itps->first)/remaining;
			if (more_temp < outp) { 
                outp=more_temp; 
                minx=x;
            }
		}
		outp*=subweight;

        /* If there's only one log-FC, we also determine which directions contribute to the combined p-value.
         * This is done by looking at the direction of the tests with p-values below that used as the combined p-value.
         * These tests must contribute because if any of them were non-significant, the combined p-value would increase.
         * Output codes are only up (1), only down (2) or mixed, i.e., both up and down (0).
         */
        if (fcn==1) {
            bool has_up=false, has_down=false;
            std::deque<std::pair<double, int> >::iterator itps=psorter.begin();
            for (int x=i; x<=minx; ++x, ++itps) {
                const double& curfc=fcs[0][itps->second];
                if (curfc > 0) { 
                    has_up=true; 
                } else if (curfc < 0) {
                    has_down=true;
                }
                if (has_up & has_down) { 
                    break; 
                }
            }

            int& current_dir=out_dir[element];
            if (has_up & !has_down) { 
                current_dir=1; 
            } else if (has_down & !has_up) { 
                current_dir=2; 
            }
        }

		// Setting it up for the next round.
        psorter.clear();
		++element;
        i=j;
	}
	
    return Rcpp::List::create(out_nwin, out_nfc, out_p, out_dir);
    END_RCPP
}

/* Computes the total cluster weight in a reasonably fast manner. */

SEXP get_cluster_weight(SEXP ids, SEXP weight) {
    BEGIN_RCPP

    Rcpp::IntegerVector _ids(ids);
    Rcpp::NumericVector _weight(weight);
    const int n=_ids.size();
	if (n!=_weight.size()) { 
        throw std::runtime_error("vector lengths are not equal"); 
    }
    int total=checkByVector(_ids.begin(), _ids.end());

    Rcpp::NumericVector output(total);
    if (total) { 
        Rcpp::NumericVector::iterator oIt=output.begin(), wIt=_weight.begin();
        *oIt=*wIt;
        ++wIt;

        for (int i=1; i<n; ++i) {
            if (_ids[i]!=_ids[i-1]) { 
                ++oIt;
            }
            (*oIt)+=*wIt;
            ++wIt;
        }
    }

    return output;
    END_RCPP
}

