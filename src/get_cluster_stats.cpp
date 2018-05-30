#include "csaw.h"
#include "utils.h"

SEXP get_cluster_stats (SEXP fcs, SEXP pvals, SEXP by, SEXP weight, SEXP fcthreshold) {
    BEGIN_RCPP

	// Checking indices.
    const Rcpp::List fclist(fcs);
	const size_t fcn=fclist.size();
    const Rcpp::NumericVector pval(pvals);
	const size_t nwin=pval.size();

	// Setting up the log-FC columns.
    std::vector<Rcpp::NumericVector> fcs(fcn);
    for (size_t i=0; i<fcn; ++i) { 
        const Rcpp::NumericVector current=fclist[i];
		if (nwin!=current.size()) { 
            throw std::runtime_error("log-FC and p-value vectors have different lengths");
        }
		fcs[i]=current;
	}

	// Setting up the remaining inputs. 
	const double fcthresh=check_numeric_scalar(fcthreshold, "log fold-change threshold");

    const Rcpp::IntegerVector clustid(by);
    size_t nclust=check_By_vector(clustid.begin(), clustid.end());
    if (clustid.size()!=nwin) {
        throw std::runtime_error("cluster ID and p-value vectors have different lengths");
    }

    const Rcpp::NumericVector winweight(weight);
	if (nwin!=winweight.size()) { 
        throw std::runtime_error("weight and p-value vectors have different lengths");
    }

	// Creating output objects to hold the results.
    Rcpp::IntegerVector out_nwin(nclust);
    std::vector<Rcpp::IntegerVector> out_nfc(fcn*2);
    for (size_t idx=0; idx<out_nfc.size(); ++idx) {
        out_nfc[idx]=Rcpp::IntegerVector(nclust); // avoid shallow copies.
    }
    
    Rcpp::NumericVector out_p(nclust);
    Rcpp::IntegerVector out_dir(fcn==1 ? nclust : 0);

    // Various temporary values.
    size_t curclust=0;
    size_t run_start=0;
    std::deque<std::pair<double, int> > psorter;

    while (run_start < nwin) {
        psorter.push_back(std::make_pair(pval[run_start], run_start));
        size_t run_end=run_start + 1;
        double subweight=winweight[run_start];
		while (run_end < nwin && clustid[run_start]==clustid[run_end]) { 
			subweight+=winweight[run_end];
            psorter.push_back(std::make_pair(pval[run_end], run_end));
			++run_end; 
		}

		// Computing the nclust number of windows, and that up/down for each fold change.
        out_nwin[curclust]=run_end - run_start;
		for (size_t fc=0; fc<fcn; ++fc) { 
			int& allup=out_nfc[fc*2][curclust];
			int& alldown=out_nfc[fc*2+1][curclust];

			for (size_t curwin=run_start; curwin<run_end; ++curwin) {  
                const double& curfc=fcs[fc][curwin];
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
        double remaining=winweight[psorter.front().second];
		double& outp=(out_p[curclust]=psorter.front().first/remaining); 
        auto itps=psorter.begin()+1;
        auto minit=psorter.begin();

        while (itps!=psorter.end()) {
	    	const int& curwin=itps->second;
			remaining+=winweight[curwin];
			const double more_temp=(itps->first)/remaining;
			if (more_temp < outp) { 
                outp=more_temp; 
                minit=itps;
            }
            ++itps;
		}
		outp*=subweight;

        /* If there's only one log-FC, we also determine which directions contribute to the combined p-value.
         * This is done by looking at the direction of the tests with p-values below that used as the combined p-value.
         * These tests must contribute because if any of them were non-significant, the combined p-value would increase.
         * Output codes are only up (1), only down (2) or mixed, i.e., both up and down (0).
         */
        if (fcn==1) {
            bool has_up=false, has_down=false;
            ++minit;

            for (auto itps=psorter.begin(); itps!=minit; ++itps) {
                const double& curfc=fcs[0][itps->second];
                if (curfc > 0) { 
                    has_up=true; 
                } else if (curfc < 0) {
                    has_down=true;
                }
                if (has_up && has_down) { 
                    break; 
                }
            }

            if (has_up && !has_down) { 
                out_dir[curclust]=1; 
            } else if (has_down && !has_up) { 
                out_dir[curclust]=2; 
            }
        }

		// Setting it up for the next round.
        psorter.clear();
		++curclust;
        run_start=run_end;
	}

    // Ensuring we return _something_ that can be put into a DataFrame with nclust rows,
    // even if there are no log-fold changes to report.
    Rcpp::RObject out_nfc_value, out_dir_value;
    if (fcn) {
        out_nfc_value=Rcpp::List(out_nfc.begin(), out_nfc.end());
        out_dir_value=out_dir;
    } else {
        out_nfc_value=Rcpp::IntegerMatrix(nclust, 0);
        out_dir_value=Rcpp::IntegerMatrix(nclust, 0);
    }

    return Rcpp::List::create(out_nwin, out_nfc_value, out_p, out_dir_value);
    END_RCPP
}
