#include "csaw.h"
#include "utils.h"

/* INTERNAL UTILITIES */

typedef std::deque<std::pair<double, int> > IndexedPValues;

class SimesPreparer {
public:
    template<class V>
    std::pair<double, size_t> operator()(IndexedPValues& psorter, const V& winweight) {
        /* Computing the (weighted) FDR and thus the (weighted) Simes value.
         * The weights are implemented as frequency weights, e.g., if you had 2
         * tests with a weight of 10 to 1, you'd consider the one with the
         * higher weight 10 more times to try to reject the global null (i.e.,
         * expanding it in-place in the sorted vector of p-values).
         */
        std::sort(psorter.begin(), psorter.end());

        double cumweight=0;
        for (auto pIt=psorter.begin(); pIt!=psorter.end(); ++pIt) {
            cumweight += winweight[pIt->second];
            (pIt->first)/=cumweight;
        }
        const double& total_weight=cumweight;

        // Backtracking to create adjusted p-values with a cumulative minimum.
        double curmin=1;
        size_t counter=psorter.size()-1, minindex=counter;
        for (auto prIt=psorter.rbegin(); prIt!=psorter.rend(); ++prIt, --counter) {
            double& current=(prIt->first);
            current*=total_weight;
            if (current < curmin) {
                curmin=current;
                minindex=counter;
            } else {
                current=curmin;
            }
        }

        return std::make_pair(curmin, minindex);
    }
};

class HolmPreparer {
public:    
    HolmPreparer(size_t mn, double mp) : min_num(std::max(size_t(1), mn)), min_prop(mp) {}

    template<class V>
    std::pair<double, size_t> operator()(IndexedPValues& psorter, const V& winweight) {
        /* Computing the (weighted) Holm correction. Weights are implemented
         * as scaling factors on the nominal type I error threshold, as described
         * by BH in the weighting paper form the Scandanavian Journal of Stats.
         */
        double total_weight=0;
        for (auto pIt=psorter.begin(); pIt!=psorter.end(); ++pIt) {
            total_weight += winweight[pIt->second];
        }

        std::sort(psorter.begin(), psorter.end());

        double remaining=total_weight;
        double cummax=0;
        for (auto pIt=psorter.begin(); pIt!=psorter.end(); ++pIt) {
            const double curweight = winweight[pIt->second];

            double& current=(pIt->first);
            current *= remaining/curweight;
            if (current > 1) {
                current=1;
            }
            if (current > cummax) {
                cummax=current;
            } else {
                current=cummax;
            }

            remaining -= curweight;
        }


        size_t index=std::max(
            min_num, 
            static_cast<size_t>(std::ceil(min_prop * static_cast<double>(psorter.size())))
        );
        index=std::min(index, psorter.size());
        return std::make_pair(psorter[index].first, index-1);
    }
private:
    const size_t min_num;
    const double min_prop;
};

class MaxedPreparer {
public:    
    MaxedPreparer(const Rcpp::NumericVector& m) : metric(m) {}

    template<class V>
    std::pair<double, size_t> operator()(IndexedPValues& psorter, const V& winweight) {
        double maxed=R_NegInf;
        auto chosen=psorter.begin();

        for (auto pIt=psorter.begin(); pIt!=psorter.end(); ++pIt) {
            const auto& current=metric[pIt->second];
            if (current > maxed) {
                maxed=current;
                chosen=pIt;
            }
        }

        // Wiping out everything that is NOT the top test with respect to some
        // independent filter statistic.
        for (auto pIt=psorter.begin(); pIt!=psorter.end(); ++pIt) {
            if (pIt!=chosen) { pIt->first=1; }
        }

        return std::make_pair(chosen->first, (chosen - psorter.begin()));
    }
private:
    const Rcpp::NumericVector& metric;
};

template<class V>
int guess_direction(const IndexedPValues& psorter, size_t minit, const V& lfc) {
    /* If there's only one log-FC, we also determine which directions
     * contribute to the combined p-value.  This is done by looking at the
     * direction of the tests with p-values below that used as the combined
     * p-value.  These tests must contribute because if any of them were
     * non-significant, the combined p-value would increase.  Output codes
     * are only up (1), only down (2) or mixed, i.e., both up and down (0).
     */
    bool has_up=false, has_down=false;

    for (size_t i=0; i<=minit; ++i) {
        const double& curfc=lfc[psorter[i].second];
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
        return 1; 
    } else if (has_down && !has_up) { 
        return 2; 
    } else {
        return 0;
    }
}

/* MAIN FUNCTION */

template<class PREP>
SEXP get_cluster_stats_internal (SEXP fcs0, SEXP pvals, SEXP by, SEXP weight, SEXP fcthreshold, PREP& preparer) {
	// Checking indices.
    const Rcpp::List fclist(fcs0);
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
    Rcpp::IntegerVector out_rep(nclust);

    // Various temporary values.
    size_t curclust=0;
    size_t run_start=0;
    std::deque<std::pair<double, int> > psorter;

    while (run_start < nwin) {
        psorter.push_back(std::make_pair(pval[run_start], run_start));
        size_t run_end=run_start + 1;
		while (run_end < nwin && clustid[run_start]==clustid[run_end]) { 
            psorter.push_back(std::make_pair(pval[run_end], run_end));
			++run_end; 
		}

        auto output=preparer(psorter, winweight);
        out_p[curclust] = output.first;
        out_nwin[curclust] = run_end - run_start;
        out_rep[curclust] = psorter[output.second].second;

		// Computing the nclust number of windows, and that up/down for each fold change.
		for (size_t fc=0; fc<fcn; ++fc) { 
			int& allup=out_nfc[fc*2][curclust];
			int& alldown=out_nfc[fc*2+1][curclust];
            const auto& curfcs=fcs[fc];

            for (auto pIt=psorter.begin(); pIt!=psorter.end(); ++pIt) {
                if (pIt->first <= fcthresh) {
                    const double& curfc=curfcs[pIt->second];
                    if (curfc > 0) { 
                        ++allup; 
                    } else if (curfc < 0) { 
                        ++alldown; 
                    }
                }
			}
		}

        if (fcn==1) {
            out_dir[curclust]=guess_direction(psorter, output.second, fcs[0]);
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

    return Rcpp::List::create(out_nwin, out_nfc_value, out_p, out_dir_value, out_rep);
}

SEXP compute_cluster_simes(SEXP fcs, SEXP pvals, SEXP by, SEXP weight, SEXP fcthreshold) {
    BEGIN_RCPP
    SimesPreparer prep;
    return get_cluster_stats_internal(fcs, pvals, by, weight, fcthreshold, prep);
    END_RCPP
}

SEXP compute_cluster_holm(SEXP fcs, SEXP pvals, SEXP by, SEXP weight, SEXP fcthreshold, SEXP min_n, SEXP min_p) {
    BEGIN_RCPP
    HolmPreparer prep(
        check_integer_scalar(min_n, "minimum number of tests"),
        check_numeric_scalar(min_n, "minimum proportion of tests")
    );
    return get_cluster_stats_internal(fcs, pvals, by, weight, fcthreshold, prep);
    END_RCPP
}

SEXP compute_cluster_maxed(SEXP fcs, SEXP pvals, SEXP by, SEXP weight, SEXP fcthreshold, SEXP metric) {
    BEGIN_RCPP
    MaxedPreparer prep(metric);
    return get_cluster_stats_internal(fcs, pvals, by, weight, fcthreshold, prep);
    END_RCPP
}
