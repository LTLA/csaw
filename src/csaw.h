#ifndef CSAW_H
#define CSAW_H

#include "Rcpp.h"

#include <deque>
#include <vector>
#include <set>
#include <map>
#include <queue>

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include <string>
#include <sstream>

extern "C" {

/* annotator.cpp */

SEXP collate_exon_data (SEXP, SEXP, SEXP, SEXP);

SEXP annotate_overlaps (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* best_in_cluster.cpp */

SEXP best_in_cluster(SEXP, SEXP, SEXP);

/* get_cluster_stats.cpp */

SEXP compute_cluster_simes(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP compute_cluster_holm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP compute_cluster_maxed(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* merge_windows.cpp */

SEXP merge_windows(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* correlate_reads.cpp */

SEXP correlate_reads(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* get_rle_counts.cpp */

SEXP get_rle_counts(SEXP, SEXP, SEXP, SEXP, SEXP);

/* get_profile.cpp */

SEXP get_profile(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* find_maxima.cpp */

SEXP find_maxima(SEXP, SEXP, SEXP, SEXP, SEXP);

/* check_bimodality.cpp */

SEXP check_bimodality(SEXP, SEXP, SEXP, SEXP, SEXP);

/* pair_reads.cpp */

SEXP extract_pair_data(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP get_leftovers(SEXP, SEXP, SEXP);

/* single_reads.cpp */

SEXP extract_single_data(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

}

#endif

