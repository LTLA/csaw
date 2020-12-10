#include "csaw.h"
#include "R_ext/Rdynload.h"
#include "R_ext/Visibility.h"

#define REGISTER(x, i) {#x, (DL_FUNC) &x, i}

extern "C" {

static const R_CallMethodDef all_call_entries[] = {
    REGISTER(annotate_overlaps, 8),
    REGISTER(best_in_cluster, 3),

    REGISTER(merge_windows, 6),
    REGISTER(correlate_reads, 6),
    REGISTER(get_rle_counts, 5), 
    REGISTER(get_profile, 6), 
    REGISTER(find_maxima, 5), 
    REGISTER(check_bimodality, 5), 

    REGISTER(extract_pair_data, 10), 
    REGISTER(get_leftovers, 3), 
    REGISTER(extract_single_data, 11), 
    {NULL, NULL, 0}
};

void attribute_visible R_init_csaw(DllInfo *dll) {
    R_registerRoutines(dll, NULL, all_call_entries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}

}
