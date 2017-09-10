#include "bam_utils.h"
#include "utils.h"

SEXP extract_single_data(SEXP bam, SEXP index, SEXP chr, SEXP start, SEXP end, 
        SEXP mapq, SEXP dedup, SEXP use_forward, SEXP use_first) {
    BEGIN_RCPP

    // Checking input values.
    Rcpp::IntegerVector _mapq(mapq);
    if (_mapq.size()!=1) { 
        throw std::runtime_error("mapping quality should be an integer scalar");
    }    
    const int minqual=_mapq[0];

    Rcpp::LogicalVector _dedup(dedup);
    if (_dedup.size()!=1) {    
        throw std::runtime_error("duplicate removal should be a logical scalar"); 
    }
    const bool rmdup=_dedup[0];

    Rcpp::LogicalVector _use_forward(use_forward);
    if (_use_forward.size()!=1) {    
        throw std::runtime_error("forward usage specification should be a logical scalar"); 
    }
    const int useforward=_use_forward[0];
    int set_flags=0, unset_flags=BAM_FSUPPLEMENTARY + BAM_FSECONDARY;
    if (useforward!=NA_LOGICAL) {
        if (useforward) {
            unset_flags+=BAM_FREVERSE;                
        } else {
            set_flags+=BAM_FREVERSE;                
        }
    }

    Rcpp::LogicalVector _use_first(use_first);
    if (_use_first.size()!=1) {
        throw std::runtime_error("first usage specification should be a logical scalar"); 
    } 
    const int usefirst=_use_first[0];
    if (usefirst!=NA_LOGICAL) {
        set_flags+=BAM_FPAIRED;
        if (usefirst) {
            set_flags+=BAM_FREAD1;
            unset_flags+=BAM_FREAD2;
        } else {
            set_flags+=BAM_FREAD2;
            unset_flags+=BAM_FREAD1;
        }
    }

    // Initializing odds and ends.
    BamFile bf(bam, index);
    BamRead br;
    BamIterator biter(bf, chr, start, end);
    std::deque<int> forward_pos, forward_len, reverse_pos, reverse_len;
    AlignData algn_data;

    while (bam_itr_next(bf.in, biter.iter, br.read) >= 0){    
//        // If we can see that it is obviously unmapped (IMPOSSIBLE for a sorted file).
//        if (((br.read -> core).flag & BAM_FUNMAP)!=0) { 
//            continue;
//        } 
        
        if (!br.is_well_mapped(minqual, rmdup)) { continue; }
        if (((br.read->core).flag & set_flags)!=set_flags) { continue; }
        if (((br.read->core).flag & unset_flags)!=0) { continue; }
        int curpos = (br.read -> core).pos + 1;
        br.extract_data(algn_data);
        
        if (algn_data.is_reverse) { 
            reverse_pos.push_back(curpos);
            reverse_len.push_back(algn_data.len);
        } else {
            forward_pos.push_back(curpos);
            forward_len.push_back(algn_data.len);
        }        
    }

    // Storing all output.
    Rcpp::List output(2);
    output[0]=Rcpp::List::create(
        Rcpp::IntegerVector(forward_pos.begin(), forward_pos.end()),
        Rcpp::IntegerVector(forward_len.begin(), forward_len.end())
    );
    output[1]=Rcpp::List::create(
        Rcpp::IntegerVector(reverse_pos.begin(), reverse_pos.end()),
        Rcpp::IntegerVector(reverse_len.begin(), reverse_len.end())
    );
 
    return output;
    END_RCPP
}

