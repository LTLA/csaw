#include "bam_utils.h"
#include "utils.h"
#include "intersector.h"

SEXP extract_single_data(SEXP bam, SEXP index, 
        SEXP chr, SEXP start, SEXP end, 
        SEXP mapq, SEXP dedup, SEXP use_forward, SEXP use_first, 
        SEXP discard_pos, SEXP discard_id) 
{
    BEGIN_RCPP

    // Checking input values.
    const int minqual=check_integer_scalar(mapq, "minimum mapping quality");
    const bool rmdup=check_logical_scalar(dedup, "duplicate removal specification");
    int set_flags=0, unset_flags=BAM_FSUPPLEMENTARY + BAM_FSECONDARY;

    const Rcpp::LogicalVector _use_forward(use_forward); // Do NOT replace with check_logical_scalar, need to check for NA!
    if (_use_forward.size()!=1) {    
        throw std::runtime_error("forward usage specification should be a logical scalar"); 
    }
    const int useforward=_use_forward[0];
    if (useforward!=NA_LOGICAL) {
        if (useforward) {
            unset_flags+=BAM_FREVERSE;                
        } else {
            set_flags+=BAM_FREVERSE;                
        }
    }

    const Rcpp::LogicalVector _use_first(use_first); // Do NOT replace with check_logical_scalar, need to check for NA!
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

    intersector discarder(discard_pos, discard_id);

    // Initializing odds and ends.
    BamFile bf(bam, index);
    BamRead br;
    BamIterator biter(bf, chr, start, end);
    std::deque<int> forward_pos, forward_len, reverse_pos, reverse_len;

    while (bam_itr_next(bf.in, biter.iter, br.read) >= 0){    
        // No need to check for unmapped directly; this should be IMPOSSIBLE for a sorted file.
        if (!br.is_well_mapped(minqual, rmdup)) { continue; }
        auto cur_flag = br.get_flag();
        if ((cur_flag & set_flags)!=set_flags) { continue; }
        if ((cur_flag & unset_flags)!=0) { continue; }

        int curpos = br.get_aln_pos() + 1;
        int curlen = br.get_aln_len();

        discarder.advance_to_start(curpos);
        if (discarder.end_is_within(curpos + curlen)) { continue; }

        if (br.is_reverse()) { 
            reverse_pos.push_back(curpos);
            reverse_len.push_back(curlen);
        } else {
            forward_pos.push_back(curpos);
            forward_len.push_back(curlen);
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

