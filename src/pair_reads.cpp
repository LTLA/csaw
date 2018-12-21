#include "bam_utils.h"
#include "utils.h"
#include "intersector.h"

struct AlignData {
    AlignData() : len(0), is_reverse(false) {}
    AlignData(int L, bool R) : len(L), is_reverse(R) {}
    int len;
    bool is_reverse;
};

struct valid_pairs {
    bool add_pair(int pos1, const AlignData& data1, int pos2, const AlignData& data2, bool isfirst1) {
        if (data2.is_reverse==data1.is_reverse) {
            return false;
        }

        int forward_pos, forward_len, reverse_pos, reverse_len;
        if (data2.is_reverse) {
            forward_pos = pos1;
            forward_len = data1.len;
            reverse_pos = pos2;
            reverse_len = data2.len;
        } else {
            forward_pos = pos2;
            forward_len = data2.len;
            reverse_pos = pos1;
            reverse_len = data1.len;
        }

        // Completely overrun fragments go into the 'unoriented' basket.
        if (forward_pos >= reverse_pos + reverse_len) {
            return false; 
        } 

        forward_pos_out.push_back(forward_pos); 
        forward_len_out.push_back(forward_len);
        reverse_pos_out.push_back(reverse_pos);
        reverse_len_out.push_back(reverse_len);
        return true;
    }

    std::deque<int> forward_pos_out, forward_len_out, reverse_pos_out, reverse_len_out;
};

/* Strolls through the file for each chromosome and accumulates paired-end statistics; 
 * forward and reverse reads (position and width), singles, unoriented, and names of
 * inter-chromosomals.
 */

SEXP extract_pair_data(SEXP bam, SEXP index, 
        SEXP chr, SEXP start, SEXP end, 
        SEXP mapq, SEXP dedup, 
        SEXP discard_pos, SEXP discard_id,
        SEXP diagnostics)
{
    BEGIN_RCPP

    // Checking input values.
    const int minqual=check_integer_scalar(mapq, "minimum mapping quality");
    const bool rmdup=check_logical_scalar(dedup, "duplicate removal specification");
    const bool storediags=check_logical_scalar(diagnostics, "diagnostics specification");
    intersector discarder(discard_pos, discard_id);

    // Initializing odds and ends.
    BamFile bf(bam, index);
    BamRead br;
    BamIterator biter(bf, chr, start, end);

    typedef std::map<std::pair<int, std::string>, AlignData> Holder;
    std::deque<Holder> all_holders(4); // four holders, one for each strand/first combination; cut down searches.
    std::pair<int, std::string> current;

    std::set<std::string> identical_pos;
    int last_identipos=-1;

    valid_pairs valid;
    int totals=0, num_singles=0, num_onemapped=0, num_unoriented=0;
    std::deque<std::string> interchr_names_1, interchr_names_2;

    while (bam_itr_next(bf.in, biter.iter, br.read) >= 0){
        auto curflag=br.get_flag();
        if ((curflag & BAM_FSECONDARY)!=0 || (curflag & BAM_FSUPPLEMENTARY)!=0) {
            continue; // These guys don't even get counted as reads.
        }

        ++totals;
        const int curpos = br.get_aln_pos() + 1; // Getting 1-indexed position.
        const AlignData algn_data(br.get_aln_len(), br.is_reverse());

        discarder.advance_to_start(curpos);
        const bool am_mapped=br.is_well_mapped(minqual, rmdup) & 
            !discarder.end_is_within(curpos + algn_data.len);

        /* Reasons to not add a read: */
       
        // If it's a singleton.
        if ((curflag & BAM_FPAIRED)==0) {
            if (storediags && am_mapped) { ++num_singles; }
            continue;
        }

        // Or, if we can see that its partner is obviously unmapped.
        if ((curflag & BAM_FMUNMAP)!=0) {
            if (storediags && am_mapped) { ++num_onemapped; }
            continue;
        }

        // Or if it's inter-chromosomal.
        const bool is_first=((curflag & BAM_FREAD1)!=0);
        if (is_first==((curflag & BAM_FREAD2)!=0)) { 
            std::stringstream err;
            err << "read '" << bam_get_qname(br.read) << "' must be either first or second in the pair";
            throw std::runtime_error(err.str()); 
        }
      
        if ((br.read -> core).mtid!=(br.read -> core).tid) { 
            if (storediags && am_mapped) { 
                auto& interchr_names=(is_first ? interchr_names_1 : interchr_names_2);
                interchr_names.push_back(bam_get_qname(br.read)); 
            } 
            continue;
        }

        /* Checking the map for its mate, adding it if it doesn't exist. */
        
        current.second.assign(bam_get_qname(br.read));
        const int mate_pos = (br.read -> core).mpos + 1; // 1-indexed position, again.
        bool mate_is_in=false;
        if (mate_pos < curpos) {
            mate_is_in=true;
        } else if (mate_pos == curpos) {
            // Identical mpos to curpos needs careful handling to figure out whether we've already seen it.
            if (curpos!=last_identipos) { 
                identical_pos.clear();
                last_identipos=curpos;
            }

            auto itip=identical_pos.lower_bound(current.second);
            if (itip!=identical_pos.end() && !(identical_pos.key_comp()(current.second, *itip))) {
                mate_is_in=true;
                identical_pos.erase(itip);
            } else {
                identical_pos.insert(itip, current.second);
            }
        }

        if (mate_is_in) {
            current.first = mate_pos;
            Holder& holder=all_holders[int(!is_first) + 2*int(bam_is_mrev(br.read))];
            Holder::iterator ith=holder.find(current);

            if (ith != holder.end()) { 
                if (!am_mapped) {
                    // Searching to pop out the mate, to reduce the size of 'holder' for the remaining searches.
                    if (storediags) { ++num_onemapped; }
                    holder.erase(ith);
                    continue;
                }

                if (!valid.add_pair(curpos, algn_data, (ith->first).first, ith->second, is_first)) {
                    ++num_unoriented;
                }
                holder.erase(ith);

            } else if (am_mapped) {
                // Only possible if the mate didn't get added because 'am_mapped' was false.
                if (storediags) { ++num_onemapped; }
            }
        } else if (am_mapped) {
            current.first = curpos;
            Holder& holder=all_holders[int(is_first) + 2*int(algn_data.is_reverse)];
            holder[current] = algn_data;
        }
    }

    // Leftovers treated as one_unmapped; marked as paired, but the mate is not in file.
    if (storediags) {
        for (auto& holder : all_holders) {
            num_onemapped+=holder.size();
        }    
    }

    // Storing all output.
    if (storediags) {
        return Rcpp::List::create(
            Rcpp::List::create(
                Rcpp::IntegerVector(valid.forward_pos_out.begin(), valid.forward_pos_out.end()),
                Rcpp::IntegerVector(valid.forward_len_out.begin(), valid.forward_len_out.end())
            ),
            Rcpp::List::create(
                Rcpp::IntegerVector(valid.reverse_pos_out.begin(), valid.reverse_pos_out.end()),
                Rcpp::IntegerVector(valid.reverse_len_out.begin(), valid.reverse_len_out.end())
            ),
            Rcpp::IntegerVector::create(totals),
            Rcpp::IntegerVector::create(num_singles),
            Rcpp::IntegerVector::create(num_unoriented),
            Rcpp::IntegerVector::create(num_onemapped),
            Rcpp::List::create(
                Rcpp::StringVector(interchr_names_1.begin(), interchr_names_1.end()),
                Rcpp::StringVector(interchr_names_2.begin(), interchr_names_2.end())
            )
        );
    } else {
        return Rcpp::List::create(
            Rcpp::List::create(
                Rcpp::IntegerVector(valid.forward_pos_out.begin(), valid.forward_pos_out.end()),
                Rcpp::IntegerVector(valid.forward_len_out.begin(), valid.forward_len_out.end())
            ),
            Rcpp::List::create(
                Rcpp::IntegerVector(valid.reverse_pos_out.begin(), valid.reverse_pos_out.end()),
                Rcpp::IntegerVector(valid.reverse_len_out.begin(), valid.reverse_len_out.end())
            )
        );
    }
    END_RCPP
}

/* Getting reads on other unprocessed chromosomes, unmapped reads. */

SEXP get_leftovers (SEXP bam, SEXP index, SEXP processed) { 
    BEGIN_RCPP
    BamFile bf(bam, index);
    BamRead br;

    Rcpp::StringVector _processed(processed);
    const int nchr=_processed.size();
    std::set<std::string> already_there;
    for (int i=0; i<nchr; ++i) {
        already_there.insert(Rcpp::as<std::string>(_processed[i]));
    }

    // Getting the reads mapped to chromosomes we didn't look at due to 'restrict'.
    int leftovers=0;
    std::set<std::string>::iterator iat;
    for (int cid=0; cid<bf.header->n_targets; ++cid) {
        iat=already_there.find(std::string(bf.header->target_name[cid]));
        if (iat!=already_there.end()) { continue; }
        BamIterator biter(bf, cid, 0, bf.header->target_len[cid]);
        while (bam_itr_next(bf.in, biter.iter, br.read) >= 0){ ++leftovers; }
    } 
    
    // Also getting the unmapped guys. 
    BamIterator biter(bf);
    while (bam_itr_next(bf.in, biter.iter, br.read) >= 0){ ++leftovers; }
    return Rcpp::IntegerVector::create(leftovers);
    END_RCPP
}

