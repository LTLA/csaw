#include "csaw.h"
#include "utils.h"

SEXP annotate_overlaps (SEXP N, SEXP Q, SEXP S, SEXP D, SEXP Gn, SEXP Gs, SEXP T, SEXP str) {
    BEGIN_RCPP

    const int nin=check_integer_scalar(N, "number of query regions");

    // Setting up overlap information.
    const Rcpp::IntegerVector query(Q), subject(S);
    const int nolaps=query.size();
	if (nolaps!=subject.size()){ 
        throw std::runtime_error("overlap vectors should have equal length"); 
    }

    const bool use_dist=(D!=R_NilValue);
    Rcpp::IntegerVector dist;
    if (use_dist) {
        dist=Rcpp::IntegerVector(D);
        if (nolaps!=dist.size()) { 
            throw std::runtime_error("overlap vectors should have the same length");
        }
    }

	// Declaring metafeatures.
    const Rcpp::StringVector gene_sym(Gs), gene_name(Gn), gene_type(T), gene_str(str);
    const int ngenes=gene_name.size();
	if (ngenes!=gene_type.size() || ngenes!=gene_sym.size() || ngenes!=gene_str.size()) { 
		throw std::runtime_error("gene data vectors should have the same length"); 
	}

    // Setting up structures for string organization.
    Rcpp::StringVector output(nin);

    typedef std::pair<int, int> overlap; // first=>subject index, second=>distance.
    std::deque<overlap> allindices;
    auto sorter = [&] (const overlap& left, const overlap& right) -> bool {
        return (gene_name[left.first] < gene_name[right.first]);
    };

    size_t counter=0;
	for (int curreg=0; curreg<nin; ++curreg) {
        allindices.clear();
        
        // Identifying all overlaps corresponding to the current region.
        // Overlaps should be sorted by query, so the loop below is valid.
        while (counter < nolaps && query[counter]==curreg) { 
            const int cur_sub=subject[counter];
    		if (cur_sub >= ngenes) { 
                throw std::runtime_error("symbol out of range for overlap index"); 
            }

            allindices.push_back(overlap(cur_sub, 0));
            if (use_dist) { 
                allindices.back().second = dist[counter];
            } 

            ++counter;
        }

		// Collapsing into a string.
		std::sort(allindices.begin(), allindices.end(), sorter);
        std::stringstream outstring;

        size_t innercount=0;
        bool empty=true;
        while (innercount < allindices.size()) { 
            bool has_prom=false, has_exon=false, has_gene=false;
            const auto ref=allindices[innercount].first;
            auto current_feature=gene_name[ref];
            int mindist=allindices[innercount].second;

            do {
                const char cur_type=Rcpp::String(gene_type[allindices[innercount].first]).get_cstring()[0];
                switch (cur_type) {
                    case 'E':
                        has_exon=true;
                        if (use_dist && mindist > allindices[innercount].second) {
                            mindist=allindices[innercount].second;
                        }
                        break;
                    case 'P':
                        has_prom=true;
                        break;
                    case 'G':
                        has_gene=true;
                        break;
                    default:
                        {
                            std::stringstream err;
                            err << "unrecognized feature type '" << cur_type << "'";
                            throw std::runtime_error(err.str().c_str());
                        }
                }
                ++innercount;
            } while (innercount < allindices.size() && gene_name[allindices[innercount].first]==current_feature);

            // Only reporting the subject if it's a full overlap or if it flanks an exon.
            if (!use_dist || has_exon) {
                if (!empty) { 
                    outstring << ",";
                }
                outstring << Rcpp::String(gene_sym[ref]).get_cstring() << ":" 
                    << Rcpp::String(gene_str[ref]).get_cstring() << ":";

                // Reporting the code for a direct overlap; reporting the distance to exon for flanking overlap.
                if (use_dist) {
                    outstring << mindist;
                } else {
                    if (has_prom) { 
                        outstring << "P";
                    }
                    if (has_exon) { 
                        outstring << "E";
                    }
                    if (has_gene && !has_exon) { // not reporting introns if it overlaps an exon.
                        outstring << "I";
                    }
                }

                empty=false;
            }

    		output[curreg] = outstring.str();
        }
	}

    return output;
    END_RCPP
}
