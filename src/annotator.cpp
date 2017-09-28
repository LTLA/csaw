#include "csaw.h"
#include "utils.h"

/* This function spits out the ID for each exon, with some degree of strand-awareness, 
 * such that the first exon in the gene is labelled as exon 1, then 2, 3, etc. They
 * are assumed to be stored in genomic order so we just label them as-is.
 */

SEXP collate_exon_data (SEXP geneid, SEXP strand, SEXP start, SEXP end) {
    BEGIN_RCPP

	// Checking inputs.
    const Rcpp::IntegerVector _geneid(geneid), _start(start), _end(end);
    const Rcpp::LogicalVector _strand(strand);
    const int n=_geneid.size();
    if (n!=_start.size() || n!=_end.size() || n!=_strand.size()) {
        throw std::runtime_error("exon data vectors should have the same length");
    }
	
	// Scanning through to determine the number of unique genes.
	int nuniq=0;
	if (n > 0) {
		nuniq=1;
		for (int x=1; x<n; ++x) {
			if (_geneid[x]!=_geneid[x-1]) {
				++nuniq;
			} else if (_strand[x]!=_strand[x-1]) {
				throw std::runtime_error("exons of the same gene should have the same strand");
			} else if (_start[x]<_start[x-1]) {
				throw std::runtime_error("exons of the same gene should be sorted by the start index");
			}
		}
	}

	// Setting up output structures.
    Rcpp::IntegerVector out_exnum(n), out_index(nuniq), out_start(nuniq), out_end(nuniq);
    int curex=0, ngene=0;
	std::deque<std::pair<int, int> > current_indices;
	while (curex<n) {
		const int& current=_geneid[curex];
	    out_index[ngene]=curex+1;
	    out_start[ngene]=_start[curex]; // Input should be sorted by start positions within each gene ID, so first element should be earliest.
		int& last_end=(out_end[ngene]=_end[curex]);
		++ngene;

		// Resorting by end positions, if the gene is on the negative strand.
		if (!_strand[curex]) { 
			do {
				current_indices.push_back(std::make_pair(_end[curex], curex));
				++curex;
			} while (curex < n && current==_geneid[curex]);
			
			const int& stretch=current_indices.size();
			std::sort(current_indices.begin(), current_indices.end()); 
			for (int counter=0; counter<stretch; ++counter) { 
				out_exnum[current_indices[counter].second]=stretch - counter;
			}
			last_end=current_indices.back().first;
			current_indices.clear();
		} else {
			int counter=1;
			do {
				out_exnum[curex]=counter;
				if (last_end < _end[curex]) { last_end=_end[curex]; }
				++counter;
				++curex;
			} while (curex < n && current==_geneid[curex]);
		}
	}

	return Rcpp::List::create(out_exnum, Rcpp::List::create(out_index, out_start, out_end));
    END_RCPP
}

/* This function collapses indices into a string. The overlaps are also assigned
 * distances in 'dists', so 'dists[start <= x < end]' will give the distances from
 * the edge of the region to the annotated bit (if the pointer is not NULL).
 * The function will return a string collating all information for that annotated
 * feature (i.e., collate exon-level overlaps to a gene-level string).
 */

struct feature_data { 
    feature_data(int ID, int Feature, int Strand, int Index) : id(ID), feature(Feature), strand(Strand), index(Index), distance(0) {}
    int id, feature, strand, index;
    double distance;

    friend bool operator< (const feature_data& left, const feature_data& right) {
        if (left.id==right.id) {
            if (left.feature==right.feature) {
                return left.index < right.index;
            }
            return left.feature < right.feature;
        } 
        return left.id < right.id;
    }
};

std::string digest2string (const std::deque<feature_data>& gelements, const Rcpp::StringVector& symbols, bool use_dist) {
	if (!gelements.size()) { return ""; }
	std::stringstream ss;
	size_t start=0, end, index=0;

	while (start < gelements.size()) {
		if (start!=0) { ss << ","; }
		ss << Rcpp::as<std::string>(symbols[gelements[start].index]) << '|'; 
		end=start+1;	
		while (end < gelements.size() && gelements[end].id==gelements[start].id) { ++end; }

		// Deciding what to print.
		if (end==start+1) {
			if (gelements[start].feature==-1) {
				ss << "I"; 
			} else {
				ss << gelements[start].feature;
			}
		} else {	
			index=start;
            // Skipping gene bodies and extra promoters.
			if (gelements[start].feature==-1) { ++index; }
			ss << gelements[index].feature;
            while (index+1 < end && gelements[index+1].feature==0) { ++index; }
			bool wasokay=false;

			// Running through and printing all stretches of contiguous exons.
			while ((++index) < end) {
                if (gelements[index].feature==gelements[index-1].feature+1) { 
					wasokay=true;
				} else {
					if (wasokay) {
						ss << '-' << gelements[index-1].feature;
						wasokay=false;
					}
					ss << ',' << gelements[index].feature;
				}
			}
			if (wasokay) { ss << '-' << gelements[index-1].feature; }
		}
		
		// Adding the strand and distance information.	
		ss << '|' << (gelements[start].strand ? '+' : '-');
		if (use_dist) {
			int lowest=gelements[start].distance;
			for (index=start+1; index < end; ++index) {
				if (lowest > gelements[index].distance) { lowest=gelements[index].distance; }
			}
			ss << '[' << lowest << ']';
		}
		start=end;
	}
	return ss.str();
}

/* The main function */

SEXP annotate_overlaps (SEXP N, SEXP fullQ, SEXP fullS, SEXP leftQ, SEXP leftS, SEXP leftDist,
		SEXP rightQ, SEXP rightS, SEXP rightDist, 
		SEXP symbol, SEXP genefeature, SEXP geneid, SEXP genestr) {

    BEGIN_RCPP
    const int nin=check_integer_scalar(N, "number of query regions");

    // Setting up overlap information.
    const Rcpp::IntegerVector _fullQ(fullQ), _fullS(fullS), 
        _leftQ(leftQ), _leftS(leftS), _leftDist(leftDist),
        _rightQ(rightQ), _rightS(rightS), _rightDist(rightDist);
	const int nfull=_fullQ.size();
	if (nfull!=_fullS.size()){ 
        throw std::runtime_error("full overlap vectors should have equal length"); 
    }
    const int nleft=_leftQ.size();
	if (nleft!=_leftS.size() || nleft!=_leftDist.size()) { 
        throw std::runtime_error("left overlap vectors should have equal length"); 
    }
	const int nright=_rightQ.size();
	if (nright!=_rightS.size() || nright!=_rightDist.size()) { 
        throw std::runtime_error("right overlap vectors should have equal length"); 
    }

	// Declaring metafeatures.
    const Rcpp::StringVector _symbol(symbol);
    const Rcpp::IntegerVector _geneid(geneid), _genefeature(genefeature), _genestr(genestr);
    const int nsym=_symbol.size();
	if (nsym!=_geneid.size() || nsym!=_genefeature.size() || nsym!=_genestr.size()) {
		throw std::runtime_error("gene data vectors should have the same length"); 
	}

    // Going through them and assembling the output vectors. 
    Rcpp::StringVector out_full(nin), out_left(nin), out_right(nin);
    
	int fullx=0, leftx=0, rightx=0;
    int* curx_p;
    const int* curn_p;
    Rcpp::IntegerVector::const_iterator cur_qIt, cur_sIt, cur_dIt;
    Rcpp::StringVector::iterator cur_oIt;
    bool use_dist=false;
    std::deque<feature_data> allindices;

	for (int curreg=0; curreg<nin; ++curreg) {
		// Adding all overlaps of each type. Assuming that findOverlaps gives ordered output, which it should.
		for (int mode=0; mode<3; ++mode) {
			if(mode==0) {
				curx_p=&fullx;
				curn_p=&nfull;
				cur_qIt=_fullQ.begin();
				cur_sIt=_fullS.begin();
				cur_oIt=out_full.begin();
                use_dist=false;
			} else if (mode==1) {
				curx_p=&leftx;
				curn_p=&nleft;
				cur_qIt=_leftQ.begin();
				cur_sIt=_leftS.begin();
				cur_dIt=_leftDist.begin();
				cur_oIt=out_left.begin();
                use_dist=true;
			} else if (mode==2) {
				curx_p=&rightx;
				curn_p=&nright;
				cur_qIt=_rightQ.begin();
				cur_sIt=_rightS.begin();
				cur_dIt=_rightDist.begin();
				cur_oIt=out_right.begin();
                use_dist=true;
			}

			// For the current region, we get everything in the current overlap.
			int& curx=*curx_p;
			const int& curn=*curn_p;
			allindices.clear();
			while (curx < curn && *(cur_qIt+curx)==curreg) {
                const int& index=*(cur_sIt+curx);
				if (index >= nsym) { 
                    throw std::runtime_error("symbol out of range for overlap index"); 
                }
				allindices.push_back(feature_data(_geneid[index], _genefeature[index], _genestr[index], index));
				if (use_dist) { 
                    // Storing distance to each overlapped feature.
                    allindices.back().distance = *(cur_dIt+curx); 
                } 
				++curx;
			}

			// Sorting by gene index, then feature index; then collapsing into a string.
			std::sort(allindices.begin(), allindices.end());
			*(cur_oIt+curreg) = digest2string(allindices, _symbol, use_dist);
		} 
	}

    return Rcpp::List::create(out_full, out_left, out_right);
    END_RCPP
}
