#include "bam_utils.h"
#include "utils.h"

BamFile::BamFile(SEXP bam, SEXP idx) {
    const char* path=check_string(bam, "BAM file path");
    const char* xpath=check_string(idx, "BAM index file path");

    in = sam_open(path, "rb");
    if (in == NULL) {
        std::stringstream err;
        err << "failed to open BAM file at '" << path << "'";
        throw std::runtime_error(err.str());
    }
    try {
        index = bam_index_load(xpath); 
        if (index==NULL) { 
            std::stringstream err;
            err << "failed to open BAM index at '" << xpath << "'";
            throw std::runtime_error(err.str());
        }
        try {
            header=sam_hdr_read(in);
        } catch (std::exception& e) {
            hts_idx_destroy(index);
            throw;
        }
    } catch (std::exception& e) {
        sam_close(in);
        throw;
    }
    bgzf_set_cache_size(in->fp.bgzf, 100*BGZF_MAX_BLOCK_SIZE);
    return;
}

BamFile::~BamFile() {
    bam_hdr_destroy(header);
    hts_idx_destroy(index); 
    sam_close(in);
    return;
}

BamRead::BamRead() {
    read=bam_init1();
    return;
}

bool BamRead::is_well_mapped(const int& minqual, const bool& rmdup) const {
    if (minqual!=NA_INTEGER && (read -> core).qual < minqual) { return false; } 
    if (rmdup && ((read -> core).flag & BAM_FDUP)!=0) { return false; }
    return true;
}

void BamRead::extract_data(AlignData& data) const {
    data.len=bam_cigar2rlen((read->core).n_cigar, bam_get_cigar(read));
    data.is_reverse=bam_is_rev(read);
    return;
}

BamRead::~BamRead() { 
    bam_destroy1(read);
    return;
}

AlignData::AlignData() : len(0), is_reverse(true) { }

BamIterator::BamIterator(const BamFile& bf) : iter(NULL) {
    iter=bam_itr_queryi(bf.index, HTS_IDX_NOCOOR, 0, 0);
    return;
}

BamIterator::BamIterator(const BamFile& bf, SEXP Chr, SEXP Start, SEXP End) : iter(NULL) {
    // Checks on inputs. Get start to 0-indexed closed, end to 0-indexed open.
    const char* chr=check_string(Chr, "chromosome name");
    int start=check_integer_scalar(Start, "start position")-1;
    int end=check_integer_scalar(End, "end position");

    // Pulling out chromsoome name.
    int cid=bam_name2id(bf.header, chr);
    if (cid==-1) {
        std::stringstream err;
        err << "reference sequence '" << chr << "' missing in BAM header";
        throw std::runtime_error(err.str());
    }

    // Checks on coordinates.
    if (start < 0) { start=0; }
    const int curlen = (bf.header->target_len)[cid];
    if (end > curlen) { end = curlen; }
    if (start > end) {
        throw std::runtime_error("invalid values for region start/end coordinates");
    }

    // Constructs the iterator.
    iter=bam_itr_queryi(bf.index, cid, start, end);
    return;
}

BamIterator::BamIterator(const BamFile& bf, int cid, int start, int end) : iter(bam_itr_queryi(bf.index, cid, start, end)) {}

BamIterator::~BamIterator() { 
    bam_itr_destroy(iter); 
    return;
}

