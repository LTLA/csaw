#ifndef BAM_UTILS_H
#define BAM_UTILS_H

#include "csaw.h"
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/bgzf.h"

struct BamFile {
    BamFile(SEXP, SEXP);
    BamFile(const BamFile&)=delete;
    BamFile& operator=(const BamFile&)=delete;
    BamFile(BamFile&&)=delete;
    BamFile& operator=(BamFile&&)=delete;
    ~BamFile();

    samFile* in;
    hts_idx_t * index;
    bam_hdr_t * header;
};
   
struct BamRead {
    BamRead();
    BamRead(const BamRead&)=delete;
    BamRead& operator=(const BamRead&)=delete;
    BamRead(BamRead&&)=delete;
    BamRead& operator=(BamRead&&)=delete;
    ~BamRead();

    uint16_t get_flag() const;
    int32_t get_aln_pos() const;
    bool is_well_mapped(int, bool) const;
    bool is_reverse() const;
    int get_aln_len() const;

    bam1_t* read;
};

struct BamIterator {
    BamIterator(const BamFile&);
    BamIterator(const BamFile&, SEXP, SEXP, SEXP);
    BamIterator(const BamFile&, int, int, int);

    BamIterator(const BamIterator&)=delete;
    BamIterator& operator=(const BamIterator&)=delete;
    BamIterator(BamIterator&&)=delete;
    BamIterator& operator=(BamIterator&&)=delete;
     ~BamIterator();

    hts_itr_t* iter;
};

#endif
