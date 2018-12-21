#include "intersector.h"

intersector::intersector(SEXP p, SEXP e) : positions(p), elements(e), index(0), num_open(0), lastpos(0) {
    // See csaw:::.setupDiscard for the expected input vectors.
    // We do some fairly careful input checks here, just in case
    // users do stupid things to the readParam() object.
    const size_t N=positions.size();
    if (N!=elements.size()) {
        throw std::runtime_error("position and element vectors should be of the same length");
    }

    if (N) {
        if (positions[0]<1) {
            throw std::runtime_error("position vector should be 1-based");
        }

        for (size_t i=1; i<N; ++i) {
            if (positions[i]<positions[i-1]) {
                throw std::runtime_error("position vector should be sorted");
            }
        }
    }

    if (N%2!=0) {
        throw std::runtime_error("each element should be present exactly twice");
    }
    const size_t nelements=N/2;

    open.resize(nelements);
    for (size_t i=0; i<N; ++i) {
        auto current=elements[i];
        if (current < 0 || current >= nelements) {
            throw std::runtime_error("element ID out of range for blacklister");
        }
        ++open[current];
    }

    for (auto n : open) {
        if (n!=2) { throw std::runtime_error("each element should be present exactly twice"); }
    }
    std::fill(open.begin(), open.end(), 0);

    return;
}

void intersector::advance_to_start(int curpos) {
    if (lastpos > curpos) {
        throw std::runtime_error("supplied base positions should not decrease");
    }

    while (index < positions.size() && curpos >= positions[index]) {
        auto cur_element=elements[index];
        auto& cur_open=open[cur_element];
        cur_open=1-cur_open;

        if (cur_open) {
            ++num_open;
        } else {
            --num_open;
        }
        ++index;
    }

    lastpos=curpos;
    return;
}

bool intersector::end_is_within(int curend) const {
    if (lastpos > curend) {
        throw std::runtime_error("end of read should not occur before the start position");
    }

    int end_index=index;
    int cur_open=num_open;
    while (end_index < positions.size() && curend > positions[end_index]) {
        /* Only considering elements that are already open. Thus, there 
         * is no need to edit 'open', as there shouldn't be any more 
         * positions for this element once it's opened _and_ closed.
         */
        if (open[elements[end_index]]) { 
            --cur_open;
        }
        ++end_index;
    }
    return (cur_open > 0);
}
