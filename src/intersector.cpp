#include "intersector.h"

intersector::intersector(SEXP p, SEXP e) : positions(p), elements(e), index(0), open(positions.size()), num_open(0) {
    // See csaw:::.setupDiscard for the expected input vectors.
    if (positions.size()!=elements.size()) {
        throw std::runtime_error("position and identity vectors should be of the same length");
    }
    return;
}

void intersector::advance_to_start(int curpos) {
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
}

bool intersector::end_is_within(int curend) const {
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
