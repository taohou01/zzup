#ifndef _CHAIN_H_
#define _CHAIN_H_

#include <vector>
#include <memory>
#include <algorithm>

#include "common.h"

typedef std::vector<Integer> Chain;

inline bool isEmptyChain(const std::shared_ptr<Chain> &c) {
    return c == nullptr || c->size() == 0;
}

void sumChain(const Chain &c1, const Chain &c2, Chain &c3);

// inline void sumChain(const std::shared_ptr<Chain> &c1, 
//     const std::shared_ptr<Chain> &c2, std::shared_ptr<Chain> &c3) { 

//     if (isEmptyChain(c1) && isEmptyChain(c2)) { c3 = nullptr; return; } 

//     if (!c3) { c3.reset(new Chain()); }

//     if (isEmptyChain(c1)) { *c3 = *c2; return; }
//     if (isEmptyChain(c2)) { *c3 = *c1; return; }

//     sumChain(*c1, *c2, *c3);
// }

// inline void sumChain(const std::shared_ptr<Chain> &c1, 
//     const std::shared_ptr<Chain> &c2, Chain &c3) { 

//     c3.clear();

//     if (c1 == nullptr) {
//         if (c2 == nullptr) { return; } 
//         else { c3 = *c2; return; }
//     } else if (c2 == nullptr) { c3 = *c1; return; }

//     sumChain(*c1, *c2, c3); 
// }

// inline void sumChain(const std::shared_ptr<Chain> &c1, std::shared_ptr<Chain> &c2) {
//     if (c1 == nullptr) { return; }
//     if (c2 == nullptr) { c2.reset(new Chain(*c1)); return; }

//     Chain c3;
//     sumChain(*c1, *c2, c3);
//     *c2 = c3;
// }

inline std::shared_ptr<Chain> sumChain(
    const std::shared_ptr<Chain> c1, const std::shared_ptr<Chain> c2) {

    if (isEmptyChain(c1)) {
        if (isEmptyChain(c2)) { return nullptr; } 
        else { return c2; }
        // else { return std::shared_ptr<Chain>(new Chain(*c2)); }
    } else if (isEmptyChain(c2)) { 
        return c1;
        // return std::shared_ptr<Chain>(new Chain(*c1));
    }

    std::shared_ptr<Chain> c3(new Chain);
    sumChain(*c1, *c2, *c3);

    return c3;
}

inline std::shared_ptr<Chain> copyChainPtr(const std::shared_ptr<Chain> c) {
    return c;
    // if (isEmptyChain(c)) { return nullptr; } 
    // else { return std::shared_ptr<Chain>(new Chain(*c)); }
}

inline bool isSimpInChain(const Integer s_id, const Chain &c) {
    return std::binary_search(c.begin(), c.end(), s_id);
}

inline bool isSimpInChain(const Integer s_id, const std::shared_ptr<Chain> &c) {
    if (c) { return std::binary_search(c->begin(), c->end(), s_id); } 
    return false;
}

#endif
