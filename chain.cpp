#include "chain.h"

void sumChain(const Chain &c1, const Chain &c2, Chain &c3) {
    c3.clear();

    Integer i = 0, j = 0;

    while (i < c1.size() && j < c2.size()) {
        if (c1.at(i) < c2.at(j)) 
        { c3.push_back(c1.at(i)); i ++; } 
        else if (c1.at(i) > c2.at(j)) 
        { c3.push_back(c2.at(j)); j ++; } 
        else 
        { i ++; j ++; }
    }

    for (; i < c1.size(); i ++) { c3.push_back(c1.at(i)); }

    for (; j < c2.size(); j ++) { c3.push_back(c2.at(j)); }
}
