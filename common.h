#ifndef _COMMON_H_
#define _COMMON_H_

#include <cmath>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>

#include <boost/functional/hash.hpp>

typedef long double Decimal;
typedef int Integer;

typedef std::vector<Integer> Simplex;

class Point {
public:
    Decimal x, y;
};

enum DpcEvtType {
    L_LOCAL_MIN,
    R_LOCAL_MIN,
    L_LOCAL_MAX,
    R_LOCAL_MAX,
    LOCAL_MIN,
    LOCAL_MAX,
    INC_CROSS,
    DEC_CROSS,
    OPPO_CROSS
};

class DpcEvent {
public:
    Decimal d;
    Decimal t;
    DpcEvtType type;
    Integer e1;
    Integer e2;
};

enum AddDelOp {
    ADD_OP = 0,
    DEL_OP = 1
};

class EdgeOp {
public:
    AddDelOp op;
    Integer e_id;
};

class SimpOp {
public:
    AddDelOp op;
    Integer s_id;
};

std::ostream& operator<<(std::ostream& os, const DpcEvent &evt);

inline Decimal eucDis(const Point &p1, const Point &p2) {
    return std::sqrt((p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y));
}

inline void lineIntersec(Decimal x1, Decimal y1, 
    Decimal x2, Decimal y2, 
    Decimal x3, Decimal y3, 
    Decimal x4, Decimal y4, 
    Decimal &r_x, Decimal &r_y) {

    Decimal denom = (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4);

    r_x = ((x1*y2 - y1*x2) * (x3-x4) - (x1-x2) * (x3*y4 - y3*x4)) / denom;
    r_y = ((x1*y2 - y1*x2) * (y3-y4) - (y1-y2) * (x3*y4 - y3*x4)) / denom;
}

inline bool isDecimalInteger(const Decimal x) {
    return x == std::floor(x);
}

inline Decimal interpVal(Decimal x1, Decimal y1, 
    Decimal x2, Decimal y2, Decimal x) {

    assert(x1 <= x && x <= x2);

    return (x-x1) / (x2-x1) * (y2-y1) + y1;
}

Integer combChoose(Integer n, Integer k);

template <class ElemType>
class VecHash { 
public:
    size_t operator()(const std::vector<ElemType>& v) const; 
};

template <class ElemType>
size_t VecHash<ElemType>
    ::operator()(const std::vector<ElemType>& v) const {

    std::size_t seed = 0;

    for (auto e : v) { boost::hash_combine(seed, e); }

    return seed;
}

template <class ElemType>
class VecEqual { 
public:
    bool operator()(const std::vector<ElemType>& v1, 
        const std::vector<ElemType>& v2) const; 
};

template <class ElemType>
bool VecEqual<ElemType>
    ::operator()(const std::vector<ElemType>& v1, 
        const std::vector<ElemType>& v2) const {

    if (v1.size() != v2.size()) { return false; }

    for (auto i = 0; i < v1.size(); i ++) {
        if (v1[i] != v2[i]) {
            return false;
        }
    }

    return true;
}

template <class T>
bool isSortedVecSubsetOf(const std::vector<T>& v1, const std::vector<T>& v2) {
    for (const auto& e : v1) {
        if ( !std::binary_search(v2.begin(), v2.end(), e) )
            { return false; }
    }

    return true;
}

typedef std::unordered_map< Simplex, Integer,
    VecHash<Integer>, VecEqual<Integer> > SimplexIdMap;

inline Integer simpDim(const Simplex &s) { return s.size()-1; }

#endif
