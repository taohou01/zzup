#ifndef _TOTAL_COMPLEX_H_
#define _TOTAL_COMPLEX_H_

#include <vector>
#include <unordered_map>

#include "common.h"
#include "chain.h"

class TotalComplex {
public:
    void init(Integer v_num, Integer dim);
    void enumComb(Simplex &simp, Integer d, Integer depth);
    void getCofaces(Integer s_id, Integer cof_dim, std::vector<Integer> &cofaces);
    void getBoundary(Integer s_id, Chain &bdry);

    std::shared_ptr<Chain> getBoundary(Integer s_id) {
        std::shared_ptr<Chain> bdry(new Chain);
        this->getBoundary(s_id, *bdry);
        return bdry;
    }

    Integer dimStart(Integer d) { return this->dim_start[d]; }
    Integer dimEnd(Integer d) { return this->dim_end[d]; }
    Integer dimSimpCnt(Integer d) { return this->dim_end[d] - this->dim_start[d]; }
    Integer edgeCnt() { return this->dimSimpCnt(1); }

    bool isSimpFaceOf(Integer s_id1, Integer s_id2) 
    { return isSortedVecSubsetOf<Integer>(this->simplices[s_id1], this->simplices[s_id2]); }

    Integer simpDim(Integer s_id) { return ::simpDim(this->simplices[s_id]); }

public:
    Integer v_num;
    std::vector<Simplex> simplices;
    std::vector<Integer> dim_start;
    std::vector<Integer> dim_end;
    SimplexIdMap simplex_id;
};

#endif
