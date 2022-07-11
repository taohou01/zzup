#include "total_complex.h"

#include <iostream>

#include "utils.h"

void TotalComplex::init(Integer v_num, Integer dim) {
    this->v_num = v_num;

    this->dim_start.resize(dim + 1);
    this->dim_end.resize(dim + 1);

    Simplex temp_simp;

    for (Integer d = 0; d <= dim; ++d) {
        this->dim_start[d] = this->simplices.size();

        temp_simp.clear();
        this->enumComb(temp_simp, d, -1);

        this->dim_end[d] = this->simplices.size();
    }
}

void TotalComplex::enumComb(Simplex &simp, Integer d, Integer depth) {
    // std::cout << d << " " << depth << std::endl;
    if (depth == d) {
        assert(std::is_sorted(simp.begin(), simp.end()));
        this->simplex_id.insert({simp, this->simplices.size()});
        this->simplices.push_back(simp);
        return;
    }

    Integer enum_start = 0;
    if (depth >= 0) {
        enum_start = simp[depth] + 1;
    }

    for (Integer i = enum_start; i + std::max(0,d-depth-1) < this->v_num; ++i) {
        simp.push_back(i);
        this->enumComb(simp, d, depth+1);
        simp.pop_back();
    }
}

void TotalComplex::getCofaces(Integer s_id, 
    Integer cof_dim, std::vector<Integer> &cofaces) {

    assert(cof_dim >= this->simplices[s_id].size());

    auto codim = cof_dim - this->simplices[s_id].size();

    Simplex cof;
    for (auto i = this->dimStart(codim); i < this->dimEnd(codim); ++i) {
        cof.clear();
        cof.insert(cof.end(), this->simplices[s_id].begin(), this->simplices[s_id].end());
        cof.insert(cof.end(), this->simplices[i].begin(), this->simplices[i].end());

        std::sort(cof.begin(), cof.end());

        Integer j;
        for (j = 1; j < cof.size() && cof[j-1] != cof[j]; ++j) ;
        if (j != cof.size()) { continue; }

        cofaces.push_back(this->simplex_id[cof]);
    }
}

void TotalComplex::getBoundary(Integer s_id, Chain &bdry) {
    bdry.clear();

    if (s_id < this->dimEnd(0)) { return; }

    Simplex &simp = this->simplices[s_id];
    Simplex bound_simp(simp.size()-1);

    for (auto i = 0; i < simp.size(); ++i) {
        bound_simp.clear();
        bound_simp.insert(bound_simp.end(), simp.begin(), simp.begin()+i);
        bound_simp.insert(bound_simp.end(), simp.begin()+i+1, simp.end());
        assert(bound_simp.size() == simp.size() - 1);
        bdry.push_back(this->simplex_id[bound_simp]);
    }

    std::sort(bdry.begin(), bdry.end());
}
