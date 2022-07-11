#ifndef _UD_H_
#define _UD_H_

#include <vector>
#include <memory>

#include "common.h"
#include "chain.h"
#include "total_complex.h"
#include "dynamic_zigzag.h"

class AliveInterval {
public:
    Integer birth;
    // Integer death;
    std::shared_ptr<Chain> left_cyc;
    std::shared_ptr<Chain> chain;
    std::shared_ptr<Chain> right_cyc;
};

class UpDownPersistence {
public:
    UpDownPersistence(TotalComplex *_tot_cplx, DynamicZigzag *_dzz) { 
        this->tot_cplx = _tot_cplx; 

        this->dzz = _dzz; 
        
        this->simp_cnt = 0;

        this->simp_ext_id.resize(this->tot_cplx->simplices.size(), -1);
        this->simp_id_by_ext_id.resize(this->tot_cplx->simplices.size(), -1);

        this->pair.resize(this->tot_cplx->simplices.size(), UNPAIRED_IND);
        // this->pivot_map.resize(this->tot_cplx->simplices.size(), -1);
        this->cycles.resize(this->tot_cplx->simplices.size());
        this->chains.resize(this->tot_cplx->simplices.size());
    }

    void addSimplex(Integer ext_s_id);

    void middle();

    void deleteSimplex(Integer ext_s_id);

private:
    std::shared_ptr<Chain> getExtChain(std::shared_ptr<Chain> chn);

    void addUpOrDownInterval(
        Integer dim, Integer birth, Integer death, 
        std::shared_ptr<Chain> cyc, 
        std::shared_ptr<Chain> chn);

    void addClosedClosedInterval(
        Integer dim, Integer birth, Integer death, 
        std::shared_ptr<Chain> left_cyc, std::shared_ptr<Chain> chn, 
        std::shared_ptr<Chain> right_cyc);

private:
    TotalComplex *tot_cplx;

    DynamicZigzag *dzz;

    Integer simp_cnt, cur_index;

    std::vector<Integer> simp_ext_id, simp_id_by_ext_id;

    std::vector<Integer> pair;
    // std::vector<Integer> pivot_map;
    std::vector<std::shared_ptr<Chain>> cycles;
    std::vector<std::shared_ptr<Chain>> chains;

    std::vector<AliveInterval> alives;

    const Integer NEG_IND = -1;
    const Integer UNPAIRED_IND = -2;
};

#endif
