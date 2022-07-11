#ifndef _DYNAMIC_ZIGZAG_H_
#define _DYNAMIC_ZIGZAG_H_

#include "common.h"

#include <vector>
#include <unordered_map>
#include <memory>

#include "chain.h"
#include "total_complex.h"

class RepInterval {
public:
    RepInterval() : b(-1), d(-1), dimension(-1) {}

    Integer dim() { return dimension; }
    Integer len() { return d-b+1; }

    std::shared_ptr<Chain>& Z(Integer i) { return this->cycles.at(i-this->b); }
    std::shared_ptr<Chain>& C(Integer i) { return this->chains.at(i-this->b+1); }


    /* Ones that read global variables 'reversed_' and 'dzz_' */

    Integer birth();
    Integer death();
    void setBirth(Integer _b);
    void setDeath(Integer _d);

    std::shared_ptr<Chain>& cyc(Integer i);
    std::shared_ptr<Chain>& chn(Integer i);

    void append();
    void prepend();
    void eraseBegin();
    void eraseEnd();

public:
    Integer b, d;
    Integer dimension;
    std::vector< std::shared_ptr<Chain> > cycles;
    std::vector< std::shared_ptr<Chain> > chains;
};

class DynamicZigzag {
public:
    DynamicZigzag(TotalComplex* _tot_cplx) : new_intv_id(0), tot_cplx(_tot_cplx) {}

    Integer m() { return this->filt.size(); }

    Integer addInterval(std::shared_ptr<RepInterval> rep);

    bool birthLessThan(Integer b1, Integer b2);

    bool deathLessThan(Integer d1, Integer d2);

    // Switch sigma_{i-1} and sigma_i
    void forwardSwitch(Integer i);

    // Switch sigma_{i-1} and sigma_i
    void backwardSwitch(Integer i);

    void outwardSwitch(
        Integer i,
        std::vector<Integer> &vanish_intvs, 
        std::vector<Integer> &new_intvs);

    void inwardContraction(
        Integer i, 
        std::vector<Integer> &vanish_intvs, 
        std::vector<Integer> &new_intvs);

    // Insert the two arrows between K_{i-1} and K_i
    void outwardExpansion(
        Integer i, 
        Integer sigma,
        std::vector<Integer> &vanish_intvs, 
        std::vector<Integer> &new_intvs);

    void getPersistence(
        std::vector< std::tuple<Integer, Integer, Integer> > &persistence);

    void printPers(const std::string outfname);

private:
    SimpOp filtOp(Integer i);

    void forwardSwitchImpl(Integer i);

    void forwardSwitchUpdate(Integer i, Integer sigma,
        const std::shared_ptr<RepInterval> &rep);

    void forwardSwitchCaseA(Integer i, Integer rid_b_i, Integer rid_b_iP1);

    void forwardSwitchCaseAUpdate(Integer i, 
        const std::shared_ptr<RepInterval> &rep_b_i, 
        const std::shared_ptr<RepInterval> &rep_b_iP1);

    void forwardSwitchCaseB(Integer i, Integer rid_d_iM1, Integer rid_d_i);

    void forwardSwitchCaseC(Integer i, Integer rid_b_i, Integer rid_d_i);

    void forwardSwitchCaseD(Integer i, Integer rid_d_iM1, Integer rid_b_iP1);

    void inwardContractionInjective(
        Integer i, 
        std::vector<Integer> &vanish_intvs, 
        std::vector<Integer> &new_intvs, 
        Integer rid_b_d_i);

    void inwardContractionSurjective(
        Integer i, 
        std::vector<Integer> &vanish_intvs, 
        std::vector<Integer> &new_intvs);

    void inwardContractionSqueezeRep(
        Integer i, std::shared_ptr<RepInterval> rep);

    void inwardContractionSqueezeRep(
        Integer i, 
        std::shared_ptr<RepInterval> rep, 
        std::shared_ptr<Chain> C_bar);

    void outwardExpansionInjective(
        Integer i, 
        Integer sigma,
        const std::vector<Integer> &sigma_rel_rids,
        std::vector<Integer> &vanish_intvs, 
        std::vector<Integer> &new_intvs);

    void outwardExpansionStretchRep(
        Integer i, Integer sigma, std::shared_ptr<RepInterval> rep);

    std::shared_ptr<RepInterval> sumRep(const std::shared_ptr<RepInterval> r1, 
        const std::shared_ptr<RepInterval> r2);

    void copyRep(const std::shared_ptr<RepInterval> from_r, 
        Integer start, Integer cnt, std::shared_ptr<RepInterval> to_r);

    void sumRepBasic(
        const std::shared_ptr<RepInterval> r1, Integer r1_start, Integer cnt,
        const std::shared_ptr<RepInterval> r2, Integer r2_start,
        std::shared_ptr<RepInterval> r3);

    bool repValid(const std::shared_ptr<RepInterval> r);

public:
    std::vector<SimpOp> filt;
    std::unordered_map< Integer, std::shared_ptr<RepInterval> > intervals;

private:
    Integer new_intv_id;
    TotalComplex* tot_cplx;
};

extern bool reversed_;
extern DynamicZigzag* dzz_;

inline SimpOp DynamicZigzag::filtOp(Integer i) {
    if (reversed_) {
        SimpOp s_op = this->filt[this->m()-i-1];
        s_op.op = (AddDelOp)(DEL_OP - s_op.op);

        // if (s_op.op == ADD_OP) { s_op.op = DEL_OP; }
        // else { s_op.op = ADD_OP; }

        return s_op;
    }

    return this->filt[i];
}

inline Integer RepInterval::birth() {
    if (reversed_) { return dzz_->m() - this->d; }
    return this->b;
}

inline Integer RepInterval::death() {
    if (reversed_) { return dzz_->m() - this->b; }
    return this->d;
}

inline void RepInterval::setBirth(Integer _b) {
    if (reversed_) { this->d = dzz_->m() - _b; return; }
    this->b = _b;
}

inline void RepInterval::setDeath(Integer _d) {
    if (reversed_) { this->b = dzz_->m() - _d; return; }
    this->d = _d;
}

inline std::shared_ptr<Chain>& RepInterval::cyc(Integer i) {
    if (reversed_) { return this->cycles.at(dzz_->m()-i - this->b); }
    return this->cycles.at(i - this->b);
}

inline std::shared_ptr<Chain>& RepInterval::chn(Integer i) {
    if (reversed_) { return this->chains.at(dzz_->m()-i - this->b); }
    return this->chains.at(i - this->b + 1);
}

inline void RepInterval::append() {
    if (reversed_) {
        this->cycles.insert(this->cycles.begin(), nullptr);
        this->chains.insert(this->chains.begin()+1, nullptr);
    } else {
        this->cycles.push_back(nullptr);
        this->chains.insert(this->chains.end()-1, nullptr);
    }

    this->setDeath(this->death()+1);
    this->cyc(this->death()) = this->cyc(this->death()-1);
}

inline void RepInterval::prepend() {
    if (reversed_) {
        this->cycles.push_back(nullptr);
        this->chains.insert(this->chains.end()-1, nullptr);
    } else {
        this->cycles.insert(this->cycles.begin(), nullptr);
        this->chains.insert(this->chains.begin()+1, nullptr);
    }

    this->setBirth(this->birth()-1);
    this->cyc(this->birth()) = this->cyc(this->birth()+1);
}

inline void RepInterval::eraseBegin() {
    if (reversed_) {
        this->cycles.erase(this->cycles.end()-1);
        this->chains.erase(this->chains.end()-2);
    } else {
        this->cycles.erase(this->cycles.begin());
        this->chains.erase(this->chains.begin()+1);
    }

    this->setBirth(this->birth()+1);
}

inline void RepInterval::eraseEnd() {
    this->chn(this->death()) = 
        sumChain( this->chn(this->death()-1), this->chn(this->death()) );

    if (reversed_) {
        this->cycles.erase(this->cycles.begin());
        this->chains.erase(this->chains.begin()+1);
    } else {
        this->cycles.erase(this->cycles.end()-1);
        this->chains.erase(this->chains.end()-2);
    }

    this->setDeath(this->death()-1);
}

#endif
