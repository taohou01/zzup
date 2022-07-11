#include "dynamic_zigzag.h"

#include <algorithm>
#include <forward_list>
#include <fstream>
#include <tuple>

#include "utils.h"

bool reversed_ = false;
DynamicZigzag* dzz_ = nullptr;

struct RepIdCmp {
    DynamicZigzag* dzz;

    bool operator()(Integer rid1, Integer rid2) const {
        auto b1 = dzz->intervals[rid1]->b;
        auto b2 = dzz->intervals[rid2]->b;
        return dzz->birthLessThan(b1, b2); 
    }
};

bool DynamicZigzag::birthLessThan(Integer b1, Integer b2) {
    assert(b1 != b2);

    if (b1 < b2) {
        if (this->filtOp(b2-1).op == ADD_OP) { return true; }
        else { return false; }
    } else {
        if (this->filtOp(b1-1).op == DEL_OP) { return true; }
        else { return false; }
    }
}

bool DynamicZigzag::deathLessThan(Integer d1, Integer d2) {
    assert(d1 != d2);

    if (d1 > d2) {
        if (this->filtOp(d2).op == DEL_OP) { return true; }
        else { return false; }
    } else {
        if (this->filtOp(d1).op == ADD_OP) { return true; }
        else { return false; }
    }
}

void DynamicZigzag::forwardSwitch(Integer i) {
    dzz_ = this;
    this->forwardSwitchImpl(i);
    std::swap(this->filt[i-1], this->filt[i]);
}

void DynamicZigzag::backwardSwitch(Integer i) {
    dzz_ = this;
    reversed_ = true;
    this->forwardSwitchImpl(this->m()-i);
    reversed_ = false;
    std::swap(this->filt[i-1], this->filt[i]);
}

void DynamicZigzag::forwardSwitchImpl(Integer i) {
    assert(i > 0 && i < this->m());
    assert(this->filtOp(i-1).op == ADD_OP && this->filtOp(i).op == ADD_OP);
    assert( !this->tot_cplx->isSimpFaceOf(this->filtOp(i-1).s_id, this->filtOp(i).s_id) );

    Integer rid_b_i = -1, rid_d_iM1 = -1, rid_b_iP1 = -1, rid_d_i = -1;

    for (const auto& it : this->intervals) {
        if (it.second->birth() == i) 
            { assert(rid_b_i < 0); rid_b_i = it.first; }
        else if (it.second->birth() == i+1) 
            { assert(rid_b_iP1 < 0); rid_b_iP1 = it.first; }

        if (it.second->death() == i-1) 
            { assert(rid_d_iM1 < 0); rid_d_iM1 = it.first; }
        else if (it.second->death() == i) 
            { assert(rid_d_i < 0); rid_d_i = it.first; }
    }

    assert((rid_b_i >= 0) != (rid_d_iM1 >= 0));
    assert((rid_b_iP1 >= 0) != (rid_d_i >= 0));

    if (rid_b_i >= 0 && rid_b_iP1 >= 0) {
        // std::cout << "  Case A" << std::endl;
        this->forwardSwitchCaseA(i, rid_b_i, rid_b_iP1);
    } else if (rid_d_iM1 >= 0 && rid_d_i >= 0) {
        // std::cout << "  Case B" << std::endl;
        this->forwardSwitchCaseB(i, rid_d_iM1, rid_d_i);
    } else if (rid_b_i >= 0 && rid_d_i >= 0) {
        // std::cout << "  Case C" << std::endl;
        this->forwardSwitchCaseC(i, rid_b_i, rid_d_i);
    } else if (rid_d_iM1 >= 0 && rid_b_iP1 >= 0) {
        // std::cout << "  Case D" << std::endl;
        this->forwardSwitchCaseD(i, rid_d_iM1 , rid_b_iP1);
    } else {
        assert(0);
    }
}

void DynamicZigzag::forwardSwitchUpdate(
    Integer i, Integer sigma,
    const std::shared_ptr<RepInterval> &rep) {

    if (rep->birth() <= i && i <= rep->death()) {
        assert(rep->birth() <= i-1 && rep->death() >= i+1);
        
        if ( isSimpInChain(sigma, rep->chn(i-1)) || 
            isSimpInChain(sigma, rep->cyc(i)) ) {

            auto new_chn = sumChain(rep->chn(i-1), rep->chn(i));

            rep->chn(i-1) = nullptr;
            rep->chn(i) = new_chn;
            rep->cyc(i) = rep->cyc(i-1);
        }
    } else {
        assert(rep->birth() > i+1 || rep->death() < i-1);
    }   
}

void DynamicZigzag::forwardSwitchCaseA(Integer i, Integer rid_b_i, Integer rid_b_iP1) {
    const auto sigma = this->filtOp(i-1).s_id;
    const auto tau = this->filtOp(i).s_id;

    for (const auto& it : this->intervals) {
        const auto& rep = it.second;
        if (rep->birth() != i && rep->birth() != i+1)
            { this->forwardSwitchUpdate(i, sigma, rep); }
    }

    auto& rep_b_i = this->intervals[rid_b_i];
    auto& rep_b_iP1 = this->intervals[rid_b_iP1];

    assert( isSimpInChain(sigma, rep_b_i->cyc(i)) );
    assert( isSimpInChain(tau, rep_b_iP1->cyc(i+1)) );

    if ( !isSimpInChain(sigma, rep_b_iP1->cyc(i+1)) ) {
        this->forwardSwitchCaseAUpdate(i, rep_b_i, rep_b_iP1);
    } else {
        assert(rep_b_i->dim() == rep_b_iP1->dim());

        const auto d1 = rep_b_i->death();
        const auto d2 = rep_b_iP1->death();

        if (this->deathLessThan(d1, d2)) {
            auto new_rep = this->sumRep(rep_b_i, rep_b_iP1);

            assert(rep_b_iP1->birth() == new_rep->birth() &&
                rep_b_iP1->death() == new_rep->death());
            assert( !isSimpInChain(sigma, new_rep->cyc(i+1)) );

            rep_b_iP1 = new_rep;
            this->forwardSwitchCaseAUpdate(i, rep_b_i, rep_b_iP1);
        } else {
            auto new_rep = this->sumRep(rep_b_i, rep_b_iP1);

            assert(new_rep->birth() == i+1 && new_rep->death() == d1);
            assert( !isSimpInChain(sigma, new_rep->cyc(i+1)) );

            new_rep->prepend();
            // new_rep->setBirth(i);
            new_rep->cyc(i) = new_rep->cyc(i+1);

            rep_b_i = new_rep;
        }
    }
}

void DynamicZigzag::forwardSwitchCaseAUpdate(Integer i, 
    const std::shared_ptr<RepInterval> &rep_b_i, 
    const std::shared_ptr<RepInterval> &rep_b_iP1) {

    rep_b_i->eraseBegin();
    rep_b_iP1->prepend();
}

void DynamicZigzag::forwardSwitchCaseB(Integer i, Integer rid_d_iM1, Integer rid_d_i) {
    const auto sigma = this->filtOp(i-1).s_id;
    const auto tau = this->filtOp(i).s_id;

    for (const auto& it : this->intervals) {
        const auto& rep = it.second;
        if (rep->death() != i-1 && rep->death() != i) 
            { this->forwardSwitchUpdate(i, sigma, rep); }
    }

    auto& zeta1 = this->intervals[rid_d_iM1];
    auto& zeta2 = this->intervals[rid_d_i];

    assert( isSimpInChain(sigma, zeta1->chn(i-1)) );
    assert( isSimpInChain(tau, zeta2->chn(i)) );

    auto Cp_iM1_sum_Cp_i = sumChain(zeta2->chn(i-1), zeta2->chn(i));

    if ( !isSimpInChain(sigma, Cp_iM1_sum_Cp_i) ) {
        zeta1->append();
        assert(zeta1->death() == i);
        // zeta1->setDeath(i);
        // zeta1->cyc(i) = zeta1->cyc(i-1);

        zeta2->eraseEnd();
        assert(zeta2->death() == i-1);
        // zeta2->setDeath(i-1);
        // zeta2->cyc(i-1) = Cp_iM1_sum_Cp_i;
    } else {
        assert(zeta1->dim() == zeta2->dim());
        const auto b1 = zeta1->birth();
        const auto b2 = zeta2->birth();

        if (this->birthLessThan(b1, b2)) {
            auto new_rep = sumRep(zeta1, zeta2);
            assert(zeta2->birth() == new_rep->birth());
            new_rep->eraseEnd();
            zeta2 = new_rep;
            assert(zeta2->death() == i-1);
            assert( !isSimpInChain(sigma, zeta2->chn(i-1)) );
            assert( isSimpInChain(tau, zeta2->chn(i-1)) );

            zeta1->append();
            assert(zeta1->death() == i);
        } else {
            zeta2->chn(i) = sumChain(zeta2->chn(i-1), zeta2->chn(i));
            zeta2->chn(i-1) = nullptr;
            zeta2->cyc(i) = zeta2->cyc(i-1);

            auto new_rep = sumRep(zeta1, zeta2);
            assert(zeta1->birth() == new_rep->birth());
            new_rep->eraseEnd();
            zeta1 = new_rep;
            assert(zeta1->death() == i-1);
            assert( !isSimpInChain(sigma, zeta1->chn(i-1)) );
            assert( isSimpInChain(tau, zeta1->chn(i-1)) );
        }
    }
}

void DynamicZigzag::forwardSwitchCaseC(Integer i, Integer rid_b_i, Integer rid_d_i) {
    assert(rid_b_i != rid_d_i);

    const auto sigma = this->filtOp(i-1).s_id;
    const auto tau = this->filtOp(i).s_id;

    for (const auto& it : this->intervals) {
        const auto& rep = it.second;
        if (rep->birth() != i && rep->death() != i)
            { this->forwardSwitchUpdate(i, sigma, rep); }
    }

    auto& zeta1 = this->intervals[rid_d_i];
    auto& zeta2 = this->intervals[rid_b_i];

    // assert( isSimpInChain(sigma, zeta1->chn(i-1)) );
    // assert( isSimpInChain(tau, zeta2->chn(i)) );

    zeta1->eraseEnd();
    assert(zeta1->death() == i-1);
    if ( isSimpInChain(sigma, zeta1->chn(i-1)) )
    { zeta1->chn(i-1) = sumChain(zeta2->cyc(i), zeta1->chn(i-1)); }
    assert( !isSimpInChain(sigma, zeta1->chn(i-1)) );

    zeta2->eraseBegin();
    assert(zeta2->birth() == i+1);
}

void DynamicZigzag::forwardSwitchCaseD(Integer i, Integer rid_d_iM1, Integer rid_b_iP1) {
    const auto sigma = this->filtOp(i-1).s_id;
    const auto tau = this->filtOp(i).s_id;

    for (const auto& it : this->intervals) {
        const auto& rep = it.second;
        if (rep->birth() != i+1 && rep->death() != i-1)
            { this->forwardSwitchUpdate(i, sigma, rep); }
    }

    auto& zeta1 = this->intervals[rid_d_iM1];
    auto& zeta2 = this->intervals[rid_b_iP1];

    if ( isSimpInChain(sigma, zeta2->cyc(i+1)) ) {
        zeta1->chn(i-1) = sumChain(zeta2->cyc(i+1), zeta1->chn(i-1));
        assert( !isSimpInChain(sigma, zeta1->chn(i-1)) );
    } else {
        zeta1->append();
        assert(zeta1->death() == i);

        zeta2->prepend();
        assert(zeta2->birth() == i);
    }
}

void DynamicZigzag::outwardSwitch(
    Integer i, 
    std::vector<Integer> &vanish_intvs, 
    std::vector<Integer> &new_intvs) {

    assert(i > 0 && i < this->m());
    assert(this->filt[i-1].op == ADD_OP && this->filt[i].op == DEL_OP);
    assert(this->filt[i-1].s_id != this->filt[i].s_id);

    vanish_intvs.clear();
    new_intvs.clear();

    std::shared_ptr<RepInterval> new_rep = nullptr;

    const auto sigma = this->filt[i-1].s_id;
    const auto tau = this->filt[i].s_id;

    for (auto& it : this->intervals) {
        auto& rep = it.second;

        if (rep->b == i && rep->d == i) { // Case A
            assert(rep->dimension > 0);
            assert( isSimpInChain(sigma, rep->Z(i)) );
            assert( isSimpInChain(tau, rep->Z(i)) );

            // -- rep->dimension;
            // rep->C(i-1).reset(new Chain {tau});
            // rep->C(i) = sumChain(rep->Z(i), rep->C(i-1));
            // rep->Z(i) = this->tot_cplx->getBoundary(tau);

            // new_rep = rep;

            new_rep = std::shared_ptr<RepInterval>(new RepInterval());
            new_rep->b = i;
            new_rep->d = i;
            new_rep->dimension = rep->dimension - 1;
            new_rep->cycles.resize(1);
            new_rep->chains.resize(2);

            new_rep->C(i-1) = std::shared_ptr<Chain>(new Chain {tau});
            new_rep->C(i) = sumChain(rep->Z(i), new_rep->C(i-1));
            new_rep->Z(i) = this->tot_cplx->getBoundary(tau);

            vanish_intvs.emplace_back(it.first);

        } else if (rep->d == i) { // Case B
            assert(rep->b < rep->d);
            rep->cycles.erase(rep->cycles.end()-1);
            rep->chains.erase(rep->chains.end()-2);
            -- rep->d;
            assert( isSimpInChain(tau, rep->Z(i-1)) );

        } else if (rep->b == i) { // Case C
            assert(rep->b < rep->d);
            rep->cycles.erase(rep->cycles.begin());
            rep->chains.erase(rep->chains.begin()+1);
            ++ rep->b;
            assert( isSimpInChain(sigma, rep->Z(i+1)) );

        } else if (rep->b < i && rep->d > i) { // Case D
            assert( !isSimpInChain(sigma, rep->Z(i)) );
            assert( !isSimpInChain(tau, rep->Z(i)) );

            if ( isSimpInChain(sigma, rep->C(i-1)) || isSimpInChain(tau, rep->C(i)) ) {
                auto C_iM1_sum_C_i = sumChain(rep->C(i-1), rep->C(i));
                bool sigma_in = isSimpInChain(sigma, C_iM1_sum_C_i);
                bool tau_in = isSimpInChain(tau, C_iM1_sum_C_i);

                if (!sigma_in) {
                    assert( !isSimpInChain(sigma, rep->Z(i+1)) );
                    rep->Z(i) = rep->Z(i+1);
                    rep->C(i-1) = C_iM1_sum_C_i;
                    rep->C(i) = nullptr;

                } else if (!tau_in) {
                    assert( !isSimpInChain(tau, rep->Z(i-1)) );
                    rep->Z(i) = rep->Z(i-1);
                    rep->C(i-1) = nullptr;
                    rep->C(i) = C_iM1_sum_C_i;

                } else {
                    rep->C(i).reset(new Chain {sigma});
                    rep->C(i-1) = sumChain(C_iM1_sum_C_i, rep->C(i));
                    rep->Z(i) = sumChain(rep->Z(i+1), this->tot_cplx->getBoundary(sigma));
                }
            }

        } else if (rep->b == i+1) { // Case E
            assert( !isSimpInChain(sigma, rep->Z(i+1)) );

            if ( !isSimpInChain(sigma, rep->C(i)) ) {
                rep->cycles.insert(rep->cycles.begin(), rep->Z(i+1));
                rep->chains.insert(rep->chains.begin()+1, nullptr);
                -- rep->b;
            } else {
                rep->cycles.insert( rep->cycles.begin(), 
                    sumChain(rep->Z(i+1), this->tot_cplx->getBoundary(sigma)) );

                std::shared_ptr<Chain> sigma_chn(new Chain {sigma});
                rep->C(i) = sumChain(rep->C(i), sigma_chn);
                rep->chains.insert(rep->chains.begin()+1, sigma_chn);

                -- rep->b;
            }

        } else if (rep->d == i-1) { // Case F
            assert( !isSimpInChain(tau, rep->Z(i-1)) );

            if ( !isSimpInChain(tau, rep->C(i-1)) ) {
                rep->cycles.emplace_back(rep->Z(i-1));
                rep->chains.insert(rep->chains.end()-1, nullptr);
                ++ rep->d;
            } else {
                rep->cycles.emplace_back( 
                    sumChain(rep->Z(i-1), this->tot_cplx->getBoundary(tau)) );

                std::shared_ptr<Chain> tau_chn(new Chain {tau});
                rep->C(i-1) = sumChain(rep->C(i-1), tau_chn);
                rep->chains.insert(rep->chains.end()-1, tau_chn);

                ++ rep->d;
            }

        } else { // Case G
            assert(rep->b > i+1 || rep->d < i-1);
        }
    }

    if (new_rep) { 
        assert(vanish_intvs.size() == 1);
        new_intvs.emplace_back(this->addInterval(new_rep)); 
    } else {
        assert(vanish_intvs.size() == 0);
    }

    for (auto rid : vanish_intvs) { 
        auto ret = this->intervals.erase(rid); 
        assert(ret == 1);
    }

    std::swap(this->filt[i-1], this->filt[i]);
}

void DynamicZigzag::inwardContraction(
    Integer i, 
    std::vector<Integer> &vanish_intvs, 
    std::vector<Integer> &new_intvs) {

    vanish_intvs.clear();
    new_intvs.clear();

    assert(i > 0 && i < this->m());
    assert(this->filt[i-1].op == ADD_OP && this->filt[i].op == DEL_OP);
    assert(this->filt[i-1].s_id == this->filt[i].s_id);

    bool i_is_birth = false, iM1_is_death = false;
    Integer rid_b_d_i;

    for (const auto& it : this->intervals) {
        if (it.second->b == i) { 
            assert(it.second->b == it.second->d); 
            i_is_birth = true;
            rid_b_d_i = it.first;
            break;

        } else if (it.second->d == i-1) { 
            iM1_is_death = true; 
            break;
        }
    }

    assert(i_is_birth != iM1_is_death);

    if (i_is_birth) { 
        this->inwardContractionInjective(i, vanish_intvs, new_intvs, rid_b_d_i); 
    } else { 
        this->inwardContractionSurjective(i, vanish_intvs, new_intvs); 
    }

    auto filt_iter_del_beg = this->filt.begin()+i-1;
    this->filt.erase(filt_iter_del_beg, filt_iter_del_beg+2);
}

void DynamicZigzag::inwardContractionInjective(
    Integer i, 
    std::vector<Integer> &vanish_intvs, 
    std::vector<Integer> &new_intvs, 
    Integer rid_b_d_i) {

    const auto sigma = this->filt[i-1].s_id;
    auto Z_tilde_i = this->intervals[rid_b_d_i]->Z(i);
    assert(isSimpInChain(sigma, Z_tilde_i));
    // auto p = this->intervals[rid_b_d_i]->dim();
    const auto p = this->tot_cplx->simpDim(sigma);

    for (const auto& it : this->intervals) {
        auto& rep = it.second;

        if (rep->d < i) {
            assert(rep->d < i-1);

        } else if (rep->b > i) {
            assert(rep->b > i+1);
            rep->b -= 2;
            rep->d -= 2;

        } else if (rep->b == i) {
            ;

        } else {
            auto C_bar = sumChain(rep->C(i-1), rep->C(i));
            if (rep->dim() == p-1 && isSimpInChain(sigma, C_bar)) 
                { C_bar = sumChain(C_bar, Z_tilde_i); } 
            this->inwardContractionSqueezeRep(i, rep, C_bar);
        }
    }

    vanish_intvs.emplace_back(rid_b_d_i);
    this->intervals.erase(rid_b_d_i);
}

void DynamicZigzag::inwardContractionSurjective(
    Integer i, 
    std::vector<Integer> &vanish_intvs, 
    std::vector<Integer> &new_intvs) {
    
    const auto sigma = this->filt[i-1].s_id;
    const auto p = this->tot_cplx->simpDim(sigma);

    Integer rid_d_iM1 = -1, rid_b_iP1 = -1;
    std::vector<Integer> sigma_rel_rids;

    for (auto& it : this->intervals) {
        const auto rid = it.first;
        auto& rep = it.second;

        if (rep->d < i-1) {
            ;

        } else if (rep->b > i+1) {
            rep->b -= 2;
            rep->d -= 2;

        } else if (rep->d == i-1) {
            rid_d_iM1 = rid;

        } else if (rep->b == i+1) {
            rid_b_iP1 = rid;

        } else {
            assert(rep->b < i && rep->d > i);
            auto C_iM1_sum_C_i = sumChain(rep->C(i-1), rep->C(i));

            if (rep->dim() == p-1 && isSimpInChain(sigma, C_iM1_sum_C_i)) 
                { sigma_rel_rids.emplace_back(rid); } 
            else 
                { this->inwardContractionSqueezeRep(i, rep, C_iM1_sum_C_i); }
        }
    }

    assert(rid_d_iM1 != -1);
    assert(this->intervals[rid_d_iM1]->dim() == p-1);
    assert(rid_b_iP1 != -1);
    assert(this->intervals[rid_b_iP1]->dim() == p-1);

    std::vector<bool> rid_valid(sigma_rel_rids.size(), true);

    for (auto j = 0; j < sigma_rel_rids.size(); ++j) {
        if (!rid_valid[j]) { continue; }
        auto& rep1 = this->intervals[sigma_rel_rids[j]];

        for (auto k = j+1; k < sigma_rel_rids.size(); ++k) {
            if (!rid_valid[k]) { continue; }
            auto& rep2 = this->intervals[sigma_rel_rids[k]];

            if (this->birthLessThan(rep1->b, rep2->b) &&
                this->deathLessThan(rep1->d, rep2->d)) {

                rep2 = sumRep(rep1, rep2);
                this->inwardContractionSqueezeRep(i, rep2);
                rid_valid[k] = false;

            } else if (this->birthLessThan(rep2->b, rep1->b) &&
                       this->deathLessThan(rep2->d, rep1->d)) {

                rep1 = sumRep(rep1, rep2);
                this->inwardContractionSqueezeRep(i, rep1);
                rid_valid[j] = false;
                break;
            }
        }
    }

    auto& rep_d_iM1 = this->intervals[rid_d_iM1]; // \zeta_*
    auto& rep_b_iP1 = this->intervals[rid_b_iP1]; // \zeta_\circ

    for (auto j = 0; j < sigma_rel_rids.size(); ++j) {
        if (!rid_valid[j]) { continue; }
        auto& rep = this->intervals[sigma_rel_rids[j]];

        assert(this->deathLessThan(i-1, rep->d));
        assert(this->birthLessThan(i+1, rep->b));

        if (this->birthLessThan(rep_d_iM1->b, rep->b)) {
            rep = sumRep(rep, rep_d_iM1);
            this->inwardContractionSqueezeRep(i, rep);
            rid_valid[j] = false;
        } else if (this->deathLessThan(rep_b_iP1->d, rep->d)) {
            rep = sumRep(rep, rep_b_iP1);
            this->inwardContractionSqueezeRep(i, rep);
            rid_valid[j] = false;
        }
    }

    // std::vector<Integer> sigma_rel_rids2;

    for (auto j = 0; j < sigma_rel_rids.size(); ++j) {
        if (rid_valid[j]) 
            { vanish_intvs.emplace_back(sigma_rel_rids[j]); }
    }

    if (vanish_intvs.empty()) {
        std::shared_ptr<RepInterval> new_rep(new RepInterval());

        new_rep->chains.emplace_back(copyChainPtr(rep_d_iM1->chains[0]));
        this->copyRep(rep_d_iM1, 0, rep_d_iM1->len(), new_rep);

        if (rep_b_iP1->len() == 1 && this->filt[rep_b_iP1->d].op == DEL_OP) {
            new_rep->chains.emplace_back(nullptr);

        } else {
            auto conn_chn = sumChain(rep_d_iM1->C(i-1), rep_b_iP1->C(i));
            conn_chn = sumChain(conn_chn, rep_b_iP1->C(i+1));
            new_rep->chains.emplace_back(conn_chn);

            if (rep_b_iP1->len() > 1) {
                this->copyRep(rep_b_iP1, 1, rep_b_iP1->len()-1, new_rep);
                new_rep->chains.emplace_back(copyChainPtr( rep_b_iP1->chains.back() ));
            }
        }

        new_rep->b = rep_d_iM1->b;
        new_rep->d = rep_b_iP1->d-2;
        new_rep->dimension = rep_d_iM1->dimension;

        new_intvs.emplace_back(this->addInterval(new_rep));

        assert(new_rep->len() > 0);
        assert(new_rep->len() == new_rep->cycles.size());
        assert(new_rep->len()+1 == new_rep->chains.size());

    } else {
        RepIdCmp rid_cmp = { this };
        std::sort(vanish_intvs.begin(), vanish_intvs.end(), rid_cmp);

        for (auto j = 0; j < vanish_intvs.size()-1; ++j) {
            auto rep1 = this->intervals[vanish_intvs[j]];
            auto rep2 = this->intervals[vanish_intvs[j+1]];

            assert(this->birthLessThan(rep1->b, rep2->b));
            assert(this->deathLessThan(rep2->d, rep1->d));

            auto new_rep = sumRep(rep1, rep2);
            this->inwardContractionSqueezeRep(i, new_rep);
            new_intvs.emplace_back(this->addInterval(new_rep));
        }

        {
            auto new_rep = sumRep(rep_d_iM1, this->intervals[vanish_intvs.back()]);
            this->inwardContractionSqueezeRep(i, new_rep);
            new_intvs.emplace_back(this->addInterval(new_rep));
        }

        {
            auto new_rep = sumRep(rep_b_iP1, this->intervals[vanish_intvs[0]]);
            this->inwardContractionSqueezeRep(i, new_rep);
            new_intvs.emplace_back(this->addInterval(new_rep));
        }
    }

    vanish_intvs.emplace_back(rid_d_iM1);
    vanish_intvs.emplace_back(rid_b_iP1);

    for (auto rid : vanish_intvs) { 
        auto ret = this->intervals.erase(rid); 
        assert(ret == 1);
    }

    assert(vanish_intvs.size()-new_intvs.size() == 1);
}

Integer DynamicZigzag::addInterval(std::shared_ptr<RepInterval> rep) {
    assert(this->intervals.find(this->new_intv_id) == this->intervals.end());
    this->intervals[this->new_intv_id] = rep;
    auto ret = this->new_intv_id;
    ++ this->new_intv_id;
    return ret;
}

void DynamicZigzag::inwardContractionSqueezeRep(
    Integer i, std::shared_ptr<RepInterval> rep) {

    auto C_iM1_sum_C_i = sumChain(rep->C(i-1), rep->C(i));
    this->inwardContractionSqueezeRep(i, rep, C_iM1_sum_C_i); 
}

void DynamicZigzag::inwardContractionSqueezeRep(
    Integer i, std::shared_ptr<RepInterval> rep, std::shared_ptr<Chain> C_bar) {

    auto sigma = this->filt[i-1].s_id;
    assert(!isSimpInChain(sigma, C_bar));

    rep->C(i+1) = sumChain(rep->C(i+1), C_bar);

    auto cyc_iter_del_beg = rep->cycles.begin()+i-rep->b;
    rep->cycles.erase(cyc_iter_del_beg, cyc_iter_del_beg+2);

    auto chn_iter_del_beg = rep->chains.begin()+i-rep->b;
    rep->chains.erase(chn_iter_del_beg, chn_iter_del_beg+2);

    rep->d -= 2;

    assert(rep->len() > 0);
    assert(rep->cycles.size() == rep->len());
}

// Insert the two arrows between K_{i-1} and K_i
void DynamicZigzag::outwardExpansion(
    Integer i, 
    Integer sigma,
    std::vector<Integer> &vanish_intvs, 
    std::vector<Integer> &new_intvs) {

    assert(i > 0 && i <= this->m());

    vanish_intvs.clear();
    new_intvs.clear();

    const auto p = this->tot_cplx->simpDim(sigma);

    std::vector<Integer> sigma_rel_rids;

    for (auto& it : this->intervals) {
        const auto rid = it.first;
        auto rep = it.second;

        if (rep->d < i-1) {
            ;

        } else if (rep->b > i-1) {
            rep->b += 2;
            rep->d += 2;

        } else {
            if (rep->dim() == p && isSimpInChain(sigma, rep->Z(i-1))) 
                { sigma_rel_rids.emplace_back(rid); } 
            else 
                { this->outwardExpansionStretchRep(i, sigma, rep); }
        }
    }

    if (sigma_rel_rids.empty()) {
        std::shared_ptr<RepInterval> new_rep(new RepInterval());

        new_rep->b = i;
        new_rep->d = i;
        new_rep->dimension = p-1;
        new_rep->cycles.emplace_back(this->tot_cplx->getBoundary(sigma));
        new_rep->chains.emplace_back(new Chain {sigma});
        new_rep->chains.emplace_back(new Chain {sigma});

        new_intvs.emplace_back(this->addInterval(new_rep));

    } else {
        this->outwardExpansionInjective(i, sigma, 
            sigma_rel_rids, vanish_intvs, new_intvs);
    }

    std::vector<SimpOp> ops { {DEL_OP, sigma}, {ADD_OP, sigma} };
    this->filt.insert(this->filt.begin()+i-1, ops.begin(), ops.end());
}

void DynamicZigzag::outwardExpansionInjective(
    Integer i, 
    Integer sigma,
    const std::vector<Integer> &sigma_rel_rids,
    std::vector<Integer> &vanish_intvs, 
    std::vector<Integer> &new_intvs) {

    std::vector<bool> rid_valid(sigma_rel_rids.size(), true);

    for (auto j = 0; j < sigma_rel_rids.size(); ++j) {
        if (!rid_valid[j]) { continue; }
        auto& rep1 = this->intervals[sigma_rel_rids[j]];

        for (auto k = j+1; k < sigma_rel_rids.size(); ++k) {
            if (!rid_valid[k]) { continue; }
            auto& rep2 = this->intervals[sigma_rel_rids[k]];

            if (this->birthLessThan(rep1->b, rep2->b) &&
                this->deathLessThan(rep1->d, rep2->d)) {

                rep2 = sumRep(rep1, rep2);
                this->outwardExpansionStretchRep(i, sigma, rep2);
                rid_valid[k] = false;

            } else if (this->birthLessThan(rep2->b, rep1->b) &&
                       this->deathLessThan(rep2->d, rep1->d)) {

                rep1 = sumRep(rep1, rep2);
                this->outwardExpansionStretchRep(i, sigma, rep1);
                rid_valid[j] = false;
                break;
            }
        }
    }

    for (auto j = 0; j < sigma_rel_rids.size(); ++j) {
        if (rid_valid[j]) 
            { vanish_intvs.emplace_back(sigma_rel_rids[j]); }
    }

    assert(!vanish_intvs.empty());

    RepIdCmp rid_cmp = { this };
    std::sort(vanish_intvs.begin(), vanish_intvs.end(), rid_cmp);

    for (auto j = 0; j < vanish_intvs.size()-1; ++j) {
        auto rep1 = this->intervals[vanish_intvs[j]];
        auto rep2 = this->intervals[vanish_intvs[j+1]];

        assert(this->birthLessThan(rep1->b, rep2->b));
        assert(this->deathLessThan(rep2->d, rep1->d));

        auto new_rep = sumRep(rep1, rep2);
        this->outwardExpansionStretchRep(i, sigma, new_rep);
        new_intvs.emplace_back(this->addInterval(new_rep));
    }

    {
        std::shared_ptr<RepInterval> new_rep(new RepInterval());
        auto rep = this->intervals[vanish_intvs[0]];

        new_rep->chains.emplace_back(copyChainPtr(rep->chains[0]));
        this->copyRep(rep, 0, i-rep->b, new_rep);
        new_rep->chains.emplace_back(nullptr);

        new_rep->b = rep->b;
        new_rep->d = i-1;
        new_rep->dimension = rep->dimension;

        new_intvs.emplace_back(this->addInterval(new_rep));
    }

    {
        std::shared_ptr<RepInterval> new_rep(new RepInterval());
        auto rep = this->intervals[vanish_intvs.back()];

        new_rep->chains.emplace_back(nullptr);
        this->copyRep(rep, i-1-rep->b, rep->d-i+2, new_rep);
        new_rep->chains.emplace_back( copyChainPtr(rep->chains.back()) );
        
        new_rep->b = i+1;
        new_rep->d = rep->d+2;
        new_rep->dimension = rep->dimension;

        new_intvs.emplace_back(this->addInterval(new_rep));
    }

    for (auto rid : vanish_intvs) { 
        auto ret = this->intervals.erase(rid); 
        assert(ret == 1);
    }

    assert(new_intvs.size()-vanish_intvs.size() == 1);
}

void DynamicZigzag::outwardExpansionStretchRep(
    Integer i, Integer sigma, std::shared_ptr<RepInterval> rep) {

    assert( !isSimpInChain(sigma, rep->Z(i-1)) );

    rep->cycles.insert(rep->cycles.begin()+i-rep->b, 2, rep->Z(i-1));
    rep->chains.insert(rep->chains.begin()+i-rep->b, 2, nullptr);

    rep->d += 2;
}

std::shared_ptr<RepInterval> DynamicZigzag::sumRep(
    const std::shared_ptr<RepInterval> r1, const std::shared_ptr<RepInterval> r2) {

    assert(r1->dim() == r2->dim());
    assert(r1->b != r2->b && r1->d != r2->d);
    assert(std::max(r1->b, r2->b) <= std::min(r1->d, r2->d));
    // assert(r3->cycles.size() == 0);
    // assert(r3->chains.size() == 0);

    std::shared_ptr<RepInterval> r3(new RepInterval());
    r3->dimension = r1->dim();

    std::shared_ptr<RepInterval> long_rep, short_rep;

    if (r1->b < r2->b) {
        long_rep = r1;
        short_rep = r2;
    } else {
        long_rep = r2;
        short_rep = r1;
    }

    if (this->filt[short_rep->b-1].op == ADD_OP) {
        r3->b = short_rep->b;
        r3->chains.push_back(nullptr);
    } else {
        r3->b = long_rep->b;

        if (long_rep->chains[0]) { 
            r3->chains.emplace_back(copyChainPtr(long_rep->chains[0]));
            // r3->chains.emplace_back(new Chain( *(long_rep->chains[0]) )); 
        } else { 
            r3->chains.push_back(nullptr); 
        }

        this->copyRep(long_rep, 0, short_rep->b-long_rep->b, r3);

        auto c = sumChain(short_rep->chains[0], long_rep->chains[short_rep->b-long_rep->b]);
        r3->chains.push_back(c);
    }


    auto intsec_d = std::min(r1->d, r2->d);
    this->sumRepBasic(
        short_rep, 0, intsec_d-short_rep->b+1,
        long_rep, short_rep->b-long_rep->b,
        r3);


    if (r1->d < r2->d) {
        short_rep = r1;
        long_rep = r2;
    } else {
        short_rep = r2;
        long_rep = r1;
    }

    if (this->filt[short_rep->d].op == DEL_OP) {
        r3->d = short_rep->d;
        r3->chains.push_back(nullptr);
    } else {
        r3->d = long_rep->d;

        auto c = sumChain(short_rep->chains[short_rep->d+1-short_rep->b], 
            long_rep->chains[short_rep->d+1-long_rep->b]);
        r3->chains.push_back(c);

        this->copyRep(long_rep, short_rep->d+1-long_rep->b, 
            long_rep->d-short_rep->d, r3);

        if (long_rep->chains.back()) { 
            r3->chains.emplace_back(copyChainPtr( long_rep->chains.back() )); 
            // r3->chains.emplace_back(new Chain( *(long_rep->chains.back()) )); 
        } else { 
            r3->chains.push_back(nullptr);
        }
    }

    return r3;
}

void DynamicZigzag::copyRep(const std::shared_ptr<RepInterval> from_r, 
    Integer start, Integer cnt, std::shared_ptr<RepInterval> to_r) {

    assert(start >= 0);
    assert(cnt > 0);
    assert(start+cnt <= from_r->cycles.size());

    to_r->cycles.emplace_back(copyChainPtr(from_r->cycles[start]));
    // to_r->cycles.emplace_back(new Chain( *(from_r->cycles[start]) ));

    for (auto i = start+1; i < start+cnt; ++i) {
        if (isEmptyChain(from_r->chains[i])) {
            to_r->cycles.emplace_back(to_r->cycles.back());
            to_r->chains.emplace_back(nullptr);
        } else {
            to_r->cycles.emplace_back(copyChainPtr(from_r->cycles[i]));
            // to_r->cycles.emplace_back(new Chain( *(from_r->cycles[i]) ));
            to_r->chains.emplace_back(copyChainPtr(from_r->chains[i]));
            // to_r->chains.emplace_back(new Chain( *(from_r->chains[i]) ));
        }
    }
}

void DynamicZigzag::sumRepBasic(
    const std::shared_ptr<RepInterval> r1, Integer r1_start, Integer cnt,
    const std::shared_ptr<RepInterval> r2, Integer r2_start,
    std::shared_ptr<RepInterval> r3) {

    assert(r1_start >= 0);
    assert(cnt > 0);
    assert(r1_start+cnt <= r1->cycles.size());

    assert(r2_start >= 0);
    assert(r2_start+cnt <= r2->cycles.size());

    assert(r1->b+r1_start == r2->b+r2_start);

    r3->cycles.emplace_back(
        sumChain(r1->cycles[r1_start], r2->cycles[r2_start]));

    for (auto i = 1; i < cnt; ++i) {
        if (isEmptyChain(r1->chains[r1_start+i]) && 
            isEmptyChain(r2->chains[r2_start+i])) {
            r3->cycles.push_back(r3->cycles.back());
            r3->chains.push_back(nullptr);
        } else {
            r3->cycles.emplace_back(
                sumChain(r1->cycles[r1_start+i], r2->cycles[r2_start+i]));

            r3->chains.emplace_back(
                sumChain(r1->chains[r1_start+i], r2->chains[r2_start+i]));
        }
    }
}

bool DynamicZigzag::repValid(const std::shared_ptr<RepInterval> r) {
    return true;
}

void DynamicZigzag::getPersistence(
    std::vector< std::tuple<Integer, Integer, Integer> > &persistence) {

    persistence.clear();

    for (auto& it : this->intervals) {
        auto rep = it.second;
        persistence.emplace_back(rep->b, rep->d, rep->dim());
    }
}

void DynamicZigzag::printPers(const std::string outfname) {
    std::vector< std::tuple<Integer, Integer, Integer> > pers;
    this->getPersistence(pers);
    std::sort(pers.begin(), pers.end());

    ::printPers(outfname, pers);
}
