#include "ud.h"

std::shared_ptr<Chain> UpDownPersistence::getExtChain(
    std::shared_ptr<Chain> chn) {

    std::shared_ptr<Chain> ext_chn(new Chain);
    for (auto i : *chn) { 
        assert(this->simp_ext_id[i] >= 0);
        ext_chn->emplace_back(this->simp_ext_id[i]); 
    }
    std::sort(ext_chn->begin(), ext_chn->end());

    return ext_chn;
}

void UpDownPersistence::addUpOrDownInterval(
    Integer dim, Integer birth, Integer death, 
    std::shared_ptr<Chain> cyc, 
    std::shared_ptr<Chain> chn) {

    assert(death < this->simp_cnt || birth > this->simp_cnt);

    std::shared_ptr<RepInterval> rep(new RepInterval);

    rep->dimension = dim;
    rep->b = birth;
    rep->d = death;

    rep->cycles.resize(rep->len(), this->getExtChain(cyc));
    rep->chains.resize(rep->len()+1, nullptr);
    if (death < this->simp_cnt)
        { rep->chains.back() = this->getExtChain(chn); }
    else
        { rep->chains[0] = this->getExtChain(chn); }

    this->dzz->addInterval(rep);
}

void UpDownPersistence::addClosedClosedInterval(
    Integer dim, Integer birth, Integer death, 
    std::shared_ptr<Chain> left_cyc, std::shared_ptr<Chain> chn, 
    std::shared_ptr<Chain> right_cyc) {

    assert(birth <= this->simp_cnt && death >= this->simp_cnt);

    std::shared_ptr<RepInterval> rep(new RepInterval);

    rep->dimension = dim;
    rep->b = birth;
    rep->d = death;

    if (isEmptyChain(chn)) {
        rep->cycles.resize(rep->len(), this->getExtChain(left_cyc));
        rep->chains.resize(rep->len()+1, nullptr);

    } else {
        assert(death > this->simp_cnt);

        rep->cycles.insert(rep->cycles.end(), 
            this->simp_cnt - rep->b + 1, this->getExtChain(left_cyc));
        rep->cycles.insert(rep->cycles.end(), 
            rep->d - this->simp_cnt, this->getExtChain(right_cyc));

        rep->chains.resize(rep->len()+1, nullptr);
        rep->C(this->simp_cnt) = this->getExtChain(chn);
    }

    this->dzz->addInterval(rep);
}

void UpDownPersistence::addSimplex(Integer ext_s_id) {
    Integer s_id = this->simp_cnt;
    this->simp_ext_id[s_id] = ext_s_id;
    this->simp_id_by_ext_id[ext_s_id] = s_id;
    ++ this->simp_cnt;

    std::shared_ptr<Chain> cyc(new Chain());
    std::shared_ptr<Chain> chn(new Chain({ s_id }));

    this->tot_cplx->getBoundary(ext_s_id, *cyc);
    for (auto &i : *cyc) { i = this->simp_id_by_ext_id[i]; }
    std::sort(cyc->begin(), cyc->end());
    assert(cyc->empty() || cyc->front() >= 0);
    assert(cyc->empty() || cyc->back() < s_id);

    while (!cyc->empty()) {
        auto i = this->pair[cyc->back()];
        if (i < 0) { break; }

        assert(this->cycles[i]->back() == cyc->back());

        cyc = sumChain(this->cycles[i], cyc);
        chn = sumChain(this->chains[i], chn);
    }

    if (cyc->empty()) {
        this->chains[s_id] = chn;
    } else {
        // std::cout 
        //    // << this->tot_cplx->simplices[ext_s_id].size()-2 << " "
        //    << cyc->back()+1 << " " << s_id << std::endl;

        this->pair[cyc->back()] = s_id;
        this->pair[s_id] = NEG_IND;
        this->chains[s_id] = chn;
        this->cycles[s_id] = cyc;
        // this->pivot_map[cyc->back()] = s_id;

        this->addUpOrDownInterval(
            this->tot_cplx->simpDim(ext_s_id)-1, 
            cyc->back()+1, s_id, cyc, chn);
    }
}

void UpDownPersistence::middle() {
    assert(this->simp_cnt == this->tot_cplx->simplices.size());
    this->cur_index = this->simp_cnt;

    for (Integer b = 0; b < this->simp_cnt; ++b) {
        if (this->pair[b] == UNPAIRED_IND) {
            alives.emplace_back();
            alives.back().birth = b+1;
            alives.back().left_cyc = this->chains[b];
            alives.back().chain = nullptr;
            alives.back().right_cyc = this->chains[b];
            // alives.back().right_cyc.reset( new Chain(*(this->chains[b])) );
        }
    }
}

void UpDownPersistence::deleteSimplex(Integer ext_s_id) {
    ++ this->cur_index;

    auto s_id = this->simp_id_by_ext_id[ext_s_id];

    Integer right_max_index = -1;
    Integer left_min_index = -1;

    for (Integer i = 0; i < alives.size(); ++i) {
        if (!std::binary_search(alives[i].right_cyc->begin(), 
            alives[i].right_cyc->end(), s_id)) 
        { continue; }

        if (alives[i].birth > this->simp_cnt) {
            if (right_max_index < 0 || 
                alives[i].birth > alives[right_max_index].birth) 
            { right_max_index = i; }
        } else {
            if (left_min_index < 0 || 
                alives[i].birth < alives[left_min_index].birth) 
            { left_min_index = i; }
        }
    }

    Integer sum_from_index = -1;

    if (right_max_index >= 0) {
        sum_from_index = right_max_index;

        this->addUpOrDownInterval(this->tot_cplx->simpDim(ext_s_id),
            alives[sum_from_index].birth, this->cur_index-1,
            alives[sum_from_index].right_cyc, alives[sum_from_index].chain);

    } else if (left_min_index >= 0) {
        sum_from_index = left_min_index;

        this->addClosedClosedInterval(this->tot_cplx->simpDim(ext_s_id),
            alives[sum_from_index].birth, this->cur_index-1,
            alives[sum_from_index].left_cyc, alives[sum_from_index].chain, 
            alives[sum_from_index].right_cyc);

    } else {
        alives.emplace_back();

        alives.back().birth = this->cur_index;
        alives.back().left_cyc = nullptr;
        alives.back().chain.reset(new Chain({ s_id }));

        alives.back().right_cyc.reset(new Chain());
        this->tot_cplx->getBoundary(ext_s_id, *(alives.back().right_cyc));
        for (auto &i : *(alives.back().right_cyc)) 
        { i = this->simp_id_by_ext_id[i]; }
        std::sort(alives.back().right_cyc->begin(), alives.back().right_cyc->end());

        return;
    }

    // std::cout 
    //     // << this->tot_cplx->simplices[ext_s_id].size()-1 << " "
    //     << alives[sum_from_index].birth << " " << this->cur_index-1 << std::endl;

    for (Integer i = 0; i < alives.size(); ++i) {
        if (i == sum_from_index) { continue; }

        if (!std::binary_search(alives[i].right_cyc->begin(), 
            alives[i].right_cyc->end(), s_id)) 
        { continue; }

        if (!isEmptyChain(alives[sum_from_index].left_cyc)) { 
            alives[i].left_cyc = 
                sumChain(alives[sum_from_index].left_cyc, alives[i].left_cyc); 
        }

        if (!isEmptyChain(alives[sum_from_index].chain)) { 
            alives[i].chain = 
                sumChain(alives[sum_from_index].chain, alives[i].chain);
        }

        if (!isEmptyChain(alives[sum_from_index].right_cyc)) { 
            alives[i].right_cyc = 
                sumChain(alives[sum_from_index].right_cyc, alives[i].right_cyc);
        }
    }

    alives.erase(alives.begin() + sum_from_index);
}
