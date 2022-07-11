#include "dpc.h"

#include <fstream>
#include <sstream>
#include <unordered_set>

#include "utils.h"
#include "ud.h"
#include "fzz.h"

struct SimpDimRevCmp {
    TotalComplex* tot_cplx;

    bool operator()(Integer s_id1, Integer s_id2) const {
        return tot_cplx->simpDim(s_id1) > tot_cplx->simpDim(s_id2);  
    }
};

void DynamicPC::init(
    const std::string &posfname, Integer _max_dim, 
    Integer _min_time, Integer _max_time) {

    this->max_dim = _max_dim;
    this->min_time = _min_time;
    this->max_time = _max_time;
    // this->num_pts = _num_pts;

    std::ifstream posfin(posfname);
    if (!posfin) {
        std::cerr << "ERR: cannot open file!" << std::endl;
        exit(-1);
    }

    posfin >> num_pts;
    char c;
    posfin.get(c);
    std::cout << "num_pts: " << num_pts << std::endl;

    /* Read in dynamic point cloud data */

    std::string line;
    std::vector<std::vector<Point>> dpoints;
    for (Integer t = 0; t <= max_time && !posfin.eof();) {
        std::getline(posfin, line);
        // if (line[0] == '#') { continue; }

        // std::cout << "[t: " << t << "]" << std::endl;

        if (t < min_time) {
            dpoints.emplace_back();
        } else {
            dpoints.emplace_back(num_pts);
        }

        for (Integer i = 0; i < num_pts; ++i) {
            std::getline(posfin, line);
            std::stringstream ss(line);

            Integer id;
            Decimal x, y;
            ss >> id >> x >> y;

            // std::cout << "  " << line << std::endl;
            // std::cout << "  " << id << " " << x << " " << y << std::endl;

            if (t >= min_time) {
                dpoints.back()[id].x = x;
                dpoints.back()[id].y = y;
            }
        }

        // trigger eof
        posfin.peek();
    
        ++t;
    }

    // std::cout << "pos loaded" << std::endl;

    // for (auto points : dpoints) {
    //     for (auto p : points) {
    //         std::cout << p.x << " " << p.y << std::endl;
    //     }

    //     std::cout << std::endl;
    // }
    // std::cout << dpoints.size() << std::endl;

    // if (dpoints.size() != max_time + 1) {
    //     std::cout << "ERR: dpoints.size() < max_time + 1" << std::endl;
    //     return -1;
    // }

    max_time = dpoints.size() - 1;


    /* Build all simplices up to the max dimension */

    tot_cplx.init(num_pts, max_dim);
    // for (auto i = 0; i <= max_dim; ++i) {
    //     std::cout << "dim " << i << " simp num: " << tot_cplx.dimSimpCnt(i) << std::endl;
    // }

    // for (auto simp : tot_cplx.simplices) {
    //     std::cout << simp << std::endl;
    // }


    /* Compute the distance-time curves for all edges */

    edge_dis.resize(tot_cplx.dimEnd(1));

    for (auto i = tot_cplx.dimStart(1); i < tot_cplx.dimEnd(1); ++i) {
        edge_dis[i-num_pts].resize(dpoints.size());

        auto v1 = tot_cplx.simplices[i][0];
        auto v2 = tot_cplx.simplices[i][1];
        // std::cout << v1 << " " << v2 << std::endl;

        for (auto t = min_time; t <= max_time; ++t) {
            edge_dis[i-num_pts][t] = eucDis(dpoints[t][v1], dpoints[t][v2]);
            // std::cout << "  t " << t << ": " << edge_dis[i-num_pts][t] << std::endl;
            // std::cout << dpoints[t][v1].x << " " << dpoints[t][v1].y << " " 
            //     << dpoints[t][v2].x << " " << dpoints[t][v2].y << " " 
            //     << edge_dis[i-num_pts][t] << std::endl;
        }
    }


    /* Open the vineyard file for write */

    std::string vines_fname;
    getFilePurename(posfname, &vines_fname);
    vines_fname += "_d_" + std::to_string(this->max_dim) + "_t_"
        + std::to_string(this->min_time) + "_" 
        + std::to_string(this->max_time) + "_vines.txt";

    this->vines_fout.open(vines_fname);
}

void DynamicPC::genEvent() {
    // Find all local max/min
    for (auto i = tot_cplx.dimStart(1); i < tot_cplx.dimEnd(1); ++i) {
        assert(edge_dis[i-num_pts][min_time] != edge_dis[i-num_pts][min_time+1]);

        if (edge_dis[i-num_pts][min_time] < edge_dis[i-num_pts][min_time+1]) {
            DpcEvent evt = {edge_dis[i-num_pts][min_time], (Decimal)min_time, L_LOCAL_MIN, i, -1};
            events.emplace_back(evt);
        }

        if (edge_dis[i-num_pts][max_time] < edge_dis[i-num_pts][max_time-1]) {
            DpcEvent evt = {edge_dis[i-num_pts][max_time], (Decimal)max_time, R_LOCAL_MIN, i, -1};
            events.emplace_back(evt);
        }

        if (edge_dis[i-num_pts][min_time] > edge_dis[i-num_pts][min_time+1]) {
            DpcEvent evt = {edge_dis[i-num_pts][min_time], (Decimal)min_time, L_LOCAL_MAX, i, -1};
            events.emplace_back(evt);
        }

        if (edge_dis[i-num_pts][max_time] > edge_dis[i-num_pts][max_time-1]) {
            DpcEvent evt = {edge_dis[i-num_pts][max_time], (Decimal)max_time, R_LOCAL_MAX, i, -1};
            events.emplace_back(evt);
        }

        for (auto t = min_time + 1; t < max_time; ++t) {
            assert(edge_dis[i-num_pts][t] != edge_dis[i-num_pts][t+1]);

            if (edge_dis[i-num_pts][t] < edge_dis[i-num_pts][t-1] && 
                edge_dis[i-num_pts][t] < edge_dis[i-num_pts][t+1]) {
                DpcEvent evt = {edge_dis[i-num_pts][t], (Decimal)t, LOCAL_MIN, i, -1};
                events.emplace_back(evt);
            }

            if (edge_dis[i-num_pts][t] > edge_dis[i-num_pts][t-1] && 
                edge_dis[i-num_pts][t] > edge_dis[i-num_pts][t+1]) {
                DpcEvent evt = {edge_dis[i-num_pts][t], (Decimal)t, LOCAL_MAX, i, -1};
                events.emplace_back(evt);
            }
        }
    }

    // Find all crossings
    for (auto t = min_time; t < max_time; ++t) {
        for (auto e1 = tot_cplx.dimStart(1); e1 < tot_cplx.dimEnd(1); ++e1) {
            for (auto e2 = e1+1; e2 < tot_cplx.dimEnd(1); ++e2) {
                assert(edge_dis[e1-num_pts][t] != edge_dis[e2-num_pts][t] && 
                    edge_dis[e1-num_pts][t+1] != edge_dis[e2-num_pts][t+1]);

                if ((edge_dis[e1-num_pts][t] < edge_dis[e2-num_pts][t]) != 
                    (edge_dis[e1-num_pts][t+1] < edge_dis[e2-num_pts][t+1])) {
                    
                    Decimal intsec_d, intsec_t;

                    lineIntersec((Decimal)t, edge_dis[e1-num_pts][t],
                        (Decimal)(t+1), edge_dis[e1-num_pts][t+1],
                        (Decimal)t, edge_dis[e2-num_pts][t],
                        (Decimal)(t+1), edge_dis[e2-num_pts][t+1],
                        intsec_t, intsec_d);

                    assert(intsec_t > t && intsec_t < t+1);

                    // std::cout << t << "," << edge_dis[e1-num_pts][t]
                    //     << " " << t+1 << "," << edge_dis[e1-num_pts][t+1]
                    //     << " " << t << "," << edge_dis[e2-num_pts][t]
                    //     << " " << t+1 << "," << edge_dis[e2-num_pts][t+1]
                    //     << " " << intsec_t << "," << intsec_d << std::endl;

                    DpcEvtType evt_type = OPPO_CROSS;

                    if (edge_dis[e1-num_pts][t] < edge_dis[e1-num_pts][t+1] && 
                        edge_dis[e2-num_pts][t] < edge_dis[e2-num_pts][t+1]) {
                        evt_type = INC_CROSS;
                    }

                    if (edge_dis[e1-num_pts][t] > edge_dis[e1-num_pts][t+1] && 
                        edge_dis[e2-num_pts][t] > edge_dis[e2-num_pts][t+1]) {
                        evt_type = DEC_CROSS;
                    }

                    DpcEvent evt = {intsec_d, intsec_t, evt_type, e1, e2};
                    events.emplace_back(evt);
                }
            }
        }
    }

    struct {
        bool operator()(const DpcEvent &evt1, const DpcEvent &evt2) const {
            return evt1.d > evt2.d; 
        }
    } DpcEventCmp;

    std::sort(events.begin(), events.end(), DpcEventCmp);

    for (auto evt_id = 0 ; evt_id < events.size(); ++evt_id) {
        // std::cout << events[evt_id] << std::endl;

        // this->printEvent(evt_id);

        assert(!(evt_id < events.size() - 1 && events[evt_id].d == events[evt_id+1].d));
    }

    std::cout << std::endl << "num of events: " << this->events.size() << std::endl;

    this->vines_fout << events[0].d << std::endl;
}

void DynamicPC::initEdgeFilt() {
    for (auto e = tot_cplx.dimStart(1); e < tot_cplx.dimEnd(1); ++e) {
        EdgeOp e_op = {ADD_OP, e};
        e_filt.emplace_back(e_op);
        e_time.emplace_back((Decimal)min_time);
    }

    for (auto e = tot_cplx.dimStart(1); e < tot_cplx.dimEnd(1); ++e) {
        EdgeOp e_op = {DEL_OP, e};
        e_filt.emplace_back(e_op);
        e_time.emplace_back((Decimal)max_time);
    }

    struct EdgeStartDisCmp{
        Integer min_time_;
        Integer num_pts_;
        std::vector<std::vector<Decimal>> *edge_dis_;

        bool operator()(EdgeOp op1, EdgeOp op2) const {
            return (*edge_dis_)[op1.e_id-num_pts_][min_time_] < 
                (*edge_dis_)[op2.e_id-num_pts_][min_time_]; 
        }
    };

    EdgeStartDisCmp edge_start_dis_cmp = {min_time, num_pts, &edge_dis};
    std::sort(e_filt.begin(), e_filt.begin() + tot_cplx.dimSimpCnt(1), edge_start_dis_cmp);

    for (auto i = 0; i < tot_cplx.dimSimpCnt(1); ++i) {
        assert(e_filt[i].op == ADD_OP);

        // std::cout << edge_dis[e_filt[i].e_id-num_pts][min_time] << std::endl;
    }

    struct EdgeEndDisCmp{
        Integer max_time_;
        Integer num_pts_;
        std::vector<std::vector<Decimal>> *edge_dis_;

        bool operator()(EdgeOp op1, EdgeOp op2) const {
            return (*edge_dis_)[op1.e_id-num_pts_][max_time_] > 
                (*edge_dis_)[op2.e_id-num_pts_][max_time_]; 
        }
    };

    EdgeEndDisCmp edge_end_dis_cmp = {max_time, num_pts, &edge_dis};
    std::sort(e_filt.begin() + tot_cplx.dimSimpCnt(1), e_filt.end(), edge_end_dis_cmp);

    std::cout << std::endl;
    for (auto i = tot_cplx.dimSimpCnt(1); i < e_filt.size(); ++i) {
        assert(e_filt[i].op == DEL_OP);

        // std::cout << edge_dis[e_filt[i].e_id-num_pts][max_time] << std::endl;
    }
}

void DynamicPC::initPersistence() {

    std::vector<Integer> ref_cnt(this->tot_cplx.simplices.size());

    for (auto dim = 2; dim <= this->max_dim; ++dim) {
        for (auto s_id = this->tot_cplx.dimStart(dim); 
            s_id < this->tot_cplx.dimEnd(dim);
            ++ s_id) {

            ref_cnt[s_id] = (dim+1) * dim / 2;
        }
    }

    std::ofstream filtfout;
    if (out_filt_pers) {
        filtfout.open(this->tmp_folder+pathSeparator()+"init_filt");
        if (!filtfout) 
            { ERR() << "open init_filt for write failed!" << std::endl; exit(-1); }

        // filtfout << this->max_dim << std::endl 
        //     << this->num_pts << std::endl
        //     << this->tot_cplx.simplices.size()*2 << std::endl;
    }

    UpDownPersistence ud_pers(&(this->tot_cplx), &(this->dzz));


    // Addition 

    SimpOp s_op;
    s_op.op = ADD_OP;

    for (auto i = this->tot_cplx.dimStart(0); 
        i < this->tot_cplx.dimEnd(0); ++i) { 
        ud_pers.addSimplex(i);

        s_op.s_id = i;
        this->dzz.filt.emplace_back(s_op);
         
        if (out_filt_pers) {
            filtfout << "i " << 
                containerToStr(this->tot_cplx.simplices[i], " ") 
                << std::endl;
        }
    }

    std::vector<Integer> cofaces;
    for (auto i = 0; i < this->tot_cplx.edgeCnt(); ++i) {
        auto &e_op = this->e_filt[i];
        assert(e_op.op == ADD_OP);

        this->simp_filt_pos.emplace_back(this->dzz.filt.size());

        // Add the edge
        ud_pers.addSimplex(e_op.e_id);

        s_op.s_id = e_op.e_id;
        this->dzz.filt.emplace_back(s_op);

        if (out_filt_pers) {
            filtfout << "i " 
                << containerToStr(this->tot_cplx.simplices[e_op.e_id], " ") 
                << std::endl;
        }


        for (auto cof_dim = 2; cof_dim <= this->max_dim; ++cof_dim) {

            // std::cout << this->tot_cplx.simplices[e_op.e_id] 
            //     << " dim " << cof_dim << ":" << std::endl;

            cofaces.clear();
            this->tot_cplx.getCofaces(e_op.e_id, cof_dim, cofaces);
            assert(cofaces.size() == combChoose(this->num_pts-2, cof_dim-1));

            for (auto cof_id : cofaces) {
                // std::cout << "  " << this->tot_cplx.simplices[cof_id] << std::endl;

                --ref_cnt[cof_id];
                if (ref_cnt[cof_id] == 0) { 
                    ud_pers.addSimplex(cof_id); 

                    s_op.s_id = cof_id;
                    this->dzz.filt.emplace_back(s_op);

                    if (out_filt_pers) {
                        filtfout << "i " 
                            << containerToStr(this->tot_cplx.simplices[cof_id], " ") 
                            << std::endl;
                    }
                }
            }
        }
    }

    for (auto c : ref_cnt) { assert(c == 0); }

    ref_cnt.clear();
    ref_cnt.shrink_to_fit();


    ud_pers.middle();


    // Deletion 

    s_op.op = DEL_OP;

    std::vector<bool> deleted(this->tot_cplx.simplices.size(), false);

    for (auto i = this->tot_cplx.edgeCnt(); i < this->e_filt.size(); ++i) {
        auto &e_op = this->e_filt[i];
        assert(e_op.op == DEL_OP);

        this->simp_filt_pos.emplace_back(this->dzz.filt.size());

        for (auto cof_dim = this->max_dim; cof_dim >= 2;  --cof_dim) {

            // std::cout << this->tot_cplx.simplices[e_op.e_id] 
            //     << " dim " << cof_dim << ":" << std::endl;

            cofaces.clear();
            this->tot_cplx.getCofaces(e_op.e_id, cof_dim, cofaces);
            // assert(cofaces.size() == combChoose(this->num_pts-2, cof_dim-1));

            for (auto cof_id : cofaces) {
                // std::cout << "  " << this->tot_cplx.simplices[cof_id] << std::endl;

                if (!deleted[cof_id]) {
                    ud_pers.deleteSimplex(cof_id);

                    s_op.s_id = cof_id;
                    this->dzz.filt.emplace_back(s_op);

                    if (out_filt_pers) {
                        filtfout << "d " 
                            << containerToStr(this->tot_cplx.simplices[cof_id], " ") 
                            << std::endl;
                    }

                    deleted[cof_id] = true;
                }
            }
        }

        ud_pers.deleteSimplex(e_op.e_id);

        s_op.s_id = e_op.e_id;
        this->dzz.filt.emplace_back(s_op);

        if (out_filt_pers) {
            filtfout << "d " 
                << containerToStr(this->tot_cplx.simplices[e_op.e_id], " ") 
                << std::endl;
        }
    }

    this->simp_filt_pos.emplace_back(this->dzz.filt.size());

    for (auto i = this->tot_cplx.dimStart(0); 
        i < this->tot_cplx.dimEnd(0); ++i) { 
        ud_pers.deleteSimplex(i); 

        s_op.s_id = i;
        this->dzz.filt.emplace_back(s_op);

        if (out_filt_pers) {
            filtfout << "d " << 
                containerToStr(this->tot_cplx.simplices[i], " ") 
                << std::endl;
        }
    }

    assert(this->dzz.filt.size() == this->tot_cplx.simplices.size()*2);
    assert(this->e_filt.size()+1 == this->simp_filt_pos.size());

    if (out_filt_pers) 
        { this->dzz.printPers(tmp_folder+pathSeparator()+"init_pers"); }

    if (!run_update) { return; }

    for (auto& it : this->dzz.intervals) {
        const auto rid = it.first;
        const auto& rep = it.second;
        
        Integer e_birth, e_death;
        Decimal birth_time, death_time;
        bool is_form_vine_intv = this->getEdgeFiltInterval(
            rid, e_birth, e_death, birth_time, death_time);

        if (is_form_vine_intv) {
            auto& vine_intv = this->vine_intervals[rid];
            vine_intv.dimension = rep->dimension;
            vine_intv.birth_e = e_birth;
            vine_intv.death_e = e_death;
            vine_intv.birth_time = birth_time;
            vine_intv.death_time = death_time;

            assert(vine_intv.birth_time == this->min_time);
            assert(vine_intv.death_time == this->max_time);

            std::ostringstream oss;
            oss << "s " << rid << " " << vine_intv.birth_time << " " 
                << vine_intv.death_time << " inf " << vine_intv.dimension;
            if (print_vines) { std::cout << oss.str() << std::endl; }
            this->vines_fout << oss.str() << std::endl;
        }
    }
}

void DynamicPC::travEvent() {
    Decimal t1_old, t2_old;

    assert(this->events.back().type == LOCAL_MIN || 
        this->events.back().type == L_LOCAL_MIN || 
        this->events.back().type == R_LOCAL_MIN);

    for (auto evt_id = 0 ; evt_id < this->events.size(); ++ evt_id) {
        const auto& evt = this->events[evt_id];

        // std::cout << "evt dis " << evt.d << " " << this->events[evt_id+1].d << std::endl;

        // this->printEvent(evt_id);
        // this->printEFilt();
        // std::cout << "e_filt len: " << this->e_filt.size() << std::endl;

        if (this->e_filt.size() > max_e_filt_len) {
            max_e_filt_len = this->e_filt.size();
        }

        Decimal dis = -1.0;
        if (evt_id < this->events.size() - 1) {
            dis = (evt.d + this->events[evt_id+1].d) / 2;
        }

        if (evt.type == LOCAL_MIN || evt.type == L_LOCAL_MIN || evt.type == R_LOCAL_MIN) {
            this->localMin(evt_id);

        } else if (evt.type == LOCAL_MAX) {
            this->localMax(evt_id);

        } else if (evt.type == INC_CROSS) {
            this->incCrossing(evt_id);
            
        } else if (evt.type == DEC_CROSS) {
            this->decCrossing(evt_id);
            
        } else if (evt.type == OPPO_CROSS) {
            this->oppoCrossing(evt_id);
            
        } else if (evt.type == L_LOCAL_MAX) {
            
            // Only for checking
            auto iter = std::upper_bound(this->e_time.begin(), this->e_time.end(), evt.t);
            Integer i = (iter - this->e_time.begin()) - 1;
            assert(i >= 0 && i < this->e_filt.size());
            assert(this->e_filt[i].e_id == evt.e1);
            assert(this->e_filt[i].op == ADD_OP);
            assert(this->e_time[i] == (Decimal)this->min_time);

            auto switch_i = this->evalFiltETime(dis, t1_old, t2_old);

            assert(switch_i < 0);
            assert(this->e_time[i] > (Decimal)this->min_time);

            std::vector<Integer> vanish_intvs;
            std::unordered_set<Integer> new_intvs;
            this->updateVines(evt_id, vanish_intvs, new_intvs);

        } else if (evt.type == R_LOCAL_MAX) {

            // Only for checking
            auto iter = std::lower_bound(this->e_time.begin(), this->e_time.end(), evt.t);
            Integer i = iter - this->e_time.begin();
            assert(i >= 0 && i < this->e_filt.size());
            assert(this->e_filt[i].e_id == evt.e1);
            assert(this->e_filt[i].op == DEL_OP);
            assert(this->e_time[i] == (Decimal)this->max_time);

            auto switch_i = this->evalFiltETime(dis, t1_old, t2_old);

            assert(switch_i < 0);
            assert(this->e_time[i] < (Decimal)this->max_time);

            std::vector<Integer> vanish_intvs;
            std::unordered_set<Integer> new_intvs;
            this->updateVines(evt_id, vanish_intvs, new_intvs);
            
        } else {
            assert(0);
        }

        assert(this->simp_filt_pos.size() == this->e_filt.size() + 1);
        assert(this->simp_filt_pos[0] == this->num_pts);
        assert( std::is_sorted(this->simp_filt_pos.begin(), 
            this->simp_filt_pos.end()) );

        // std::cout << std::endl << "max_e_filt_len: " << max_e_filt_len << std::endl;
        // std::cout << "max_simp_filt_len: " << max_simp_filt_len << std::endl;
        // if (run_update) { std::cout << "up time: " << this->up_timing.count() << std::endl; }
        // if (run_fzz) { std::cout << "fzz time: " << this->fzz_timing.count() << std::endl; }
    }

    assert(this->e_filt.size() == 0 && this->e_time.size() == 0);
    assert(this->vine_intervals.size() == 0);
}

void DynamicPC::localMin(Integer evt_id) {

    Decimal t1_old, t2_old;
    const auto& evt = this->events[evt_id];

    Decimal dis = -1.0;
    if (evt_id < this->events.size() - 1) {
        dis = (evt.d + this->events[evt_id+1].d) / 2;
    }

    std::vector<Decimal>::iterator iter;

    if (evt.type == L_LOCAL_MIN) {
        iter = std::upper_bound(this->e_time.begin(), this->e_time.end(), evt.t);
    } else {
        iter = std::lower_bound(this->e_time.begin(), this->e_time.end(), evt.t);
    }

    auto i = iter - this->e_time.begin();

    assert(this->e_filt[i-1].e_id == evt.e1 && this->e_filt[i].e_id == evt.e1 && 
        this->e_filt[i-1].op == ADD_OP && this->e_filt[i].op == DEL_OP);

    std::ostringstream oss;
    oss << "[E]: " << evt.d << " erase " << i-1;

    if (print_e_filt_op) { std::cout << std::endl << oss.str() << std::endl; }
    if (out_filt_pers) { oper_fout << std::endl << oss.str() << std::endl; }

    std::vector<Integer> vanish_intvs;
    std::unordered_set<Integer> new_intvs;

    if (exec_simp_filt_op) { this->localMinExec(i-1, vanish_intvs, new_intvs); }

    this->e_filt.erase(this->e_filt.begin()+i-1, this->e_filt.begin()+i+1);
    this->e_time.erase(this->e_time.begin()+i-1, this->e_time.begin()+i+1);

    if (evt_id < this->events.size() - 1) {
        auto switch_i = this->evalFiltETime(dis, t1_old, t2_old);
        assert(switch_i < 0);
    }

    assert(isFiltValid());

    if (evt_id < this->events.size() - 1) {
        this->updateVines(evt_id, vanish_intvs, new_intvs);

    } else if (run_update) {
        for (auto rid : vanish_intvs) {
            if (this->vine_intervals.count(rid) > 0) 
            { this->endVineInterval(rid, this->events[evt_id].d); }
        }

        assert(this->dzz.intervals.size() == this->num_pts);

        for (auto& it : this->dzz.intervals) {
            const auto rid = it.first;
            const auto& rep = it.second;
            assert(rep->dimension == 0);

            if (this->vine_intervals.count(rid) > 0) {
                assert(this->vine_intervals[rid].dimension == 0);

                std::ostringstream oss;
                oss << "e " << rid << " " << this->min_time
                    << " " << this->max_time << " " << 0;
                if (print_vines) { std::cout << oss.str() << std::endl; }
                this->vines_fout << oss.str() << std::endl;

                this->vine_intervals.erase(rid);
            }
        }
    }
}

void DynamicPC::localMax(Integer evt_id) {

    Decimal t1_old, t2_old;
    const auto& evt = this->events[evt_id];

    Decimal dis = -1.0;
    if (evt_id < this->events.size() - 1) {
        dis = (evt.d + this->events[evt_id+1].d) / 2;
    }

    auto lt = this->descendFindTime(evt.e1, 'l', (Integer)evt.t - 1, dis);
    auto rt = this->descendFindTime(evt.e1, 'r', (Integer)evt.t + 1, dis);
    
    auto switch_i = this->evalFiltETime(dis, t1_old, t2_old);
    assert(switch_i < 0);

    auto iter = std::lower_bound(this->e_time.begin(), this->e_time.end(), evt.t);
    auto i = iter - this->e_time.begin();
    assert(*(iter-1) < lt && *iter > rt);

    this->e_time.insert(iter, 2, -1.0);
    this->e_time[i] = lt;
    this->e_time[i+1] = rt;

    this->e_filt.insert(this->e_filt.begin()+i, 2, EdgeOp());
    this->e_filt[i].op = DEL_OP;
    this->e_filt[i].e_id = evt.e1;
    this->e_filt[i+1].op = ADD_OP;
    this->e_filt[i+1].e_id = evt.e1;

    std::ostringstream oss;
    oss << "[E]: " << evt.d << " insert " << i-1 << " " << evt.e1;

    if (print_e_filt_op) { std::cout << std::endl << oss.str() << std::endl; }
    if (out_filt_pers) { oper_fout << std::endl << oss.str() << std::endl; }

    std::vector<Integer> vanish_intvs;
    std::unordered_set<Integer> new_intvs;

    if (exec_simp_filt_op) { this->localMaxExec(i, evt.e1, vanish_intvs, new_intvs); }

    assert(isFiltValid());

    this->updateVines(evt_id, vanish_intvs, new_intvs);
}

void DynamicPC::incCrossing(Integer evt_id) {

    Decimal t1_old, t2_old;
    const auto& evt = this->events[evt_id];

    Decimal dis = -1.0;
    if (evt_id < this->events.size() - 1) {
        dis = (evt.d + this->events[evt_id+1].d) / 2;
    }

    auto switch_i = this->evalFiltETime(dis, t1_old, t2_old);

    assert(switch_i >= 0);
    assert(this->e_time[switch_i+1] < evt.t && evt.t < t1_old);
    assert(this->e_filt[switch_i].op == DEL_OP && 
        this->e_filt[switch_i+1].op == DEL_OP);
    assert((evt.e1 == this->e_filt[switch_i].e_id && 
        evt.e2 == this->e_filt[switch_i+1].e_id) ||
        (evt.e2 == this->e_filt[switch_i].e_id && 
        evt.e1 == this->e_filt[switch_i+1].e_id));

    std::ostringstream oss;
    oss << std::endl << "[E]: " 
            << evt.d << " switch " << switch_i 
            // << " " << evt.e1 << " " << evt.e2 
            << " inc";

    if (print_e_filt_op) { std::cout << std::endl << oss.str() << std::endl; }
    if (out_filt_pers) { oper_fout << std::endl << oss.str() << std::endl; }

    if (exec_simp_filt_op) { this->incCrossingExec(switch_i); }
    
    std::swap(this->e_filt[switch_i], this->e_filt[switch_i+1]);

    assert(isFiltValid());

    std::vector<Integer> vanish_intvs;
    std::unordered_set<Integer> new_intvs;
    this->updateVines(evt_id, vanish_intvs, new_intvs);
}

void DynamicPC::decCrossing(Integer evt_id) {

    Decimal t1_old, t2_old;
    const auto& evt = this->events[evt_id];

    Decimal dis = -1.0;
    if (evt_id < this->events.size() - 1) {
        dis = (evt.d + this->events[evt_id+1].d) / 2;
    }

    auto switch_i = this->evalFiltETime(dis, t1_old, t2_old);

    assert(switch_i >= 0);
    assert(t2_old < evt.t && evt.t < this->e_time[switch_i]);
    assert(this->e_filt[switch_i].op == ADD_OP && 
        this->e_filt[switch_i+1].op == ADD_OP);
    assert((evt.e1 == this->e_filt[switch_i].e_id && 
        evt.e2 == this->e_filt[switch_i+1].e_id) ||
        (evt.e2 == this->e_filt[switch_i].e_id && 
        evt.e1 == this->e_filt[switch_i+1].e_id));

    std::ostringstream oss;
    oss << "[E]: " << evt.d << " switch " << switch_i << " dec";

    if (print_e_filt_op) { std::cout << std::endl << oss.str() << std::endl; }
    if (out_filt_pers) { oper_fout << std::endl << oss.str() << std::endl; }

    if (exec_simp_filt_op) { this->decCrossingExec(switch_i); }
    
    std::swap(this->e_filt[switch_i], this->e_filt[switch_i+1]);

    assert(isFiltValid());

    std::vector<Integer> vanish_intvs;
    std::unordered_set<Integer> new_intvs;
    this->updateVines(evt_id, vanish_intvs, new_intvs);
}

void DynamicPC::oppoCrossing(Integer evt_id) {

    Decimal t1_old, t2_old;
    const auto& evt = this->events[evt_id];

    Decimal dis = -1.0;
    if (evt_id < this->events.size() - 1) {
        dis = (evt.d + this->events[evt_id+1].d) / 2;
    }

    auto switch_i = this->evalFiltETime(dis, t1_old, t2_old);

    assert(switch_i >= 0);
    assert(t1_old < evt.t && evt.t < t2_old &&
        this->e_time[switch_i] < evt.t &&
        evt.t < this->e_time[switch_i+1]);
    assert(this->e_filt[switch_i].op == ADD_OP && 
        this->e_filt[switch_i+1].op == DEL_OP);
    assert((evt.e1 == this->e_filt[switch_i].e_id && 
        evt.e2 == this->e_filt[switch_i+1].e_id) ||
        (evt.e2 == this->e_filt[switch_i].e_id && 
        evt.e1 == this->e_filt[switch_i+1].e_id));

    std::ostringstream oss;
    oss << "[E]: " << evt.d << " switch " << switch_i << " oppo";

    if (print_e_filt_op) { std::cout << std::endl << oss.str() << std::endl; }
    if (out_filt_pers) { oper_fout << std::endl << oss.str() << std::endl; }

    std::vector<Integer> vanish_intvs;
    std::unordered_set<Integer> new_intvs;

    if (exec_simp_filt_op) 
        { this->oppoCrossingExec(switch_i, vanish_intvs, new_intvs); }
    
    std::swap(this->e_filt[switch_i], this->e_filt[switch_i+1]);

    assert(isFiltValid());

    this->updateVines(evt_id, vanish_intvs, new_intvs);
}

void DynamicPC::endVineInterval(const Integer rid, const Decimal evt_dis) {
    const auto& edge_intv = this->vine_intervals.at(rid);

    auto b = this->descendFindTimeNoDir(
        edge_intv.birth_e, edge_intv.birth_time, evt_dis);
    auto d = this->descendFindTimeNoDir(
        edge_intv.death_e, edge_intv.death_time, evt_dis);

    std::ostringstream oss;
    oss << "e " << rid << " " << b << " " << d << " " << evt_dis;
    if (print_vines) { std::cout << oss.str() << std::endl; }
    this->vines_fout << oss.str() << std::endl;

    this->vine_intervals.erase(rid);
}

void DynamicPC::updateVines(
    const Decimal evt_id,
    const std::vector<Integer> &vanish_intvs, 
    const std::unordered_set<Integer> &new_intvs) {

    if (!run_update) { return; }

    assert(evt_id < this->events.size() - 1);

    const Decimal evt_dis = this->events[evt_id].d;

    Decimal intm_dis = -1.0;
    if (evt_id < this->events.size() - 1) {
        intm_dis = (evt_dis + this->events[evt_id+1].d) / 2;
    }

    for (auto rid : vanish_intvs) {
        if (this->vine_intervals.count(rid) > 0) {
            this->endVineInterval(rid, evt_dis);
        }
    }

    for (auto& it : this->dzz.intervals) {
        const auto rid = it.first;
        const auto& rep = it.second;
        
        Integer e_birth, e_death;
        Decimal birth_time, death_time;
        bool is_form_vine_intv = this->getEdgeFiltInterval(
            rid, e_birth, e_death, birth_time, death_time);

        if (is_form_vine_intv) {

            bool is_new_vine_intv = (this->vine_intervals.count(rid) == 0);

            // std::cout << "is_form_vine_intv " << e_birth << " " << e_death 
            //     << " " << birth_time << " " << death_time << std::endl;
            // std::cout << (e_birth >= 0 ? this->edgeDisInterp(e_birth, birth_time) : -1.0)
            //     << " " << (e_death >= 0 ?  this->edgeDisInterp(e_death, death_time) : -1.0)
            //     << std::endl;

            if (is_new_vine_intv) {
                // insert the new vine interval
                this->vine_intervals[rid].dimension = rep->dimension;
                auto b = this->ascendFindTimeNoDir(e_birth, birth_time, evt_dis);
                auto d = this->ascendFindTimeNoDir(e_death, death_time, evt_dis);

                {
                    std::ostringstream oss;
                    oss << "s " << rid << " " << b << " " << d 
                        << " " << evt_dis << " " << rep->dimension;
                    if (print_vines) { std::cout << oss.str() << std::endl; }
                    this->vines_fout << oss.str() << std::endl;
                }
            }

            auto& vine_intv = this->vine_intervals[rid];
            assert(vine_intv.dimension == rep->dimension);
            vine_intv.birth_e = e_birth;
            vine_intv.death_e = e_death;
            vine_intv.birth_time = birth_time;
            vine_intv.death_time = death_time;

            {
                std::ostringstream oss;
                oss << "c " << rid << " " << vine_intv.birth_time
                    << " " << vine_intv.birth_time << " " << intm_dis;
                if (print_vines) { std::cout << oss.str() << std::endl; }
                this->vines_fout << oss.str()  << std::endl;
            }

        } else {
            if (this->vine_intervals.count(rid) > 0) // vine interval ends
                { this->endVineInterval(rid, evt_dis); }
        }
    }
}

bool DynamicPC::getEdgeFiltInterval(
    const Integer rid, 
    Integer &birth_e_id, Integer &death_e_id,
    Decimal &birth_time, Decimal &death_time) {

    Integer birth_pos = -1, death_pos = -1;

    const auto b = this->dzz.intervals[rid]->b;
    const auto d = this->dzz.intervals[rid]->d;

    auto b_iter = std::upper_bound(
        this->simp_filt_pos.begin(), this->simp_filt_pos.end(), b-1);

    assert(b_iter != this->simp_filt_pos.end());

    if (b_iter == this->simp_filt_pos.begin()) {
        birth_e_id = this->LEFT_MIN_E;
        birth_time = this->min_time;
    } else {
        birth_pos = b_iter - this->simp_filt_pos.begin() - 1;
        assert(birth_pos >= 0 && birth_pos < this->e_filt.size());

        birth_e_id = this->e_filt[birth_pos].e_id;
        birth_time = this->e_time[birth_pos];
    }

    auto d_iter = std::upper_bound(
        this->simp_filt_pos.begin(), this->simp_filt_pos.end(), d);

    assert(d_iter != this->simp_filt_pos.begin());
    
    if (d_iter == this->simp_filt_pos.end()) {
        death_e_id = this->RIGHT_MAX_E;
        death_time = this->max_time;
    } else {
        death_pos = d_iter - this->simp_filt_pos.begin() - 1;
        assert(death_pos >= 0 && death_pos < this->e_filt.size());

        death_e_id = this->e_filt[death_pos].e_id;
        death_time = this->e_time[death_pos];
    }

    if (birth_time == this->min_time || death_time == this->max_time) 
        { return birth_time < death_time; }

    assert(birth_pos <= death_pos);
    return birth_pos < death_pos;
}

void DynamicPC::addVanishNewIntervals(
    std::vector<Integer> &vanish_intvs, 
    std::unordered_set<Integer> &new_intvs,
    const std::vector<Integer> &vanish_intvs_local, 
    const std::vector<Integer> &new_intvs_local) {

    for (auto rid : vanish_intvs_local) {
        if (new_intvs.count(rid) > 0)
            { new_intvs.erase(rid); }
        else
            { vanish_intvs.emplace_back(rid); }
    }

    for (auto rid : new_intvs_local) {
        assert(new_intvs.count(rid) == 0);
        new_intvs.insert(rid);
    }
}

void DynamicPC::localMinExec(
    Integer i, 
    std::vector<Integer> &vanish_intvs, 
    std::unordered_set<Integer> &new_intvs) {

    assert(i+2 < this->simp_filt_pos.size());
    assert(this->simp_filt_pos[i+1] - this->simp_filt_pos[i] == 
        this->simp_filt_pos[i+2] - this->simp_filt_pos[i+1]);

    auto s_cnt = this->simp_filt_pos[i+1] - this->simp_filt_pos[i];

    std::unordered_map<Integer, Integer> dest_pos_map;
    for (auto j = this->simp_filt_pos[i+1]; j < this->simp_filt_pos[i+2]; ++j) {
        const auto& s_op = this->dzz.filt[j];
        assert(s_op.op == DEL_OP);

        dest_pos_map[s_op.s_id] = s_cnt-1 - (j-this->simp_filt_pos[i+1]);
    }

    std::vector<Integer> cur_pos_by_dest_pos(s_cnt, -1);
    for (auto j = this->simp_filt_pos[i]; j < this->simp_filt_pos[i+1]; ++j) {
        const auto& s_op = this->dzz.filt[j];
        assert(s_op.op == ADD_OP);

        assert(dest_pos_map.find(s_op.s_id) != dest_pos_map.end());
        cur_pos_by_dest_pos.at(dest_pos_map[s_op.s_id]) = j;
    }

    const auto start_pos = this->simp_filt_pos[i];
    for (auto dest_pos = s_cnt-1; dest_pos >= 0; --dest_pos) {
        assert(cur_pos_by_dest_pos[dest_pos] <= dest_pos+start_pos);

        for (auto k = cur_pos_by_dest_pos[dest_pos]; k < dest_pos+start_pos; ++k) {
            auto cur_sid = this->dzz.filt[k].s_id;
            assert(dest_pos_map[cur_sid] == dest_pos);

            cur_pos_by_dest_pos[dest_pos] ++;

            auto next_sid = this->dzz.filt[k+1].s_id;
            cur_pos_by_dest_pos[dest_pos_map[next_sid]] --;

            this->forwardSwitch(k+1);
        }
    }

    std::vector<Integer> vanish_intvs_local, new_intvs_local;
    for (auto j = s_cnt-1; j >= 0; --j) { 
        this->inwardContraction(start_pos+j+1, 
            vanish_intvs_local, new_intvs_local); 

        this->addVanishNewIntervals(vanish_intvs, new_intvs, 
            vanish_intvs_local, new_intvs_local);
    }


    for (auto j = i+2; j < this->simp_filt_pos.size(); ++j) 
        { this->simp_filt_pos[j] -= 2 * s_cnt; }

    this->simp_filt_pos.erase(
        this->simp_filt_pos.begin()+i, 
        this->simp_filt_pos.begin()+i+2);

    assert(i == 0 || this->simp_filt_pos[i-1] < this->simp_filt_pos[i]);
}

// inserts the two edge operations before the i-th edge operation
void DynamicPC::localMaxExec(
    Integer i, Integer e_id, 
    std::vector<Integer> &vanish_intvs, 
    std::unordered_set<Integer> &new_intvs) {

    const auto c = this->simp_filt_pos[i];
    std::vector<Integer> rel_simps;

    {
        std::unordered_set<Integer> cur_simps;
        for (auto j = 0; j < c; ++j) {
            const auto& s_op = this->dzz.filt[j];

            if (s_op.op == ADD_OP) {
                assert(cur_simps.count(s_op.s_id) == 0);
                cur_simps.insert(s_op.s_id);
            } else {
                assert(cur_simps.count(s_op.s_id) == 1);
                cur_simps.erase(s_op.s_id);
            }
        }

        for (const auto s_id : cur_simps) {
            if (this->tot_cplx.isSimpFaceOf(e_id, s_id)) 
                { rel_simps.emplace_back(s_id); }
        }
    }

    assert(rel_simps.size() > 0);

    SimpDimRevCmp cmp { &(this->tot_cplx) };
    std::sort(rel_simps.begin(), rel_simps.end(), cmp);

    std::vector<Integer> vanish_intvs_local, new_intvs_local;
    for (auto j = 0; j < rel_simps.size(); ++j) {
        assert( j == 0 || this->tot_cplx.simpDim(rel_simps.at(j-1)) 
            >= this->tot_cplx.simpDim(rel_simps.at(j)) );

        this->outwardExpansion(c+j+1, rel_simps[j], 
            vanish_intvs_local, new_intvs_local);

        this->addVanishNewIntervals(vanish_intvs, new_intvs, 
            vanish_intvs_local, new_intvs_local);
    }


    for (auto j = i; j < this->simp_filt_pos.size(); ++j) 
        { this->simp_filt_pos[j] += 2 * rel_simps.size(); }

    this->simp_filt_pos.insert(this->simp_filt_pos.begin() + i, 2, c);
    this->simp_filt_pos[i+1] = c + rel_simps.size();
}

void DynamicPC::incCrossingExec(Integer i) {
    assert(i+2 < this->simp_filt_pos.size());

    const auto e1 = this->e_filt[i].e_id;
    const auto e2 = this->e_filt[i+1].e_id;

    auto next_start = this->simp_filt_pos[i+2];

    for (auto j = this->simp_filt_pos[i+1]-1; 
        j >= this->simp_filt_pos[i]; --j) {

        const auto s_op = this->dzz.filt[j];
        assert(this->tot_cplx.isSimpFaceOf(e1, s_op.s_id));

        if (!this->tot_cplx.isSimpFaceOf(e2, s_op.s_id)) {
            assert(j < next_start - 1);

            for (auto k = j; k < next_start - 1; ++k) 
                { this->backwardSwitch(k+1); }

            -- next_start;
        }
    }

    assert(this->simp_filt_pos[i] < next_start && 
        next_start < this->simp_filt_pos[i+2]);
    this->simp_filt_pos[i+1] = next_start;
}

void DynamicPC::decCrossingExec(Integer i) {
    assert(i+2 < this->simp_filt_pos.size());

    const auto e1 = this->e_filt[i].e_id;
    const auto e2 = this->e_filt[i+1].e_id;

    // start of e1 after switch
    auto next_start = this->simp_filt_pos[i]; 

    for (auto j = this->simp_filt_pos[i+1]; 
        j < this->simp_filt_pos[i+2]; ++j) {

        const auto s_op = this->dzz.filt[j];
        assert(this->tot_cplx.isSimpFaceOf(e2, s_op.s_id));

        if (!this->tot_cplx.isSimpFaceOf(e1, s_op.s_id)) {
            assert(j > next_start);

            for (auto k = j; k > next_start; --k) 
                { this->forwardSwitch(k); }

            assert(s_op.s_id == this->dzz.filt[next_start].s_id);
            ++ next_start;
        }
    }

    assert(this->simp_filt_pos[i] < next_start && 
        next_start < this->simp_filt_pos[i+2]);
    this->simp_filt_pos[i+1] = next_start;
}

void DynamicPC::oppoCrossingExec(
    Integer i, 
    std::vector<Integer> &vanish_intvs, 
    std::unordered_set<Integer> &new_intvs) {

    assert(i+2 < this->simp_filt_pos.size());

    const auto e1 = this->e_filt[i].e_id;
    const auto e2 = this->e_filt[i+1].e_id;

    auto j = this->simp_filt_pos[i+1];
    auto end = this->simp_filt_pos[i+2]-1;
    auto next_start = this->simp_filt_pos[i];
    Integer num_inw_contra = 0;
    std::vector<Integer> vanish_intvs_local, new_intvs_local;

    while (j <= end) {
        const auto s_op = this->dzz.filt[j];
        assert(this->tot_cplx.isSimpFaceOf(e2, s_op.s_id));

        auto k = j;
        for (; k > next_start; --k) { 
            assert(this->dzz.filt[k].s_id == s_op.s_id);
            assert(k-1 >= this->simp_filt_pos[i]);
            assert(this->tot_cplx.isSimpFaceOf(e1, this->dzz.filt[k-1].s_id));

            if (this->dzz.filt[k].s_id == this->dzz.filt[k-1].s_id) {
                this->inwardContraction(k, vanish_intvs_local, new_intvs_local);

                this->addVanishNewIntervals(vanish_intvs, new_intvs, 
                    vanish_intvs_local, new_intvs_local);

                ++ num_inw_contra;
                break;
            }

            this->outwardSwitch(k, vanish_intvs_local, new_intvs_local);
            this->addVanishNewIntervals(vanish_intvs, new_intvs, 
                vanish_intvs_local, new_intvs_local);
        }

        if (k == next_start) {
            ++ next_start;
            ++ j;
        } else {
            j -= 1;
            end -= 2;
        }
    }

    if (num_inw_contra > 0) {
        for (j = i+2; j < this->simp_filt_pos.size(); ++j) 
            { this->simp_filt_pos[j] -= 2 * num_inw_contra; }
    }

    assert(this->simp_filt_pos[i] < next_start && 
        next_start < this->simp_filt_pos[i+2]);
    this->simp_filt_pos[i+1] = next_start;
}

void DynamicPC::forwardSwitch(Integer i) {

    std::ostringstream oss;
    oss << this->oper_cnt << ": fw_sw " << i;
    if (print_simp_filt_op) { std::cout << oss.str() << this->oper_end_str; }
    if (out_filt_pers) { oper_fout << oss.str() << std::endl; }

    if (!run_update) {
        assert(i > 0 && i < this->dzz.m());
        assert(this->dzz.filt[i-1].op == ADD_OP && this->dzz.filt[i].op == ADD_OP);
        assert( !this->tot_cplx.isSimpFaceOf(
            this->dzz.filt[i-1].s_id, this->dzz.filt[i].s_id) );

        std::swap(this->dzz.filt[i-1], this->dzz.filt[i]);

    } else {
        auto start = std::chrono::high_resolution_clock::now();
        this->dzz.forwardSwitch(i);
        auto end = std::chrono::high_resolution_clock::now();
        this->up_timing += end - start;
    }

    this->runFZZ();

    ++ this->fw_sw_cnt;
    ++ this->oper_cnt;

    this->afterSimplexOper();
}

void DynamicPC::backwardSwitch(Integer i) {

    std::ostringstream oss;
    oss << this->oper_cnt << ": bw_sw " << i;
    if (print_simp_filt_op) { std::cout << oss.str() << this->oper_end_str; }
    if (out_filt_pers) { oper_fout << oss.str() << std::endl; }

    if (!run_update) {
        assert(i > 0 && i < this->dzz.m());
        assert(this->dzz.filt[i-1].op == DEL_OP && this->dzz.filt[i].op == DEL_OP);
        assert( !this->tot_cplx.isSimpFaceOf(
            this->dzz.filt[i].s_id, this->dzz.filt[i-1].s_id) );

        std::swap(this->dzz.filt[i-1], this->dzz.filt[i]);

    } else {
        auto start = std::chrono::high_resolution_clock::now();
        this->dzz.backwardSwitch(i);
        auto end = std::chrono::high_resolution_clock::now();
        this->up_timing += end - start;
    }

    this->runFZZ();

    ++ this->bw_sw_cnt;
    ++ this->oper_cnt;

    this->afterSimplexOper();
}

void DynamicPC::outwardSwitch(
    Integer i, 
    std::vector<Integer> &vanish_intvs, 
    std::vector<Integer> &new_intvs) {

    std::ostringstream oss;
    oss << this->oper_cnt << ": ow_sw " << i;
    if (print_simp_filt_op) { std::cout << oss.str() << this->oper_end_str; }
    if (out_filt_pers) { oper_fout << oss.str() << std::endl; }

    if (!run_update) {
        assert(i > 0 && i < this->dzz.m());
        assert(this->dzz.filt[i-1].op == ADD_OP && this->dzz.filt[i].op == DEL_OP);
        assert(this->dzz.filt[i-1].s_id != this->dzz.filt[i].s_id);

        std::swap(this->dzz.filt[i-1], this->dzz.filt[i]);

    } else {
        auto start = std::chrono::high_resolution_clock::now();
        this->dzz.outwardSwitch(i, vanish_intvs, new_intvs);
        auto end = std::chrono::high_resolution_clock::now();
        this->up_timing += end - start;
    }

    this->runFZZ();

    ++ this->ow_sw_cnt;
    ++ this->oper_cnt;

    this->afterSimplexOper();
}

void DynamicPC::inwardContraction(
    Integer i, 
    std::vector<Integer> &vanish_intvs, 
    std::vector<Integer> &new_intvs) {

    std::ostringstream oss;
    oss << this->oper_cnt << ": iw_con " << i;
    if (print_simp_filt_op) { std::cout << oss.str() << this->oper_end_str; }
    if (out_filt_pers) { oper_fout << oss.str() << std::endl; }

    if (!run_update) {
        assert(i > 0 && i < this->dzz.m());
        assert(this->dzz.filt[i-1].op == ADD_OP && this->dzz.filt[i].op == DEL_OP);
        assert(this->dzz.filt[i-1].s_id == this->dzz.filt[i].s_id);

        auto filt_iter_del_beg = this->dzz.filt.begin()+i-1;
        this->dzz.filt.erase(filt_iter_del_beg, filt_iter_del_beg+2);

    } else {
        auto start = std::chrono::high_resolution_clock::now();
        this->dzz.inwardContraction(i, vanish_intvs, new_intvs);
        auto end = std::chrono::high_resolution_clock::now();
        this->up_timing += end - start;
    }

    this->runFZZ();

    ++ this->iw_con_cnt;
    ++ this->oper_cnt;

    this->afterSimplexOper();
}

void DynamicPC::outwardExpansion(
    Integer i, 
    Integer sigma,
    std::vector<Integer> &vanish_intvs, 
    std::vector<Integer> &new_intvs) {

    std::ostringstream oss;
    oss << this->oper_cnt << ": ow_exp " << i << ", " 
        << containerToStr(this->tot_cplx.simplices[sigma], " ");
    if (print_simp_filt_op) { std::cout << oss.str() << this->oper_end_str; }
    if (out_filt_pers) { oper_fout << oss.str() << std::endl; }

    if (!run_update) {
        assert(i > 0 && i <= this->dzz.m());

        std::vector<SimpOp> ops { {DEL_OP, sigma}, {ADD_OP, sigma} };
        this->dzz.filt.insert(this->dzz.filt.begin()+i-1, ops.begin(), ops.end());

    } else {
        auto start = std::chrono::high_resolution_clock::now();
        this->dzz.outwardExpansion(i, sigma, vanish_intvs, new_intvs);
        auto end = std::chrono::high_resolution_clock::now();
        this->up_timing += end - start;
    }

    this->runFZZ();

    if (this->dzz.filt.size() > this->max_simp_filt_len) 
        { this->max_simp_filt_len = this->dzz.filt.size(); }

    ++ this->ow_exp_cnt;
    ++ this->oper_cnt;

    this->afterSimplexOper();
}

void DynamicPC::afterSimplexOper() {
    // if (print_simp_filt_op) {
    //     std::cout << "simp filt len: " << this->dzz.filt.size() << std::endl;
    //     std::cout << "max simp filt len: " << this->max_simp_filt_len << std::endl;
    //     if (run_update) { std::cout << "up time: " << this->up_timing.count() << std::endl; }
    //     if (run_fzz) { std::cout << "fzz time: " << this->fzz_timing.count() << std::endl; }
    // }
}

bool DynamicPC::isFiltValid() {
    assert(this->simp_filt_pos.size() == this->e_filt.size()+1);

    for (auto j = 0; j < this->e_filt.size(); ++j) {
        const auto start = this->simp_filt_pos[j];
        const auto end = this->simp_filt_pos[j+1];

        assert(start < end);

        for (auto k = start; k < end; ++k) {
            assert(this->dzz.filt[k].op == this->e_filt[j].op);
            assert(this->tot_cplx.isSimpFaceOf(
                this->e_filt[j].e_id, this->dzz.filt[k].s_id));
        }

        if (this->e_filt[j].op == ADD_OP) 
            { assert(this->dzz.filt[start].s_id == this->e_filt[j].e_id); }
        else
            { assert(this->dzz.filt[end-1].s_id == this->e_filt[j].e_id); }
    }

    return true;
}

void DynamicPC::runFZZ() {

    if (!run_fzz) { return; }

    std::vector<Simplex> filt_simp;
    std::vector<AddDelOp> filt_op;
    for (const auto e : this->dzz.filt) {
        filt_simp.emplace_back(this->tot_cplx.simplices[e.s_id]);
        filt_op.emplace_back(e.op);
    }

    // std::cout << std::endl << "runFZZ filt len: " << filt_simp.size() << std::endl;

    FastZigzagPHAT fzz;
    auto start = std::chrono::high_resolution_clock::now();
    fzz.compute(filt_simp, filt_op);
    auto end = std::chrono::high_resolution_clock::now();
    this->fzz_timing += end - start;

    if (!run_update || !check_barcode) { return; }

    std::sort(fzz.persistence.begin(), fzz.persistence.end());

    std::vector< std::tuple<Integer, Integer, Integer> > dzz_pers;
    this->dzz.getPersistence(dzz_pers);
    std::sort(dzz_pers.begin(), dzz_pers.end());

    if (dzz_pers == fzz.persistence) { return; }
    
    printPers(
        tmp_folder+pathSeparator()+std::to_string(oper_cnt)+"_dzz_pers",
        dzz_pers);

    printPers(
        tmp_folder+pathSeparator()+std::to_string(oper_cnt)+"_fzz_pers",
        fzz.persistence);
    
    assert(false);
}

Integer DynamicPC::evalFiltETime(const Decimal dis, Decimal &t1_old, Decimal &t2_old) {

    char direction;
    Integer start_t;
    Integer switch_i = -1;
    // Decimal temp_t1_old = -1, temp_t2_old = -1;

    for (auto i = 0; i < this->e_time.size(); ++i) {
        const EdgeOp &e_op = this->e_filt[i];

        if (switch_i < 0) {
            t1_old = t2_old;
            t2_old = this->e_time[i];
        }

        if (this->e_time[i] == (Decimal)this->min_time) {
            if (dis >= this->edgeDisSample(e_op.e_id, this->min_time)) {
                continue;
            }

            assert(e_op.op == ADD_OP);
            assert(this->edgeDisSample(e_op.e_id, this->min_time) >
                this->edgeDisSample(e_op.e_id, this->min_time+1));
 
            direction = 'r';
            start_t = this->min_time + 1;

        } else if (this->e_time[i] == (Decimal)this->max_time) {
            if (dis >= this->edgeDisSample(e_op.e_id, this->max_time)) {
                continue;
            }

            assert(e_op.op == DEL_OP);
            assert(this->edgeDisSample(e_op.e_id, this->max_time) >
                this->edgeDisSample(e_op.e_id, this->max_time-1));

            direction = 'l';
            start_t = this->max_time - 1;

        } else {
            assert(this->edgeDisInterp(e_op.e_id, this->e_time[i]) > dis);

            Integer lt, rt;
            if (isDecimalInteger(this->e_time[i])) {
                lt = (Integer)this->e_time[i] - 1;
                rt = (Integer)this->e_time[i] + 1;
            } else {
                lt = (Integer)std::floor(this->e_time[i]);
                rt = lt + 1;
            }

            if (this->edgeDisSample(e_op.e_id, lt) < 
                this->edgeDisInterp(e_op.e_id, this->e_time[i])) {

                direction = 'l';
                start_t = lt;
            } else {
                direction = 'r';
                start_t = rt;
            } 
        }

        // if (direction == 'l') {
        //     Integer t;

        //     for (t = start_t; t >= this->min_time && 
        //         this->edgeDisSample(e_op.e_id, t) > dis; --t) {

        //         if (this->edgeDisSample(e_op.e_id, t) > 
        //             this->edgeDisSample(e_op.e_id, t+1)) {

        //             std::cout << "ERR: this->edgeDisSample(e_op.e_id, t) > "
        //                 "this->edgeDisSample(e_op.e_id, t+1)" << std::endl;
        //             exit(-1);
        //         }
        //     }

        //     if (t == this->min_time) {
        //         std::cout << "ERR: t == this->min_time" << std::endl;
        //         exit(-1);
        //     }

        //     this->e_time[i] = interpVal(this->edgeDisSample(e_op.e_id, t), (Decimal)t,
        //         this->edgeDisSample(e_op.e_id, t+1), (Decimal)(t+1), dis);
        // } else {

        //     Integer t;

        //     for (t = start_t; t <= this->max_time && 
        //         this->edgeDisSample(e_op.e_id, t) > dis; ++t) {

        //         if (this->edgeDisSample(e_op.e_id, t) > 
        //             this->edgeDisSample(e_op.e_id, t-1)) {

        //             std::cout << "ERR: this->edgeDisSample(e_op.e_id, t) > "
        //                 "this->edgeDisSample(e_op.e_id, t-1)" << std::endl;
        //             exit(-1);
        //         }
        //     }

        //     if (t == this->max_time) {
        //         std::cout << "ERR: t == this->max_time" << std::endl;
        //         exit(-1);
        //     }

        //     this->e_time[i] = interpVal(this->edgeDisSample(e_op.e_id, t), (Decimal)t,
        //         this->edgeDisSample(e_op.e_id, t-1), (Decimal)(t-1), dis);
        // }

        // std::cout << i << ": ";
        // this->printEdgeOp(e_op);
        // std::cout << std::endl;

        this->e_time[i] = this->descendFindTime(e_op.e_id, direction, start_t, dis);

        if (i > 0 && this->e_time[i] < this->e_time[i-1]) {
            assert(switch_i < 0);

            switch_i = i-1;
            std::swap(this->e_time[i-1], this->e_time[i]);
        }
    }

    assert(std::is_sorted(this->e_time.begin(), this->e_time.end()));

    return switch_i;
}

Decimal DynamicPC::descendFindTimeNoDir(
    const Integer e_id, const Decimal start_t, const Decimal dis) {

    if (e_id == this->LEFT_MIN_E) { return this->min_time; }
    if (e_id == this->RIGHT_MAX_E) { return this->max_time; }

    if (std::floor(start_t) == start_t) { // start_t is an integer
        auto t = (Integer) start_t;
        
        if (t == this->min_time) {
            if (this->edgeDisSample(e_id, t) <= dis) { return t; }
            else { assert(false); }

            // assert(this->edgeDisSample(e_id, t) > 
            //     this->edgeDisSample(e_id, t+1));

            // return descendFindTime(e_id, 'r', t+1, dis);
        }

        if (t == this->max_time) {
            if (this->edgeDisSample(e_id, t) <= dis) { return t; }
            else { assert(false); }

            // assert(this->edgeDisSample(e_id, t) > 
            //     this->edgeDisSample(e_id, t-1));

            // return descendFindTime(e_id, 'l', t-1, dis);
        }

        assert(this->edgeDisSample(e_id, t) >= dis);

        assert( 
            (this->edgeDisSample(e_id, t-1) < this->edgeDisSample(e_id, t))
            == 
            (this->edgeDisSample(e_id, t) < this->edgeDisSample(e_id, t+1)) );

        if (this->edgeDisSample(e_id, t-1) < this->edgeDisSample(e_id, t)) 
            { return descendFindTime(e_id, 'l', t-1, dis); }
        else 
            { return descendFindTime(e_id, 'r', t+1, dis); }
    }

    {
        assert(this->edgeDisInterp(e_id, start_t) >= dis);

        auto t = std::floor(start_t);

        if (this->edgeDisSample(e_id, t) < this->edgeDisInterp(e_id, start_t)) 
            { return descendFindTime(e_id, 'l', t, dis); }
        else 
            { return descendFindTime(e_id, 'r', t+1, dis); }
    }
}

Decimal DynamicPC::descendFindTime(
    const Integer e_id, const char direction, 
    const Integer start_t, const Decimal dis) {

    if (direction == 'l') {
        Integer t;

        for (t = start_t; t >= this->min_time && 
            this->edgeDisSample(e_id, t) > dis; --t) {

            assert(this->edgeDisSample(e_id, t) <
                this->edgeDisSample(e_id, t+1));
        }

        assert(t >= this->min_time);
        assert(this->edgeDisSample(e_id, t+1) >= dis);

        return interpVal(this->edgeDisSample(e_id, t), (Decimal)t,
            this->edgeDisSample(e_id, t+1), (Decimal)(t+1), dis);
    } else {

        Integer t;

        for (t = start_t; t <= this->max_time && 
            this->edgeDisSample(e_id, t) > dis; ++t) {

            assert(this->edgeDisSample(e_id, t) < 
                this->edgeDisSample(e_id, t-1));
        }

        assert(t <= this->max_time);
        assert(this->edgeDisSample(e_id, t-1) >= dis);

        return interpVal(this->edgeDisSample(e_id, t), (Decimal)t,
            this->edgeDisSample(e_id, t-1), (Decimal)(t-1), dis);
    }
}

Decimal DynamicPC::ascendFindTimeNoDir(
    const Integer e_id, const Decimal start_t, const Decimal dis) {

    if (e_id == this->LEFT_MIN_E) { return this->min_time; }
    if (e_id == this->RIGHT_MAX_E) { return this->max_time; }

    if (std::floor(start_t) == start_t) { // start_t is an integer
        auto t = (Integer) start_t;
        assert(this->edgeDisSample(e_id, t) <= dis);
        
        if (t == this->min_time || t == this->max_time) { return t; }

        assert( 
            (this->edgeDisSample(e_id, t-1) < this->edgeDisSample(e_id, t))
            == 
            (this->edgeDisSample(e_id, t) < this->edgeDisSample(e_id, t+1)) );

        if (this->edgeDisSample(e_id, t-1) > this->edgeDisSample(e_id, t)) 
            { return ascendFindTime(e_id, 'l', t-1, dis); }
        else 
            { return ascendFindTime(e_id, 'r', t+1, dis); }
    }

    {
        assert(this->edgeDisInterp(e_id, start_t) <= dis);

        auto t = std::floor(start_t);

        if (this->edgeDisSample(e_id, t) > this->edgeDisInterp(e_id, start_t)) 
            { return ascendFindTime(e_id, 'l', t, dis); }
        else 
            { return ascendFindTime(e_id, 'r', t+1, dis); }
    }
}

Decimal DynamicPC::ascendFindTime(
    const Integer e_id, const char direction, 
    const Integer start_t, const Decimal dis) {

    if (direction == 'l') {
        Integer t;

        for (t = start_t; t >= this->min_time && 
            this->edgeDisSample(e_id, t) < dis; --t) {

            assert(this->edgeDisSample(e_id, t) >
                this->edgeDisSample(e_id, t+1));
        }

        assert(t >= this->min_time);
        assert(this->edgeDisSample(e_id, t+1) <= dis);

        return interpVal(this->edgeDisSample(e_id, t+1), (Decimal)(t+1),
            this->edgeDisSample(e_id, t), (Decimal)t, dis);
    } else {

        Integer t;

        for (t = start_t; t <= this->max_time && 
            this->edgeDisSample(e_id, t) < dis; ++t) {

            assert(this->edgeDisSample(e_id, t) >
                this->edgeDisSample(e_id, t-1));
        }

        assert(t <= this->max_time);
        assert(this->edgeDisSample(e_id, t-1) <= dis);

        return interpVal(this->edgeDisSample(e_id, t-1), (Decimal)(t-1),
            this->edgeDisSample(e_id, t), (Decimal)t, dis);
    }
}
