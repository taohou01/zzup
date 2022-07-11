#ifndef _DPC_H_
#define _DPC_H_

#include <iomanip>
#include <unordered_set>
#include <unordered_map>
#include <chrono>

#include "common.h"
#include "utils.h"
#include "total_complex.h"
#include "dynamic_zigzag.h"

class VineInterval {
public:
    Integer dimension;
    Integer birth_e;
    Decimal birth_time;
    Integer death_e;
    Decimal death_time;
};

class DynamicPC {
public:
    DynamicPC() : dzz(&tot_cplx), up_timing(0.0), fzz_timing(0.0) {
        if (out_filt_pers) { oper_fout.open(tmp_folder + pathSeparator() + "oper"); }
    }

    void init(
        const std::string &posfname, Integer _max_dim, 
        Integer _min_time, Integer _max_time);

    void genEvent();

    void initEdgeFilt();

    void initPersistence();

    void travEvent();


    void localMin(Integer evt_id);

    void localMax(Integer evt_id);

    void incCrossing(Integer evt_id);

    void decCrossing(Integer evt_id);

    void oppoCrossing(Integer evt_id);


    Integer evalFiltETime(const Decimal dis, Decimal &t1_old, Decimal &t2_old);

    void addVanishNewIntervals(
        std::vector<Integer> &vanish_intvs, 
        std::unordered_set<Integer> &new_intvs,
        const std::vector<Integer> &vanish_intvs_local, 
        const std::vector<Integer> &new_intvs_local);


    void localMinExec(
        Integer i, 
        std::vector<Integer> &vanish_intvs, 
        std::unordered_set<Integer> &new_intvs);

    void localMaxExec(
        Integer i, Integer e_id,
        std::vector<Integer> &vanish_intvs, 
        std::unordered_set<Integer> &new_intvs);

    void incCrossingExec(Integer i);

    void decCrossingExec(Integer i);

    void oppoCrossingExec(
        Integer i,
        std::vector<Integer> &vanish_intvs, 
        std::unordered_set<Integer> &new_intvs);


    bool getEdgeFiltInterval(
        const Integer rid, 
        Integer &birth_e_id, Integer &death_e_id,
        Decimal &birth_time, Decimal &death_time);

    void updateVines(const Decimal evt_id,
        const std::vector<Integer> &vanish_intvs, 
        const std::unordered_set<Integer> &new_intvs);

    void endVineInterval(const Integer rid, const Decimal evt_dis);


    Decimal descendFindTimeNoDir(
        const Integer e_id, const Decimal start_t, const Decimal dis);

    Decimal descendFindTime(
        const Integer e_id, const char direction, 
        const Integer start_t, const Decimal dis);

    Decimal ascendFindTimeNoDir(
        const Integer e_id, const Decimal start_t, const Decimal dis);

    Decimal ascendFindTime(
        const Integer e_id, const char direction, 
        const Integer start_t, const Decimal dis);


    Decimal edgeDisSample(const Integer e_id, const Integer t) { 
        assert(e_id >= this->tot_cplx.dimStart(1) && e_id < this->tot_cplx.dimEnd(1));
        assert(t >= this->min_time && t <= this->max_time);
        return this->edge_dis.at(e_id - this->num_pts).at(t); 
    }

    Decimal edgeDisInterp(const Integer e_id, const Decimal t) {
        assert(e_id >= this->tot_cplx.dimStart(1) && e_id < this->tot_cplx.dimEnd(1));
        assert(t >= this->min_time && t <= this->max_time);

        if (std::floor(t) == t) { return edgeDisSample(e_id, (Integer)t); }

        auto lt = static_cast<Integer>(std::floor(t));
        auto lval = this->edgeDisSample(e_id, lt);
        auto rval = this->edgeDisSample(e_id, lt+1);

        return lval + (t - lt) * (rval - lval);
    }

    void printEvent(const Integer evt_id) {
        const char* evt_name[] = 
        {"l_min", "r_min", "l_max", "r_max", "min", "max", "inc", "dec", "oppo"};

        const auto& evt = this->events[evt_id];

        std::cout << std::setprecision(std::numeric_limits<Decimal>::digits10) 
            << evt_id << " " << evt.d  << " " 
            << evt.t << " " << evt_name[evt.type] 
            << " (" << tot_cplx.simplices[evt.e1][0] << " "
            << tot_cplx.simplices[evt.e1][1] << ")";

        if (evt.e2 >= 0) {
            std::cout << " (" << tot_cplx.simplices[evt.e2][0] << " "
                << tot_cplx.simplices[evt.e2][1] << ")";
        }

        std::cout << std::endl;
    }

    void printEdgeOp(const EdgeOp e_op) {
        if (e_op.op == ADD_OP) {
            std::cout << "a";
        } else {
            std::cout << "d";
        }

        std::cout << " (" << tot_cplx.simplices[e_op.e_id][0] << " "
            << tot_cplx.simplices[e_op.e_id][1] << ")";
    }

    void printEFilt() {
        for (auto i = 0; i < this->e_filt.size(); ++i) {
            std::cout << i << ": " << this->e_time[i] << " ";
            this->printEdgeOp(this->e_filt[i]);
            std::cout << std::endl;
        }

        std::cout << std::endl;
    }

    bool isFiltValid();

    void runFZZ();

    void forwardSwitch(Integer i);

    void backwardSwitch(Integer i);

    void outwardSwitch(
        Integer i,
        std::vector<Integer> &vanish_intvs, 
        std::vector<Integer> &new_intvs);

    void inwardContraction(
        Integer i, 
        std::vector<Integer> &vanish_intvs, 
        std::vector<Integer> &new_intvs);

    void outwardExpansion(
        Integer i, 
        Integer sigma,
        std::vector<Integer> &vanish_intvs, 
        std::vector<Integer> &new_intvs);
    
    void afterSimplexOper();

    // Decimal getEdgeTimeByPos(Integer e_pos) {
    //     if (e_pos == this->LEFT_MIN_E) { return this->min_time; }
    //     if (e_pos == this->RIGHT_MAX_E) { return this->max_time; }
    //     return this->e_time.at(e_pos);
    // }

public:
    Integer max_dim;
    Integer min_time, max_time;
    Integer num_pts;
    std::vector<std::vector<Decimal>> edge_dis;

    TotalComplex tot_cplx;

    std::vector<DpcEvent> events;

    std::vector<EdgeOp> e_filt;
    std::vector<Decimal> e_time;
    std::vector<Integer> simp_filt_pos;

    // std::vector<Integer> e_filt_pos;

    DynamicZigzag dzz;

    Integer oper_cnt = 0;

    Integer fw_sw_cnt = 0;
    Integer bw_sw_cnt = 0;
    Integer ow_sw_cnt = 0;
    Integer iw_con_cnt = 0;
    Integer ow_exp_cnt = 0;

    Integer max_e_filt_len = 0;
    Integer max_simp_filt_len = 0;

    std::unordered_map<Integer, VineInterval> vine_intervals;

    std::ofstream oper_fout;

    std::ofstream vines_fout;

    // The time before any edges are added
    const Integer LEFT_MIN_E = -1;
    // The time after all edges are deleted
    const Integer RIGHT_MAX_E = -2;
    
    bool print_e_filt_op = false;

    bool print_simp_filt_op = false;

    // const bool print_timing_per_filt_op = true;

    const bool print_vines = false;

    const bool out_filt_pers = false;

    const bool exec_simp_filt_op = true;

    const bool run_update = true;

    bool run_fzz = false;

    const bool check_barcode = false;

    // const bool check_filt = true;

    const std::string tmp_folder { "tmp_filt_pers" };
    
    const std::string oper_end_str { "\n" };
    // const std::string oper_end_str { " | \r" };

    std::chrono::duration<double> up_timing;
    std::chrono::duration<double> fzz_timing;
};

#endif
