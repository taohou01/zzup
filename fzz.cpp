#include "fzz.h"

// phat headers
// wrapper algorithm that computes the persistence pairs of a given boundary matrix using a specified algorithm
#include <phat/compute_persistence_pairs.h>

// main data structure (choice affects performance)
#include <phat/representations/vector_vector.h>
#include <phat/representations/bit_tree_pivot_column.h>

// algorithm (choice affects performance)
#include <phat/algorithms/standard_reduction.h>
#include <phat/algorithms/chunk_reduction.h>
#include <phat/algorithms/row_reduction.h>
#include <phat/algorithms/twist_reduction.h>

void getBoundaryChainPhat(const SimplexIdMap &id_map, 
    const Simplex &simp, std::vector<phat::index> &bound_c) {

    bound_c.clear();

    if (simp.size() <= 1) { return; }

    bound_c.reserve(simp.size());

    Simplex bound_simp(simp.begin()+1, simp.end());
    bound_c.push_back(id_map.at(bound_simp));
    // std::cout << "  " << bound_simp << endl;

    for (Integer i = 0; i < simp.size()-1; ++i) {
        bound_simp[i] = simp[i];
        // std::cout << "  " << bound_simp << endl;
        bound_c.push_back(id_map.at(bound_simp));
    }

    std::sort(bound_c.begin(), bound_c.end());
}

inline Integer getDim(const std::vector<phat::index> &bound_c) {
    if (bound_c.size() == 0) { return 0; }
    return bound_c.size() - 1;
}

void FastZigzagPHAT::compute(
    const std::vector<Simplex> &filt_simp, 
    const std::vector<AddDelOp> &filt_op) {

    auto filt_len = filt_simp.size();
    assert(filt_len % 2 == 0);
    simp_num = filt_len / 2;

    std::vector<phat::index> bound_c;
    // phat::boundary_matrix< phat::vector_vector > bound_chains;
    phat::boundary_matrix< phat::bit_tree_pivot_column > bound_chains;
    bound_chains.set_num_cols(filt_len + 1);

    // Add the Omega vertex for the coning
    bound_chains.set_col(0, bound_c);
    bound_chains.set_dim(0, 0);

    orig_f_add_id.reserve(simp_num);
    orig_f_del_id.reserve(simp_num);

    std::vector<Integer> del_ids;
    del_ids.reserve(simp_num);

    SimplexIdMap *p_id_map = new SimplexIdMap();
    SimplexIdMap &id_map = *p_id_map;

    Integer orig_f_id = 0;
    char op;
    Simplex simp;
    Integer s_id = 1;
    Integer death;

    for (auto i = 0; i < filt_simp.size(); ++i) {
        const Simplex &simp = filt_simp[i];

        if (filt_op[i] == ADD_OP) {
            getBoundaryChainPhat(id_map, simp, bound_c);
            bound_chains.set_col(s_id, bound_c);
            bound_chains.set_dim(s_id, getDim(bound_c));

            // assert(s_id == bound_chains.size()-1);
            id_map[simp] = s_id;
            orig_f_add_id.push_back(orig_f_id);
            s_id ++;
        } else {
            del_ids.push_back(id_map[simp]);
            orig_f_del_id.push_back(orig_f_id);
        }

        orig_f_id ++;
    }

    assert(del_ids.size() == s_id-1);
    delete p_id_map;

    assert(simp_num == del_ids.size());

    std::vector<Integer> cone_sid(simp_num+1);
    Integer dim;

    for (auto del_id_it = del_ids.rbegin(); del_id_it != del_ids.rend(); ++del_id_it) {
        bound_c.clear();
        bound_c.push_back(*del_id_it);

        // std::cout << std::endl << bound_chains[*del_id_it] << std::endl;

        std::vector<phat::index> orig_bound_c;
        bound_chains.get_col(*del_id_it, orig_bound_c);

        if (orig_bound_c.size() == 0) {
            bound_c.push_back(0);
        } else {
            for (auto bsimp : orig_bound_c) {
                // assert(cone_sid[bsimp] >= 0);
                bound_c.push_back(cone_sid[bsimp]);
            }
        }

        std::sort(bound_c.begin(), bound_c.end());
        // std::cout << s_id << ": " << *del_id_it << " " << bound_c << std::endl;

        bound_chains.set_col(s_id, bound_c);
        bound_chains.set_dim(s_id, getDim(bound_c));

        cone_sid[*del_id_it] = s_id;

        s_id ++;
    }

    phat::persistence_pairs pairs;
    phat::compute_persistence_pairs< phat::twist_reduction >( pairs, bound_chains );

    for (phat::index idx = 0; idx < pairs.get_num_pairs(); idx++) {
            Integer b = pairs.get_pair(idx).first;
            Integer d = pairs.get_pair(idx).second - 1;
            Integer p = bound_chains.get_dim(b);

            if (d < simp_num) { mapOrdIntv(b, d); } 
            else { mapRelExtIntv(p, b, d); }

            this->persistence.emplace_back(b, d, p);
    }
}
