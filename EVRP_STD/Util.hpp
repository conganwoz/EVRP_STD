//
//  Util.hpp
//  EVRP_STD
//
//  Created by MAC on 5/15/20.
//  Copyright © 2020 MAC. All rights reserved.
//

#ifndef Util_hpp
#define Util_hpp

#include <stdio.h>
#include "Solution.hpp"
class Util {
public:
    int Best_Route_Education[200];
    double compute_H_value(Solution *Sols, int num_sols);
    double compute_AVG_q(Solution *Sols, int num_sols);
    double compute_alpha(double H, double avg_q);
    bool check_selected_customer(int *seq_customer, int begin, int end, int value);
    void Permutation_Order_1(Solution sol1, Solution sol2, int *child_1, int num_node);
    void Cycle(Solution sol1, Solution sol2, int * child, int num_node);
    double compute_cost_sol(int *sol_seq, double **Distances, int num_c, int num_v);
    bool choose_better_sol(int *sol_seq1, int *sol_seq2, double **Distances, int num_c, int num_v);
    
    int Rand(int i, int j);
    
    // Education phase
    void set_init_best(int *sol, int num_node);
    bool find_through_station(int *sol, bool *through_station, int num_c, int num_v, double max_eng, double eng_consum, double **Distances, int **best_stat, double **best_stat_distances);
    double dist_comsum(double distance, double eng_consum);
    
    // Operators
    void two_exchange_education(int *sol, double ** Distances, int num_c, int num_v);
    void or_exhcange(int *sol, double **Distances, int num_c, int num_v);
    void cross_exchange(int *sol, double ** Distances, int num_c, int num_v);
    
    bool Education(int *seq_node, bool * though_station, double ** Distances, int num_c, int num_v, double max_eng, double eng_consum, int **best_stat, double **best_stat_distances);
    
    // save sol to pool
    void save_sol_to_pool(int *seq, int * pool_seq, int num_node);
    
    //optimize local search for best sol
    void Interchange10APR(int *seq, double **Distances, int num_c, int num_v, int comb[][2], int loop);
    void Interchange20APR(int *seq, double **Distances, int num_c, int num_v);
    void Interchange11APR(int *seq, double **Distances, int num_c, int num_v);
    void Interchange21APR(int *seq, double **Distances, int num_c, int num_v);
    void Interchange22APR(int *seq, double **Distances, int num_c, int num_v);
    void local_search(int *seq, double **Distances, int num_c, int num_v, int num);
};
#endif /* Util_hpp */
