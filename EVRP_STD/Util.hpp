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
#include "Move.hpp"
class Util {
public:
    int Best_Route_Education[200];
    int Best_Tabu_Search[200];
    bool FOUND_NEW_BEST = false;
    double FIT_BEST = 0.0;
    int tabu[200][200];
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
    //void update_Best_Search(int **plitted_route, int R1, int R2, int idx1, int idx2, int num_c, int num_v);
    void Interchange10APR(int *seq, double **Distances, int num_c, int num_v, int comb[][2], int loop);
    void Interchange20APR(int *seq, double **Distances, int num_c, int num_v, int comb[][2], int loop);
    void Interchange11APR(int *seq, double **Distances, int num_c, int num_v, int comb[][2], int loop);
    void Interchange21APR(int *seq, double **Distances, int num_c, int num_v, int comb[][2], int loop);
    void Interchange22APR(int *seq, double **Distances, int num_c, int num_v, int comb[][2], int loop);
    void local_search(int *seq, double **Distances, int num_c, int num_v, int num);
    void Tabu_search(int *seq, double **Distances, int num_c, int num_v, double *Demands, double init_finess, int **List_Nearest_Cus, double init_cost, double max_cap, double ALPHA, double max_eng, double eng_consum, int **Best_Stat, double ** Best_Stat_Distances);
    Move CallEvaluate(int *seq, double fitness_Zt, int num_c, int num_v, int **List_Nearest_Cus, double init_cost, double ** Distances, double *Demands, double max_cap, double ALPHA, double max_eng, double eng_consum, int **Best_Stat, double ** Best_Stat_Distances, int IT);
    bool validate_through_stat(int *seq, int begin, int num_c, int num_v, double max_eng, double eng_consum, double **Distances, int **Best_Stat, double **Best_Stat_Distances);
};
#endif /* Util_hpp */
