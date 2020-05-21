//
//  Solution.hpp
//  EVRP_STD
//
//  Created by MAC on 5/14/20.
//  Copyright © 2020 MAC. All rights reserved.
//

#ifndef Solution_hpp
#define Solution_hpp

#include <stdio.h>
class Solution {
public:
    int *seq_node; // 0 -> 1 -> ... 24
    bool *is_though_stat; // 0 -> 1 -> ... 0
    int num_node; // num_current node in seq_node
    
    double fitness = 0.0;
    double cost = 0.0;
    bool is_feasible = false;
    double over_capacity = 0.0;
    int num_c = 0;
    int num_v = 0;
    
    // meta data
    double *F_CAP;
    double *B_CAP;
    double *F_ENG;
    double *B_ENG;
    
    // method
    void streamlining_seq(int *seq_customers, bool *is_though_station);
    void init_mem_space(int num_vhs, int num_customers);
    void init_data();
    void random_init();
    void Potvin_init(int * list_seeds, double **Distances, double *Cap, double max_cap);
    double compute_cost(double ** Distances, double ** Best_Station_Distances);
    double compute_over_cap(double *Demands, double max_capacity_vh);
    void Potvin_init_test(int * list_seeds, double **Distances, double *Demands, double max_cap);
    void reset_in_select(bool *list);
    
};
#endif /* Solution_hpp */
