//
//  Solution.cpp
//  EVRP_STD
//
//  Created by MAC on 5/14/20.
//  Copyright © 2020 MAC. All rights reserved.
//

#include "Solution.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void Solution::init_mem_space(int num_vhs, int num_customers)
{
    seq_node = (int *)malloc((num_vhs + num_customers + 1) * sizeof(int));
    is_though_stat = (bool *)malloc((num_vhs + num_customers + 1) * sizeof(bool));
    num_node = num_vhs + num_customers + 1;
    
    fitness = 0.0;
    cost = 0.0;
    is_feasible = true;
    over_capacity = 0.0;
    num_c = num_customers;
    num_v = num_vhs;
    
    // init for meta data later
}

void Solution::streamlining_seq(int *seq_customers, bool *is_though_station)
{
    int num_item_new_seq = 1;
    seq_node[0] = 0;
    is_though_stat[0] = is_though_station[0];
    for(int s = 1; s < num_c + num_v + 1; s++)
    {
        int prev_node = seq_node[num_item_new_seq - 1];
        int current_node = seq_customers[s];
        if(prev_node != 0 && prev_node < num_c)
        {
            seq_node[num_item_new_seq] = seq_customers[s];
            is_though_stat[num_item_new_seq] = is_though_station[s];
            num_item_new_seq++;
        } else
        {
            if(current_node >= num_c || current_node == 0) continue;
            else
            {
                seq_node[num_item_new_seq] = seq_customers[s];
                is_though_stat[num_item_new_seq] = is_though_station[s];
                num_item_new_seq++;
            }
        }
    }
    num_node = num_item_new_seq;
}

void Solution::random_init()
{
    // random init
    int j, k, rand_through;
    int num_custom_each_route = (int)((num_c - 1) / num_v) + 2;
    
    int temp_route[num_c + num_v + 1];
    bool temp_through_stat[num_c + num_v + 1];
    int rand_j = (int)(rand() % (num_c - 1)) + 1;
    k = 1;
    int route_num = num_c - 1;
    for(j = rand_j; j < num_c; j++)
    {
        if(k % num_custom_each_route == 0)
        {
            route_num++;
            temp_route[k] = route_num;
            rand_through = (int)(rand() % 2);
            if(rand_through == 0)
                temp_through_stat[k] = false;
            else temp_through_stat[k] = true;
            k++;
        }
        temp_route[k] = j;
        rand_through = (int)(rand() % 2);
        if (rand_through == 0)
            temp_through_stat[k] = false;
        else
            temp_through_stat[k] = true;
        k++;
    }
    
    for(j = 1; j < rand_j; j++)
    {
        if(k % num_custom_each_route == 0)
        {
            route_num++;
            temp_route[k] = route_num;
            rand_through = (int)(rand() % 2);
            if(rand_through == 0)
                temp_through_stat[k] = false;
            else temp_through_stat[k] = true;
            k++;
        }
        temp_route[k] = j;
        rand_through = (int)(rand() % 2);
        if(rand_through == 0)
            temp_through_stat[k] = false;
        else temp_through_stat[k] = true;
        k++;
    }
    
    route_num++;
    temp_route[k] = route_num;
    
    // set though station for depot
    rand_through = (int)(rand() % 2);
    if (rand_through == 0)
        temp_through_stat[0] = false;
    else temp_through_stat[0] = true;
    
    temp_route[0] = 0;
    temp_route[num_c + num_v] = 0;
    streamlining_seq(temp_route, temp_through_stat);
}














void Solution::Potvin_init(int *list_seeds, double ** Distances)
{
    int Routes[num_v][num_c]; // a -> 1 -> 2 -> 3 -> ...
    double best_cost_route[num_v][num_c];
    int best_route[num_c];
    int seq_selected[num_c];
    int num_selected = 0;
    double cost_route[num_v];
    double gap[num_c];
    int after_index[num_v][num_c];
    
    // init seq_selected data
    for(int i = 0; i < num_c; i++)
    {
        seq_selected[i] = i;
    }
    // init Routes data
    for(int m = 0; m < num_v; m++)
    {
        for(int n = 0; n < num_c; n++)
        {
            Routes[m][n] = -1;
        }
    }
    
    // insert seed customer to routes
    for(int s = 0; s < num_v; s++)
    {
        int seed = list_seeds[s];
        Routes[s][0] = 2;
        Routes[s][1] = 0;
        Routes[s][2] = seed;
        // put seed index to tail of arr
        //list_seeds[s] = list_seeds[num_c - 1 - num_selected];
        int index_seed = -1;
        for(int k = 0 ; k < num_c - num_selected; k++)
        {
            if(seq_selected[k] == seed) {
                index_seed = k;
                break;
            }
        }
        
        seq_selected[index_seed] = seq_selected[num_c - 1 - num_selected];
        seq_selected[num_c - 1 - num_selected] = seed;
        //list_seeds[num_c - 1 - num_selected] = seed;
        num_selected++;
        cost_route[s] = 2 * Distances[seed][0];
    }
    
    while(num_selected < num_c - 1)
    {
        // calculate gap of unrounded customers
        for(int k = 1; k < num_c - num_selected; k++)
        {
            int cus = seq_selected[k];
            
            for(int s = 0; s < num_v; s++)
            {
                double best_cost = 100000;
                int num_curr_cus = Routes[s][0];
                
                double current_route_cost = cost_route[s];
                for(int i = 1; i <= num_curr_cus; i++)
                {
                    int prev_node = Routes[s][i];
                    if(prev_node == -1) prev_node = 0;
                    int next_node = Routes[s][i + 1];
                    if(next_node == -1) next_node = 0;
                    
                    double new_cost = current_route_cost + (Distances[cus][prev_node] + Distances[cus][next_node] - Distances[prev_node][next_node]);
                    if(new_cost < best_cost){
                        best_cost = new_cost;
                        after_index[s][cus] = i;
                    }
                }
                best_cost_route[s][cus] = best_cost;
                
            }
        }
        
        // calculate gap
        for(int i = 1; i < num_c - num_selected; i++)
        {
            int cus = seq_selected[i];
            double best_cost = best_cost_route[0][cus];
            best_route[cus] = 0;
            double sum_cost = 0.0;
            for(int s = 0; s < num_v; s++)
            {
                double cost_route = best_cost_route[s][cus];
                sum_cost += cost_route;
                if(cost_route < best_cost)
                {
                    best_cost = cost_route;
                    best_route[cus] = s;
                }
            }
            gap[cus] = sum_cost - num_v * best_cost;
        }
        
        // find max gap
        
        int best_max_gap = -1000000000;
        int index_best_gap = 1;
        for(int i = 1; i < num_c - num_selected; i++)
        {
            int cus = seq_selected[i];
            int gap_cus = gap[cus];

            if(best_max_gap < gap_cus){
                best_max_gap = gap_cus;
                index_best_gap = i;
            }
        }
        
        
        int best_cus = seq_selected[index_best_gap];
        int best_route_c = best_route[best_cus];
        seq_selected[index_best_gap] = seq_selected[num_c - 1 - num_selected];
        seq_selected[num_c - 1 - num_selected] = best_cus;
//        printf("\n----------\n");
//        for(int m = 0 ; m < num_c; m++)
//        {
//            if(m == num_c - num_selected) printf("|| %d -> ", seq_selected[m]);
//            else printf("%d -> ", seq_selected[m]);
//        }
        
        num_selected++;
        // insert
        int after = after_index[best_route_c][best_cus];
        int prev_node = Routes[best_route_c][after];
        if(prev_node == -1) prev_node = 0;
        int next_node = Routes[best_route_c][after + 1];
        if(next_node == -1) next_node = 0;
        //// update cost
        cost_route[best_route_c] = cost_route[best_route_c] + (Distances[best_cus][prev_node] + Distances[best_cus][next_node] - Distances[prev_node][next_node]);
        
        int curr_num_best_route = Routes[best_route_c][0];
        for(int k = curr_num_best_route + 1; k > after + 1; k--)
        {
            Routes[best_route_c][k] = Routes[best_route_c][k - 1];
        }
        Routes[best_route_c][after + 1] = best_cus;
        Routes[best_route_c][0] = curr_num_best_route + 1;
    }
    
    // concat to solution
    int current_node = 1;
    int current_route = num_c;
    seq_node[0] = 0;
    for(int i = 0; i < num_v; i++)
    {
        int num_node = Routes[i][0];
        for(int j = 2; j <= num_node; j++)
        {
            seq_node[current_node] = Routes[i][j];
            current_node++;
        }
        seq_node[current_node] = current_route;
        current_route++;
        current_node++;
    }
}

double Solution::compute_cost(double **Distances, double **Best_Station_Distances)
{
    double sum_cost = 0.0;
    for(int i = 0; i < num_c + num_v - 1; i++)
    {
        // current node
        int node = seq_node[i];
        if(node >= num_c) node = 0;
        
        // next node
        int next_node = seq_node[i + 1];
        if(next_node >= num_c) next_node = 0;
        
        if(node != next_node)
        {
            if(is_though_stat[i])
            {
                sum_cost += Best_Station_Distances[node][next_node];
            } else sum_cost += Distances[node][next_node];
        }
    }
    cost = sum_cost;
    return sum_cost;
}

double Solution::compute_over_cap(double *Demands, double max_capacity_vh)
{
    int i;
    int begin_route_num = num_c;
    double sum_capacity_over = 0.0;
    double sum_capacity_route = 0.0;
    
    for(i = 1; i < num_c + num_v; i++)
    {
        int node = seq_node[i];
        if(node < begin_route_num)
        {
            sum_capacity_route = Demands[node];
        } else{
            double temp_over = max_capacity_vh - sum_capacity_route;
            sum_capacity_route = 0.0;
            if(temp_over < 0)
            {
                sum_capacity_over += temp_over;
            }
        }
    }
    if(sum_capacity_route < 0.0)
    {
        is_feasible = false;
        over_capacity = -sum_capacity_route;
        
    }
    
    return sum_capacity_over;
}
