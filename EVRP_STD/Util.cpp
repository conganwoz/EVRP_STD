//
//  Util.cpp
//  EVRP_STD
//
//  Created by MAC on 5/15/20.
//  Copyright © 2020 MAC. All rights reserved.
//

#include "Util.hpp"
#include "Move.hpp"
#include <stdlib.h>
double Util::compute_H_value(Solution *Sols, int num_sols)
{
    int i;
    bool check = true;
    double h_s_worst = -100000.0;
    double h_s_best = 100000.0;
    for(i = 0; i < num_sols; i++)
    {
        if(Sols[i].is_feasible)
        {
            check = false;
            if(h_s_best > Sols[i].cost)
                h_s_best = Sols[i].cost;
        } else
        {
            if(h_s_worst < Sols[i].cost)
            {
                h_s_worst = Sols[i].cost;
            }
        }
    }
    
    if(check)
        return h_s_worst;
    else
        return h_s_best;
}

int Util::Rand(int i, int j)
{
    int iInterval = j - i + 1;
    return rand() % iInterval + i;
}

double Util::compute_AVG_q(Solution *Sols, int num_sols)
{
    int i;
    double sum_q = 0.0;
    for(i = 0; i < num_sols; i++)
    {
        sum_q += Sols[i].over_capacity;
    }
    return sum_q / num_sols;
}

double Util::compute_alpha(double H, double avg_q)
{
    return H / avg_q;
}

bool Util::check_selected_customer(int *seq_customer, int begin, int end, int value)
{
    int i;
    for(i = begin; i <= end; i++)
    {
        if(seq_customer[i] == value)
            return true;
    }
    return false;
}

void Util::Permutation_Order_1(Solution sol1, Solution sol2, int *child, int num_node)
{
    int i, begin_rand, end_rand;
    int seq_choosen[num_node];
    for(i = 0; i < num_node; i++)
    {
        seq_choosen[i] = i;
    }
    
    begin_rand = (int)(rand() % (num_node - 1)) + 1;
    int temp = seq_choosen[begin_rand];
    seq_choosen[begin_rand] = seq_choosen[num_node - 1];
    seq_choosen[num_node - 1] = temp;
    begin_rand = seq_choosen[num_node - 1];
    
    end_rand = (int)(rand() % (num_node - 2)) + 1;
    temp = seq_choosen[end_rand];
    seq_choosen[end_rand] = seq_choosen[num_node - 2];
    seq_choosen[num_node - 2] = temp;
    end_rand = seq_choosen[num_node - 2];
    
    if(begin_rand > end_rand)
    {
        temp = begin_rand;
        begin_rand = end_rand;
        end_rand = temp;
    }
    
    for(i = begin_rand; i <= end_rand; i++)
    {
        child[i] = sol1.seq_node[i];
    }
    
    int j = 1;
    for(i = 1; i < begin_rand; i++)
    {
        for(; j < num_node; j++)
        {
            if(!check_selected_customer(sol1.seq_node, begin_rand, end_rand, sol2.seq_node[j]))
            {
                child[i] = sol2.seq_node[j];
                j++;
                break;
            }
        }
    }
    
    for(i = end_rand + 1; i < num_node; i++)
    {
        for(; j < num_node; j++)
        {
            if(!check_selected_customer(sol1.seq_node, begin_rand, end_rand, sol2.seq_node[j])){
                child[i] = sol2.seq_node[j];
                j++;
                break;
            }
        }
    }
    child[0] = 0;
}

void Util::Permutation_Order_1_v1(Solution sol1, Solution sol2, int *child, int num_node)
{
    int i, begin_rand, end_rand;
    int seq_choosen[num_node];
    for(i = 0; i < num_node; i++)
    {
        seq_choosen[i] = i;
    }
    
    begin_rand = (int)(rand() % (num_node - 1)) + 1;
    int temp = seq_choosen[begin_rand];
    seq_choosen[begin_rand] = seq_choosen[num_node - 1];
    seq_choosen[num_node - 1] = temp;
    begin_rand = seq_choosen[num_node - 1];
    
    end_rand = (int)(rand() % (num_node - 2)) + 1;
    temp = seq_choosen[end_rand];
    seq_choosen[end_rand] = seq_choosen[num_node - 2];
    seq_choosen[num_node - 2] = temp;
    end_rand = seq_choosen[num_node - 2];
    
    if(begin_rand > end_rand)
    {
        temp = begin_rand;
        begin_rand = end_rand;
        end_rand = temp;
    }
    
    for(i = begin_rand; i <= end_rand; i++)
    {
        child[i] = sol1.seq_node[i];
    }
    
    // create rotate seq
    int temp_sol2[200];
    int index = 0;
    for(int j = end_rand + 1; j < num_node; j++)
    {
        index++;
        temp_sol2[index] = sol2.seq_node[j];
    }
    
    for(int j = 1; j <= end_rand; j++)
    {
        index++;
        temp_sol2[index] = sol2.seq_node[j];
    }
    
    temp_sol2[0] = 0;
    
    printf("\n INSPECR Order 1 v1:\n");
    for(int k = 0; k < num_node; k++)
        printf("%d -> ", temp_sol2[k]);
    
    int j = 1;
    for(i = 1; i < begin_rand; i++)
    {
        for(; j < num_node; j++)
        {
            int node_2 = temp_sol2[j];
            if(!check_selected_customer(sol1.seq_node, begin_rand, end_rand, node_2))
            {
                child[i] = node_2;
                j++;
                break;
            }
        }
    }
    
    for(i = end_rand + 1; i < num_node; i++)
    {
        for(; j < num_node; j++)
        {
            int node_2 = temp_sol2[j];
            if(!check_selected_customer(sol1.seq_node, begin_rand, end_rand, node_2))
            {
                child[i] = node_2;
                j++;
                break;
            }
        }
    }
    child[0] = 0;
}



void Util::Cycle(Solution sol1, Solution sol2, int *child, int num_node)
{
//    printf("\n--------\n");
//    for(int m = 0; m < num_node; m++)
//    {
//        printf("%d ->", sol1.seq_node[m]);
//    }
//    printf("\n");
//    for(int n = 0; n < num_node; n++)
//    {
//        printf("%d ->", sol2.seq_node[n]);
//    }
//    printf("\n");
    bool is_selected[num_node];
    for(int k = 0; k < num_node; k++)
    {
        is_selected[k] = false;
    }
    
    int num_selected = 0;
    bool is_get_from_p1 = true;
    int idx_choose_from_p1 = 0;
//    int prev_index = 0;
    int init_index = 0;
    
    while (num_selected <= num_node) {
        //printf("\n %d", idx_choose_from_p1);
        if(is_get_from_p1)
        {
            child[idx_choose_from_p1] = sol1.seq_node[idx_choose_from_p1];
        } else {
            child[idx_choose_from_p1] = sol2.seq_node[idx_choose_from_p1];
        }
        
        is_selected[idx_choose_from_p1] = true;
        num_selected++;
        
        int get_from_p2 = sol2.seq_node[idx_choose_from_p1];
        if(get_from_p2 == sol1.seq_node[init_index])
        {
            is_get_from_p1 = !is_get_from_p1;
            for(int k = 0; k < num_node; k++)
            {
                if(!is_selected[k])
                {
                    init_index = k;
                    idx_choose_from_p1 = k;
                    break;
                }
            }
            continue;
        }
        
        // search on p1
        for(int s = 0; s < num_node; s++)
        {
            if(get_from_p2 == sol1.seq_node[s])
            {
                idx_choose_from_p1 = s;
                break;
            }
        }
    }
}

double Util::compute_cost_sol(int *sol_seq, double **Distances, int num_c, int num_v)
{
    double cost = 0.0;
    for(int i = 0; i < num_c + num_v - 1; i++)
    {
        int node = sol_seq[i];
        if(node >= num_c) node = 0;
        
        int next_node = sol_seq[i + 1];
        if(next_node >= num_c) next_node = 0;
        
        if(node != next_node) cost += Distances[node][next_node];
    }
    return cost;
}

bool Util::choose_better_sol(int *sol_seq1, int *sol_seq2, double **Distances, int num_c, int num_v)
{
    double cost1 = compute_cost_sol(sol_seq1, Distances, num_c, num_v);
    double cost2 = compute_cost_sol(sol_seq2, Distances, num_c, num_v);
    
    if(cost1 > cost2) return true;
    else return  false;
}

void Util::set_init_best(int *sol, int num_node)
{
    int i;
    for(i = 0; i < num_node; i++)
    {
        Best_Route_Education[i] = sol[i];
    }
}

void Util::two_exchange_education(int *sol, double **Distances, int num_c, int num_v)
{
    int begin_route_num = num_c;
    int s, begin = 1, end = 1, i, j;
    int parent[num_c + num_v];
    for(s = 0; s < num_c + num_v; s++)
    {
        parent[s] = sol[s];
        Best_Route_Education[s] = sol[s];
    }
    double best_cost = compute_cost_sol(sol, Distances, num_c, num_v);
    double init_cost = best_cost;
    while(parent[begin] < (num_c + num_v))
    {
        begin = end + 1;
        while (parent[begin] >= begin_route_num) {
            begin++;
        }
        end = begin + 1;
        while (parent[end] < begin_route_num) {
            end++;
        }
        
        if((parent[begin] > num_c + num_v) || (parent[end] > num_c + num_v)) break;
        if(begin >= num_c + num_v || end >= num_c + num_v) break;
        
        for(i = begin; i < end - 2; i++)
        {
            for(j = i + 2; j < end - 1; j++)
            {
                int i_1 = i + 1;
                int j_1 = j + 1;
                int node_i = parent[i];
                int node_i_1 = parent[i_1];
                int node_j = parent[j];
                int node_j_1 = parent[j_1];
                double route_cost = init_cost + ((Distances[node_i][node_j] + Distances[node_i_1][node_j_1]) - (Distances[node_i][node_i_1] + Distances[node_j][node_j_1]));
                
                if(route_cost < best_cost)
                {
                    best_cost = route_cost;
                    for(int t = 0; t < num_c + num_v; t++)
                    {
                        Best_Route_Education[t] = parent[t];
                    }
                    
                    int temp_node = Best_Route_Education[i+1];
                    Best_Route_Education[i + 1] = Best_Route_Education[j];
                    Best_Route_Education[j] = temp_node;
                    // revert route seq
                    for(int m = i + 2; m < (i + j + 1) / 2; m++)
                    {
                        int temp_node = Best_Route_Education[m];
                        Best_Route_Education[m] = Best_Route_Education[j - 1 - (m - i - 2)];
                        Best_Route_Education[j - 1 - (m - i - 2)] =temp_node;
                    }
                }
            }
        }
        
        if(end == num_c + num_v) break;
    }
    
    // update route
    for(int s = 0; s < num_c + num_v; s++)
    {
        sol[s] = Best_Route_Education[s];
    }
}

void Util::or_exhcange(int *sol, double **Distances, int num_c, int num_v)
{
    int begin_route_num = num_c;
    int s, begin = 0, end = 0, i, j;
    int parent[num_c + num_v];
    for(s = 0; s < num_c + num_v; s++)
    {
        parent[s] = sol[s];
        Best_Route_Education[s] = sol[s];
    }
    
    double best_cost = compute_cost_sol(sol, Distances, num_c, num_v);
    double init_cost = best_cost;
    while (parent[begin] < begin_route_num) {
        begin = end + 1;
        while (parent[begin] >= begin_route_num) {
            begin++;
        }
        end = begin + 1;
        while (parent[end] < begin_route_num) {
            end++;
        }
        
        if((parent[begin] > num_c + num_v) || (parent[end] > num_c + num_v)) break;
        
        for(i = begin; i < end - 5; i++)
        {
            for(j = i + 4; j < end; j++)
            {
                int j_1 = j + 1;
                int i1_sub_1 = i - 1;
                int i1 = i;
                int i2 = i1 + 1;
                int i2_add_1 = i2 + 1;
                
                int node_i1 = parent[i1];
                int node_i1_sub_1 = parent[i1_sub_1];
                int node_i2 = parent[i2];
                int node_i2_add_1 = parent[i2_add_1];
                int node_j = parent[j];
                int node_j_add_1 = parent[j_1];
                
                double new_cost = init_cost + ((Distances[node_i1_sub_1][node_i2_add_1] + Distances[node_i2][node_j] + Distances[node_i1][node_j_add_1]) - (Distances[node_i1_sub_1][node_i1] + Distances[node_i2][node_i2_add_1] + Distances[node_j_add_1][node_j]));
                
                if(new_cost < best_cost)
                {
                    best_cost = new_cost;
                    // copy to best_route_education
                    int k = 0;
                    for(k = 0; k <= i1_sub_1; k++)
                    {
                        Best_Route_Education[k] = sol[k];
                    }
                    Best_Route_Education[k] = sol[i2_add_1];
                    for(int k1 = i2_add_1 + 1; k1 <= j; k1++)
                    {
                        k++;
                        Best_Route_Education[k] = sol[k1];
                    }
                    for(int k1 = i1; k1 <= i2; k1++)
                    {
                        k++;
                        Best_Route_Education[k] = sol[k1];
                    }
                    k++;
                    Best_Route_Education[k] = sol[j_1];
                    for(int k1 = j_1 + 1; k1 < num_c + num_v; k1++)
                    {
                        k++;
                        Best_Route_Education[k] = sol[k1];
                    }
                    Best_Route_Education[0] = 0;
                }
            }
        }
        if(end == num_c + num_v) break;
    }
    // update route
    for(int s = 0; s < num_c + num_v; s++)
    {
        sol[s] = Best_Route_Education[s];
    }
}

void Util::cross_exchange(int *sol, double **Distances, int num_c, int num_v)
{
    int temp_route_split[num_v][num_c + num_v];
    for(int i = 0; i < num_v; i++)
    {
        for(int j = 0; j < num_v + num_c; j++)
        {
            temp_route_split[i][j] = num_c + num_v + 1;
        }
    }
    
    int s;
    for(s = 0; s < num_c + num_v; s++)
    {
        Best_Route_Education[s] = sol[s];
    }
    
    double best_cost = compute_cost_sol(sol, Distances, num_c, num_v);
    double init_cost = best_cost;
    
    // split route
    int index_route = 0;
    int index_node = 0;
    for(int i = 1; i < num_c + num_v; i++)
    {
        int node = sol[i];
        if(node < num_c)
        {
            index_node++;
            temp_route_split[index_route][index_node] = node;
        } else
        {
            index_node++;
            temp_route_split[index_route][index_node] = node;
            temp_route_split[index_route][0] = index_node;
            if(i < num_c + num_v - 1)
            {
                index_route++;
                index_node = 0;
                if(index_route > num_v)
                    return;
            }
        }
    }
    
    // cross exchange
    for(int i = 0; i < index_route; i++)
    {
        for(int j = i + 1; j < index_route; j++)
        {
            int numc_route_i = temp_route_split[i][0];
            int numc_route_j = temp_route_split[j][0];
            
            if(numc_route_i >= 5 && numc_route_j >= 5)
            {
                for(int X1 = 1; X1 <= numc_route_i - 4; X1++)
                {
                    for(int Y1 = X1 + 2; Y1 <= numc_route_i - 2; Y1++)
                    {
                        int X1_P = X1 + 1;
                        int Y1_P = Y1 + 1;
                        
                        for(int X2 = 1; X2 < numc_route_j - 4; X2++)
                        {
                            for(int Y2 = X2 + 2; Y2 < numc_route_j - 2; Y2++)
                            {
                                int X2_P = X2 + 1;
                                int Y2_P = Y2 + 1;
                                
                                
                                int node_X1 = temp_route_split[i][X1];
                                int node_X1_P = temp_route_split[i][X1_P];
                                int node_Y1 = temp_route_split[i][Y1];
                                int node_Y1_P = temp_route_split[i][Y1_P];
                                
                                int node_X2 = temp_route_split[j][X2];
                                int node_X2_P = temp_route_split[j][X2_P];
                                int node_Y2 = temp_route_split[j][Y2];
                                int node_Y2_P = temp_route_split[j][Y2_P];
                                
                                double new_cost = init_cost + ((Distances[node_X1][node_X2_P] + Distances[node_Y2][node_Y1_P] + Distances[node_X2][node_X1_P] + Distances[node_Y1][node_Y2_P]) - (Distances[node_X1][node_X1_P] + Distances[node_Y1][node_Y1_P] + Distances[node_X2][node_X2_P] + Distances[node_Y2][node_Y2_P]));
                                
                                if(new_cost < init_cost)
                                {
                                    // build route
                                    int index = 0;
                                    Best_Route_Education[0] = 0;
                                    for(int i1 = 1; i1 <= X1; i1++)
                                    {
                                        index++;
                                        if(index > num_c + num_v) break;
                                        Best_Route_Education[index] = temp_route_split[i][i1];
                                    }
                                    for(int i1 = X2_P; i1 <= Y2; i1++)
                                    {
                                        index++;
                                        if(index > num_c + num_v) break;
                                        Best_Route_Education[index] = temp_route_split[j][i1];
                                    }
                                    for(int i1 = Y1_P; i1 <= numc_route_i; i1++)
                                    {
                                        index++;
                                        if(index > num_c + num_v) break;
                                        Best_Route_Education[index] = temp_route_split[i][i1];
                                    }
                                    for(int i1 = 1; i1 <= X2; i1++)
                                    {
                                        index++;
                                        if(index > num_c + num_v) break;
                                        Best_Route_Education[index] = temp_route_split[j][i1];
                                    }
                                    for(int i1 = X1_P; i1 <= Y1; i1++)
                                    {
                                        index++;
                                        if(index > num_c + num_v) break;
                                        Best_Route_Education[index] = temp_route_split[i][i1];
                                    }
                                    for(int i1 = Y2_P; i1 <= numc_route_j; i1++)
                                    {
                                        index++;
                                        if(index > num_c + num_v) break;
                                        Best_Route_Education[index] = temp_route_split[j][i1];
                                    }
                                    
                                    for(int i2 = 0; i2 <= index_route; i2++)
                                    {
                                        if(i2 != i && i2 != j)
                                        {
                                            int numc = temp_route_split[i2][0];
                                            for(int i1 = 1; i1 <= numc; i1++)
                                            {
                                                index++;
                                                if(index > num_c + num_v) break;
                                                Best_Route_Education[index] = temp_route_split[i2][i1];
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    // update route
    for(int s = 0; s < num_c + num_v; s++)
    {
        sol[s] = Best_Route_Education[s];
    }
}

double Util::dist_comsum(double distance, double eng_consum)
{
    return distance * eng_consum;
}

bool Util::find_through_station(int * sol, bool *through_stat, int num_c, int num_v, double max_eng, double eng_consum, double ** Distances,int **Best_Stat, double **Best_Stat_Distances)
{    
    bool is_back_track[num_c + num_v];
    for(int i = 0; i < num_c + num_v; i++)
    {
        is_back_track[i] = false;
    }
    int j;
    int first_route_num = num_c;
    int num_node = num_c + num_v;
    double available_eng = max_eng;
    bool is_feasible = true;
    
    for(j = 0; j < num_node - 1;j++)
    {
        int node = sol[j];
        if(node >= first_route_num)
            node = 0;
        
        int next_node = sol[j + 1];
        if(next_node >= first_route_num)
            next_node = 0;
        
        if(node == 0 || node >= first_route_num)
        {
            available_eng = max_eng;
        }
        
        if(node != next_node && node >= 0 && next_node >= 0)
        {
            if(is_feasible)
            {
                double dist = Distances[node][next_node];
                if(available_eng > dist_comsum(dist, eng_consum))
                {
                    available_eng -= dist_comsum(dist, eng_consum);
                    through_stat[j] = false;
                    if(available_eng < 0) is_feasible = false;
                } else
                {
                    int best_stat = Best_Stat[node][next_node];
                    if(available_eng > dist_comsum(Distances[best_stat][node], eng_consum) && max_eng > dist_comsum(Distances[best_stat][next_node], eng_consum)){
                        available_eng = max_eng - dist_comsum(Distances[best_stat][next_node], eng_consum);
                        if(available_eng < 0)
                            is_feasible = false;
                        through_stat[j] = true;
                    } else {
                        // try to back track 1 time to find better sol
                        int prev_node = sol[j - 1];
                        if(prev_node >= num_c) prev_node = 0;
                        double prev_avail = available_eng + dist_comsum(Distances[prev_node][node], eng_consum);

                        through_stat[j - 1] = true;
                        int best_prev_stat = Best_Stat[prev_node][node];
                        if(prev_avail > dist_comsum(Distances[best_prev_stat][prev_node], eng_consum) && max_eng > dist_comsum(Distances[best_prev_stat][node], eng_consum) && !is_back_track[j])
                        {
                            available_eng = max_eng - dist_comsum(Distances[best_prev_stat][node], eng_consum);
                            is_back_track[j] = true;
                            j--;
                        }
                        else {
                            is_feasible = false;
                            through_stat[j] = false;
                        }
//                        is_feasible = false;
//                        through_stat[j] = false;
                    }
                }
            }else
            {
                through_stat[j] = false;
            }
        }
    }
    return is_feasible;
}

bool Util::Education(int *seq_node, bool *though_station, double **Distances, int num_c, int num_v, double max_eng, double eng_consum, int **best_stat, double **best_stat_distances)
{
    set_init_best(seq_node, num_c + num_v);
    two_exchange_education(seq_node, Distances, num_c, num_v);
    or_exhcange(seq_node, Distances, num_c, num_v);
    //cross_exchange(seq_node, Distances, num_c, num_v);
    //return find_through_station(seq_node, though_station, num_c, num_v, max_eng, eng_consum, Distances, best_stat, best_stat_distances);
    return false;
}

void Util::save_sol_to_pool(int *seq, int *pool_seq, int num_node)
{
    for(int i = 0; i < num_node; i++)
    {
        pool_seq[i] = seq[i];
    }
}

void Util::Interchange10APR(int *sol, double **Distances, int num_c, int num_v, int comb[][2], int loop)
{
    int temp_route_split[num_v][num_c + num_v];
    for(int i = 0 ; i < num_v; i++)
    {
        for(int j = 0; j < num_v + num_c; j++)
        {
            temp_route_split[i][j] = num_c + num_v + 1;
        }
    }
    
    int s;
    for(s = 0; s < num_c + num_v; s++)
    {
        Best_Route_Education[s] = sol[s];
    }
    
    double best_cost = compute_cost_sol(sol, Distances, num_c, num_v);
    double init_cost = best_cost;
    
    // split route
    int index_route = 0;
    int index_node = 0;
    for(int i = 1; i < num_c + num_v; i++)
    {
        int node = sol[i];
        if(node < num_c)
        {
            index_node++;
            temp_route_split[index_route][index_node] = node;
        } else
        {
            index_node++;
            temp_route_split[index_route][index_node] = node;
            temp_route_split[index_route][0] = index_node;
            if(i < num_c + num_v - 1)
            {
                index_route++;
                index_node = 0;
                if(index_route > num_v)
                    return;
            }
        }
    }
    int i, j, k, temp;
    int ar[1000];
    for(i = 0; i < loop; i++)
        ar[i] = i;
    
    for(i = 0; i < loop - 1; i++)
    {
        j = Rand(i, loop - 1);
        temp = ar[i];
        ar[i] = ar[j];
        ar[j] = temp;
    }
    int best_index_1 = -1;
    int best_index_2 = -1;
    int best_route_1 = -1;
    int best_route_2 = -1;

    for(i = 0; i < loop; i++)
    {
        j = comb[ar[i]][0];
        k = comb[ar[i]][1];
        if(temp_route_split[j][0] > 1 && temp_route_split[k][0] > 1)
        {
            int num_node = temp_route_split[j][0];
            // search for the best position
            for(int s = 1; s < num_node; s++)
            {
                int node_route_1 = temp_route_split[j][s];
                if(temp_route_split[j][s] < num_c)
                {
                    int num_node_2 = temp_route_split[k][0];
                    for(int g = 1; g < num_node_2 - 1; g++)
                    {
                        int index_after = temp_route_split[k][g];
                        int next_node = temp_route_split[k][g + 1];
                        if(next_node > num_c) next_node = 0;
                        double new_cost = init_cost + (Distances[index_after][node_route_1] + Distances[node_route_1][next_node] - Distances[index_after][next_node]);
                        if(new_cost < best_cost)
                        {
                            best_index_1 = s;
                            best_index_2 = g;
                            best_route_1 = j;
                            best_route_2 = k;
                        }
                    }
                }
            }
        }
    }
    
    if(best_index_1 != -1 && best_index_2 != -1)
    {
        int val1 = temp_route_split[best_route_1][best_index_1];
        int val2 = temp_route_split[best_route_2][best_index_2];
        // change route
        int curr_idx = 0;
        sol[curr_idx] = 0;
        int num_node_1 = temp_route_split[best_route_1][0];
        int num_node_2 = temp_route_split[best_route_2][0];
        for(int i1 = 1; i1 <= num_node_1; i1++)
        {
            if(temp_route_split[best_route_1][i1] != val1)
            {
                curr_idx++;
                sol[curr_idx] = temp_route_split[best_route_1][i1];
            }
        }
        
        for(int i1 = 1; i1 <= best_index_2; i1++)
        {
            curr_idx++;
            sol[curr_idx] = temp_route_split[best_route_2][i1];
        }
        curr_idx++;
        sol[curr_idx] = val1;
        for(int i1 = best_route_2 + 1; i1 <= num_node_2; i1++)
        {
            curr_idx++;
            sol[curr_idx] = temp_route_split[best_route_2][i1];
        }
        
        for(int i2 = 0; i2 <= index_route; i2++)
        {
            if(i2 != best_route_1 && i2 != best_route_2 )
            {
                int n_node = temp_route_split[i2][0];
                for(int i1 = 1; i1 <= n_node; i1++)
                {
                    curr_idx++;
                    sol[curr_idx] = temp_route_split[i2][i1];
                }
            }
        }
        
    }
    
}

void Util::Interchange20APR(int *sol, double **Distances, int num_c, int num_v, int (*comb)[2], int loop)
{
    
    // split route
    int temp_route_split[num_v][num_c + num_v];
    for(int i = 0; i < num_v; i++)
    {
        for(int j = 0; j < num_v + num_c; j++)
        {
            temp_route_split[i][j] = num_c + num_v + 1;
        }
    }
    int s;
    for(s = 0; s < num_c + num_v; s++)
    {
        Best_Route_Education[s] = sol[s];
    }
    
    double best_cost = compute_cost_sol(sol, Distances, num_c, num_v);
    double init_cost = best_cost;
    
    // split route
    int index_route = 0;
    int index_node = 0;
    for(int i = 1; i < num_c + num_v; i++)
    {
        int node = sol[i];
        if(node < num_c)
        {
            index_node++;
            temp_route_split[index_route][index_node] = node;
        } else
        {
            index_node++;
            temp_route_split[index_route][index_node] = node;
            temp_route_split[index_route][0] = index_node;
            if(i < num_c + num_v - 1)
            {
                index_route++;
                index_node = 0;
                if(index_route > num_v)
                    return;
            }
        }
    }
    
    // finish split route
    int i, j, k, temp, cont;
    
    int ar[1000];
    for(i = 0; i < loop; i++)
        ar[i] = i;
    for(i = 0;i < loop - 1; i++)
    {
        j = Rand(i, loop - 1);
        temp = ar[i];
        ar[i] = ar[j];
        ar[j] = temp;
    }
    
    int best_index_1 = -1;
    int best_index_2 = -1;
    int best_route_1 = -1;
    int best_route_2 = -1;
    
    double best_gain = 0.0;
    for(i = 0; i < loop; i++)
    {
        j = comb[ar[i]][0];
        k = comb[ar[i]][1];
        if(temp_route_split[j][0] > 3 && temp_route_split[k][0] > 1 && j != k)
        {
            // search for position
            int numc_route_1 = temp_route_split[j][0];
            int numc_route_2 = temp_route_split[k][0];
            for(int s  = 1; s < numc_route_1 - 2; s++)
            {
                int nodeR1 = temp_route_split[j][s];
                int nextNodeR1 = temp_route_split[j][s + 1];
                int prev1 = 0;
                if(s > 1) prev1 = temp_route_split[j][s - 1];
                int nextR1 = temp_route_split[j][s + 2];
                if(nextR1 > num_c) nextR1 = 0;
                double reduce = Distances[prev1][nodeR1] + Distances[nextNodeR1][nextR1] - Distances[prev1][nextR1];
                for(int g = 1; g < numc_route_2 - 1; g++)
                {
                    int nodeR2 = temp_route_split[k][g];
                    int nextNodeR2 = temp_route_split[k][g + 1];
                    int prev2 = 0;
                    if(g > 1) prev2 = temp_route_split[k][g - 1];
                    int nextR2 = temp_route_split[k][g + 2];
                    double increase = Distances[nodeR2][nodeR1] + Distances[nextNodeR1][nextNodeR2] - Distances[nodeR2][nextNodeR2];
                    
                    double gain = reduce - increase;
                    if(gain > 0.0 && gain > best_gain)
                    {
                        best_gain = gain;
                        best_index_1 = s;
                        best_index_2 = g;
                        best_route_1 = j;
                        best_route_2 = k;
                    }
                }
            }
        }
    }
    
    // update best route to sol
    if(best_index_1 != -1 && best_index_2 != -1)
    {
//        printf("\nInterchange 20: %d - %d - %d - %d", best_route_1, best_route_2, best_index_1, best_index_2);
//            printf("\ndata route split\n");
//            for(int k = 0 ; k <= index_route; k++)
//            {
//                printf("\n");
//                int num_route_node = temp_route_split[k][0];
//                for(int m = 0; m <= num_route_node; m++)
//                {
//                    printf("%d -> ", temp_route_split[k][m]);
//                }
//            }
        int nodeR1 = temp_route_split[best_route_1][best_index_1];
        int nextR1 = temp_route_split[best_route_1][best_index_1 + 1];
        int nodeR2 = temp_route_split[best_route_2][best_index_2];
        int nextR2 = temp_route_split[best_route_2][best_index_2 + 1];
        // change route
        int curr_idx = 0;
        sol[curr_idx] = 0;
        int num_node_1 = temp_route_split[best_route_1][0];
        int num_node_2 = temp_route_split[best_route_2][0];
        for(int i1 = 1; i1 <= num_node_1; i1++)
        {
            if(temp_route_split[best_route_1][i1] == nodeR1 || temp_route_split[best_route_1][i1] == nextR1) continue;
            else{
                curr_idx++;
                sol[curr_idx] = temp_route_split[best_route_1][i1];
            }
        }
        
        for(int i1 = 1; i1 <= best_index_2; i1++)
        {
            curr_idx++;
            sol[curr_idx] = temp_route_split[best_route_2][i1];
        }
        curr_idx++;
        sol[curr_idx] = nodeR1;
        curr_idx++;
        sol[curr_idx] = nextR1;
        for(int i1 = best_index_2 + 1; i1 <= num_node_2; i1++)
        {
            curr_idx++;
            sol[curr_idx] = temp_route_split[best_route_2][i1];
        }
        // copy all other routes
        for(int i2 = 0; i2 <= index_route; i2++)
        {
            if(i2 != best_route_1 && i2 != best_route_2)
            {
                int n_node = temp_route_split[i2][0];
                for(int i1 = 1; i1 <= n_node; i1++)
                {
                    curr_idx++;
                    sol[curr_idx] = temp_route_split[i2][i1];
                }
            }
        }
    }
}

void Util::Interchange11APR(int *sol, double **Distances, int num_c, int num_v, int (*comb)[2], int loop)
{
    int temp_route_split[num_v][num_c + num_v];
    for(int i = 0; i < num_v; i++)
    {
        for(int j = 0; j < num_v + num_c; j++)
        {
            temp_route_split[i][j] = num_c + num_v + 1;
        }
    }
    
    int s;
    for(s = 0; s < num_c + num_v; s++)
    {
        Best_Route_Education[s] = sol[s];
    }
    
    double best_cost = compute_cost_sol(sol, Distances, num_c, num_v);
    double init_cost = best_cost;
    
    //split route
    int index_route = 0;
    int index_node = 0;
    for(int i = 1; i < num_c + num_v; i++)
    {
        int node = sol[i];
        if(node < num_c)
        {
            index_node++;
            temp_route_split[index_route][index_node] = node;
        } else
        {
            index_node++;
            temp_route_split[index_route][index_node] = node;
            temp_route_split[index_route][0] = index_node;
            if(i < num_c + num_v - 1)
            {
                index_route++;
                index_node = 0;
                if(index_route > num_v)
                    return;
            }
        }
    }
    
    int i, j, k, temp;
    int ar[1000];
    for(i = 0; i < loop; i++)
        ar[i] = i;
    for(i = 0; i < loop - 1; i++)
    {
        j = Rand(i, loop - 1);
        temp = ar[i];
        ar[i] = ar[j];
        ar[j] = temp;
    }
    
    int best_index_1 = -1;
    int best_index_2 = -1;
    int best_route_1 = -1;
    int best_route_2 = -1;
            printf("\ndata route split\n");
            for(int k = 0 ; k <= index_route; k++)
            {
                printf("\n");
                int num_route_node = temp_route_split[k][0];
                for(int m = 0; m <= num_route_node; m++)
                {
                    printf("%d -> ", temp_route_split[k][m]);
                }
            }
    double best_gain = 0.0;
    for(i = 0; i < loop; i++)
    {
        j = comb[ar[i]][0];
        k = comb[ar[i]][1];
        if(temp_route_split[j][0] > 1 && temp_route_split[k][0] > 1)
        {
            int numc_route_1 = temp_route_split[j][0];
            int numc_route_2 = temp_route_split[k][0];
            for(int s = 1; s < numc_route_1 - 1; s++)
            {
                int nodeR1 = temp_route_split[j][s];
                int nextNodeR1 = temp_route_split[j][s + 1];
                int prev1 = 0;
                if(s > 1) prev1 = temp_route_split[j][s - 1];
                int nextR1 = temp_route_split[j][s + 2];
                if(nextR1 > num_c) nextR1 = 0;
                
                for(int g = 1; g < numc_route_2 - 1; g++)
                {
                    int nodeR2 = temp_route_split[k][g];
                    int nextNodeR2 = temp_route_split[k][g + 1];
                    int prev2 = 0;
                    if(g > 1) prev2 = temp_route_split[k][g - 1];
                    int nextR2 = temp_route_split[k][g + 2];
                    double reduce1 = ((Distances[prev1][nodeR1] + Distances[nodeR1][nextNodeR1]) - (Distances[prev1][nodeR2] + Distances[nodeR2][nextNodeR1]));
                    double reduce2 = ((Distances[prev2][nodeR2] + Distances[nodeR2][nextNodeR2]) - (Distances[prev2][nodeR1] + Distances[nodeR1][nextNodeR2]));
                    
                    double gain = reduce1 + reduce2;
                    if(gain > 0 && gain > best_gain)
                    {
                        best_gain = gain;
                        best_index_1 = s;
                        best_index_2 = g;
                        best_route_1 = j;
                        best_route_2 = k;
                    }
                    
                }
            }
        }
    }
    
    // update route to sol
    if(best_index_1 != -1 && best_index_2 != -1)
    {
//        printf("\nInterchange 20: %d - %d - %d - %d", best_route_1, best_route_2, best_index_1, best_index_2);
//        printf("\ndata route split\n");
//        for(int k = 0 ; k <= index_route; k++)
//        {
//            printf("\n");
//            int num_route_node = temp_route_split[k][0];
//            for(int m = 0; m <= num_route_node; m++)
//            {
//                printf("%d -> ", temp_route_split[k][m]);
//            }
//        }
        int curr_idx = 0;
        sol[curr_idx] = 0;
        int num_node_1 = temp_route_split[best_route_1][0];
        int num_node_2 = temp_route_split[best_route_2][0];
        for(int i1 = 1; i1 <= num_node_1; i1++)
        {
            if(i1 == best_index_1)
            {
                curr_idx++;
                sol[curr_idx] = temp_route_split[best_route_2][best_index_2];
            } else {
                curr_idx++;
                sol[curr_idx] = temp_route_split[best_route_1][i1];
            }
        }
        
        for(int i1 = 1; i1 <= num_node_2; i1++)
        {
            if(i1 == best_index_2)
            {
                curr_idx++;
                sol[curr_idx] = temp_route_split[best_route_1][best_index_1];
            } else {
                curr_idx++;
                sol[curr_idx] = temp_route_split[best_route_2][i1];
            }
        }
        
        for(int i2 = 0 ; i2 <= index_route; i2++)
        {
            if(i2 != best_route_1 && i2 != best_route_2)
            {
                int n_node = temp_route_split[i2][0];
                for(int i1 = 1; i1 <= n_node; i1++)
                {
                    curr_idx++;
                    sol[curr_idx] = temp_route_split[i2][i1];
                }
            }
        }
    }
}

void Util::Interchange21APR(int *sol, double **Distances, int num_c, int num_v, int (*comb)[2], int loop)
{
    int temp_route_split[num_v][num_c + num_v];
    for(int i = 0; i < num_v; i++)
    {
        for(int j = 0; j < num_v + num_c; j++)
        {
            temp_route_split[i][j] = num_c + num_v + 1;
        }
    }
    int s;
    for(s = 0; s < num_c + num_v; s++)
    {
        Best_Route_Education[s] = sol[s];
    }
    double best_cost = compute_cost_sol(sol, Distances, num_c, num_v);
    double init_cost = best_cost;
    //split route
    int index_route = 0;
    int index_node = 0;
    for(int i = 1; i < num_c + num_v; i++)
    {
        int node = sol[i];
        if(node < num_c)
        {
            index_node++;
            temp_route_split[index_route][index_node] = node;
        } else
        {
            index_node++;
            temp_route_split[index_route][index_node] = node;
            temp_route_split[index_route][0] = index_node;
            if(i < num_c + num_v - 1)
            {
                index_route++;
                index_node = 0;
                if(index_route > num_v)
                    return;
            }
        }
    }
    
    int i, j, k, temp;
    int ar[1000];
    for(i = 0; i < loop; i++)
        ar[i] = i;
    for(i = 0; i < loop - 1; i++)
    {
        j = Rand(i, loop - 1);
        temp = ar[i];
        ar[i] = ar[j];
        ar[j] = temp;
    }
    
    int best_index_1 = -1;
    int best_index_2 = -1;
    int best_route_1 = -1;
    int best_route_2 = -1;
    
    double best_gain = 0.0;
    printf("\ndata route split\n");
    for(int k = 0 ; k <= index_route; k++)
    {
        printf("\n");
        int num_route_node = temp_route_split[k][0];
        for(int m = 0; m <= num_route_node; m++)
        {
            printf("%d -> ", temp_route_split[k][m]);
        }
    }
    for(i = 0; i < loop; i++)
    {
        j = comb[ar[i]][0];
        k = comb[ar[i]][1];
        if(temp_route_split[j][0] > 3 && temp_route_split[k][0] > 1)
        {
            int numc_route_1 = temp_route_split[j][0];
            int numc_route_2 = temp_route_split[k][0];
            for(int s = 1; s < numc_route_1 - 2; s++)
            {
                int nodeR1 = temp_route_split[j][s];
                int nextNodeR1 = temp_route_split[j][s + 1];
                int prev1 = 0;
                if(s > 1) prev1 = temp_route_split[j][s - 1];
                int nextR1 = temp_route_split[j][s + 2];
                if(nextR1 > num_c) nextR1 = 0;
                
                for(int g = 1; g < numc_route_2 - 1; g++)
                {
                    int nodeR2 = temp_route_split[k][g];
                    int nextNodeR2 = temp_route_split[k][g + 1];
                    int prev2 = 0;
                    if(g > 1) prev2 = temp_route_split[k][g - 1];
                    int nextR2 = temp_route_split[k][g + 2];
                    
                    double reduce1 = ((Distances[prev1][nodeR1] + Distances[nodeR1][nextNodeR1] + Distances[nextNodeR1][nextR1]) - (Distances[prev1][nodeR2] + Distances[nodeR2][nextR1]));
                    double reduce2 = ((Distances[prev2][nodeR2] + Distances[nodeR2][nextNodeR2]) - (Distances[prev2][nodeR1] + Distances[nodeR1][nextNodeR1] + Distances[nextNodeR1][nextNodeR2]));
                    
                    double gain = reduce1 + reduce2;
                    if(gain > 0 && gain > best_gain)
                    {
                        best_gain = gain;
                        best_index_1 = s;
                        best_index_2 = g;
                        best_route_1 = j;
                        best_route_2 = k;
                    }
                }
            }
        }
    }
    
    if(best_index_1 != -1 && best_index_2 != -1)
    {
        int curr_idx = 0;
        sol[curr_idx] = 0;
        int num_node_1 = temp_route_split[best_route_1][0];
        int num_node_2 = temp_route_split[best_route_2][0];
        
        int nodeR1 = temp_route_split[best_route_1][best_index_1];
        int nextR1 = temp_route_split[best_route_1][best_index_1 + 1];
        int nodeR2 = temp_route_split[best_route_2][best_index_2];
        int nextR2 = temp_route_split[best_route_2][best_index_2 + 1];
        for(int i1 = 1; i1 <= num_node_1; i1++)
        {
            if(temp_route_split[best_route_1][i1] == nodeR1)
            {
                curr_idx++;
                sol[curr_idx] = nodeR2;
            }
            else if(temp_route_split[best_route_1][i1] == nextR1) continue;
            else {
                curr_idx++;
                sol[curr_idx] = temp_route_split[best_route_1][i1];
            }
        }
        
        for(int i1 = 1; i1 < best_index_2; i1++)
        {
            curr_idx++;
            sol[curr_idx] = temp_route_split[best_route_2][i1];
        }
        curr_idx++;
        sol[curr_idx] = nodeR1;
        curr_idx++;
        sol[curr_idx] = nextR1;
        for(int i1 = best_index_2 + 1; i1 <= num_node_2; i1++)
        {
            curr_idx++;
            sol[curr_idx] = temp_route_split[best_route_2][i1];
        }
        
        for(int i2 = 0; i2 <= index_route; i2++)
        {
            if(i2 != best_route_1 && i2 != best_route_2)
            {
                int n_node = temp_route_split[i2][0];
                for(int i1 = 1; i1 <= n_node; i1++)
                {
                    curr_idx++;
                    sol[curr_idx] = temp_route_split[i2][i1];
                }
            }
        }
    }
    
            printf("\nInterchange 20: %d - %d - %d - %d", best_route_1, best_route_2, best_index_1, best_index_2);
            printf("\ndata route split\n");
            for(int k = 0 ; k <= index_route; k++)
            {
                printf("\n");
                int num_route_node = temp_route_split[k][0];
                for(int m = 0; m <= num_route_node; m++)
                {
                    printf("%d -> ", temp_route_split[k][m]);
                }
            }
}

void Util::Interchange22APR(int *sol, double **Distances, int num_c, int num_v, int (*comb)[2], int loop)
{
    int temp_route_split[num_v][num_c + num_v];
    for(int i = 0; i < num_v; i++)
    {
        for(int j = 0; j < num_v + num_c; j++)
        {
            temp_route_split[i][j] = num_c + num_v + 1;
        }
    }
    int s;
    for(s = 0; s < num_c + num_v; s++)
    {
        Best_Route_Education[s] = sol[s];
    }
    double best_cost = compute_cost_sol(sol, Distances, num_c, num_v);
    double init_cost = best_cost;
    
    //split route
    int index_route = 0;
    int index_node = 0;
    for(int i = 1; i < num_c + num_v; i++)
    {
        int node = sol[i];
        if(node < num_c)
        {
            index_node++;
            temp_route_split[index_route][index_node] = node;
        } else
        {
            index_node++;
            temp_route_split[index_route][index_node] = node;
            temp_route_split[index_route][0] = index_node;
            if(i < num_c + num_v - 1)
            {
                index_route++;
                index_node = 0;
                if(index_route > num_v)
                    return;
            }
        }
    }
    
    int i, j, k, temp;
    int ar[1000];
    for(i = 0; i < loop; i++)
        ar[i] = i;
    for(i = 0; i < loop - 1; i++)
    {
        j = Rand(i, loop - 1);
        temp = ar[i];
        ar[i] = ar[j];
        ar[j] = temp;
    }
    
    int best_index_1 = -1;
    int best_index_2 = -1;
    int best_route_1 = -1;
    int best_route_2 = -1;
    double best_gain = 0.0;
    
    for(i = 0; i < loop; i++)
    {
        for(i = 0; i < loop; i++)
        {
            j = comb[ar[i]][0];
            k = comb[ar[i]][1];
            if(temp_route_split[j][0] > 3 && temp_route_split[k][0] > 3)
            {
                int numc_route_1 = temp_route_split[j][0];
                int numc_route_2 = temp_route_split[k][0];
                for(int s = 1; s < numc_route_1 - 2; s++)
                {
                    int nodeR1 = temp_route_split[j][s];
                    int nextNodeR1 = temp_route_split[j][s + 1];
                    int prev1 = 0;
                    if(s > 1) prev1 = temp_route_split[j][s - 1];
                    int nextR1 = temp_route_split[j][s + 2];
                    if(nextR1 > num_c) nextR1 = 0;
                    for(int g = 1; g < numc_route_2 - 2; g++)
                    {
                        int nodeR2 = temp_route_split[k][g];
                        int nextNodeR2 = temp_route_split[k][g + 1];
                        int prev2 = 0;
                        if(g > 1) prev2 = temp_route_split[k][g - 1];
                        int nextR2 = temp_route_split[k][g + 2];
                        
                        double reduce1 = ((Distances[prev1][nodeR1] + Distances[nextNodeR1][nextR1]) - (Distances[prev1][nodeR2] + Distances[nextNodeR2][nextR1]));
                        double reduce2 = ((Distances[prev2][nodeR2] + Distances[nextNodeR2][nextR2]) - (Distances[prev2][nodeR1] + Distances[nextNodeR1][nextR2]));
                        
                        double gain = reduce1 + reduce2;
                        if(gain > 0 && gain > best_gain)
                        {
                            best_gain = gain;
                            best_index_1 = s;
                            best_index_2 = g;
                            best_route_1 = j;
                            best_route_2 = k;
                        }
                    }
                }
            }
        }
    }
    
    if(best_index_1 != -1 && best_index_2 != -1)
    {
        int curr_idx = 0;
        sol[curr_idx] = 0;
        int num_node_1 = temp_route_split[best_route_1][0];
        int num_node_2 = temp_route_split[best_route_2][0];
        
        int nodeR1 = temp_route_split[best_route_1][best_index_1];
        int nextR1 = temp_route_split[best_route_1][best_index_1 + 1];
        int nodeR2 = temp_route_split[best_route_2][best_index_2];
        int nextR2 = temp_route_split[best_route_2][best_index_2 + 1];
        
        for(int i1 = 1; i1 <= num_node_1; i1++)
        {
            if(temp_route_split[best_route_1][i1] == nodeR1)
            {
                curr_idx++;
                sol[curr_idx] = nodeR2;
                continue;
            }
            if(temp_route_split[best_route_1][i1] == nextR1)
            {
                curr_idx++;
                sol[curr_idx] = nextR2;
                continue;
            }
            curr_idx++;
            sol[curr_idx] = temp_route_split[best_route_1][i1];
        }
        
        for(int i1 = 1; i1 <= num_node_2; i1++)
        {
            if(temp_route_split[best_route_2][i1] == nodeR2)
            {
                curr_idx++;
                sol[curr_idx] = nodeR1;
                continue;
            }
            if(temp_route_split[best_route_2][i1] == nextR2)
            {
                curr_idx++;
                sol[curr_idx] = nextR1;
                continue;
            }
            curr_idx++;
            sol[curr_idx] = temp_route_split[best_route_2][i1];
        }
        for(int i2 = 0; i2 <= index_route; i2++)
        {
            if(i2 != best_route_1 && i2 != best_route_2)
            {
                int n_node = temp_route_split[i2][0];
                for(int i1 = 1; i1 <= n_node; i1++)
                {
                    curr_idx++;
                    sol[curr_idx] = temp_route_split[i2][i1];
                }
            }
        }
    }
}

void Util::local_search(int *seq, double **Distances, int num_c, int num_v, int num)
{
    int ar[50];
    int i, j, k, temp, loop;
    for(i = 0; i < num; i++)
        ar[i] = i;
    for(i = 0; i < num - 1; i++)
    {
        j = (int)(rand() % 2);
        temp = ar[i];
        ar[i] = ar[j];
        ar[j] = temp;
    }
    
    int comb[1000][2];
    k = -1;
    for(i = 0; i < num_v - 1; i++)
        for(j = i + 1; j < num_v; j++)
        {
            comb[++k][0] = i;
            comb[k][1] = j;
            comb[++k][0] = j;
            comb[k][1] = i;
        }
    loop = num_v * (num_v - 1);
    printf("\nTest LOCAL SEARCH\n");
    for(i = 0; i < num; i++)
    {
        switch (ar[i]) {
            case 0:
                Interchange10APR(seq, Distances, num_c, num_v, comb, loop);
                //printf("\nin interchange_10\n");
                break;
            case 1:
                Interchange20APR(seq, Distances, num_c, num_v, comb, loop);
                //printf("\nin interchange_20\n");
                break;
            case 2:
                Interchange11APR(seq, Distances, num_c, num_v, comb, loop);
                //printf("\nin interchange_11\n");
                break;
            case 3:
                Interchange21APR(seq, Distances, num_c, num_v, comb, loop);
                //printf("\nin interchange_21\n");
                break;
            case 4:
                Interchange22APR(seq, Distances, num_c, num_v, comb, loop);
                //printf("\nin interchange_22\n");
                break;
            default:
                break;
        }
    }
    printf("\n END LOCAL SEARCH\n");
}

// CODE FOR TABU SEARCH

int find_index(int *seq, int num_node, int cus)
{
    for(int i = 0; i < num_node; i++)
    {
        if(seq[i] == cus) return i;
    }
    return  -1;
}

bool Util::validate_through_stat(int *seq, int begin, int num_c, int num_v, double max_eng, double eng_consum, double **Distances, int **Best_Stat, double **Best_Stat_Distances, int idx_j, int idx_i, int node_j, int node_i){
    int current_idx = begin;
    double available_eng = max_eng;
    bool is_begin = true;
    int temp_route[num_c + num_v];
    int index = begin - 1;
    for(int s = begin; s < num_c + num_v; s++)
    {
        if(s == idx_i)
        {
            continue;
        } else if(s == idx_j)
        {
            index++;
            temp_route[index] = seq[s];
            index++;
            temp_route[index] = node_i;
        } else
        {
            index++;
            temp_route[index] = seq[s];
        }
    }
    while ((temp_route[current_idx] != 0 && temp_route[current_idx] < num_c) || (is_begin && temp_route[current_idx] == 0)) {
        is_begin = false;
        int current_node = temp_route[current_idx];
        if(current_node >= num_c) current_node = 0;
        int next_node = temp_route[current_idx + 1];
        if(next_node >= num_c) next_node = 0;
        if(current_node != next_node)
        {
            if(current_node == 0 || current_node >= num_c) available_eng = max_eng;
        }
        double dist = Distances[current_node][next_node];
        if(available_eng > dist_comsum(dist, eng_consum))
        {
            available_eng -= dist_comsum(dist, eng_consum);
            if(available_eng < 0) return false;
        } else
        {
            int best_stat = Best_Stat[current_node][next_node];
            if(available_eng > dist_comsum(Distances[best_stat][current_node], eng_consum) && max_eng > dist_comsum(Distances[best_stat][next_node], eng_consum))
            {
                available_eng = max_eng - dist_comsum(Distances[best_stat][next_node], eng_consum);
                if(available_eng < 0) return false;
            } else {
                // try to back track 1 time to find better sol
                if(current_idx == begin) return false;
                else {
                    int prev_node = temp_route[current_idx - 1];
                    if(prev_node >= num_c) prev_node = 0;
                    double prev_avail_eng = available_eng + dist_comsum(Distances[prev_node][current_node], eng_consum);
                    
                    int best_prev_stat = Best_Stat[prev_node][current_node];
                    if(prev_avail_eng > dist_comsum(Distances[best_prev_stat][best_prev_stat], eng_consum) && max_eng > dist_comsum(Distances[best_prev_stat][current_node], eng_consum)){
                        available_eng = max_eng - dist_comsum(Distances[best_prev_stat][current_node], eng_consum);
                    } else return false;
                }
            }
        }
        current_idx++;
    }
    return  true;
}

Move Util::CallEvaluate(int *seq, double fitness_Zt, int num_c, int num_v, int **List_Nearest_Cus, double init_cost, double ** Distances, double *Demands, double max_cap, double ALPHA, double max_eng, double eng_consum, int **Best_Stat, double ** Best_Stat_Distances, int IT)
{
    FOUND_NEW_BEST = false;
    Move best_move;
    best_move.cus_1 = -1;
    best_move.cus_2 = -1;
    best_move.varFitness = 10000000.0;
    for(int i = 1; i < num_c; i++)
    {
        int idx_i = find_index(seq, num_c + num_v, i);
        for(int k = 0; k < 10; k++)
        {
            Move move;
            int j = List_Nearest_Cus[i][k];
            int idx_j = find_index(seq, num_c + num_v, j);
            if(idx_i == idx_j + 1) continue;
            
            int prev_node_i = seq[idx_i - 1];
            int next_node_i = seq[idx_i + 1];
            int prev_node_j = seq[idx_j - 1];
            int next_node_j = seq[idx_j + 1];
            
            if(prev_node_i >= num_c) prev_node_i = 0;
            if(next_node_i >= num_c) next_node_i = 0;
            if(prev_node_j >= num_c) prev_node_j = 0;
            if(next_node_j >= num_c) next_node_j = 0;
            
            move.varCost = Distances[prev_node_i][next_node_i] + Distances[j][i] + Distances[i][next_node_j] - (Distances[prev_node_i][i] + Distances[i][next_node_i] + Distances[j][next_node_j]);
            //printf("\nVar_cost: %lf", move.varCost);
            //move.varVioCap =
            double sum_cap_j = Demands[j];
            int before_j = idx_j - 1;
            int after_j = idx_j + 1;
            bool is_on_the_route = false;
            if(before_j >= 0)
            {
                while(seq[before_j] < num_c && seq[before_j] != 0 && before_j >= 0)
                {
                    if(idx_i == before_j) is_on_the_route = true;
                    int node = seq[before_j];
                    sum_cap_j += Demands[node];
                    before_j--;
                }
            } else before_j = 0;
            
            while (seq[after_j] < num_c && seq[after_j] != 0) {
                if(after_j == idx_i) is_on_the_route = true;
                sum_cap_j += Demands[seq[after_j]];
                after_j++;
            }
            
            if(is_on_the_route) {
                move.varVioCap = 0.0;
            } else
            {
                // find over cap of route i
                double sum_cap_i = Demands[i];
                int before_i = idx_i - 1;
                int after_i = idx_i + 1;
                if(before_i >= 0)
                {
                    while (seq[before_i] < num_c && seq[before_i] != 0 && before_i >= 0)
                    {
                        int node = seq[before_i];
                        sum_cap_i += Demands[node];
                        before_i--;
                    }
                } else before_i = 0;
                while (seq[after_i] < num_c && seq[after_i] != 0)
                {
                    sum_cap_i += Demands[seq[after_i]];
                    after_i++;
                }
                
                double over_cap_j = sum_cap_j - max_cap;
                if(over_cap_j < 0.0) over_cap_j = 0.0;
                
                double over_cap_i = sum_cap_i - max_cap;
                if(over_cap_i < 0.0) over_cap_i = 0.0;
                
                double over_cap_i_new = (sum_cap_i - Demands[i]) - max_cap;
                if(over_cap_i_new < 0.0) over_cap_i_new = 0.0;
                
                double over_cap_j_new = (sum_cap_j + Demands[i]) - max_cap;
                if(over_cap_j_new < 0.0) over_cap_j_new = 0.0;
                
                move.varVioCap = (over_cap_i_new + over_cap_j_new) - (over_cap_i + over_cap_j);
            }
            
            move.varFitness = move.varCost + move.varVioCap * ALPHA;
            if(before_j == -1) before_j = 0;
            // validate through stat
            if(move.varVioCap > 0 || !validate_through_stat(seq, before_j, num_c, num_v, max_eng, eng_consum, Distances, Best_Stat, Best_Stat_Distances, idx_j, idx_i, j, i))
            {
                move.feasible = false;
            }
            else {
                move.feasible = true;
            }
            
            move.tabu = true;
            bool FOUND = false;
            // check aspiration
            if(move.feasible && (fitness_Zt + move.varFitness < FIT_BEST)) {
                move.tabu = false;
                FIT_BEST = fitness_Zt + move.varFitness;
                FOUND = true;
            } else if (tabu[i][j] < IT)
            {
                move.tabu = false;
            }
            
            if(!move.tabu) {
                if(FOUND_NEW_BEST)
                {
                    if(FOUND){
                        if(move.varFitness < best_move.varFitness){
                            best_move.cus_1 = i;
                            best_move.cus_2 = j;
                            best_move.varFitness = move.varFitness;
                            best_move.varVioCap = move.varVioCap;
                            best_move.varCost = move.varCost;
                        }
                    }
                } else {
                    if(FOUND) {
                        FOUND_NEW_BEST = true;
                        best_move.cus_1 = i;
                        best_move.cus_2 = j;
                        best_move.varFitness = move.varFitness;
                        best_move.varVioCap = move.varVioCap;
                        best_move.varCost = move.varCost;
                    } else if(move.varFitness < best_move.varFitness) {
                        best_move.cus_1 = i;
                        best_move.cus_2 = j;
                        best_move.varFitness = move.varFitness;
                        best_move.varVioCap = move.varVioCap;
                        best_move.varCost = move.varCost;
                    }
                }
            }
        }
    }
    return best_move;
}

void update_seq(int *origin, int i, int j, int num_node){
    int index = -1;
    int temp[num_node];
    for(int k = 0; k < num_node; k++){
        if(origin[k] == i) {
            continue;
        } else if(origin[k] == j) {
            index++;
            temp[index] = j;
            index++;
            temp[index] = i;
        } else {
            index++;
            temp[index] = origin[k];
        }
    }
    
    for(int k = 0; k < num_node; k++)
    {
        origin[k] = temp[k];
    }
}

void Util::Tabu_search(int *seq, double **Distances, int num_c, int num_v, double *Demands, double init_fitness, int **List_Nearest_Cus, double init_cost, double max_cap, double ALPHA, double max_eng, double eng_consum, int **Best_Stat, double ** Best_Stat_Distances)
{
    int tabu_num = (int)((num_c + num_v)*0.1);
    // reset tabu counter
    for(int i = 0; i < num_c; i++)
    {
        for(int j = 0; j < num_c; j++)
        {
            tabu[i][j] = -1;
        }
    }
    // init best tabu
    for(int i = 0; i < num_c + num_v; i++)
    {
        Best_Tabu_Search[i] = seq[i];
    }
    double finess_Zt = init_fitness;
    FIT_BEST = init_fitness;
    for(int it = 0; it < 150; it++)
    {
        Move best_move = CallEvaluate(seq, finess_Zt, num_c, num_v, List_Nearest_Cus, init_cost, Distances, Demands, max_cap, ALPHA, max_eng, eng_consum, Best_Stat, Best_Stat_Distances, it);
        if(best_move.cus_1 != -1)
        {
            int i = best_move.cus_1;
            int j = best_move.cus_2;
            tabu[i][j] = it + tabu_num;
            tabu[j][i]= it + tabu_num;
            
            
            
            finess_Zt = finess_Zt + best_move.varFitness;
            // create new route
            update_seq(seq, i, j, num_c + num_v);
            if(FOUND_NEW_BEST) {
                for(int s = 0; s < num_c + num_v; s++)
                {
                    Best_Tabu_Search[s] = seq[s];
                }
            }
        }
    }
    // update best search
    for(int i = 0; i < num_c + num_v; i++)
    {
        seq[i] = Best_Tabu_Search[i];
    }
}



// LOCAL SEARCH FIX

bool Util::validate_partial_route(int *route, int idx_j, int node_i, double **Distances, int **Best_Stat, double **Best_Stat_Distances, double max_eng, double eng_consum, int num_c)
{
    int num_node = route[0];
    int temp[200];
    temp[0] = 0;
    int curr_idx = 0;
    for(int i = 1; i <= num_node; i++)
    {
        if(i == idx_j)
        {
            curr_idx++;
            temp[curr_idx] = route[i];
            curr_idx++;
            temp[curr_idx] = node_i;
        } else
        {
            curr_idx++;
            temp[curr_idx] = route[i];
        }
    }
    
    curr_idx = 0;
    double avail_eng = max_eng;
    while(curr_idx < num_node + 1)
    {
        int current_node = temp[curr_idx];
        int next_node = temp[curr_idx + 1];
        if(next_node > num_c) next_node = 0;
        
        if(current_node != next_node)
        {
            double dist = Distances[current_node][next_node];
            if(avail_eng > dist_comsum(dist, eng_consum))
                avail_eng -= dist_comsum(dist, eng_consum);
            else
            {
                int best_stat = Best_Stat[current_node][next_node];
                if(avail_eng > dist_comsum(Distances[best_stat][current_node], eng_consum) && max_eng > dist_comsum(Distances[best_stat][next_node], eng_consum))
                {
                    avail_eng = max_eng - dist_comsum(Distances[best_stat][next_node], eng_consum);
                    if(avail_eng < 0.0) return false;
                } else
                {
                    // try to back track 1 time to find better sol
                    if(curr_idx == 0) return false;
                    else
                    {
                        int prev_node = temp[curr_idx - 1];
                        double prev_avail_eng = avail_eng + dist_comsum(Distances[prev_node][current_node], eng_consum);
                        
                        int best_prev_stat = Best_Stat[prev_node][current_node];
                        if(prev_avail_eng > dist_comsum(Distances[best_prev_stat][current_node], eng_consum) && max_eng > dist_comsum(Distances[best_prev_stat][current_node], eng_consum))
                        {
                            avail_eng = max_eng - dist_comsum(Distances[best_prev_stat][current_node], eng_consum);
                        } else return false;
                    }
                }
            }
        }
        curr_idx++;
    }
    return true;
}

void Util::local_search_FIX(Solution sol, double **Distances, int num_c, int num_v, double *Demands, double max_cap, double max_eng, double eng_consum, int **Best_Stat, double **Best_Stat_Distances)
{
    int num = 5;
    int ar[50];
    int i, j, k, temp, loop;
    
    for(i = 0; i < num; i++)
        ar[i] = i;
    
    for(i = 0; i < num - 1; i++)
    {
        j = (int)(rand() % 2);
        temp = ar[i];
        ar[i] = ar[j];
        ar[j] = temp;
    }
    
    int comb[1000][2];
    k = -1;
    for(i = 0; i < num_v - 1; i++)
        for(j = i + 1; j < num_v; j++)
        {
            comb[++k][0] = i;
            comb[k][1] = j;
            comb[++k][0] = j;
            comb[k][1] = i;
        }
    
    loop = num_v * (num_v - 1);
    
    for(i = 0; i < num; i++)
    {
        switch (ar[i]) {
            case 0:
                printf("\nin optimise10\n");
                Interchange10FIX(sol, Distances, num_c, num_v, Demands, max_cap, max_eng, eng_consum, Best_Stat, Best_Stat_Distances, loop, comb);
                break;
            case 1:
                printf("\nin optimise20\n");
                Interchange20FIX(sol, Distances, num_c, num_v, Demands, max_cap, max_eng, eng_consum, Best_Stat, Best_Stat_Distances, loop, comb);
                break;
            case 2:
                printf("\nin optimise11\n");
                Interchange11FIX(sol, Distances, num_c, num_v, Demands, max_cap, max_eng, eng_consum, Best_Stat, Best_Stat_Distances, loop, comb);
                break;
            case 3:
                printf("\nin optimise21\n");
                Interchange21FIX(sol, Distances, num_c, num_v, Demands, max_cap, max_eng, eng_consum, Best_Stat, Best_Stat_Distances, loop, comb);
                break;
            case 4:
                printf("\nin optimise21\n");
                Interchange22FIX(sol, Distances, num_c, num_v, Demands, max_cap, max_eng, eng_consum, Best_Stat, Best_Stat_Distances, loop, comb);
                break;
            default:
                break;
        }
    }
}

void Util::Interchange10FIX(Solution sol, double **Distances, int num_c, int num_v, double *Demands, double max_cap, double max_eng, double eng_consum, int **Best_Stat, double **Best_Stat_Distances, int loop, int comb[][2])
{
    int temp_route_split[num_v][num_c];
    double route_capacity[num_v];
    for(int i = 0; i < num_v; i++)
    {
        for(int j = 0; j < num_c + num_v; j++)
        {
            temp_route_split[i][j] = num_c + num_v + 1;
        }
    }
    
    for(int i = 0; i < num_v; i++)
        route_capacity[i] = 0.0;
    
    double best_cost = sol.cost;
    
    // SPLIT ROUTE
    int index_route = 0;
    int index_node = 0;
    for(int i = 1; i < num_c + num_v; i++)
    {
        int node = sol.seq_node[i];
        if(node < num_c)
        {
            index_node++;
            temp_route_split[index_route][index_node] = node;
            route_capacity[index_route] += Demands[node];
        } else
        {
            index_node++;
            temp_route_split[index_route][index_node] = node;
            
            temp_route_split[index_route][0] = index_node;
            
            if(i < num_c + num_v - 1)
            {
                index_route++;
                index_node = 0;
                if(index_route > num_v) return;
            }
        }
    }
    
    int i, j, temp;
    int route_1, route_2;
    int ar[1000];
    for(i = 0; i < loop; i++)
        ar[i] = i;
    
    for(i = 0; i < loop - 1; i++)
    {
        j = Rand(i, loop - 1);
        temp = ar[i];
        ar[i] = ar[j];
        ar[j] = temp;
    }
    
    int Best_Index_1 = -1;
    int Best_Index_2 = -1;
    int Best_R_Idx_1 = -1;
    int Best_R_Idx_2 = -1;
    double best_delta_cost = 10000.0;
    
    for(i = 0; i < loop; i++)
    {
        route_1 = comb[ar[i]][0];
        route_2 = comb[ar[i]][1];
        int num_node_1 = temp_route_split[route_1][0];
        int num_node_2 = temp_route_split[route_2][0];
        
        if(num_node_1 > 2 && num_node_2 > 2)
        {
            // SEARCH
            for(int it1 = 1; it1 < num_node_1; it1++)
            {
                int node_1 = temp_route_split[route_1][it1];
                int prev_node_1 = temp_route_split[route_1][it1 - 1];
                if(it1 == 1) prev_node_1 = 0;
                int next_node_1 = temp_route_split[route_1][it1 + 1];
                if(next_node_1 >= num_c) next_node_1 = 0;
                for(int it2 = 1; it2 < num_node_2; it2++)
                {
                    int node_2 = temp_route_split[route_2][it2];
                    int next_node_2 = temp_route_split[route_2][it2 + 1];
                    if(next_node_2 >= num_c) next_node_2 = 0;
                    double delta_cost = (Distances[prev_node_1][next_node_1] + Distances[node_2][node_1] + Distances[node_1][next_node_2]) - (Distances[prev_node_1][node_1] + Distances[node_1][next_node_1] + Distances[node_2][next_node_2]);
                    
                    if(delta_cost < 0.0 && delta_cost < best_delta_cost && (route_capacity[route_2] + Demands[node_1] <= max_cap) && validate_partial_route(temp_route_split[route_2], it2, node_1, Distances, Best_Stat, Best_Stat_Distances, max_eng, eng_consum, num_c))
                    {
                        printf("\nFOUND BETTER: %lf\n", delta_cost);
                        best_delta_cost = delta_cost;
                        Best_Index_1 = it1;
                        Best_Index_2 = it2;
                        Best_R_Idx_1 = route_1;
                        Best_R_Idx_2 = route_2;
                    }
                }
            }
        }
    }
    // update route
    if(Best_Index_1 != -1 && Best_Index_2 != -1)
    {
        int val1 = temp_route_split[Best_R_Idx_1][Best_Index_1];
        
        int cur_idx = 0;
        sol.seq_node[cur_idx] = 0;
        int num_node_1 = temp_route_split[Best_R_Idx_1][0];
        int num_node_2 = temp_route_split[Best_R_Idx_2][0];
        for(int i1 = 1; i1 <= num_node_1; i1++)
        {
            if(temp_route_split[Best_R_Idx_1][i1] != val1)
            {
                cur_idx++;
                sol.seq_node[cur_idx] = temp_route_split[Best_R_Idx_1][i1];
            }
        }
        for(int i1 = 1; i1 <= Best_Index_2; i1++)
        {
            cur_idx++;
            sol.seq_node[cur_idx] = temp_route_split[Best_R_Idx_2][i1];
        }
        cur_idx++;
        sol.seq_node[cur_idx] = val1;
        for(int i1 = Best_Index_2 + 1; i1 <= num_node_2; i1++)
        {
            cur_idx++;
            sol.seq_node[cur_idx] = temp_route_split[Best_R_Idx_2][i1];
        }
        
        for(int i2 = 0; i2 <= index_route; i2++)
        {
            if(i2 != Best_R_Idx_1 && i2 != Best_R_Idx_2)
            {
                int n_node = temp_route_split[i2][0];
                for(int i1 = 1; i1 <= n_node; i1++)
                {
                    cur_idx++;
                    sol.seq_node[cur_idx] = temp_route_split[i2][i1];
                }
            }
        }
    }
}

bool Util::validate_partial_route20(int *route, int idx_j, int idx_i, int node_j, int node_i, int next_i, double **Distances, int **Best_Stat, double **Best_Stat_Distances, double max_eng, double eng_consum, int numc)
{
    int num_node = route[0];
    int temp[200];
    temp[0] = 0;
    int curr_idx = 0;
    for(int i = 1; i < num_node; i++)
    {
        if(i == idx_j)
        {
            curr_idx++;
            temp[curr_idx] = route[i];
            curr_idx++;
            temp[curr_idx] = node_i;
            curr_idx++;
            temp[curr_idx] = next_i;
        } else
        {
            curr_idx++;
            temp[curr_idx] = route[i];
        }
    }
    
    curr_idx = 0;
    double avail_eng = max_eng;
    while (curr_idx < num_node + 2)
    {
        int current_node = temp[curr_idx];
        int next_node = temp[curr_idx + 1];
        if(next_node > numc) next_node = 0;
        
        if(current_node != next_node)
        {
            //double dist = dist_comsum(Distances[current_node][next_node], eng_consum);
            double dist = Distances[current_node][next_node];
            if(avail_eng > dist_comsum(dist, eng_consum))
                avail_eng -= dist_comsum(dist, eng_consum);
            else
            {
                int best_stat = Best_Stat[current_node][next_node];
                if(avail_eng > dist_comsum(Distances[best_stat][current_node], eng_consum) && max_eng > dist_comsum(Distances[best_stat][next_node], eng_consum))
                {
                    avail_eng = max_eng - dist_comsum(Distances[best_stat][next_node], eng_consum);
                    if(avail_eng < 0.0) return false;
                } else
                {
                    // try to back track 1 time to find better sol
                    if(curr_idx == 0) return  false;
                    else
                    {
                        int prev_node = temp[curr_idx - 1];
                        double prev_avail_eng = avail_eng + dist_comsum(Distances[prev_node][current_node], eng_consum);
                        int best_prev_stat = Best_Stat[prev_node][current_node];
                        if(prev_avail_eng > dist_comsum(Distances[best_prev_stat][current_node], eng_consum) && max_eng > dist_comsum(Distances[best_prev_stat][current_node], eng_consum))
                        {
                            avail_eng = max_eng - dist_comsum(Distances[best_prev_stat][current_node], eng_consum);
                        } else return false;
                    }
                }
            }
        }
        curr_idx++;
    }
    return true;
}

void Util::Interchange20FIX(Solution sol, double **Distances, int num_c, int num_v, double *Demands, double max_cap, double max_eng, double eng_consum, int **Best_Stat, double **Best_Stat_Distances, int loop, int (*comb)[2])
{
    int temp_route_split[num_v][num_c];
    double route_capacity[num_v];
    for(int i = 0; i < num_v; i++)
    {
        for(int j = 0; j < num_c + num_v; j++)
        {
            temp_route_split[i][j] = num_c + num_v + 1;
        }
    }
    
    for(int i = 0; i < num_v; i++)
        route_capacity[i] = 0.0;
    
    // SPLIT ROUTE
    int index_route = 0;
    int index_node = 0;
    for(int i = 1; i < num_c + num_v; i++)
    {
        int node = sol.seq_node[i];
        if(node < num_c)
        {
            index_node++;
            temp_route_split[index_route][index_node] = node;
            route_capacity[index_route] += Demands[node];
        } else
        {
            index_node++;
            temp_route_split[index_route][index_node] = node;
            
            temp_route_split[index_route][0] = index_node;
            
            if(i < num_c + num_v - 1)
            {
                index_route++;
                index_node = 0;
                if(index_route > num_v) return;
            }
        }
    }
    
    int i, j, temp;
    int route_1, route_2;
    int ar[1000];
    for(i = 0; i < loop; i++)
        ar[i] = i;
    
    for(i = 0; i < loop - 1; i++)
    {
        j = Rand(i, loop - 1);
        temp = ar[i];
        ar[i] = ar[j];
        ar[j] = temp;
    }
    
    int Best_Index_1 = -1;
    int Best_Index_2 = -1;
    int Best_R_Idx_1 = -1;
    int Best_R_Idx_2 = -1;
    double best_delta_cost = 10000.0;
    
    for(i = 0; i < loop; i++)
    {
        route_1 = comb[ar[i]][0];
        route_2 = comb[ar[i]][1];
        int num_node_1 = temp_route_split[route_1][0];
        int num_node_2 = temp_route_split[route_2][0];
        
        if(num_node_1 > 3 && num_node_2 > 3)
        {
            for(int it1 = 1; it1 < num_node_1 - 1; it1++)
            {
                int node_1 = temp_route_split[route_1][it1];
                int prev_node_1 = temp_route_split[route_1][it1 - 1];
                if(it1 == 1) prev_node_1 = 0;
                int next_node_1 = temp_route_split[route_1][it1 + 1];
                int more_node_1 = temp_route_split[route_1][it1 + 1];
                if(more_node_1 >= num_c) more_node_1 = 0;
                for(int it2 = 1; it2 < num_node_2; it2++)
                {
                    int node_2 = temp_route_split[route_2][it2];
                    int next_node_2 = temp_route_split[route_2][it2 + 1];
                    if(next_node_2 >= num_c) next_node_2 = 0;
                    
                    double delta_cost = (Distances[prev_node_1][more_node_1] + Distances[node_2][node_1] + Distances[node_1][next_node_1] + Distances[next_node_1][next_node_2]) - (Distances[prev_node_1][node_1] + Distances[node_1][next_node_1] + Distances[next_node_1][more_node_1] + Distances[node_2][next_node_2]);
                    
                    if(delta_cost < 0.0 && delta_cost < best_delta_cost && (route_capacity[route_2] + Demands[node_1] + Demands[next_node_1] < max_cap) && validate_partial_route20(temp_route_split[route_2], it2, it1, node_2, node_1, next_node_1, Distances, Best_Stat, Best_Stat_Distances, max_eng, eng_consum, num_c))
                    {
                        
                        best_delta_cost = delta_cost;
                        Best_Index_1 = it1;
                        Best_Index_2 = it2;
                        Best_R_Idx_1 = route_1;
                        Best_R_Idx_2 = route_2;
                    }
                }
            }
        }
    }
    
    // update route
    if(Best_Index_1 != -1 && Best_Index_2 != -1)
    {
        int val1 = temp_route_split[Best_R_Idx_1][Best_Index_1];
        int next_val = temp_route_split[Best_R_Idx_1][Best_Index_1 + 1];
        int cur_idx = 0;
        sol.seq_node[cur_idx] = 0;
        int num_node_1 = temp_route_split[Best_R_Idx_1][0];
        int num_node_2 = temp_route_split[Best_R_Idx_2][0];
        
        for(int i1 = 1; i1 <= num_node_1; i1++)
        {
            if(temp_route_split[Best_R_Idx_1][i1] != val1 && temp_route_split[Best_R_Idx_1][i1] != next_val)
            {
                cur_idx++;
                sol.seq_node[cur_idx] = temp_route_split[Best_R_Idx_1][i1];
            }
        }
        
        for(int i1 = 1; i1 <= Best_Index_2; i1++)
        {
            cur_idx++;
            sol.seq_node[cur_idx] = temp_route_split[Best_R_Idx_2][i1];
        }
        cur_idx++;
        sol.seq_node[cur_idx] = val1;
        cur_idx++;
        sol.seq_node[cur_idx] = next_val;
        
        for(int i1 = Best_Index_2 + 1; i1 < num_node_2; i1++)
        {
            cur_idx++;
            sol.seq_node[cur_idx] = temp_route_split[Best_R_Idx_2][i1];
        }
        
        for(int i2 = 0; i2 <= index_route; i2++)
        {
            if(i2 != Best_R_Idx_1 && i2 != Best_R_Idx_2)
            {
                int n_node = temp_route_split[i2][0];
                for(int i1 = 1; i1 <= n_node; i1++)
                {
                    cur_idx++;
                    sol.seq_node[cur_idx] = temp_route_split[i2][i1];
                }
            }
        }
    }
}

bool Util::validate_partial_route11(int *route, int idx_j, int idx_i, int node_j, int node_i, int next_i, double **Distances, int **Best_Stat, double **Best_Stat_Distances, double max_eng, double eng_consum, int numc)
{
    int num_node = route[0];
    int temp[200];
    temp[0] = 0;
    int curr_idx = 0;
    for(int i = 1; i < num_node; i++)
    {
        if(i == idx_j)
        {
            curr_idx++;
            temp[curr_idx] = node_i;
        } else
        {
            curr_idx++;
            temp[curr_idx] = route[i];
        }
    }
    
    curr_idx = 0;
    double avail_eng = max_eng;
    
    while (curr_idx < num_node) {
        int current_node = temp[curr_idx];
        int next_node = temp[curr_idx + 1];
        if(next_node > numc) next_node = 0;
        
        if(current_node != next_node)
        {
            double dist = Distances[current_node][next_node];
            if(avail_eng > dist_comsum(dist, eng_consum))
                avail_eng -= dist_comsum(dist, eng_consum);
            else
            {
                int best_stat = Best_Stat[current_node][next_node];
                if(avail_eng > dist_comsum(Distances[best_stat][current_node], eng_consum) && max_eng > dist_comsum(Distances[best_stat][next_node], eng_consum))
                {
                    avail_eng = max_eng - dist_comsum(Distances[best_stat][next_node], eng_consum);
                    if(avail_eng < 0.0) return false;
                } else
                {
                    // try to back track 1 time to find better sol
                    if(curr_idx == 0) return false;
                    else
                    {
                        int prev_node = temp[curr_idx - 1];
                        double prev_avail_eng = avail_eng + dist_comsum(Distances[prev_node][current_node], eng_consum);
                        
                        int best_prev_stat = Best_Stat[prev_node][current_node];
                        if(prev_avail_eng > dist_comsum(Distances[best_prev_stat][current_node], eng_consum) && max_eng > dist_comsum(Distances[best_prev_stat][current_node], eng_consum))
                        {
                            avail_eng = max_eng - dist_comsum(Distances[best_prev_stat][current_node], eng_consum);
                        } else return false;
                    }
                }
            }
        }
        curr_idx++;
    }
    return true;
}

void Util::Interchange11FIX(Solution sol, double **Distances, int num_c, int num_v, double *Demands, double max_cap, double max_eng, double eng_consum, int **Best_Stat, double **Best_Stat_Distances, int loop, int (*comb)[2])
{
    int temp_route_split[num_v][num_c];
    double route_capacity[num_v];
    for(int i = 0; i < num_v; i++)
    {
        for(int j = 0; j < num_c + num_v; j++)
        {
            temp_route_split[i][j] = num_c + num_v + 1;
        }
    }
    for(int i = 0; i < num_v; i++)
        route_capacity[i] = 0.0;
    
    double best_cost = sol.cost;
    
    // SPLIT ROUTE
    int index_route = 0;
    int index_node = 0;
    
    for(int i = 1; i < num_c + num_v; i++)
    {
        int node = sol.seq_node[i];
        if(node < num_c)
        {
            index_node++;
            temp_route_split[index_route][index_node] = node;
            route_capacity[index_route] += Demands[node];
        } else
        {
            index_node++;
            temp_route_split[index_route][index_node] = node;
            
            temp_route_split[index_route][0] = index_node;
            
            if(i < num_c + num_v - 1)
            {
                index_route++;
                index_node = 0;
                if(index_route > num_v) return;
            }
        }
    }
    
    int i, j, temp;
    int route_1, route_2;
    int ar[1000];
    for(i = 0; i < loop; i++)
        ar[i] = i;
    
    for(i = 0; i < loop - 1; i++)
    {
        j = Rand(i, loop - 1);
        temp = ar[i];
        ar[i] = ar[j];
        ar[j] = temp;
    }
    
    int Best_Index_1 = -1;
    int Best_Index_2 = -1;
    int Best_R_Idx_1 = -1;
    int Best_R_Idx_2 = -1;
    double best_delta_cost = 10000.0;
    
    for(i = 0; i < loop; i++)
    {
        route_1 = comb[ar[i]][0];
        route_2 = comb[ar[i]][1];
        int num_node_1 = temp_route_split[route_1][0];
        int num_node_2 = temp_route_split[route_2][0];
        
        if(num_node_1 > 2 && num_node_2 > 2)
        {
            // SEARCH
            for(int it1 = 1; it1 < num_node_1 ; it1++)
            {
                int node_1 = temp_route_split[route_1][it1];
                int prev_1 = temp_route_split[route_1][it1 - 1];
                if(it1 == 1) prev_1 = 0;
                int next_1 = temp_route_split[route_1][it1 + 1];
                if(next_1 >= num_c) next_1 = 0;
                for(int it2 = 1; it2 < num_node_2; it2++)
                {
                    int node_2 = temp_route_split[route_2][it2];
                    int next_2 = temp_route_split[route_2][it2 + 1];
                    if(next_2 >= num_c) next_2 = 0;
                    int prev_2 = temp_route_split[route_2][it2 - 1];
                    if(it2 == 1) prev_2 = 0;
                    
                    double delta_cost = (Distances[prev_2][node_1] + Distances[node_1][next_2] + Distances[prev_1][node_2] + Distances[node_2][next_1]) - (Distances[prev_1][node_1] + Distances[node_1][next_1] + Distances[prev_2][node_2] + Distances[node_2][next_2]);
                    if(delta_cost < 0.0 && delta_cost < best_delta_cost && (route_capacity[route_1] + Demands[node_2] < max_cap) && (route_capacity[route_2] + Demands[node_1]< max_cap) && validate_partial_route11(temp_route_split[route_2], it2, it1, node_2, node_1, next_1, Distances, Best_Stat, Best_Stat_Distances, max_eng, eng_consum, num_c))
                    {
                       printf("\nFOUND BETTER: %lf\n", delta_cost);
                        best_delta_cost = delta_cost;
                        Best_Index_1 = it1;
                        Best_Index_2 = it2;
                        Best_R_Idx_1 = route_1;
                        Best_R_Idx_2 = route_2;
                    }
                }
            }
        }
    }
}

void Util::Interchange21FIX(Solution sol, double **Distances, int num_c, int num_v, double *Demands, double max_cap, double max_eng, double eng_consum, int **Best_Stat, double **Best_Stat_Distances, int loop, int (*comb)[2])
{
    int temp_route_split[num_v][num_c];
    double route_capacity[200];
    for(int i = 0; i < num_v; i++)
    {
        for(int j = 0; j < num_c + num_v; j++)
        {
            temp_route_split[i][j] = num_c + num_v + 1;
        }
    }
    
    for(int i = 0; i < num_v; i++)
        route_capacity[i] = 0.0;
    
    double best_cost = sol.cost;
    
    // SPLIT ROUTE
    int index_route = 0;
    int index_node = 0;
    for(int i = 1; i < num_c + num_v; i++)
    {
        int node = sol.seq_node[i];
        if(node < num_c)
        {
            index_node++;
            temp_route_split[index_route][index_node] = node;
            route_capacity[index_route] += Demands[node];
        } else
        {
            index_node++;
            temp_route_split[index_route][index_node] = node;
            
            temp_route_split[index_route][0] = index_node;
            
            if(i < num_c + num_v - 1)
            {
                index_route++;
                index_node = 0;
                if(index_route > num_v) return;
            }
        }
    }
    
    int i, j, temp;
    int route_1, route_2;
    int ar[1000];
    for(i = 0; i < loop; i++)
        ar[i] = i;
    
    for(i = 0; i < loop - 1; i++)
    {
        j = Rand(i, loop - 1);
        temp = ar[i];
        ar[i] = ar[j];
        ar[j] = temp;
    }
    
    int Best_Index_1 = -1;
    int Best_Index_2 = -1;
    int Best_R_Idx_1 = -1;
    int Best_R_Idx_2 = -1;
    double best_delta_cost = 10000.0;
    
    for(i = 0; i < loop; i++)
    {
        route_1 = comb[ar[i]][0];
        route_2 = comb[ar[i]][1];
        int num_node_1 = temp_route_split[route_1][0];
        int num_node_2 = temp_route_split[route_2][0];
        
        if(num_node_1 > 3 && num_node_2 > 3)
        {
            // SEARCH
            for(int it1 = 1; it1 < num_node_1 - 2; it1++)
            {
                int node_1 = temp_route_split[route_1][it1];
                int prev_node_1 = temp_route_split[route_1][it1 - 1];
                if(it1 == 1) prev_node_1 = 0;
                int next_node_1 = temp_route_split[route_1][it1 + 1];
                int more_node_1 = temp_route_split[route_1][it1 + 2];
                if(more_node_1 >= num_c) more_node_1 = 0;
                
                for(int it2 = 1; it2 < num_node_2; it2++)
                {
                    int node_2 = temp_route_split[route_2][it2];
                    int next_node_2 = temp_route_split[route_2][it2 + 1];
                    if(next_node_2 >= num_c) next_node_2 = 0;
                    int prev_node_2 = temp_route_split[route_2][it2 - 1];
                    if(it2 == 1) prev_node_2 = 0;
                    
                    double detal_cost = (Distances[prev_node_2][node_1] + Distances[node_1][next_node_1] + Distances[next_node_1][next_node_2] + Distances[prev_node_1][node_2] + Distances[node_2][next_node_1]) - (Distances[prev_node_1][node_1] + Distances[node_1][next_node_1] + Distances[next_node_1][more_node_1] + Distances[prev_node_2][node_2] + Distances[node_1][next_node_2]);
                    
                    if(detal_cost < -5.0 && detal_cost < best_delta_cost && (route_capacity[route_1] - Demands[node_1] - Demands[next_node_1] + Demands[node_2] < max_cap) && (route_capacity[route_2] - Demands[node_2] + Demands[node_1] + Demands[next_node_1] < max_cap))
                    {
                        best_delta_cost = detal_cost;
                        Best_Index_1 = it1;
                        Best_Index_2 = it2;
                        Best_R_Idx_1 = route_1;
                        Best_R_Idx_2 = route_2;
                    }
                }
            }
        }
    }
    
    // update route
    if(Best_Index_1 != -1 && Best_Index_2 != 2)
    {
         int curr_idx = 0;
        sol.seq_node[0] = 0;
        int num_node_1 = temp_route_split[Best_R_Idx_1][0];
        int num_node_2 = temp_route_split[Best_Index_2][0];
        
        int nodeR1 = temp_route_split[Best_R_Idx_1][Best_Index_1];
        int next_node_R1 = temp_route_split[Best_R_Idx_1][Best_Index_1 + 1];
        int nodeR2 = temp_route_split[Best_R_Idx_2][Best_Index_2];
        int next_node_R2 = temp_route_split[Best_R_Idx_2][Best_Index_2 + 1];
        
        for(int i1 = 1; i1 <= num_node_1; i1++)
        {
            if(temp_route_split[Best_R_Idx_1][i1] == nodeR1)
            {
                curr_idx++;
                sol.seq_node[curr_idx] = nodeR2;
            } else if(temp_route_split[Best_R_Idx_1][i1] == next_node_R1) continue;
            else
            {
                curr_idx++;
                sol.seq_node[curr_idx] = temp_route_split[Best_R_Idx_1][i1];
            }
        }
        
        for(int i1 = 1; i1 < Best_Index_2; i1++)
        {
            curr_idx++;
            sol.seq_node[curr_idx] = temp_route_split[Best_R_Idx_2][i1];
        }
        curr_idx++;
        sol.seq_node[curr_idx] = nodeR1;
        curr_idx++;
        sol.seq_node[curr_idx] = next_node_R1;
        
        for(int i1 = Best_Index_2 + 1; i1 <- num_node_2; i1++)
        {
            curr_idx++;
            sol.seq_node[curr_idx] = temp_route_split[Best_R_Idx_2][i1];
        }
        
        for(int i2 = 0; i2 <= index_route; i2++)
        {
            if(i2 != Best_R_Idx_2 && i2 != Best_R_Idx_1)
            {
                int n_node = temp_route_split[i2][0];
                for(int i1 = 1; i1 <= n_node; i1++)
                {
                    curr_idx++;
                    sol.seq_node[curr_idx] = temp_route_split[i2][i1];
                }
            }
        }
    }
}

void Util::Interchange22FIX(Solution sol, double **Distances, int num_c, int num_v, double *Demands, double max_cap, double max_eng, double eng_consum, int **Best_Stat, double **Best_Stat_Distances, int loop, int (*comb)[2])
{
    int temp_route_split[num_v][num_c];
    double route_capacity[100];
    for(int i = 0; i < num_v; i++)
    {
        for(int j = 0; j < num_c + num_v; j++)
        {
            temp_route_split[i][j] = num_c + num_v + 1;
        }
    }
    
    for(int i = 0; i < num_v; i++)
        route_capacity[i] = 0.0;
    
    // SPLIT ROUTE
    int index_route = 0;
    int index_node = 0;
    for(int i = 1; i < num_c + num_v; i++)
    {
        int node = sol.seq_node[i];
        if(node < num_c)
        {
            index_node++;
            temp_route_split[index_route][index_node] = node;
            route_capacity[index_route] += Demands[node];
        } else
        {
            index_node++;
            temp_route_split[index_route][index_node] = node;
            
            temp_route_split[index_route][0] = index_node;
            
            if(i < num_c + num_v - 1)
            {
                index_route++;
                index_node = 0;
                if(index_route > num_v) return;
            }
        }
    }
    
    int i, j, temp;
    int route_1, route_2;
    int ar[1000];
    for(i = 0; i < loop; i++)
        ar[i] = i;
    
    for(i = 0; i < loop - 1; i++)
    {
        j = Rand(i, loop - 1);
        temp = ar[i];
        ar[i] = ar[j];
        ar[j] = temp;
    }
    
    int Best_Index_1 = -1;
    int Best_Index_2 = -1;
    int Best_R_Idx_1 = -1;
    int Best_R_Idx_2 = -1;
    double best_delta_cost = 10000.0;
    
    for(i = 0; i < loop; i++)
    {
        route_1 = comb[ar[i]][0];
        route_2 = comb[ar[i]][1];
        int num_node_1 = temp_route_split[route_1][0];
        int num_node_2 = temp_route_split[route_2][0];
        
        if(num_node_1 > 3 && num_node_2 > 3)
        {
            for(int it1 = 1; it1 < num_node_1 - 2; it1++)
            {
                int node_1 = temp_route_split[route_1][it1];
                int prev_1 = temp_route_split[route_1][it1 - 1];
                if(it1 == 1) prev_1 = 0;
                int next_1 = temp_route_split[route_1][it1 + 1];
                int more_1 = temp_route_split[route_1][it1 + 2];
                if(more_1 >= num_c) more_1 = 0;
                
                for(int it2 = 1; it2 < num_node_2 - 2; it2++)
                {
                    int node_2 = temp_route_split[route_2][it2];
                    int prev_2 = temp_route_split[route_2][it2 - 1];
                    if(it2 == 1) prev_2 = 0;
                    int next_2 = temp_route_split[route_2][it2 + 1];
                    int more_2 = temp_route_split[route_2][it2 + 2];
                    if(more_2 >= num_c) more_2 = 0;
                    
                    double delta_cost = (Distances[prev_1][node_2] + Distances[node_2][next_2] + Distances[next_2][more_2] + Distances[prev_2][node_1] + Distances[node_1][next_1] + Distances[next_1][next_2]) - (Distances[prev_1][node_1] + Distances[node_1][next_1] + Distances[next_1][more_1] + Distances[prev_2][node_2] + Distances[node_2][next_2] + Distances[next_2][more_2]);
                    
                    if(delta_cost < best_delta_cost && delta_cost < 0.0 && (route_capacity[route_1] - Demands[node_1] - Demands[next_1] + Demands[node_2] + Demands[next_2] < max_cap) && (route_capacity[route_2] - Demands[node_2] - Demands[next_2] + Demands[node_1] + Demands[next_1] < max_cap))
                    {
                        printf("\nFOUND BETTER: %lf\n", delta_cost);
                        best_delta_cost = delta_cost;
                        Best_Index_1 = it1;
                        Best_Index_2 = it2;
                        Best_R_Idx_1 = route_1;
                        Best_R_Idx_2 = route_2;
                    }
                }
            }
        }
    }
    
    // update route
    if(Best_Index_1 != -1 && Best_Index_2 != -1)
    {
        int cur_idx = 0;
        sol.seq_node[cur_idx] = 0;
        int num_node_1 = temp_route_split[Best_R_Idx_1][0];
        int num_node_2 = temp_route_split[Best_R_Idx_2][0];
        
        int nodeR1 = temp_route_split[Best_R_Idx_1][Best_Index_1];
        int nodeR2 = temp_route_split[Best_R_Idx_2][Best_Index_2];
        int nextR1 = temp_route_split[Best_R_Idx_1][Best_Index_1 + 1];
        int nextR2 = temp_route_split[Best_R_Idx_2][Best_Index_2 + 1];
        
        for(int i1 = 1; i1 <= num_node_1; i1++)
        {
            if(temp_route_split[Best_R_Idx_1][i1] == nodeR1)
            {
                cur_idx++;
                sol.seq_node[cur_idx] = nodeR2;
            } else if(temp_route_split[Best_R_Idx_1][i1] == nextR1) {
                cur_idx++;
                sol.seq_node[cur_idx] = nextR2;
            } else {
                cur_idx++;
                sol.seq_node[cur_idx] = temp_route_split[Best_R_Idx_1][i1];
            }
        }
        
        for(int i1 = 1;  i1 <= num_node_2; i1++)
        {
            if(temp_route_split[Best_R_Idx_2][i1] == nodeR2)
            {
                cur_idx++;
                sol.seq_node[cur_idx] = nodeR1;
            } else if(temp_route_split[Best_R_Idx_2][i1] == nextR2){
                cur_idx++;
                sol.seq_node[cur_idx] = nextR1;
            } else {
                cur_idx++;
                sol.seq_node[cur_idx] = temp_route_split[Best_R_Idx_2][i1];
            }
        }
        
        for(int i2 = 0; i2 <= index_route; i2++)
        {
            if(i2 != Best_R_Idx_1 && i2 != Best_R_Idx_1)
            {
                int n_node = temp_route_split[i2][0];
                for(int i1 = 1; i1 <= n_node; i1++)
                {
                    cur_idx++;
                    sol.seq_node[cur_idx] = temp_route_split[i2][i1];
                }
            }
        }
    }
}

void Util::Two_Exchange(int *sol, double **Distances, int num_c, int num_v, double *Demands, double init_finess, int **List_Nearest_Cus, double init_cost, double max_cap, double ALPHA, double max_eng, double eng_consum, int **Best_Stat, double **Best_Stat_Distances)
{
    Move best_move;
    best_move.varCost = 100000.0;
    best_move.cus_1 = -1;
    best_move.cus_2 = -1;
    int begin_route_num = num_c;
    int s, begin = 0, end = 0, i, j;
    int parent[num_c + num_v];
    for(s = 0; s < num_c + num_v; s++)
    {
        parent[s] = sol[s];
        Best_Route_Education[s] = sol[s];
    }
    while (parent[begin] < (num_c + num_v))
    {
        begin = end + 1;
        while (parent[begin] >= begin_route_num)
            begin++;
        end = begin + 1;
        while (parent[end] < begin_route_num)
            end++;
        
        if((parent[begin] > num_c + num_v) || (parent[end] > num_c + num_v)) break;
        if(begin >= num_c + num_v || end >= num_c + num_v) break;
        
        for(i = begin; i < end - 2; i++)
            for(j = i + 2; j < end - 1; j++)
            {
                int i_1 = i + 1;
                int j_1 = j + 1;
                int node_i = parent[i];
                int node_i_1 = parent[i_1];
                int node_j = parent[j];
                int node_j_1 = parent[j_1];
                double var_cost = (Distances[node_i][node_j] + Distances[node_i_1][node_j_1] - (Distances[node_i][node_i_1] + Distances[node_j][node_j_1]));
                if(var_cost < 0.0 && var_cost < best_move.varCost && validate_two_exchange(sol, begin, end, i, j, max_eng, eng_consum, Distances, Best_Stat, Best_Stat_Distances, num_c, num_v))
                {
                    best_move.varCost = var_cost;
                    best_move.cus_1 = i;
                    best_move.cus_2 = j;
                }
            }
        if(end == num_c + num_v) break;
    }
    
    // update route
    if(best_move.cus_1 != -1 &&best_move.cus_2 != -1)
    {
        int cus_1 = best_move.cus_1;
        int cus_2 = best_move.cus_2;
        for(int it = 0; it <= cus_1; it++)
            sol[it] = parent[it];
        
        for(int it = cus_2; it < num_c + num_v; it++)
            sol[it] = parent[it];
        
        int index = cus_1;
        for(int it = cus_2 - 1; it > cus_1; it--)
        {
            index++;
            sol[index] = parent[it];
        }
    }
}

bool Util::validate_two_exchange(int *seq, int begin, int end, int i, int j, double max_eng, double eng_consum, double **Distances, int **Best_Stat, double **Best_Stat_Distances, int num_c, int num_v)
{
    int new_seq[num_c + num_v];
    bool though_stat[num_c + num_v];
    for(int it = 0; it <= i; it++)
        new_seq[it] = seq[it];
    for(int it = j; it < num_c + num_v; it++)
        new_seq[it] = seq[it];
    
    int index = i;
    for(int it = j - 1; it > i; it--)
    {
        index++;
        new_seq[index] = seq[it];
    }
    
    return find_through_station(new_seq, though_stat, num_c, num_v, max_eng, eng_consum, Distances, Best_Stat, Best_Stat_Distances);
}











void Util::Or_Exchange(int *sol, double **Distances, int num_c, int num_v, double *Demands, double init_finess, int **List_Nearest_Cus, double init_cost, double max_cap, double ALPHA, double max_eng, double eng_consum, int **Best_Stat, double **Best_Stat_Distances)
{
    Move best_move;
    best_move.varCost = 100000.0;
    best_move.cus_1 = -1;
    best_move.cus_2 = -1;
    int begin_route_num = num_c;
    int s, begin = 0, end = 0, i, j;
    int parent[num_c + num_v];
    for(s = 0; s < num_c + num_v; s++)
    {
        parent[s] = sol[s];
        Best_Route_Education[s] = sol[s];
    }
    
    while (parent[begin] < (num_c + num_v)) {
        begin = end + 1;
        while (parent[begin] >= begin_route_num)
            begin++;
        end = begin + 1;
        while (parent[end] < begin_route_num)
            end++;
        
        if((parent[begin] > num_c + num_v) || (parent[end] > num_c + num_v)) break;
        if(begin >= num_c + num_v || end >= num_c + num_v) break;
        for(i = begin + 1; i < end - 5;i++)
        {
            for(j = i + 3; j < end - 1; j++)
            {
                int j_add_1 = j + 1;
                int i1_sub_1 = i - 1;
                int i1 = i;
                int i2 = i1 + 1;
                int i2_add_1 = i2 + 1;
                
                int node_j = parent[j];
                int node_j_add_1 = parent[j_add_1];
                int node_i1 = parent[i1];
                int node_i1_sub_1 = parent[i1_sub_1];
                int node_i2 = parent[i2];
                int node_i2_add_1 = parent[i2_add_1];
                
                double var_cost = (Distances[node_i1_sub_1][node_i2_add_1] + Distances[node_j_add_1][node_i2] + Distances[node_i2][node_i1] + Distances[node_j][node_i1]) - (Distances[node_i1_sub_1][node_i1] + Distances[node_i1][node_i2] + Distances[node_i2][node_i2_add_1] + Distances[node_j][node_j_add_1]);
                
                if(var_cost < 0.0 && var_cost < best_move.varCost && validate_or_exchange(sol, i, j, max_eng, eng_consum, Distances, Best_Stat, Best_Stat_Distances, num_c, num_v))
                {
                    best_move.varCost = var_cost;
                    best_move.cus_1 = i1;
                    best_move.cus_2 = j;
                }
            }
        }
    }
    
    // update best move
    if(best_move.cus_1 != -1 && best_move.cus_2 != -1)
    {
        int i1 = best_move.cus_1;
        int i1_sub_1 = i1 - 1;
        int i2 = i1 + 1;
        int i2_add_1 = i2 + 1;
        int j = best_move.cus_2;
        int j_add_1 = j + 1;
        
        int k = 0;
        for(k = 0; k <= i1_sub_1; k++)
        {
            sol[k] = parent[k];
        }
        
        sol[k] = parent[i2_add_1];
        for(int k1 = i2_add_1 + 1; k1 <= j; k1++)
        {
            k++;
            sol[k] = parent[k1];
        }
        for(int k1 = i1; k1 <= i2; k1++)
        {
            k++;
            sol[k] = parent[k1];
        }
        for(int k1 = j_add_1; k1 < num_c + num_v; k1++)
        {
            k++;
            sol[k] = parent[k1];
        }
    }
}

bool Util::validate_or_exchange(int *seq, int i, int j, double max_eng, double eng_consum, double **Distances, int **Best_Stat, double **Best_Stat_Distances, int num_c, int num_v)
{
    int i1_sub_1 = i - 1;
    int i1 = i;
    int i2 = i1 + 1;
    int i2_add_1 = i2 + 1;
    int j_add_1 = j + 1;
    int new_seq[num_c + num_v];
    bool though_stat[num_c + num_v];
    int k = 0;
    for(k = 0; k <= i1_sub_1; k++)
    {
        new_seq[k] = seq[k];
    }
    
    new_seq[k] = seq[i2_add_1];
    for(int k1 = i2_add_1 + 1; k1 <= j; k1++)
    {
        k++;
        new_seq[k] = seq[k1];
    }
    for(int k1 = i1; k1 <= i2; k1++)
    {
        k++;
        new_seq[k] = seq[k1];
    }
    for(int k1 = j_add_1; k1 < num_c + num_v; k1++)
    {
        k++;
        new_seq[k] = seq[k1];
    }
    
    return  find_through_station(new_seq, though_stat, num_c, num_v, max_eng, eng_consum, Distances, Best_Stat, Best_Stat_Distances);
}


void Util::Cross_Exchange(int *sol, double **Distances, int num_c, int num_v, double *Demands, double init_finess, int **List_Nearest_Cus, double init_cost, double max_cap, double ALPHA, double max_eng, double eng_consum, int **Best_Stat, double **Best_Stat_Distances)
{
    int parent[num_c + num_v];
    bool through_stat[num_c + num_v];
    double route_cap[num_v];
    for(int s = 0; s < num_v; s++)
        route_cap[s] = 0.0;
    double best_var_cost = 10000.0;
    int temp_route_split[num_v][num_c + num_v];
    for(int i = 0; i < num_v; i++)
        for(int j = 0; j < num_v + num_c; j++)
            temp_route_split[i][j] = num_c + num_v + 1;
    
    // split route
    int index_route = 0;
    int index_node = 0;
    for(int i = 1; i < num_c + num_v; i++)
    {
        int node = sol[i];
        if(node < num_c)
        {
            index_node++;
            temp_route_split[index_route][index_node] = node;
            route_cap[index_route] += Demands[node];
        }
        else
        {
            index_node++;
            temp_route_split[index_route][index_node] = node;
            temp_route_split[index_route][0] = index_node;
            if(i < num_c + num_v - 1)
            {
                index_route++;
                index_node = 0;
                if(index_route > num_v)
                    return;
            }
        }
    }
    
    // cross exchange
    for(int i = 0; i < index_route; i++)
        for(int j = i + 1; j < index_route; j++)
        {
            int numc_route_i = temp_route_split[i][0];
            int numc_route_j = temp_route_split[j][0];
            
            if(numc_route_i >= 5 && numc_route_j >= 5)
                for(int X1 = 1; X1 <= numc_route_i - 4; X1++)
                    for(int Y1 = X1 + 2; Y1 <= numc_route_i - 2; Y1++)
                    {
                        int X1_P = X1 + 1;
                        int Y1_P = Y1 + 1;
                        
                        for(int X2 = 1; X2 < numc_route_j - 4; X2++)
                            for(int Y2 = X2 + 2; Y2 < numc_route_j - 2; Y2++)
                            {
                                int X2_P = X2 + 1;
                                int Y2_P = Y2 + 1;
                                
                                int node_X1 = temp_route_split[i][X1];
                                int node_X1_P = temp_route_split[i][X1_P];
                                int node_Y1 = temp_route_split[i][Y1];
                                int node_Y1_P = temp_route_split[i][Y1_P];
                                
                                int node_X2 = temp_route_split[j][X2];
                                int node_X2_P = temp_route_split[j][X2_P];
                                int node_Y2 = temp_route_split[j][Y2];
                                int node_Y2_P = temp_route_split[j][Y2_P];
                                
                                double var_cost = ((Distances[node_X1][node_X2_P] + Distances[node_Y2][node_Y1_P] + Distances[node_X2][node_X1_P] + Distances[node_Y1][node_Y2_P]) - (Distances[node_X1][node_X1_P] + Distances[node_Y1][node_Y1_P] + Distances[node_X2][node_X2_P] + Distances[node_Y2][node_Y2_P]));
                                
                                double demand_i = 0.0;
                                for(int s = X1_P; s <= Y1;s++)
                                    demand_i += Demands[temp_route_split[i][s]];
                                
                                double demand_j = 0.0;
                                for(int s = X2_P; s <= Y2; s++)
                                    demand_j += Demands[temp_route_split[j][s]];
                                
                                if(var_cost < 0.0 && var_cost < best_var_cost && (route_cap[i] + demand_j - demand_i) < max_cap && (route_cap[j] + demand_i - demand_j) < max_cap)
                                {
                                    best_var_cost = var_cost;
                                    // build route
                                    int index = 0;
                                    parent[0] = 0;
                                    for(int i1 = 1; i1 <= X1; i1++)
                                    {
                                        index++;
                                        if(index > num_c + num_v) break;
                                        parent[index] = temp_route_split[i][i1];
                                    }
                                    for(int i1 = X2_P; i1 <= Y2; i1++)
                                    {
                                        index++;
                                        if(index > num_c + num_v) break;
                                        parent[index] = temp_route_split[j][i1];
                                    }
                                    for(int i1 = Y1_P; i1 <= numc_route_i; i1++)
                                    {
                                        index++;
                                        if(index > num_c + num_v) break;
                                        parent[index] = temp_route_split[i][i1];
                                    }
                                    for(int i1 = 1; i1 <= X2; i1++)
                                    {
                                        index++;
                                        if(index > num_c + num_v) break;
                                        parent[index] = temp_route_split[j][i1];
                                    }
                                    for(int i1 = X1_P; i1 <= Y1; i1++)
                                    {
                                        index++;
                                        if(index > num_c + num_v) break;
                                        parent[index] = temp_route_split[i][i1];
                                    }
                                    for(int i1 = Y2_P; i1 <= numc_route_j; i1++)
                                    {
                                        index++;
                                        if(index > num_c + num_v) break;
                                        parent[index] = temp_route_split[j][i1];
                                    }
                                    
                                    for(int i2 = 0; i2 <= index_route; i2++)
                                    {
                                        if(i2 != i && i2 != j)
                                        {
                                            int numc = temp_route_split[i2][0];
                                            for(int i1 = 1; i1 <= numc; i1++)
                                            {
                                                index++;
                                                if(index > num_c + num_v) break;
                                                parent[index] = temp_route_split[i2][i1];
                                            }
                                        }
                                    }
                                    
                                    if(find_through_station(parent, through_stat, num_c, num_v, max_eng, eng_consum, Distances, Best_Stat, Best_Stat_Distances))
                                    {
                                        printf("\nRoute Cap i: %lf", (route_cap[i] + (Demands[X2_P] + Demands[Y2]) - (Demands[X1_P] + Demands[Y1])));
                                        printf("\nRoute Cap j: %lf", (route_cap[j] + (Demands[X1_P] + Demands[Y1]) - (Demands[X2_P] + Demands[Y2])));
                                        printf("\nMax Cap: %lf", max_cap);
                                        for(int s = 0; s < num_c + num_v; s++)
                                            sol[s] = parent[s];
                                        return;
                                    }
                                }
                            }
                    }
        }
}


bool check_move_feasible(double *Demands, double max_capacity_vh, int num_c, int num_v, int *seq_node)
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
            sum_capacity_route += Demands[node];
        } else
        {
            double temp_over = max_capacity_vh - sum_capacity_route;
            sum_capacity_route = 0.0;
            if(temp_over < 0)
            {
                sum_capacity_over += temp_over;
            }
        }
    }
    if(sum_capacity_over < 0.0)
    {
        return false;
    }
    
    return true;
}


void Util::Tabu_search_cvrp(int *seq, double **Distances, int num_c, int num_v, double *Demands, double init_finess, int **List_Nearest_Cus, double init_cost, double max_cap, double ALPHA, double max_eng, double eng_consum, int **Best_Stat, double **Best_Stat_Distances)
{
    int tabu_num = (int)((num_c + num_v)*0.1);
    // reset tabu counter
    for(int i = 0; i < num_c; i++)
    {
        for(int j = 0; j < num_c; j++)
        {
            tabu[i][j] = -1;
        }
    }
    
    // init best tabu
    for(int i = 0; i < num_c + num_v; i++)
    {
        Best_Tabu_Search[i] = seq[i];
    }
    
    double finess_Zt = init_finess;
    FIT_BEST = init_finess;
    for(int it = 0; it < 100; it++)
    {
        Move best_move = CallEvaluate_cvrp(seq, finess_Zt, num_c, num_v, List_Nearest_Cus, init_cost, Distances, Demands, max_cap, ALPHA_TABU, max_eng, eng_consum, Best_Stat, Best_Stat_Distances, it);
        if(best_move.cus_1 != -1)
        {
            int i = best_move.cus_1;
            int j = best_move.cus_2;
            tabu[i][j] = it + tabu_num;
            tabu[j][i]= it + tabu_num;
            
            finess_Zt = finess_Zt + best_move.varFitness;
            // create new route
            update_seq(seq, i, j, num_c + num_v);
            if(FOUND_NEW_BEST) {
                if(check_move_feasible(Demands, max_cap, num_c, num_v, seq))
                {
                    for(int s = 0; s < num_c + num_v; s++)
                    {
                        Best_Tabu_Search[s] = seq[s];
                    }
                }
            }
        }
    }
    // update best search
    for(int i = 0; i < num_c + num_v; i++)
    {
        seq[i] = Best_Tabu_Search[i];
    }
}

Move Util::CallEvaluate_cvrp(int *seq, double fitness_Zt, int num_c, int num_v, int **List_Nearest_Cus, double init_cost, double **Distances, double *Demands, double max_cap, double ALPHA, double max_eng, double eng_consum, int **Best_Stat, double **Best_Stat_Distances, int IT)
{
    FOUND_NEW_BEST = false;
    Move best_move;
    best_move.cus_1 = -1;
    best_move.cus_2 = -1;
    best_move.varFitness = 10000000.0;
    for(int i = 1; i < num_c; i++)
    {
        int idx_i = find_index(seq, num_c + num_v, i);
        for(int k = 0; k < 10; k++)
        {
            Move move;
            int j = List_Nearest_Cus[i][k];
            int idx_j = find_index(seq, num_c + num_v, j);
            if(idx_i == idx_j + 1) continue;
            
            int prev_node_i = seq[idx_i - 1];
            int next_node_i = seq[idx_i + 1];
            int prev_node_j = seq[idx_j - 1];
            int next_node_j = seq[idx_j + 1];
            
            if(prev_node_i >= num_c) prev_node_i = 0;
            if(next_node_i >= num_c) next_node_i = 0;
            if(prev_node_j >= num_c) prev_node_j = 0;
            if(next_node_j >= num_c) next_node_j = 0;
            
            move.varCost = Distances[prev_node_i][next_node_i] + Distances[j][i] + Distances[i][next_node_j] - (Distances[prev_node_i][i] + Distances[i][next_node_i] + Distances[j][next_node_j]);
            //printf("\nVar_cost: %lf", move.varCost);
            //move.varVioCap =
            double sum_cap_j = Demands[j];
            int before_j = idx_j - 1;
            int after_j = idx_j + 1;
            bool is_on_the_route = false;
            if(i == 4 && j == 8 && IT == 4)
                printf("fix bugs");
            if(before_j >= 0)
            {
                while(seq[before_j] < num_c && seq[before_j] != 0 && before_j >= 0)
                {
                    if(idx_i == before_j) is_on_the_route = true;
                    int node = seq[before_j];
                    sum_cap_j += Demands[node];
                    before_j--;
                }
            } else before_j = 0;
            
            while (seq[after_j] < num_c && seq[after_j] != 0) {
                if(after_j == idx_i) is_on_the_route = true;
                sum_cap_j += Demands[seq[after_j]];
                after_j++;
            }
            
            if(is_on_the_route) {
                move.varVioCap = 0.0;
            } else
            {
                // find over cap of route i
                double sum_cap_i = Demands[i];
                int before_i = idx_i - 1;
                int after_i = idx_i + 1;
                if(before_i >= 0)
                {
                    while (seq[before_i] < num_c && seq[before_i] != 0 && before_i >= 0)
                    {
                        int node = seq[before_i];
                        sum_cap_i += Demands[node];
                        before_i--;
                    }
                } else before_i = 0;
                while (seq[after_i] < num_c && seq[after_i] != 0)
                {
                    sum_cap_i += Demands[seq[after_i]];
                    after_i++;
                }
                
                double over_cap_j = sum_cap_j - max_cap;
                if(over_cap_j < 0.0) over_cap_j = 0.0;
                
                double over_cap_i = sum_cap_i - max_cap;
                if(over_cap_i < 0.0) over_cap_i = 0.0;
                
                double over_cap_i_new = (sum_cap_i - Demands[i]) - max_cap;
                if(over_cap_i_new < 0.0) over_cap_i_new = 0.0;
                
                double over_cap_j_new = (sum_cap_j + Demands[i]) - max_cap;
                if(over_cap_j_new < 0.0) over_cap_j_new = 0.0;
                
                move.varVioCap = (over_cap_i_new + over_cap_j_new) - (over_cap_i + over_cap_j);
            }
            
            move.varFitness = move.varCost + move.varVioCap * ALPHA;
            if(before_j == -1) before_j = 0;
            // validate through stat
            if(move.varVioCap > 0)
                move.feasible = false;
            else
                move.feasible = true;
            
            move.tabu = true;
            bool FOUND = false;
            // check aspiration
            if(move.feasible && (fitness_Zt + move.varFitness < FIT_BEST)) {
                move.tabu = false;
                FIT_BEST = fitness_Zt + move.varFitness;
                FOUND = true;
            } else if (tabu[i][j] < IT)
                move.tabu = false;
            
            if(!move.tabu) {
                if(FOUND_NEW_BEST)
                {
                    if(FOUND){
                        if(move.varFitness < best_move.varFitness){
                            best_move.cus_1 = i;
                            best_move.cus_2 = j;
                            best_move.varFitness = move.varFitness;
                            best_move.varVioCap = move.varVioCap;
                            best_move.varCost = move.varCost;
                            best_move.feasible = move.feasible;
                        }
                    }
                } else {
                    if(FOUND) {
                        FOUND_NEW_BEST = true;
                        best_move.cus_1 = i;
                        best_move.cus_2 = j;
                        best_move.varFitness = move.varFitness;
                        best_move.varVioCap = move.varVioCap;
                        best_move.varCost = move.varCost;
                        best_move.feasible = move.feasible;
                    } else if(move.varFitness < best_move.varFitness) {
                        best_move.cus_1 = i;
                        best_move.cus_2 = j;
                        best_move.varFitness = move.varFitness;
                        best_move.varVioCap = move.varVioCap;
                        best_move.varCost = move.varCost;
                        best_move.feasible = move.feasible;
                    }
                }
            }
        }
    }
    if(!best_move.feasible || best_move.cus_1 == -1)
    {
        ALPHA_TABU = ALPHA_TABU * (1.5);
    } else if(best_move.feasible)
    {
        ALPHA_TABU = ALPHA_TABU / (1.5);
    }
    return best_move;
}
