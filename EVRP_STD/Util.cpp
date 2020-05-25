//
//  Util.cpp
//  EVRP_STD
//
//  Created by MAC on 5/15/20.
//  Copyright © 2020 MAC. All rights reserved.
//

#include "Util.hpp"
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
                double route_cost = init_cost + ((Distances[i][j] + Distances[i_1][j_1]) - (Distances[parent[i]][parent[i_1]] + Distances[parent[j]][parent[j_1]]));
                
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
                
                double new_cost = init_cost + ((Distances[i1_sub_1][i2_add_1] + Distances[i2][j] + Distances[i1][j_1]) - (Distances[i1_sub_1][i1] + Distances[i2][i2_add_1] + Distances[j_1][j]));
                
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
                                
                                double new_cost = init_cost + ((Distances[X1][X2_P] + Distances[Y2][Y1_P] + Distances[X2][X1_P] + Distances[Y1][Y2_P]) - (Distances[X1][X1_P] + Distances[Y1][Y1_P] + Distances[X2][X2_P] + Distances[Y2][Y2_P]));
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
//                printf("\nInterchange 20: %d - %d - %d - %d", best_route_1, best_route_2, best_index_1, best_index_2);
//                printf("\ndata route split\n");
//                for(int k = 0 ; k <= index_route; k++)
//                {
//                    printf("\n");
//                    int num_route_node = temp_route_split[k][0];
//                    for(int m = 0; m <= num_route_node; m++)
//                    {
//                        printf("%d -> ", temp_route_split[k][m]);
//                    }
//                }
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
    for(i = 0; i < num; i++)
    {
        switch (ar[i]) {
            case 0:
                Interchange10APR(seq, Distances, num_c, num_v, comb, loop);
                break;
            case 1:
                Interchange20APR(seq, Distances, num_c, num_v, comb, loop);
                break;
            case 2:
                Interchange11APR(seq, Distances, num_c, num_v, comb, loop);
                break;
            case 3:
                Interchange21APR(seq, Distances, num_c, num_v, comb, loop);
                break;
            case 4:
                Interchange22APR(seq, Distances, num_c, num_v, comb, loop);
                break;
            default:
                break;
        }
    }
}
