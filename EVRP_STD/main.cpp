//
//  main.cpp
//  EVRP_STD
//
//  Created by MAC on 5/14/20.
//  Copyright © 2020 MAC. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <chrono>
#include <iostream>
#include <cstdlib>
#include "Solution.hpp"
#include "Util.hpp"
using namespace std::chrono;
using namespace std;

// -- data problem --
int DIMENTION = 0;
int NUM_STATIONS = 0;
int NUM_CUSTOMERS = 0;
int NUM_VEHICLES = 0;

double MAX_CAPACITY_VH = 0.0;
double MAX_ENERGY_VH = 0.0;
double ENG_CONSUMTION = 0.0;
double OPTIMAL_VALUE = 0.0;

double **Distances;
double *Demands;
double **Coords;

int **Best_Station;
double **Best_Station_Distances;

Solution *Solutions;
int **Route_In_Pool;
bool **Though_Stattion_In_Pool;
bool *Is_Feasible_Pool;
double Alpha, Beta, H;

// --- meta data ---
//int NUM_SOL = 100;
#define NUM_SOL 100
double H_VAL = 0.0;
double ALPHA = 0.0;

int *list_farthest_seeds;
int *list_nearest_seeds;
int *list_rand_seeds;


// Util init
Util U;

// data for create pool
double virtual_fitness[NUM_SOL];
int virtual_index[NUM_SOL];
double Roulette_Wheel_Arr[NUM_SOL];
bool Is_In_Pool[NUM_SOL];
int Parent_Pool[NUM_SOL];

// -------------------------------------------------

double dist_comsum(double distance)
{
    return distance * ENG_CONSUMTION;
}

void create_space_mem()
{
    int i = 0;
    Coords = (double **)malloc(DIMENTION * sizeof(double *));
    Distances = (double **)malloc(DIMENTION * sizeof(double *));
    Demands = (double *)malloc(NUM_CUSTOMERS * sizeof(double));
    
    Best_Station = (int **)malloc(NUM_CUSTOMERS * sizeof(int *));
    Best_Station_Distances = (double **)malloc(NUM_CUSTOMERS * sizeof(double *));
    
    Solutions = (Solution *)malloc(NUM_SOL *sizeof(Solution));
    Route_In_Pool = (int **)malloc((int)(NUM_SOL * 0.7) * sizeof(int*));
    Though_Stattion_In_Pool = (bool **)malloc((int)(NUM_SOL * 0.7) * sizeof(bool *));
    Is_Feasible_Pool = (bool*)malloc((int)(NUM_SOL * 0.7) * sizeof(bool));
    
    for(i = 0; i < DIMENTION; i++){
        Coords[i] = (double *)malloc(2 * sizeof(double));
        Distances[i] = (double *)malloc(DIMENTION * sizeof(double));
        if(i < NUM_CUSTOMERS)
        {
            Best_Station[i] = (int *)malloc(NUM_CUSTOMERS * sizeof(int));
            Best_Station_Distances[i] = (double *)malloc(NUM_CUSTOMERS * sizeof(double));
            
        }
    }
    
    // init mem space  for each solutions
    for(i = 0 ; i < NUM_SOL; i++)
    {
        Solutions[i].init_mem_space(NUM_VEHICLES, NUM_CUSTOMERS);
        if(i < (int)(NUM_SOL * 0.7))
        {
            Route_In_Pool[i] = (int *)malloc((NUM_CUSTOMERS + NUM_VEHICLES) * sizeof(int));
            Though_Stattion_In_Pool[i] = (bool *)malloc((NUM_CUSTOMERS + NUM_VEHICLES) * sizeof(bool));
        }
    }
    
    // init space for seeds data
    list_farthest_seeds = (int *)malloc(NUM_VEHICLES * sizeof(int));
    list_nearest_seeds = (int *)malloc(NUM_VEHICLES * sizeof(int));
    list_rand_seeds = (int *)malloc(NUM_VEHICLES * sizeof(int));
    
}

void read_file(char *file_src)
{
    int i;
    FILE *infile;
    infile = fopen(file_src, "r");
    if(infile == NULL)
    {
        printf("READ FILE ERROR\n");
        exit(-1);
    } else{
        printf("READ FILE SUCCESS\n");
    }
    
    fscanf(infile, "%lf\n", &OPTIMAL_VALUE);
    fscanf(infile, "%d\n", &NUM_VEHICLES);
    fscanf(infile, "%d\n", &DIMENTION);
    fscanf(infile, "%d\n", &NUM_STATIONS);
    fscanf(infile, "%lf\n", &MAX_CAPACITY_VH);
    fscanf(infile, "%lf\n", &MAX_ENERGY_VH);
    fscanf(infile, "%lf\n", &ENG_CONSUMTION);
    
    NUM_CUSTOMERS = DIMENTION;
    DIMENTION = NUM_CUSTOMERS + NUM_STATIONS;
    
    printf("data: dimention: %d - num_customer: %d - capacity_vh: %lf - energy_vh: %lf - energy_consumtion: %lf - num_vehicles: %d - num_station: %d\n", DIMENTION, NUM_CUSTOMERS, MAX_CAPACITY_VH, MAX_ENERGY_VH, ENG_CONSUMTION, NUM_VEHICLES, NUM_STATIONS);
    
    create_space_mem();
    
    double index, x, y;
    double demand = 0.0;
    for(i = 0; i < NUM_CUSTOMERS; i++)
    {
        fscanf(infile, "%lf %lf %lf\n", &index, &x, &y);
        Coords[i][0] = x;
        Coords[i][1] = y;
    }
    
    for(i = NUM_CUSTOMERS; i < DIMENTION; i++)
    {
        fscanf(infile, "%lf %lf %lf\n", &index, &x, &y);
        Coords[i][0] = x;
        Coords[i][1] = y;
    }
    
    for(i = 0; i < NUM_CUSTOMERS; i++)
    {
        fscanf(infile, "%lf %lf\n", &index, &demand);
        Demands[i] = demand;
    }
    
    fclose(infile);
    printf("\n");
}

void find_best_stattion(int i, int j) {
    double min_distance = 10000000;
    int k;
    int best_stat = -1;
    for(k = NUM_CUSTOMERS; k < DIMENTION; k++)
    {
        if(Distances[i][k] + Distances[j][k] < min_distance)
        {
            min_distance = Distances[i][k] + Distances[j][k];
            best_stat = k;
        }
    }
    
    Best_Station_Distances[i][j] = min_distance;
    Best_Station_Distances[j][i] = min_distance;
    
    Best_Station[i][j] = best_stat;
    Best_Station[j][i] = best_stat;
}

bool find_exist_seed(int index, int num_curr_seeds, int *arr)
{
    for(int i = 0 ; i < num_curr_seeds; i++)
    {
        if(arr[i] == index) return true;
    }
    return false;
}

void reset_seeds(int * list) {
    for(int i = 0; i < NUM_VEHICLES; i++)
    {
        list[i] = -1;
    }
}

// finding seeds

void find_farthest_seed(){
    reset_seeds(list_farthest_seeds);
    int i = 0;
    for(i = 0 ; i < NUM_VEHICLES; i++){
        double max_distance = -10000000;
        int index_seed = -1;
        for(int j = 1; j < NUM_CUSTOMERS; j++)
        {
            if(!find_exist_seed(j, i, list_farthest_seeds))
            {
                double distance_depot = Distances[0][j];
                if(distance_depot > max_distance){
                    max_distance = distance_depot;
                    index_seed = j;
                }
            }
        }
        list_farthest_seeds[i] = index_seed;
    }
}

void find_nearest_seed(){
    reset_seeds(list_nearest_seeds);
    int i = 0;
    for(i = 0 ; i < NUM_VEHICLES; i++){
        double min_distance = 10000000;
        int index_seed = -1;
        for(int j = 1; j < NUM_CUSTOMERS; j++)
        {
            if(!find_exist_seed(j, i, list_nearest_seeds))
            {
                double distance_depot = Distances[0][j];
                if(distance_depot < min_distance)
                {
                    min_distance = distance_depot;
                    index_seed = j;
                }
            }
        }
        list_nearest_seeds[i] = index_seed;
    }
}

void prepare_data()
{
    int i, j;
    for(i = 0; i < DIMENTION; i++)
    {
        for(j = i + 1; j < DIMENTION; j++)
        {
            double distance = sqrt((Coords[i][0] - Coords[j][0]) * (Coords[i][0] - Coords[j][0]) + (Coords[i][1] - Coords[j][1]) * (Coords[i][1] - Coords[j][1]));
            Distances[i][j] = distance;
            Distances[j][i] = distance;
        }
    }
    
    for(i = 0 ; i < NUM_CUSTOMERS; i++)
    {
        for(j = i + 1; j < NUM_CUSTOMERS; j++)
        {
            find_best_stattion(i, j);
        }
    }
    
    // finding seeds
    find_farthest_seed();
    find_nearest_seed();
}

void random_select_seeds()
{
    // clone list index
    int list_cus_index[NUM_CUSTOMERS];
    for(int i = 0; i < NUM_CUSTOMERS; i++)
    {
        list_cus_index[i] = i;
    }
    
    for(int i = 0 ; i < NUM_VEHICLES; i++)
    {
        int rand_val = (int)(rand() % (NUM_CUSTOMERS - 1 - i)) + 1;
        int index = list_cus_index[rand_val];
        list_cus_index[rand_val] = list_cus_index[NUM_CUSTOMERS - 1 - i];
        list_cus_index[NUM_CUSTOMERS - 1 - i] = index;
    }
    
    for(int i = 0; i < NUM_VEHICLES; i++)
    {
        list_rand_seeds[i] = list_cus_index[NUM_CUSTOMERS - 1 - i];
    }
}

void init_population()
{
    // random init
    int i = 0;
    
    // init Potvin route with farthest seeds
    Solutions[i].Potvin_init_test(list_farthest_seeds, Distances, Demands, MAX_CAPACITY_VH);
    Solutions[i].is_feasible = U.find_through_station(Solutions[i].seq_node, Solutions[i].is_though_stat, NUM_CUSTOMERS, NUM_VEHICLES, MAX_ENERGY_VH, ENG_CONSUMTION, Distances, Best_Station, Best_Station_Distances);
    Solutions[i].compute_cost(Distances, Best_Station_Distances);
    Solutions[i].compute_over_cap(Demands, MAX_CAPACITY_VH);
    
    // init Potvin route with nearest seeds
    i++;
    Solutions[i].Potvin_init_test(list_nearest_seeds, Distances, Demands, MAX_CAPACITY_VH);
    Solutions[i].is_feasible = U.find_through_station(Solutions[i].seq_node, Solutions[i].is_though_stat, NUM_CUSTOMERS, NUM_VEHICLES, MAX_ENERGY_VH, ENG_CONSUMTION, Distances, Best_Station, Best_Station_Distances);
    Solutions[i].compute_cost(Distances, Best_Station_Distances);
    Solutions[i].compute_over_cap(Demands, MAX_CAPACITY_VH);
    i++;
    for(; i < NUM_SOL; i++)
    {
        random_select_seeds();
        Solutions[i].Potvin_init_test(list_rand_seeds, Distances, Demands, MAX_CAPACITY_VH);
        Solutions[i].is_feasible = U.find_through_station(Solutions[i].seq_node, Solutions[i].is_though_stat, NUM_CUSTOMERS, NUM_VEHICLES, MAX_ENERGY_VH, ENG_CONSUMTION, Distances, Best_Station, Best_Station_Distances);
        Solutions[i].compute_cost(Distances, Best_Station_Distances);
        Solutions[i].compute_over_cap(Demands, MAX_CAPACITY_VH);
    }
}

// prepare meta data
void compute_meta_data()
{
    int i;
    double H = U.compute_H_value(Solutions, NUM_SOL);
    double AVG_q = U.compute_AVG_q(Solutions, NUM_SOL);
    double alpha = 0.0;
    if(AVG_q == 0.0) alpha = 0;
    else alpha = H / AVG_q;
    
    H_VAL = H;
    ALPHA = alpha;
    for(i = 0; i < NUM_SOL; i++)
    {
        Solutions[i].fitness = Solutions[i].cost + alpha * Solutions[i].over_capacity;
    }
}

// === CREATE POOL FOR MATING ===

int compare(const void* a, const void* b)
{
    const double* x = (double*) a;
    const double* y = (double*) b;
    
    if(*x > *y) return 1;
    else if (*x < *y) return -1;
    
    return 0;
}

int compare_sols_fitness(const void *ptr1, const void *ptr2)
{
    Solution m1 = *(Solution*)ptr1;
    Solution m2 = *(Solution*)ptr2;
    
    if(m1.fitness > m2.fitness) return 1;
    else if(m1.fitness < m2.fitness) return -1;
    
    return  0;
}

void sort_sols()
{
    qsort(Solutions, NUM_SOL, sizeof(Solution), compare_sols_fitness);
}

void build_Roulette_wheel_arr()
{
    int i;
    sort_sols();
    Roulette_Wheel_Arr[0] = Solutions[NUM_SOL - 1].fitness;
    for(i = 1; i < 100; i++)
    {
        Roulette_Wheel_Arr[i] = Roulette_Wheel_Arr[i - 1] + Solutions[NUM_SOL - 1 - i].fitness;
    }
}

// call after Roulette_Wheel_Arr builded
void reset_in_pool()
{
    int i;
    for(i = 0; i < NUM_SOL; i++)
    {
        Is_In_Pool[i] = false;
    }
}

void select_parent_to_pool_distinct()
{
    reset_in_pool();
    int i, j, m;
    int total = (int)Roulette_Wheel_Arr[NUM_SOL - 1];
    
    for(i = 0; i < (int)(NUM_SOL * (0.7)); i++)
    {
        int rand_wheel = (int)(rand() % total);
        for(j = 0; j < NUM_SOL; j++)
        {
            if(Roulette_Wheel_Arr[j] > rand_wheel)
            {
                if(Is_In_Pool[j])
                {
                    for(m = 0; m < NUM_SOL; m++)
                    {
                        if(!Is_In_Pool[m])
                        {
                            Parent_Pool[i] = m;
                            Is_In_Pool[m] = true;
                            break;
                        }
                    }
                    
                } else
                {
                    Parent_Pool[i] = j;
                    Is_In_Pool[j] = true;
                }
                break;
            }
        }
    }
}


// Begin Cross over
void cross_over(){
    // create child 1 and 2 space
    int i;
    int num_node = NUM_CUSTOMERS + NUM_VEHICLES;
    int child_1[num_node];
    int child_2[num_node];
    int child_3[num_node];
    int child_4[num_node];
    child_1[0] = 0;
    child_2[0] = 0;
    child_3[0] = 0;
    child_4[0] = 0;
    
    // child solutions
    Solution Child1;
    Solution Child2;
    Solution Child3;
    Solution Child4;
    
    Child1.init_mem_space(NUM_VEHICLES, NUM_CUSTOMERS);
    Child2.init_mem_space(NUM_VEHICLES, NUM_CUSTOMERS);
    Child3.init_mem_space(NUM_VEHICLES, NUM_CUSTOMERS);
    Child4.init_mem_space(NUM_VEHICLES, NUM_CUSTOMERS);
    
    Child1.seq_node[0] = 0;
    Child2.seq_node[0] = 0;
    Child3.seq_node[0] = 0;
    Child4.seq_node[0] = 0;
    
    
    int seq_choosen[num_node];
    for(i = 0; i < num_node; i++)
    {
        seq_choosen[i] = i;
    }
    
    int count_pool = -1;
    
    for(int s = 0; s < (int)(NUM_SOL * 0.7) / 2; s++)
    {
        // select parent
        int parent_1 = (int)(rand() % ((int)(NUM_SOL * 0.7) - s * 2));
        int temp = Parent_Pool[parent_1];
        Parent_Pool[parent_1] = Parent_Pool[(int)(NUM_SOL * 0.7) - 1 - s * 2];
        Parent_Pool[(int)(NUM_SOL * 0.7) - 1 - s * 2] = temp;
        parent_1 = temp;
        
        int parent_2 = (int)(rand() % ((int)(NUM_SOL * 0.7) - s * 2 - 1));
        temp = Parent_Pool[parent_2];
        Parent_Pool[parent_2] = Parent_Pool[(int)(NUM_SOL * 0.7) - 2 - s * 2];
        Parent_Pool[(int)(NUM_SOL * 0.7) - 2 - s * 2] = temp;
        parent_2 = temp;
        
        /* ====== in permutation order 1 ====== */
        /// CHILD 1
        U.Permutation_Order_1(Solutions[parent_1], Solutions[parent_2], child_1, num_node);
        // CHILD 2
        U.Permutation_Order_1(Solutions[parent_2], Solutions[parent_1], child_2, num_node);
        
        /* ====== in cycle ====== */
        // CHILD 3
        U.Cycle(Solutions[parent_1], Solutions[parent_2], child_3, num_node);
        //CHILD 4
        U.Cycle(Solutions[parent_2], Solutions[parent_1], child_4, num_node);
        
        bool is_choose_child_1 = U.choose_better_sol(child_1, child_2, Distances, NUM_CUSTOMERS, NUM_VEHICLES);
        bool is_choose_child_3 = U.choose_better_sol(child_3, child_4, Distances, NUM_CUSTOMERS, NUM_VEHICLES);
        
        // Education 2 better sol
        count_pool++;
        if(is_choose_child_1){
            Is_Feasible_Pool[count_pool] = U.Education(child_1, Though_Stattion_In_Pool[count_pool], Distances, NUM_CUSTOMERS, NUM_VEHICLES, MAX_ENERGY_VH, ENG_CONSUMTION, Best_Station, Best_Station_Distances);
            
            U.save_sol_to_pool(child_1, Route_In_Pool[count_pool], NUM_CUSTOMERS + NUM_VEHICLES);
        } else
        {
            Is_Feasible_Pool[count_pool] = U.Education(child_2, Though_Stattion_In_Pool[count_pool], Distances, NUM_CUSTOMERS, NUM_VEHICLES, MAX_ENERGY_VH, ENG_CONSUMTION, Best_Station, Best_Station_Distances);
            
            U.save_sol_to_pool(child_2, Route_In_Pool[count_pool], NUM_CUSTOMERS + NUM_VEHICLES);
        }
        // save to Pool
        
        count_pool++;
        if(is_choose_child_3)
        {
            Is_Feasible_Pool[count_pool] = U.Education(child_3, Though_Stattion_In_Pool[count_pool], Distances, NUM_CUSTOMERS, NUM_VEHICLES, MAX_ENERGY_VH, ENG_CONSUMTION, Best_Station, Best_Station_Distances);
            
            U.save_sol_to_pool(child_3, Route_In_Pool[count_pool], NUM_CUSTOMERS + NUM_VEHICLES);
        } else
        {
            Is_Feasible_Pool[count_pool] = U.Education(child_4, Though_Stattion_In_Pool[count_pool], Distances, NUM_CUSTOMERS, NUM_VEHICLES, MAX_ENERGY_VH, ENG_CONSUMTION, Best_Station, Best_Station_Distances);
            
            U.save_sol_to_pool(child_4, Route_In_Pool[count_pool], NUM_CUSTOMERS + NUM_VEHICLES);
        }
    }
    // create new pool
    int index_pool = -1;
    for(int i = ((int)(NUM_SOL * 0.3)); i < NUM_SOL; i++)
    {
        index_pool++;
        Solutions[index_pool].is_feasible = Is_Feasible_Pool[index_pool];
        for(int j = 0; j < NUM_CUSTOMERS + NUM_VEHICLES; j++)
        {
            Solutions[i].seq_node[j] = Route_In_Pool[index_pool][j];
            Solutions[i].is_though_stat[j] = Though_Stattion_In_Pool[index_pool][j];
            Solutions[i].compute_cost(Distances, Best_Station_Distances);
            Solutions[i].compute_over_cap(Demands, MAX_CAPACITY_VH);
        }
    }
}

void tune_result(){
    for(int i = 0; i < NUM_SOL; i++)
    {
        Solutions[i].is_feasible =  U.find_through_station(Solutions[i].seq_node, Solutions[i].is_though_stat, NUM_CUSTOMERS, NUM_VEHICLES, MAX_ENERGY_VH, ENG_CONSUMTION, Distances, Best_Station, Best_Station_Distances);
        Solutions[i].compute_over_cap(Demands, MAX_CAPACITY_VH);
        Solutions[i].compute_cost(Distances, Best_Station_Distances);
        //printf("\nis_feasible: %d - over_cap: %lf", Solutions[i].is_feasible, Solutions[i].over_capacity);
    }
}

int main(int argc, const char * argv[]) {
    auto start = high_resolution_clock::now();
    FILE *fp;
    fp = fopen("result.txt", "a");
    srand((unsigned)time(NULL));
    read_file((char *)"E-n22-k4.evrp");
    prepare_data();
    init_population();
    compute_meta_data();
    
    build_Roulette_wheel_arr();
    select_parent_to_pool_distinct();
    cross_over();
    tune_result();
    compute_meta_data();
    
//    for(int i = 0; i < 100; i++)
//    {
//        build_Roulette_wheel_arr();
//        select_parent_to_pool_distinct();
//        cross_over();
//        tune_result();
//        compute_meta_data();
//    }
//
//    auto stop = high_resolution_clock::now();
//    auto duration = duration_cast<microseconds>(stop - start);
//
//    int best_sol = -1;
//    double best_cost = 10000000;
//    for(int i = 0; i < NUM_SOL; i++)
//    {
//        if(Solutions[i].cost < best_cost && Solutions[i].is_feasible)
//        {
//            best_sol = i;
//        }
//    }
//
//    printf("\n BEST SOL: %d", best_sol);
//    fprintf(fp, "\ntime_run: %llu milliseconds", duration.count() / 1000);
    fclose(fp);
    return 0;
}
