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
#include "Customer.hpp"
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
double best_sol = 1000000.0;

// --- meta data ---
//int NUM_SOL = 100;
#define NUM_SOL 100
double H_VAL = 0.0;
double ALPHA = 0.0;

int *list_farthest_seeds;
int *list_nearest_seeds;
int *list_rand_seeds;

// data for best solution
double Best_Cost = 10000000.0;
int Index_Best = -1;

// Util init
Util U;

// data for create pool
double virtual_fitness[NUM_SOL];
int virtual_index[NUM_SOL];
double Roulette_Wheel_Arr[NUM_SOL];
bool Is_In_Pool[NUM_SOL];
int Parent_Pool[NUM_SOL];

//Solution Parent_Pool_Sol[NUM_SOL];
int Index_Best_Sols[(int)(NUM_SOL * 0.3)];

// data for nearest customers

Customer *List_Customers;
int **List_Nearest_Cus;

/* -------------------------------------------------------------------------------------------------------------------------------- */

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
    List_Nearest_Cus = (int **)malloc(NUM_CUSTOMERS * sizeof(int*));
    
    List_Customers = (Customer*) malloc(NUM_CUSTOMERS * sizeof(Customer));
    
    for(i = 0; i < DIMENTION; i++){
        Coords[i] = (double *)malloc(2 * sizeof(double));
        Distances[i] = (double *)malloc(DIMENTION * sizeof(double));
        if(i < NUM_CUSTOMERS)
        {
            Best_Station[i] = (int *)malloc(NUM_CUSTOMERS * sizeof(int));
            Best_Station_Distances[i] = (double *)malloc(NUM_CUSTOMERS * sizeof(double));
            List_Nearest_Cus[i] = (int *)malloc(10 * sizeof(int));
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
        List_Customers[i].x = x;
        List_Customers[i].y = y;
        List_Customers[i].id = i;
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


void mutation()
{
    int i, begin_rand, end_rand;
    int num_node = NUM_CUSTOMERS + NUM_VEHICLES;
    int seq_choosen[NUM_CUSTOMERS + NUM_VEHICLES];
    int seq[num_node];
    for(i = 0; i < NUM_CUSTOMERS + NUM_VEHICLES; i++)
    {
        seq_choosen[i] = i;
    }
    
    int parent = (int)(rand() % (int)(NUM_SOL * 0.3)) + 1;
    for(int i = 0; i < num_node; i++)
        seq[i] = Solutions[parent].seq_node[i];
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
        int temp = begin_rand;
        begin_rand = end_rand;
        end_rand = temp;
    }
    
//    temp = seq[begin_rand];
//    seq[begin_rand] = seq[end_rand];
//    seq[end_rand] = temp;
    
    for(int i = 0; i < (end_rand - begin_rand)/2; i++)
    {
        int temp = seq[begin_rand + i];
        seq[begin_rand + i] = seq[end_rand - i];
        seq[end_rand - i] = temp;
    }
    
    Solution new_sol;
    new_sol.init_mem_space(NUM_VEHICLES, NUM_CUSTOMERS);
    for(int i = 0; i < num_node; i++)
        new_sol.seq_node[i] = seq[i];
    
    new_sol.is_feasible = U.find_through_station(new_sol.seq_node, new_sol.is_though_stat, NUM_CUSTOMERS, NUM_VEHICLES, MAX_ENERGY_VH, ENG_CONSUMTION, Distances, Best_Station, Best_Station_Distances);
    new_sol.compute_cost(Distances, Best_Station_Distances);
    new_sol.compute_over_cap(Demands, MAX_CAPACITY_VH);
    if(Solutions[parent].over_capacity > 0.0) new_sol.is_feasible = false;
    new_sol.fitness = new_sol.cost + ALPHA * new_sol.over_capacity;
    Solutions[parent] = new_sol;
}


// compare Customer's distances
int compare_customer_distances(const void *ptr1, const void *ptr2)
{
    Customer c1 = *(Customer*)ptr1;
    Customer c2 = *(Customer*)ptr2;
    
    if(c1.temp_distance > c2.temp_distance) return 1;
    else if(c1.temp_distance < c2.temp_distance) return  -1;
    
    return 0;
}

// find list nearest customers for each customer
void find_nearest_cus() {
    for(int i = 1; i < NUM_CUSTOMERS; i++)
    {
        // compute temp distance
        for(int j = 0; j < NUM_CUSTOMERS; j++)
        {
            List_Customers[j].temp_distance = Distances[j][i];
            List_Customers[j].id = j;
        }
        // sort list customers
        qsort(List_Customers, NUM_CUSTOMERS, sizeof(Customer), compare_customer_distances);
        // copy id to list nearest cus
        int index = -1;
        for(int k = 0; k < 10; k++)
        {
            index++;
            if(List_Customers[index].id == i || List_Customers[index].id == 0) index++;
            List_Nearest_Cus[i][k] = List_Customers[index].id;
        }
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
    find_nearest_cus();
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
    Solutions[i].is_feasible = true;
    Solutions[i].Potvin_init_test(list_farthest_seeds, Distances, Demands, MAX_CAPACITY_VH);
    Solutions[i].is_feasible = U.find_through_station(Solutions[i].seq_node, Solutions[i].is_though_stat, NUM_CUSTOMERS, NUM_VEHICLES, MAX_ENERGY_VH, ENG_CONSUMTION, Distances, Best_Station, Best_Station_Distances);
    Solutions[i].compute_cost(Distances, Best_Station_Distances);
    Solutions[i].compute_over_cap(Demands, MAX_CAPACITY_VH);
    if(Solutions[i].over_capacity > 0.0) Solutions[i].is_feasible = false;
    
    // init Potvin route with nearest seeds
    i++;
    Solutions[i].is_feasible = true;
    Solutions[i].Potvin_init_test(list_nearest_seeds, Distances, Demands, MAX_CAPACITY_VH);
    Solutions[i].is_feasible = U.find_through_station(Solutions[i].seq_node, Solutions[i].is_though_stat, NUM_CUSTOMERS, NUM_VEHICLES, MAX_ENERGY_VH, ENG_CONSUMTION, Distances, Best_Station, Best_Station_Distances);
    Solutions[i].compute_cost(Distances, Best_Station_Distances);
    Solutions[i].compute_over_cap(Demands, MAX_CAPACITY_VH);
    if(Solutions[i].over_capacity > 0.0) Solutions[i].is_feasible = false;
    i++;
    for(; i < NUM_SOL; i++)
    {
        Solutions[i].is_feasible = true;
        random_select_seeds();
        Solutions[i].Potvin_init_test(list_rand_seeds, Distances, Demands, MAX_CAPACITY_VH);
        Solutions[i].is_feasible = U.find_through_station(Solutions[i].seq_node, Solutions[i].is_though_stat, NUM_CUSTOMERS, NUM_VEHICLES, MAX_ENERGY_VH, ENG_CONSUMTION, Distances, Best_Station, Best_Station_Distances);
        Solutions[i].compute_cost(Distances, Best_Station_Distances);
        Solutions[i].compute_over_cap(Demands, MAX_CAPACITY_VH);
        if(Solutions[i].over_capacity > 0.0) Solutions[i].is_feasible = false;
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

int compare_sols_cost(const void *ptr1, const void *ptr2)
{
    Solution m1 = *(Solution*)ptr1;
    Solution m2 = *(Solution*)ptr2;
    
    if(m1.cost  > m2.cost) return 1;
    else if(m1.cost < m2.cost) return -1;
    
    return 0;
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
    for(i = 1; i < NUM_SOL; i++)
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
        Parent_Pool[i] = -1;
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

// select best feasible sol
void select_best_feasible_sol()
{
    int latest_idx = -1;
    for(int i = 0; i < NUM_SOL; i++)
    {
        if(Solutions[i].is_feasible)
        {
            latest_idx++;
            Solution temp = Solutions[latest_idx];
            Solutions[latest_idx] = Solutions[i];
            Solutions[i] = temp;
        }
    }
    qsort(Solutions, latest_idx + 1, sizeof(Solution), compare_sols_cost);
    qsort(&Solutions[latest_idx + 1], NUM_SOL - latest_idx - 1, sizeof(Solution), compare_sols_fitness);
}

// inspect sol
void inspect_sol(Solution sol, int index) {
    printf("\n %d - COST: %lf - OVER_CAP: %lf - IS_FEASIBLE: %d - FITNESS: %lf - ALPHA: %lf",index, sol.cost, sol.over_capacity, sol.is_feasible, sol.fitness, ALPHA);
}

// Begin Cross over
void cross_over(){
    Solution Parent_Pool_Sol[NUM_SOL];
    // create child 1 and 2 space
    int num_node = NUM_CUSTOMERS + NUM_VEHICLES;

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
        U.Permutation_Order_1(Solutions[parent_1], Solutions[parent_2], Child1.seq_node, num_node);
        U.Permutation_Order_1(Solutions[parent_2], Solutions[parent_1], Child2.seq_node, num_node);
        
        /* ====== in cycle ====== */
        U.Cycle(Solutions[parent_1], Solutions[parent_2], Child3.seq_node, num_node);
        U.Cycle(Solutions[parent_2], Solutions[parent_1], Child4.seq_node, num_node);

        // compute meta data for 4 childs
        
        //CHILD 1
        Child1.is_feasible = U.find_through_station(Child1.seq_node, Child1.is_though_stat, NUM_CUSTOMERS, NUM_VEHICLES, MAX_ENERGY_VH, ENG_CONSUMTION, Distances, Best_Station, Best_Station_Distances);
        Child1.compute_cost(Distances, Best_Station_Distances);
        Child1.compute_over_cap(Demands, MAX_CAPACITY_VH);
        if(Child1.over_capacity > 0.0) Child1.is_feasible = false;
        Child1.fitness = Child1.cost + ALPHA * Child1.over_capacity;
        
        // CHILD 2
        Child2.is_feasible = U.find_through_station(Child2.seq_node, Child2.is_though_stat, NUM_CUSTOMERS, NUM_VEHICLES, MAX_ENERGY_VH, ENG_CONSUMTION, Distances, Best_Station, Best_Station_Distances);
        Child2.compute_cost(Distances, Best_Station_Distances);
        Child2.compute_over_cap(Demands, MAX_CAPACITY_VH);
        if(Child2.over_capacity > 0.0) Child2.is_feasible = false;
        Child2.fitness = Child2.cost + ALPHA * Child2.over_capacity;
        
        // CHILD 3
        Child3.is_feasible = U.find_through_station(Child3.seq_node, Child3.is_though_stat, NUM_CUSTOMERS, NUM_VEHICLES, MAX_ENERGY_VH, ENG_CONSUMTION, Distances, Best_Station, Best_Station_Distances);
        Child3.compute_cost(Distances, Best_Station_Distances);
        Child3.compute_over_cap(Demands, MAX_CAPACITY_VH);
        if(Child3.over_capacity > 0.0) Child3.is_feasible = false;
        Child3.fitness = Child3.cost + ALPHA * Child3.over_capacity;
        
        // CHILD 4
        Child4.is_feasible = U.find_through_station(Child4.seq_node, Child4.is_though_stat, NUM_CUSTOMERS, NUM_VEHICLES, MAX_ENERGY_VH, ENG_CONSUMTION, Distances, Best_Station, Best_Station_Distances);
        Child4.compute_cost(Distances, Best_Station_Distances);
        Child4.compute_over_cap(Demands, MAX_CAPACITY_VH);
        if(Child4.over_capacity > 0.0) Child4.is_feasible = false;
        Child4.fitness = Child4.cost + ALPHA * Child4.over_capacity;
        
        bool is_choose_child_1 = false;
        bool is_choose_child_3 = false;
        if(Child1.fitness < Child2.fitness) is_choose_child_1 = true;
        if(Child3.fitness < Child4.fitness) is_choose_child_3 = true;
        
        // Education 2 better sol
        count_pool++;
        if(is_choose_child_1) {
            Parent_Pool_Sol[count_pool].init_mem_space(NUM_VEHICLES, NUM_CUSTOMERS);
            Parent_Pool_Sol[count_pool].copy_seq_node(Child1);
            Parent_Pool_Sol[count_pool].cost = Child1.cost;
            Parent_Pool_Sol[count_pool].fitness = Child1.fitness;
            
            
            
            U.Tabu_search(Parent_Pool_Sol[count_pool].seq_node, Distances, NUM_CUSTOMERS, NUM_VEHICLES, Demands, Parent_Pool_Sol[count_pool].fitness, List_Nearest_Cus, Parent_Pool_Sol[count_pool].cost, MAX_CAPACITY_VH, ALPHA, MAX_ENERGY_VH, ENG_CONSUMTION, Best_Station, Best_Station_Distances);
            
            Parent_Pool_Sol[count_pool].is_feasible = U.find_through_station(Parent_Pool_Sol[count_pool].seq_node, Parent_Pool_Sol[count_pool].is_though_stat, NUM_CUSTOMERS, NUM_VEHICLES, MAX_ENERGY_VH, ENG_CONSUMTION, Distances, Best_Station, Best_Station_Distances);
            Parent_Pool_Sol[count_pool].compute_cost(Distances, Best_Station_Distances);
            Parent_Pool_Sol[count_pool].compute_over_cap(Demands, MAX_CAPACITY_VH);
            if(Parent_Pool_Sol[count_pool].over_capacity > 0.0) Parent_Pool_Sol[count_pool].is_feasible = false;
            Parent_Pool_Sol[count_pool].fitness = Parent_Pool_Sol[count_pool].cost + ALPHA * Parent_Pool_Sol[count_pool].over_capacity;
        } else {
            Parent_Pool_Sol[count_pool].init_mem_space(NUM_VEHICLES, NUM_CUSTOMERS);
            Parent_Pool_Sol[count_pool].copy_seq_node(Child2);
            Parent_Pool_Sol[count_pool].cost = Child2.cost;
            Parent_Pool_Sol[count_pool].fitness = Child2.fitness;
            
            
            
            U.Tabu_search(Parent_Pool_Sol[count_pool].seq_node, Distances, NUM_CUSTOMERS, NUM_VEHICLES, Demands, Parent_Pool_Sol[count_pool].fitness, List_Nearest_Cus, Parent_Pool_Sol[count_pool].cost, MAX_CAPACITY_VH, ALPHA, MAX_ENERGY_VH, ENG_CONSUMTION, Best_Station, Best_Station_Distances);
            
            Parent_Pool_Sol[count_pool].is_feasible = U.find_through_station(Parent_Pool_Sol[count_pool].seq_node, Parent_Pool_Sol[count_pool].is_though_stat, NUM_CUSTOMERS, NUM_VEHICLES, MAX_ENERGY_VH, ENG_CONSUMTION, Distances, Best_Station, Best_Station_Distances);
            Parent_Pool_Sol[count_pool].compute_cost(Distances, Best_Station_Distances);
            Parent_Pool_Sol[count_pool].compute_over_cap(Demands, MAX_CAPACITY_VH);
            if(Parent_Pool_Sol[count_pool].over_capacity > 0.0) Parent_Pool_Sol[count_pool].is_feasible = false;
            Parent_Pool_Sol[count_pool].fitness = Parent_Pool_Sol[count_pool].cost + ALPHA * Parent_Pool_Sol[count_pool].over_capacity;
        }
        
        count_pool++;
        if(is_choose_child_3) {
            Parent_Pool_Sol[count_pool].init_mem_space(NUM_VEHICLES, NUM_CUSTOMERS);
            Parent_Pool_Sol[count_pool].copy_seq_node(Child3);
            Parent_Pool_Sol[count_pool].cost = Child3.cost;
            Parent_Pool_Sol[count_pool].fitness = Child3.fitness;
            
            
            
            U.Tabu_search(Parent_Pool_Sol[count_pool].seq_node, Distances, NUM_CUSTOMERS, NUM_VEHICLES, Demands, Parent_Pool_Sol[count_pool].fitness, List_Nearest_Cus, Parent_Pool_Sol[count_pool].cost, MAX_CAPACITY_VH, ALPHA, MAX_ENERGY_VH, ENG_CONSUMTION, Best_Station, Best_Station_Distances);
            
            Parent_Pool_Sol[count_pool].is_feasible = U.find_through_station(Parent_Pool_Sol[count_pool].seq_node, Parent_Pool_Sol[count_pool].is_though_stat, NUM_CUSTOMERS, NUM_VEHICLES, MAX_ENERGY_VH, ENG_CONSUMTION, Distances, Best_Station, Best_Station_Distances);
            Parent_Pool_Sol[count_pool].compute_cost(Distances, Best_Station_Distances);
            Parent_Pool_Sol[count_pool].compute_over_cap(Demands, MAX_CAPACITY_VH);
            if(Parent_Pool_Sol[count_pool].over_capacity > 0.0) Parent_Pool_Sol[count_pool].is_feasible = false;
            Parent_Pool_Sol[count_pool].fitness = Parent_Pool_Sol[count_pool].cost + ALPHA * Parent_Pool_Sol[count_pool].over_capacity;
        }
        else {
            Parent_Pool_Sol[count_pool].init_mem_space(NUM_VEHICLES, NUM_CUSTOMERS);
            Parent_Pool_Sol[count_pool].copy_seq_node(Child4);
            Parent_Pool_Sol[count_pool].cost = Child4.cost;
            Parent_Pool_Sol[count_pool].fitness = Child4.fitness;
            
            
            
            U.Tabu_search(Parent_Pool_Sol[count_pool].seq_node, Distances, NUM_CUSTOMERS, NUM_VEHICLES, Demands, Parent_Pool_Sol[count_pool].fitness, List_Nearest_Cus, Parent_Pool_Sol[count_pool].cost, MAX_CAPACITY_VH, ALPHA, MAX_ENERGY_VH, ENG_CONSUMTION, Best_Station, Best_Station_Distances);
            
            Parent_Pool_Sol[count_pool].is_feasible = U.find_through_station(Parent_Pool_Sol[count_pool].seq_node, Parent_Pool_Sol[count_pool].is_though_stat, NUM_CUSTOMERS, NUM_VEHICLES, MAX_ENERGY_VH, ENG_CONSUMTION, Distances, Best_Station, Best_Station_Distances);
            Parent_Pool_Sol[count_pool].compute_cost(Distances, Best_Station_Distances);
            Parent_Pool_Sol[count_pool].compute_over_cap(Demands, MAX_CAPACITY_VH);
            if(Parent_Pool_Sol[count_pool].over_capacity > 0.0) Parent_Pool_Sol[count_pool].is_feasible = false;
            Parent_Pool_Sol[count_pool].fitness = Parent_Pool_Sol[count_pool].cost + ALPHA * Parent_Pool_Sol[count_pool].over_capacity;
        }
    }
//    printf("\nParent Pool\n");
//    for(int s = 0; s < NUM_CUSTOMERS + NUM_VEHICLES; s++)
//    {
//        printf("\n\n\n\n=====");
//        if(Parent_Pool_Sol[s].is_feasible) inspect_sol(Parent_Pool_Sol[s]);
//        Parent_Pool_Sol[s].is_feasible = U.find_through_station(Parent_Pool_Sol[s].seq_node, Parent_Pool_Sol[s].is_though_stat, NUM_CUSTOMERS, NUM_VEHICLES, MAX_ENERGY_VH, ENG_CONSUMTION, Distances, Best_Station, Best_Station_Distances);
//        Parent_Pool_Sol[s].compute_cost(Distances, Best_Station_Distances);
//        Parent_Pool_Sol[s].compute_over_cap(Demands, MAX_CAPACITY_VH);
//        if(Parent_Pool_Sol[s].over_capacity > 0.0) Parent_Pool_Sol[s].is_feasible = false;
//        Parent_Pool_Sol[s].fitness = Parent_Pool_Sol[s].cost + ALPHA * Parent_Pool_Sol[s].over_capacity;
//        inspect_sol(Parent_Pool_Sol[s]);
//        printf("\n\n\n------");
//    }
    // create new pool
    select_best_feasible_sol();
    int index_pool = -1;
    for(int i = ((int)(NUM_SOL * 0.3)); i < NUM_SOL; i++)
    {
        index_pool++;
        Solutions[i].init_mem_space(NUM_VEHICLES, NUM_CUSTOMERS);
        Solutions[i].copy_seq_node(Parent_Pool_Sol[index_pool]);
        //inspect_sol(Parent_Pool_Sol[index_pool]);
        Solutions[i].is_feasible = U.find_through_station(Solutions[i].seq_node, Solutions[i].is_though_stat, NUM_CUSTOMERS, NUM_VEHICLES, MAX_ENERGY_VH, ENG_CONSUMTION, Distances, Best_Station, Best_Station_Distances);
        Solutions[i].compute_cost(Distances, Best_Station_Distances);
        Solutions[i].compute_over_cap(Demands, MAX_CAPACITY_VH);
        if(Solutions[i].over_capacity > 0.0) Solutions[i].is_feasible = false;
        Solutions[i].fitness = Solutions[i].cost + ALPHA * Solutions[i].over_capacity;
        //inspect_sol(Solutions[i]);
    }
}



void tune_result(){
    for(int i = 0; i < NUM_SOL; i++)
    {
        Solutions[i].is_feasible =  U.find_through_station(Solutions[i].seq_node, Solutions[i].is_though_stat, NUM_CUSTOMERS, NUM_VEHICLES, MAX_ENERGY_VH, ENG_CONSUMTION, Distances, Best_Station, Best_Station_Distances);
        Solutions[i].compute_over_cap(Demands, MAX_CAPACITY_VH);
        if(Solutions[i].over_capacity > 0.0) Solutions[i].is_feasible = false;
        Solutions[i].compute_cost(Distances, Best_Station_Distances);
    }
}

void local_optimize(int *seq, double init_fitness, double init_cost)
{
    int list_rand[3] = {0,1,2};
    int rand1 = rand() % 3;
    int temp = list_rand[0];
    list_rand[0] = list_rand[rand1];
    list_rand[rand1] = temp;
    int rand2 = rand() % 3;
    temp = list_rand[0];
    list_rand[0] = list_rand[rand2];
    list_rand[rand2] = temp;
    
    for(int i = 0; i < 3; i++)
    {
        switch (i) {
            case 0:
                U.Two_Exchange(seq, Distances, NUM_CUSTOMERS, NUM_VEHICLES, Demands, init_fitness, List_Nearest_Cus, init_cost, MAX_CAPACITY_VH, ALPHA, MAX_ENERGY_VH, ENG_CONSUMTION, Best_Station, Best_Station_Distances);
                break;
            case 1:
                U.Or_Exchange(seq, Distances, NUM_CUSTOMERS, NUM_VEHICLES, Demands, init_fitness, List_Nearest_Cus, init_cost, MAX_CAPACITY_VH, ALPHA, MAX_ENERGY_VH, ENG_CONSUMTION, Best_Station, Best_Station_Distances);
                break;
            case 2:
                U.Cross_Exchange(seq, Distances, NUM_CUSTOMERS, NUM_VEHICLES, Demands, init_fitness, List_Nearest_Cus, init_cost, MAX_CAPACITY_VH, ALPHA, MAX_ENERGY_VH, ENG_CONSUMTION, Best_Station, Best_Station_Distances);
            default:
                break;
        }
    }
}


void create_space_mem_cvrp()
{
    int i = 0;
    Distances = (double **)malloc(DIMENTION * sizeof(double *));
    List_Nearest_Cus = (int **)malloc(DIMENTION * sizeof(int*));
    for(i = 0; i < DIMENTION; i++)
    {
        Distances[i] = (double *)malloc(DIMENTION * sizeof(double));
        List_Nearest_Cus[i] = (int *)malloc(10 * sizeof(int));
    }
}

void read_baygn29k4(char *file_src)
{
    int i, j;
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
    fscanf(infile, "%lf\n", &MAX_CAPACITY_VH);
    NUM_CUSTOMERS = DIMENTION;
    
    printf("data: dimention: %d - capacity_vh: %lf - num_vehicles: %d - num_customer: %d\n", DIMENTION, MAX_CAPACITY_VH, NUM_VEHICLES, NUM_CUSTOMERS);
    
    create_space_mem();
    
    double weight = 0.0;
    for(i = 0; i < DIMENTION; i++)
    {
        for(j = i + 1; j < DIMENTION; j++)
        {
            fscanf(infile,"%lf",&weight);
            Distances[i][j] = weight;
            Distances[j][i] = weight;
            Distances[i][i] = 0.0;
        }
    }
    double index;
    double demand = 0.0;
    for(i = 0; i < DIMENTION; i++)
    {
        fscanf(infile, "%lf %lf\n", &index, &demand);
        Demands[i] = demand;
    }
    
    find_nearest_cus();
    fclose(infile);
}

void Tabu_search_Test()
{
    FILE *fp;
    fp = fopen("./Result/tabu_search_baygn_29_k4.txt", "a");
    
    
    Solution sol;
    sol.init_mem_space(NUM_VEHICLES, NUM_CUSTOMERS);
    sol.is_feasible = true;
    random_select_seeds();
    //find_farthest_seed();
    
    sol.Potvin_init_test(list_rand_seeds, Distances, Demands, MAX_CAPACITY_VH);
    for(int s = 0; s < NUM_CUSTOMERS + NUM_VEHICLES; s++)
        sol.is_though_stat[s]= false;
    sol.compute_cost(Distances, Best_Station_Distances);
    sol.compute_over_cap(Demands, MAX_CAPACITY_VH);
    if(sol.over_capacity > 0.0) sol.is_feasible = false;
    sol.fitness = sol.cost + sol.over_capacity;
    
    inspect_sol(sol, 1);
    printf("\n BEGIN ROUTE: ");
    for(int s = 0; s < NUM_CUSTOMERS + NUM_VEHICLES; s++)
        printf("%d -> ", sol.seq_node[s]);
    
    U.Tabu_search_cvrp(sol.seq_node, Distances, NUM_CUSTOMERS, NUM_VEHICLES, Demands, sol.fitness, List_Nearest_Cus, sol.cost, MAX_CAPACITY_VH, 1.0, MAX_ENERGY_VH, ENG_CONSUMTION, Best_Station, Best_Station_Distances);
    
    sol.is_feasible = true;
    sol.compute_cost(Distances, Best_Station_Distances);
    sol.compute_over_cap(Demands, MAX_CAPACITY_VH);
    if(sol.over_capacity > 0.0) sol.is_feasible = false;
    sol.fitness = sol.cost + sol.over_capacity;
    
    inspect_sol(sol, 1);
    
    printf("\n");
    for(int i = 0; i < DIMENTION + NUM_VEHICLES; i++)
        printf("%d -> ", sol.seq_node[i]);
    
     fprintf(fp, "\nis_feasible: %d - cost: %lf - over_cap: %lf - Fitness: %lf\n", sol.is_feasible, sol.cost, sol.over_capacity, sol.fitness);
    for(int s = 0; s < NUM_CUSTOMERS + NUM_VEHICLES; s++)
        fprintf(fp, "%d -> ", sol.seq_node[s]);
     fprintf(fp, "\n");
    for(int s = 0; s < NUM_CUSTOMERS + NUM_VEHICLES;s++)
        fprintf(fp, "%d -> ", sol.is_though_stat[s]);
    fclose(fp);
    
}

int main(int argc, const char * argv[]) {
    srand((unsigned)time(NULL));
    read_baygn29k4((char *)"./Data/CVRP/bayg-n29-k4.vrp");
    Tabu_search_Test();
    return 0;
}
