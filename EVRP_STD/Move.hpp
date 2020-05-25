//
//  Move.hpp
//  EVRP_STD
//
//  Created by MAC on 5/25/20.
//  Copyright © 2020 MAC. All rights reserved.
//

#ifndef Move_hpp
#define Move_hpp

#include <stdio.h>
class Move {
public:
    int cus_1;
    int cus_2;
    double varCost;
    double varVioCap;
    double varFitness;
    bool feasible = false;
    bool tabu = false;
};
#endif /* Move_hpp */
