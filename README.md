## EVRP_STD
This project solve the CVRP and EVRP problem
it cans find the best result as possible compare to the optimal result in each instance problem

## Algorithm
GA(Gennetic Algorithm) + Tabu search(in education phase) to reachieve the best result

### Genetic algorithms(GA)
GA is a metaheuristic inspired by the process of natural selection that belongs to the larger class of evolutionary algorithms (EA). Genetic algorithms are commonly used to generate high-quality solutions to optimization and search problems by relying on biologically inspired operators such as mutation, crossover and selection

### Tabu search
Tabu search enhances the performance of local search by relaxing its basic rule. First, at each step worsening moves can be accepted if no improving move is available (like when the search is stuck at a strict local minimum). In addition, prohibitions (henceforth the term tabu) are introduced to discourage the search from coming back to previously-visited solutions. Tabu search is apply in step of choosing mutation, crossover and selection of GA

## Launch
- Open project in xcode and run
- The program will solve all the problems store in E-n*.evrp files
