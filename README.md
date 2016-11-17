# Evaluation Ant Colony Optimization
An advanced ant colony optimization for solving travelling salesman problem.

## Background
[Travelling salesman problem](https://en.wikipedia.org/wiki/Travelling_salesman_problem)
[Ant colony optimization algorithms](https://en.wikipedia.org/wiki/Ant_colony_optimization_algorithms)

## Introduction
Suggested a new idea to improve the original ant colony optimization.
Split the TSP map into multiple smaller maps, running solution search on each sub-maps.</br>
After search on each sub-map get converged, combine the ant from different sub-maps and use [genetic algorithm](https://en.wikipedia.org/wiki/Genetic_algorithm) to make the ants search and evaluation by themselves.

## Sample result
This is the sample result of TSP225 problem.</br>
Result of sub-maps:</br>
![submaps_img](https://github.com/lightbluex/Evaluation-Ant-Colony-Optimization/blob/master/sample_result/tsp225_submaps.jpg?raw=true "sub-maps")
</br>
Final result:</br>
![final_result_img](https://github.com/lightbluex/Evaluation-Ant-Colony-Optimization/blob/master/sample_result/tsp225_final_result.jpg?raw=true "sub-maps")

## Getting Started

1. Download TSP problems from [TSP Library](http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/) to your local.

2. Change the paths of TSP in ```AS.c``` to your local path.

3. Run AS.c.
