#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <limits>
#include <fstream>

const int N = 100; 
const int M = 30;  
const int MAX_ITER = 1000;
const double ALPHA = 2.0;
const double BETA = 4.0;
const double RHO = 0.4;
const int STEP = 20;
const int MAX_DISTANCE = 50;
const int MIN_DISTANCE = 5;

double distances[N][N]; 
std::vector<int> bestPath(N); 

double pheromones[N][N];
double bestCost = std::numeric_limits<double>::infinity();

void initializeDistances() {
    srand(time(0));
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            distances[i][j] = distances[j][i] = MIN_DISTANCE + rand() % (MAX_DISTANCE - MIN_DISTANCE + 1);
        }
        distances[i][i] = 0;
    }
}

double greedySolution(std::vector<int>& path) {
    std::vector<bool> visited(N, false);
    path.clear();
    path.push_back(0);
    visited[0] = true;
    double totalCost = 0;

    for (int step = 1; step < N; ++step) {
        int current = path.back();
        int next = -1;
        double minDistance = std::numeric_limits<double>::infinity();

        for (int j = 0; j < N; ++j) {
            if (!visited[j] && distances[current][j] < minDistance) {
                next = j;
                minDistance = distances[current][j];
            }
        }

        path.push_back(next);
        visited[next] = true;
        totalCost += minDistance;
    }

    totalCost += distances[path.back()][path[0]];
    return totalCost;
}

void initializePheromones() {
    std::vector<int> path;
    double Lmin = greedySolution(path);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            pheromones[i][j] = 1.0 / (N * Lmin);
        }
    }
}

std::vector<int> generateAntSolution() {
    std::vector<int> path;
    std::vector<bool> visited(N, false);
    int start = rand() % N;
    path.push_back(start);
    visited[start] = true;

    for (int step = 1; step < N; step++) {
        int current = path.back();
        double totalProb = 0;
        std::vector<double> probabilities(N, 0.0);

        for (int j = 0; j < N; ++j) {
            if (!visited[j]) {
                probabilities[j] = pow(pheromones[current][j], ALPHA) * pow(1.0 / distances[current][j], BETA);
                totalProb += probabilities[j];  
            }
        }
        for (int j = 0; j < N; ++j) {
            if (!visited[j]) {
                probabilities[j] /= totalProb;  
            }
        }

        double threshold = ((double)rand() / RAND_MAX);
        double cumulative = 0;
        int next = -1;

        for (int j = 0; j < N; ++j) {
            if (!visited[j]) {
                cumulative += probabilities[j];
                if (cumulative >= threshold) {
                    next = j;
                    break;
                }
            }
        }

        path.push_back(next);
        visited[next] = true;
    }

    return path;
}

double computeCost(const std::vector<int>& path) {
    double totalCost = 0;
    for (int i = 0; i < N - 1; ++i) {
        totalCost += distances[path[i]][path[i + 1]];
    }
    totalCost += distances[path.back()][path[0]];
    return totalCost;
}

void updatePheromones(const std::vector<std::vector<int>>& solutions, const std::vector<double>& costs) {
    std::vector<int> greedyPath;
    double Lmin = greedySolution(greedyPath);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            pheromones[i][j] *= (1 - RHO);
        }
    }

    for (int k = 0; k < M; ++k) {
        for (int i = 0; i < N - 1; ++i) {
            pheromones[solutions[k][i]][solutions[k][i + 1]] += Lmin / costs[k];
        }
        pheromones[solutions[k].back()][solutions[k][0]] += Lmin / costs[k];
    }
}

int main() {
    initializeDistances();
    initializePheromones();

    std::ofstream results("results.txt");

    for (int iter = 0; iter <= MAX_ITER; iter += STEP) {
        std::vector<std::vector<int>> solutions(M);
        std::vector<double> costs(M);

        for (int k = 0; k < M; ++k) {
            solutions[k] = generateAntSolution();
            costs[k] = computeCost(solutions[k]);
            if (costs[k] < bestCost) {
                bestCost = costs[k];
                bestPath = solutions[k];
            }
        }

        updatePheromones(solutions, costs);

            results << "Iteration " << iter << ": Best Cost = " << bestCost << std::endl;
            std::cout << "Iteration " << iter << ": Best Cost = " << bestCost << std::endl;
    }
}