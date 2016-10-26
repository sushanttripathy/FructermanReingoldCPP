//
// Created by sushant on 10/24/16.
//

#ifndef FR_FRUCTERMANREINGOLD_H
#define FR_FRUCTERMANREINGOLD_H

#include <cstdint>
#include <vector>
#include <map>
#include <random>
#include <math.h>
#include <iostream>
#include <fstream>


#define MIN_DIFF 0.001
#define MIN_DISTANCE MIN_DIFF*1.4142

struct Edge {
    int64_t source, target;
    float weight;

    Edge(int64_t source = 0, int64_t target = 0, float weight = 1.0);
};


struct Displacements_2d {
    float x_disp, y_disp;

    Displacements_2d operator+(const Displacements_2d &operand) const;

    Displacements_2d operator-(const Displacements_2d &operand) const;

    Displacements_2d &operator+=(const Displacements_2d &operand);

    Displacements_2d &operator-=(const Displacements_2d &operand);

    Displacements_2d operator*(float operand) const;

    Displacements_2d operator*(double operand) const;

    Displacements_2d(float x_disp = 0, float y_disp = 0);

    float getMagnitude() const;
};

struct Coords_2d {
    float x, y;

    Coords_2d(float x = 0, float y = 0);

    Coords_2d operator+(const Displacements_2d &operand) const;

    Coords_2d &operator+=(const Displacements_2d &operand);

    Displacements_2d operator-(const Coords_2d &operand) const;
};


class FructermanReingold {
public:

    FructermanReingold(int width, int height);

    FructermanReingold(int width, int height, std::vector<Edge> edge_list, bool bidirectional = true);

    void addEdge(int64_t source, int64_t target, float weight, bool bidirectional = true);

    void addNodeLabel(int64_t node_id, std::string node_label);

    void addNodeWeight(int64_t node_id, float node_weight);

    void randomGraphInitialize();

    void setIgnoreEdgeWeights(bool value);

    void setVerbose(int level);

    void run(long num_iterations, float init_temperature = -1.0, float temperature_decay_constant = 0.9);

    void saveCSV(std::string file_path);

protected:
    int width, height;
    float k; //Optimal pairwise distance
    float k_squared;
    //std::vector <Edge> edge_list;

    std::map<int64_t, std::map<int64_t, float>> edges;
    std::map<int64_t, std::string> node_labels;
    std::map<int64_t, Coords_2d> node_coords;
    std::map<int64_t, Displacements_2d> node_displacements;
    std::map<int64_t, float> node_weights;

    bool weigh_repulsive_force;
    bool weigh_attractive_force;
    int verbose;
};


#endif //FR_FRUCTERMANREINGOLD_H
