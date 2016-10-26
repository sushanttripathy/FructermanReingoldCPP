//
// Created by sushant on 10/24/16.
//

#include "fructermanreingold.h"


Edge::Edge(int64_t source, int64_t target, float weight) {
    this->source = source;
    this->target = target;
    this->weight = weight;
}

Coords_2d::Coords_2d(float x, float y) {
    this->x = x;
    this->y = y;
}

Displacements_2d Coords_2d::operator-(const Coords_2d &operand) const {
    float delta_x = this->x - operand.x;
    float delta_y = this->y - operand.y;
    return Displacements_2d(delta_x, delta_y);
}

Coords_2d Coords_2d::operator+(const Displacements_2d &operand) const {
    float new_x = this->x + operand.x_disp;
    float new_y = this->y + operand.y_disp;

    return Coords_2d(new_x, new_y);
}

Coords_2d &Coords_2d::operator+=(const Displacements_2d &operand) {
    this->x += operand.x_disp;
    this->y += operand.y_disp;
    return *this;
}

Displacements_2d::Displacements_2d(float x_disp, float y_disp) {

    this->x_disp = x_disp;
    this->y_disp = y_disp;
}

Displacements_2d Displacements_2d::operator+(const Displacements_2d &a) const {
    float new_x_disp = this->x_disp + a.x_disp;
    float new_y_disp = this->y_disp + a.y_disp;
    return Displacements_2d(new_x_disp, new_y_disp);
}

Displacements_2d &Displacements_2d::operator+=(const Displacements_2d &operand) {
    this->x_disp += operand.x_disp;
    this->y_disp += operand.y_disp;
    return *this;
}

Displacements_2d Displacements_2d::operator-(const Displacements_2d &a) const {
    float new_x_disp = this->x_disp - a.x_disp;
    float new_y_disp = this->y_disp - a.y_disp;
    return Displacements_2d(new_x_disp, new_y_disp);
}

Displacements_2d &Displacements_2d::operator-=(const Displacements_2d &operand) {
    this->x_disp -= operand.x_disp;
    this->y_disp -= operand.y_disp;
    return *this;
}

Displacements_2d Displacements_2d::operator*(float operand) const {
    float new_x_disp = this->x_disp * operand;
    float new_y_disp = this->y_disp * operand;
    return Displacements_2d(new_x_disp, new_y_disp);
}

Displacements_2d Displacements_2d::operator*(double operand) const {
    float new_x_disp = this->x_disp * operand;
    float new_y_disp = this->y_disp * operand;
    return Displacements_2d(new_x_disp, new_y_disp);
}

float Displacements_2d::getMagnitude() const {
    return sqrt(this->x_disp * this->x_disp + this->y_disp * this->y_disp);
}


FructermanReingold::FructermanReingold(int width, int height) {
    this->width = width;
    this->height = height;
    this->weigh_attractive_force = false;
    this->weigh_repulsive_force = false;
    this->verbose = 0;
}

FructermanReingold::FructermanReingold(int width, int height, std::vector<Edge> edge_list, bool bidirectional) {
    this->width = width;
    this->height = height;
    this->weigh_attractive_force = false;
    this->weigh_repulsive_force = false;
    this->verbose = 0;
}

void FructermanReingold::addEdge(int64_t source, int64_t target, float weight, bool bidirectional) {
    if (bidirectional) {
        if (source < target)
            this->edges[source][target] = weight;
        else
            this->edges[target][source] = weight;
    } else {
        this->edges[source][target] = weight;
    }

    if (weight < 1.0) {
        this->weigh_attractive_force = true;
    }

    if (this->node_coords.find(source) == this->node_coords.end()) {
        this->node_coords[source] = Coords_2d(0, 0);
    }
    if (this->node_coords.find(target) == this->node_coords.end()) {
        this->node_coords[target] = Coords_2d(0, 0);
    }
}

void FructermanReingold::addNodeLabel(int64_t node_id, std::string node_label) {
    this->node_labels[node_id] = node_label;
}


void FructermanReingold::randomGraphInitialize() {
    std::mt19937 rng;
    rng.seed(std::random_device()());
    std::uniform_int_distribution<std::mt19937::result_type> dist_width(0, this->width), dist_height(0, this->height);

    for (std::map<int64_t, Coords_2d>::iterator map_it = this->node_coords.begin();
         map_it != this->node_coords.end(); ++map_it) {
        //std::cout << map_it->first << std::endl;
        map_it->second.x = ((float) dist_width(rng) - this->width / 2.0);
        map_it->second.y = ((float) dist_height(rng) - this->height / 2.0);
    }

    if (this->verbose > 1) {
        for (std::map<int64_t, Coords_2d>::iterator map_it = this->node_coords.begin();
             map_it != this->node_coords.end(); ++map_it) {
            std::cout << map_it->first << " " << map_it->second.x << "," << map_it->second.y << std::endl;

        }
    }

}

void FructermanReingold::run(long num_iterations, float init_temperature, float temperature_decay_constant) {
    this->randomGraphInitialize();

    if (this->node_weights.size() == this->node_coords.size()) {
        this->weigh_repulsive_force &= true;
    }

    if (this->node_coords.size())
        this->k = sqrt(
                float(this->width * this->height) / float(this->node_coords.size()));
    else
        this->k = sqrt(float(this->width * this->height)); // Assume only one node
    this->k_squared = this->k * this->k;

    float temperature = float(std::min(this->width, this->height)) /
                        2.0;//sqrt(this->node_coords.size());//float(std::min(this->width, this->height)) / 2.0;
    if (init_temperature != -1.0)
        temperature = init_temperature;


    float bounds_x_min = -1.0 * this->width / 2.0, bounds_x_max = this->width / 2.0;
    float bounds_y_min = -1.0 * this->height / 2.0, bounds_y_max = this->height / 2.0;

    float linear_cooling_factor = temperature / float(num_iterations);

    for (long i = 0; i < num_iterations; ++i) {

        //Initialize the displacement vectors and add pairwise repulsive forces

        for (std::map<int64_t, Coords_2d>::iterator map_it = this->node_coords.begin();
             map_it != this->node_coords.end(); ++map_it) {
            this->node_displacements[map_it->first] = Displacements_2d(0, 0); //Initialize the displacement vector
            for (std::map<int64_t, Coords_2d>::iterator map_it2 = this->node_coords.begin();
                 map_it2 != this->node_coords.end(); ++map_it2) {

                if (map_it->first != map_it2->first) {
                    Displacements_2d delta = map_it->second - map_it2->second;
                    float delta_mag = 0;


                    delta_mag = delta.getMagnitude();

                    if (delta_mag < MIN_DISTANCE) {
                        delta_mag = MIN_DISTANCE;
                        if (delta.x_disp > 0) {
                            delta.x_disp = MIN_DIFF;
                        } else {
                            delta.x_disp = -1.0 * MIN_DIFF;
                        }
                        if (delta.y_disp > 0) {
                            delta.y_disp = MIN_DIFF;
                        } else {
                            delta.y_disp = -1.0 * MIN_DIFF;
                        }
                    }

                    float repulsive_force = this->k_squared / delta_mag;

                    if (this->weigh_repulsive_force) {
                        repulsive_force *= sqrt(this->node_weights[map_it->first] * this->node_weights[map_it2->first]);
                    }

                    Displacements_2d new_disp;
                    new_disp.x_disp = (delta.x_disp / delta_mag) * repulsive_force;
                    new_disp.y_disp = (delta.y_disp / delta_mag) * repulsive_force;

                    if (this->verbose > 3) {
                        std::cout << " Node " << map_it->first << " delta_x_disp " << new_disp.x_disp
                                  << " delta_y_disp " << new_disp.y_disp << std::endl;
                    }

                    this->node_displacements[map_it->first] += new_disp;
                }

            }

        }


        //Add edge dependent attractive forces
        for (std::map<int64_t, std::map<int64_t, float>>::iterator map_it = this->edges.begin();
             map_it != this->edges.end(); ++map_it) {
            for (std::map<int64_t, float>::iterator map_it2 = map_it->second.begin();
                 map_it2 != map_it->second.end(); ++map_it2) {
                int64_t source = map_it->first, target = map_it2->first;
                float weight = map_it2->second;

                Displacements_2d delta = this->node_coords[source] - this->node_coords[target];
                float delta_mag = 0;

                delta_mag = delta.getMagnitude();

                if (delta_mag < MIN_DISTANCE) {
                    delta_mag = MIN_DISTANCE;
                    if (delta.x_disp > 0) {
                        delta.x_disp = MIN_DIFF;
                    } else {
                        delta.x_disp = -1.0 * MIN_DIFF;
                    }
                    if (delta.y_disp > 0) {
                        delta.y_disp = MIN_DIFF;
                    } else {
                        delta.y_disp = -1.0 * MIN_DIFF;
                    }
                }

                float attractive_force = delta_mag * delta_mag / this->k;
                if (this->weigh_attractive_force) {
                    attractive_force *= weight;//
                }
                Displacements_2d new_disp;
                new_disp.x_disp = (delta.x_disp / delta_mag) * attractive_force;
                new_disp.y_disp = (delta.y_disp / delta_mag) * attractive_force;
                this->node_displacements[source] -= new_disp;
                this->node_displacements[target] += new_disp;
            }
        }

        //Finalize new positions at the end of this iteration
        for (std::map<int64_t, Coords_2d>::iterator map_it = this->node_coords.begin();
             map_it != this->node_coords.end(); ++map_it) {
            //Scaling to temperature
            map_it->second += (this->node_displacements[map_it->first] *
                               ((1.0 / this->node_displacements[map_it->first].getMagnitude()) *
                                std::min(this->node_displacements[map_it->first].getMagnitude(), temperature)));
            if (this->verbose > 2) {
                Displacements_2d d = (this->node_displacements[map_it->first] *
                                      ((1.0 / this->node_displacements[map_it->first].getMagnitude()) *
                                       std::min(this->node_displacements[map_it->first].getMagnitude(), temperature)));
                std::cout << " Unclipped displacement for " << map_it->first << " " << d.x_disp << " " << d.y_disp
                          << std::endl;
            }
            //Clipping to bounds
            map_it->second.x = std::min(bounds_x_max,
                                        std::max(bounds_x_min, map_it->second.x));
            map_it->second.y = std::min(bounds_y_max,
                                        std::max(bounds_y_min, map_it->second.y));


        }

        temperature *= temperature_decay_constant;//Cooling //-= linear_cooling_factor;//*
        if (this->verbose) {
            std::cout << "Finished iteration " << i;
            if (this->verbose > 1) {
                std::cout << " Current temperature " << temperature;
            }
            std::cout << std::endl;
        }
    }

    if (this->verbose > 1) {
        for (std::map<int64_t, Coords_2d>::iterator map_it = this->node_coords.begin();
             map_it != this->node_coords.end(); ++map_it) {
            std::cout << map_it->first << " " << map_it->second.x << "," << map_it->second.y << std::endl;

        }
    }
}

void FructermanReingold::saveCSV(std::string file_path) {
    std::ofstream outfile(file_path);
    if (outfile.is_open()) {
        outfile << "node id" << ",";
        if (this->node_labels.size()) {
            outfile << "node label" << ",";
        }
        outfile << "node_x,node_y" << std::endl;

        for (std::map<int64_t, Coords_2d>::iterator map_it = this->node_coords.begin();
             map_it != this->node_coords.end(); ++map_it) {

            outfile << map_it->first << ",";

            if (this->node_labels.size()) {
                if (this->node_labels.find(map_it->first) != this->node_labels.end()) {
                    outfile << this->node_labels[map_it->first] << ",";
                } else {
                    outfile << ",";
                }
            }
            outfile << map_it->second.x << "," << map_it->second.y << std::endl;

        }

        outfile.close();
    }
}


void FructermanReingold::addNodeWeight(int64_t node_id, float node_weight) {
    this->node_weights[node_id] = node_weight;
    if (node_weight < 1.0) {
        this->weigh_repulsive_force = true;
    }
}

void FructermanReingold::setIgnoreEdgeWeights(bool value) {
    this->weigh_attractive_force = !value;
}

void FructermanReingold::setVerbose(int level) {
    this->verbose = level;
}