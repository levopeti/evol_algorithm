//
// Created by biot on 2018.04.06..
//

#ifndef EVOLUTIONARY_ALGORITHM_SWARM_H
#define EVOLUTIONARY_ALGORITHM_SWARM_H
#include <vector>
#include "Chromosome.h"

class Gene;
class Chromosome;
using namespace std;


class Swarm_Chromosome : public Chromosome
{
private:
    double personal_best;
    vector<Gene> personal_best_chrom;
    vector<double> velocity;
public:

    Swarm_Chromosome();

    Swarm_Chromosome(const vector<double>& size_of_feature, const int& number_of_feature, const vector<double>& delta, const int& number_of_gene);

    ~Swarm_Chromosome();

    Swarm_Chromosome& operator=(const Swarm_Chromosome& chrom);

    //const double& getij_feature(const int& i, const int& j) const;
    //
    //void setij_feature(const int& i, const int& j, const double& feature);

    void set_personal_best();

    double get_personal_best();

    void update_PSO(const double& inertia, const double& phi_p, const double& phi_g, Swarm_Chromosome& global_best_chrom, int& number_of_fit_func_call);

    void update_SSO(const double& inertia, const double& phi_p, const double& phi_g, Swarm_Chromosome& global_best_chrom, const vector<double>& C_SSO, int& number_of_fit_func_call);

    vector<Swarm_Chromosome>  update_IWO(const int& max_iteration, int& iteration, const int& max_number_of_seeds, const int& min_number_of_seeds, const double& exponent, const double& sigma_initial, const double& sigma_final, const double& fit_max, const double& fit_min, int& number_of_fit_func_call);

    void update_velocity(const double inertia, const double phi_p, const double phi_g, Swarm_Chromosome& global_best_chrom);

    void write_chrom(string file_name) const;
};
#endif //EVOLUTIONARY_ALGORITHM_SWARM_H
