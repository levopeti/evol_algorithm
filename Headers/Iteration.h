//
// Created by biot on 2018.04.06..
//

#ifndef EVOLUTIONARY_ALGORITHM_ITERATION_H
#define EVOLUTIONARY_ALGORITHM_ITERATION_H
//class Chromosome;
class Swarm_Chromosome;
using namespace std;
#include <vector>

class Iteration
{
private:
    vector<Swarm_Chromosome> chromosome;
public:
    vector<double> size_of_feature;
    int number_of_feature;
    vector<double> delta;
    int number_of_gene;
    int number_of_chrom;
    int mode;
    vector<double> C_SSO;

public:
    Iteration();

    Iteration(const vector<double>& size_of_feature, const int& number_of_feature, const vector<double>& delta, const int& number_of_gene, const int& number_of_chrom, const vector<double>& C_SSO, const int& mode);

    ~Iteration();

    const Swarm_Chromosome& geti(const int& i) const;

    void seti(const int& i, const Swarm_Chromosome& chrom);

    void update(const double& inertia, const double& phi_p, const double& phi_g, const int& max_iteration, int& iteration, const int& max_number_of_seeds, const int& min_number_of_seeds, const double& exponent, const double& sigma_initial, const double& sigma_final, Swarm_Chromosome& global_best_chrom, double& step_size, const int& number_of_steps, bool& Lamarck, int& number_of_fit_func_call);

    void update_global_best(double& global_best, Swarm_Chromosome& global_best_chrom, int& counter_till_end);

    void sorting(int& number_of_fit_func_call);

    void save_fit(string& fit_ness_fn);
};
#endif //EVOLUTIONARY_ALGORITHM_ITERATION_H
