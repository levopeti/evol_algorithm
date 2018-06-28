//
// Created by biot on 2018.04.06..
//

#ifndef EVOLUTIONARY_ALGORITHM_POPULATION_H
#define EVOLUTIONARY_ALGORITHM_POPULATION_H
class Chromosome;
using namespace std;
#include <vector>
#include <string>

class Population
{
private:
    vector<Chromosome> chromosome;
public:
    vector<double> size_of_feature;
    int number_of_feature;
    vector<double> delta;
    int number_of_gene;
    int number_of_chrom;
    int number_of_clones;
    int mode;
    double pmut;

    Population();

    Population(const vector<double>& size_of_feature, const int& number_of_feature, const vector<double>& delta, const int& number_of_gene, const int& number_of_chrom, const int& number_of_clones, const int& mode, const double& pmut);

    ~Population();

    Population& operator=(const Population& pop);

    const Chromosome& geti(const int& i) const;

    void seti(const int& i, const Chromosome& chrom);

    //sorting of the chromosomes by fit values in the old population
    void sorting(int& number_of_fit_func_call);

    //sorting of the chromosomes by fit values in the new population
    //void sorting_new();

    Population new_pop(const vector<double>& size_of_feature, const int& number_of_feature, const vector<double>& delta, const int& number_of_gene, const int& number_of_chrom, const int& number_of_clones, const int& mode, const double& pmut, const double& CR, double& step_size, const int& number_of_steps, bool& Lamarck, int& number_of_fit_func_call);

    //print out the best fitness value in a file
    void save_fit(string& fit_ness_fn);
};
#endif //EVOLUTIONARY_ALGORITHM_POPULATION_H
