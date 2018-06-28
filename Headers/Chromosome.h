//
// Created by biot on 2018.04.06..
//

#ifndef EVOLUTIONARY_ALGORITHM_CHROMOSOME_H
#define EVOLUTIONARY_ALGORITHM_CHROMOSOME_H
class Gene;
using namespace std;
#include <vector>
#include <string>

class Chromosome
{
private:
    vector<Gene> gene;
    double fitness_value;
public:
    vector<double> size_of_feature;
    int number_of_feature;
    vector<double> delta;
    int number_of_gene;

    Chromosome();

    Chromosome(const vector<double>& size_of_feature, const int& number_of_feature, const vector<double>& delta, const int& number_of_gene);

    ~Chromosome();

    Chromosome& operator=(const Chromosome& chrom);

    //get the i.th gene
    const Gene& geti(const int& i) const;

    //set the i.th gene
    void seti(const int& i, const Gene& gene);

    //get the j.th feature of i.th gene
    const double getij_feature(const int& i, const int& j) const;

    //set the j.th feature of i.th gene
    void setij_feature(const int& i, const int& j, const double& feature);

    //get the fitness of the i.th gene
    const double get_fit() const;

    //set the fitness of the i.th gene
    void set_fit(double fit);

    //Fitness function
    virtual double fit_func(int& number_of_fit_func_call) const;

    //GA Mutation
    void GAmutacio(const double& p);

    //BEA Mutation
    void BEAmutacio(const int& noumber_of_clones, int& number_of_fit_func_call);

    //DE Mutation
    void DEmutation(const Chromosome& rand_chrom_1, const Chromosome& rand_chrom_2, const Chromosome& rand_chrom_3);

    //Gene transmission
    void gen_trans(const Chromosome& chrom);

    //Single crossover
    Chromosome operator+(const Chromosome& chrom) const;

    //DE Recombination
    const Chromosome DErecombination(const Chromosome& donor_chrom, const double& CR) const;

    //Gradient descent search
    void gradient_descent(double& step_size, const int& number_of_steps, int& number_of_fit_func_call);

    //write the chromosome into a file with file_name
    void write_chrom(string file_name) const;
};
#endif //EVOLUTIONARY_ALGORITHM_CHROMOSOME_H
