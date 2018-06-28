//
// Created by biot on 2018.04.06..
//

#ifndef EVOLUTIONARY_ALGORITHM_GENE_H
#define EVOLUTIONARY_ALGORITHM_GENE_H
using namespace std;
#include <vector>

class Gene
{
private:
    vector<double> feature;
public:

    vector<double> size_of_feature;
    int number_of_feature;
    vector<double> delta;
    int number_of_gene;

    Gene();

    Gene(const vector<double>& size_of_feature, const int& number_of_feature, const vector<double>& delta, const int& number_of_gene);

    ~Gene();

    Gene& operator=(const Gene& gene);

    double geti(const int& i) const;

    void seti(const int& i, const double& feature);
};
#endif //EVOLUTIONARY_ALGORITHM_GENE_H
