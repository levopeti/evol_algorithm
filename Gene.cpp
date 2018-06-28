//
// Created by biot on 2018.04.06..
//

#include "Gene.h"
#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>

using namespace std;

#define PI 3.14159

Gene::Gene()
{}

Gene::Gene(const vector<double>& size_of_feature, const int& number_of_feature, const vector<double>& delta, const int& number_of_gene)
{
    this->size_of_feature = size_of_feature;
    this->number_of_feature = number_of_feature;
    this->delta = delta;
    this->number_of_gene = number_of_gene;

    for (int i = 0; i < this->number_of_feature; i++)
    {
        //cout << "gene1" << endl;
        feature.push_back(((double)rand() / RAND_MAX) * this->size_of_feature[i] + this->delta[i]);
        //cout << "feature: "<< number_of_feature << endl;
    }
}

Gene::~Gene()
{}

Gene& Gene::operator=(const Gene& gene)
{
    size_of_feature = gene.size_of_feature;
    number_of_feature = gene.number_of_feature;
    delta = gene.delta;
    number_of_gene = gene.number_of_gene;

    if (feature.size() == 0)
    {
        for (int i = 0; i<number_of_feature; ++i)
        {
            double new_feature;
            new_feature = gene.geti(i);
            feature.push_back(new_feature);
        }
    }
    else
    {
        for (int i = 0; i<number_of_feature; ++i)
        {
            //cout << "gene2" << endl;
            feature[i] = gene.geti(i);
            //cout << "gene3" << endl;

        }
    }

    return *this;
}

//get the i.th feature
double Gene::geti(const int& i) const
{
    return feature[i];
}

//set the i.th feature
void Gene::seti(const int& i, const double& feature)
{
    if (feature - this->delta[i] > this->size_of_feature[i] )
    {
        this->feature[i] = this->size_of_feature[i] + this->delta[i];
    }
    else if(feature - this->delta[i] < 0)
    {
        this->feature[i] = this->delta[i];
    }
    else
    {
        this->feature[i] = feature;
    }

}