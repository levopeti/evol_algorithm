//
// Created by biot on 2018.04.06..
//

#include "Swarm.h"
#include "Chromosome.h"
#include <iostream>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <fstream>
#include <time.h>
#include "Gene.h"
#include <numeric>

using namespace std;

#define PI 3.14159

Swarm_Chromosome::Swarm_Chromosome()
{}

Swarm_Chromosome::Swarm_Chromosome(const vector<double>& size_of_feature, const int& number_of_feature, const vector<double>& delta, const int& number_of_gene) : Chromosome(size_of_feature, number_of_feature, delta, number_of_gene)
{
    personal_best = FLT_MAX;

    for (int i = 0; i < number_of_gene; ++i)
    {
        this->personal_best_chrom.push_back(this->geti(i));
    }

    for (int i = 0; i < number_of_gene; ++i)
    {
        for (int j = 0; j < number_of_feature; ++j)
        {
            double value_0 = .0;
            double value = ((((double)rand() / RAND_MAX) * 2 * this->size_of_feature[j]) - this->size_of_feature[j]);
            this->velocity.push_back(value_0);
        }
    }

    /*double sum = 0.;

    for (int i = 0; i < number_of_feature * number_of_gene; ++i)
    {
        sum += velocity[i] * velocity[i];
    }

    double norm = sqrt(sum);

    for (int i = 0; i < number_of_gene; ++i)
    {
        for (int j = 0; j < number_of_feature; ++j)
        {
            int index = i * number_of_feature + j;
            velocity[index] *= 5 / norm;
        }
    }*/

}

Swarm_Chromosome::~Swarm_Chromosome()
{}

Swarm_Chromosome& Swarm_Chromosome::operator=(const Swarm_Chromosome& chrom)
{
    this->size_of_feature = chrom.size_of_feature;
    this->number_of_feature = chrom.number_of_feature;
    this->delta = chrom.delta;
    this->number_of_gene = chrom.number_of_gene;
    this->personal_best = chrom.personal_best;

    double fit_value;
    fit_value = chrom.get_fit();
    this->set_fit(fit_value);


    for (int i = 0; i<number_of_gene; ++i)
    {
        for (int j = 0; j < number_of_feature; ++j)
        {
            double value;
            value = chrom.getij_feature(i, j);
            this->setij_feature(i, j, value);
        }
    }
    return *this;
}

void Swarm_Chromosome::set_personal_best()
{
    if (this->get_fit() < personal_best)
    {
        personal_best = this->get_fit();

        for (int i = 0; i < number_of_gene; ++i)
        {
            for (int j = 0; j < number_of_feature; ++j)
            {
                double value = this->getij_feature(i, j);
                personal_best_chrom[i].seti(j, value);
            }
        }
    }
}

double Swarm_Chromosome::get_personal_best()
{
    return this->personal_best;
}

void Swarm_Chromosome::update_PSO(const double& inertia, const double& phi_p, const double& phi_g, Swarm_Chromosome& global_best_chrom, int& number_of_fit_func_call)
{
    for (int i = 0; i<number_of_gene; ++i)
    {
        for (int j = 0; j < number_of_feature; ++j)
        {
            int index = i * number_of_feature + j;
            double feature = this->getij_feature(i, j);
            double velocity = this->velocity[index];
            this->setij_feature(i, j, feature + velocity);
        }
    }
    double sum = 0.;

    /*for (int i = 0; i < number_of_feature * number_of_gene; ++i)
    {
        sum += velocity[i] * velocity[i];
    }

    double norm = sqrt(sum);*/


    this->set_fit(this->fit_func(number_of_fit_func_call));
    this->set_personal_best();
    this->update_velocity(inertia, phi_p, phi_g, global_best_chrom);
    //int sum = accumulate(velocity.begin(), velocity.end(), 0);
    //cout << "sum: " << sum << endl;
}

void Swarm_Chromosome::update_SSO(const double& inertia, const double& phi_p, const double& phi_g, Swarm_Chromosome& global_best_chrom, const vector<double>& C_SSO, int& number_of_fit_func_call)
{
    double feature;

    for (int i = 0; i<number_of_gene; ++i)
    {
        for (int j = 0; j < number_of_feature; ++j)
        {
            double random = (double)rand() / RAND_MAX;

            // previous value
            if (random <= C_SSO[0])
            {
                feature = this->getij_feature(i, j);
                this->setij_feature(i, j, feature);
            }

            // personal best value
            if (random <= C_SSO[1] && random > C_SSO[0])
            {
                feature = personal_best_chrom[i].geti(j);
                this->setij_feature(i, j, feature);
            }

            // global best value
            if (random <= C_SSO[2] && random > C_SSO[1])
            {
                feature = global_best_chrom.getij_feature(i, j);
                this->setij_feature(i, j, feature);
            }

            // random value
            if (random > C_SSO[2])
            {
                feature = ((double)rand() / RAND_MAX) * this->size_of_feature[j] + this->delta[j];
                this->setij_feature(i, j, feature);
            }
        }
    }

    this->set_fit(this->fit_func(number_of_fit_func_call));
    this->set_personal_best();
}

vector<Swarm_Chromosome> Swarm_Chromosome::update_IWO(const int& max_iteration, int& iteration, const int& max_number_of_seeds, const int& min_number_of_seeds, const double& exponent, const double& sigma_initial, const double& sigma_final, const double& fit_max, const double& fit_min, int& number_of_fit_func_call)
{
    double sigma;
    vector<Swarm_Chromosome> offsprings;

    if (max_iteration - iteration > 0)
    {
        sigma = (pow((((double)max_iteration - (double)iteration) / (double)max_iteration), exponent)) * (sigma_initial - sigma_final) + sigma_final;
    }
    else
    {
        sigma = sigma_final;
    }

    double ratio = (fit_max - this->get_fit()) / (fit_max - fit_min);
    int number_of_seeds = (int)(min_number_of_seeds + (max_number_of_seeds - min_number_of_seeds) * ratio);

    for (int k = 0; k < number_of_seeds; k++)
    {
        offsprings.push_back(*this);

        for (int i = 0; i < number_of_gene; ++i)
        {
            for (int j = 0; j < number_of_feature; ++j)
            {
                double random = ((double)rand() / RAND_MAX) * 2 * sigma - sigma;
                double feature = this->getij_feature(i, j);
                offsprings[k].setij_feature(i, j, feature + random);
            }
        }
        offsprings[k].set_fit(offsprings[k].fit_func(number_of_fit_func_call));
        offsprings[k].set_personal_best();
    }

    return offsprings;
}

void Swarm_Chromosome::update_velocity(const double inertia, const double phi_p, const double phi_g, Swarm_Chromosome& global_best_chrom)
{
    for (int i = 0; i < number_of_gene; ++i)
    {
        for (int j = 0; j < number_of_feature; ++j)
        {
            int index = i * number_of_feature + j;
            double r_p = ((double)rand() / RAND_MAX);
            double r_g = ((double)rand() / RAND_MAX);
            double value = inertia * velocity[index] + phi_p * r_p * (personal_best_chrom[i].geti(j) - this->getij_feature(i, j)) + phi_g * r_g * (global_best_chrom.getij_feature(i, j) - this->getij_feature(i, j));
            velocity[index] = value;
        }
    }

    double sum = 0.;

    for (int i = 0; i < number_of_feature * number_of_gene; ++i)
    {
        sum += velocity[i] * velocity[i];
    }

    double norm = sqrt(sum);

    for (int i = 0; i < number_of_gene; ++i)
    {
        for (int j = 0; j < number_of_feature; ++j)
        {
            int index = i * number_of_feature + j;
            velocity[index] *= (5 / norm);
        }
    }

//    double sum1 = 0.;
//    for (int i = 0; i < number_of_feature * number_of_gene; ++i)
//    {
//        sum1 += velocity[i] * velocity[i];
//    }
//
//    double norm1 = sqrt(sum1);
}


//write the chromosome into a file with file_name
void Swarm_Chromosome::write_chrom(string file_name) const
{
    ofstream file(file_name.c_str(), ios::out);

    file << 'x' << '\t' << 'y' << '\t' << 'z' << '\t' << endl;

    for (int i = 0; i<number_of_gene; ++i)
    {
        for (int j = 0; j<number_of_feature; ++j)
        {
            file << this->getij_feature(i, j) << '\t';
        }
        file << endl;
    }

    file << endl;
    file.close();
}