//
// Created by biot on 2018.04.06..
//

#include "Population.h"
#include <time.h>
#include "Gene.h"
#include "Chromosome.h"
#include <cfloat>
#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <fstream>
#include "Iteration.h"
#include "Swarm.h"

using namespace std;

#define PI 3.14159

Iteration::Iteration()
{}

Iteration::Iteration(const vector<double>& size_of_feature, const int& number_of_feature, const vector<double>& delta, const int& number_of_gene, const int& number_of_chrom, const vector<double>& C_SSO, const int& mode)
{
    this->size_of_feature = size_of_feature;
    this->number_of_feature = number_of_feature;
    this->delta = delta;
    this->number_of_gene = number_of_gene;
    this->number_of_chrom = number_of_chrom;
    this->mode = mode;
    this->C_SSO = C_SSO;

    for (int i = 0; i < this->number_of_chrom; ++i)
    {
        chromosome.push_back(Swarm_Chromosome(this->size_of_feature, this->number_of_feature, this->delta, this->number_of_gene));
    }
}

Iteration::~Iteration()
{
}

//get the i.th chromosome
const Swarm_Chromosome& Iteration::geti(const int& i) const
{
    return chromosome[i];
}

//set the i.th chromosome
void Iteration::seti(const int& i, const Swarm_Chromosome& chrom)
{
    this->chromosome[i] = chrom;
}

void Iteration::update(const double& inertia, const double& phi_p, const double& phi_g, const int& max_iteration, int& iteration, const int& max_number_of_seeds, const int& min_number_of_seeds, const double& exponent, const double& sigma_initial, const double& sigma_final, Swarm_Chromosome& global_best_chrom, double& step_size, const int& number_of_steps, bool& Lamarck, int& number_of_fit_func_call)
{
    //PSO
    if (mode == 0)
    {
        for (int j = 0; j < number_of_chrom; ++j)
        {
            chromosome[j].update_PSO(inertia, phi_p, phi_g, global_best_chrom, number_of_fit_func_call);
        }
    }

    //SSO
    if (mode == 1)
    {
        for (int j = 0; j < number_of_chrom; ++j)
        {
            chromosome[j].update_SSO(inertia, phi_p, phi_g, global_best_chrom, C_SSO, number_of_fit_func_call);
        }
    }

    //IWO
    if (mode == 2)
    {
        this->sorting(number_of_fit_func_call);

        const double fit_min = chromosome[0].get_fit();
        const double fit_max = chromosome[number_of_chrom - 1].get_fit();

        for (int j = 0; j < number_of_chrom; ++j)
        {
            vector<Swarm_Chromosome> offsprings;
            offsprings = chromosome[j].update_IWO(max_iteration, iteration, max_number_of_seeds, min_number_of_seeds, exponent, sigma_initial, sigma_final, fit_max, fit_min, number_of_fit_func_call);
            chromosome.insert(chromosome.end(), offsprings.begin(), offsprings.end());
        }

        this->sorting(number_of_fit_func_call);

        while (chromosome.size() != number_of_chrom)
        {
            chromosome.pop_back();
        }
    }

    //Gradient descent
    if (Lamarck)
    {
        for (int i = 0; i < number_of_chrom; ++i)
        {
            chromosome[i].gradient_descent(step_size, number_of_steps, number_of_fit_func_call);
            chromosome[i].set_personal_best();
        }
        //cout << "feature: " << chromosome[0].getij_feature(0, 0) << endl;
    }

    this->sorting(number_of_fit_func_call);
}

void Iteration::update_global_best(double& global_best, Swarm_Chromosome& global_best_chrom, int& counter_till_end)
{
    counter_till_end++;

    for (int j = 0; j<number_of_chrom; ++j)
    {
        if (chromosome[j].get_personal_best() < global_best)
        {
            global_best = chromosome[j].get_personal_best();
            global_best_chrom = chromosome[j];
            counter_till_end = 0;
        }
    }
}

//sorting of the chromosomes by fit values in the old population
void Iteration::sorting(int& number_of_fit_func_call)
{
    vector<Swarm_Chromosome> sorted;
    vector<double> fit_array;
    int best_fit_index = 0;
    double fit_best;
    int fit = 0;
    double value;

    //fill the fitness array
    for (int i = 0; i<chromosome.size(); ++i)
    {
        value = chromosome[i].get_fit();

        if (value == 0)
        {
            value = chromosome[i].fit_func(number_of_fit_func_call);
        }

        chromosome[i].set_fit(value);
        chromosome[i].set_personal_best();
        fit_array.push_back(value);
    }

    //sorting of the chromosomes
    for (int j = 0; j<chromosome.size(); ++j)
    {
        fit_best = FLT_MAX;

        for (int i = 0; i<chromosome.size(); ++i)
        {
            if (fit_array[i]<fit_best && fit_array[i] >= 0)
            {
                fit = i;
                fit_best = fit_array[i];
            }
        }

        fit_array[fit] = -1;
        sorted.push_back(chromosome[fit]);
    }

    chromosome = sorted;
}

//print out the best fitness value in a file
void Iteration::save_fit(string& fit_ness_fn)
{
    ofstream file(fit_ness_fn.c_str(), ios::out | std::ios::app);
    file << chromosome[0].get_fit() << endl;
    file.close();
}