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

using namespace std;

#define PI 3.14159

Population::Population()
{}

Population::Population(const vector<double>& size_of_feature, const int& number_of_feature, const vector<double>& delta, const int& number_of_gene, const int& number_of_chrom, const int& number_of_clones, const int& mode, const double& pmut)
{
    this->size_of_feature = size_of_feature;
    this->number_of_feature = number_of_feature;
    this->delta = delta;
    this->number_of_gene = number_of_gene;
    this->number_of_chrom = number_of_chrom;
    this->number_of_clones = number_of_clones;
    this->mode = mode;
    this->pmut = pmut;

    for (int i = 0; i < this->number_of_chrom; ++i)
    {
        chromosome.push_back(Chromosome(this->size_of_feature, this->number_of_feature, this->delta, this->number_of_gene));
    }
}

Population::~Population()
{}

Population& Population::operator=(const Population& pop)
{
    size_of_feature = pop.size_of_feature;
    number_of_feature = pop.number_of_feature;
    delta = pop.delta;
    number_of_gene = pop.number_of_gene;
    number_of_chrom = pop.number_of_chrom;
    number_of_clones = pop.number_of_clones;
    mode = pop.mode;
    pmut = pop.pmut;

    if (chromosome.size() == 0)
    {
        for (int i = 0; i<number_of_feature; ++i)
        {
            chromosome.push_back(pop.geti(i));
        }
    }
    else
    {
        for (int i = 0; i<number_of_chrom; ++i)
        {
            chromosome[i] = pop.geti(i);
        }
    }

    return *this;
}

//get the i.th chromosome
const Chromosome& Population::geti(const int& i) const
{
    return chromosome[i];
}

//set the i.th chromosome
void Population::seti(const int& i, const Chromosome& chrom)
{
    this->chromosome[i] = chrom;
}

//create a new population
Population Population::new_pop(const vector<double>& size_of_feature, const int& number_of_feature, const vector<double>& delta, const int& number_of_gene, const int& number_of_chrom, const int& number_of_clones, const int& mode, const double& pmut, const double& CR, double& step_size, const int& number_of_steps, bool& Lamarck, int& number_of_fit_func_call)
{
    Population new_pop(size_of_feature, number_of_feature, delta, number_of_gene, number_of_chrom, number_of_clones, mode, pmut);

    //GA
    if (mode == 0)
    {
        //save the best part of population
        for (int j = 0; j<number_of_chrom / 4; ++j)
        {
            new_pop.seti(j, chromosome[j]);
        }

        //Crossover
        for (int j = number_of_chrom / 4; j < number_of_chrom * 3 / 4; ++j)
        {
            double fit_value;
            int random1 = (int)rand() % (number_of_chrom / 4);
            int random2 = (int)rand() % (number_of_chrom / 4);

            Chromosome new_chrom(size_of_feature, number_of_feature, delta, number_of_gene);
            new_chrom = new_pop.geti(random1) + new_pop.geti(random2);
            fit_value = new_chrom.fit_func(number_of_fit_func_call);
            new_chrom.set_fit(fit_value);
            new_pop.seti(j, new_chrom);
        }

        new_pop.sorting(number_of_fit_func_call);

        //GA Mutation
        for (int i = 1; i<number_of_chrom; ++i)
        {
            double fit_value;
            Chromosome new_chrom(size_of_feature, number_of_feature, delta, number_of_gene);
            new_chrom = new_pop.geti(i);
            new_chrom.GAmutacio(pmut);
            fit_value = new_chrom.fit_func(number_of_fit_func_call);
            new_chrom.set_fit(fit_value);
            new_pop.seti(i, new_chrom);
        }
    }

    //BEA
    if (mode == 1)
    {
        //save the best part of population
        for (int j = 0; j<number_of_chrom / 4; ++j)
        {
            new_pop.seti(j, chromosome[j]);
        }

        //Gene Transmission
        for (int j = number_of_chrom / 4; j < number_of_chrom * 3 / 4; ++j)
        {
            int random = (int)rand() % (number_of_chrom / 4);
            Chromosome new_chrom(size_of_feature, number_of_feature, delta, number_of_gene);
            new_chrom = new_pop.geti(random);
            new_chrom.gen_trans(new_pop.geti(random));
            new_pop.seti(j, new_chrom);
        }

        new_pop.sorting(number_of_fit_func_call);

        //BEA Mutation
        for (int i = 0; i<number_of_chrom; ++i)
        {
            double fit_value;
            Chromosome new_chrom(size_of_feature, number_of_feature, delta, number_of_gene);
            new_chrom = new_pop.geti(i);
            new_chrom.BEAmutacio(number_of_clones, number_of_fit_func_call);
//            fit_value = new_chrom.fit_func(number_of_fit_func_call);
//            new_chrom.set_fit(fit_value);
            new_pop.seti(i, new_chrom);
        }
    }

    //DE
    if (mode == 2)
    {
        for (int i = 0; i < number_of_chrom; ++i)
        {
            double fit_value;

            //DE Mutation
            int random1 = (int)rand() % (number_of_chrom);
            int random2 = (int)rand() % (number_of_chrom);
            int random3 = (int)rand() % (number_of_chrom);

            Chromosome donor_chrom(size_of_feature, number_of_feature, delta, number_of_gene);
            donor_chrom.DEmutation(chromosome[random1], chromosome[random2], chromosome[random3]);

            //Recombination
            Chromosome trial_chrom(size_of_feature, number_of_feature, delta, number_of_gene);
            trial_chrom = chromosome[i].DErecombination(donor_chrom, CR);
            fit_value = trial_chrom.fit_func(number_of_fit_func_call);
            trial_chrom.set_fit(fit_value);

            //Selection
            if (chromosome[i].get_fit()>=trial_chrom.get_fit())
            {
                new_pop.seti(i, trial_chrom);
            }
            else
            {
                new_pop.seti(i, chromosome[i]);
            }
        }
    }

    new_pop.sorting(number_of_fit_func_call);

    //Gradient descent
    if (Lamarck)
    {
        for (int i = 0; i < number_of_chrom; ++i)
        {
            Chromosome new_chrom(size_of_feature, number_of_feature, delta, number_of_gene);

            if (mode == -1)
            {
                new_chrom = this->geti(i);
            }
            else
            {
                new_chrom = new_pop.geti(i);
            }

            new_chrom.gradient_descent(step_size, number_of_steps, number_of_fit_func_call);

//            if (mode == -1)
//            {
//                double fit_value;
//                fit_value = new_chrom.fit_func(number_of_fit_func_call);
//                new_chrom.set_fit(fit_value);
//            }

            new_pop.seti(i, new_chrom);
        }
    }

    new_pop.sorting(number_of_fit_func_call);
    return new_pop;
}

//sorting of the chromosomes by fit values in the old population
void Population::sorting(int& number_of_fit_func_call)
{
    vector<Chromosome> sorted;
    vector<double> fit_array;
    int best_fit_index = 0;
    double fit_best;
    int fit = 0;
    double value;

    //fill the fitness array
    for (int i = 0; i<number_of_chrom; ++i)
    {
        value = chromosome[i].get_fit();

        if (value < 0.1)
        {
            value = chromosome[i].fit_func(number_of_fit_func_call);
        }

        chromosome[i].set_fit(value);
        fit_array.push_back(value);
    }

    //sorting of the chromosomes
    for (int j = 0; j<number_of_chrom; ++j)
    {
        fit_best = FLT_MAX;

        for (int i = 0; i<number_of_chrom; ++i)
        {
            if (fit_array[i]<fit_best && fit_array[i]>=0)
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
void Population::save_fit(string& fit_ness_fn)
{
    ofstream file(fit_ness_fn.c_str(), ios::out | std::ios::app);
    file << chromosome[0].get_fit() << endl;
    file.close();
}

////sorting of the chromosomes by fit values in the new population
//void Population::sorting_new()
//{
//	vector<Chromosome> sorted;
//	vector<double> fit_array;
//	int best_fit_index = 0;
//	double fit_best;
//	int fit = 0;
//	double value;
//
//	//fill the fitness array
//	for (int i = 0; i<number_of_chrom; ++i)
//	{
//		value = chromosome[i].fit_func();
//		chromosome[i].set_fit(value);
//		fit_array.push_back(value);
//	}
//
//	//sorting of the chromosomes
//	for (int j = 0; j<number_of_chrom; ++j)
//	{
//		fit_best = FLT_MAX;
//
//		for (int i = 0; i<number_of_chrom; ++i)
//		{
//			if (fit_array[i]<fit_best && fit_array[i] >= 0)
//			{
//				fit = i;
//				fit_best = fit_array[i];
//			}
//		}
//		fit_array[fit] = -1;
//		sorted.push_back(chromosome[fit]);
//	}
//
//	chromosome = sorted;
//
//	//print out the best fitness value in a file
//	ofstream file(fit_ness_fn.c_str(), ios::out | std::ios::app);
//	file << chromosome[0].get_fit() << endl;
//	file.close();
//}

////index of chromosome with the best fitness
//int Population::best(double& best_fit_value)
//{
//	vector<double> fit_array;
//	int best_fit_index = 0;
//	double fit_best;
//
//	//fill the fitness array
//	for (int i = 0; i<this->number_of_chrom; ++i)
//	{
//		fit_array.push_back(chromosome[i].fit_func());
//	}
//
//	fit_best = FLT_MAX;
//	//cout << "fit_best: " << fit_best << endl;
//
//	for (int i = 0; i<this->number_of_chrom; ++i)
//	{
//		if (fit_array[i]<fit_best)
//		{
//			best_fit_index = i;
//			fit_best = fit_array[i];
//		}
//	}
//	best_fit_value = fit_best;
//
//	//cout << "best_fit_index: " << best_fit_index << endl;
//	//cout << "fit_best: " << fit_best << endl << endl;
//
//	return best_fit_index;
//}