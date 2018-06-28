//
// Created by biot on 2018.04.06..
//

#include "Chromosome.h"
#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <fstream>
#include <time.h>
#include "Gene.h"

#include <string>


using namespace std;

#define PI 3.14159

Chromosome::Chromosome()
{}

Chromosome::Chromosome(const vector<double>& size_of_feature, const int& number_of_feature, const vector<double>& delta, const int& number_of_gene)
{
    this->size_of_feature = size_of_feature;
    this->number_of_feature = number_of_feature;
    this->delta = delta;
    this->number_of_gene = number_of_gene;
    this->fitness_value = 0;


    for (int i = 0; i < this->number_of_gene; ++i)
    {
        //cout << "chrom1" << endl;
        gene.push_back(Gene(this->size_of_feature, this->number_of_feature, this->delta, this->number_of_gene));
        //cout << "chrom12" << endl;
    }
}

Chromosome::~Chromosome()
{}

Chromosome& Chromosome::operator=(const Chromosome& chrom)
{
    this->size_of_feature = chrom.size_of_feature;
    this->number_of_feature = chrom.number_of_feature;
    this->delta = chrom.delta;
    this->number_of_gene = chrom.number_of_gene;
    this->fitness_value = chrom.fitness_value;
    //cout << "chr123" << endl;
    //cout << "number_of_gene: " << number_of_gene << endl;

    if (gene.size() == 0)
    {
        for (int i = 0; i<number_of_gene; ++i)
        {
            //cout << "chrom2" << endl;
            gene.push_back(chrom.geti(i));
            //cout << "chrom22" << endl;
        }
    }
    else
    {
        for (int i = 0; i<number_of_gene; ++i)
        {
            Gene new_gene(size_of_feature, number_of_feature, delta, number_of_gene);
            new_gene = chrom.geti(i);
            gene[i] = new_gene;
            //cout << "chrom32" << endl;
        }
    }

    return *this;
}

//get the i.th gene
const Gene& Chromosome::geti(const int& i) const
{
    return gene[i];
}

//set the i.th gene
void Chromosome::seti(const int& i, const Gene& gene)
{
    //cout << "chrom4" << endl;
    this->gene[i] = gene;
    //cout << "chrom42" << endl;
}

//get the j.th feature of i.th gene
const double Chromosome::getij_feature(const int& i, const int& j) const
{
    return gene[i].geti(j);
}

//set the j.th feature of i.th gene
void Chromosome::setij_feature(const int& i, const int& j, const double& feature)
{
    this->gene[i].seti(j, feature);
}

//get the fitness of the i.th gene
const double Chromosome::get_fit() const
{
    return this->fitness_value;
}

//set the fitness of the i.th gene
void Chromosome::set_fit(double fit)
{
    if(fit >= 0)
        this->fitness_value = fit;
}

//Fitness function
double Chromosome::fit_func(int& number_of_fit_func_call) const
{
    double fit = 0;


//    //Rastrigin function
//    int A = 10;
//
//    fit = A * number_of_gene;
//
//    for (int i = 0; i < number_of_gene; i++)
//    {
//        fit += (gene[i].geti(0) * gene[i].geti(0)) - A * cos(2 * PI * gene[i].geti(0));
//    }


//    //variance of the distance from the centre of the points
//    vector<double> avg(number_of_feature, 0);
//    double avg_dist_from_centre = 0;
//
//    //average of the features
//    for (int i = 0; i<number_of_feature; ++i)
//    {
//        for (int j = 0; j < number_of_gene; ++j)
//        {
//            avg[i] += gene[j].geti(i);
//        }
//    }
//
//    for (int i = 0; i < number_of_feature; ++i)
//    {
//        avg[i] = avg[i] / number_of_gene;
//    }
//
//    //average distance from the centrum
//    for (int i = 0; i<number_of_gene; ++i)
//    {
//        double dist = 0;
//
//        for (int j = 0; j < number_of_feature; ++j)
//        {
//            dist += (avg[j] - gene[i].geti(j)) * (avg[j] - gene[i].geti(j));
//        }
//        avg_dist_from_centre += sqrt(dist);
//    }
//    avg_dist_from_centre = avg_dist_from_centre / number_of_gene;
//
//    //summary of the distence
//    for (int i = 0; i<number_of_gene; ++i)
//    {
//        double dist = 0;
//
//        for (int j = 0; j < number_of_feature; ++j)
//        {
//            dist += (avg[j] - gene[i].geti(j)) * (avg[j] - gene[i].geti(j));
//        }
//        fit += abs(sqrt(dist) - avg_dist_from_centre);
//    }


    //ofstream file("/home/biot/projects/szakdolgozat/LOG/parameters.csv", std::ios::app);
    ofstream file("/root/inspiron/szakdoga/LOG/parameters.csv", std::ios::app);

    for (int i = 0; i<number_of_gene; ++i)
    {
        file << gene[i].geti(0) << ',';
    }
    file << endl;
    file.close();

    //system("/usr/bin/python /home/biot/projects/szakdolgozat/Python_for_EA/CNN/cnn.py");
    system("/usr/bin/python3 /root/inspiron/szakdoga/Python_for_EA/CNN/cnn.py");

    ifstream infile;
    //infile.open("/home/biot/projects/szakdolgozat/LOG/fitness.csv");
    infile.open("/root/inspiron/szakdoga/LOG/fitness.csv");
    string lastLine;

    if(infile.is_open())
    {
        getline(infile, lastLine);
    }
    infile.close();

    //cout << "The last line : " << lastLine << '\n';
    //cout << "The last line length : " << lastLine.length() << '\n';
    fit = std::stod(lastLine);
    //cout << "fit : " << fit << '\n';
    fit = 110 - (fit * 100);

    number_of_fit_func_call++;
    return fit;
}

//GA Mutation
void Chromosome::GAmutacio(const double& p)
{
    for (int i = 0; i<number_of_gene; ++i)
    {
        //cout << "chrom5" << endl;
        if (p>(rand() % 10000) / 10000)
            gene[i] = Gene(size_of_feature, number_of_feature, delta, number_of_gene);
        //cout << "chrom52" << endl;
    }
}

//BEA Mutation
void Chromosome::BEAmutacio(const int& number_of_clones, int& number_of_fit_func_call)
{
    vector<Chromosome> clones;
    int random = (int)rand() % (number_of_gene);

    //which feature has changed
    vector<int> random_place(number_of_gene, 0);

    for (int i = 0; i < number_of_clones; ++i)
    {
        clones.push_back(*this);
    }

    for (int i = 0; i < number_of_gene; ++i)
    {
        if (random_place[random] == 0)
        {
            for (int k = 1; k < number_of_clones; ++k)
            {
                clones[k].seti(i, Gene(size_of_feature, number_of_feature, delta, number_of_gene));
            }

            for (int k = 1; k < number_of_clones; ++k)
            {
                double fitbest = clones[0].get_fit();
                double fit = clones[k].fit_func(number_of_fit_func_call);
                clones[k].set_fit(fit);

                if (fit < fitbest)
                {
                    clones[0] = clones[k];
                }
            }

            for (int k = 1; k < number_of_clones; ++k)
            {
                clones[k] = clones[0];
            }

            random_place[random] = 1;
            random = (int)rand() % (number_of_gene);
        }
        else
        {
            //if random_place[random] == 1  is true, increase the random by one (in range)
            random = (random + 1) % number_of_gene;
            i--;
        }
    }

    *this = clones[0];
}

//DE Mutation
void Chromosome::DEmutation(const Chromosome& rand_chrom_1, const Chromosome& rand_chrom_2, const Chromosome& rand_chrom_3)
{

    const double F = 0.2;

    for (int i = 0; i < number_of_gene; ++i)
    {
        for (int j = 0; j < number_of_feature; ++j)
        {
            double feature = rand_chrom_1.getij_feature(i, j) + F * (rand_chrom_2.getij_feature(i, j) - rand_chrom_3.getij_feature(i, j));

            this->setij_feature(i, j, feature);
        }
    }
}

//Gene transmission
void Chromosome::gen_trans(const Chromosome& chrom)
{
    int random = ((int)rand() % number_of_gene);
    gene[random] = chrom.geti(random);
}

//Single crossover
Chromosome Chromosome::operator+(const Chromosome& chrom) const
{
    Chromosome new_chrom(this->size_of_feature, this->number_of_feature, this->delta, this->number_of_gene);
    int point = (int)rand() % number_of_gene;

    for (int i = 0; i<point; ++i)
    {
        new_chrom.seti(i, chrom.geti(i));
    }
    return new_chrom;
}

//DE Recombination
const Chromosome Chromosome::DErecombination(const Chromosome& donor_chrom, const double& CR) const
{
    int I_rand = ((int)rand() % number_of_gene * number_of_feature);
    Chromosome trial_chrom(this->size_of_feature, this->number_of_feature, this->delta, this->number_of_gene);

    for (int i = 0; i < number_of_gene; ++i)
    {
        for (int j = 0; j < number_of_feature; ++j)
        {
            double random = (double)rand() / RAND_MAX;

            if (random <= CR || ((i * number_of_gene) + j) == I_rand)
            {
                double feature = donor_chrom.getij_feature(i, j);
                trial_chrom.setij_feature(i, j, feature);
            }
            else
            {
                trial_chrom.setij_feature(i, j, this->getij_feature(i, j));
            }
        }
    }

    return trial_chrom;
}

//Gradient descent search
void Chromosome::gradient_descent(double& step_size, const int& number_of_steps, int& number_of_fit_func_call)
{
    //which feature has changed
    vector<int> random_place(number_of_gene * number_of_feature, 0);
    //delta is a step size for the derivation
    double delta;
    //fit values of the directions
    double positive_fitvalue;
    double negative_fitvalue;
    double feature;
    double original_fitvalue;
    double gradient;
    int random = (int)rand() % (number_of_gene * number_of_feature);

    for (int i = 0; i < number_of_gene; ++i)
    {
        for (int j = 0; j < number_of_feature; ++j)
        {
            if (random_place[random] == 0)
            {
                //delta = size_of_feature[random % number_of_feature] / 1000;
                delta = step_size;

                int index_gene = random / number_of_feature;
                int index_feature = random % number_of_feature;

                for (int k = 0; k < number_of_steps; ++k)
                {
                    feature = this->getij_feature(index_gene, index_feature);
                    original_fitvalue = this->get_fit();

                    this->setij_feature(index_gene, index_feature, feature + delta);
                    positive_fitvalue = this->fit_func(number_of_fit_func_call);

                    this->setij_feature(index_gene, index_feature, feature - delta);
                    negative_fitvalue = this->fit_func(number_of_fit_func_call);

                    if (positive_fitvalue <= negative_fitvalue && positive_fitvalue < original_fitvalue)
                    {
                        gradient = original_fitvalue - positive_fitvalue;
                        this->setij_feature(index_gene, index_feature, feature + delta);
                        //this->setij_feature(index_gene, index_feature, feature + (gradient * step_size));
                        this->set_fit(positive_fitvalue);
                    }
                    else if(positive_fitvalue > negative_fitvalue && negative_fitvalue < original_fitvalue)
                    {
                        gradient = original_fitvalue - negative_fitvalue;
                        this->setij_feature(index_gene, index_feature, feature - delta);
                        //this->setij_feature(index_gene, index_feature, feature - (gradient * step_size));
                        this->set_fit(negative_fitvalue);
                    }
                    else
                    {
                        this->setij_feature(index_gene, index_feature, feature);
                    }
                }

                random_place[random] = 1;
                random = (int)rand() % (number_of_gene * number_of_feature);
            }
            else
            {
                //if random_place[random] == 1  is true, increase the random by one (in range)
                random = (random + 1) % (number_of_gene * number_of_feature);
                j--;
            }
        }
    }
}

//write the chromosome into a file with file_name
void Chromosome::write_chrom(string file_name) const
{
    ofstream file(file_name.c_str(), ios::out);

    file << 'x'<< '\t' << 'y' << '\t' << 'z' << '\t' << endl;

    for (int i = 0; i<number_of_gene; ++i)
    {
        for (int j = 0; j<number_of_feature; ++j)
        {
            file << gene[i].geti(j) << '\t';
        }
        file << endl;
    }

    file << endl;
    file.close();
}