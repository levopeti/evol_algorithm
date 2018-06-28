// Evolutionary algorithm.cpp : Defines the entry point for the console application.
//

#undef main
#include <iostream>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <fstream>
#include "Gene.h"
#include "Chromosome.h"
#include "Population.h"
#include "Swarm.h"
#include "Iteration.h"
#include "report.h"
#include "Grid_Search.h"
#include "python_test.h"
#include <string>
#include <sys/stat.h>
//#include <Python.h>

using namespace std;

//size of area
#define A  500
#define R  5.12
#define POLYGON 0
#define CIRCLE 1
#define EVOLUTION 1
#define SWARM 1 - EVOLUTION

// GE:    0 = GA,  1 = BEA,  2 = DE
// SWARM: 0 = PSO, 1 = SSO, 2 = IWO
// -1 = only gradient descent

bool Grid_Search = false;
bool Lamarck = false;
const int mode = 0;

clock_t t_start, t_end;
double running_time;
int number_of_fit_func_call = 0;
const int patience = 5;
int LIMIT = 5;
int id_of_search = 0;
int iteration_block = 1;

ofstream report;

// Base directory name
//string base_dir_name = "/home/biot/projects/szakdolgozat/LOG/";
string base_dir_name = "/root/inspiron/szakdoga/LOG/";
string grid_base_dir_name = base_dir_name;

// Directory name for chromosomes
string chrom_dir_name;

// Fitness value file name
string fit_ness_fn;
string g_b_fit_ness_fn;
string step_size_fn;

// Report file name
string report_fn;

// Grid search report file name
string grid_report_fn;

int main(int argc, wchar_t *argv[])
{

    make_grid_base_dir_name(grid_base_dir_name);
    grid_report_fn = grid_base_dir_name + "grid_report.dat";
    mkdir(grid_base_dir_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    init_make_grid_search_report(grid_report_fn);
    //////////////////////////////////////////////////////////////////////////
    // Common parameters
    //////////////////////////////////////////////////////////////////////////

//    const int number_of_feature = 3;
//
//    //2*k|keZ  number of point
//    const int number_of_gene = 1000;
//
//    //4*k|keZ  number of instances
//    const int number_of_chrom = 1;
//
//    const double features[number_of_feature] = { A, A, A };
//    const double deltas[number_of_feature] = { 0, 0, 0 };

//    const int number_of_feature = 1;
//
//    //2*k|keZ  number of point
//    const int number_of_gene = 10;
//
//    //4*k|keZ  number of instances
//    const int number_of_chrom = 160;
//
//    const double features[number_of_feature] = {10.24};
//    const double deltas[number_of_feature] = {-5.12};

    const int number_of_feature = 1;

    //2*k|keZ  number of point
    const int number_of_gene = 11;

    //4*k|keZ  number of instances
    const int number_of_chrom = 8;

    const double features[number_of_feature] = { 100};
    const double deltas[number_of_feature] = { 0};

    vector<double> size_of_feature;
    for (int i = 0; i < number_of_feature; i++)
    {
        size_of_feature.push_back(features[i]);
    }

    vector<double> delta;
    for (int i = 0; i < number_of_feature; i++)
    {
        delta.push_back(deltas[i]);
    }


    //////////////////////////////////////////////////////////////////////////
    // Parameters of Gradient descent
    //////////////////////////////////////////////////////////////////////////


    double init_step_size = 16;
    int number_of_steps = 1;
    double decay = 0.5;

    //////////////////////////////////////////////////////////////////////////
    // Parameters of GA
    //////////////////////////////////////////////////////////////////////////

    //probability of GAmutation
    double pmut = 0.7;


    //////////////////////////////////////////////////////////////////////////
    // Parameters of BA
    //////////////////////////////////////////////////////////////////////////

    //number of clones in the bacterial mutation
    int number_of_clones = 2;


    //////////////////////////////////////////////////////////////////////////
    // Parameters of DE
    //////////////////////////////////////////////////////////////////////////

    //probability of DE recombination
    double CR = 0.3;


    ///////////////////////////////////////////////////////////////////////////
    // Parameters os PSO
    ///////////////////////////////////////////////////////////////////////////

    double inertia = 1;
    double phi_p = 2;
    double phi_g = 6;


    //////////////////////////////////////////////////////////////////////////
    // Parameters os SSO
    //////////////////////////////////////////////////////////////////////////

    vector<double> C_SSO;
    double C_w = 0.25;
    double C_p = 0.5;
    double C_g = 0.9;
    C_SSO.push_back(C_w);
    C_SSO.push_back(C_p);
    C_SSO.push_back(C_g);


    ///////////////////////////////////////////////////////////////////////////
    // Parameters os IWO
    ///////////////////////////////////////////////////////////////////////////

    int max_iteration = 10;
    int max_number_of_seeds = 7;
    int min_number_of_seeds = 1;
    double exponent = 2;
    double sigma_initial = 32;
    double sigma_final = 8;

    ///////////////////////////////////////////////////////////////////////////

    srand(time(NULL));
    bool no_grid_search = false;

    double init_pmut;
    int init_number_of_clones;
    double init_CR;
    double init_init_step_size;
    int init_number_of_steps;
    double init_decay;
    double init_inertia;
    double init_phi_p;
    double init_phi_g;
    vector<double> init_C_SSO;
    init_C_SSO.push_back(C_w);
    init_C_SSO.push_back(C_p);
    init_C_SSO.push_back(C_g);

    int init_max_iteration;
    int init_max_number_of_seeds;
    int init_min_number_of_seeds;
    double init_exponent;
    double init_sigma_initial;
    double init_sigma_final;

    if (Grid_Search)
    {
        //minimum values of parameters to grid search
        //GA
        init_pmut = pmut = 0;


        //BEA
        init_number_of_clones = number_of_clones = 1;  //minimum 1

        //DE
        init_CR = CR = 0;

        //PSO
        init_inertia = inertia = 0.1;
        init_phi_p = phi_p = 0.1;
        init_phi_g = phi_g = 0.1;


        //SSO
        init_C_SSO[0] = C_SSO[0] = 0;
        init_C_SSO[1] = C_SSO[1] = 0.1;
        init_C_SSO[2] = C_SSO[2] = 0.2;

        //IWO
        init_max_iteration = max_iteration = 5000;
        init_max_number_of_seeds = max_number_of_seeds = 1;
        init_min_number_of_seeds = min_number_of_seeds = 1;
        init_exponent = exponent = 2;
        init_sigma_initial = sigma_initial = 50;
        init_sigma_final = sigma_final = 0;

        //Lamarck
        init_init_step_size = init_step_size = 1;
        init_number_of_steps = number_of_steps = 1;
        init_decay = decay = 0.1;
    }
    else
    {
        Grid_Search = true;
        no_grid_search = true;
    }

    while(Grid_Search)
    {
        //cout << "1: " << init_step_size << endl;
        //label, hogy be ért-e a limit alá!
        //////////////////////////////////////////////////////////////////////////
        // Initialize
        //////////////////////////////////////////////////////////////////////////

        int counter_till_end = 0;
        number_of_fit_func_call = 0;
        t_start = clock();
        int generation = 0;
        int limit_label = 0;
        double global_best = FLT_MAX;
        double step_size = init_step_size;

        //GE
        Population old_pop(size_of_feature, number_of_feature, delta, number_of_gene, number_of_chrom, number_of_clones, mode, pmut);
        old_pop.sorting(number_of_fit_func_call);
        Population new_pop(size_of_feature, number_of_feature, delta, number_of_gene, number_of_chrom, number_of_clones, mode, pmut);

        //Swarm
        Swarm_Chromosome global_best_chrom(size_of_feature, number_of_feature, delta, number_of_gene);
        Iteration Iteration(size_of_feature, number_of_feature, delta, number_of_gene, number_of_chrom, C_SSO, mode);
        Iteration.sorting(number_of_fit_func_call);

        _init_(grid_report_fn,grid_base_dir_name, chrom_dir_name, fit_ness_fn, g_b_fit_ness_fn, step_size_fn, report_fn, mode, Lamarck, EVOLUTION, pmut, CR, number_of_clones, init_step_size, number_of_steps, decay,
               inertia, phi_p, phi_g, C_SSO, max_iteration, max_number_of_seeds, min_number_of_seeds, exponent, sigma_initial, sigma_final, id_of_search);

        //////////////////////////////////////////////////////////////////////////
        while (counter_till_end < patience && global_best >= LIMIT)
        {
            //cout << endl;

            for (int i = 0; i < iteration_block; ++i)
            {
                generation++;
                //iteration_block++;

                if (EVOLUTION)
                {
                    new_pop = old_pop.new_pop(size_of_feature, number_of_feature, delta, number_of_gene, number_of_chrom, number_of_clones, mode, pmut, CR, step_size, number_of_steps, Lamarck,
                                              number_of_fit_func_call);
                    old_pop = new_pop;
                    old_pop.save_fit(fit_ness_fn);
                    double new_glob_best = old_pop.geti(0).get_fit();
                    update_glob_best(new_glob_best, global_best, counter_till_end);
                    global_best = new_glob_best;

                    //cout << "generation: " << generation << endl;
                    //cout << "best_fit_value: " << old_pop.geti(0).get_fit() << endl << endl;
                }

                if (SWARM)
                {
                    Iteration.update(inertia, phi_p, phi_g, max_iteration, generation, max_number_of_seeds, min_number_of_seeds, exponent, sigma_initial, sigma_final, global_best_chrom, step_size,
                                     number_of_steps, Lamarck, number_of_fit_func_call);
                    Iteration.update_global_best(global_best, global_best_chrom, counter_till_end);
                    Iteration.save_fit(fit_ness_fn);
                    save_gb_fit(g_b_fit_ness_fn, global_best);

                    //cout << "iteration: " << generation << endl;
                    //cout << "best_fit_value: " << global_best << endl << endl;
                }

                if (Lamarck)
                {
                    save_step_size(step_size_fn, step_size);
                }
            }

            if (counter_till_end > 2 && Lamarck)
            {
                step_size *= decay;
            }

            if (EVOLUTION)
            {
                old_pop.geti(0).write_chrom(chrom_dir_name + "/" + to_string(generation) + ".dat");
            }

            if (SWARM)
            {
                global_best_chrom.write_chrom(chrom_dir_name + "/" + to_string(generation) + ".dat");
            }
        }

        t_end = clock();
        running_time = t_end - t_start;
        running_time = ((float)running_time)/CLOCKS_PER_SEC;

        make_report_file(report_fn, size_of_feature, number_of_feature, delta, number_of_gene, number_of_chrom, number_of_clones, mode, pmut, CR, step_size, number_of_steps, Lamarck, inertia,
                         phi_p, phi_g, max_iteration, generation, max_number_of_seeds, min_number_of_seeds, exponent, sigma_initial, sigma_final, EVOLUTION, running_time, counter_till_end, C_SSO,
                         init_step_size, global_best, number_of_fit_func_call);

        if (global_best < LIMIT)
        {
            limit_label = 1;
        }

        make_grid_search_report(grid_report_fn, number_of_clones, pmut, CR, step_size, number_of_steps, inertia,
                                phi_p, phi_g, max_iteration, generation, max_number_of_seeds, min_number_of_seeds, exponent, sigma_initial, sigma_final, running_time, C_SSO,
                                init_step_size, global_best, limit_label, id_of_search, number_of_fit_func_call);

        Parameter_variation(mode, Lamarck, EVOLUTION, Grid_Search,
                            number_of_clones, init_number_of_clones, pmut, init_pmut, CR, init_CR, number_of_steps, init_number_of_steps, inertia, init_inertia, phi_p, init_phi_p, phi_g, init_phi_g,
                            max_iteration, init_max_iteration, max_number_of_seeds, init_max_number_of_seeds, min_number_of_seeds, init_min_number_of_seeds, exponent, init_exponent,
                            sigma_initial, init_sigma_initial, sigma_final, init_sigma_final, C_SSO, init_C_SSO, init_step_size, init_init_step_size, decay, init_decay);

        if (no_grid_search)
        {
            Grid_Search = false;
        }
        cout << "id_of_search: " << id_of_search << endl;
        id_of_search++;
    }

 return 0;
}