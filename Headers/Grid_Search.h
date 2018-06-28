//
// Created by biot on 2018.04.07..
//

#ifndef EVOLUTIONARY_ALGORITHM_GRID_SEARCH_H
#define EVOLUTIONARY_ALGORITHM_GRID_SEARCH_H
void Parameter_variation(const int mode, bool Lamarck, int Evol, bool& Grid_Search,
                         int& number_of_clones, int& init_number_of_clones, double& pmut, double& init_pmut, double& CR, double& init_CR, int& number_of_steps, int& init_number_of_steps,
                         double& inertia, double& init_inertia, double& phi_p, double& init_phi_p, double& phi_g, double& init_phi_g, int& max_iteration, int& init_max_iteration, int& max_number_of_seeds,
                         int& init_max_number_of_seeds, int& min_number_of_seeds, int& init_min_number_of_seeds, double& exponent, double& init_exponent, double& sigma_initial, double& init_sigma_initial,
                         double& sigma_final,double& init_sigma_final, vector<double>& C_SSO, vector<double>& init_C_SSO, double& init_step_size, double& init_init_step_size, double& decay, double& init_decay)
{
    //GA
    double pmut_max = 1;
    double pmut_step = 0.05;

    //BEA
    int number_of_clones_max = 6;
    int number_of_clones_step = 1;

    //DE
    double CR_max = 1;
    double CR_step = 0.1;

    //////////////////////////////////////////////////////////////
    //PSO
    double inertia_max = 1;
    double inertia_step = 0.3;

    double phi_p_max = 1;
    double phi_p_step = 0.3;

    double phi_g_max = 1;
    double phi_g_step = 0.3;

    //SSO
    double C_w_max = C_SSO[1];
    double C_w_step = 0.1;

    double C_p_max = C_SSO[2];
    double C_p_step = 0.1;

    double C_g_max = 1;
    double C_g_step = 0.1;

    //IWO
    int max_iteration_max = 5000;
    int max_iteration_step = 5000;

    int max_number_of_seeds_max = 16;
    int max_number_of_seeds_step = 5;

    int min_number_of_seeds_max = max_number_of_seeds;
    int min_number_of_seeds_step = 5;

    double exponent_max = 2;
    double exponent_step = 0.5;

    double sigma_initial_max = 200;
    double sigma_initial_step = 50;

    double sigma_final_max = 10; //sigma_initial;
    double sigma_final_step = 10;

    //////////////////////////////////////////////////////////////
    //Lamarck
    double init_step_size_max = 1;
    double init_step_size_step = 1;

    int number_of_steps_max = 1;
    int number_of_steps_step = 1;

    double decay_max = 0.1;
    double decay_step = 0.05;


    if (Evol)
    {
        switch (mode)
        {
            case 0://GA

                if (floorf(pmut * 1000) < floorf(pmut_max * 1000))
                {
                    pmut += pmut_step;
                }
                else  if (floorf(pmut * 1000) >= floorf(pmut_max * 1000) && Lamarck)
                {
                    pmut = init_pmut;

                    if (floorf(init_step_size * 1000) < floorf(init_step_size_max * 1000))
                    {
                        init_step_size += init_step_size_step;
                    }
                    else if (floorf(init_step_size * 1000) >= floorf(init_step_size_max * 1000))
                    {
                        init_step_size = init_init_step_size;

                        if (number_of_steps < number_of_steps_max)
                        {
                            number_of_steps += number_of_steps_step;
                        }
                        else if (number_of_steps >= number_of_steps_max)
                        {
                            number_of_steps = init_number_of_steps;

                            if (floorf(decay * 1000) < floorf(decay_max * 1000))
                            {
                                decay += decay_step;
                            }
                            else if (floorf(decay * 1000) >= floorf(decay_max * 1000))
                            {
                                Grid_Search = false;
                            }
                        }
                    }
                }
                else
                {
                    Grid_Search = false;
                }

                break;
            case 1://BEA

                if (number_of_clones < number_of_clones_max)
                {
                    number_of_clones += number_of_clones_step;
                }
                else  if (floorf(number_of_clones * 1000) >= floorf(number_of_clones * 1000) && Lamarck)
                {
                    number_of_clones = init_number_of_clones;

                    if (floorf(init_step_size * 1000) < floorf(init_step_size_max * 1000))
                    {
                        init_step_size += init_step_size_step;
                    }
                    else if (floorf(init_step_size * 1000) >= floorf(init_step_size_max * 1000))
                    {
                        init_step_size = init_init_step_size;

                        if (number_of_steps < number_of_steps_max)
                        {
                            number_of_steps += number_of_steps_step;
                        }
                        else if (number_of_steps >= number_of_steps_max)
                        {
                            number_of_steps = init_number_of_steps;

                            if (floorf(decay * 1000) < floorf(decay_max * 1000))
                            {
                                decay += decay_step;
                            }
                            else if (floorf(decay * 1000) >= floorf(decay_max * 1000))
                            {
                                Grid_Search = false;
                            }
                        }
                    }
                }
                else
                {
                    Grid_Search = false;
                }

                break;
            case 2://DE

                if (floorf(CR * 1000) < floorf(CR_max * 1000))
                {
                    CR += CR_step;
                }
                else  if (floorf(CR * 1000) >= floorf(CR_max * 1000) && Lamarck)
                {
                    CR = init_CR;

                    if (floorf(init_step_size * 1000) < floorf(init_step_size_max * 1000))
                    {
                        init_step_size += init_step_size_step;
                    }
                    else if (floorf(init_step_size * 1000) >= floorf(init_step_size_max * 1000))
                    {
                        init_step_size = init_init_step_size;

                        if (number_of_steps < number_of_steps_max)
                        {
                            number_of_steps += number_of_steps_step;
                        }
                        else if (number_of_steps >= number_of_steps_max)
                        {
                            number_of_steps = init_number_of_steps;

                            if (floorf(decay * 1000) < floorf(decay_max * 1000))
                            {
                                decay += decay_step;
                            }
                            else if (floorf(decay * 1000) >= floorf(decay_max * 1000))
                            {
                                Grid_Search = false;
                            }
                        }
                    }
                }
                else
                {
                    Grid_Search = false;
                }

                break;
            case -1:
                if (Lamarck)
                {
                    if (floorf(init_step_size * 1000) < floorf(init_step_size_max * 1000))
                    {
                        init_step_size += init_step_size_step;
                    }
                    else if (floorf(init_step_size * 1000) >= floorf(init_step_size_max * 1000))
                    {
                        init_step_size = init_init_step_size;

                        if (number_of_steps < number_of_steps_max)
                        {
                            number_of_steps += number_of_steps_step;
                        }
                        else if (number_of_steps >= number_of_steps_max)
                        {
                            number_of_steps = init_number_of_steps;

                            if (floorf(decay * 1000) < floorf(decay_max * 1000))
                            {
                                decay += decay_step;
                            }
                            else if (floorf(decay * 1000) >= floorf(decay_max * 1000))
                            {
                                Grid_Search = false;
                            }
                        }
                    }
                }
                break;
        }

    }
    else
    {
        switch (mode)
        {
            case 0://PSO
                if (floorf(inertia * 1000) < floorf(inertia_max * 1000))
                {
                    inertia += inertia_step;
                }
                else if (floorf(inertia * 1000) >= floorf(inertia_max * 1000))
                {
                    inertia = init_inertia;

                    if (floorf(phi_p * 1000) < floorf(phi_p_max * 1000))
                    {
                        phi_p += phi_p_step;
                    }
                    else if (floorf(phi_p * 1000) >= floorf(phi_p_max * 1000))
                    {
                        phi_p = init_phi_p;

                        if (floorf(phi_g * 1000) < floorf(phi_g_max * 1000))
                        {
                            phi_g += phi_g_step;
                        }
                        else  if (floorf(phi_g * 1000) >= floorf(phi_g * 1000) && Lamarck)
                        {
                            phi_g = init_phi_g;

                            if (floorf(init_step_size * 1000) < floorf(init_step_size_max * 1000))
                            {
                                init_step_size += init_step_size_step;
                            }
                            else if (floorf(init_step_size * 1000) >= floorf(init_step_size_max * 1000))
                            {
                                init_step_size = init_init_step_size;

                                if (number_of_steps < number_of_steps_max)
                                {
                                    number_of_steps += number_of_steps_step;
                                }
                                else if (number_of_steps >= number_of_steps_max)
                                {
                                    number_of_steps = init_number_of_steps;

                                    if (floorf(decay * 1000) < floorf(decay_max * 1000))
                                    {
                                        decay += decay_step;
                                    }
                                    else if (floorf(decay * 1000) >= floorf(decay_max * 1000))
                                    {
                                        Grid_Search = false;
                                    }
                                }
                            }
                        }
                        else
                        {
                            Grid_Search = false;
                        }
                    }
                }
                break;
            case 1://SSO
                if (floorf(C_SSO[0] * 1000) < floorf(C_w_max * 1000))
                {
                    C_SSO[0] += C_w_step;
                }
                else if (floorf(C_SSO[0] * 1000) >= floorf(C_w_max * 1000))
                {
                    C_SSO[0] = init_C_SSO[0];

                    if (floorf(C_SSO[1] * 1000) < floorf(C_p_max * 1000))
                    {
                        C_SSO[1] += C_p_step;
                    }
                    else if (floorf(C_SSO[1] * 1000) >= floorf(C_p_max * 1000))
                    {
                        C_SSO[1] = init_C_SSO[1];

                        if (floorf(C_SSO[2] * 1000) < floorf(C_g_max * 1000))
                        {
                            C_SSO[2] += C_g_step;
                        }
                        else if (floorf(C_SSO[2] * 1000) >= floorf(C_g_max * 1000) && Lamarck)
                        {
                            C_SSO[2] = init_C_SSO[2];

                            if (floorf(init_step_size * 1000) < floorf(init_step_size_max * 1000))
                            {
                                init_step_size += init_step_size_step;
                            }
                            else if (floorf(init_step_size * 1000) >= floorf(init_step_size_max * 1000))
                            {
                                init_step_size = init_init_step_size;

                                if (number_of_steps < number_of_steps_max)
                                {
                                    number_of_steps += number_of_steps_step;
                                }
                                else if (number_of_steps >= number_of_steps_max)
                                {
                                    number_of_steps = init_number_of_steps;

                                    if (floorf(decay * 1000) < floorf(decay_max * 1000))
                                    {
                                        decay += decay_step;
                                    }
                                    else if (floorf(decay * 1000) >= floorf(decay_max * 1000))
                                    {
                                        Grid_Search = false;
                                    }
                                }
                            }
                        }
                        else
                        {
                            Grid_Search = false;
                        }
                    }
                }
                break;
            case 2://IWO
                if (floorf(max_iteration * 1000) < floorf(max_iteration_max * 1000))
                {
                    max_iteration += max_iteration_step;
                    }
                    else if (floorf(max_iteration * 1000) >= floorf(max_iteration_max * 1000))
                    {
                        max_iteration = init_max_iteration;

                        if (floorf(min_number_of_seeds * 1000) < floorf(min_number_of_seeds_max * 1000))
                        {
                            min_number_of_seeds += min_number_of_seeds_step;
                        }
                        else if (floorf(min_number_of_seeds * 1000) >= floorf(min_number_of_seeds_max * 1000))
                        {
                            min_number_of_seeds = init_min_number_of_seeds;

                            if (floorf(max_number_of_seeds * 1000) < floorf(max_number_of_seeds_max * 1000))
                            {
                                max_number_of_seeds += max_number_of_seeds_step;
                            }
                            else if (floorf(max_number_of_seeds * 1000) >= floorf(max_number_of_seeds_max * 1000))
                            {
                                max_number_of_seeds = init_max_number_of_seeds;

                                if (floorf(exponent * 1000) < floorf(exponent_max * 1000))
                                {
                                    exponent += exponent_step;
                                }
                                else if (floorf(exponent * 1000) >= floorf(exponent_max * 1000))
                                {
                                    exponent = init_exponent;

                                    if (floorf(sigma_final * 1000) < floorf(sigma_final_max * 1000))
                                    {
                                            sigma_final += sigma_final_step;
                                    }
                                    else if (floorf(sigma_final * 1000) >= floorf(sigma_final_max * 1000)) {
                                        sigma_final = init_sigma_final;

                                        if (floorf(sigma_initial * 1000) < floorf(sigma_initial_max * 1000))
                                        {
                                            sigma_initial += sigma_initial_step;
                                        }
                                        else if (floorf(sigma_initial * 1000) >= floorf(sigma_initial_max * 1000) && Lamarck) {
                                            sigma_initial = init_sigma_initial;

                                            if (floorf(init_step_size * 1000) < floorf(init_step_size_max * 1000))
                                            {
                                                init_step_size += init_step_size_step;
                                            }
                                            else if (floorf(init_step_size * 1000) >= floorf(init_step_size_max * 1000))
                                            {
                                                init_step_size = init_init_step_size;

                                                if (number_of_steps < number_of_steps_max)
                                                {
                                                    number_of_steps += number_of_steps_step;
                                                }
                                                else if (number_of_steps >= number_of_steps_max) {
                                                    number_of_steps = init_number_of_steps;

                                                    if (floorf(decay * 1000) < floorf(decay_max * 1000))
                                                    {
                                                        decay += decay_step;
                                                    }
                                                    else if (floorf(decay * 1000) >= floorf(decay_max * 1000))
                                                    {
                                                        Grid_Search = false;
                                                    }
                                                }
                                            }
                                        }
                                        else
                                        {
                                            Grid_Search = false;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    break;
            case -1:
                        if (Lamarck)
                        {
                            if (floorf(init_step_size * 1000) < floorf(init_step_size_max * 1000))
                            {
                                init_step_size += init_step_size_step;
                            }
                            else if (floorf(init_step_size * 1000) >= floorf(init_step_size_max * 1000))
                            {
                                init_step_size = init_init_step_size;

                                if (number_of_steps < number_of_steps_max)
                                {
                                    number_of_steps += number_of_steps_step;
                                }
                                else if (number_of_steps >= number_of_steps_max)
                                {
                                    number_of_steps = init_number_of_steps;

                                    if (floorf(decay * 1000) < floorf(decay_max * 1000))
                                    {
                                        decay += decay_step;
                                    }
                                    else if (floorf(decay * 1000) >= floorf(decay_max * 1000))
                                    {
                                        Grid_Search = false;
                                    }
                                }
                            }
                        }
                break;
        }

    }
}
#endif //EVOLUTIONARY_ALGORITHM_GRID_SEARCH_H
