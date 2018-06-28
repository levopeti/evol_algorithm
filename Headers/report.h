//
// Created by biot on 2018.04.06..
//

#ifndef EVOLUTIONARY_ALGORITHM_REPORT_H
#define EVOLUTIONARY_ALGORITHM_REPORT_H
#include <string>
#include <sys/stat.h>
#include <iomanip>
#include <sstream>


//write the chromosome into a file with file_name
void save_gb_fit(string file_name, double& value)
{
    ofstream file(file_name.c_str(), ios::out | std::ios::app);
    file << value << endl;
    file.close();
}

//write the step size into a file with file_name
void save_step_size(string file_name, double& value)
{
    ofstream file(file_name.c_str(), ios::out | std::ios::app);
    file << value << endl;
    file.close();
}

//write the head of the grid report file
void init_make_grid_search_report(string file_name)
{
    ofstream file(file_name.c_str(), ios::out);// | std::ios::app);
    file << "id_of_search" << '\t' << "limit_label" << '\t' << "global_best" << '\t' << "init_step_size" << '\t'
         << "C_w" << '\t' << "C_p" << '\t' << "C_g" << '\t' << "running_time" << '\t' << "sigma_final" <<
         '\t' << "sigma_initial" << '\t' << "exponent" << '\t' << "min_number_of_seeds" << '\t' << "max_number_of_seeds"
         << '\t' << "iteration" << '\t' << "max_iteration" << '\t' << "phi_g" << '\t' << "phi_p" << '\t' <<
         "inertia" << '\t' << "number_of_steps" << '\t' << "step_size" << '\t' << "CR" << '\t' << "pmut" << '\t'
         << "number_of_clones" << '\t' << "number_of_fit_func_call" << endl;
    file.close();
}

//write the finish values of the algorithm into a file with file_name
void make_grid_search_report(string file_name, const int& number_of_clones, const double& pmut, const double& CR, double& step_size, const int& number_of_steps, const double& inertia,
                             const double& phi_p, const double& phi_g, const int& max_iteration, int& iteration, const int& max_number_of_seeds, const int& min_number_of_seeds, const double& exponent,
                             const double& sigma_initial, const double& sigma_final, double running_time, vector<double> C_SSO, double& init_step_size, double global_best, int limit_label, int id_of_search,
                             int number_of_fit_func_call)
{
    ofstream file(file_name.c_str(), ios::out | std::ios::app);

    file << id_of_search << '\t' << limit_label << '\t' << global_best << '\t' << init_step_size << '\t' << C_SSO[0] << '\t' << C_SSO[1] << '\t' << C_SSO[2] << '\t' << running_time << '\t' << sigma_final <<
         '\t' << sigma_initial << '\t' << exponent << '\t' << min_number_of_seeds << '\t' << max_number_of_seeds << '\t' << iteration << '\t' << max_iteration << '\t' << phi_g << '\t' << phi_p << '\t' <<
         inertia << '\t' << number_of_steps << '\t' << step_size << '\t' << CR << '\t' << pmut << '\t' << number_of_clones << '\t' << number_of_fit_func_call << endl;

    file.close();
}

void make_grid_base_dir_name(string& file_name)
{
    // current date/time based on current system
    time_t now = time(0);

    // convert now to string form
    tm *ltm = localtime(&now);
    stringstream dt;
    dt << (1900 + ltm->tm_year) << '_' << 1 + ltm->tm_mon<< '_' <<ltm->tm_mday<< '_' << ltm->tm_hour << '_' << ltm->tm_min << '/';
    file_name = file_name + dt.str();
}

stringstream _mode_to_file(int mode, bool Lamarck, int Evol,
                           double pmut, double CR, int number_of_clones, double init_step_size, int number_of_steps, double decay,
                           double inertia, double phi_p, double phi_g, vector<double> C_SSO, int max_iteration, int max_number_of_seeds, int min_number_of_seeds, double exponent,
                           double sigma_initial, double sigma_final)
{
    stringstream mode_s;
    mode_s << fixed << setprecision(2);

    if (Evol)
    {
        switch (mode)
        {
            case 0:
                mode_s << "GA--pmut:" << pmut;
                break;
            case 1:
                mode_s << "BEA--n_o_c:" << number_of_clones;
                break;
            case 2:
                mode_s << "DE--CR:" << CR;
                break;
            case -1:
                mode_s << "GD";
                break;
        }

    }
    else
    {
        switch (mode)
        {
            case 0:
                mode_s << "PSO--iner:" << inertia << "--phi_p:" << phi_p << "--phi_g:" << phi_g;
                break;
            case 1:
                mode_s << "SSO--C_w:" << C_SSO[0] << "--C_p:" << C_SSO[1] << "--C_g:" << C_SSO[2];
                break;
            case 2:
                mode_s << "IWO--max_it:" << max_iteration << "--max_n_o_s:" << max_number_of_seeds << "--min_n_o_s:" << min_number_of_seeds << "--exp:" << exponent <<
                       "--sig_init:" << sigma_initial << "--sig_fin:" << sigma_final;
                break;
            case -1:
                mode_s << "GD";
                break;
        }
    }

    if (Lamarck)
    {
        mode_s << "--GD:i_s_s:" << init_step_size << "--n_o_s:" << number_of_steps << "--dec:" << decay;
    }

    return mode_s;
}

stringstream _mode_to_report(const vector<double>& size_of_feature, const int& number_of_feature, const vector<double>& delta, const int& number_of_gene, const int& number_of_chrom,
                             const int& number_of_clones, const int& mode, const double& pmut, const double& CR, double& step_size, const int& number_of_steps, bool& Lamarck, const double& inertia,
                             const double& phi_p, const double& phi_g, const int& max_iteration, int& iteration, const int& max_number_of_seeds, const int& min_number_of_seeds, const double& exponent,
                             const double& sigma_initial, const double& sigma_final, int Evol, vector<double> C_SSO, double& init_step_size)
{
    stringstream mode_s;

    if (Evol)
    {
        switch (mode)
        {
            case 0:
                mode_s << fixed << setprecision(2) << "Genetic Algorithm \nnumber of chromosomes: " << number_of_chrom << "\nnumber of genes: " << number_of_gene << "\nnumber of features: "<<
                       number_of_feature  << "\nprobability of mutation: " << pmut;
                break;
            case 1:
                mode_s << fixed << setprecision(2) << "Bacterial Algorithm \nnumber of chromosomes: " << number_of_chrom << "\nnumber of genes: " << number_of_gene << "\nnumber of features: " <<
                       number_of_feature << "\nnumber of clones: " << number_of_clones;
                break;
            case 2:
                mode_s << fixed << setprecision(2) << "Deferential Evolution \nnumber of chromosomes: " << number_of_chrom << "\nnumber of genes: " << number_of_gene << "\nnumber of features: " <<
                       number_of_feature << "\nprobability of recombination: " << CR;
                break;
            case -1:
                mode_s << fixed << setprecision(2) << "Only Gradient Descent \nnumber of chromosomes: " << number_of_chrom << "\nnumber of genes: " << number_of_gene << "\nnumber of features: " <<
                       number_of_feature << "\ninit step size: " << init_step_size << "\nend step size: " << scientific << step_size << "\nnumber of steps: " << fixed << setprecision(2) << number_of_steps;
                break;
        }

    }
    else
    {
        switch (mode)
        {
            case 0:
                mode_s << fixed << setprecision(2) << "Particle Swarm Optimization \nnumber of chromosomes: " << number_of_chrom << "\nnumber of genes: " << number_of_gene << "\nnumber of features: " <<
                       number_of_feature << "\ninertia: " << inertia << "\nphi p: " << phi_p << "\nphi g: " << phi_g;
                break;
            case 1:
                mode_s << fixed << setprecision(2) << "Simplified Swarm Optimization \nnumber of chromosomes: " << number_of_chrom << "\nnumber of genes: " << number_of_gene << "\nnumber of features: " <<
                       number_of_feature << "\nC w: " << C_SSO[0] << "\nC p: " << C_SSO[1] << "\nC g: " << C_SSO[2];
                break;
            case 2:
                mode_s << fixed << setprecision(2) << "Invasive Weed Optimization \nnumber of chromosomes: " << number_of_chrom << "\nnumber of genes: " << number_of_gene << "\nnumber of features: " <<
                       number_of_feature << "\nmax iteration: " << max_iteration << "\nmax number of seeds: " << max_number_of_seeds << "\nmin number of seeds: " << min_number_of_seeds << "\nexponent: " <<
                       exponent << "\nsigma initial: " << sigma_initial << "\nsigma_final: " << sigma_final;
                break;
            case -1:
                mode_s << fixed << setprecision(2) << "Only Gradient Descent \nnumber of chromosomes: " << number_of_chrom << "\nnumber of genes: " << number_of_gene << "\nnumber of features: " <<
                       number_of_feature << "\ninit step size: " << init_step_size << "\nend step size: " << scientific << step_size << "\nnumber of steps: " << fixed << setprecision(2) << number_of_steps;
                break;
        }
    }

    if (Lamarck)
    {
        if (mode != -1)
        {
            mode_s << fixed << setprecision(2) << "\n\nGradient Descent"  << "\ninit step size: " << init_step_size << "\nend step size: " << scientific << step_size <<
                   "\nnumber of steps: " <<
                   fixed << setprecision(2) << number_of_steps;
        }
    }

    return mode_s;
}

void _init_(string& grid_report_fn, string& base_dir_name, string& chrom_dir_name, string& fit_ness_fn, string& g_b_fit_ness_fn, string& step_size_fn, string& report_fn, int mode, bool Lamarck, int Evol,
            double pmut, double CR, int number_of_clones, double init_step_size, int number_of_steps, double decay, double inertia, double phi_p, double phi_g, vector<double> C_SSO, int max_iteration,
            int max_number_of_seeds, int min_number_of_seeds, double exponent, double sigma_initial, double sigma_final, int id_of_search)
{
    stringstream mode_s;
    mode_s << "__" << id_of_search << "__";
    stringstream mode_add = _mode_to_file(mode, Lamarck, Evol, pmut, CR, number_of_clones, init_step_size, number_of_steps, decay,
                                        inertia, phi_p, phi_g, C_SSO, max_iteration, max_number_of_seeds, min_number_of_seeds, exponent, sigma_initial, sigma_final);
    mode_s << mode_add.str();
    string current_dir = base_dir_name + mode_s.str();
    chrom_dir_name = current_dir + "/Chromosomes/";
    mkdir(current_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(chrom_dir_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    fit_ness_fn = current_dir + "/fit_ness.dat";
    g_b_fit_ness_fn = current_dir + "/g_best_fit_ness.dat";
    step_size_fn = current_dir + "/gd_step_size.dat";
    report_fn = current_dir + "/report.dat";
}

void make_report_file(string file_name, const vector<double>& size_of_feature, const int& number_of_feature, const vector<double>& delta, const int& number_of_gene, const int& number_of_chrom,
                      const int& number_of_clones, const int& mode, const double& pmut, const double& CR, double& step_size, const int& number_of_steps, bool& Lamarck, const double& inertia,
                      const double& phi_p, const double& phi_g, const int& max_iteration, int& iteration, const int& max_number_of_seeds, const int& min_number_of_seeds, const double& exponent,
                      const double& sigma_initial, const double& sigma_final, int Evol, double running_time, int counter_till_end, vector<double> C_SSO, double& init_step_size, double global_best,
                      int number_of_fit_func_call)
{
    ofstream file(file_name.c_str(), ios::out);
    stringstream mode_s = _mode_to_report(size_of_feature, number_of_feature, delta, number_of_gene, number_of_chrom, number_of_clones, mode, pmut, CR, step_size, number_of_steps, Lamarck, inertia, phi_p,
                                          phi_g, max_iteration, iteration, max_number_of_seeds, min_number_of_seeds, exponent, sigma_initial, sigma_final, Evol, C_SSO, init_step_size);
    mode_s << fixed << setprecision(2) << "\nrunning time: " << running_time << " sec" << "\nnumber of iteration: " << iteration << scientific << "\nbest fitness value: " << global_best
           << "\nnumber_of_fit_func_call: " << number_of_fit_func_call;
    string report = mode_s.str();
    file << report << endl;

    file.close();
}

void update_glob_best(double& new_global_best, double& old_global_best, int& counter_till_end)
{
    counter_till_end++;

    if (new_global_best < old_global_best)
    {
        counter_till_end = 0;
    }
}

#endif //EVOLUTIONARY_ALGORITHM_REPORT_H
