#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <filesystem>
#include <chrono>
#include <set>

#include "constants.hpp"
#include "integration.hpp"
#include "interpolation.hpp"
#include "file_read.hpp"
#include "file_write.hpp"
#include "process_file.hpp"
#include "utilities.hpp"
#include "runtime_variables.hpp"
#include "linearize.hpp"
#include "energy_grid.hpp"

#include <iostream>
#include <vector>

namespace {
    int log_depth = 0;
    struct LogIndent {
        LogIndent()  { ++log_depth; }
        ~LogIndent() { --log_depth; }
    };
    std::string log_prefix() { return std::string(log_depth * 2, ' '); }
}

void read_cdf_file__(std::string const & file_path, std::vector<double> &grid){
    LogIndent _li;
    std::cout << log_prefix() << "[read_cdf_file] Reading: " << file_path << std::endl;
    std::fstream file(file_path);
    std::string line;
    while (getline(file, line))
    {
        // right trim the line, should only be one value per line
        line.erase(line.find_last_not_of(" \n\r\t")+1);
        grid.push_back(std::stod(line));
    }
    std::cout << log_prefix() << "[read_cdf_file] Read " << grid.size() << " points" << std::endl;
}

// Class constructor
DistData::DistData(TslFileData& file_data) : tsl_data(file_data){
    LogIndent _li;
    beta_grid = tsl_data.return_betas();
    calculation_alphas = tsl_data.return_scaled_alphas();
    calculation_betas = tsl_data.return_full_scaled_betas();
    calculation_half_betas = tsl_data.return_scaled_betas();
    calculation_tsl_vals = tsl_data.return_tsl_vals();
    if (tsl_data.lasym != 0){throw std::runtime_error("Asymmetric TSLs given.");}

    set_interp_integration_schemes__();
    set_initial_cdf_grids__();
    if (!silence){
        std::cout << log_prefix() << "[DistData] Initialized: " << calculation_alphas.size() << " alphas, "
                  << calculation_half_betas.size() << " half-betas, "
                  << calculation_betas.size() << " full betas" << std::endl;
    }
}

void DistData::set_interp_integration_schemes__(){
    LogIndent _li;
    /// NOTE: Due to the handling of TSL data, these schemes may be unreliable.
    /// Alpha interpolation and integration will be invalid as long as lln is set in the tsl data (ln(s) instead of s is stored) and interpolation cannot be altered
    /// Some considerations have been taken but not guaranteed to be valid
    /// Beta interpolation and integration will be invalid as long as lasym is set in the tsl data (symmetric for of s in stored)
    /// All beta interpolation should be handled in file_read.cpp and integration should be done with lin-lin due to linearization routines
    alpha_interpolation_scheme = tsl_data.return_alpha_schemes();
    alpha_integration_scheme = alpha_interpolation_scheme;
    beta_interpolation_scheme = tsl_data.return_beta_schemes();
    /// NOTE: Beta arrays are linearized
    // beta_integration_scheme = beta_interpolation_scheme;
    beta_integration_scheme = 2;
    if (!silence){
        std::cout << log_prefix() << "[set_interp_integration_schemes] alpha_interp=" << alpha_interpolation_scheme
                  << " beta_interp=" << beta_interpolation_scheme
                  << " beta_integ=" << beta_integration_scheme << std::endl;
    }
}

void DistData::set_initial_cdf_grids__(){
    LogIndent _li;
    if (linearize_cdfs){
        initialization_beta_cdf_grid.push_back(0.0001);
        initialization_beta_cdf_grid.push_back(0.9999);
        initialization_alpha_cdf_grid.push_back(0.0001);
        initialization_alpha_cdf_grid.push_back(0.9999);
        if (!silence){std::cout << log_prefix() << "[set_initial_cdf_grids] mode=linearize (2-point seed)" << std::endl;}
    }
    else if (use_sigmoid_cdfs)
    {
        initialization_beta_cdf_grid = sigmoid_space(0, 1, num_beta_cdf_points + 2, beta_cdf_extent);
        initialization_alpha_cdf_grid = sigmoid_space(0, 1, num_alpha_cdf_points + 2, alpha_cdf_extent);

        // Trim off the 0 and 1 for the cdf grids
        // Including them causes trouble with the fitting across multiple temperatures as finding fit points at 0 and 1 cause
        // large spikes as the edges of the fit values
        initialization_beta_cdf_grid.erase(initialization_beta_cdf_grid.begin());
        initialization_beta_cdf_grid.erase(initialization_beta_cdf_grid.end()-1);
        initialization_alpha_cdf_grid.erase(initialization_alpha_cdf_grid.begin());
        initialization_alpha_cdf_grid.erase(initialization_alpha_cdf_grid.end()-1);
        if (!silence){
            std::cout << log_prefix() << "[set_initial_cdf_grids] mode=sigmoid beta_pts=" << initialization_beta_cdf_grid.size()
                      << " alpha_pts=" << initialization_alpha_cdf_grid.size() << std::endl;
        }
    }
    else
    {
        read_cdf_file__(alpha_cdf_grid_loc, initialization_alpha_cdf_grid);
        read_cdf_file__(beta_cdf_grid_loc, initialization_beta_cdf_grid);
        if (!silence){
            std::cout << log_prefix() << "[set_initial_cdf_grids] mode=file beta_pts=" << initialization_beta_cdf_grid.size()
                      << " alpha_pts=" << initialization_alpha_cdf_grid.size() << std::endl;
        }
    }

}

std::pair<double, bool> DistData::return_arbitrary_TSL_val(double alpha, double beta)
{
    // NOTE: Not logged -- called in the innermost interpolation hot path (millions of calls)
    beta = abs(beta);
    bool off_data = (
        alpha > calculation_alphas.back() ||
        beta  > calculation_half_betas.back()
        );
    if (off_data){
        return std::make_pair(tsl_data.return_asym_SCT(alpha, beta), true);
    }

    int alpha_lo_insert = std::lower_bound(calculation_alphas.begin()+1, calculation_alphas.end(), alpha) - calculation_alphas.begin();
    int beta_lo_insert = std::lower_bound(calculation_half_betas.begin()+1, calculation_half_betas.end(), beta) - calculation_half_betas.begin();
    double& f11 = calculation_tsl_vals[beta_lo_insert-1][alpha_lo_insert-1];
    double& f12 = calculation_tsl_vals[beta_lo_insert-1][alpha_lo_insert];
    double& f21 = calculation_tsl_vals[beta_lo_insert][alpha_lo_insert-1];
    double& f22 = calculation_tsl_vals[beta_lo_insert][alpha_lo_insert];
    bool below_cutoff = (f11<sct_cutoff || f12<sct_cutoff || f21<sct_cutoff || f22<sct_cutoff);
    if (below_cutoff){
        return std::make_pair(tsl_data.return_asym_SCT(alpha, beta), true);
    }
    double& b_l = calculation_half_betas[beta_lo_insert-1];
    double& b_u = calculation_half_betas[beta_lo_insert];
    double& a_l = calculation_alphas[alpha_lo_insert-1];
    double& a_u = calculation_alphas[alpha_lo_insert];
    return std::make_pair(bi_interp(b_l, b_u, a_l, a_u,
                                     f11, f12, f21, f22,
                                     beta, alpha,
                                     beta_interpolation_scheme, alpha_interpolation_scheme),
                                     false);
}

std::vector<double> DistData::get_viable_betas__(double const& inc_energy){
    LogIndent _li;
    double b_min = tsl_data.calculate_beta_min(inc_energy);
    double b_max = tsl_data.calculate_beta_max(inc_energy);
    auto lo_insert = std::lower_bound(calculation_betas.begin(), calculation_betas.end(), b_min) - calculation_betas.begin();
    auto hi_insert = std::lower_bound(calculation_betas.begin(), calculation_betas.end(), b_max) - calculation_betas.begin();
    std::vector<double> result;
    result.reserve(hi_insert - lo_insert + 2);
    result.push_back(b_min);
    for (int i=lo_insert; i<hi_insert; i++){
        result.push_back(calculation_betas[i]);
    }
    result.push_back(b_max);
    if (!silence){
        std::cout << log_prefix() << "[get_viable_betas] E=" << inc_energy
                  << " b_min=" << b_min << " b_max=" << b_max << " n=" << result.size() << std::endl;
    }
    return result;
}

std::vector<double> DistData::get_viable_alphas__(double const &inc_energy, double const &beta){
    LogIndent _li;
    auto [a_min, a_max] = tsl_data.calculate_alpha_extrema(inc_energy, beta);
    auto lo_insert = std::lower_bound(calculation_alphas.begin(), calculation_alphas.end(), a_min) - calculation_alphas.begin();
    auto hi_insert = std::lower_bound(calculation_alphas.begin(), calculation_alphas.end(), a_max) - calculation_alphas.begin();
    std::vector<double> result;
    result.reserve(hi_insert - lo_insert + 2);
    result.push_back(a_min);
    for (int i=lo_insert; i<hi_insert; i++){
        result.push_back(calculation_alphas[i]);
    }
    result.push_back(a_max);
    if (!silence){
        std::cout << log_prefix() << "[get_viable_alphas] E=" << inc_energy << " beta=" << beta
                  << " a_min=" << a_min << " a_max=" << a_max << " n=" << result.size() << std::endl;
    }
    return result;
}

std::pair<std::vector<double>, std::vector<bool>> DistData::get_alpha_line__(std::vector<double> const& alpha_vals, double const& beta){
    // NOTE: Not logged -- thin loop over return_arbitrary_TSL_val, in hot path
    std::vector<double> vals;
    std::vector<bool> truthy;
    vals.reserve(alpha_vals.size());
    truthy.reserve(alpha_vals.size());
    for(double alpha:alpha_vals){
        auto [val, thruth] = return_arbitrary_TSL_val(alpha, beta);
        vals.push_back(val);
        truthy.push_back(thruth);
    }
    return std::make_pair(vals, truthy);
}

double DistData::integrate_alpha_line__(std::vector<double> const& alpha_vals, std::vector<double> const& vals, std::vector<bool> const& truthy, double const& beta){
    // NOTE: Not logged -- integration kernel, in hot path
    double integral = 0;
    for(int i=0; i<alpha_vals.size()-1; i++){
        if (truthy[i] || truthy[i+1]){
            integral += tsl_data.return_asym_SCT_alpha_integral(alpha_vals[i], alpha_vals[i+1], beta);
        }
        else{
            /// NOTE: This function will probably need to be moved inside tsl_data in order to properly integrate
            integral += ENDF_integrate(alpha_vals[i], alpha_vals[i+1], vals[i], vals[i+1], alpha_integration_scheme);
        }
    }
    return integral;
}

double DistData::get_beta_pdf_val__(double const& inc_energy, double const& beta){
    LogIndent _li;
    std::vector<double> alpha_vals = get_viable_alphas__(inc_energy, beta);
    auto [vals, truthy] = get_alpha_line__(alpha_vals, beta);
    // Multiplying by exp(-b/2) makes asymmetric
    double result = exp(-beta/2) * integrate_alpha_line__(alpha_vals, vals, truthy, beta);
    // NOTE: Each log line here is one linearizer callback evaluation -- call count is diagnostic
    if (!silence){
        std::cout << log_prefix() << "[get_beta_pdf_val] E=" << inc_energy << " beta=" << beta << " pdf=" << result << std::endl;
    }
    return result;
}

std::pair<std::vector<double>, std::vector<double>> DistData::return_beta_pdf(double const& inc_energy){
    LogIndent _li;
    std::vector<double> beta_vals = get_viable_betas__(inc_energy);
    std::vector<double> beta_pdf;
    beta_pdf.reserve(beta_vals.size());
    for (double beta: beta_vals){
        beta_pdf.push_back(get_beta_pdf_val__(inc_energy, beta));
    }
    if (!silence){
        std::cout << log_prefix() << "[return_beta_pdf] E=" << inc_energy << " n_betas=" << beta_vals.size() << std::endl;
    }
    return std::make_pair(beta_vals, beta_pdf);
}

std::pair<std::vector<double>, std::vector<double>> DistData::return_linearized_beta_pdf(double const& inc_energy){
    LogIndent _li;
    auto [beta_vals, beta_pdf] = return_beta_pdf(inc_energy);
    auto get_new_beta_pdf_val = [this, inc_energy](double beta) {return get_beta_pdf_val__(inc_energy, beta);};
    size_t pre_size = beta_vals.size();
    auto t0 = std::chrono::steady_clock::now();
    linearize(beta_vals, beta_pdf, get_new_beta_pdf_val, 1e-15, 1e-3, ToleranceCondition::Adaptive);
    auto t1 = std::chrono::steady_clock::now();
    if (!silence){
        std::cout << log_prefix() << "[return_linearized_beta_pdf] E=" << inc_energy
                  << " n: " << pre_size << " -> " << beta_vals.size()
                  << " (" << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << " ms)" << std::endl;
    }
    return std::make_pair(beta_vals, beta_pdf);
}

double DistData::return_ii_xs_value(double const& inc_energy){
    LogIndent _li;
    auto [beta_vals, beta_pdf] = return_linearized_beta_pdf(inc_energy);
    double xs = ((tsl_data.a0*boltz*tsl_data.temp*tsl_data.bound_xs)/(4*inc_energy))*ENDF_integrate_vector(beta_vals, beta_pdf, beta_integration_scheme);
    if (!silence){
        std::cout << log_prefix() << "[return_ii_xs_value] E=" << inc_energy << " xs=" << xs << std::endl;
    }
    return xs;
}

std::vector<double> DistData::return_ii_xs_vector(std::vector<double> const& inc_energies){
    LogIndent _li;
    std::vector<double> result(inc_energies.size());
    for (int i=0; i<inc_energies.size(); i++){
        if (!silence){
            std::cout << log_prefix() << "[return_ii_xs_vector] " << i+1 << "/" << inc_energies.size()
                      << " E=" << inc_energies[i] << std::endl;
        }
        result[i] = return_ii_xs_value(inc_energies[i]);
    }
    return result;
}

std::pair<std::vector<double>, std::vector<double>> DistData::return_linearized_ii_xs(){
    LogIndent _li;
    std::vector<double> energies = logspace(e_min, tsl_data.e_max, num_energies); // Don't know what the best choice for starting grid would be or how many points
    std::vector<double> xs = return_ii_xs_vector(energies);
    auto get_new_xs = [&](double x) {return return_ii_xs_value(x);};
    size_t pre_size = energies.size();
    auto t0 = std::chrono::steady_clock::now();
    linearize(energies, xs, get_new_xs, 1e-15, 1e-3, ToleranceCondition::Adaptive);
    auto t1 = std::chrono::steady_clock::now();
    if (!silence){
        std::cout << log_prefix() << "[return_linearized_ii_xs] n: " << pre_size << " -> " << energies.size()
                  << " (" << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << " ms)" << std::endl;
    }
    return std::make_pair(energies, xs);
}

std::pair<std::vector<double>, std::vector<double>> DistData::return_final_ii_xs(){
    LogIndent _li;
    if (use_external_energy_grid || use_internal_energy_grid)
    {
        if (!silence){
            std::cout << log_prefix() << "[return_final_ii_xs] mode=" << (use_external_energy_grid ? "external_grid" : "internal_grid") << std::endl;
        }
        std::vector<double> energies = return_incident_energy_grid();
        std::vector<double> xs = return_ii_xs_vector(energies);
        return std::make_pair(energies, xs);
    }
    else
    {
        if (!silence){std::cout << log_prefix() << "[return_final_ii_xs] mode=linearize" << std::endl;}
        return return_linearized_ii_xs();
    }
}

void DistData::get_beta_sampling_dists__(){
    LogIndent _li;
    calculation_beta_vals.reserve(incident_energy_grid.size());
    calculation_beta_cdfs.reserve(incident_energy_grid.size());
    beta_vals.reserve(incident_energy_grid.size());
    int energy_idx = 0;
    for (auto inc_energy: incident_energy_grid){
        if (!silence){
            std::cout << log_prefix() << "[get_beta_sampling_dists] energy " << energy_idx+1 << "/" << incident_energy_grid.size()
                      << " E=" << inc_energy << std::endl;
        }
        auto [vals, pdf] = return_linearized_beta_pdf(inc_energy);
        calculation_beta_vals.push_back(vals);
        calculation_beta_cdfs.push_back(pdf_to_cdf(vals, pdf));
        beta_vals.push_back(fit_cdf(vals, calculation_beta_cdfs.back(), initialization_beta_cdf_grid));
        beta_cdf_grid.push_back(initialization_beta_cdf_grid);
        energy_idx++;
    }
}

void DistData::linearize_beta_sampling_dists__(){
    LogIndent _li;
    for (int i = 0; i < calculation_beta_cdfs.size(); i++){
        // Lambda function to get a new beta value that corresponds to a desired cdf value
       auto get_new_cdf_point = [&](double desired_cdf_val){
            int cdf_insert = std::lower_bound(calculation_beta_cdfs[i].begin()+1, calculation_beta_cdfs[i].end(), desired_cdf_val) - calculation_beta_cdfs[i].begin();
            return ENDF_interp(calculation_beta_cdfs[i][cdf_insert],
                               calculation_beta_cdfs[i][cdf_insert+1],
                               calculation_beta_vals[i][cdf_insert],
                               calculation_beta_vals[i][cdf_insert+1],
                               desired_cdf_val);
        };
        size_t pre_size = beta_vals[i].size();
        auto t0 = std::chrono::steady_clock::now();
        linearize(beta_cdf_grid[i], beta_vals[i], get_new_cdf_point, 1e-15, 1e-3, ToleranceCondition::Adaptive);
        auto t1 = std::chrono::steady_clock::now();
        if (!silence){
            std::cout << log_prefix() << "[linearize_beta_sampling_dists] i=" << i << "/" << calculation_beta_cdfs.size()
                      << " n: " << pre_size << " -> " << beta_vals[i].size()
                      << " (" << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << " ms)" << std::endl;
        }
    }
}

std::pair<std::vector<double>, std::vector<double>> DistData::return_linearized_alpha_pdf(double const& beta){
    LogIndent _li;
    std::vector<double> a_vals = calculation_alphas;
    auto [alpha_pdf, truthy] = get_alpha_line__(a_vals, beta);
    auto get_new_alpha_pdf_val = [this, beta](double alpha) {return return_arbitrary_TSL_val(alpha, beta).first;};
    size_t pre_size = a_vals.size();
    auto t0 = std::chrono::steady_clock::now();
    linearize(a_vals, alpha_pdf, get_new_alpha_pdf_val, 1e-15, 1e-3, ToleranceCondition::Adaptive);
    auto t1 = std::chrono::steady_clock::now();
    if (!silence){
        std::cout << log_prefix() << "[return_linearized_alpha_pdf] beta=" << beta
                  << " n: " << pre_size << " -> " << a_vals.size()
                  << " (" << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << " ms)" << std::endl;
    }
    return std::make_pair(a_vals, alpha_pdf);
}

std::pair<std::vector<double>, std::vector<double>> DistData::return_viable_linearized_alpha_pdf(double const& inc_energy, double const& beta){
    LogIndent _li;
    std::vector<double> a_vals = get_viable_alphas__(inc_energy, beta);
    auto [alpha_pdf, truthy] = get_alpha_line__(a_vals, beta);
    auto get_new_alpha_pdf_val = [this, beta](double alpha) {return return_arbitrary_TSL_val(alpha, beta).first;};
    size_t pre_size = a_vals.size();
    auto t0 = std::chrono::steady_clock::now();
    linearize(a_vals, alpha_pdf, get_new_alpha_pdf_val, 1e-15, 1e-3, ToleranceCondition::Adaptive);
    auto t1 = std::chrono::steady_clock::now();
    if (!silence){
        std::cout << log_prefix() << "[return_viable_linearized_alpha_pdf] E=" << inc_energy << " beta=" << beta
                  << " n: " << pre_size << " -> " << a_vals.size()
                  << " (" << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << " ms)" << std::endl;
    }
    return std::make_pair(a_vals, alpha_pdf);
}

void DistData::get_alpha_sampling_dists__(){
    LogIndent _li;
    /// NOTE: Beta grid for storage is determined at initialization and is set to the leapr input grid and is not used here
    /// NOTE: Alpha distributions are symmetric about b=0 so only do the positive half
    /// NOTE: Ensure that you are using the scaled half betas for calculation
    calculation_alpha_vals.reserve(beta_grid.size());
    calculation_alpha_cdfs.reserve(beta_grid.size());
    alpha_vals.reserve(beta_grid.size());
    int beta_idx = 0;
    for (auto beta: calculation_half_betas){
        if (!silence){
            std::cout << log_prefix() << "[get_alpha_sampling_dists] beta " << beta_idx+1 << "/" << calculation_half_betas.size()
                      << " = " << beta << std::endl;
        }
        auto [vals, pdf] = return_linearized_alpha_pdf(beta);
        calculation_alpha_vals.push_back(vals);
        calculation_alpha_cdfs.push_back(pdf_to_cdf(vals, pdf));
        alpha_vals.push_back(fit_cdf(vals, calculation_alpha_cdfs.back(), initialization_alpha_cdf_grid));
        alpha_cdf_grid.push_back(initialization_alpha_cdf_grid);
        beta_idx++;
    }
}

void DistData::linearize_alpha_sampling_dists__(){
    LogIndent _li;
    for (int i = 0; i < calculation_alpha_cdfs.size(); i++){
        // Lambda function to get a new beta value that corresponds to a desired cdf value
       auto get_new_cdf_point = [&](double desired_cdf_val){
            int cdf_insert = std::lower_bound(calculation_alpha_cdfs[i].begin()+1, calculation_alpha_cdfs[i].end(), desired_cdf_val) - calculation_alpha_cdfs[i].begin();
            return ENDF_interp(calculation_alpha_cdfs[i][cdf_insert],
                               calculation_alpha_cdfs[i][cdf_insert+1],
                               calculation_alpha_vals[i][cdf_insert],
                               calculation_alpha_vals[i][cdf_insert+1],
                               desired_cdf_val);
        };
        size_t pre_size = alpha_vals[i].size();
        auto t0 = std::chrono::steady_clock::now();
        linearize(alpha_cdf_grid[i], alpha_vals[i], get_new_cdf_point, 1e-15, 1e-3, ToleranceCondition::Adaptive);
        auto t1 = std::chrono::steady_clock::now();
        if (!silence){
            std::cout << log_prefix() << "[linearize_alpha_sampling_dists] i=" << i << "/" << calculation_alpha_cdfs.size()
                      << " n: " << pre_size << " -> " << alpha_vals[i].size()
                      << " (" << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << " ms)" << std::endl;
        }
    }
}

void DistData::calculate_sampling_dists(){
    LogIndent _li;
    auto total_t0 = std::chrono::steady_clock::now();

    if (!silence){std::cout << log_prefix() << "[calculate_sampling_dists] Computing cross sections..." << std::endl;}
    auto xs_t0 = std::chrono::steady_clock::now();
    auto [ene, xs] = return_final_ii_xs();
    incident_energy_grid = ene;
    cross_section = xs;
    if (!silence){
        std::cout << log_prefix() << "[calculate_sampling_dists] XS done: n_energies=" << incident_energy_grid.size()
                  << " (" << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - xs_t0).count() << " ms)" << std::endl;
    }

    if (!silence){std::cout << log_prefix() << "[calculate_sampling_dists] Computing beta sampling distributions..." << std::endl;}
    auto beta_t0 = std::chrono::steady_clock::now();
    get_beta_sampling_dists__();
    if (!silence){
        std::cout << log_prefix() << "[calculate_sampling_dists] Beta sampling done ("
                  << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - beta_t0).count() << " ms)" << std::endl;
    }

    if (!silence){std::cout << log_prefix() << "[calculate_sampling_dists] Computing alpha sampling distributions..." << std::endl;}
    auto alpha_t0 = std::chrono::steady_clock::now();
    get_alpha_sampling_dists__();
    if (!silence){
        std::cout << log_prefix() << "[calculate_sampling_dists] Alpha sampling done ("
                  << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - alpha_t0).count() << " ms)" << std::endl;
    }

    if (linearize_cdfs){
        if (!silence){std::cout << log_prefix() << "Linearizing sampling distributions" << std::endl;}
        auto lin_beta_t0 = std::chrono::steady_clock::now();
        linearize_beta_sampling_dists__();
        if (!silence){
            std::cout << log_prefix() << "[calculate_sampling_dists] Beta CDF linearization done ("
                      << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - lin_beta_t0).count() << " ms)" << std::endl;
        }
        auto lin_alpha_t0 = std::chrono::steady_clock::now();
        linearize_alpha_sampling_dists__();
        if (!silence){
            std::cout << log_prefix() << "[calculate_sampling_dists] Alpha CDF linearization done ("
                      << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - lin_alpha_t0).count() << " ms)" << std::endl;
        }
    }

    if (!silence){
        std::cout << log_prefix() << "[calculate_sampling_dists] Total time: "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - total_t0).count() << " ms" << std::endl;
    }
}

double DistData::beta_to_outgoing_energy(double const& inc_energy, double const& beta){
    LogIndent _li;
    double result = inc_energy + tsl_data.temp*boltz*beta;
    if (!silence){
        std::cout << log_prefix() << "[beta_to_outgoing_energy] E_in=" << inc_energy << " beta=" << beta << " E_out=" << result << std::endl;
    }
    return result;
}

std::vector<double> DistData::betas_to_outgoing_energies(double const& inc_energy, std::vector<double> const& betas){
    LogIndent _li;
    if (!silence){
        std::cout << log_prefix() << "[betas_to_outgoing_energies] E_in=" << inc_energy << " n_betas=" << betas.size() << std::endl;
    }
    std::vector<double> outs;
    outs.reserve(betas.size());
    for (double beta:betas){
        outs.push_back(beta_to_outgoing_energy(inc_energy, beta));
    }
    return outs;
}

double DistData::return_beta(double const& inc_energy, double const& out_energy){
    LogIndent _li;
    double result = (out_energy - inc_energy)/(boltz*tsl_data.temp);
    if (!silence){
        std::cout << log_prefix() << "[return_beta] E_in=" << inc_energy << " E_out=" << out_energy << " beta=" << result << std::endl;
    }
    return result;
}

double DistData::alpha_to_scattering_angle(double const& inc_energy, double const& out_energy, double const& alpha){
    LogIndent _li;
    double t1 = inc_energy + out_energy;
    double t2 = alpha*tsl_data.a0*boltz*tsl_data.temp;
    double t3 = 2*sqrt(inc_energy*out_energy);
    double result = (t1 - t2)/t3;
    if (!silence){
        std::cout << log_prefix() << "[alpha_to_scattering_angle] E_in=" << inc_energy << " E_out=" << out_energy
                  << " alpha=" << alpha << " cos_theta=" << result << std::endl;
    }
    return result;
}

std::vector<double> DistData::alphas_to_scatting_angles(double const& inc_energy, double const& out_energy, std::vector<double> const& alphas){
    LogIndent _li;
    if (!silence){
        std::cout << log_prefix() << "[alphas_to_scatting_angles] E_in=" << inc_energy << " E_out=" << out_energy
                  << " n_alphas=" << alphas.size() << std::endl;
    }
    std::vector<double> angles;
    angles.reserve(alphas.size());
    for (double alpha:alphas){
        angles.push_back(alpha_to_scattering_angle(inc_energy, out_energy, alpha));
    }
    return angles;
}
