/**
 * rap_mb.cc
 * 
 * 
*/
#include <iostream>
#include <vector>
#include <unordered_map>
#include <random>
#include <math.h>
#include <omp.h>
#include <arkode/arkode_erkstep.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>

// Some Fundamental Constants
const double pi = M_PI;
const double c = 299792458.0;               // speed of light in m/s
const double hbar = 1.054571817E-34;        // reduced Planck's constant in Js
const double epsilon_0 = 8.8541878128E-12; // vacuum electric permittivity in F/m 

// a data structure to hold bloch params
struct bloch_params{
    double E0; 
    double w_x, w_y;
    double roc_zx, roc_zy;
    double omega_10;
    double mu_10;
    double vx, vy, vz;
};


class rap_mb {
    // A C++ class for Rapid Adiabatic Passage in a molecular beam interacting with CW laser.

public:
    rap_mb(const std::unordered_map<std::string, double> &input_dict): input_dict(input_dict){}
    const std::unordered_map<std::string, double> &input_dict;
    // field params
    const double wl_nm = input_dict.at("wavelength_nm");
    const double r = input_dict.at("beam_radius"); 
    const double f = input_dict.at("focal_length");
    const double P = input_dict.at("power");
    const double zi = input_dict.at("zi"); 
    // beam params
    const double v = input_dict.at("velocity");
    const double mu_10 = input_dict.at("dipole_mom");
    const double dv = input_dict.at("delta_v");
    const double divg = input_dict.at("angular_divg");
    // simulations params
    const double U0 = input_dict.at("U0");
    const double V0 = input_dict.at("V0");
    const double W0 = input_dict.at("W0");
    const double Ns = input_dict.at("N_samples");
    const double rel_tol = input_dict.at("rel_tol");    // relative error
    const double abs_tol = input_dict.at("abs_tol");    // absolute error
    const double dt_max = input_dict.at("dt_max");      // maximum size of dt    
    const int max_tsteps = input_dict.at("max_tsteps"); // maximum number of time steps
    // Calculating parameter for optical-bloch equations.
    double omega_10 = 2.0*pi*c/(wl_nm*1e-9);  
    double w_x0 = f*wl_nm*1.0e-7/(r*pi);
    double w_y0 = r;
    double zr_x = pi*pow(w_x0, 2.0)/(wl_nm*1.0e-7);
    double zr_y = pi*pow(w_y0, 2.0)/(wl_nm*1.0e-7);
    double w_x = w_x0*sqrt(1.0+pow((zi/zr_x), 2.0));
    double w_y = r;//w_y0*sqrt(1.0+pow((zi/zr_y), 2.0));
    double roc_zx = zi + pow(zr_x, 2.0)/zi;
    double roc_zy = zi + pow(zr_y, 2.0)/zi;
    double E0 = sqrt(4.0*P/(epsilon_0*c*pi*w_x*w_y*1e-4));

    // random sampling of velocities
    std::vector<double> sample_velocities(){
        // Random number generator
        std::random_device rd;
        std::mt19937 gen_rnd(rd());
        double v_fwhm = v*dv/100.0;
        std::normal_distribution<double> vel_dist(v, v_fwhm/2.0);
        std::normal_distribution<double> ang_dist(0.0, divg);
        double v0 = vel_dist(gen_rnd);
        double alpha = ang_dist(gen_rnd)*pi/180.0;
        double beta = ang_dist(gen_rnd)*pi/180.0;
        double gamma = ang_dist(gen_rnd)*pi/180.0;
        double vx = v0*cos(alpha);
        double vy = v0*sin(beta);
        double vz = v0*sin(gamma);
        std::vector<double> vel_xyz = {vx, vy, vz};
        return vel_xyz;
    };

    // A function for the RHS of Optical-Bloch equations
    static int bloch_equation(sunrealtype t, N_Vector psi, N_Vector dpsi_dt, void *params) {
        double *psi_data = NV_DATA_S(psi);
        double *dpsi_dt_data = NV_DATA_S(dpsi_dt);
        // extract params
        bloch_params *params_ptr = (bloch_params *) params;
        // Calculating detuning and rabi frequency
        double E = params_ptr->E0*exp(-pow(params_ptr->vx*t/(params_ptr->w_x*1.0e-2), 2.0) 
                                      - pow(params_ptr->vy*t/(params_ptr->w_y*1.0e-2), 2.0));
        double omega = E*params_ptr->mu_10/hbar;
        double delta = params_ptr->vz*params_ptr->omega_10/c -  pow(params_ptr->vx,2)*params_ptr->omega_10*t/((params_ptr->roc_zx*1.0e-2)*c);
        dpsi_dt_data[0] = -delta*psi_data[1];
        dpsi_dt_data[1] = delta*psi_data[0] + omega*psi_data[2];
        dpsi_dt_data[2] = -omega*psi_data[1];
        return 0;
    };

    // A function that solves the IVP and returns W
    double ivp_solve(){
        SUNContext sunctx;
        SUNContext_Create(SUN_COMM_NULL, &sunctx);
        //initiating and setting bloch params
        bloch_params params;
        params.E0 = E0; 
        params.w_x = w_x;
        params.w_y = w_y;
        params.roc_zx = roc_zx;
        params.roc_zy = roc_zy;
        params.omega_10 = omega_10;
        params.mu_10 = mu_10;
        std::vector<double> vel_xyz = sample_velocities();
        params.vx = vel_xyz[0];
        params.vy = vel_xyz[1];
        params.vz = vel_xyz[2];
        // calculating time limits
        double T = 2.0*(w_x*1.0e-2)/vel_xyz[0];
        double ti = -1.50*T;
        double tf = 1.50*T;
        // initiating the state vector and set psi0
        N_Vector psi0 = N_VNew_Serial(3, sunctx); 
        NV_Ith_S(psi0, 0) = U0;
        NV_Ith_S(psi0, 1) = V0;
        NV_Ith_S(psi0, 2) = W0;
        void* arkode_mem = ERKStepCreate(bloch_equation, ti, psi0, sunctx);
        // Passing bloch params
        ERKStepSetUserData(arkode_mem, &params); 
        // setting tolerances
        ERKStepSStolerances(arkode_mem, rel_tol, abs_tol);
        // setting max stepsize
        ERKStepSetMaxStep(arkode_mem, dt_max);
        //setting max number of steps
        ERKStepSetMaxNumSteps(arkode_mem, max_tsteps);
        // Setting Butcher Table of the method
        ERKStepSetTableNum(arkode_mem, ARKODE_DORMAND_PRINCE_7_4_5);
        // time iteration variable for the solver
        double t = 0.0;
        ERKStepEvolve(arkode_mem, tf, psi0, &t, ARK_NORMAL);
        // Extracting psi at t=tf, to get final W
        sunrealtype *psi = N_VGetArrayPointer(psi0);
        double W = psi[2];
        // cleaning up memory
        N_VDestroy(psi0);
        ERKStepFree(&arkode_mem);
        return W;
    };

    double monte_carlo(){
        double W_avg = 0.0;
        // parallelizing the Monte Carlo runs
        const int nthr = 4; //omp_get_num_threads();
        omp_set_num_threads(nthr);
        #pragma omp parallel for reduction(+:W_avg)
        for (int i=0; i<int(Ns); i++) {
            W_avg += ivp_solve();
        };
        W_avg = W_avg/Ns;
        // calculate final population in the excited state.
        double n1_avg = (W_avg+1.0)/2.0;
        return n1_avg;
    };
};

