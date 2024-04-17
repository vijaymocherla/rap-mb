# include "rap_mb.cc"

int main() {
    // Define the input dictionary
    std::unordered_map<std::string, double> input_dict = {
        {"wavelength_nm", 4666.0},
        {"beam_radius", 0.20},
        {"focal_length", 15.0},
        {"power", 0.020},
        {"zi", 10.0},
        {"velocity", 1200.0},
        {"dipole_mom", 4.77e-31},
        {"delta_v", 0.1},
        {"angular_divg", 0.01},
        {"N_samples", 1000.0},
        {"N_tstep", 100.0},
        {"U0", 0.0},
        {"V0", 0.0},
        {"W0", -1.0}
    };
    // Create an instance of the rap_mb class
    std::cout << "Input dictionary created" << std::endl;
    rap_mb rap(input_dict);
    std::cout << "Initiated rap" << std::endl;
    // Call the monte_carlo function to get the final population in the excited state
    double n1_avg = rap.monte_carlo();
    // Print the result
    std::cout << "Final population in the excited state: " << n1_avg << std::endl;

    return 0;
}
