#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "rap_mb.cc"

namespace py = pybind11;


std::unordered_map<std::string, double> convert_dict_to_map(py::dict dictionary)
{
    std::unordered_map<std::string, double> result;
    for (std::pair<py::handle, py::handle> item : dictionary)
    {
        auto key = item.first.cast<std::string>();
        auto value = item.second.cast<double>();
        //cout << key << " : " << value;
        result[key] = value;
    }
    return result;
}

double calculate_ensemble_average(py::dict dictionary){
    // Create an instance of the rap_mb class
    std::unordered_map<std::string, double> input_dict = convert_dict_to_map(dictionary);
    rap_mb rap(input_dict);
    // Call the monte_carlo function to get the final population in the excited state
    double n1_avg = rap.monte_carlo();
    // Print the result
    return n1_avg;
}

PYBIND11_MODULE(_rap_mb, m) {
    m.doc() = "Rap MB module";
    m.def("calculate_ensemble_average", &calculate_ensemble_average, "A function that calculates the final population in the excited state");
};
