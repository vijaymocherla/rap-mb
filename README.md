# RAP-MB

This repo contains a program for simulating Rapid Adiabatic Passage (RAP) in Molecular Beam (MB) experiments. The interaction between the molecular beam and a CW laser is modeled using the model of a two-level system (TLS) intereacting with radiation and, its time evolution is described by the optical-bloch equations. Further, the Monte Carlo method is used to sample molecular velocities, to calculate an ensemble average of the excited state population over many different trajectories of molecules passing through the laser field.   

To optimize computational efficiency and ease of use, `rap-mb` is written in C++ with a python interface. For better performance, we use the `ARKODE` module within `SUNDIALS` to numerically solve the Optical-Bloch equations. Further, the Monte-Carlo trajectories are parallelised using `OpenMP`. 


In the following jupyter notebooks we illustrate the usuage of this program with a few examples.
- [jcp-methane.ipynb](notebooks/jcp-methane.ipynb)
- [carbon-monoxide.ipynb](notebooks/carbon-monoxide.ipynb)


**References**
- Chadwick, Helen, P. Morten Hundt, Maarten E. van Reijzen, Bruce L. Yoder, and Rainer D. Beck. "Quantum state specific reactant preparation in a molecular beam by rapid adiabatic passage." The Journal of chemical physics 140, no. 3 (2014).
- Reynolds, Daniel R., David J. Gardner, Carol S. Woodward, and Rujeko Chinomona. "ARKODE: A flexible IVP solver infrastructure for one-step methods." ACM Transactions on Mathematical Software 49, no. 2 (2023): 1-26.

## Installation
To use the `rap-mb` program, clone this repository and the compile the code. 
> Note: Make sure to install `SUNDIALS` before you do this. 

```
git clone github.com/vijaymocherla/rap-mb
cd rap_mb
bash compile.sh
``` 

## SUNDIALS

SUNDIALS is a SUite of Nonlinear and DIfferential/ALgebraic equation Solvers developed and maintained by LLNL. 

Here's a breif walkthrough on building SUNDIALS from source. A more detailed guide can be found [here](https://sundials.readthedocs.io/en/latest/sundials/Install_link.html#building-and-installing-with-cmake). 

### Dependencies

- Make sure you have install `intel-oneapi`'s base  and HPC toolkits. This should basically install the intel compilers (icx, icpx and ifx) and, with them the libraries math kernel libraries (MKL, LAPACK, BLAS etc.) and, parallel processing libraries (MPI, openMP, SYCL etc.). Make sure to intialise (activate) the intel oneapi environment. 

```
source /opt/intel/oneapi/setvars.sh
```

> Note: One could also install SUNDIALS using the GNU equivalents: gcc, g++ and gfortran, along with openBLAS, openMP, mpich.

- cmake and cmake-gui are needed for this installation. This can be done as follows:
    ```
    sudo apt install cmake
    sudo apt install cmake-curses-gui
    ```
### Building from Source

- Obtaining and building SUNDIALS: 

Git clone the SUNDIALS repo. (Or alternatively, it can be downloaded from https://computing.llnl.gov/projects/sundials/sundials-software)

```
git clone https://github.com/LLNL/sundials
```
- Go the `sundials` directory and build the package. 
    - First, set the environment variables for intel compiler (use `icx`, `icpx`, `ifx`) for `cmake` and generate the initial configuration. 
        ```
        cd sundials
        mkdir BUILDIR
        cd BUILDIR
        export CC=icx
        export CXX=icpx
        export FC=ifx
        cmake ..
        ```
    - Then following that, use the cmake-gui to change and adjust the configuration. Check the following options (set them to `ON`): 
        ```
        BUILD_FORTRAN_MODULE_INTERFACE
        ENABLE_ONEMKL
        ENABLE_LAPACK
        ENABLE_MPI
        ENABLE_OPENMP
        ```
    - Use the configure option [c] to build this configuration. You may run into an error for MKL Configuration and would have to manually provide the CMake config files for MKL. Usually, these can be found in the folder `/opt/intel/oneapi/mkl/latest/lib/cmake/mkl`, just look in the folder where you have oneapi installed. Use the exit options [e] and, you can set the `MKL_DIR` option to this folder in the cmake-gui window.
    - Now, use [c] to start building the configuration again. Now, you would have to choose the math libraries to link to the SUNDIALS installation. You should check the following options (set them to `ON`):
    ```
    BLACS
    BLAS95
    CDFT
    LAPACK95
    SCALAPACK
    ``` 
    - Configure these options with [c] and generate this configuration with [g] in `cmake-gui` and quite the window.

    - Now, make the package and install it.
        ```
        make -j 4
        make install 
        ```

