This repository contains scripts for benchmarking dynamical systems on networks. Usage Instructions are found below. Details on our methodology are reported in our [paper]((https://arxiv.org/abs/2012.12696) on [NetworkDynamics.jl](https://github.com/PIK-ICoN/NetworkDynamics.jl). Below you see plots of competing runtimes simulating a Kuramoto network of 10, 100 and 1000 oscillators.

![Work precision diagram of a network with 10 Kuramoto oscillators.](https://github.com/PIK-ICoN/NetworkDynamicsBenchmarks/blob/main/utils/plotting/WPD10.png?raw=true)
![Work precision diagram of a network with 100 Kuramoto oscillators.](https://github.com/PIK-ICoN/NetworkDynamicsBenchmarks/blob/main/utils/plotting/WPD100.png?raw=true)
![Work precision diagram of a network with 1000 Kuramoto oscillators.](https://github.com/PIK-ICoN/NetworkDynamicsBenchmarks/blob/main/utils/plotting/WPD1000.png?raw=true)

In the concrete implementations we tried to strike a balance between usability (writing code should feel "simple") and performance. We compare two implementations in Python (SciPy and JiTCODE, the later JIT-compiles a C function for the ODE), two in Julia (NetworkDynamics.jl and incidence matrix multiplication with SparseArrays.jl) and one in Fortran (similar approach as for SparseArrays.jl).

The JITCODE (and SciPy) version is derived from this publication and does not make use of the symmetry of the sine function 𝑠𝑖𝑛(𝑦−𝑥)=−𝑠𝑖𝑛(𝑥−𝑦), which means that the right hand side function evaluates roughly twice as many floating point operations.

All programs employ Dormand-Prince's 5th order Runge-Kutta method. For Julia we use DifferentialEquations.jl to solve the ODE.

When startup and JIT-compile time is taking into account we see a different picture. The following plot shows the time until the first trajectory is solved.

![System size vs. total time](https://github.com/PIK-ICoN/NetworkDynamicsBenchmarks/blob/main/utils/plotting/size_vs_jit.png?raw=true)

Fortran isn't shown in this plot since that program doesn't have a JIT compile stage. In fact the whole Fortan benchmark (300 integrations) finishes roughly in the time it takes to startup a Julia session and import all required libraries (Julia 1.5.1).  If you know how to improve the benchmarks or want to add your favourite network dynamics solution in yet another language, we are happy about contributions. Just open a pull request.


# Instructions

## High level Requirements

 - sqlite
 - make
 - python3
 - Julia 1.5

## Setup python environment
We suggest a virtual env for python

```
python3 -m venv .env
```

activate and install the requirements

```
source .env/bin/activate
pip install -r src/python/requirements.pip
```

## Setup Julia enviroment

The julia env should be activated automatically by the scripts.

## Compile fortran program

If you would like to run the fortran benchmarks, you must have a Fortran
2003 compiler. In the example we use GFortran. If you do not wish to run
the fortran benchmark, bypass this section.

The fortran program uses FPM, the fortran package manager. Install from
https://github.com/fortran-lang/fpm

Using this tool, you can build the fortran program with fpm

```
cd src/fortran
fpm build --flag -ldl --flag -lpthread
```

now, set KURABENCH_FEXE to the newly built binary

```
cd ../../
export
KURABENCH_FEXE=`pwd`/src/fortran/build/<compiler>_<uuid>/app/FKuraBenchmark
```

where compiler and uuid are the details of your compilation process
resulting from the fpm build command.

## Prepare Database

In these benchmarks, we use an SQLite database to store the setup of
each system (node parameters, edge parameters, system parameters,
connectivity) and also the experiment log and trajectory outputs.

You may place your database at a location of your choosing. Each program
uses the environmental variable KURABENCH_DB to identify the location of
the database. In this case, we will place the database in the `data`
directory.

```
export KURABENCH_DB=`pwd`/data/kurabench.db
make database
```

## Generate systems

Be sure that your python virtual env is activated

```
make systems
```

## Generate reference trajectories

The "ground-truth" trajectories are generated using the (implicit) radau
solver at low tolerance settings. Julia is used to run this generation
step. Some timing information will be printed to stdout during the
program execution. Compute the trajectories with:

```
make reference
```

## Generate benchmarks

Now, generate any benchmarks of interest by running

```
make <X>
```

where `<X>` is one of

   networkdynamics
   diffeq
   jitcode
   scipy
   fortran

The resulting experiment summary is stored in the database in the
experiments table. The trajectory of each state is stored in the
`states` table, indexed by the experiment id found in the `experiments`
table.

The `experiments` table includes the timing and `err_v_ref`, the average
deviation (across all states) of the experiment trajectory and the
reference (radau-computed) trajectory for that system and initial
condition.
