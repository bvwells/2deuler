# 2deuler
[![Build Status](https://travis-ci.org/bvwells/2deuler.svg?branch=master)](https://travis-ci.org/bvwells/2deuler)

Two dimensional Compressible Euler Equations solved with moving mesh approach described in the PhD thesis

*A moving mesh finite element method for the numerical solution of partial differential equations and systems.*

which can be found [here][1].

The Compressible Euler equations are described by the non-linear system of partial differential equations

```
ρ_t + (ρu)_x + (ρv)_y = 0
(ρu)_t + (ρu^2 + P)_x + (ρuv)_y = 0
(ρv)_t + (ρuv)_x + (ρv^2+P)_y = 0
E_t + (u(E+P))_x + (v(E+P))_y = 0
```

where ```ρ```, ```u```, ```u```, ```P``` and ```E``` are the density, ```x``` and ```y``` velocity components, pressure and energy of the gas being modelled. The equation of state for an ideal gas is used which is given by

```
E=P/(1-γ) + 0.5ρ(u^2 + v^2)
```
where ```γ``` is the ratio of specific heats for the gas.

## Numerical Solution

The two-dimensional Euler Equations are solved using a moving mesh 
method which uses the monitor function ```M=u(x,t)``` in the moving mesh 
equations for the mesh velocity. The mesh is advanced forwards in time 
using a forward Euler time-stepping scheme. The solution to the Euler
equations are obtained by solving the Arbitrary Lagrangian Eulerian (ALE)
formation of the equations with a finite volume method. All the moving mesh
equations are solved using linear finite elements.

## Building and Developing

Developing locally requires Docker for Windows. Run the command

```
docker build -t 2deuler .
```

to build the docker image which pulls the gcc base image containing gfortran and maps the source code into the container.

Then run image with the command:

```
docker run -i -t -v /f/git/src/github.com/bvwells/2deuler:/app 2deuler
```

This command maps the local workspace into the running image so any changes made in the running image will be reflected on the local workspace.

Within the running image generate the make files for the release version by running the command:

```
cmake .
```

To build the debug version of the code run the command:

```
cmake -DCMAKE_BUILD_TYPE=Debug
```

Build the executable by running the command:

```
make
```

## Running

The program takes the file [variables.data](./variables.data) as input to the simulation. The program can be run from the base of the repo with the command:

```
./bin/2deuler.exe
```

The program outputs the mesh and solution over time into the files ```SolutionXXX.m```. The variables for the solution are written to the file ```variables.m```.

## Plotting Solution

The output from the simulation can be plotted in [Octave](https://www.gnu.org/software/octave/) by running the plotting file
[plot_solution.m](./plot_solution.m) in the root of the repo.

[1]: http://www.reading.ac.uk/nmsruntime/saveasdialog.aspx?lID=24080&sID=90294
