# ParametricSurfaces.jl

*A simple package to generate some basic parametric shapes for testing Nystrom methods* 

![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)<!--
![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-retired-orange.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-archived-red.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-dormant-blue.svg) -->

![CI](https://github.com/IntegralEquations/ParametricSurfaces/workflows/CI/badge.svg?branch=master)
[![codecov.io](http://codecov.io/github/IntegralEquations/ParametricSurfaces.jl/coverage.svg?branch=master)](http://codecov.io/github/maltezfaria/ParametricSurfaces.jl?branch=master)

## Installation
Install from the Pkg REPL:
```
pkg> add https://github.com/IntegralEquations/ParametricSurfaces
```

## Usage

The most basic use case is to create some simple shape, and then compute a quadrature of the shape:
```julia
    using ParametricSurfaces
    # create a geometry 
    geo = Bean(center=(0,0,0),paxis=(1,1,1))
    # generate a tensor quadrature with 10 points per patch per direction
    quad = TensorQuadrature((10,10),geo)
    # plot 
    using Plots
    pyplot()
    plot(quad)
```
You can also refine globally by splitting:
```julia
    refine!(geo)
    quad = TensorQuadrature((10,10),geo)
    plot(quad)
```
or split a given part:
```julia
    geo = Bean()
    refine!(geo.parts[1])
    quad = TensorQuadrature((10,10),geo)
    plot(quad)
```

See the [tests](./test/runtests.jl) for more.

## Integration with GMSH

There is currently some experimental integration with GMSH
```julia
    gmsh.initialize()
    gmsh.open("./meshes/halfmodel.stp")
    body = GmshParametricBody(3,1,skip=[41])
    # refine!(body)
    quad = TensorQuadrature((10,10),body)
    pyplot()
    plot(quad)
    gmsh.finalize()
```
If everything went right you should a mesh of the generated quadrature of the halfmodel.

:warning: There appears to be a problem with surfaces of type `Plane` in gmsh. Need to investigate this further.
