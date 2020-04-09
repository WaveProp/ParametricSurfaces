# ParametricSurfaces.jl

*A simple package to generate some basic parametric shapes for testing Nystrom methods* 

![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)<!--
![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-retired-orange.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-archived-red.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-dormant-blue.svg) -->
[![Build Status](https://travis-ci.com/maltezfaria/ParametricSurfaces.jl.svg?branch=master)](https://travis-ci.com/maltezfaria/ParametricSurfaces.jl)
[![codecov.io](http://codecov.io/github/maltezfaria/ParametricSurfaces.jl/coverage.svg?branch=master)](http://codecov.io/github/maltezfaria/ParametricSurfaces.jl?branch=master)

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
