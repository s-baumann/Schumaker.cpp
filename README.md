# SchumakerSpline

C++ code to create a shape preserving spline. This is guaranteed to be monotonic and concave or convex if the data is monotonic and concave or convex. It does not use any optimisation and is therefore quick and smoothly converges to a fixed point in economic dynamics problems including value function iteration. This package has the same functionality as the R package called [schumaker](https://cran.r-project.org/web/packages/schumaker/index.html) and the Julia package called [SchumakerSpline.jl](https://github.com/s-baumann/SchumakerSpline.jl).

Note that the main.cpp file only demonstrates the use of the spline and has no generalisable functionality. The two important files that you might want to include in a project are Schumaker.hpp and Schumaker.cpp.

Copyright Stuart Baumann 2024 and under MIT license (https://en.wikipedia.org/wiki/MIT_License).