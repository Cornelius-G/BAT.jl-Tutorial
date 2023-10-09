# <img style="height:2.8em;" alt="BAT.jl" src="https://bat.github.io/BAT.jl/dev/assets/logo.svg"/> &nbsp;  Tutorial for the PUNCH Young Academy  

If you have access to the PUNCH jupyterhub service, you can log in to https://jupyterhub.uni-muenster.de and start the Julia-1.9.3 kernel.

If you want to install on your machine follow the steps below:
1. Download the Julia 1.9.3 release for your system from https://julialang.org/downloads/
2. Unpack/Install Julia according to the platform specific instructions: https://julialang.org/downloads/platform/
3. Launch Julia 
4. Copy and paste the following snippet to install all required packages
   ```julia
   using Pkg
   pkg"add BAT, Distributions, IntervalSets, ValueShapes, Plots, ArraysOfArrays, StatsBase, LinearAlgebra, DensityInterface, Optim, EmpiricalDistributions"
   ```
You should now have everything installed and should be able to run the tutorial notebook.
