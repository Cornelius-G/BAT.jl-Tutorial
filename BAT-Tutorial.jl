#   BAT.jl Tutorial - Poisson Counting Experiment
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

#   Note: This tutorial is written in the style of a cloze test. Parts that have
#   to be filled out by you are marked #=Please Fill=#.
# 
#   For additional help, you can take a look at the BAT.jl tutorial
#   (https://bat.github.io/BAT.jl/stable/tutorial/) or the BAT.jl documentation
#   (https://bat.github.io/BAT.jl/stable/).

#   To get started, we first need to import some packages we will use in this
#   tutorial:

using BAT
using Distributions 
using IntervalSets
using ValueShapes
using Plots
using ArraysOfArrays
using StatsBase 
using LinearAlgebra
using DensityInterface

#   Situation
#   ===========
# 
#   In this tutorial, we want to determine the properties of a radioactive
#   singal source in the presence of background from natural sources of
#   radioactivity. We assume to have one signal source S and only one source of
#   background B. All measurements are taken for the same duration, so we do not
#   distinguish between counts and rates.

#   Part 1: A first background-only measurement
#   =============================================
# 
#   We start by using our detector without the signal source installed. This
#   measurement yields a number of k_B=10 counts. Using this measurement we want
#   to gain information about the event rate of the natural radioactive
#   background.
# 
#   Task:
#   –––––––
# 
#   We want to perform a Bayesian analysis to estimate the event rate of the
#   natural background \lambda_b using a Poisson model.
# 
#   We start by defining the corresponding log-likelihood function, using the
#   function logpdf() and the type Poisson provided by the package
#   Distributions.jl
#   (https://juliastats.github.io/Distributions.jl/latest/univariate/).

# Number of observed background events
kb = 10

likelihood_B = let k = kb
    logfuncdensity(function (params) # this function is part of the package DensityInterface.jl and is used to define the log-likelihood
        return logpdf(Poisson(params.λb), k) # poisson log-likelihood
    end)
end;

#   Next, we define the Prior distribution using the distprod function. We use a
#   flat prior between 0 and 30:

prior_B = distprod(
    λb = 0..30., # This is a shorthand notation for Uniform(0, 30)
)

#   With the likelihood and prior, we can now define the PosteriorDensity() for
#   the background-only scenario:

posterior_B = PosteriorMeasure(likelihood_B, prior_B)

#   We can now explore the posterior distribution, for example by first
#   searching for the mode, i.e. the value of \lambda_b with the highest
#   probability density:

using Optim # we need to explicitly load the Optim package to use the bat_findmode function
mode_B = bat_findmode(posterior_B)

#   The value of the mode can be obtained as:

mode_B.result

#   In a next step, we want to actually sample the posterior distribution to
#   check other numerical properties and create plots. For this, we first
#   specify some settings for the sampling. We choose the MetropolisHastings()
#   algorithm and use 4 Markov chains with 10^5 steps each:

algorithm = MCMCSampling(mcalg = MetropolisHastings(), nsteps = 10^5, nchains = 4)

#   Using the function bat_sample we can then sample the posterior distribution:

samples_B = bat_sample(posterior_B, algorithm).result

#   Let us now look at some numerical properties of the posterior distribution
#   obtained from the samples:

bat_report(samples_B)

#   The individual numerical values shown above can be also obtained using the
#   functions mean, mode, std, ... on the samples object, e.g.:

mean(samples_B).λb

#   Let us now also take a look at the resulting disribution for the background
#   rate using the plot() function:

plot(samples_B, :λb)

#   We can also plot the prior distribution on top of the posterior to visualize
#   the knowledge update we gained by analyzing the data:

plot(samples_B, :λb)

samples_prior_B = rand(prior_B, 100_000).λb
plot!(samples_prior_B, st=:stephist, normalize=true, lw=2, label="Prior")

#   Additional Tasks:
#   –––––––––––––––––––
# 
#   Play around with the settings of the sampling algorithm. What happens if
#   you, e.g.
# 
#     •  change the number of chains ?
# 
#     •  change the number of steps per chain ?
# 
#     •  pick a different sampling algorithm from the List of algorithms
#        (https://bat.github.io/BAT.jl/stable/list_of_algorithms/#Sampling-algorithms)
#        ?

#   Part 2: Second background-only measurement
#   ============================================
# 
#   A second measurement of the natural background is perfomed (for the same
#   duration) and yields a number of k_B=8 counts. Therefore, we want to update
#   our estimation of the background rate using this new knowledge together with
#   the pervious results.
# 
#   Task:
#   –––––––
# 
#   We want to perform an analysis of the new measurement similar to the first
#   one. We want to use the posterior distribution of the previous background
#   measurement as a prior for this analysis. For this, we need to convert our
#   previous samples into a StatsBase histogram
#   (http://juliastats.github.io/StatsBase.jl/latest/empirical/#Histograms-1).

#   We build a StatsBase histogram containing the previous posterior
#   distribution:

samples_flat = samples_B.v.λb
weights = FrequencyWeights(samples_B.weight)

posterior_hist_B1 = fit(Histogram, samples_flat, weights, nbins = 400)

#   This histogram can now be used as a prior by converting it into a univariate
#   binned distribution (UvBinnedDist) using the package
#   EmpiricalDistributions.jl:

using EmpiricalDistributions

prior_B2 = distprod(
    λb = UvBinnedDist(posterior_hist_B1)
)

#   Now we define the new log-likelihood function and build the new posterior
#   distribution:

kb2 = 8

likelihood_B2 = let k = kb2
    logfuncdensity(function (params)
        return #=Please Fill=#
    end)
end;

posterior_B2 = #=Please Fill=#

#   And generate samples: (Note: By default, BAT.jl internally performes a phase
#   space transformation according to the prior distribution in order to make
#   the sampling more efficient. As this is currently not supported for priors
#   containing UvBinnedDists, we need to explicitly disable this transformation
#   using trafo=DoNotTransform().)

samples_B2 = #=Please Fill=#

#   Let us now visualize both the prior (i.e. the posterior of the first
#   analysis) and the updated posterior:

plot(samples_B2, :λb)
plot!(rand(prior_B2, 1_000_000).λb, st=:stephist, lw=2, normalize=true, label="Prior")

bat_report(samples_B2)

#   Part 3: Signal + Background Measurement
#   =========================================
# 
#   Now we include the radioactive source into the experimental setup. We repeat
#   the measurement and now observe k_{S+B}=12 counts. With this measurement and
#   our prior knowledge about the background we are now able to estimate the
#   signal rate \lambda_s.
# 
#   Task
#   ––––––
# 
#   Perform a third analysis using a poisson model with the combined singal +
#   background rate. Use the known information about the background from the
#   previous tasks as a prior for the backround and choose a suitable prior for
#   the signal.

#   Define the likelihood for the signal + background model:

# Number of observed events
kSB = 12

likelihood_SB = let k = kSB
    logfuncdensity(function (params)
        return #=Please Fill=# # poisson log-likelihood for b+s 
    end)
end;

#   Define the prior for both the signal and backgound parameters. Choose a flat
#   prior in the range (0, 30) for the signal rate. Remember to use the
#   knowledge from the previous tasks for the prior of the background rate.

posterior_hist_B2 = fit(Histogram, samples_B2.v.λb, FrequencyWeights(samples_B2.weight), nbins = 400)

prior_SB = NamedTupleDist( 
    #=Please Fill=#
)

#   Define the posterior for the signal + backround analysis:

posterior_SB = #=Please Fill=#

#   Generate samples for the signal + background model: (Note: Again use
#   trafo=DoNotTransform().)

samples_SB = #=Please Fill=#

#   To visualize an overview of the results for both prameters simply use
#   plot(samples_SB): You can also indicate the global mode adding the keyword
#   globalmode=true.

plot(samples_SB, globalmode=true)

#   Alsoe, print some statistics of the samples:

#=Please Fill=#

#   Part 3 b) - Combined analyses
#   ===============================
# 
#   Let us now perform the same analysis but this time using all three
#   measurments (first background measuerment, second background measurement,
#   and the signal + backround measurements) at the same time instead of
#   inserting them sequentially.

prior_SB_2 = NamedTupleDist( 
    #=Please Fill=#
)


likelihood_SB_2 = let k = [10, 8, 12]
    logfuncdensity(function (params)
        b1 = #=Please Fill=#
        b2 = #=Please Fill=#
        sb = #=Please Fill=#
        return b1 + b2 + sb  # poisson log-likelihood for b+s 
    end)
end;


posterior_SB_2 = PosteriorMeasure(likelihood_SB_2, prior_SB_2);

samples_SB_2 = #=Please Fill=#
bat_report(samples_SB_2)

#   4. Error propagation
#   ======================

#   Finally, we want to caluclate the cross section of the signal process. The
#   rate of measured events in the detector of a couting experiment can be
#   written as
# 
#   \frac{\mathrm d N}{\mathrm d t} = \epsilon \cdot σ \cdot L
# 
#   ,
# 
#   with the Luminosity L and the efficiency of the detector \epsilon. The
#   signal cross section is therefore given as:
# 
#   σ_S = \frac{λ_s}{\epsilon \cdot L}.
#   -------------------------------------
# 
#   For this experiment we assume a luminosity of L = 1.1 (neglecting units).
# 
#   As a final result we want to obtain either a measurement or an upper limit
#   on the signal cross section.

#   Task 4 a) Known efficiency with gaussian uncertainty
#   ––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 
#   The detector efficiency has been measured to be \epsilon = 0.1 \pm 0.02,
#   assuming the uncertainties to follow a normal distribution.
# 
#   We want to calculate the signal cross section σ_S using the equation above.
# 
#   Use the Distrtibutions.jl
#   (https://juliastats.github.io/Distributions.jl/latest/univariate/) package
#   and rand() to obain samples for \epsilon. In order to obtain the same number
#   of unweigthed samples of \lambda_S the function bat_sample(samples,
#   RandResampling(nsamples=nsamples)) can be used for resampling the posterior
#   samples.
# 
#   The function broadcast() or the .operator (e.g. a .+ b) might be useful for
#   element-wise operation when handeling the samples.

#   Define the luminosity and the efficiency:

nsamples = 100_000

L = 1.1
ϵ = rand(Normal(0.1, 0.02), nsamples)

#   Let us quickly plot the efficiency. Use a StatsBase histogram
#   (http://juliastats.github.io/StatsBase.jl/latest/empirical/#Histograms-1) to
#   visualize the distribution. (Hint: The plot recipes can also be used for the
#   StatsBase histograms)

hist_ϵ = fit(Histogram, ϵ, nbins=200) 

plot(hist_ϵ, st=:steps, normalize=true, lw=2, xlabel="\$\\epsilon\$", ylabel="\$p(\\epsilon)\$")

#   Now, let's get unweighted samples for the signal rate and calculate the
#   cross section distribution:

resampled_SB = #=Please Fill=#
λs = resampled_SB.v.λs

σS = (λs)./(ϵ*L);

#   Plot the distribution of the signal cross section. (Hint: use again a
#   histogram)

hist_σS = fit(Histogram, σS, nbins=300)
bat_hist = BAT.MarginalDist(UvBinnedDist(hist_σS))
plot(bat_hist, st=:smallest_intervals, xlabel="\$\\sigma_S\$", ylabel="\$p(\\sigma_S)\$")

#   Did we measure the signal process? Or can we only set an upper limit?

#   Task 4 b) Binomial analysis of calibration measurement to determine
#  efficiency
#   ––––––––––––––––––––––––––
# 
#   We now want to perform the same analysis as in task 4 a) but for the case
#   that the detector efficiency \epsilon is not yet known.
# 
#   Instead, \epsilon is to be determined using a calibration measurement with a
#   source for which the signal rate is known. Then it is possible to calculate
#   the efficiency from the expected number of counts and the number of counts
#   actually measured with the detector.
# 
#   For our example, the number of expected events is assumed to be
#   N_\text{expected} = 200. The detector measures only N_\text{measured} = 21
#   events.
# 
#   Task: Implement a binomial model using the Binomial(n,p) function of the
#   Distrtibutions.jl
#   (https://juliastats.github.io/Distributions.jl/latest/univariate/) package
#   and determine the distribution of the detector efficiency. Afterwards,
#   repeat the steps from 4 a) using the obtained distribution of the efficiency
#   to calculate the signal cross section .

#   Define the binomial likelihood:

n_expected = 200
n_measured = 21

likelihood_binomial = let n = n_expected, k = n_measured
    logfuncdensity(function (params)
        return #=Please Fill=#
    end)
end

#   Define the prior (flat) for the efficiency and define the posterior:

prior_binomial = distprod(
    #=Please Fill=#
)

posterior_binomial = PosteriorMeasure(likelihood_binomial, prior_binomial)

#   Generate the samples:

samples_binomial = #=Please Fill=#

#   Plot the distribution of the efficiency: (What do you observe when comparing
#   to the efficiency used in 4 a) ?)

plot(samples_binomial, :ε)

#   Print some statistics of the samples:

#=Please Fill=#

#   Calculate the cross section by sampling the same number of events from the
#   efficiency and from the signal samples: (Hint: proceed like in 4 a) for the
#   samples of \lambda_SB )

nsamples = 100_000
resampled_SB = bat_sample(samples_SB, RandResampling(nsamples=nsamples)).result
resampled_binomial = bat_sample(samples_binomial, RandResampling(nsamples=nsamples)).result

@assert length(resampled_SB) == length(resampled_binomial)

λ_SB = resampled_SB.v.λs
ϵ   = resampled_binomial.v.ε
σS = (λ_SB)./(ϵ*L)

#   Use a StatsBasehistogram to visualize the cross section distribution. From
#   the plot, determine the 95% upper limit on the cross section.

hist_σ = fit(Histogram, σS, nbins=200)
bat_hist = BAT.MarginalDist(UvBinnedDist(hist_σS))
plot(bat_hist, st=:smallest_intervals, xlabel="\$\\sigma_S\$", ylabel="\$p(\\sigma_S)\$")

#   Compare the distribution of the signal cross section to the distribution
#   from 4 a).