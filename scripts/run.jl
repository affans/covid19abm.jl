using Distributed
using Random
using DelimitedFiles
using Distributions
using Base.Filesystem
using DataFrames
using CSV
using ClusterManagers
using Query
using Statistics
using UnicodePlots
using Dates
using ProgressMeter

using covid19abm
const cv=covid19abm

addprocs(SlurmManager(512), N=16) 

@everywhere using covid19abm
@everywhere using ProgressMeter

function run(beta) 
    nsims = 1000
    p = cv.ModelParameters()
    
    cd = pmap(1:nsims) do x 
        vals = main(x, p)        
    end    

    for i = 1:nsims
        dts = cd[i]
        insertcols!(dts, 1, :sim => i)    
        insertcols!(dts, 1, :time => 1:500)
    end

    d = vcat([cd[i] for i = 1:length(cd)]...)
    ## create yearly average
    dd = d |> @groupby(_.time) |>
        @map({time=key(_), cnt=length(_),
              sus=mean(_.sus), lat=mean(_.lat), inf=mean(_.inf), iso=mean(_.iso), infiso=mean(_.infiso), 
              hos=mean(_.hos), icu=mean(_.icu), ded=mean(_.ded), rec=mean(_.rec)}) |> DataFrame
    return dd
end


function calibrate(beta)
    p = cv.ModelParameters()
    p.β = beta 
    p.calibration = true
    nsims = 1000
    vals = zeros(Int64, nsims)
    println("calibrating with beta: $(p.β)")
    for i = 1:nsims
        a = main(i, p)
        vals[i] = a.lat[end] - 1 ## minus becuase the initial latent guy ends up in latent by the end cuz of swapupdate()
    end
    println("mean R0: $(mean(vals)) with std: $(std(vals))")
    return vals
end

function calibrate_robustness(beta)
    #once a beta is found based on nsims simulations, 
    # see how robust it is. run calibration with same beta 100 times 
    # to see the variation in R0 produced. 
    # means = zeros(Float64, 100)
    # for i = 1:100 
    #     vals = calibrate(beta)
    #     means[i] = mean(vals)
    # end
    cd = pmap(1:100) do x 
        vals = calibrate(beta)
        return mean(vals)
    end
    return cd
end