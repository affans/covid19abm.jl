using Distributed
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
    myp = ModelParameters()
    myp.β = beta
    @everywhere reset_params($myp)
    nsims = 1000      
    cd = pmap(1:nsims) do x 
        vals = main(x)        
    end    
    for i = 1:nsims
        dts = cd[i]
        insertcols!(dts, 1, :sim => i)    
        insertcols!(dts, 1, :time => 1:myp.modeltime)
    end

    d = vcat([cd[i] for i = 1:length(cd)]...)
    ## create yearly average
    dd = d |> @groupby(_.time) |>
        @map({time=key(_), cnt=length(_),
              sus=mean(_.sus), 
              lat=mean(_.lat), 
              mild=mean(_.mild), 
              miso=mean(_.miso), 
              inf=mean(_.inf), 
              iiso=mean(_.iiso), 
              hos=mean(_.hos), 
              icu=mean(_.icu), 
              rec=mean(_.rec), 
              ded=mean(_.ded), 
              lat_inc=mean(_.lat_inc), 
              mild_inc=mean(_.mild_inc), 
              miso_inc=mean(_.miso_inc), 
              inf_inc=mean(_.inf_inc),
              iiso_inc=mean(_.iiso_inc),
              hos_inc=mean(_.hos_inc),
              icu_inc=mean(_.icu_inc),
              rec_inc=mean(_.rec_inc),
              ded_inc=mean(_.rec_inc)
              }) |> DataFrame
    return dd
end

function calibrate(beta)
    myp = ModelParameters()
    myp.β = beta
    myp.h = 0 
    myp.f = 0
    myp.calibration = true
    @everywhere reset_params($myp)
    
    nsims = 1000
    vals = zeros(Int64, nsims)
    println("calibrating with beta: $(cv.p.β)")
    cd = pmap(1:nsims) do i 
        a = main(i)
        val = a.lat[end] - 1 ## minus becuase the initial latent guy ends up in latent by the end cuz of swapupdate()
        return val
    end
    println("mean R0: $(mean(cd)) with std: $(std(cd))")
    return cd
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