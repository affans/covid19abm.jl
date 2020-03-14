using Distributed
using Base.Filesystem
using DataFrames
using CSV
using Query
using Statistics
using UnicodePlots
using ClusterManagers
using Dates

## load the packages by covid19abm
using Parameters, Distributions, StatsBase, StaticArrays, Random, Match, DataFrames
using covid19abm
const cv=covid19abm

#addprocs(30, exeflags="--project=.")
#@everywhere using covid19abm

addprocs(SlurmManager(512), N=16, topology=:master_worker, exeflags="--project=.")
@everywhere using covid19abm

function run(myp::ModelParameters, folderprefix="./")
    nsims = 1000      ## robustness function gives us stable point r0 for 1000 sims
    println("starting $nsims simulations...")
    println("save folder set to $(pwd())...")
    println("running with parameters...")
    dump(myp)
    # will return 5 dataframes. 1 total, 4 age-specific 
    cd = pmap(1:nsims) do x         
        reset_params(myp)  ## reset parameters for each simulation. 
        hmatrix, ags = main(x)        
        all = _collectdf(hmatrix)
        spl = _splitstate(hmatrix, ags)
        ag1 = _collectdf(spl[1])
        ag2 = _collectdf(spl[2])
        ag3 = _collectdf(spl[3])
        ag4 = _collectdf(spl[4])
        return (a=all, g1=ag1, g2=ag2, g3=ag3, g4=ag4)
    end    
    for i = 1:nsims
         dts = cd[i]
         for dt in dts 
            insertcols!(dt, 1, :sim => i)    
            insertcols!(dt, 1, :time => 1:myp.modeltime)
         end         
    end
    println("simulations finished")
    println("total size of simulations: $(Base.summarysize(cd))")

    ## stack the sims together
    all = vcat([cd[i].a  for i = 1:nsims]...)
    ag1 = vcat([cd[i].g1 for i = 1:nsims]...)
    ag2 = vcat([cd[i].g2 for i = 1:nsims]...)
    ag3 = vcat([cd[i].g3 for i = 1:nsims]...)
    ag4 = vcat([cd[i].g4 for i = 1:nsims]...)
    mydfs = Dict("all" => all, "ag1" => ag1, "ag2" => ag2, "ag3" => ag3, "ag4" => ag4)

    ## save at the simulation and time level
    ## to ignore for now: miso, iiso, mild, ded
    c1 = Symbol.((:LAT, :INF, :HOS, :ICU, :DED), :_INC)
    c2 = Symbol.((:LAT, :INF, :HOS, :ICU, :DED), :_PREV)
    for (k, df) in mydfs
        println("saving dataframe sim level: $k")
        # simulation level, save file per health status, per age group
        for c in vcat(c1..., c2...)
            udf = unstack(df, :time, :sim, c) 
            fn = (lowercase(string("$(folderprefix)/simlevel_", c, "_", k, ".dat")))
            CSV.write(fn, udf)
        end
        println("saving dataframe time level: $k")
        # time level, save file per age group
        yaf = compute_yearly_average(df)
        fn = string("$(folderprefix)/timelevel_", k, ".dat")        
        CSV.write(fn, yaf)
    end
    return mydfs
end

function compute_yearly_average(df)
    ya = df |> @groupby(_.time) |> @map({time=key(_), cnt=length(_),
              sus_prev=mean(_.SUS_PREV), 
              lat_prev=mean(_.LAT_PREV), 
              mild_prev=mean(_.MILD_PREV), 
              miso_prev=mean(_.MISO_PREV), 
              inf_prev=mean(_.INF_PREV), 
              iiso_prev=mean(_.IISO_PREV), 
              hos_prev=mean(_.HOS_PREV), 
              icu_prev=mean(_.ICU_PREV), 
              rec_prev=mean(_.REC_PREV), 
              ded_prev=mean(_.DED_PREV), 
              sus_inc=mean(_.SUS_INC),
              lat_inc=mean(_.LAT_INC), 
              mild_inc=mean(_.MILD_INC), 
              miso_inc=mean(_.MISO_INC), 
              inf_inc=mean(_.INF_INC),
              iiso_inc=mean(_.IISO_INC),
              hos_inc=mean(_.HOS_INC),
              icu_inc=mean(_.ICU_INC),
              rec_inc=mean(_.REC_INC),
              ded_inc=mean(_.DED_INC)
              }) |> DataFrame
    return ya
end

function run_scenarios()
    myp = covid19abm.ModelParameters()
    start = time()
    betas = [0.035] 
    prov = [:ontario, :alberta, :bc, :manitoba, :newbruns, :newfdland, :nwterrito, :novasco, :nunavut, :pei, :quebec, :saskat, :yukon]
    fs = (0.05, 0.1, 0.2)
    τs = (1, 2)      
    ts = length(betas) * length(fs) * length(τs) 
    #pr = Progress(ts, dt=1, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:yellow)
    for p in prov, r in betas
        myp.β = r
        myp.prov = p        
        ## run with no isolation
        myp.fmild = 0 
        myp.τmild = 0                 
        prefix = savestr(myp)
        run(myp, prefix)
        
        ## run with isolation 
        for f in fs, τ in τs    
            myp.τmild  = τ
            myp.fmild = f            
            prefix = savestr(myp)
            run(myp, prefix)
        end
    end  
    elapsed = time() - start
    println("done, time: $elapsed")
end

function savestr(p::ModelParameters)
    datestr = (Dates.format(Dates.now(), dateformat"mmdd_HHMM"))
    ## setup folder name based on model parameters
    taustr = replace(string(p.τmild), "." => "")
    fstr = replace(string(p.fmild), "." => "")
    rstr = replace(string(p.β), "." => "")
    prov = replace(string(p.prov), "." => "")
    fldrname = "/data/covid19abm/simresults/$prov/b$rstr/tau$(taustr)_f$(fstr)/"
    mkpath(fldrname)
end

function calibrate(beta, nsims, prov=:ontario)
    println("did you remember to turn off hospitalization/icu for the initial infected case? Check move_to_inf()")
    myp = ModelParameters()
    myp.β = beta
    myp.prov = prov
    myp.calibration = true
    vals = zeros(Int64, nsims)
    println("calibrating with beta: $beta, total sims: $nsims, province: $prov")
    cd = pmap(1:nsims) do i 
        reset_params(myp)  
        h, ags = main(i) ## gets the entire model. 
        val = sum(_get_column_incidence(h, covid19abm.LAT))            
        val = val - 1 ## minus becuase the initial latent guy ends up in latent by the end cuz of swapupdate()
        return val
    end
    println("mean R0: $(mean(cd)) with std: $(std(cd))")
    return mean(cd)
end

function calibrate_robustness(beta, prov=:ontario)
    #[:ontario, :alberta, :bc, :manitoba, :newbruns, :newfdland, :nwterrito, :novasco, :nunavut, :pei, :quebec, :saskat, :yukon]
    # once a beta is found based on nsims simulations, 
    # see how robust it is. run calibration with same beta 100 times 
    # to see the variation in R0 produced. 
    nsims = [500, 1000, 2000]
    reps = 5
    means = zeros(Float64, reps, length(nsims))
    for (i, ns) in enumerate(nsims)
        cd = map(1:reps) do x 
            println("iter: $x, sims: $ns")
            mval = calibrate(beta, ns, prov)         
            return mval
        end
        means[:, i] = cd
    end
    # for i in 2:nworkers()
    #     ## mf defined as: @everywhere mg() = covid19abm.p.β     
    #     rpr = remotecall_fetch(mf,  i+1).prov
    #     rpr != prov && error("province didn't get set in the remote workers")
    # end
    return means
end
#provs = []
#for x in [:ontario, :alberta, :bc, :manitoba, :newbruns]
#    m = calibrate_robustness(0.101, x)
#    push!(provs, m)
#end

#@everywhere mg() = covid19abm.p.β    
#@everywhere mf(b) = covid19abm.p.β = b
# function hmatrix_analysis()
#     ## function takes the output of the model 
#     ## and splits it into an array of 2d matrices for plotting purposes. 
#     myp = cv.ModelParameters()
#     myp.β = 0.038
#     reset_params($myp)
#     cd, ags = main(1)

#     #mats = Array{Int64, 3}(undef, 500) 
#     mats = Array{Int64, 3}(undef, 100,100,500)
#     for i = 1:500
#         col = cd[:, i]
#         mats[:, :, i] = reshape(col, (100, 100))
#     end

#     colorarray = IndirectArray(mats, [colorant"green", colorant"teal", colorant"blue", colorant"azure", 
#     colorant"red", colorant"lightsalmon", colorant"magenta", colorant"magenta", colorant"grey60", colorant"black"])
#     save("myanimation.gif", colorarray; fps=5)

#     idxarray = rand(1:3, 100, 100, 50);
#     colorarray = IndirectArray(idxarray, [colorant"red", colorant"purple", colorant"indigo"])

#     return mats
# end

# function makeheatmaps(bufs)
#     heatmaps = map(bufs) do buf
#         heatmap(
#             buf, scale_plot = false, show_axis = false)

#     end
#     scene = hbox(map(i-> vbox(heatmaps[i, :]), 1:size(bufs, 1))...)
#     scene, last.(heatmaps)
# end
#datarows = 500; datacols = 500
#plotrows = 4; plotcols = 4
#bufs = [fill(0.0f0, datarows, datacols) for _ in 1:plotrows, _ in 1:plotcols]

# http://juliaplots.org/MakieReferenceImages/gallery//lots_of_heatmaps/index.html
# http://juliaplots.org/MakieReferenceImages/gallery//chess_game/index.html
