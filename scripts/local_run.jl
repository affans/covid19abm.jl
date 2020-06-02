using Distributed
using Base.Filesystem
using DataFrames
using CSV
using Query
using Statistics
using UnicodePlots
using ClusterManagers
using Dates
using DelimitedFiles

## load the packages by covid19abm
using Parameters, Distributions, StatsBase, StaticArrays, Random, Match, DataFrames
using covid19abm
const cv=covid19abm

#addprocs(30, exeflags="--project=.")
#@everywhere using covid19abm

addprocs(SlurmManager(288), N=9, topology=:master_worker, exeflags="--project=.")
@everywhere using covid19abm

function run(myp::ModelParameters, nsims=500, folderprefix="./")
    println("starting $nsims simulations...\nsave folder set to $(folderprefix)")
    #dump(myp)
    myp.calibration && error("can not run simulation, calibration is on.")
    # will return 6 dataframes. 1 total, 4 age-specific 
    cdr = pmap(1:nsims) do x                 
            cv.runsim(x, myp)
    end      

    println("simulations finished")
    println("total size of simulation dataframes: $(Base.summarysize(cdr))")
    ## write the infectors 
    writedlm("$(folderprefix)/infectors.dat", [cdr[i].infectors for i = 1:nsims])    

    ## write contact numbers
    writedlm("$(folderprefix)/ctnumbers.dat", [cdr[i].ct_numbers for i = 1:nsims])    
    ## stack the sims together
    allag = vcat([cdr[i].a  for i = 1:nsims]...)
    ag1 = vcat([cdr[i].g1 for i = 1:nsims]...)
    ag2 = vcat([cdr[i].g2 for i = 1:nsims]...)
    ag3 = vcat([cdr[i].g3 for i = 1:nsims]...)
    ag4 = vcat([cdr[i].g4 for i = 1:nsims]...)
    ag5 = vcat([cdr[i].g5 for i = 1:nsims]...)
    #mydfs = Dict("all" => all, "ag1" => ag1, "ag2" => ag2, "ag3" => ag3, "ag4" => ag4, "ag5" => ag5)
    mydfs = Dict("all" => allag)
    
    ## save at the simulation and time level
    ## to ignore for now: miso, iiso, mild 
    #c1 = Symbol.((:LAT, :ASYMP, :MILD, :MISO, :INF, :IISO, :HOS, :ICU, :DED), :_INC)
    #c2 = Symbol.((:LAT, :ASYMP, :MILD, :MISO, :INF, :IISO, :HOS, :ICU, :DED), :_PREV)
    c1 = Symbol.((:LAT, ), :_INC)
    c2 = Symbol.((:LAT, ), :_PREV)
    for (k, df) in mydfs
        println("saving dataframe sim level: $k")
        # simulation level, save file per health status, per age group
        for c in vcat(c1..., c2...)
            udf = unstack(df, :time, :sim, c) 
            fn = string("$(folderprefix)/simlevel_", lowercase(string(c)), "_", k, ".dat")
            CSV.write(fn, udf)
        end
        println("saving dataframe time level: $k")
        # time level, save file per age group
        yaf = compute_yearly_average(df)       
        fn = string("$(folderprefix)/timelevel_", k, ".dat")   
        CSV.write(fn, yaf)       

        # write the secondary (Ro) values 
        secvals = unstack(df, :time, :sim, :secondarys)
        fn = string("$(folderprefix)/secondarys_", k, ".dat")   
        CSV.write(fn, secvals)       
    end
    return mydfs
end

function compute_yearly_average(df)
    ya = df |> @groupby(_.time) |> @map({time=key(_), cnt=length(_),
              sus_prev=mean(_.SUS_PREV), 
              lat_prev=mean(_.LAT_PREV), 
              pre_prev=mean(_.PRE_PREV), 
              asymp_prev=mean(_.ASYMP_PREV), 
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
              pre_inc=mean(_.PRE_INC), 
              asymp_inc=mean(_.ASYMP_INC), 
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

function savestr(p::ModelParameters, custominsert="simresults", customstart="")
    datestr = (Dates.format(Dates.now(), dateformat"mmdd_HHMM"))
    ## setup folder name based on model parameters
    taustr = replace(string(p.τmild), "." => "")
    fstr = replace(string(p.fmild), "." => "")
    rstr = replace(string(p.β), "." => "")
    prov = replace(string(p.prov), "." => "")
    eldr = replace(string(p.eldq), "." => "")
    eldqag = replace(string(p.eldqag), "." => "")     
    fpreiso = replace(string(p.fpreiso), "." => "")
    tpreiso = replace(string(p.tpreiso), "." => "")
    fsev = replace(string(p.fsevere), "." => "")    
    frelasymp = replace(string(p.frelasymp), "." => "")
    fasymp = replace(string(p.fasymp), "." => "") 
    fasympiso = replace(string(p.fasympiso), "." => "")
    if p.ctstrat == 1 
        strat = "strat1sc" 
    elseif p.ctstrat == 3
        strat = "strat3s" * string(p.strat3qdays)
    else 
        strat = "strat0"
    end
    pct = replace(string(p.fctcapture), "." => "")
    cct = replace(string(p.fcontactst), "." => "")
    idt = replace(string(p.cidtime), "." => "") 
    tback = replace(string(p.cdaysback), "." => "")     
    fldrname = "/data/covid19abm/$(custominsert)/$(customstart)_$(prov)_$(strat)_pct$(pct)_cct$(cct)_idt$(idt)_tback$(tback)_fsev$(fsev)_tau$(taustr)_fmild$(fstr)_q$(eldr)_qag$(eldqag)_relasymp$(frelasymp)_asymp$(fasymp)_tpreiso$(tpreiso)_preiso$(fpreiso)_asympiso$(fasympiso)/"
    mkpath(fldrname)
end

function _calibrate(nsims, myp::ModelParameters)
    myp.calibration != true && error("calibration parameter not turned on")
    vals = zeros(Int64, nsims)
    println("calibrating with beta: $(myp.β), total sims: $nsims, province: $(myp.prov)")
    println("calibration parameters:")
    dump(myp)
    cdr = pmap(1:nsims) do i 
        h = main(myp) ## gets the entire model. 
        val = sum(_get_column_incidence(h, covid19abm.LAT))            
        return val
    end
    return mean(cdr), std(cdr)
end

function calibrate(beta, nsims, init_inf, size, prov=:ontario)
    myp = ModelParameters() # set up default parameters 
    myp.β = beta
    myp.prov = prov
    myp.popsize = size
    myp.modeltime = 30
    myp.calibration = true
    myp.initialinf = init_inf
    m, sd = _calibrate(nsims, myp)
    println("mean R0: $(m) with std: $(sd)")
    myp.calibration = false       
    return m
end

function calibrate_robustness(beta, reps, prov=:ontario)
    #[:ontario, :alberta, :bc, :manitoba, :newbruns, :newfdland, :nwterrito, :novasco, :nunavut, :pei, :quebec, :saskat, :yukon]
    # once a beta is found based on nsims simulations, 
    # see how robust it is. run calibration with same beta 100 times 
    # to see the variation in R0 produced. 
    nsims = [500, 1000]
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

