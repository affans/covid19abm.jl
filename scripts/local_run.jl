using Distributed
using Base.Filesystem
using DataFrames
using CSV
using Query
using Statistics
using UnicodePlots
using covid19abm
const cv=covid19abm

addprocs(4, exeflags="--project=.")
@everywhere using covid19abm

function run(beta) 
    myp = cv.ModelParameters()
    myp.β = beta
    @everywhere reset_params($myp)
    nsims = 4     
    println("starting $nsims simulations...")
    dump(myp)
    # will return 5 dataframes. 1 total, 4 age-specific 
    cd = pmap(1:nsims) do x         
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

    ## stack the sims together
    all = vcat([cd[i].a  for i = 1:nsims]...)
    ag1 = vcat([cd[i].g1 for i = 1:nsims]...)
    ag2 = vcat([cd[i].g2 for i = 1:nsims]...)
    ag3 = vcat([cd[i].g3 for i = 1:nsims]...)
    ag4 = vcat([cd[i].g4 for i = 1:nsims]...)
    mydfs = Dict("all" => all, "ag1" => ag1, "ag2" => ag2, "ag3" => ag3, "ag4" => ag4)

    ## save at the simulation and time level
    c1 = Symbol.((:SUS, :LAT, :MILD, :MISO, :INF, :IISO, :HOS, :ICU), :_INC)
    c2 = Symbol.((:SUS, :LAT, :MILD, :MISO, :INF, :IISO, :HOS, :ICU), :_PREV)
    for (k, df) in mydfs
        println("saving dataframe sim level: $k")
        # simulation level, save file per health status, per age group
        for c in vcat(c1..., c2...)
            udf = unstack(df, :time, :sim, c) 
            #fn = (lowercase(string("simlevel_", c, "_", k, "_.dat")))
            #CSV.write(fn, udf)
        end
        println("saving dataframe time level: $k")
        # time level, save file per age group
        yaf = compute_yearly_average(df)
        fn = string("timelevel_", k, ".dat")
        println("typeof: $(typeof(yaf))")
        CSV.write(fn, yaf)
    end
    return cd
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

function calibrate(beta)
    myp = ModelParameters()
    myp.β = beta
    myp.calibration = true
    @everywhere reset_params($myp)
    
    nsims = 16
    vals = zeros(Int64, nsims)
    println("calibrating with beta: $(cv.p.β)")
    cd = pmap(1:nsims) do i 
        h, ags = main(i) ## gets the entire model. 
        val = sum(_get_column_incidence(h, covid19abm.LAT))            
        val = val - 1 ## minus becuase the initial latent guy ends up in latent by the end cuz of swapupdate()
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




function split_model(pic)
    ## function takes the output of the model 
    ## and splits it into an array of 2d matrices for plotting purposes. 
    # go through each column
    mats = Array{Array{Int64, 2}, 1}(undef, 500) 
    for i = 1:500
        col = pic[:, i]
        mats[i] = reshape(col, (100, 100))
    end
    return mats
end

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
