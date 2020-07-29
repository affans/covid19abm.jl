module covid19abm
using Parameters, Distributions, StatsBase, StaticArrays, Random, Match, DataFrames

@enum HEALTH SUS LAT PRE ASYMP MILD MISO INF IISO HOS ICU REC DED UNDEF

Base.@kwdef mutable struct Human
    idx::Int64 = 0 
    health::HEALTH = SUS
    swap::HEALTH = UNDEF
    sickfrom::HEALTH = UNDEF
    nextday_meetcnt::Int16 = 0 ## how many contacts for a single day
    age::Int16   = 0    # in years. don't really need this but left it incase needed later
    ag::Int16   = 0
    tis::Int16   = 0   # time in state 
    exp::Int16   = 0   # max statetime
    dur::NTuple{4, Int8} = (0, 0, 0, 0)   # Order: (latents, asymps, pres, infs) TURN TO NAMED TUPS LATER
    doi::Int16   = 999   # day of infection.
    iso::Bool = false  ## isolated (limited contacts)
    isovia::Symbol = :null ## isolated via quarantine (:qu), preiso (:pi), intervention measure (:im), or contact tracing (:ct)    
    tracing::Bool = false ## are we tracing contacts for this individual?
    tracestart::Int16 = -1 ## when to start tracing, based on values sampled for x.dur
    traceend::Int16 = -1 ## when to end tracing
    tracedby::UInt32 = 0 ## is the individual traced? property represents the index of the infectious person 
    tracedxp::Int16 = 0 ## the trace is killed after tracedxp amount of days
end

## default system parameters
@with_kw mutable struct ModelParameters @deftype Float64    ## use @with_kw from Parameters
    β = 0.0       
    seasonal::Bool = false ## seasonal betas or not
    popsize::Int64 = 10000
    prov::Symbol = :newyork 
    calibration::Bool = false 
    modeltime::Int64 = 300
    initialinf::Int64 = 1
    initialhi::Int64 = 0 ## initial herd immunity, inserts number of REC individuals
    τmild::Int64 = 0 ## days before they self-isolate for mild cases
    fmild::Float64 = 0.0  ## percent of people practice self-isolation
    fsevere::Float64 = 0.0 # fixed at 0.80
    eldq::Float64 = 0.0 ## complete isolation of elderly
    eldqag::Int8 = 5 ## default age group, if quarantined(isolated) is ag 5. 
    fasymp::Float64 = 0.50 ## NOT USED ANYMORE ## percent going to asymp (may 10, removed all tests and references)
    fpre::Float64 = 1.0 ## NOT USED ANYMORE (percent going to presymptomatic)
    fpreiso::Float64 = 0.0 ## percent that is isolated at the presymptomatic stage
    tpreiso::Int64 = 0## preiso is only turned on at this time. 
    frelasymp::Float64 = 0.11 ## relative transmission of asymptomatic
    ctstrat::Int8 = 0 ## strategy 
    fctcapture::Float16 = 0.0 ## how many symptomatic people identified
    fcontactst::Float16 = 0.0 ## fraction of contacts being isolated/quarantined
    cidtime::Int8 = 0  ## time to identification (for CT) post symptom onset, NOT USED ANYMORE
    cdaysback::Int8 = 0 ## number of days to go back and collect contacts
    strat3qdays::Int8 = 0
end

Base.@kwdef mutable struct ct_data_collect
    total_symp_id::Int64 = 0  # total symptomatic identified
    totaltrace::Int64 = 0     # total contacts traced
    totalisolated::Int64 = 0  # total number of people isolated
    iso_sus::Int64 = 0        # total susceptible isolated IN DAYS
    iso_lat::Int64 = 0        # total latent isolated IN DAYS
    iso_asymp::Int64 = 0      # total asymp isolated IN DAYS
    iso_symp::Int64 = 0       # total symp (mild, inf) isolated IN DAYS
end

Base.show(io::IO, ::MIME"text/plain", z::Human) = dump(z)

## constants 
const humans = Array{Human}(undef, 0) 
const p = ModelParameters()  ## setup default parameters
const agebraks = @SVector [0:4, 5:19, 20:49, 50:64, 65:99]
const BETAS = Array{Float64, 1}(undef, 0) ## to hold betas (whether fixed or seasonal), array will get resized
const ct_data = ct_data_collect()
export ModelParameters, HEALTH, Human, humans, BETAS

function runsim(simnum, ip::ModelParameters)
    Random.seed!(simnum*726)

    # function runs the `main` function, and collects the data as dataframes. 
    hmatrix = main(ip)            
    # get infectors counters
    infectors = _count_infectors()

    ct_numbers = (ct_data.total_symp_id, ct_data.totaltrace, ct_data.totalisolated, 
                    ct_data.iso_sus, ct_data.iso_lat, ct_data.iso_asymp, ct_data.iso_symp)
    
    # get simulation age groups
    ags = [x.ag for x in humans] # store a vector of the age group distribution 
    all = _collectdf(hmatrix)
    spl = _splitstate(hmatrix, ags)
    ag1 = _collectdf(spl[1])
    ag2 = _collectdf(spl[2])
    ag3 = _collectdf(spl[3])
    ag4 = _collectdf(spl[4])
    ag5 = _collectdf(spl[5])
    insertcols!(all, 1, :sim => simnum); insertcols!(ag1, 1, :sim => simnum); insertcols!(ag2, 1, :sim => simnum); 
    insertcols!(ag3, 1, :sim => simnum); insertcols!(ag4, 1, :sim => simnum); insertcols!(ag5, 1, :sim => simnum); 
    return (a=all, g1=ag1, g2=ag2, g3=ag3, g4=ag4, g5=ag5, infectors=infectors, ct_numbers=ct_numbers)
end
export runsim

function main(ip::ModelParameters)
    ## datacollection            
    # matrix to collect model state for every time step

    # reset the parameters for the simulation scenario
    reset_params(ip)

    p.popsize == 0 && error("no population size given")
    
    hmatrix = zeros(Int16, p.popsize, p.modeltime)
    initialize() # initialize population
    # insert initial infected agents into the model
    # and setup the right swap function. 
    if p.calibration 
        insert_infected(PRE, p.initialinf, 4)
    else 
        insert_infected(PRE, p.initialinf, 4)  
        insert_infected(REC, p.initialhi, 4)
    end    
    
    ## save the preisolation isolation parameters
    _fpreiso = p.fpreiso
    p.fpreiso = 0

    # split population in agegroups 
    grps = get_ag_dist()

    # start the time loop
    for st = 1:p.modeltime
        # start of day
        if st == p.tpreiso ## time to introduce testing
            p.fpreiso = _fpreiso
        end
        _get_model_state(st, hmatrix) ## this datacollection needs to be at the start of the for loop
        dyntrans(st, grps)
        sw = time_update()
        # end of day
    end
    return hmatrix ## return the model state as well as the age groups. 
end
export main

reset_params() = reset_params(ModelParameters())
function reset_params(ip::ModelParameters)
    # the p is a global const
    # the ip is an incoming different instance of parameters 
    # copy the values from ip to p. 
    for x in propertynames(p)
        setfield!(p, x, getfield(ip, x))
    end

    # reset the contact tracing data collection structure
    for x in propertynames(ct_data)
        setfield!(ct_data, x, 0)
    end

    # resize and update the BETAS constant array
    init_betas()

    # resize the human array to change population size
    resize!(humans, p.popsize)
end
export reset_params

function _model_check() 
    ## checks model parameters before running 
    ## move code from here to unit tests.
    (p.fctcapture > 0 && p.fpreiso > 0) && error("Can not do contact tracing and ID/ISO of pre at the same time.")
end

## Data Collection/ Model State functions
function _get_model_state(st, hmatrix)
    # collects the model state (i.e. agent status at time st)
    for i=1:length(humans)
        hmatrix[i, st] = Int(humans[i].health)
    end    
end
export _get_model_state

function _collectdf(hmatrix)
    ## takes the output of the humans x time matrix and processes it into a dataframe
    #_names_inci = Symbol.(["lat_inc", "mild_inc", "miso_inc", "inf_inc", "iiso_inc", "hos_inc", "icu_inc", "rec_inc", "ded_inc"])    
    #_names_prev = Symbol.(["sus", "lat", "mild", "miso", "inf", "iiso", "hos", "icu", "rec", "ded"])
    mdf_inc, mdf_prev = _get_incidence_and_prev(hmatrix)
    mdf = hcat(mdf_inc, mdf_prev)    
    _names_inc = Symbol.(string.((Symbol.(instances(HEALTH)[1:end - 1])), "_INC"))
    _names_prev = Symbol.(string.((Symbol.(instances(HEALTH)[1:end - 1])), "_PREV"))
    _names = vcat(_names_inc..., _names_prev...)
    datf = DataFrame(mdf, _names)

    lat_inc = mdf_inc[:, 2]
    inf_prv = vec(sum(mdf_prev[:, 3:8], dims=2))
    secondarys = lat_inc ./ inf_prv 
    replace!(secondarys, NaN => 0)
    replace!(secondarys, Inf => 0)
    insertcols!(datf, 1, :time => 1:p.modeltime) ## add a time column to the resulting dataframe
    insertcols!(datf, 2, :secondarys => secondarys) ## add a time column to the resulting dataframe
    return datf
end

function _splitstate(hmatrix, ags)
    #split the full hmatrix into 4 age groups based on ags (the array of age group of each agent)
    #sizes = [length(findall(x -> x == i, ags)) for i = 1:4]
    matx = []#Array{Array{Int64, 2}, 1}(undef, 4)
    for i = 1:length(agebraks)
        idx = findall(x -> x == i, ags)
        push!(matx, view(hmatrix, idx, :))
    end
    return matx
end
export _splitstate

function _get_incidence_and_prev(hmatrix)
    cols = instances(HEALTH)[1:end - 1] ## don't care about the UNDEF health status
    inc = zeros(Int64, p.modeltime, length(cols))
    pre = zeros(Int64, p.modeltime, length(cols))
    for i = 1:length(cols)
        inc[:, i] = _get_column_incidence(hmatrix, cols[i])
        pre[:, i] = _get_column_prevalence(hmatrix, cols[i])
    end
    return inc, pre
end

function _get_column_incidence(hmatrix, hcol)
    inth = Int(hcol)
    timevec = zeros(Int64, p.modeltime)
    for r in eachrow(hmatrix)
        idx = findfirst(x -> x == inth, r)
        if idx !== nothing 
            timevec[idx] += 1
        end
    end
    return timevec
end

function _get_column_prevalence(hmatrix, hcol)
    inth = Int(hcol)
    timevec = zeros(Int64, p.modeltime)
    for (i, c) in enumerate(eachcol(hmatrix))
        idx = findall(x -> x == inth, c)
        if idx !== nothing
            ps = length(c[idx])    
            timevec[i] = ps    
        end
    end
    return timevec
end

function _count_infectors()     
    pre_ctr = asymp_ctr = mild_ctr = inf_ctr = 0
    for x in humans 
        if x.health != SUS ## meaning they got sick at some point
            if x.sickfrom == PRE
                pre_ctr += 1
            elseif x.sickfrom == ASYMP
                asymp_ctr += 1
            elseif x.sickfrom == MILD || x.sickfrom == MISO 
                mild_ctr += 1 
            elseif x.sickfrom == INF || x.sickfrom == IISO 
                inf_ctr += 1 
            else 
                error("sickfrom not set right: $(x.sickfrom)")
            end
        end
    end
    return (pre_ctr, asymp_ctr, mild_ctr, inf_ctr)
end

export _collectdf, _get_incidence_and_prev, _get_column_incidence, _get_column_prevalence, _count_infectors

## initialization functions 
function get_province_ag(prov) 
    ret = @match prov begin        
        :alberta => Distributions.Categorical(@SVector [0.0655, 0.1851, 0.4331, 0.1933, 0.1230])
        :bc => Distributions.Categorical(@SVector [0.0475, 0.1570, 0.3905, 0.2223, 0.1827])
        :canada => Distributions.Categorical(@SVector [0.0540, 0.1697, 0.3915, 0.2159, 0.1689])
        :manitoba => Distributions.Categorical(@SVector [0.0634, 0.1918, 0.3899, 0.1993, 0.1556])
        :newbruns => Distributions.Categorical(@SVector [0.0460, 0.1563, 0.3565, 0.2421, 0.1991])
        :newfdland => Distributions.Categorical(@SVector [0.0430, 0.1526, 0.3642, 0.2458, 0.1944])
        :nwterrito => Distributions.Categorical(@SVector [0.0747, 0.2026, 0.4511, 0.1946, 0.0770])
        :novasco => Distributions.Categorical(@SVector [0.0455, 0.1549, 0.3601, 0.2405, 0.1990])
        :nunavut => Distributions.Categorical(@SVector [0.1157, 0.2968, 0.4321, 0.1174, 0.0380])
        :ontario => Distributions.Categorical(@SVector [0.0519, 0.1727, 0.3930, 0.2150, 0.1674])
        :pei => Distributions.Categorical(@SVector [0.0490, 0.1702, 0.3540, 0.2329, 0.1939])
        :quebec => Distributions.Categorical(@SVector [0.0545, 0.1615, 0.3782, 0.2227, 0.1831])
        :saskat => Distributions.Categorical(@SVector [0.0666, 0.1914, 0.3871, 0.1997, 0.1552])
        :yukon => Distributions.Categorical(@SVector [0.0597, 0.1694, 0.4179, 0.2343, 0.1187])
        :newyork   => Distributions.Categorical(@SVector [0.064000, 0.163000, 0.448000, 0.181000, 0.144000])
        _ => error("shame for not knowing your canadian provinces and territories")
    end       
    return ret  
end
export get_province_ag

function initialize() 
    agedist = get_province_ag(p.prov)
    for i = 1:p.popsize 
        humans[i] = Human()              ## create an empty human       
        x = humans[i]
        x.idx = i 
        x.ag = rand(agedist)
        x.age = rand(agebraks[x.ag]) 
        x.exp = 999  ## susceptible people don't expire.
        x.dur = sample_epi_durations() # sample epi periods   
        if rand() < p.eldq && x.ag == p.eldqag   ## check if elderly need to be quarantined.
            x.iso = true   
            x.isovia = :qu         
        end
        # initialize the next day counts (this is important in initialization since dyntrans runs first)
        get_nextday_counts(x)
    end
end
export initialize

function init_betas() 
    if p.seasonal  
        tmp = p.β .* td_seasonality()
    else 
        tmp = p.β .* ones(Float64, p.modeltime)
    end
    resize!(BETAS, length(tmp))
    for i = 1:length(tmp)
        BETAS[i] = tmp[i]
    end
end

function td_seasonality()
    ## returns a vector of seasonal oscillations
    t = 1:p.modeltime
    a0 = 6.261
    a1 = -11.81
    b1 = 1.817
    w = 0.022 #0.01815    
    temp = @. a0 + a1*cos((80-t)*w) + b1*sin((80-t)*w)  #100
    #temp = @. a0 + a1*cos((80-t+150)*w) + b1*sin((80-t+150)*w)  #100
    temp = (temp .- 2.5*minimum(temp))./(maximum(temp) .- minimum(temp)); # normalize  @2
    return temp
end

function get_ag_dist() 
    # splits the initialized human pop into its age groups
    grps =  map(x -> findall(y -> y.ag == x, humans), 1:length(agebraks)) 
    return grps
end

function insert_infected(health, num, ag) 
    ## inserts a number of infected people in the population randomly
    ## this function should resemble move_to_inf()
    ## when health == INF, the ag matters. 
    ## for  health == REC, just use all humans 
    if health == INF
        l = findall(x -> x.health == SUS && x.ag == ag, humans)
    else
        l = findall(x -> x.health == SUS, humans)
    end
    # we can throw an error here because this function is run right at the start
    length(l) == 0 && error("No susceptibles in the population to change status")
    num > length(l) && error("No susceptibles in the population to change status")

    h = sample(l, num; replace = false)
    @inbounds for i in h 
        x = humans[i]
        if health == PRE 
            move_to_pre(x) ## the swap may be asymp, mild, or severe, but we can force severe in the time_update function
        elseif health == LAT 
            move_to_latent(x)
            #if the person is moving to latent, let them expire immediately. 
            #this removes the model artifact of "dipping" at the start of simulations.
            x.tis = x.exp  
        elseif health == INF
            move_to_infsimple(x)
        elseif health == REC 
            move_to_recovered(x)
        else 
            error("can not insert human of health $(health)")
        end       
        x.sickfrom = INF # this will add +1 to the INF count in _count_infectors()... keeps the logic simple in that function.    
    end
        
    return h
end
export insert_infected

function time_update()
    # counters to calculate incidence
    lat=0; pre=0; asymp=0; mild=0; miso=0; inf=0; infiso=0; hos=0; icu=0; rec=0; ded=0;
    for x in humans 
        x.tis += 1 
        x.doi += 1 # increase day of infection. variable is garbage until person is latent
        if x.tis >= x.exp             
            @match Symbol(x.swap) begin
                :LAT  => begin move_to_latent(x); lat += 1; end
                :PRE  => begin move_to_pre(x); pre += 1; end
                :ASYMP => begin move_to_asymp(x); asymp += 1; end
                :MILD => begin move_to_mild(x); mild += 1; end
                :MISO => begin move_to_miso(x); miso += 1; end
                :INF  => begin move_to_inf(x); inf +=1; end    
                :IISO => begin move_to_iiso(x); infiso += 1; end
                :HOS  => begin move_to_hospicu(x); hos += 1; end 
                :ICU  => begin move_to_hospicu(x); icu += 1; end
                :REC  => begin move_to_recovered(x); rec += 1; end
                :DED  => begin move_to_dead(x); ded += 1; end
                _    => error("swap expired, but no swap set.")
            end
        end
        # run covid-19 functions for other integrated dynamics. 
        ct_dynamics(x)

        # get the meet counts for the next day 
        get_nextday_counts(x)
    end
    return (lat, mild, miso, inf, infiso, hos, icu, rec, ded)
end
export time_update

@inline _set_isolation(x::Human, iso) = _set_isolation(x, iso, x.isovia)
@inline function _set_isolation(x::Human, iso, via)
    # a helper setter function to not overwrite the isovia property. 
    # a person could be isolated in susceptible/latent phase through contact tracing
    # --> in which case it will follow through the natural history of disease 
    # --> if the person remains susceptible, then iso = off
    # a person could be isolated in presymptomatic phase through fpreiso
    # --> if x.iso == true from CT and x.isovia == :ct, do not overwrite
    # a person could be isolated in mild/severe phase through fmild, fsevere
    # --> if x.iso == true from CT and x.isovia == :ct, do not overwrite
    # --> if x.iso == true from PRE and x.isovia == :pi, do not overwrite
    x.iso = iso 
    x.isovia == :null && (x.isovia = via)
end

function sample_epi_durations()
    # when a person is sick, samples the 
    lat_dist = truncated(LogNormal(log(5.2), 0.1), 4, 7) # truncated between 4 and 7
    pre_dist = truncated(Gamma(1.058, 5/2.3), 0.8, 3)#truncated between 0.8 and 3
    asy_dist = Gamma(5, 1)
    inf_dist = Gamma((3.2)^2/3.7, 3.7/3.2)

    latents = Int.(round.(rand(lat_dist)))
    pres = Int.(round.(rand(pre_dist)))
    latents = latents - pres # ofcourse substract from latents, the presymp periods
    asymps = Int.(ceil.(rand(asy_dist)))
    infs = Int.(ceil.(rand(inf_dist)))
    return (latents, asymps, pres, infs)
end

function move_to_latent(x::Human)
    ## transfers human h to the incubation period and samples the duration
    x.health = LAT
    x.doi = 0 ## day of infection is reset when person becomes latent
    x.tis = 0   # reset time in state 
    x.exp = x.dur[1] # get the latent period
    # the swap to asymptomatic is based on age group.
    # ask seyed for the references
    asymp_pcts = (0.25, 0.25, 0.14, 0.07, 0.07)    
    x.swap = rand() < asymp_pcts[x.ag] ? ASYMP : PRE 
    ## in calibration mode, latent people never become infectious.
    if p.calibration 
        x.swap = LAT 
        x.exp = 999
    end
end
export move_to_latent

function move_to_asymp(x::Human)
    ## transfers human h to the asymptomatic stage 
    x.health = ASYMP     
    x.tis = 0 
    x.exp = x.dur[2] # get the presymptomatic period
    x.swap = REC 
    # x.iso property remains from either the latent or presymptomatic class
    # if x.iso is true, the asymptomatic individual has limited contacts
end
export move_to_asymp

function move_to_pre(x::Human)
    θ = (0.95, 0.9, 0.85, 0.6, 0.2)  # percentage of sick individuals going to mild infection stage
    x.health = PRE
    x.tis = 0   # reset time in state 
    x.exp = x.dur[3] # get the presymptomatic period

    if rand() < θ[x.ag]
        x.swap = MILD
    else 
        x.swap = INF
    end
    # calculate whether person is isolated
    rand() < p.fpreiso && _set_isolation(x, true, :pi)
end
export move_to_pre

function move_to_mild(x::Human)
    ## transfers human h to the mild infection stage for γ days
    x.health = MILD     
    x.tis = 0 
    x.exp = x.dur[4]
    x.swap = REC 
    # x.iso property remains from either the latent or presymptomatic class
    # if x.iso is true, staying in MILD is same as MISO since contacts will be limited. 
    # we still need the separation of MILD, MISO because if x.iso is false, then here we have to determine 
    # how many days as full contacts before self-isolation
    # NOTE: if need to count non-isolated mild people, this is overestimate as isolated people should really be in MISO all the time
    #   and not go through the mild compartment 
    if x.iso || rand() < p.fmild
        x.swap = MISO  
        x.exp = p.τmild
    end
end
export move_to_mild

function move_to_miso(x::Human)
    ## transfers human h to the mild isolated infection stage for γ days
    x.health = MISO
    x.swap = REC
    x.tis = 0 
    x.exp = x.dur[4] - p.τmild  ## since tau amount of days was already spent as infectious
    _set_isolation(x, true, :mi) 
end
export move_to_miso

function move_to_infsimple(x::Human)
    ## transfers human h to the severe infection stage for γ days 
    ## simplified function for calibration/general purposes
    x.health = INF
    x.tis = 0 
    x.exp = x.dur[4]
    x.swap = REC 
    _set_isolation(x, false, :null) 
end

function move_to_inf(x::Human)
    ## transfers human h to the severe infection stage for γ days
    ## for swap, check if person will be hospitalized, selfiso, die, or recover
 
    # h = prob of hospital, c = prob of icu AFTER hospital    
    h = (rand(Uniform(0.02, 0.03)), rand(Uniform(0.02, 0.03)), rand(Uniform(0.28, 0.34)), rand(Uniform(0.28, 0.34)), rand(Uniform(0.60, 0.68)))
    c = (rand(Uniform(0.01, 0.015)), rand(Uniform(0.01, 0.015)), rand(Uniform(0.03, 0.05)), rand(Uniform(0.05, 0.20)), rand(Uniform(0.05, 0.15))) 
    mh = [0.01/5, 0.01/5, 0.0135/3, 0.01225/1.5, 0.04/2]     # death rate for severe cases.
    
    if p.calibration
        h =  (0, 0, 0, 0, 0)
        c =  (0, 0, 0, 0, 0)
        mh = (0, 0, 0, 0, 0)
    end

    time_to_hospital = Int(round(rand(Uniform(2, 5)))) # duration symptom onset to hospitalization
   
    x.health = INF
    x.swap = UNDEF
    x.tis = 0 
    if rand() < h[x.ag]     # going to hospital or ICU but will spend delta time transmissing the disease with full contacts 
        x.exp = time_to_hospital    
        x.swap = rand() < c[x.ag] ? ICU : HOS        
    else ## no hospital for this lucky (but severe) individual 
        if rand() < mh[x.ag]
            x.exp = x.dur[4]  
            x.swap = DED
        else 
            x.exp = x.dur[4]  
            x.swap = REC
            if x.iso || rand() < p.fsevere 
                x.exp = 1  ## 1 day isolation for severe cases     
                x.swap = IISO
            end  
        end
    end
    ## before returning, check if swap is set 
    x.swap == UNDEF && error("agent I -> ?")
end

function move_to_iiso(x::Human)
    ## transfers human h to the sever isolated infection stage for γ days
    x.health = IISO   
    x.swap = REC
    x.tis = 0     ## reset time in state 
    x.exp = x.dur[4] - 1  ## since 1 day was spent as infectious
    _set_isolation(x, true, :mi)
end 

function move_to_hospicu(x::Human)   
    ## to do: don't sample from distribution unless needed
    ## to do: check whether the parameters seyed gave are for mean/std or shape/scale
    #mh = [0.001, 0.001, 0.00135, 0.01225, 0.04]
    #mc = [0.002, 0.002, 0.0027, 0.0245, 0.08] 
    mh = [0.01*2, 0.01*2, 0.0135, 0.01225, 0.03*2]
    mc = 2*mh
    #mh = [0, 0, 0, 0.01225*2, 0.04*2]
    #mc = [0, 0, 0, 0.01225*3.5, 0.04*3.5]

    psiH = Int(round(rand(truncated(Gamma(4.5, 2.75), 8, 17))))
    psiC = Int(round(rand(truncated(Gamma(4.5, 2.75), 8, 17)))) + 2
    muH = Int(round(rand(truncated(Gamma(5.3, 2.1), 9, 15))))
    muC = Int(round(rand(truncated(Gamma(5.3, 2.1), 9, 15)))) + 2

    swaphealth = x.swap 
    x.health = swaphealth ## swap either to HOS or ICU
    x.swap = UNDEF
    x.tis = 0
    _set_isolation(x, true) # do not set the isovia property here.  

    if swaphealth == HOS 
        if rand() < mh[x.ag] ## person will die in the hospital 
            x.exp = muH 
            x.swap = DED
        else 
            x.exp = psiH 
            x.swap = REC
        end        
    end

    if swaphealth == ICU         
        if rand() < mc[x.ag] ## person will die in the ICU 
            x.exp = muC
            x.swap = DED
        else 
            x.exp = psiC
            x.swap = REC
        end
    end 
    ## before returning, check if swap is set 
    x.swap == UNDEF && error("agent H -> ?")    
end

function move_to_dead(h::Human)
    # no level of alchemy will bring someone back to life. 
    h.health = DED
    h.swap = UNDEF
    h.tis = 0 
    h.exp = 999 ## stay recovered indefinitely
    h.iso = true # a dead person is isolated
    _set_isolation(h, true)  # do not set the isovia property here.  
    # isolation property has no effect in contact dynamics anyways (unless x == SUS)
end

function move_to_recovered(h::Human)
    h.health = REC
    h.swap = UNDEF
    h.tis = 0 
    h.exp = 999 ## stay recovered indefinitely    
    # do not set isolation property for recovered. 
    # if they were isolated through contact tracing, they will unisolate automatically. 
    # to do: check other pathways. 
    # _set_isolation(h, false)  # do not set the isovia property here.      
end

function apply_ct_strategy(y::Human)
    iso = false  # to collect data at the end of the function
    yhealth = y.health
    if p.ctstrat == 1 
        # in strategy 1, all traced individuals are isolated for 14 days. 
        _set_isolation(y, true, :ct)
        iso = true
        y.tracedxp = 14 ## trace isolation will last for 14 days before expiry                
        ct_data.totalisolated += 1  ## update counter               
    end
    if p.ctstrat == 2 ## not really useful
        if y.health in (PRE, ASYMP, MILD, MISO, INF, IISO)
            _set_isolation(y, true, :ct)
            iso = true
            y.tracedxp = 14 ## trace isolation will last for 14 days before expiry                
            ct_data.totalisolated += 1  ## update counter
        else
            _set_isolation(y, false)
            # kill the trace
            y.tracedby = 0 
            y.tracedxp = 0 
        end
    end
    if p.ctstrat == 3
         # in strategy 3, all traced individuals are isolated for only 4 days. 
         _set_isolation(y, true, :ct)
         iso = true
         y.tracedxp = p.strat3qdays
         ct_data.totalisolated += 1  ## update counter 
    end
    # count data at time of infection... 
    # this is technical a bug... for example, a susceptible individual may be quarantined for 14 days 
    # but because they still have contacts, they may become sick and the category of the count should change. 
end

function ct_dynamics(x::Human)
    # main function for ct dynamics 
    # turns tracing on if person is infectious and in the right time window
    # applies isolation to all traced contacts
    # turns of isolation for susceptibles > 14 days
    # order of if statements matter here 
    (p.ctstrat == 0) && (return)
    xh = x.health 
    xs = x.swap
    dur = x.dur
    doi = x.doi
    
    ## person is newly infectious, and will show symptoms. 
    ## initialize tracestart and traceend times
    if xh == LAT && xs != ASYMP && doi == 0 
        if rand() < p.fctcapture 
            #delta = latent + presymp + time to identification time
            delta = dur[1] + dur[3] + Int(round(rand(Gamma(3.2, 1)))) #see chads email for reference
            q = delta - p.cdaysback  ## how many days to trace prior to identification
            x.tracestart = max(0, q) # minimum of zero since delta < p.cdaysback
            x.traceend = delta
            (x.tracestart > x.traceend) && error("tracestart < traceend")
            ct_data.total_symp_id += 1 ## update ct data ctr: total number of symptomatic identified
        end
    end

    ## turn on trace when day of infection == tracestart    
    if doi == x.tracestart 
        x.tracing = true 
    end 
    if doi == x.traceend 
        ## have to do more work here
        x.tracing = false 
        _set_isolation(x, true, :ct)
        alltraced = findall(y -> y.tracedby == x.idx, humans)
        for i in alltraced
            y = humans[i]
            yhealth = y.health
            apply_ct_strategy(y)
        end
    end

    if x.tracedxp > 0  # person was traced, isolated, and will be for x days
        x.iso == false && error("a traced person should always be isolated until the trace is killed")
        # collect stats on the traced individuals         
        if xh == INF || xh == MILD || xh == PRE
            ct_data.iso_symp += 1
        elseif xh == LAT 
            ct_data.iso_lat += 1
        elseif xh == ASYMP 
            ct_data.iso_asymp += 1
        else        
            ct_data.iso_sus += 1 
        end

        x.tracedxp -= 1
        if x.tracedxp == 0 
            # check whether isolation is turned on/off based on ctstrat
            if p.ctstrat == 3 
                rn = rand()
                # ask seyed for refernece, jeff
                if (x.health == LAT && rn < 0.05) || (x.health == PRE && rn < 0.95) || (x.health == ASYMP && rn < 0.70) || x.health in (MILD, INF, MISO, IISO)
                    # if strategy 3, only those that are tested positive are furthered isolated for 14 days.
                    _set_isolation(x, true) ## isovia should already be set to :ct 
                    x.tracedxp = 14 ## trace isolation will last for 14 days before expiry                
                    #ct_data.totalisolated += 1  ## update counter  (BUG: double counting)
                else 
                    _set_isolation(x, false) ## isovia should already be set to :ct 
                    x.isovia = :null # not isolated via contact tracing anymore
                    x.tracedby = 0 # the trace is killed
                end
            else
                # else their trace is killed
                _set_isolation(x, false)
                x.isovia = :null # not isolated via contact tracing anymore
                x.tracedby = 0 # the trace is killed
            end            
        end
    end
    return
end
export ct_dynamics

@inline function _get_betavalue(sys_time, xhealth) 
    #bf = p.β ## baseline PRE
    length(BETAS) == 0 && return 0
    bf = BETAS[sys_time]
    # values coming from FRASER Figure 2... relative tranmissibilities of different stages.
    if xhealth == ASYMP
        bf = bf * p.frelasymp
    elseif xhealth == MILD || xhealth == MISO 
        bf = bf * 0.44
    elseif xhealth == INF || xhealth == IISO 
        bf = bf * 0.89
    end
    return bf
end
export _get_betavalue

@inline function get_nextday_counts(x::Human)
    # get all people to meet and their daily contacts to recieve
    # we can sample this at the start of the simulation to avoid everyday    
    cnt = 0
    ag = x.ag
    #if person is isolated, they can recieve only 3 maximum contacts
    if x.iso 
        cnt = rand() < 0.5 ? 0 : rand(1:3)
    else 
        cnt = rand(nbs[ag])  # expensive operation, try to optimize
    end
    if x.health == DED 
        cnt = 0 
    end
    x.nextday_meetcnt = cnt
    return cnt
end

function dyntrans(sys_time, grps)
    totalmet = 0 # count the total number of contacts (total for day, for all INF contacts)
    totalinf = 0 # count number of new infected 
    ## find all the people infectious
    infs = findall(x -> x.health in (PRE, ASYMP, MILD, MISO, INF, IISO), humans)
    
    # go through every infectious person
    for xid in infs 
        x = humans[xid]
        xhealth = x.health
        cnts = x.nextday_meetcnt
        if cnts > 0  
            gpw = Int.(round.(cm[x.ag]*cnts)) # split the counts over age groups
            for (i, g) in enumerate(gpw)
                # sample the people from each group
                meet = rand(grps[i], g)
                # go through each person
                for j in meet 
                    y = humans[j]
                    ycnt = y.nextday_meetcnt             
                    
                    if ycnt > 0 
                        y.nextday_meetcnt = y.nextday_meetcnt - 1 # remove a contact
                        totalmet += 1
                        # there is a contact to recieve
                        # tracing dynamics
                        
                        if x.tracing  
                            if y.tracedby == 0 && rand() < p.fcontactst
                                y.tracedby = x.idx
                                ct_data.totaltrace += 1 
                            end
                        end
                        
                        # tranmission dynamics
                        if  y.health == SUS && y.swap == UNDEF                  
                            beta = _get_betavalue(sys_time, xhealth)
                            if rand() < beta
                                totalinf += 1
                                y.swap = LAT
                                y.exp = y.tis   ## force the move to latent in the next time step.
                                y.sickfrom = xhealth ## stores the infector's status to the infectee's sickfrom                            
                            end  
                        end    
                    end                    
                end
            end
        end        
    end
    return totalmet, totalinf
end
export dyntrans

function contact_matrix() 
    # regular contacts, just with 5 age groups. 
    #  0-4, 5-19, 20-49, 50-64, 65+
    CM = Array{Array{Float64, 1}, 1}(undef, 5)
    CM[1] = [0.2287, 0.1839, 0.4219, 0.1116, 0.0539]
    CM[2] = [0.0276, 0.5964, 0.2878, 0.0591, 0.0291]
    CM[3] = [0.0376, 0.1454, 0.6253, 0.1423, 0.0494]
    CM[4] = [0.0242, 0.1094, 0.4867, 0.2723, 0.1074]
    CM[5] = [0.0207, 0.1083, 0.4071, 0.2193, 0.2446]
    return CM
end
# 
# calibrate for 2.7 r0
# 20% selfisolation, tau 1 and 2.

function negative_binomials() 
    ## the means/sd here are calculated using _calc_avgag
    means = [10.21, 16.793, 13.7950, 11.2669, 8.0027]
    sd = [7.65, 11.7201, 10.5045, 9.5935, 6.9638]
    totalbraks = length(means)
    nbinoms = Vector{NegativeBinomial{Float64}}(undef, totalbraks)
    for i = 1:totalbraks
        p = 1 - (sd[i]^2-means[i])/(sd[i]^2)
        r = means[i]^2/(sd[i]^2-means[i])
        nbinoms[i] =  NegativeBinomial(r, p)
    end
    return nbinoms   
end
const nbs = negative_binomials()
const cm = contact_matrix()
export negative_binomials, contact_matrix, nbs, cm


## internal functions to do intermediate calculations
function _calc_avgag(lb, hb) 
    ## internal function to calculate the mean/sd of the negative binomials
    ## returns a vector of sampled number of contacts between age group lb to age group hb
    dists = _negative_binomials_15ag()[lb:hb]
    totalcon = Vector{Int64}(undef, 0)
    for d in dists 
        append!(totalcon, rand(d, 10000))
    end    
    return totalcon
end
export _calc_avgag

function _negative_binomials_15ag()
    ## negative binomials 15 agegroups
    AgeMean = Vector{Float64}(undef, 15)
    AgeSD = Vector{Float64}(undef, 15)
    #0-4, 5-9, 10-14, 15-19, 20-24, 25-29, 30-34, 35-39, 40-44, 45-49, 50-54, 55-59, 60-64, 65-69, 70+
    AgeMean = [10.21, 14.81, 18.22, 17.58, 13.57, 13.57, 14.14, 14.14, 13.83, 13.83, 12.3, 12.3, 9.21, 9.21, 6.89]
    AgeSD = [7.65, 10.09, 12.27, 12.03, 10.6, 10.6, 10.15, 10.15, 10.86, 10.86, 10.23, 10.23, 7.96, 7.96, 5.83]

    nbinoms = Vector{NegativeBinomial{Float64}}(undef, 15)
    for i = 1:15
        p = 1 - (AgeSD[i]^2-AgeMean[i])/(AgeSD[i]^2)
        r = AgeMean[i]^2/(AgeSD[i]^2-AgeMean[i])
        nbinoms[i] =  NegativeBinomial(r, p)
    end
    return nbinoms    
end
export negative_binomials


## references: 
# critical care capacity in Canada https://www.ncbi.nlm.nih.gov/pubmed/25888116
end # module end