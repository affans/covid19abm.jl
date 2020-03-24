module covid19abm
using Parameters, Distributions, StatsBase, StaticArrays, Random, Match, DataFrames

@enum HEALTH SUS LAT MILD MISO INF IISO HOS ICU REC DED UNDEF

Base.@kwdef mutable struct Human
    idx::Int64 = 0 
    health::HEALTH = SUS
    swap::HEALTH = UNDEF
    age::Int64   = 0    # in years. don't really need this but left it incase needed later
    ag::Int64    = 0
    tis::Int64   = 0   # time in state 
    exp::Int64   = 0    # max statetime
    iso::Bool = false  ## isolation
    icu::Bool = false ## is going to be in icu?
end

## default system parameters
@with_kw mutable struct ModelParameters @deftype Float64    ## use @with_kw from Parameters
    β = 0.0       
    prov::Symbol = :ontario 
    τmild::Int64 = 1 ## days before they self-isolate for mild cases
    fmild::Float64 = 0.05  ## percent of people practice self-isolation
    fsevere::Float64 = 0.80 # fixed at 0.80
    eldq::Float64 = 0.0 ## complete isolation of elderly
    calibration::Bool = false 
    modeltime::Int64 = 500    
    beta_timedepedent::Bool = false
end

Base.show(io::IO, ::MIME"text/plain", z::Human) = dump(z)

## constants 
const HSIZE = 10000
const humans = Array{Human}(undef, HSIZE) # run 100,000
const p = ModelParameters()  ## setup default parameters
const agebraks = @SVector [0:4, 5:19, 20:49, 50:64, 65:99]

export ModelParameters, HEALTH, Human, humans

function main(ip::ModelParameters)
    #Random.seed!(sim*726)
    ## datacollection            
    # matrix to collect model state for every time step

    # reset the parameters for the simulation scenario
    reset_params(ip)

    hmatrix = zeros(Int64, HSIZE, p.modeltime)
    initialize() # initialize population
    ags = [x.ag for x in humans] # store a vector of the age group distribution 
    # insert initial infected agents into the model
    # and setup the right swap function. 
    if p.calibration 
        # TO DO: checks whether mild, severe, and eldq parameters are set properly. 
        swapupdate = time_update_cal
        insert_infected(1, 4)  ## function will never isolation nor put in hospital/icu 
    else 
        swapupdate = time_update
        insert_infected(5, 4)  
    end    
    
    if p.beta_timedepedent
        tmp = td_beta.(1:p.modeltime)
        betavalues = (tmp .- minimum(tmp))./(maximum(tmp) .- 4*minimum(tmp));
    else
        betavalues = [p.β for _=1:p.modeltime]
    end
    
    # start the time loop
    for st = 1:p.modeltime
        # start of day
        _get_model_state(st, hmatrix) ## this datacollection needs to be at the start of the for loop
        totalinf = dyntrans(st, betavalues)  ## pass the transmission function time value.
        sw = swapupdate()
        # end of day
    end
    return hmatrix, ags ## return the model state as well as the age groups. 
end
export main

function reset_params(ip::ModelParameters)
    # the p is a global const
    # the ip is an incoming different instance of parameters 
    # copy the values from ip to p. 
    for x in propertynames(p)
        setfield!(p, x, getfield(ip, x))
    end
end
export reset_params

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
export _collectdf, _get_incidence_and_prev, _get_column_incidence, _get_column_prevalence

## initialization functions 
function get_province_ag(prov) 
    println("new func")
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
    @inbounds for i = 1:length(humans) 
        humans[i] = Human()              ## create an empty human       
        x = humans[i]
        x.idx = i 
        x.ag = rand(agedist)
        x.age = rand(agebraks[x.ag]) 
        x.exp = 999  ## susceptible people don't expire.
        if x.age >= 60  ## check if elderly need to be quarantined.
            if rand() < p.eldq
                x.iso = true
            end            
        end
    end
end
export initialize

function insert_infected(num, ag) 
    ## inserts a number of infected people in the population randomly
    ## this function should resemble move_to_inf()
    l = findall(x -> x.health == SUS && x.ag == ag, humans)
    if length(l) > 0 && num < length(l)
        h = sample(l, num; replace = false)
        @inbounds for i in h 
            x = humans[i]
            x.health = INF
            x.swap = REC
            x.tis = 0 
            x.iso = false 
            x.exp = 5  ## change if needed.  
            if rand() < p.fsevere    
                x.exp = 1  ## 1 day isolation for severe cases     
                x.swap = IISO
            end          
        end
    end    
    return h
end
export insert_infected

function time_update()
    # counters to calculate incidence
    lat=0; mild=0; miso=0; inf=0; infiso=0; hos=0; icu=0; rec=0; ded=0;
    for x in humans 
        x.tis += 1 
        if x.tis >= x.exp             
            @match Symbol(x.swap) begin
                :LAT  => begin move_to_latent(x); lat += 1; end
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
    end
    return (lat, mild, miso, inf, infiso, hos, icu, rec, ded)
end

function time_update_cal()
    ## calibration specific time_update function 
    ## keeps everyone in latent (so that we can count the number of secondary infections)
    a = 0
    for x in humans 
        x.tis += 1 
        if x.tis >= x.exp     
            move_to_latent(x)  ## after expiry from latent, they renew in latent.    
            a += 1           
        end
    end
    return (a, 0, 0, 0, 0, 0, 0, 0, 0) ## needs to return the same length tuple as time_update    
end
export time_update, time_update_cal

function move_to_latent(x::Human)
    ## transfers human h to the incubation period and samples the duration
    σ = LogNormal(log(5.2), 0.1) # duration of incubation period 
    θ = (0.8, 0.8, 0.8, 0.4, 0.2)  # percentage of sick individuals going to mild infection stage
    x.health = LAT
    x.swap = rand() < θ[x.ag] ? MILD : INF  # check whether person will INF or MILD after
    x.tis = 0   # reset time in state 
    x.exp = Int(round(rand(σ)))
    x.iso = false # starting at latent, person is not isolated from the community   
    x.icu = false # reset the icu so it dosn't pop up as a bug later.
end
export move_to_latent

function move_to_mild(x::Human)
    ## transfers human h to the mild infection stage for γ days
    γ = 5 # duration symptom onset to recovery, assumed fixed, based on serial interval... sampling creates a problem negative numbers
    x.health = MILD     
    x.tis = 0 
    x.exp = γ
    x.swap = REC 
    x.iso = false
    if rand() < p.fmild 
        x.swap = MISO  
        x.exp = p.τmild
    end
end
export move_to_mild

function move_to_miso(x::Human)
    ## transfers human h to the mild isolated infection stage for γ days
    γ = 5 # duration symptom onset to recovery, assumed fixed, based on serial interval... sampling creates a problem negative numbers
    oldexp = x.exp  ## tau amount of days spent in isolation as infectious already.
    x.health = MISO
    x.swap = REC
    x.tis = 0 
    x.exp = γ - oldexp  ## since tau amount of days was already spent as infectious but it dosn't really matter
    x.iso = true  ## self isolated from the community
end
export move_to_miso

function move_to_inf(x::Human)
    ## transfers human h to the mild isolated infection stage for γ days
    ## for swap, check if person will be hospitalized, selfiso, or recover
    ## if changing this function, also check if insert_infected needs to be changed. 

    # h = prob of hospital, c = prob of icu AFTER hospital    
    h = (rand(Uniform(0.02, 0.03)), rand(Uniform(0.02, 0.03)), rand(Uniform(0.28, 0.34)), rand(Uniform(0.28, 0.34)), rand(Uniform(0.60, 0.68)))
    c = (rand(Uniform(0.01, 0.015)), rand(Uniform(0.01, 0.015)), rand(Uniform(0.03, 0.05)), rand(Uniform(0.05, 0.1)), rand(Uniform(0.05, 0.15))) 
     
    δ = Int(round(rand(Uniform(2, 5)))) # duration symptom onset to hospitalization
    γ = 5 # duration symptom onset to recovery, assumed fixed, based on serial interval... sampling creates a problem negative numbers

    # death rate for severe cases.
    mh = [0.01/5, 0.01/5, 0.0135/3, 0.01225/1.5, 0.04/2]

    #mh = [0.01*3, 0.01*3, 0.0135*2, 0.01225, 0.03*2]
    
    x.health = INF
    x.swap = UNDEF
    x.tis = 0 
    x.iso = false # person is not isolated while infectious. 
    if rand() < h[x.ag]     # going to hospital or ICU but will spend delta time transmissing the disease with full contacts 
        x.exp = δ     
        x.swap = rand() < c[x.ag] ? ICU : HOS        
    else ## no hospital for this lucky (but severe) individual 
        if rand() < mh[x.ag]
            x.swap = DED
            x.iso = true  ## self isolated
            x.tis = 0     ## reset time in state 
            x.exp = γ   ## since tau amount of days was already spent as infectious 
        else 
            x.exp = γ  # as in the other functions.  
            x.swap = REC
            if rand() < p.fsevere 
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
    γ = 5 # duration symptom onset to recovery, assumed fixed, based on serial interval... sampling creates a problem negative numbers
    oldexp = x.exp
    x.health = IISO   
    x.swap = REC
    x.iso = true  ## self isolated
    x.tis = 0     ## reset time in state 
    x.exp = γ - oldexp  ## since tau amount of days was already spent as infectious 
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
    x.iso = true  ## hospitalized patients are isolated by default.

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
end

function move_to_recovered(h::Human)
    h.health = REC
    h.swap = UNDEF
    h.tis = 0 
    h.exp = 999 ## stay recovered indefinitely
    h.iso = false ## a recovered person has ability to meet others
end

function td_beta(t)
    ## returns a time-dependent beta value
    a0 = 6.261
    a1 = -11.81
    b1 = 1.817
    w = 0.01815
    temp =  a0 + a1*cos(t*w) + b1*sin(t*w)
    return temp
end


function dyntrans(t, betavalues)
    infs = findall(x -> x.health in (INF, MILD, MISO, IISO), humans)
    # if p.calibration 
    #     length(infs) > 1 && error("more than one infected person: length: $(length(infs))")
    # end
    totalinf = 0
    tomeet = map(1:length(agebraks)) do grp
        findall(x -> x.iso == false && x.ag == grp , humans)    
    end
    #length(tomeet) <= 5 && return totalinf
    for xid in infs
        x = humans[xid]
        ag = x.ag   
        if (x.health == MISO || x.health == IISO)
            cnt = rand(1:3)
        else 
            cnt = rand(nbs[ag])  ## get number of contacts/day
        end
        cnt >= length(tomeet[ag]) && error("error here")
        
        # distribute cnt_meet to different groups based on contact matrix. 
        # these are not probabilities, but proportions. be careful. 
        # going from cnt to the gpw array might remove a contact or two due to rounding. 
        gpw = Int.(round.(cm[ag]*cnt)) 
        #println("cnt: $cnt, gpw: $gpw")
        # let's stratify the human population in it's age groups. 
        # this could be optimized by putting it outside the contact_dynamic2 function and passed in as arguments               
        # enumerate over the 15 groups and randomly select contacts from each group
        for (i, g) in enumerate(gpw)
            if length(tomeet[i]) > 0 
                meet = rand(tomeet[i], g)    # sample 'g' number of people from this group  with replacement
            
                for j in meet
                    y = humans[j]
                    bf = betavalues[t]#p.β
                    if (x.health == MILD || x.health == MISO)
                        bf = bf * (1 - 0.5)  ## reduction factor 0.5
                    end
                    if y.health == SUS && rand() < bf
                        totalinf += 1
                        y.swap = LAT
                        y.exp = y.tis ## force the move to latent in the next time step.
                    end  
                end
            end            
        end
    end
    return totalinf
end
export dyntrans

## old contact matrix
# function contact_matrix()
#     CM = Array{Array{Float64, 1}, 1}(undef, 4)
#     CM[1]=[0.5712, 0.3214, 0.0722, 0.0353]
#     CM[2]=[0.1830, 0.6253, 0.1423, 0.0494]
#     CM[3]=[0.1336, 0.4867, 0.2723, 0.1074]    
#     CM[4]=[0.1290, 0.4071, 0.2193, 0.2446]
#     return CM
# end

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