module covid19abm
using Parameters, Distributions, StatsBase, StaticArrays, Random, Match, DataFrames

@enum HEALTH UNDEF SUS LAT INF INFISO HOS REC DED
Base.@kwdef mutable struct Human
    idx::Int64 = 0 
    health::HEALTH = SUS
    swap::HEALTH = UNDEF
    age::Int64   = 0    # in years. 
    ag::Int64    = 0
    tis::Int64   = 0   # time in state 
    exp::Int64   = 0    # max statetime
    iso::Bool = false  ## isolation
    icu::Bool = false ## is going to be in icu?
end

@with_kw mutable struct ModelParameters @deftype Float64    ## use @with_kw from Parameters
    β = 0.0 
    I₀::Int64 = 1    ## initial infected
    prov::Symbol = :ontario 
    #σ::Distribution = LogNormal(log(5.2), 0.1)  ## distribution for incubation period
    σ::Int64 = 6
    γ::Int64 = 5 ## 5 days infectious period for non hospitalized
    δ::Int64 = 4 ## 4 days infectious period for those going to hospital
    τ::Int64 = 1 ## days before they self-isolate 
    h::Float64 = 0.025 ## proporiton requiring hospital
    c::Float64 = 0.07  ## proporiton requiring ICU
    f::Float64 = 0.05  ## percent of people practice self-isolation
    psiC::Int64 = 13   ## days to recover from ICU
    psiH::Int64 = 10   ## days to recover from hospitalization
    muC::Int64 = 0     ## days to death from hospital (if death happening)
    muH::Int64 = 0     ## days to death from hospital (if death happening)
    mh::Float64 = 0.0  ## probability of death in hospital
    mc::Float64 = 0.0  ## probability of death in hospital
    ## internal parameters
    calibration::Bool = false 
end

Base.show(io::IO, ::MIME"text/plain", z::Human) = dump(z)

## constants
const humans = Array{Human}(undef, 10000)
const p = ModelParameters()
const maxtime = 500
const agebraks = @SVector [0:19, 20:49, 50:64, 65:99]

export ModelParameters, HEALTH, Human, humans

function main(sim)
    #Random.seed!(sim*726)
    modeltime = 500 #maxtime

    ## datacollection 
    _names = Symbol.(["sus", "lat", "inf", "iso", "infiso", "hos", "icu", "ded", "rec"])
    #prev = DataFrame([Int64 for i = 1:length(_names)], _names, maxtime) 
    dat = DataFrame([zeros(Int64, modeltime) for i = 1:length(_names)], _names) 
    initialize() # initialize population
    e = insert_exposed()

    swapupdate = time_update
    if p.calibration 
        move_to_infected(humans[e[1]])
        swapupdate = time_update_cal
    end
    dump(p)
    #println("swapupdate func: $(swapupdate)")
    #_modelstate() 
    for st = 1:modeltime
        # start of day
        collect_data(st, dat)
        dyntrans()
        swapupdate()
        # end of day
    end
    #_modelstate() 
    return dat
end
export main

function reset_params(ip::ModelParameters)
    # the p is a global const
    # the ip is an incoming DIFFERENT instance of parameters 
    # copy the values. 
    for x in propertynames(p)
        setfield!(p, x, getfield(ip, x))
    end
end
export reset_params

function _modelstate() 
    ## prints a string with model data
    hstr = ""
    sstr = ""
    for y in instances(HEALTH)        
        hcnt = length(findall(x -> x.health == y , humans))
        scnt = length(findall(x -> x.swap == y , humans))
        hstr = string(hstr, "$(string(y)): $hcnt ") 
        sstr = string(sstr, "$(string(y)): $scnt ") 
    end
    println("class:$hstr \nswaps: $sstr")
end
export _modelstate

function collect_data(t, dat)
    dat[t, :sus] = length(findall(x -> x.health == SUS, humans))
    dat[t, :lat] = length(findall(x -> x.health == LAT, humans))
    dat[t, :inf] = length(findall(x -> x.health == INF, humans))
    dat[t, :iso] = length(findall(x -> x.health == INFISO, humans))
    dat[t, :hos] = length(findall(x -> x.health == HOS && x.icu == false, humans))
    dat[t, :icu] = length(findall(x -> x.health == HOS && x.icu == true, humans))
    dat[t, :ded] = length(findall(x -> x.health == DED, humans))
    dat[t, :rec] = length(findall(x -> x.health == REC, humans))
end
export collect_data

## initialization functions 
function get_province_ag(prov) 
    ret = @match prov begin
        :ontario   => Categorical(@SVector [0.190674655, 0.41015843, 0.224428346, 0.17473857])   
        :alberta   => Categorical(@SVector [0.209225217, 0.457056611, 0.203940396, 0.129777775])
        :bc        => Categorical(@SVector [0.172844215, 0.406074014, 0.231167323, 0.189914448])
        :manitoba  => Categorical(@SVector [0.215274172, 0.410804507, 0.209933322, 0.163988])
        :newbruns  => Categorical(@SVector [0.171227734, 0.370378096, 0.251498338, 0.206895832])
        :newfdland => Categorical(@SVector [0.166635101, 0.377318807, 0.254632629, 0.201413463])
        :nwterrito => Categorical(@SVector [0.231552163, 0.479643766, 0.206870229, 0.081933842])
        :novasco   => Categorical(@SVector [0.169896365, 0.373800103, 0.249685273, 0.206618259])
        :nunavut   => Categorical(@SVector [0.35, 0.478, 0.13, .042])
        :pei       => Categorical(@SVector [0.187067395, 0.36856102, 0.242440801, 0.201930783])
        :quebec    => Categorical(@SVector [0.179873623, 0.395567459, 0.232995696, 0.191563221])
        :saskat    => Categorical(@SVector [0.216249874, 0.408915409, 0.210946403, 0.163888315])
        :yukon     => Categorical(@SVector [0.19054588, 0.43875311, 0.246012001, 0.124689009])
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
        x.exp = maxtime + 1  ## susceptible people don't expire. 
    end
end
export initialize

insert_exposed() = insert_exposed(p.I₀)
function insert_exposed(num) 
    ## inserts a number of infected people in the population randomly
    l = findall(x -> x.health == SUS, humans)
    if length(l) > 0 && num < length(l)
        h = sample(l, num; replace = false)
        @inbounds for i in h 
            move_to_latent(humans[i])            
        end
    end    
    return h
end
export insert_exposed

function time_update()
    for x in humans 
        x.tis += 1 
        if x.tis > x.exp             
            @match Symbol(x.swap) begin
                :LAT  => move_to_latent(x)                
                :INF  => move_to_infected(x)    
                :INFISO => move_to_infiso(x)            
                :HOS  => move_to_hospicu(x)
                :REC  => move_to_recovered(x)
                :DED  => move_to_dead(x)
                _    => error("swap expired, but no swap set.")
            end
        end
    end
end

function time_update_cal()
    ## calibration specific time_update function 
    ## keeps everyone in latent (so that we can count the number of secondary infections)
    for x in humans 
        x.tis += 1 
        if x.tis > x.exp     
            move_to_latent(x)  ## after expiry from latent, they renew in latent.               
        end
    end
end
export time_update, time_update_cal, move_to_latent, move_to_hospicu, move_to_infected, move_to_recovered, move_to_dead

function move_to_latent(h::Human)
    ## human h is now in incubation period.
    h.health = LAT
    h.swap = INF
    h.tis = 0    ## reset time in state 
    h.exp = p.σ    
    h.iso = false    
    h.icu = false # reset the icu so it dosn't pop up as a bug later.
end

function move_to_infected(h::Human)
    ## human h is now in infectious period.
    ## for swap, check if person will be hospitalized, selfiso, or recover
    h.health = INF
    h.swap = UNDEF
    h.tis = 0    ## reset time in state 
    h.iso = false ## person is not isolated while infectious. 
    if rand() < p.h     ## going to hospital but will spend time transmissing the disease 
        h.exp = p.δ            
        h.swap = HOS  ## per will go to hospital after delta amount of time.       
    else ## no hospital for this lucky individual ... will be transmitting to others. 
        h.exp = p.γ
        h.swap = REC
        if rand() < p.f      
            h.exp = p.τ      
            h.swap = INFISO
        end
    end
    ## before returning, check if swap is set 
    h.swap == UNDEF && error("agent I -> ?")
end

function move_to_infiso(h::Human)
    ## human h is now in self-isolation period.
    h.health = INFISO 
    h.swap = REC
    h.iso = true  ## self isolated
    h.tis = 0     ## reset time in state 
    h.exp = p.γ ##p.γ - p.τ  ## since tau amount of days was already spent as infectious 
    ## but who cares.. self-isolated people will eventually recover and have no input in dynamics. 
end 

function move_to_hospicu(h::Human)
    # check for death, rec, icu here here. 
    h.health = HOS
    h.swap = UNDEF
    h.tis = 0
    h.iso = true  ## hospitalized patients are isolated by default.
    ## check for icu 
    ## Q: why is stay in ICU shorter than stay in hospital. 
    ## we should just have a distribution that assigns days in hospital whether ICU or not. 
    if rand() < p.c  
        h.icu = true
        if rand() < p.mc ## person will die in the ICU 
            h.exp = p.muC
            h.swap = DED
        else 
            h.exp = p.psiC
            h.swap = REC
        end
    else 
        h.icu = false
        ## person stays in hospital
        if rand() < p.mh ## person will die in the hospital 
            h.exp = p.muH 
            h.swap = DED
        else 
            h.exp = p.psiH 
            h.swap = REC
        end        
    end    
     ## before returning, check if swap is set 
     h.swap == UNDEF && error("agent H -> ?")
end

function move_to_dead(h::Human)
    # no level of alchemy will bring someone back to life. 
    h.health = DED
    h.swap = UNDEF
    h.tis = 0 
    h.exp = maxtime + 1 ## stay recovered indefinitely
    h.iso = true # a dead person is isolated
end

function move_to_recovered(h::Human)
    h.health = REC
    h.swap = UNDEF
    h.tis = 0 
    h.exp = maxtime + 1 ## stay recovered indefinitely
    h.iso = false ## a recovered person has ability to meet others
end


function dyntrans()
    infs = findall(x -> x.health == INF && x.iso == false, humans)
    tomeet = findall(x -> x.iso == false, humans) #sus, lat, inf, rec 
    totalinf = 0
    
    tomeet = map(1:4) do grp
        findall(x -> x.iso == false && x.ag == grp , humans)    
    end
    #length(tomeet) <= 5 && return totalinf
    for xid in infs
        x = humans[xid]
        ag = x.ag   
        cnt = rand(nbs[ag])  ## get number of contacts/day
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
            meet = rand(tomeet[i], g)    # sample 'g' number of people from this group 
            #println("... meeting grp: $g, meet: $meet")
            for j in meet
                y = humans[j]
                if y.health == SUS && rand() < p.β
                    totalinf += 1
                    y.swap = LAT
                    y.exp = y.tis ## force the move to latent in the next time step.
                end  
            end
        end
    end
    return totalinf
end
export dyntrans

function contact_matrix()
    CM = Array{Array{Float64, 1}, 1}(undef, 4)
    CM[1]=[0.5712, 0.3214, 0.0722, 0.0353]
    CM[2]=[0.1830, 0.6253, 0.1423, 0.0494]
    CM[3]=[0.1336, 0.4867, 0.2723, 0.1074]    
    CM[4]=[0.1290, 0.4071, 0.2193, 0.2446]
    return CM
end
function negative_binomials() 
    ## the means/sd here are calculated using _calc_avgag
    means = [15.30295, 13.7950, 11.2669, 8.0027]
    sd = [11.1901, 10.5045, 9.5935, 6.9638]
    nbinoms = Vector{NegativeBinomial{Float64}}(undef, 4)
    for i = 1:4
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

## deprecated function 
function transmission() 
    ## basic homogeneous transmission dynamics for testing purposes
    ## each infected person meets 5 contacts from the susceptible pool. 
    ## beta value is fixed
    infs = findall(x -> x.health == INF, humans)
    sus = findall(x -> x.health == SUS, humans)
    length(sus) <= 5 && return 
    for i in infs 
        x = humans[i] 
        cnt = rand(1:5)         
        smeet = sample(sus, cnt; replace = false)  ## will throw error if not enough susceptibles. 
        for j in smeet 
            if rand() < 0.90 
                humans[j].swap = LAT
                humans[j].exp = humans[j].tis ## force the move to latent in the next time step.
            end
        end
    end
end
export transmission

## references: 
# critical care capacity in Canada https://www.ncbi.nlm.nih.gov/pubmed/25888116
end # module end