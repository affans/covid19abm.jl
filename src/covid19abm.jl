module covid19abm
using Parameters, Distributions, StatsBase, StaticArrays, Random, Match, DataFrames

@enum HEALTH SUS LAT PRE ASYMP MILD MISO INF IISO HOS ICU REC DED UNDEF

Base.@kwdef mutable struct Human
    idx::Int64 = 0 
    health::HEALTH = SUS
    swap::HEALTH = UNDEF
    sickfrom::HEALTH = UNDEF
    age::Int64   = 0    # in years. don't really need this but left it incase needed later
    ag::Int64    = 0
    tis::Int64   = 0   # time in state 
    exp::Int64   = 0    # max statetime
    iso::Bool = false  ## isolated (limited contacts)
    icu::Bool = false ## is going to be in icu?
    bhpcat::Symbol = :none
end

## default system parameters
@with_kw mutable struct ModelParameters @deftype Float64    ## use @with_kw from Parameters
    β = 0.0       
    prov::Symbol = :ontario 
    calibration::Bool = false 
    modeltime::Int64 = 500
    initialinf::Int64 = 1
    bhpinitialcat::Symbol = :mp ## the initial person is introduced in this category
    τmild::Int64 = 0 ## days before they self-isolate for mild cases
    fmild::Float64 = 0.0  ## percent of people practice self-isolation
    fsevere::Float64 = 0.0 # fixed at 0.80
    eldq::Float64 = 0.0 ## complete isolation of elderly
    fasymp::Float64 = 0.50 ## percent going to asymp
    fpre::Float64 = 1.0 ## percent going to presymptomatic
    fpreiso::Float64 = 0.0 ## percent that is isolated at the presymptomatic stage
    tpreiso::Int64 = 0## preiso is only turned on at this time. 
end

Base.show(io::IO, ::MIME"text/plain", z::Human) = dump(z)

## constants 
const HSIZE = 816
const humans = Array{Human}(undef, HSIZE) # run 100,000
const p = ModelParameters()  ## setup default parameters
const agebraks = @SVector [0:4, 5:19, 20:49, 50:64, 65:99]
const bhpdist = @SVector [8, 164, 149, 34, 34, 161, 80, 80, 60, 9, 19, 18]
const bhpcats = @SVector [:c1, :c2, :c3, :c4, :c5, :c6, :c7, :c8, :c9, :c10, :c11, :c12]

export ModelParameters, HEALTH, Human, humans

function runsim(simnum, ip::ModelParameters)
    # function runs the `main` function, and collects the data as dataframes. 
    hmatrix = main(ip)            
    # get infectors counters
    infectors = _count_infectors()
    # get simulation age groups
    ags = [x.ag for x in humans] # store a vector of the age group distribution 

    all = _collectdf(hmatrix)
    spl = _splitstate(hmatrix, ags)
    g1 = _collectdf(spl[1])
    g2 = _collectdf(spl[2])
    g3 = _collectdf(spl[3])
    g4 = _collectdf(spl[4])
    g5 = _collectdf(spl[5])
    g6 = _collectdf(spl[6])
    g7 = _collectdf(spl[7])
    g8 = _collectdf(spl[8])
    g9 = _collectdf(spl[9])
    g10 = _collectdf(spl[10])
    g11 = _collectdf(spl[11])
    g12 = _collectdf(spl[12])
    insertcols!(all, 1, :sim => simnum);    
    insertcols!(g1, 1, :sim => simnum); 
    insertcols!(g2, 1, :sim => simnum); 
    insertcols!(g3, 1, :sim => simnum); 
    insertcols!(g4, 1, :sim => simnum); 
    insertcols!(g5, 1, :sim => simnum); 
    insertcols!(g6, 1, :sim => simnum); 
    insertcols!(g7, 1, :sim => simnum); 
    insertcols!(g8, 1, :sim => simnum); 
    insertcols!(g9, 1, :sim => simnum); 
    insertcols!(g10, 1, :sim => simnum); 
    insertcols!(g11, 1, :sim => simnum); 
    insertcols!(g12, 1, :sim => simnum); 
    
    return (a=all, g1=g1, g2=g2, g3=g3, g4=g4, g5=g5, g6=g6, g7=g7, g8=g8, g9=g9, g10=g10, g11=g11, g12=g12, infectors=infectors)
end
export runsim

function main(ip::ModelParameters)
    #Random.seed!(sim*726)
    ## datacollection            
    # matrix to collect model state for every time step

    # reset the parameters for the simulation scenario
    reset_params(ip)

    hmatrix = zeros(Int64, HSIZE, p.modeltime)
    initialize() # initialize population
    # insert initial infected agents into the model
    # and setup the right swap function. 
    # [:mp, :po, :im, :en, :ma, :pr, :si, :wa]
    if p.calibration 
        insert_infected(PRE, p.initialinf, rand(bhpcats)) ## for calibration, insert in random 
    else 
        insert_infected(LAT, p.initialinf, p.bhpinitialcat) ## for main simulations, insert in category
    end    
    
    ## save the preisolation isolation parameters
    _fpreiso = p.fpreiso
    p.fpreiso = 0

    # start the time loop
    for st = 1:p.modeltime
        # start of day
        if st == p.tpreiso ## time to introduce testing
            p.fpreiso = _fpreiso
        end
        _get_model_state(st, hmatrix) ## this datacollection needs to be at the start of the for loop
        dyntrans()
        sw = time_update()
        # end of day
    end
    return hmatrix ## return the model state as well as the age groups. 
end
export main

reset_params_default() = reset_params(ModelParameters())
function reset_params(ip::ModelParameters)
    # the p is a global const
    # the ip is an incoming different instance of parameters 
    # copy the values from ip to p. 
    for x in propertynames(p)
        setfield!(p, x, getfield(ip, x))
    end
end
export reset_params, reset_params_default

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
    insertcols!(datf, 1, :time => 1:p.modeltime) ## add a time column to the resulting dataframe
    return datf
end

function _splitstate(hmatrix, ags)
    #split the full hmatrix into 4 age groups based on ags (the array of age group of each agent)
    #sizes = [length(findall(x -> x == i, ags)) for i = 1:4]
    matx = []#Array{Array{Int64, 2}, 1}(undef, 4)
    for ct in bhpcats
        idx = findall(x -> x.bhpcat == ct, humans)
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
    # alternative: vcat(collect.(Iterators.repeated.(bhpcats,bhpdist))...)
    bhpsymbols = vcat([[bhpcats[i] for j = 1:bhpdist[i]] for i = 1:length(bhpdist)]...)
    shuffle!(bhpsymbols)
    agedist = get_province_ag(p.prov)
    @inbounds for i = 1:HSIZE 
        humans[i] = Human()              ## create an empty human       
        x = humans[i]
        x.idx = i 
        x.ag = rand(agedist)
        x.age = rand(agebraks[x.ag]) 
        x.exp = 999  ## susceptible people don't expire.
        if rand() < p.eldq && x.age >= 60  ## check if elderly need to be quarantined.
            x.iso = true            
        end
        ## for bhp: need to assign according to vector
        x.bhpcat = bhpsymbols[i]        
    end
end
export initialize

function insert_infected(health, num, ag) 
    ## inserts a number of infected people in the population randomly
    ## this function should resemble move_to_inf()
    l = findall(x -> x.health == SUS && x.bhpcat == ag, humans)
    if length(l) > 0 && num < length(l)
        h = sample(l, num; replace = false)
        @inbounds for i in h 
            x = humans[i]
            if health == PRE 
                move_to_pre(x) ## the swap may be asymp, mild, or severe, but we can force severe in the time_update function
            elseif health == LAT 
                move_to_latent(x)
            elseif health == INF
                move_to_infsimple(x)
            else 
                error("can not insert human of health $(health)")
            end       
            x.sickfrom = INF # this will add +1 to the INF count in _count_infectors()... keeps the logic simple in that function.    
        end
    end    
    return h
end
export insert_infected

function time_update()
    # counters to calculate incidence
    lat=0; pre=0; asymp=0; mild=0; miso=0; inf=0; infiso=0; hos=0; icu=0; rec=0; ded=0;
    for x in humans 
        x.tis += 1 
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
    end
    return (lat, mild, miso, inf, infiso, hos, icu, rec, ded)
end
export time_update

function move_to_latent(x::Human)
    ## transfers human h to the incubation period and samples the duration
    σ = LogNormal(log(5.2), 0.1) # duration of incubation period 
    θ = (0.8, 0.8, 0.8, 0.4, 0.2)  # percentage of sick individuals going to mild infection stage
    x.health = LAT
    x.tis = 0   # reset time in state 

    ## set up swap
    if rand() < p.fpre
        x.swap = PRE
        x.exp = Int(round(rand(σ))) - 1
    else 
        if rand() < θ[x.ag]
            if rand() < p.fasymp 
                x.swap = ASYMP
            else 
                x.swap = MILD 
            end        
        else 
            x.swap = INF
        end
        x.exp = Int(round(rand(σ)))
    end
    if p.calibration ## in calibration mode, latent people never become infectious.
        x.swap = LAT 
        x.exp = 999
    end
    x.iso == true && error("how did isolated person become infected")
    x.iso = false # starting at latent, person is not isolated from the community   
    x.icu = false # reset the icu so it dosn't pop up as a bug later.
end
export move_to_latent

function move_to_pre(x::Human)
    θ = (0.8, 0.8, 0.8, 0.4, 0.2)  # percentage of sick individuals going to mild infection stage
    x.health = PRE
    x.tis = 0   # reset time in state 
    x.exp = 1   # one day in presymptomatic

    # calculate whether person is isolated
    if rand() < p.fpreiso
        x.iso = true # starting at latent, person is not isolated from the community   
    else 
        x.iso = false
    end    

    # the reason why we can't move to MILD/MISO and INF/IISO here is that, 
    # particularly, for inf we need to check for hospitalization
    if rand() < θ[x.ag]
        if rand() < p.fasymp 
            x.swap = ASYMP
        else 
            x.swap = MILD 
        end        
    else 
        x.swap = INF
    end

end
export move_to_pre

function move_to_asymp(x::Human)
     ## transfers human h to the asymptomatic stage 
     γ = 5 
     x.health = ASYMP     
     x.tis = 0 
     x.exp = γ
     x.swap = REC 
     # x.iso property remains from either the latent or presymptomatic class
     # if x.iso is true, the asymptomatic individual has limited contacts
end
export move_to_asymp

function move_to_mild(x::Human)
    ## transfers human h to the mild infection stage for γ days
    γ = 5 # duration symptom onset to recovery, assumed fixed, based on serial interval... sampling creates a problem negative numbers
    x.health = MILD     
    x.tis = 0 
    x.exp = γ
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
    γ = 5 # duration symptom onset to recovery, assumed fixed, based on serial interval... sampling creates a problem negative numbers
    oldexp = x.exp  ## tau amount of days spent in isolation as infectious already.
    x.health = MISO
    x.swap = REC
    x.tis = 0 
    x.exp = γ - oldexp  ## since tau amount of days was already spent as infectious but it dosn't really matter
    x.iso = true  ## isolated from the community ... may already be set from before. still need it here if not set before
end
export move_to_miso

function move_to_infsimple(x::Human)
    ## transfers human h to the severe infection stage for γ days 
    ## simplified function for calibration/general purposes
    x.health = INF
    x.swap = REC
    x.tis = 0 
    x.iso = false 
    x.exp = 5  ## change if needed.   
end

function move_to_inf(x::Human)
    ## transfers human h to the severe infection stage for γ days
    ## for swap, check if person will be hospitalized, selfiso, die, or recover
 
    # h = prob of hospital, c = prob of icu AFTER hospital    
    h = (rand(Uniform(0.02, 0.03)), rand(Uniform(0.02, 0.03)), rand(Uniform(0.28, 0.34)), rand(Uniform(0.28, 0.34)), rand(Uniform(0.60, 0.68)))
    c = (rand(Uniform(0.01, 0.015)), rand(Uniform(0.01, 0.015)), rand(Uniform(0.03, 0.05)), rand(Uniform(0.05, 0.1)), rand(Uniform(0.05, 0.15))) 
    mh = [0.01/5, 0.01/5, 0.0135/3, 0.01225/1.5, 0.04/2]     # death rate for severe cases.
    
    if p.calibration
        h =  (0, 0, 0, 0, 0)
        c =  (0, 0, 0, 0, 0)
        mh = (0, 0, 0, 0, 0)
    end

    δ = Int(round(rand(Uniform(2, 5)))) # duration symptom onset to hospitalization
    γ = 5 # duration symptom onset to recovery, assumed fixed, based on serial interval... sampling creates a problem negative numbers

    x.health = INF
    x.swap = UNDEF
    x.tis = 0 
    if rand() < h[x.ag]     # going to hospital or ICU but will spend delta time transmissing the disease with full contacts 
        x.exp = δ     
        x.swap = rand() < c[x.ag] ? ICU : HOS        
    else ## no hospital for this lucky (but severe) individual 
        if rand() < mh[x.ag]
            x.exp = γ   
            x.swap = DED
        else 
            x.exp = γ  # as in the other functions.  
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
    γ = 5 # duration symptom onset to recovery, assumed fixed, based on serial interval... sampling creates a problem negative numbers
    oldexp = x.exp
    x.health = IISO   
    x.swap = REC
    x.iso = true  ## isolated from the community ... may already be set from before.
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

function dyntrans()
    totalinf = 0
    ## find all the people infectious
    infs = findall(x -> x.health in (PRE, ASYMP, MILD, MISO, INF, IISO), humans)
    # tomeet = map(1:length(agebraks)) do grp ## will also meet dead people, but ignore for now because it's such a small group
    #     findall(x -> x.ag == grp , humans)    
    # end
    tomeet = map(bhpcats) do y
        findall(x -> x.bhpcat == y, humans)
    end ## the order of tomeet is according to the vector bhpcats 
    ## which is different than the order in the contact distribution dicts, but it's okay since 
    ## we get the proper key/value from the dict (where the value IS in order... ie nbs[:mp] gives the number of contacts in order of bhpcats)
    #length(tomeet) <= 5 && return totalinf
    nbss = portteams_binomials()
    for xid in infs
        x = humans[xid]
        ag = x.bhpcat  
        ih = x.health
        if x.iso ## isolated infectious person has limited contacts
            cnt = rand(1:3)
        else 
            cnt = rand(nbss[ag])  ## get number of contacts/day
        end
        #cnt >= length(tomeet[ag]) && error("error here")
        
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
                ## remember, two infected people may meet the same susceptible so its possible that disease can be transferred "twice"
                for j in meet
                    y = humans[j]
                    if y.health == SUS && y.swap == UNDEF
                        bf = p.β ## baseline PRE
                        # values coming from FRASER Figure 2... relative tranmissibilities of different stages.
                        if ih == ASYMP
                            bf = bf * 0.11
                        elseif ih == MILD || ih == MISO 
                            bf = bf * 0.44
                        elseif ih == INF || ih == IISO 
                            bf = bf * 0.89
                        end
                        if rand() < bf
                            totalinf += 1
                            y.swap = LAT
                            y.exp = y.tis   ## force the move to latent in the next time step.
                            y.sickfrom = ih ## stores the infector's status to the infectee's sickfrom
                        end  
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

function bhp_matrix()
    #  0-4, 5-19, 20-49, 50-64, 65+
    dd = Dict{Symbol, Array{Float64, 1}}()
    CM = Array{Array{Float64, 1}, 1}(undef, 8)    
    CM[1] = [0.6444, 0.1, 0.0185, 0.0185, 0.1, 0.1, 0.0185, 0]
    CM[2] = [0.0063, 0.7704, 0.0012, 0.0023, 0.1481, 0.0469, 0.0012, 0.0235]
    CM[3] = [0.0833, 0.0833, 0.5, 0.0833, 0.0833, 0.0833, 0.0833, 0]
    CM[4] = [0.05, 0.05, 0.025, 0.75, 0.075, 0.025, 0.025, 0]
    CM[5] = [0.0063, 0.1467, 0.0012, 0.0023, 0.6844, 0.1114, 0.0012, 0.0465]
    CM[6] = [0.0112, 0.0826, 0.0021, 0.0041, 0.1979, 0.6587, 0.0021, 0.0413]
    CM[7] = [0.0312, 0.0625, 0.0312, 0.0312, 0.0625, 0.0625, 0.7188, 0]
    CM[8] = [0, 0.1242, 0, 0.0062, 0.4969, 0.1242, 0, 0.2484]
    for i = 1:8
        push!(dd, bhpcats[i] => CM[i])
    end 
    return dd
end


function portteam_matrix() 
    dd = Dict{Symbol, Array{Float64, 1}}()
    CM = Array{Array{Float64, 1}, 1}(undef, length(bhpcats))    
    CM[1] = [0.4706, 0.1059, 0.1059, 0.0471, 0.0235, 0.0941, 0.0235, 0.0235, 0.0235, 0.0235, 0.0294, 0.0294]
    CM[2] = [0.004, 0.874, 0.0002, 0.0132, 0, 0.0612, 0.0092, 0, 0.0193, 0, 0.0189, 0]
    CM[3] = [0.0043, 0.0002, 0.8963, 0, 0.0142, 0.034, 0, 0.0099, 0.0208, 0, 0, 0.0203]
    CM[4] = [0.0173, 0.1296, 0, 0.6479, 0, 0.108, 0.0259, 0, 0.013, 0, 0.0583, 0]
    CM[5] = [0.0118, 0, 0.177, 0, 0.59, 0.0885, 0, 0.0354, 0.0177, 0, 0, 0.0796]
    CM[6] = [0.0072, 0.1248, 0.0644, 0.0224, 0.0134, 0.6783, 0.026, 0.026, 0.034, 0.0027, 0.0004, 0.0004]
    CM[7] = [0.0048, 0.0504, 0, 0.0144, 0, 0.0696, 0.8391, 0, 0.0072, 0, 0.0144, 0]
    CM[8] = [0.0048, 0, 0.0504, 0, 0.0144, 0.0696, 0, 0.8391, 0.0072, 0, 0, 0.0144]
    CM[9] = [0.0029, 0.064, 0.064, 0.0044, 0.0044, 0.0552, 0.0044, 0.0044, 0.7922, 0.0044, 0, 0]
    CM[10] = [0.0357, 0, 0, 0, 0, 0.0536, 0, 0, 0.0536, 0.8571, 0, 0]
    CM[11] = [0.0182, 0.3127, 0, 0.0982, 0, 0.0036, 0.0436, 0, 0, 0, 0.5236, 0]
    CM[12] = [0.0119, 0, 0.2043, 0, 0.0641, 0.0024, 0, 0.0285, 0, 0, 0, 0.6888]
    for i = 1:length(bhpcats)
        push!(dd, bhpcats[i] => CM[i])
    end 
    return dd
end
const cm = portteam_matrix()


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

function bhp_binomials()
    dd = Dict{Symbol, DiscreteUniform}()
    means = [12, 12, 9, 9, 12, 9, 11, 50]
    lb = [11, 11, 7, 7, 11, 8, 9, 45]
    hb = [13, 13, 10, 10, 13, 10, 12, 55]
    for i = 1:8
        push!(dd, bhpcats[i] => DiscreteUniform(lb[i], hb[i]))
    end
    return dd
end
const nbs = bhp_binomials()
export negative_binomials, contact_matrix, nbs, cm

function portteams_binomials()
    dd = Dict{Symbol, DiscreteUniform}()
    #means = [12, 12, 9, 9, 12, 9, 11, 50]

    lb = [20, 27, 27, 13, 9, 13, 9, 9, 22, 12, 14, 22]
    hb = [22, 29, 29, 15, 11, 15, 11, 11, 24, 14, 16, 24]
    for i = 1:length(bhpcats)
        push!(dd, bhpcats[i] => DiscreteUniform(lb[i], hb[i]))
    end
    return dd
end


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