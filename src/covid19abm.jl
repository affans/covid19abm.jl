module covid19abm
using Parameters, Distributions, StatsBase, StaticArrays, Random

@enum HEALTH UNDEF=0 SUSC=1 EXP=2 INFMILD=4 INFHOSP=5 INFICU=6
Base.@kwdef mutable struct Human
    idx::Int64 = 0 
    health::HEALTH = SUSC
    swap::HEALTH = UNDEF
    age::Int64   = 0    # in years. 
    ag::Int64    = 0
    tis::Int64   = 0   # time in state 
    exp::Int64   = 0    # max statetime
    iso_time::Int64 = 0 
    is_iso::Bool = false
end

@with_kw mutable struct ModelParameters @deftype Float64    ## use @with_kw from Parameters
    st::Int64 = 500 ## sim time in days
    β = 0.0 
    I₀::Int64 = 1    ## initial infected
    σ::Distribution = LogNormal(log(5.2), 0.1)  ## distribution for incubation period
    γ::Int64 = 5 ## 5 days infectious period
    q::Int64 = 0.05 ## percentage of immediate self-isolation
end
Base.show(io::IO, ::MIME"text/plain", z::Human) = dump(z)

## constants
const agedist =  Categorical(@SVector [0.190674655, 0.41015843, 0.224428346, 0.17473857])
const agebraks = @SVector [0:19, 20:49, 50:64, 65:99]
const humans = Array{Human}(undef, 1000)
const p = ModelParameters()
export ModelParameters, HEALTH, Human, humans

## initialization functions 
function initialize() 
    @inbounds for i = 1:length(humans) 
        x = humans[i]
        x = Human()              ## create an empty human
        x.idx = i 
        x.ag = rand(agedist)
        x.age = rand(agebraks[x.ag])     
    end
end
export initialize

function insert_exposed(num) 
    ## inserts a number of infected people in the population randomly
    l = findall(x -> x.health == SUSC, humans)
    if length(l) > 0 
        h = sample(l, num; replace = false)
        @inbounds for i in h 
            move_to_latent(humans[i])            
        end
    end    
end
export insert_exposed

function move_to_latent(h::Human)
    ## make the i'th human exposed
    h.health = EXP
    h.swap = UNDEF
    h.exp = Int(ceil(rand(p.σ)))
    h.tis = 0    
end

function move_to_infected(h::Human)
    ## to do: determine if infection is mild, hosp, or ICU
    h.health = INFMILD
    h.swap = UNDEF 
    h.exp = p.γ
    h.tis = 0    
    
    ## determine qurantine levels.
    # determine, f, tau, hosp, icu.
    error("not implemented")
    
end

## references: 
# critical care capacity in Canada https://www.ncbi.nlm.nih.gov/pubmed/25888116
