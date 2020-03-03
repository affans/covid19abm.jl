module covid19abm
using Parameters, Distributions, StatsBase, StaticArrays, Random

@enum HEALTH UNDEF=0 SUSC=1 INF=2 REC=3
mutable struct Human
    idx::Int64
    health::HEALTH
    swap::HEALTH       
    age::Int64      # in years. 
    tis::Int64      # time in state 
    exp::Int64       # max statetime
    Human() = new(0, UNDEF, UNDEF, 0, 0, 0)
end

@with_kw mutable struct ModelParameters @deftype Float64    
    st::Int64 = 500 ## sim time in days
    β = 0.0 
    ii::Int64 = 20    ## initial infected
end
export ModelParameters, Human

end # module
