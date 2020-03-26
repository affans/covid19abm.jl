using covid19abm
const cv = covid19abm
using Test

@testset "parameters" begin 
    ip = ModelParameters()
    mod_ip = cv.p ## get the module parameters    
    ## set random parameters 
    ip.β = 99
    ip.prov = :noknownlocation 
    ip.τmild = 1 ## days before they self-isolate for mild cases
    ip.fmild = 0.05  ## percent of people practice self-isolation
    ip.fsevere = 0.80 # fixed at 0.80, always within 1 day. 
    ip.eldq = 0.0 ## percentage quarantine of 60+ individuals, removes from chain of transmission
    ip.calibration = false 
    ip.modeltime = 500    
    reset_params(ip)
    for x in propertynames(cv.p)
        @test getfield(ip, x) == getfield(mod_ip, x)
    end
    reset_params_default() ## check if it goes back to default
    ip = ModelParameters()
    for x in propertynames(cv.p)
        @test getfield(ip, x) == getfield(mod_ip, x)
    end        
end

@testset "demographics" begin
    reset_params_default()
    initialize()

    ## age groups
    inmodel = (1, 2, 3, 4, 5) ## possible age groups in model
    ags = []  
    for x in cv.humans
        push!(ags, x.ag)
        @test x.ag in inmodel
        @test x.exp == 999
        @test x.health == cv.SUS
        @test x.swap == cv.UNDEF
        @test x.iso == false
    end
    @test length(unique(ags)) == length(inmodel) # check if all age groups are sampled
    
    ## insert_infected
    for ag in inmodel 
        initialize() # reset population
        insert_infected(1, ag) # 1 infected in age group ag
        @test length(findall(x -> x.health == cv.INF && x.ag == ag, cv.humans)) == 1
    end
    cv.p.fsevere = 1.0 
    initialize() # reset population
    insert_infected(1, 1) # 1 infected in age group ag
    @test length(findall(x -> x.health == cv.INF && x.swap == cv.IISO, cv.humans)) == 1
end
