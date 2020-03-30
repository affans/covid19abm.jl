using covid19abm
const cv = covid19abm
using Test

const infectiousstates = (cv.LAT, cv.MILD, cv.MISO, cv.INF, cv.IISO, cv.HOS, cv.ICU)


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
    cv.reset_params_default()
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
        @test x.tis == 0
        @test x.exp == 999
    end
    @test length(unique(ags)) == length(inmodel) # check if all age groups are sampled
    
    ## insert_infected check if infected person is added in the correct age group
    for ag in inmodel 
        initialize() # reset population
        insert_infected(1, ag) # 1 infected in age group ag
        @test length(findall(x -> x.health == cv.INF && x.ag == ag, cv.humans)) == 1
    end
    ## check if the initial infected person is NOT IISO 
    cv.p.fsevere = 0.0 
    initialize() # reset population
    insert_infected(1, 1) # 1 infected in age group ag
    @test length(findall(x -> x.health == cv.INF && x.swap == cv.REC, cv.humans)) == 1
    ## check if the initial infected person is IISO 
    cv.p.fsevere = 1.0 
    initialize() # reset population
    insert_infected(1, 1) # 1 infected in age group ag
    @test length(findall(x -> x.health == cv.INF && x.swap == cv.IISO, cv.humans)) == 1
end

@testset "transitions" begin
    cv.reset_params_default()
    initialize()
    
    ## check if tis is up by one
    time_update()
    for x in humans 
        @test x.tis == 1  
        @test x.health == cv.SUS
        @test x.swap == cv.UNDEF ## shouldn't set a swap 
    end

    #check if it goes through all the move compartments 
    
    for h in 2:9  ## latent to ded, ignore susceptible
        initialize()
        rh = cv.HEALTH(h)
        for x in humans 
            x.swap = rh
            x.exp = 0 ## to force the swap
        end
        time_update() ## since tis > exp 
        for x in humans[1:5]
            @test x.health == rh
            if rh ∈ infectiousstates 
                @test x.swap != cv.UNDEF
            end
        end    
    end

    # to do: check each compartment separately.
    # to do: follow the movement of a single human
end


@testset "dyntrans" begin
    cv.reset_params_default()
    initialize()
    
    # since beta = 0 default, and everyone sus
    totalinf = dyntrans()  
    @test totalinf == 0

    # check with only a single infected person 
    insert_infected(1, 1) # 1 infected in age group ag
    totalinf = dyntrans()  
    @test totalinf == 0 # still zero cuz beta = 0

    # now change beta 
    cv.reset_params(ModelParameters(β = 1.0))
    totalinf = dyntrans()  
    @test totalinf > 0 ## actually may still be zero because of stochasticity but very unliekly 

    ## somehow check the transmission reduction and number of contacts
end

@testset "model" begin
    # to do, run the model and test total number of infections 
end
