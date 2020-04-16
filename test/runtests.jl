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
end

@testset "init" begin
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
        @test x.isovia == :null
        @test x.tis == 0
        @test x.exp == 999
        @test x.sickfrom == cv.UNDEF
    end
    @test length(unique(ags)) == length(inmodel) # check if all age groups are sampled
    
    ## insert_infected check if infected person is added in the correct age group
    for ag in inmodel 
        initialize() # reset population
        insert_infected(cv.INF, 1, ag) # 1 infected in age group ag
        @test length(findall(x -> x.health == cv.INF && x.ag == ag, cv.humans)) == 1
    end
    ## check if the initial infected person is NOT IISO 
    cv.p.fsevere = 0.0 
    initialize() # reset population
    insert_infected(cv.INF, 1, 1) # 1 infected in age group ag
    @test length(findall(x -> x.health == cv.INF && x.swap == cv.REC, cv.humans)) == 1
    ## check if the initial infected person is REC (since simpleinf function only puts them to REC)
    cv.p.fsevere = 1.0 
    initialize() # reset population
    insert_infected(cv.INF, 1, 1) # 1 infected in age group ag
    @test length(findall(x -> x.health == cv.INF && x.swap == cv.REC, cv.humans)) == 1
end

@testset "transitions" begin
    cv.reset_params_default()
    initialize()
    
    ## check if time in state is up by one
    cv.time_update()
    for x in cv.humans 
        @test x.tis == 1  
        @test x.exp == 999
        @test x.health == cv.SUS
        @test x.swap == cv.UNDEF ## shouldn't set a swap 
    end

    for x in cv.humans 
        x.exp = 0 ## this will trigger a swap
    end
    @test_throws ErrorException("swap expired, but no swap set.") cv.time_update()

    #check if it goes through all the move compartments and see if health/swap changes
    
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
            if rh ∈ infectiousstates  ## for all states, except rec/ded there should be a swap
                @test x.swap != cv.UNDEF
            end
        end    
    end

    
    # CHECKING LATENT
    initialize()
    x = humans[1]
    x.ag = 1 ## move to the first age group manually.
    myp = cv.ModelParameters()   
    myp.fpre = 0.0  ## turn off presymptomatic 
    cv.reset_params(myp)
    move_to_latent(x)
    @test x.swap ∈ (cv.ASYMP, cv.MILD, cv.INF)
    @test x.iso == false

    myp.fpre = 1.0  ## turn on presymptomatic 
    cv.reset_params(myp)
    move_to_latent(x)
    @test x.swap == cv.PRE
    @test x.iso == false # isolation should not be turned on in latent at all.

    myp.fpre = 0.0  ## turn off presymptomatic but turn on asymp 
    myp.fasymp = 1.0
    cv.reset_params(myp)
    move_to_latent(x)
    @test x.swap ∈ (cv.ASYMP, cv.INF)
    @test x.iso == false

    # CHECKING PRE
    initialize()
    x = humans[1]
    cv.reset_params(myp)
    move_to_pre(x) 
    @test x.health == cv.PRE 
    @test x.swap ∈ (cv.ASYMP, cv.MILD, cv.INF)
   
    myp.fpreiso = 0.0
    cv.reset_params(myp)
    move_to_pre(x) 
    @test x.health == cv.PRE 
    @test x.swap ∈ (cv.ASYMP, cv.MILD, cv.INF)
    
    ## check if individual moves through mild, miso through fmild, tmild parameters
    initialize()
    x = humans[1]
    x.ag = 1 ## move to the first age group manually.
    x.iso = false ## turn this off so we can test the effect of fmild, tmild
    myp.fmild = 0.0
    cv.reset_params(myp)
    move_to_mild(x)
    @test x.swap == cv.REC

    myp.fmild = 1.0
    myp.τmild = 1
    cv.reset_params(myp)
    move_to_mild(x)
    @test x.swap == cv.MISO
    @test x.exp == 1
    
    time_update() 
    @test x.health == cv.MISO
    @test x.swap == cv.REC
    
    ## todo: check if x.iso and x.isovia property are set correctly. 
    # 1) check if they are isolated through quarantine
    myp = cv.ModelParameters()   
    myp.eldq = 1.0  ## turn on eldq
    cv.reset_params(myp)
    cv.initialize()    
    for x in cv.humans
        if x.age >= 60 
            @test x.iso == true 
            @test x.isovia == :qu
        else 
            @test x.iso == false 
            @test x.isovia == :null
        end
    end
    # 2) check if they are isolated through presymptomatic capture 
    # no need to test for mild/inf movement since it's a simple assignment in code
    # more important to check through `main` when all the function dynamics are happening. 
    myp = cv.ModelParameters()   
    myp.eldq = 0.0  ## turn off eldq
    myp.fsevere = 0.0 ## turn on fsevere
    myp.fmild = 0.0
    myp.fpreiso = 1.0 
    cv.reset_params(myp)
    cv.initialize()  
    for x in cv.humans 
        cv.move_to_pre(x)
        @test x.iso == true && x.isovia == :pi
    end
end


@testset "tranmission" begin
    cv.reset_params_default()
    initialize()
    
    # since beta = 0 default, and everyone sus
    totalinf = dyntrans(1)  
    @test totalinf == 0

    # check with only a single infected person 
    insert_infected(cv.INF, 1, 1) # 1 infected in age group ag
    totalinf = dyntrans(1)  
    @test totalinf == 0 # still zero cuz beta = 0

    # now change beta 
    cv.reset_params(ModelParameters(β = 1.0))
    totalinf = dyntrans(1)  
    @test totalinf > 0 ## actually may still be zero because of stochasticity but very unliekly 

    ## somehow check the transmission reduction and number of contacts
end

@testset "calibration" begin
    # to do, run the model and test total number of infections 
    myp = ModelParameters()
    myp.β = 1.0
    myp.prov = :ontario
    myp.calibration = true
    myp.fmild = 0.0 
    myp.fsevere = 0.0
    myp.fpreiso = 0.0
    myp.fasymp = 0.5
    myp.initialinf = 1
    cv.reset_params(myp)
    cv.initialize()
    cv.insert_infected(cv.PRE, 1, 4)
    # find the single insert_presymptomatic person
    h = findall(x -> x.health == cv.PRE, cv.humans)
    x = humans[h[1]]
    @test length(h) == 1
    @test x.ag == 4
    @test x.swap ∈ (cv.ASYMP, cv.MILD, cv.INF) ## always true for calibration 
    for i = 1:20 ## run for 20 days 
        cv.dyntrans(i)
        cv.time_update()
    end
    @test x.health == cv.REC  ## make sure the initial guy recovered
    @test x.exp == 999
    ## everyone should really be in latent (or susceptible) except the recovered guy
    all = findall(x -> x.health ∈ (cv.PRE, cv.ASYMP, cv.MILD, cv.MISO, cv.INF, cv.IISO, cv.HOS, cv.ICU, cv.DED), cv.humans)
    @test length(all) == 0
    ## to do, make sure everyone stays latent
    
end

@testset "contact trace" begin
    myp = ModelParameters()
    myp.β = 1.0 
    cv.reset_params(myp)
    cv.initialize()
    hdx = rand(1:10000, 100) ## sample 100 humans instead of all 10000 
    for i in hdx 
        move_to_pre(humans[i])
        @test humans[i].tracing == true
    end

    ## test the contact_tracing() function 
    cv.initialize()    
    tracer = cv.humans[1]
    move_to_pre(tracer) ## move random human to presymptomatic
    @test tracer.tis == 0 ## these tests are not really needed, but good to verify again
    @test tracer.exp == 1
    @test tracer.tracing == true
    cv.dyntrans(1) ## go through a tranmission cycle
    alltraced = findall(x -> x.tracedby > 0, cv.humans) 
    @test length(alltraced) == 0## since we havn't turned fctcapture > 0
    myp.fctcapture = 1.0 
    cv.reset_params(myp)
    cv.dyntrans(1) ## go through a tranmission cycle
    alltraced = findall(x -> x.tracedby > 0, cv.humans) ## since we havn't turned fctcapture > 0
    @test length(alltraced) > 0
    cv.contact_tracing(1)
    for i in alltraced
        y = cv.humans[i]
        @test y.tracedby == 1 ## since we used the first human as the tracing contact
        @test y.tracedxp == 1 + 14 
        @test y.iso == false ## the first human is not in INF stage yet.. still in presymp
        @test y.isovia == :null 
    end
    cv.move_to_inf(tracer) ## can use time_update to move because may go to mild/asymp, tis = 0
    cv.time_update() #will move, tis = 1
    @test tracer.health == cv.INF && tracer.tis == 1  ## this condition is used inside contact_tracing
    cv.contact_tracing(3)  ## technically time is now "3" after three time_updates
    for i in alltraced
        y = cv.humans[i]
        @test y.tracedby == 1 ## since we used the first human as the tracing contact
        @test y.tracedxp == 1 + 14 
        @test y.iso == true ## the first human is not in INF stage yet.. still in presymp
        @test y.isovia == :ct 
    end
end

@testset "main run" begin 
    ## run model with high beta 
    myp = cv.ModelParameters()
    myp.β = 0.0525 
    myp.prov = :newyork    
    ## run empty scenario
    ## this wont return calibrated scenario since fasymp = 0
    myp.τmild = 0
    myp.fmild = 0.0
    myp.fsevere = 0.0
    myp.eldq = 0.0  
    myp.fasymp = 0.5
    myp.fpre = 1.0
    myp.fpreiso = 0.0 
    myp.tpreiso = 0
    myp.fctcapture = 0.0
    cv.runsim(1, myp) # warm up the functions
    println("time with contact tracing off:")
    @time results = cv.main(myp)

    myp.fctcapture = 1.0
    cv.runsim(1, myp)
    println("time with contact tracing on:")
    @time results = cv.main(myp)
    

end


