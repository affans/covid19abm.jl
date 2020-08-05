using covid19abm
const cv = covid19abm
using Test

const infectiousstates = (cv.LAT, cv.MILD, cv.MISO, cv.INF, cv.IISO, cv.HOS, cv.ICU)


# @testset "parameters" begin 
#     ip = ModelParameters()
#     mod_ip = cv.p ## get the module parameters    
#     ## set random parameters 
#     ip.β = 99
#     ip.prov = :noknownlocation 
#     ip.τmild = 1 ## days before they self-isolate for mild cases
#     ip.fmild = 0.05  ## percent of people practice self-isolation
#     ip.fsevere = 0.80 # fixed at 0.80, always within 1 day. 
#     ip.eldq = 0.0 ## percentage quarantine of 60+ individuals, removes from chain of transmission
#     ip.calibration = false 
#     ip.modeltime = 500    
#     ip.popsize = 20000
#     cv.reset_params(ip)
#     for x in propertynames(cv.p)
#         @test getfield(ip, x) == getfield(mod_ip, x)
#     end
# end

# @testset "init" begin
#     cv.reset_params()
#     cv.initialize()

#     # test if human array is resized to the correct length
#     @test length(cv.humans) == cv.p.popsize

#     ## age groups
#     inmodel = (1, 2, 3, 4, 5) ## possible age groups in model
#     ags = []  
#     for x in cv.humans
#         push!(ags, x.ag)
#     end
#     @test length(unique(ags)) == length(inmodel) # check if all age groups are sampled

#     ## initial human states (check only subset of humans)
#     for i in rand(1:100, 100)
#         x = cv.humans[i]
#         @test x.health == cv.SUS
#         @test x.swap == cv.UNDEF
#         @test x.sickfrom == cv.UNDEF
#         @test x.nextday_meetcnt >= 0 
#         @test x.tis == 0 
#         @test x.exp == 999
        
#         @test x.iso == false
#         @test x.isovia == :null

#         @test minimum(x.dur) != 0 #check if durations are properly set 
#         @test x.doi == 999 
#         @test x.tracing == false
#         @test x.tracestart == -1 
#         @test x.traceend == -1 
#         @test x.tracedxp == 0 
#     end
   
#     ## insert_infected check if infected person is added in the correct age group
#     for ag in inmodel 
#         cv.initialize() # reset population
#         cv.insert_infected(cv.INF, 1, ag) # 1 infected in age group ag
#         @test length(findall(x -> x.health == cv.INF && x.ag == ag, cv.humans)) == 1
#         cv.insert_infected(cv.REC, 1, ag) # 1 infected in age group ag
#         #@test length(findall(x -> x.health == cv.REC, cv.humans)) == 1
#     end

#     ## check if the initial infected person is NOT IISO 
#     cv.p.fsevere = 1.0 
#     cv.initialize() # reset population
#     cv.insert_infected(cv.INF, 1, 1) # 1 infected in age group ag
#     @test length(findall(x -> x.health == cv.INF && x.swap == cv.REC, cv.humans)) == 1

#     ## check if beta values are time dependent or NOT
#     cv.p.seasonal = true 
#     cv.p.β = 1.0
#     svalues = cv.td_seasonality()
#     cv.init_betas()
#     @test sum(svalues) == sum(cv.BETAS)

#     cv.p.seasonal = false 
#     cv.p.β = 1.0
#     cv.init_betas()
#     @test cv.p.modeltime == sum(cv.BETAS)
# end

# @testset "transitions" begin
#     cv.reset_params()
#     cv.initialize()
    
#     ## check if time in state is up by one
#     cv.time_update()
#     for i in rand(1:100, 100)
#         x = cv.humans[i]
#         @test x.tis == 1  
#         @test x.exp == 999
#         @test x.health == cv.SUS
#         @test x.swap == cv.UNDEF ## shouldn't set a swap 
#     end

#     for x in cv.humans 
#         x.exp = 0 ## this will trigger a swap, but no swap is set
#     end
#     @test_throws ErrorException("swap expired, but no swap set.") cv.time_update()

#     #check if it goes through all the move compartments and see if health/swap changes
    
#     for h in 2:9  ## latent to ded, ignore susceptible
#         cv.initialize()
#         rh = cv.HEALTH(h)
#         randhumans = rand(1:100, 100)
#         for i in randhumans 
#             x = cv.humans[i]
#             x.swap = rh
#             x.exp = 0 ## to force the swap
#         end
#         cv.time_update() ## since tis > exp 
#         for i in randhumans 
#             x = cv.humans[i]
#             @test x.health == rh
#             if rh ∈ infectiousstates  ## for all states, except rec/ded there should be a swap
#                 @test x.swap != cv.UNDEF
#             end
#         end    
#     end

#     # check latent function
#     cv.initialize()
#     x = cv.humans[1]
#     x.ag = 1 ## move to the first age group manually.
#     myp = cv.ModelParameters()   
#     cv.reset_params(myp)
#     # check basic variables    
#     cv.move_to_latent(x)
#     @test x.swap ∈ (cv.ASYMP, cv.PRE)
#     @test x.iso == false
#     @test x.doi == 0 
#     @test x.tis == 0 
#     @test x.exp == x.dur[1]

#     # check split between asymp/pre
#     # split to asymptomatic is now fixed in the model... not a parameter, see old commit for test code

#     # CHECKING PRE
#     cv.initialize()
#     x = cv.humans[1]
#     cv.reset_params(myp)
#     cv.move_to_pre(x) 
#     @test x.health == cv.PRE 
#     @test x.swap ∈ (cv.MILD, cv.INF)
   
#     ## check if individual moves through mild, miso through fmild, tmild parameters
#     cv.initialize()
#     x = cv.humans[1]
#     x.ag = 1 ## move to the first age group manually.
#     x.iso = false ## turn this off so we can test the effect of fmild, tmild
#     myp.fmild = 0.0
#     cv.reset_params(myp)
#     cv.move_to_mild(x)
#     @test x.swap == cv.REC

#     myp.fmild = 1.0
#     myp.τmild = 1
#     cv.reset_params(myp)
#     cv.move_to_mild(x)
#     @test x.swap == cv.MISO
#     @test x.exp == 1
    
#     cv.time_update() 
#     @test x.health == cv.MISO
#     @test x.swap == cv.REC
    
#     # check if x.iso and x.isovia property are set correctly. 
#     # 1) check if they are isolated through quarantine
#     myp = cv.ModelParameters()   
#     myp.eldq = 1.0  ## turn on eldq
#     myp.eldqag = 4
#     cv.reset_params(myp)
#     cv.initialize()    
#     for x in cv.humans
#         if x.ag  == myp.eldqag 
#             @test x.iso == true 
#             @test x.isovia == :qu
#         else 
#             @test x.iso == false 
#             @test x.isovia == :null
#         end
#     end
#     # 2) check if they are isolated through presymptomatic capture 
#     # no need to test for mild/inf movement since it's a simple assignment in code
#     # more important to check through `main` when all the function dynamics are happening. 
#     myp = cv.ModelParameters()   
#     myp.eldq = 0.0  ## turn off eldq
#     myp.fsevere = 0.0 ## turn on fsevere
#     myp.fmild = 0.0
#     myp.fpreiso = 1.0 
#     cv.reset_params(myp)
#     cv.initialize()  
#     for x in cv.humans 
#         cv.move_to_pre(x)
#         @test x.iso == true && x.isovia == :pi
#     end
#     # testing contact tracing isolation below

#     cv.initialize()
#     x = cv.humans[1]
#     x.iso = true
#     for i = 1:100 # stress test the function 
#         @test cv.get_nextday_counts(x) <= 3
#     end

#     x.health = cv.DED
#     for i = 1:100 # stress test the function 
#         @test cv.get_nextday_counts(x) == 0
#     end
# end

# @testset "tranmission" begin
#     cv.reset_params()
#     cv.initialize()
#     grps = cv.get_ag_dist()
    
#     # get an array of everyones meet counts 
#     ev_meet_cnts = [x.nextday_meetcnt for x in cv.humans]

#     # since beta = 0 default, and everyone sus so no inf to loop over
#     tm, ti = cv.dyntrans(1, grps)  
#     @test ti == 0
#     @test tm == 0

#     # insert some infectious
#     cv.insert_infected(cv.PRE, 1, 1) # 1 infected in age group ag
#     tm, ti = cv.dyntrans(1, grps) 
#     @test ti == 0 # still zero cuz beta = 0
#     #@test tm > 0 # count still actually be zero but very unlikely. 
#     updated_ev_meet_cnts = [x.nextday_meetcnt for x in cv.humans]
#     lc_by_one = [x.nextday_meetcnt for x in cv.humans]
#     @test tm == sum(ev_meet_cnts - lc_by_one) # test algo. ev_meet_cnts is the original counts. lc_by_one should be counts -1 (except that single infectious person could meet the same person twice)
#     # but now some contacts should be -1 


#     # now change beta 
#     cv.reset_params(cv.ModelParameters(β = 1.0)) # to account for relative transmission    
#     tm, ti = cv.dyntrans(1, grps)  
#     #@test ti > 0 ## actually may still be zero because of stochasticity but very unliekly 

#     cv.reset_params(cv.ModelParameters(β = 1.0, seasonal = false, frelasymp=0.11))   
#     @test cv._get_betavalue(1, cv.ASYMP) == 1.0*0.11

#     cv.reset_params(cv.ModelParameters(β = 1.0, seasonal = false, frelasymp=0.11))   
#     @test cv._get_betavalue(1, cv.MILD) == 1.0*0.44

#     cv.reset_params(cv.ModelParameters(β = 1.0, seasonal = false, frelasymp=0.11))   
#     @test cv._get_betavalue(1, cv.INF) == 1.0*0.89
# end

# @testset "calibration" begin
#     # to do, run the model and test total number of infections 
#     myp = cv.ModelParameters()
#     myp.β = 1.0
#     myp.prov = :ontario
#     myp.calibration = true
#     myp.fmild = 0.0 
#     myp.fsevere = 0.0
#     myp.fpreiso = 0.0
#     myp.initialinf = 1
#     cv.reset_params(myp)
#     cv.initialize()
#     grps = cv.get_ag_dist()
#     cv.insert_infected(cv.PRE, 1, 4)
#     # find the single insert_presymptomatic person
#     h = findall(x -> x.health == cv.PRE, cv.humans)
#     x = cv.humans[h[1]]
#     @test length(h) == 1
#     # @test x.ag == 4  #not true anymore
#     @test x.swap ∈ (cv.ASYMP, cv.MILD, cv.INF) ## always true for calibration 
#     for i = 1:20 ## run for 20 days 
#         cv.dyntrans(i, grps)
#         cv.time_update()
#     end
#     @test x.health == cv.REC  ## make sure the initial guy recovered
#     @test x.exp == 999
#     ## everyone should really be in latent (or susceptible) except the recovered guy
#     all = findall(x -> x.health ∈ (cv.PRE, cv.ASYMP, cv.MILD, cv.MISO, cv.INF, cv.IISO, cv.HOS, cv.ICU, cv.DED), cv.humans)
#     @test length(all) == 0
#     ## to do, make sure everyone stays latent
    
# end

@testset "contact trace" begin
    myp = cv.ModelParameters()          
    myp.ctstrat = 1  # turn on contact tracing
    myp.fctcapture = 1.0 # force everyone to get captured.
    myp.fcontactst = 1.0
    myp.cdaysback = 3 # go back 3 days from time of identification
    myp.mild_props = (0, 0, 0, 0, 0)
    cv.reset_params(myp)
    cv.initialize()
    x = cv.humans[1]

    # for dyntrans function 
    grps = cv.get_ag_dist()

    # test default values. important. 
    @test x.tracestart == -1 
    @test x.traceend == -1

    # test whether tracestart and end are set correctly
    cv.move_to_latent(x)
    @test x.doi == 0 
    @test x.tis == 0 
    @test x.exp == x.dur[1]
    
    x.swap = cv.ASYMP 
    @test x.swap == cv.ASYMP 
    cv.ct_dynamics(x) 
    # test default values. important. 
    @test x.tracestart == -1 
    @test x.traceend == -1

    x.swap = cv.PRE # force to PRE instead of ASYMP (ASYMP can work too. )
    @test x.swap == cv.PRE 
    cv.ct_dynamics(x) # since ct_dynamics will run when move_to_latent will run in time update func
    @test x.tracestart > 0
    @test x.traceend > 0 
    @test x.traceend > x.dur[1] + x.dur[3] ## bigger than since day to identification is randomly sampled
    @test x.traceend <=  x.dur[1] + x.dur[3] + x.dur[4]
    
    # save the tracestart and end in local variables. 
    trace_srt = x.tracestart
    trace_end = x.traceend
    println("t start: $trace_srt t end: $trace_end")

    # go through entire inf period of the person
    # first, lets make sure the person goes to 
    totalinfperiod = x.dur[1] + x.dur[3] + x.dur[4] # we can do this because we forced x to go to PRE
    @test totalinfperiod > 0 
    println("lat, pre, symp durations: $(x.dur[1]) + $(x.dur[3]) + $(x.dur[4])")
   
    totaltracedays = 0 # how how many days are true, should match cdaysback
    for i = 1:(trace_end + 2) # plus an 2 days to check if tracing turns off
        #println("doi: $0")
        println("tis: $(x.tis), doi: $(x.doi), health: $(x.health),")
        cv.dyntrans(1, grps)
        cv.time_update()
        println("...dyntrans(); timeupdate() => health: $(x.health), tis: $(x.tis), doi: $(x.doi)")
        @test x.doi == i # should increase with the counter.
        
        # make sure the tracestart/end only are set once when LAT and doi==0
        if x.health == cv.LAT
            @test x.tracestart == trace_srt
            @test x.traceend == trace_end
            @test sum(x.met_contacts) == 0
        end

        if x.health == cv.ASYMP || x.health == cv.MILD
            error("should never happen since forcing x to PRE/INF")
        end

        traced_cts = findall(x -> x == true, x.met_contacts)
        println("...total traced: $(length(traced_cts))")
        # when doi == tracestart, check if x.tracing is set true 
        # and remains true until traceend
        if x.doi >= x.tracestart && x.doi < x.traceend
            @test x.tracing == true 
            totaltracedays += 1  
            for i in traced_cts
                y = humans[i] 
                @test y.tracedxp == 0  # tracedxp is set to 14 on day of identification
            end
        end

        # day of identification, check if trace is off. 
        if x.doi == x.traceend
            @test x.tracing == false          
            for i in traced_cts
                y = humans[i] 
                @test y.tracedxp == 13  # tracedxp is set to 14 on day of identification (one day is already passed in time_update so we check for 13)
            end 
        end        
    end
    @test totaltracedays == cv.p.cdaysback

    traced_cts = findall(x -> x == true, x.met_contacts)
    for i in traced_cts
        y = humans[i] 
        @test y.iso == true 
        @test y.isovia == :ct
    end

    ## have to check whether the trace goes away. 
    ## reset the tracedxp to 14 days (since it's at 11 based on the above loop)
    for i in traced_cts
        y = humans[i] 
        y.tracedxp = 14
    end

    for t = 1:13
        cv.time_update()
        for i in traced_cts
            y = humans[i] 
            @test y.tracedxp == 14 - t     
            @test y.iso == true 
            @test y.isovia == :null       
        end    
    end
    # do t = 14 manually 
    cv.time_update()
    for i in traced_cts
        y = humans[i] 
        @test y.tracedxp == 0
        @test y.iso == false
        @test y.isovia == :null        
    end     

    ## reset everything ... check if cdaysback is a larger number trace start should be zero
end

@testset "main run" begin 
    ## run model with high beta 
    myp = cv.ModelParameters()
    myp.β = 0.0525 
    myp.prov = :newyork    
    ## run empty scenario
    
    myp.τmild = 0
    myp.fmild = 0.0
    myp.fsevere = 0.0
    myp.eldq = 0.0  
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
    
#     function string_or_char(s::SubString)
#         length(s) == 1 && return first(s)
#         return String(s)
#     end
#     arr = Union{Char,String}[]
# for i in split_string
#     push!(arr, string_or_char(i))
# end

    # function addday(n)
    #     strtorep = "---" * "\u00B7"
    #     str = ""
    #     for i = 1:n 
    #         str =str * strtorep 
    #     end
    #     str = str[1:end-1]
    #     str = " " * str * " "
    #     return str
    # end
    

    # function printascii_epi(x) 
    #     fs = "" 
    #     fs = fs * "L" * addday(x.dur[1])
     
    #     fs = fs * "P" * addday(x.dur[3])
    #     fs = fs * "S" * addday(x.dur[4])
    #     return fs
    # end
   

end


