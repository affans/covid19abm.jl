# COVID-19 Agent Based Model
Affan Shoukat, 2020  
Center for Infectious Disease Modelling and Analysis  
Yale University, New Haven, USA

### Model details: 
A stochastic, age-stratified agent-based computational model for the transmission dynamics of COVID-19. The computational model simulates autonomous agents (representing individuals in a human population) and their interactions within a constrained virtual environment. Agents follow the natural history of disease, including epidemiological stages of susceptible, infected and incubating, asymptomatic, presymptomatic, and symptomatic with either mild, severe, or critical illness, recovered, and dead. 

Model features include: 
 - Age structured with realistic contact dynamics
 - Asymptomatic, Presymptomatic transmission
 - Isolation of mild/severe cases 
 - Isolation of presymptomatic individuals 
 - Highly flexible contact tracing dynamics
 - Shelter-in-place mechanics for specific 
 - Calibration mode
 - Unit Tested with over 100,000 tests. 

### How to download and run
The model is a self-contained (non-registered) Julia package. To use the model, add the package by typing 

    ] add https://github.com/affans/covid19abm.jl 
    
in the Julia REPL. 

The model is run by 1) instantiating a `ModelParameters` object with the desired parameter values and 2) running the `main` function. See the snippet below. 

    julia> using covid19abm

    julia> myp = ModelParameters()
    ModelParameters
        β: Float64 0.0             # transmission value
        prov: Symbol ontario       # demographics (see available demographics in get_province_ag() function)
        calibration: Bool false    # turn calibration mode on/off
        modeltime: Int64 500       # model run time (in days)
        initialinf: Int64 1        # initial infected to seed the model
        τmild: Int64 0             # time to isolation for mild individuals
        fmild: Float64 0.0         # fraction of mild individuals to be isolated
        fsevere: Float64 0.0       # fraction of severe individuals to be isolated
        eldq: Float64 0.0          # fraction of population isolated at start of simulation
        eldqag: Int8 5             # the age group isolated at start of simulation (in combination with eldq)
        fasymp: Float64 0.5        # proportion of infected individuals that will be asymptomatic
        fpre: Float64 1.0          # (not used)
        fpreiso: Float64 0.0       # proportion of presymptomatics that are isolated 
        tpreiso: Int64 0           # time at which presymptomatic isolation is turned on
        frelasymp: Float64 0.11    # the relative transmission of asymptomatic individuals
        ctstrat: Int8 0            # contact tracing strategy 1 or 2 
        fctcapture: Float16 Float16(0.0)  # fraction of symptomatic individuals identified for contact tracing
        fcontactst: Float16 Float16(0.0)  # fraction of the contacts that are traced
        cidtime: Int8 0            # how many days post symptom onset is contact tracing turned on
        cdaysback: Int8 0          # how many days to trace back

    julia> main(myp)
    10000×500 Array{Int64,2}:

The function`main` returns a `10000 x modeltime` matrix where each element is one of the following model states: `(SUS LAT MILD MISO INF IISO HOS ICU REC DED)` coded from `1 ... 9`.  This matrix can be used to calculate the incidence and prevalence of the different model states. See functions `_get_column_incidence` and `_get_column_prevalence`. 

To evaluate different scenarios set the desired model parameters. Here are the description of the model parameters.

Since the model is stochastic, many realizations are required to get an accurate picture of the results. We recommend running this in an embarrassingly parallel manner. This essentially means running the `main` function repeatedly, saving the results for each replicate. This can be done very easily using Julia's Distributed library. Simply `addprocs` or using `ClusterManagers` to connect to a cluster to launch `n` number of parallel workers. Then use `pmap` to run the main function on each worker, saving the results to be processed later. 

See the function `run()` in `scripts/local_run.jl` for an example of using the parallel programming library.  

    
### Published Studies and Reproducibility 
In order to reproduce the results for published studies, `clone` the repository rather than adding it as a package to Julia. Once the repository is cloned, the version of the code that was used to produced the results can be pulled by the relevant tag. See below

> Projecting Critical Care Demand for COVID-19 Outbreaks in Canada. Affan Shoukat, Chad R. Wells, Joanne M. Langley, Burton H. Singer, Alison P. Galvani, Seyed M. Moghadas, 2020. Canadian Medical Association Journal (in review).

`git checkout v2CAN` 

> The implications of silent transmission for the control of COVID-19 outbreaks. Seyed M. Moghadas, Meagan C. Fitzpatrick, Pratha Sah, Abhishek Pandey, Affan Shoukat, Burton H. Singer, Alison P. Galvani. 2020. Proceedings of the National Academy of Science 
 
`git checkout v1PNAS_RoleOfPre` 


### Contribute and Development
PRs are welcome. Message me to get an explanation of how the model works. If submitting PR, please add relevant tests and make sure all tests pass. 

```
(covid19abm) pkg> test
   Testing covid19abm
 Resolving package versions...
Test Summary: | Pass  Total
parameters    |   20     20
Test Summary: |  Pass  Total
init          | 90108  90108
Test Summary: |  Pass  Total
transitions   | 70085  70085
Test Summary: | Pass  Total
tranmission   |    3      3
Test Summary: | Pass  Total
calibration   |    6      6
Test Summary: | Pass  Total
contact trace |   16     16
time with contact tracing off:
  2.134143 seconds (26.56 M allocations: 652.800 MiB, 8.63% gc time)
time with contact tracing on:
  2.259248 seconds (26.41 M allocations: 633.013 MiB, 12.05% gc time)
Test Summary: |
main run      | No tests
   Testing covid19abm tests passed
``` 
The model runtime should be less than 5 seconds. 