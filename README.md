# COVID-19 Agent Based Model
Affan Shoukat and colleagues, 2020  
Center for Infectious Disease Modelling and Analysis  
Yale University, New Haven, USA

### Model details: 
A stochastic, age-stratified agent-based computational model for the transmission dynamics of COVID-19. The computational model simulates autonomous agents (representing individuals in a human population) and their interactions within a constrained virtual environment. We accounted for the various epidemiological statuses of individuals, including susceptible, infected and incubating, infectious and symptomatic with either mild, severe, or critical illness, recovered, and dead. 

The latest model iteration also includes dynamics for presymptomatic and asymptomatic transmission. 

### How to download and run
The model is a self-contained (non-registered) Julia package. To use the model, add the package by typing 

    ] add https://github.com/affans/covid19abm.jl#v2CAN` 
    
in the Julia REPL. If trying reproduce CMAJ results, please see instructions below. 

The model is run by 1) instantiating a `ModelParameters` with the desired parameter values and 2) running the `main` function. See the snippet below.

    julia> using covid19abm

    julia> myp = ModelParameters()
    ModelParameters
    β: Float64 0.0
    prov: Symbol ontario
    τmild: Int64 1
    fmild: Float64 0.05
    fsevere: Float64 0.8
    eldq: Float64 0.0
    calibration: Bool false
    modeltime: Int64 500

    julia> hmatrix, ags = main(myp)
    new func
    ([0 0 … 0 0; 0 0 … 0 0; … ; 0 0 … 0 0; 0 0 … 0 0], [3, 4, 4, 4, 5, 3, 5, 2, 2, 5  …  2, 5, 4, 3, 5, 2, 4, 4, 4, 3])

The function`main` returns two objects: A matrix of size `(10,000, modeltime)` where each element is one of the following model states: `(SUS LAT MILD MISO INF IISO HOS ICU REC DED)` coded from `1 ... 9`.  This matrix can be used to calculate the incidence and prevalence of the different model states. The second object is how the model distributed the age groups over `10,000` individuals. 

To evaluate different scenarios set the desired model parameters. Here are the description of the model parameters.

    β: Float64 0.0           ## transmission probability
    prov: Symbol ontario     ## province or geographic location 
    τmild: Int64 1           ## time to self-isolate for mild symptomatic
    fmild: Float64 0.05      ## fraction of mild symptomatic who self-isolate
    fsevere: Float64 0.8     ## fraction of severe symptomatic who self-isolate
    eldq: Float64 0.0        ## fraction of individuals (60+ years old) completely quarantined at start of simulation
    calibration: Bool false  ## turn calibration mode on or off
    modeltime: Int64 500     ## the simulation run time in days

Since the model is stochastic, many realizations are required to get an accurate picture of the results. We recommend running this in an embarrassingly parallel manner. This essentially means running the `main` function repeatedly, saving the results for each replicate. This can be done very easily using Julia's Distributed library. Simply `addprocs` or using `ClusterManagers` to connect to a cluster to launch `n` number of parallel workers. Then use `pmap` to run the main function on each worker, saving the results to be processed later. See the function `run()` in `scripts/local_run.jl` for an example of using `addprocs`. 

#### CMAJ results reproducibility 
The model was used to project hospital and ICU capacity in Canada. The results are published in the journal: Canadian Medical Association Journal (CMAJ). For scientific reproducibility of the results in the manuscript, please add the following version of the code (and follow similar steps as above to run the model)

    ] add https://github.com/affans/covid19abm.jl#v2CAN 
    
### Published Studies
> Projecting Critical Care Demand for COVID-19 Outbreaks in Canada, Affan Shoukat, Chad R. Wells, Joanne M. Langley, Burton H. Singer, Alison P. Galvani, Seyed M. Moghadas, 2020. Canadian Medical Association Journal (in review).

### Contribute and Development
PRs are welcome. Message me to get an explanation of how the model works. 