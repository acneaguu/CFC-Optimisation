%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%README:
%%This is the main script for the optimisation unit. In this script, the
%%required variables are initialised, desired functions are called and the
%%optimisation is started and ended. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

%%Start timer
total_execution_time = tic;

%%Load MATPOWER constants for convenience in the struct CONSTANTS
define_constants_struct();

%%For reproducibility (needed for PS algorithm)
rng default  

%%-------------------------------------------------------------------------
%%%Optimisation problem specification and settings 
%%Optimisation containts the optimisation problem parameters
global Optimisation ff_par Systemdata;
%%Description of variables to optimise
Optimisation.Nturbines = 13;                    %Number of turbines:
                                                %-13 for stringlevel, 
                                                %-91 for turbinelevel
Optimisation.Npv = 4;                           %Number of pv generator strings
Optimisation.Ntr = 2;                           %Number of transformers with discrete tap positions
Optimisation.Ntaps = [17;17];                   %Number of tap positions per transformer                                             
                                                %(must have dimension of Ntr and separate by ;)
Optimisation.Nr = 1;                            %Number of discrete reactors

Optimisation.Nvars = Optimisation.Nturbines + Optimisation.Npv + ...
    Optimisation.Ntr + Optimisation.Nr;         %Number of optimisation variables
Optimisation.which_discrete = [18:20];          %Indeces of the discrete variables
% Optimisation.steps =[0.0168235 0.0168235 1];  %Steps of the discrete variables
logic_optvars();                                %Logic vectors for optimisation vector
if Optimisation.Nturbines == 13
    initialise_systemdata(system_13_350MVA);    %Initialise strinlevel topology
elseif Optimisation.Nturbines == 91
    initialise_systemdata(system_91_100MVA);    %Initialise turbinelevel topology
end

%%Optimisation run settings
initialise_optimisation_weights();  %Sets the weights of the different 
                                    %constraints and objectives
Optimisation.Ncases = 25;            %Number of evaluated time instances
Optimisation.Nruns = 5;             %Number of runs per case
Optimisation.Neval = 500*35;        %Max allowed function evaluations
Optimisation.Populationsize = 35;   %Size of the population
Optimisation.algorithm = 4;         %1 for ga, 2 for pso, 3 for cdeepso 
                                    %4 for MVMO_SHM
Optimisation.print_progress = 1;    %Plots runs in command window
Optimisation.print_interval = 2000; %Interval of printed steps
Optimisation.print_pfresults = 1;   %Plots powerflow results of optimal solution

%%-------------------------------------------------------------------------
%%Settings to plot the power flow and store the results of the optimisation
plot = 0;
store_results = 0;

%%-------------------------------------------------------------------------
%Results struct consits of the results of each optimal powerflow
%%variables containing the best solutions at all evaluated optimisation
%%runs (Fbest), a matrix containing the best solution at each optimisation
%%run (Xbest), the progress of the best fitness value of each run
%%(Fit_progress), the accuracy of Qpcc (Qaccuracy) and the values of the OF
%%paramers at each run (Ploss, tchanges and rchanges)
global Results;
Optimisation.t = 1;
initialise_results_struct();

%%-------------------------------------------------------------------------
%%Fitness evaluation function
switch Optimisation.algorithm
    case {1,2,3}
    fun = @(X)fitness_eval(X);
    case 4
    fun = str2func('fitness_eval');
end

%%Parameters for the different algorithms
switch Optimisation.algorithm
    case 1
    options = optimoptions('ga', 'FunctionTolerance', 1e-9, ...
    'MaxGenerations',Optimisation.Neval/Optimisation.Populationsize,...
    'PopulationSize',Optimisation.Populationsize);
    case 2
    options=optimoptions('particleswarm','MaxIterations',...
        Optimisation.Neval,'SwarmSize',Optimisation.Populationsize);
    case 3
    initialise_cdeepso(); %This function initialises the CDEEPSO settings
    
    case 4 
    initialise_mvmoshm(); %This function initialises the MVMO-SHM settings
end
%%-------------------------------------------------------------------------
%% run optimisation
global Keeptrack FCount;    %Some global vars to keep track of the calls of 
                            %the fitness evaluation funtion 
%%-------------------------------------------------------------------------
%%Setpoint at PCC given by TSO
global Qref;    
Qref.setpoint =  [-0.286; -0.143; 0; 0.143; 0.286]; %in p.u. of baseMVA
Qref.tolerance = 6.25/Systemdata.mpc.baseMVA;

%%-------------------------------------------------------------------------
%%Define the testcase
 v = [4.5 4.5 4.5 4.5 4.5 5 5 5 5 5 7 7 7 7 7 12 12 12 12 12 15 15 15 15 15]';

cases(:,1) = v;
cases(:,2) =repmat(Qref.setpoint,5,1);

if Optimisation.Npv > 0
    irradiance = [50 50 50 50 50 340 340 340 340 340 680 680 680 680 680 ...
        510 510 510 510 510 170 170 170 170 170];
    cases(:,3) = irradiance;
end
%%-------------------------------------------------------------------------
%%Run different cases
    for j = 2:Optimisation.Ncases+1
        %%-------------Prepare for optimisation:---------------------
        %%Set j for internal use
        Optimisation.t = j;
        
        %%Initialise the Results struct with NaNs for each case
        initialise_results_struct(); 
        
        %%Compute the allowed range of Qpcc w.r.t. the setpoints
        qpcc_limits(cases(j-1,2)); 
        
        %%Compute the reactive power generation per string depending on the
        %%windspeed
        if Optimisation.Npv > 0
            [Qmin_wtg, Qmax_wtg, Qmin_pvg, Qmax_pvg] = generate_case(cases(j-1,1),cases(j-1,3));
        else
            [Qmin_wtg, Qmax_wtg, Qmin_pvg, Qmax_pvg] = generate_case(cases(j-1,1));          
        end
        
        %%Update boundaries lb/ub
        [lb, ub]= boundary_initialise(Qmin_wtg, Qmax_wtg, Qmin_pvg, Qmax_pvg);
 
        %%-------------Start optimisation----------------------------------
        %%Case duration timer
        start_case = tic;
        
        %%Run a case multiple times
        for i = 1:Optimisation.Nruns
        tic;
        fprintf('************* Case %d, Run %d *************\n',j-1, i);

        %%Reinitialise fitness evaluation counter
        FCount = 0;

        %%Case of the different algorithms
        switch Optimisation.algorithm
            case 1
                X = ga(fun,Optimisation.Nvars,[],[],[],[],lb,ub,[],options);
            case 2
                X = particleswarm(fun,Optimisation.Nvars,lb,ub,options);
            case 3
                ff_par.fitEval = 0;
                ff_par.bestFitEval = 0;
                [gbestfit, X] = CDEEPSO_algorithm(fun,lb,ub);
            case 4
                [gbestfit, X] = mvmo_ceno(fun,lb,ub);
        end
        
        %%-------------Store the results into the results struct-----------
        %%Store the best solution and fitness of this run
        Results(j).Xbest(i+1,:) = round_discrete_vars(X);
        switch Optimisation.algorithm
            case {1,2}
                Results(j).Fbest(i+1) = Keeptrack.FitBest(end);
            case {3,4}
                Results(j).Fbest(i+1) = gbestfit;
        end

        %%Compute the Results(j) of the different OF parameters and Qpcc using the
        %%final solution and store them in results
        [Results(j).Ploss(i+1), Results(j).Tap_changes(i+1), Results(j).Reactors_changes(i+1)...
            ,Results(j).extremeness_setpoints(i+1), Results(j).total_cost_per_run(i+1), Results(j).Qaccuracy(i+1)]...
            = compute_results(Results(j).Xbest(i+1,:),cases(j-1,2));

        %%Initilise matrix with FitBest progress at each iteration
        if i == 1
            Results(j).Fit_progress = NaN * zeros(Optimisation.Nruns+1,FCount);
            Results(j).Violation_composition_progress = NaN * zeros(FCount,3,Optimisation.Nruns+1);
        elseif FCount > size(Results(j).Fit_progress,2)
            dis = FCount-size(Results(j).Fit_progress,2);
            Results(j).Fit_progress(:,end+1:FCount) = repmat(Results(j).Fit_progress(:,end),1,dis);
            Results(j).Violation_composition_progress(end+1:FCount,:,:) = ...
                 repmat(Results(j).Violation_composition_progress(end,:,:),dis,1,1);
        end
        
        %%-------------Results struct used for performance evaluation------
        %%Store the progress of FitBest of this iteration
        Results(j).Fit_progress(i+1,:) = Keeptrack.FitBest;
        Results(j).Violation_composition_progress(:,:,i+1) = Keeptrack.violation_composition;
        Results(j).runtime(i,1) = toc;
        
        %%Print the runtime of a run
        fprintf('Case %2d, Run %2d: %2f seconds \n',j-1,i,Results(j).runtime(i,1));
        
        %%Plot if desired
        if plot == 1
            animated_plot_fitness(Keeptrack.SolBest,Keeptrack.FitBest);
        end

        end
        
        %%Calculate best/worst/mean of Ploss
        MaxPloss = Systemdata.mpc.baseMVA;
        Results(j).Ploss_best = min(Results(j).Ploss);
        Results(j).Ploss_worst = max(Results(j).Ploss(Results(j).Ploss < MaxPloss));
        Results(j).Ploss_mean = mean(Results(j).Ploss(Results(j).Ploss < MaxPloss));
        
        %%Save the best fitness and solution 
        Results(j).Times_converged = sum(Results(j).Fbest<=1e3);
        best_index = find(Results(j).Fbest == min(Results(j).Fbest),1);
        Results(j).best_run_fitness = min(Results(j).Fbest);
        Results(j).best_run_solution = Results(j).Xbest(best_index,:);
        
        %%Calculates consistency performance
        Results(j).avg_fitness = mean(Results(j).Fbest(2:end));
        Results(j).std_fitness = std(Results(j).Fbest(2:end));
        Results(j).std_solution = std(Results(j).Xbest(2:end,:));
        
        %%Calculate cost per case
        Results(j).total_cost_per_case = mean(Results(j).total_cost_per_run(2:end));
        
        %%Compute the average runtime
        Results(j).avg_runtime = mean(Results(j).runtime(:,1));
        Results(j).total_runtime = toc(start_case);
        
        %%Print the total case runtime
        fprintf('Case %2d, Total Runtime: %2f seconds \n',j-1,Results(j).total_runtime);
    end

%Total costs of optimisation
total_cost = 0;
for j = 2:Optimisation.Ncases+1
    total_cost = total_cost + Results(j).total_cost_per_case;
end

%%-------------------------------------------------------------------------
%%Save the result if desired
if store_results == 1
    savedata
end

%%Save and print the total execution time
total_execution_time = toc(total_execution_time);
fprintf('Total Execution time: %2f seconds \n',total_execution_time);
