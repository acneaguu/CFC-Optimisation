%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%README:
%%This function creates the Result struct and initialises the variables
%%with NaNs. This is done for the number of runs per case + 1. This
%%function should be called for each case.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initialise_results_struct
global Optimisation Results; 

%%Total runs is nruns + 1 (init). The initial run is used to specify the
%%initial topology status. 
total_runs = Optimisation.Nruns;  

%%Check whether it is the first case. If true, then initialise initial tap
%%positions (tap 9) and reactor status (on). 
if Optimisation.t ==1
    Results(Optimisation.t).best_run_solution = ...
        NaN * zeros(1,Optimisation.Nvars);
    
    %%Initialise first tap positions for Ntransformers
    if Optimisation.Ntr ~= 0
        Results(Optimisation.t).best_run_solution(1,Optimisation.tr_pos) = ...
        repmat(1.01,1,Optimisation.Ntr); %1.01 corresponds to the standard tap position
    end
    
    %Initialise initial reactor status for Nreactors
    if Optimisation.Nr ~= 0
        Results(Optimisation.t).best_run_solution(1,Optimisation.r_pos) = ...
        ones(1,Optimisation.Nr); 
    end
else
    %%Initialise the result structs:
    
    %%Fitness and solution vector
    Results(Optimisation.t).Fbest = NaN * zeros(total_runs,1);
    Results(Optimisation.t).Xbest = NaN * zeros(total_runs,Optimisation.Nvars);

    %%Ploss
    Results(Optimisation.t).Ploss = NaN * zeros(total_runs,1);
    Results(Optimisation.t).Ploss_best = NaN;
    Results(Optimisation.t).Ploss_worst = NaN;
    Results(Optimisation.t).Ploss_mean = NaN;

    %%Transformer and reactor status
    Results(Optimisation.t).Tap_changes = NaN * zeros(total_runs,1);
    Results(Optimisation.t).Reactors_changes = NaN * zeros(total_runs,1);
    
    %%Qdistance
    Results(Optimisation.t).extremeness_setpoints = NaN * zeros(total_runs,1);
    
    %%Total cost of optimisation
    Results(Optimisation.t).total_cost_per_run = NaN * zeros(total_runs,1);

    %%Accuracy of Qpcc w.r.t Qref
    Results(Optimisation.t).Qaccuracy = NaN * zeros(total_runs,1);

    %%Times converged and average runtime
    Results(Optimisation.t).Times_converged = NaN * zeros(total_runs,1);
    Results(Optimisation.t).avg_runtime  = NaN * zeros(total_runs,1);
end
end