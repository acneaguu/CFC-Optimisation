%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%README:
%%This function is used to determine the fitness of a given set of
%%solutions specified by 'Xin'. The fitness consists of two parts: the
%%objective function and the constraint penalty. The first part is the
%%objective function which is the main factor for the optimisation. This
%%term is relevant when a set of solutions is feasible i.e. do not violate
%%the set constraints. The constraint penalty should penalise constraint
%%violations. The more constraints are violated, the higher the penalty
%%should be and thus the fitness is worse. 

%%'Xin' can be both a vector (only one particle) or a matrix (where each
%%particle is represented by a different row. This function returns the
%%fitness, the OF, the composition of the violated constraints and the
%%(rounded) control variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F,OF,g,Xout] = fitness_eval(Xin)
global Keeptrack mpopt Systemdata PFresults Optimisation FCount 

%%Initialise outputs
Xout = NaN * ones(size(Xin));
F = NaN * ones(size(Xin,1),1);

%%Determine number of inputs i.e. number of different solutions
NXin = size(Xin,1);

%%For each proposed solution:
for np = 1:NXin
    %%Increment the function evaluation counter
    FCount = FCount+1;
    %% Run powerflow
    %%Round discrete variables of Xin
    Xout(np,:) = round_discrete_vars(Xin(np,:));
    
    %%Change topology according to solutions
    update_casefile(Xout(np,:),1);
    
    %%Run pf on the system
    PFresults = runpf(Systemdata.mpc,mpopt);
    
    %%If the power flow converged:
    if PFresults.success == 1
        %------------------------------------------------------------------
        %%CONSTRAINTS:
        [g, total_violations,composition] = compute_violation_constraints_v3();
        %------------------------------------------------------------------
        %%Objective function:
        OF = compute_costs_v2(Xout(np,:));
        %------------------------------------------------------------------
    
    
        %%Consider only the OF if there are no violations and otherwise
        %%consider only the constraint violations.
        if total_violations == 0
            %%Feasible solution
            F = OF;                    
        else
            %%Infeasible solution
            F = total_violations*1e20;
        end
    else
        %%Big penalty if powerflow runs are unsuccesful
        g = 100*ones(2*(Systemdata.Nbus+Systemdata.Nbranch)+1,1);
        composition = [2*Optimisation.p1*Systemdata.Nbus,Optimisation.p2,...
            2*Optimisation.p3*Systemdata.Nbranch];
        OF = 1e50;
        F = 1e50; 
    end
    
    %%Update the fitness progress archive:
    %%Keeptrack.Fitness keeps track of fitness of every particle
    Keeptrack.Fitness(FCount) = F; 
    
    %%Keeps track of the solutions corresponding to the fitnesses
    Keeptrack.solution(FCount,:)= Xout(np,:);
    
    %%Keeps track of overall best fitnesses and solutions. A solution is
    %%considered better if its fitness is smaller or equal to the previous
    %%fitness.
    if FCount > 1
        if Keeptrack.Fitness(FCount) <= Keeptrack.FitBest(FCount-1)
            Keeptrack.FitBest(FCount) = Keeptrack.Fitness(FCount);
            Keeptrack.SolBest(FCount,:) = Keeptrack.solution(FCount,:); 
            Keeptrack.violation_composition(FCount,:) = composition;
        else
            Keeptrack.FitBest(FCount) = Keeptrack.FitBest(FCount-1);
            Keeptrack.SolBest(FCount,:) = Keeptrack.SolBest(FCount-1,:);
            Keeptrack.violation_composition(FCount,:) = ...
                Keeptrack.violation_composition(FCount-1,:);
        end
    else
        Keeptrack.FitBest(FCount) = Keeptrack.Fitness(FCount);
        Keeptrack.SolBest(FCount,:) = Keeptrack.solution(FCount,:); 
        Keeptrack.violation_composition(FCount,:) = composition;
    end
end

%%For MVMO-SHM: update the counter of the MVMO function with the value of
%%FCount i.e. the internal timer.
if Optimisation.algorithm == 4
global proc %#ok<TLEV>
    proc.i_eval = FCount;
    if proc.i_eval>=proc.n_eval 
        proc.finish=1;  
    end
end

%%Prints the progress if desired at a specified interval
if Optimisation.print_progress == 1 
    if (FCount == 1) || (FCount == Optimisation.Neval) || mod(FCount,Optimisation.print_interval) == 0
      fprintf('Neval: %7d,   fitness: %12.7E \n',...
          FCount, Keeptrack.FitBest(FCount))
    end
end
end
