%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%README:
%%By: Dr. Jose L. Rueda
%%Date: 02.03.2012
%%Edited: Alex Neagu & Jinhan Bai
%%Date: 11/05/2020

%%This function is an edited version of an example made by Dr. J. Rueda.
%%This version is adapted to initialise the needed parameters for the
%%MVMO-SHM algorithm using the main. Call this function before running the
%%MVMO-SHM algorithm.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initialise_mvmoshm()
global ps proc parameter %%vars of MVMO-SHM
global Optimisation

ps.D=Optimisation.Nvars;                     %Dimension of optimization problem      
proc.n_eval = Optimisation.Neval;            %Maximum allowed evaluations

%%Strategic parameters MVMO
parameter.n_par=Optimisation.Populationsize; %Number of particles 
parameter.n_tosave=3;                        %Archive size
parameter.fs_factor_start=1;                 %Initial fs-factor 
parameter.fs_factor_end= 1;                  %Final fs-factor
parameter.ratio_gute_max=0.3;                %Initial portion of good particles    
parameter.ratio_gute_min=0.3;                %Final portion of good particles
parameter.n_random_ini =round(0.54*ps.D/1.0);%Initial number of variables selected for mutation 
parameter.n_random_last=round(0.31*ps.D/1.0);%Final number of variables selected for mutation 
parameter.local_prob= 0;                     %Probability value between 0 and 1. 
                                             %Set to 0 to deactivate local search
parameter.ratio_local = 0.09;                %Percentage of total runs without local search possibility
end
