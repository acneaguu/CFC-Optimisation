%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Leonel Carvalho, PhD (email: leonel.m.carvalho@inescporto.pt)
% Carolina Marcelino, PhD (email: carolimarc@cos.ufrj.br)
% April 02, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright (c) 2018, Leonel Carvalho and Carolina Marcelino
%All rights reserved.
%
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are
%met:
%
%    * Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in
%      the documentation and/or other materials provided with the 
%      distribution
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%POSSIBILITY OF SUCH DAMAGE.
%CITE "Applying C-DEEPSO to solve Large Scale Global Optimization Problems"
%Proc. IEEE on Congress on Evolutionary Computation (CEC), 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C-DEEPSO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ gbestval, gbest ] = CDEEPSO_algorithm(fun,lb,ub)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETTING PARAMETERS
% Global variable
global cdeepso_par;
global ff_par;

cdeepso_par.fun = fun;

% Weights matrix
% 1 - inertia
% 2 - memory
% 3 - cooperation
% 4 - perturbation
weights = rand( cdeepso_par.popSize, 5 );
weights( :, 6 ) = 2 * rand( cdeepso_par.popSize, 1 );
if cdeepso_par.typeCDEEPSO == 1
    % Sign of the perception term
    signPerception = zeros( cdeepso_par.popSize, ff_par.D );
else
    % Sign of the perception term
    signPerception = ones( cdeepso_par.popSize, ff_par.D );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RANDOMLY INITIALIZE CURRENT population
% Particles' lower bounds
Xmin = lb;
% Particles' upper bounds
Xmax = ub;
Vmin = -Xmax + Xmin;
Vmax = -Vmin;
pos = zeros( cdeepso_par.popSize, ff_par.D );
vel = zeros( cdeepso_par.popSize, ff_par.D );
for i = 1 : cdeepso_par.popSize
    pos( i, : ) = Xmin + ( Xmax - Xmin ) .* rand( 1, ff_par.D );
    vel( i, : ) = Vmin + ( Vmax - Vmin ) .* rand( 1, ff_par.D );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EVALUATE the CURRENT population
[ fit ] = FITNESS_FUNCTION( cdeepso_par.popSize, pos );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UPDATE GLOBAL BEST
[ gbestval, gbestid ] = min( fit );
gbest = pos( gbestid, : );
memGbestval = zeros( 1, ( cdeepso_par.maxGen + 1 ) );
memGbestval( 1 ) = gbestval;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Memory of the CDEEPSO
if cdeepso_par.memGBestMaxSize > 0
    cdeepso_par.memGBestSize = 1;
    memGBest( cdeepso_par.memGBestSize, : ) = gbest;
    memGBestFit( 1, cdeepso_par.memGBestSize ) = gbestval;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UPDATE INDIVIDUAL BEST
% Individual best position ever of the particles of CURRENT population
myBestPos = pos;
% Fitness of the individual best position ever of the particles of CURRENT population
myBestPosFit = fit;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE generation counter
countGen = 0;
countGenWoChangeBest = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOOP until termination criteria isn't met
while countGen < cdeepso_par.maxGen && countGenWoChangeBest <= cdeepso_par.maxGenWoChangeBest && ff_par.fitEval <= cdeepso_par.maxFitEval
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % UPDATE generation counter
    countGen = countGen + 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % UPDATE MEMORY according to DEEPSO strategy
    if cdeepso_par.strategyCDEEPSO == 1 || cdeepso_par.strategyCDEEPSO == 3
        tmpMemGBestSize = cdeepso_par.popSize;
        tmpMemGBestFit = fit;
        tmpMemGBest = pos;
    end
    if cdeepso_par.strategyCDEEPSO == 2 || cdeepso_par.strategyCDEEPSO == 4
        tmpMemGBestSize = cdeepso_par.memGBestSize;
        tmpMemGBestFit = memGBestFit;
        tmpMemGBest = memGBest;
    end
    if cdeepso_par.strategyCDEEPSO == 5
        tmpMemGBestSize = cdeepso_par.memGBestSize + cdeepso_par.popSize;
        tmpMemGBestFit = cat( 2, memGBestFit, fit );
        tmpMemGBest = cat( 1, memGBest, pos );
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if cdeepso_par.typeCDEEPSO == 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % UPDATE PERSONAL BEST Prof. Vladimiro
        for i = 1 : cdeepso_par.popSize
            [ tmpMyBestPos, tmpSignPerception ] = DEEPSO_COMPUTE_NEW_PERSONAL_BEST_PROF_VLAD( ff_par.D, fit( i ), tmpMemGBestSize, tmpMemGBest, tmpMemGBestFit );
            myBestPos( i, : ) = tmpMyBestPos;
            signPerception( i, : ) = tmpSignPerception;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % UPDATE PERSONAL BEST DE
        for i = 1 : cdeepso_par.popSize
            [ tmpMyBestPos, tmpMyBestPosFit ] = CDEEPSO_COMPUTE_NEW_PERSONAL_BEST( ff_par.D, weights( i, 5 ), weights( i, 6 ), ...
                myBestPos( i, : ), myBestPosFit( i ), gbest, tmpMemGBestSize, tmpMemGBest, tmpMemGBestFit, Xmin, Xmax );
            myBestPos( i, : ) = tmpMyBestPos;
            myBestPosFit( i ) = tmpMyBestPosFit;
            if myBestPosFit( i ) < fit( i )
                pos( i, : ) = myBestPos( i, : );
                fit( i ) = myBestPosFit( i );
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % UPDATE GLOBAL BEST
        [ tmpgbestval, gbestid ] = min( fit );
        if tmpgbestval < gbestval
            gbestval = tmpgbestval;
            gbest = pos( gbestid, : );
            % UPDATE MEMORY DEEPSO
            if cdeepso_par.memGBestMaxSize > 0
                if cdeepso_par.memGBestSize < cdeepso_par.memGBestMaxSize
                    cdeepso_par.memGBestSize = cdeepso_par.memGBestSize + 1;
                    memGBest( cdeepso_par.memGBestSize, : ) = gbest;
                    memGBestFit( 1, cdeepso_par.memGBestSize ) = gbestval;
                else
                    [ ~, tmpgworstid ] = max( memGBestFit );
                    memGBest( tmpgworstid, : ) = gbest;
                    memGBestFit( 1, tmpgworstid ) = gbestval;
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COPY CURRENT population
    copyPos = pos;
    copyVel = vel;
    copyWeights = weights;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % APPLY DEEPSO movement rule
    for i = 1 : cdeepso_par.popSize 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MUTATE WEIGHTS of the particles of the COPIED population
        copyWeights( i, : ) = MUTATE_WEIGHTS( weights( i, : ), cdeepso_par.mutationRate );
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % COMPUTE NEW VELOCITY for the particles of the COPIED population
        copyVel( i, : ) = CDEEPSO_COMPUTE_NEW_VEL( ff_par.D, copyPos( i, : ), myBestPos( i, : ), gbest, ...
            copyVel( i, : ), Vmin, Vmax, copyWeights( i, : ), cdeepso_par.communicationProbability, signPerception( i, : ) ) ;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % COMPUTE NEW POSITION for the particles of the COPIED population
        [ copyPos( i, : ), copyVel( i, : ) ] = COMPUTE_NEW_POS( copyPos( i, : ), copyVel( i, : ) );
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % COMPUTE NEW VELOCITY for the particles of the CURRENT population
        vel( i, : ) = CDEEPSO_COMPUTE_NEW_VEL( ff_par.D, pos( i, : ), myBestPos( i, : ), gbest, ...
            vel( i, : ), Vmin, Vmax, weights( i, : ), cdeepso_par.communicationProbability, signPerception( i, : ) );
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % COMPUTE NEW POSITION for the particles of the CURRENT population
        [ pos( i, : ), vel( i, : ) ] = COMPUTE_NEW_POS( pos( i, : ), vel( i, : ) );
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ENFORCE search space limits of the COPIED population
    [ copyPos, copyVel ] = ENFORCE_POS_LIMITS( ff_par.D, cdeepso_par.popSize, copyPos, Xmin, Xmax, copyVel, Vmin, Vmax );
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ENFORCE search space limits of the CURRENT population
    [ pos, vel ] = ENFORCE_POS_LIMITS( ff_par.D, cdeepso_par.popSize, pos, Xmin, Xmax, vel, Vmin, Vmax );
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EVALUATE the COPIED population
    [ copyFit ] = FITNESS_FUNCTION( cdeepso_par.popSize, copyPos );
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EVALUATE the CURRENT population
    [ fit ] = FITNESS_FUNCTION( cdeepso_par.popSize, pos );
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CREATE NEW population to replace CURRENT population
    selParNewSwarm = ( copyFit < fit );
    for i = 1 : cdeepso_par.popSize
        if selParNewSwarm( i )
            fit( i ) = copyFit( i );
            pos( i, : ) = copyPos( i, : );
            vel( i, : ) = copyVel( i, : );
            weights( i, : ) = copyWeights( i, : );
        end
        if fit( i ) < myBestPosFit( i )
            myBestPos( i, : ) = pos( i, : );
            myBestPosFit( i ) = fit( i );
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % UPDATE GLOBAL BEST
    [ tmpgbestval, gbestid ] = min( fit );
    if tmpgbestval < gbestval
        gbestval = tmpgbestval;
        gbest = pos( gbestid, : );
        % UPDATE MEMORY DEEPSO
        if cdeepso_par.memGBestMaxSize > 0
            if cdeepso_par.memGBestSize < cdeepso_par.memGBestMaxSize
                cdeepso_par.memGBestSize = cdeepso_par.memGBestSize + 1;
                memGBest( cdeepso_par.memGBestSize, : ) = gbest;
                memGBestFit( 1, cdeepso_par.memGBestSize ) = gbestval;
            else
                [ ~, tmpgworstid ] = max( memGBestFit );
                memGBest( tmpgworstid, : ) = gbest;
                memGBestFit( 1, tmpgworstid ) = gbestval;
            end
        end
        countGenWoChangeBest = 0;
    else
        countGenWoChangeBest = countGenWoChangeBest + 1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SAVE fitness
    memGbestval( countGen + 1 ) = gbestval;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PRINT results
    if rem( countGen, cdeepso_par.printConvergenceResults ) == 0 || countGen == 1
        fprintf('Gen: %-8d Best Fit: %.3e\n', countGen, gbestval );
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
if rem( countGen, cdeepso_par.printConvergenceResults ) ~= 0
    fprintf('Gen: %-8d Best Fit: %.3e\n', countGen, gbestval );
end
fprintf('\n');
if cdeepso_par.printConvergenceChart == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PRINT final results
    x = 0 : 1 : countGen;
    memGbestval = memGbestval( 1 : countGen + 1 );
    plot( x, memGbestval);
    xlabel('generation');
    ylabel('fitness');
    grid on;
    print('Convergencia','-dpng');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
end