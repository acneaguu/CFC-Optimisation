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
% C-DEEPSO_COMPUTE_NEW_PERSONAL_BEST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ tmpMyBestPos, tmpMyBestPosFit ] = CDEEPSO_COMPUTE_NEW_PERSONAL_BEST( D, CR, F, myBestPos, myBestPosFit, gbest, numGBestSaved, memGBest, memGBestFit , ...
    Xmin, Xmax )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select subset of particles to sample myBestPos from
% Get the index of the best particles ever visited that have a fitness less
% than or equal to the fitness of particle i
global cdeepso_par;
tmpMyBestPos = myBestPos;
tmpMyBestPosFit = myBestPosFit;
tmpMemoryVect = zeros( 1, numGBestSaved );
tmpMemoryVectSize = 0;
for i = 1 : numGBestSaved
    if ( memGBestFit( 1, i ) < myBestPosFit ) && ( ~isequal( myBestPos, memGBest( i, :) ) );
        tmpMemoryVectSize = tmpMemoryVectSize + 1;
        tmpMemoryVect( 1, tmpMemoryVectSize ) = i;
    end
end
tmpMemoryVect = tmpMemoryVect( 1, 1:tmpMemoryVectSize );
if cdeepso_par.typeCDEEPSO == 2
    % DE/Rand/1/Bin
    if tmpMemoryVectSize >= 3
        tmpIndexMemoryVect = randsample( tmpMemoryVect, 3, false );
        tmpMyBestPos = memGBest( tmpIndexMemoryVect( 1 ), : ) + F * ( memGBest( tmpIndexMemoryVect( 2 ), : ) - memGBest( tmpIndexMemoryVect( 3 ), : ) );
        tmpIndexD = randsample( D, 1 );
        tmpRand = rand( 1, D );
        for i = 1 : D
            if ~( ( tmpRand( i ) < CR ) || ( i == tmpIndexD ) )
                tmpMyBestPos( i ) = myBestPos( i );
            end
            % check pos limits
            if tmpMyBestPos( i ) < Xmin( i )
                tmpMyBestPos( i ) = Xmin( i );
            elseif tmpMyBestPos( i ) > Xmax( i )
                tmpMyBestPos( i ) = Xmax( i );
            end
        end
        % select the position to use as memory
        tmpMyBestPosFit = FITNESS_FUNCTION( 1, tmpMyBestPos );
        if( tmpMyBestPosFit > myBestPosFit )
            tmpMyBestPos = myBestPos;
            tmpMyBestPosFit = myBestPosFit;
        end
    end
elseif cdeepso_par.typeCDEEPSO == 3
    % DE/Best/1/Bin
    if tmpMemoryVectSize >= 2
        tmpIndexMemoryVect = randsample( tmpMemoryVect, 2, false );
        tmpMyBestPos = gbest + F * ( memGBest( tmpIndexMemoryVect( 1 ), : ) - memGBest( tmpIndexMemoryVect( 2 ), : ) );
        tmpIndexD = randsample( D, 1 );
        tmpRand = rand( 1, D );
        for i = 1 : D
            if ~( ( tmpRand( i ) < CR ) || ( i == tmpIndexD ) )
                tmpMyBestPos( i ) = gbest( i );
            end
            % check pos limits
            if tmpMyBestPos( i ) < Xmin( i )
                tmpMyBestPos( i ) = Xmin( i );
            elseif tmpMyBestPos( i ) > Xmax( i )
                tmpMyBestPos( i ) = Xmax( i );
            end
        end
        % select the position to use as memory
        tmpMyBestPosFit = FITNESS_FUNCTION( 1, tmpMyBestPos );
        if( tmpMyBestPosFit > myBestPosFit )
            tmpMyBestPos = myBestPos;
            tmpMyBestPosFit = myBestPosFit;
        end      
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
