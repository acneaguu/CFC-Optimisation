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
%CITE "Solving security constrained optimal power flow problems: a hybrid
%evolutionary approach". Applied Intelligence, Springer, pp.: 1-19. 2018
%and
%"Applying C-DEEPSO to solve Large Scale Global Optimization Problems"
%Proc. IEEE on Congress on Evolutionary Computation (CEC), 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MUTATE_WEIGHTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ mutated_Weights ] = MUTATE_WEIGHTS( weights, mutationRate )
% Mutate weights & check weights limits
mutated_Weights = weights;
[ ~, numWeights ] = size( weights );
if numWeights > 4
    for i = 1 : 5
        mutated_Weights( i ) = weights( i ) + normrnd( 0, 1 ) * mutationRate;
        if  mutated_Weights( i ) > 1
            mutated_Weights( i ) = 1;
        elseif  mutated_Weights( i ) < 0
            mutated_Weights( i ) = 0;
        end
    end
    mutated_Weights( 6 ) = weights( 6 ) + normrnd( 0, 1 ) * mutationRate;
    if  mutated_Weights( 6 ) > 2
        mutated_Weights( 6 ) = 2;
    elseif  mutated_Weights( 6 ) < 0
        mutated_Weights( 6 ) = 0;
    end
else
    for i = 1 : 4
        mutated_Weights( i ) = weights( i ) + normrnd( 0, 1 ) * mutationRate;
        if  mutated_Weights( i ) > 1
            mutated_Weights( i ) = 1;
        elseif  mutated_Weights( i ) < 0
            mutated_Weights( i ) = 0;
        end
    end
end
end
