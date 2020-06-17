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
% C-DEEPSO_COMPUTE_NEW_VEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ new_vel ] = CDEEPSO_COMPUTE_NEW_VEL( D, pos, myBestPos, gbest, vel, Vmin, Vmax, weights, communicationProbability, signPerception )
% Computes new velocity according to the movement rule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute inertial term
inertiaTerm = weights( 1 ) * vel;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute memory term
memoryTerm = zeros( 1, D );
for i = 1 : D
    memoryTerm( 1, i ) = signPerception( 1, i ) * weights( 2 ) * ( myBestPos( 1, i ) - pos( 1, i ) );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute cooperation term
% Sample normally distributed number to perturbate the best position
cooperationTerm = weights( 3 ) * ( gbest * ( 1 + weights( 4 ) * normrnd( 0, 1 ) ) - pos );
communicationProbabilityMatrix = rand( 1, D ) < communicationProbability;
cooperationTerm = cooperationTerm .* communicationProbabilityMatrix;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute velocity
new_vel = inertiaTerm + memoryTerm + cooperationTerm;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check velocity limits
new_vel = ( new_vel > Vmax ) .* Vmax + ( new_vel <= Vmax ) .* new_vel;
new_vel = ( new_vel < Vmin ) .* Vmin + ( new_vel >= Vmin ) .* new_vel;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end