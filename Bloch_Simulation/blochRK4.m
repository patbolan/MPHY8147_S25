%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% blochRK4
% Calculates the time evolution of magnetization using the Bloch equations.
%
% Inputs:
%  	B1x,B1y - The x and y components of the complex RF amplitude in the
%       rotating frame. A purely amplitude modulated pulse will have all B1Y=0.
%		Units are rad/s [real].
%	deltaT - time step, in seconds [real]
% 	Omega - The off-resonance frequency, in rad/s/ [real]
%	R1 - Spin-lattice relaxation rate (1/T1), seconds^-1 [real]
%	R2 - Spin-spin relaxation rate (1/T2), seconds^-1 [real]
%	Minit - initial magnetization, normalized to 1 (ie equilibirum
%       magnetization M0=1) [Row vector of 3 reals]
%
% Outputs:
%	Mnext - The resultant magnetization. [Row vector of 3 reals]
%
% This implementation uses a "classic" or "original" 4th order Runge-Kutta 
% method to calculate the derivatives. There is no adaptive step size, and
% no checking that the step is small enough. For a specific case, just half
% the step size and see if the result is different (which is what RKF45
% does as part of the algorithm, with a performance penalty). 
%
% This is Kutta's original 4th order version, which is well-described in 
% many texts and websites as well as in: 
%   Fehlberg E, "Low-order classical Runge-Kutta formulas with stepsize 
%   control and their application to some heat transfer problems",  
%   NASA document 19690021375 , NASA-TR-R-315
%
% To calculate this over an RF pulse or pulse sequence, call it iteratively
% for each timestep and for each off-resonance component
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2017 Patrick Bolan and Michael Garwood, University of Minnesota
% We release this code into the public domain, without conditions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Mnext = blochRK4(Minit, B1x, B1y, Omega, R1, R2, deltaT)

% Fourth order Runge-Kutta. The derivatives are iteratively calculated 4 
% times, each with refined estimates. The final result is a weighted sum of
% the intermediate estimates. 
KK = zeros(4,3); % 4 iterations, each holding a 3D vector [Mx, My, Mz]

KK(1,:) = deltaT .* bloch_eqn(B1x, B1y, R1, R2, Omega, Minit);
KK(2,:) = deltaT .* bloch_eqn(B1x, B1y, R1, R2, Omega, Minit + KK(1,:)./2);
KK(3,:) = deltaT .* bloch_eqn(B1x, B1y, R1, R2, Omega, Minit + KK(2,:)./2);
KK(4,:) = deltaT .* bloch_eqn(B1x, B1y, R1, R2, Omega, Minit + KK(3,:));

% Take weighted sum of all K values
Mnext = Minit + 1/6 .* (KK(1,:) + 2.*KK(2,:) + 2.*KK(3,:) + KK(4,:) );

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following function implements the the bloch equations, calculating
% dM/dt from the current value of M considering B1 fields, off-resonance,
% and relaxation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dMdt = bloch_eqn(b1x, b1y, R1, R2, Omega, M)
dMdt = [0 0 0]; % Initize the 3D vector

% The three components of dM/dt: 
dMdt(1) = (Omega * M(2)) - (b1y * M(3)) - (M(1) * R2);
dMdt(2) = (b1x * M(3)) - (Omega * M(1)) - (M(2) * R2);
dMdt(3) = (b1y * M(1)) - (b1x * M(2)) + ((1.0 - M(3)) * R1);

return;




