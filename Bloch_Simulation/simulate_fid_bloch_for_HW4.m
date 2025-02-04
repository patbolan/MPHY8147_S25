% Simulation of the fid experiment, in rotating frame, using a numerical
% solver for the bloch equations

% Specify the RF pulse
Tp = 0.005; % a 100 us pulse. 
% Note if you use a longer pulse you can see off-resonance and relaxation
% effects
flipAngle = 90; % Specified flip angle, degrees

% Configure the the time axis to simulate
dT = 0.00001; % Time resolution. Needs to be finer than the RF pulsewidth
timeax = 0:dT:Tp; % Time axis for simulation
Nt = size(timeax,2); % number of timepoints

% Add T1 and T2 (set them to ~1000s to disable relaxation)
T1 = 1E10; % typ ~1s
T2 = 1E10; % typ is 20 ms

% Consider a single isochromat, starting at equilibrium (Mz=1)
offset = 30 * 2 * pi; % Off-resonance term
M0 = [0 0 1];

% Define the B1 function. i want a -90y pulse, which will rotate M down to
% the x-axis (M=[1 0 0]). Calculate the B1max needed, in rad/s. Note 
%   500 Hz in 1ms gives 180 degrees, which is how I calibrate
B1max = flipAngle/180 * 500 * (0.001/Tp) * 2 * pi;
B1y = zeros(Nt,1); % Set B1 to all zeros, same size as time axis
B1y(timeax<=Tp) = -B1max; % From t=[0,Tp] set B1y to -B1max (ie square pulse)
B1x = B1y.*0; % Set the B1x component to zero

% Define the ADC trace. Turn it on (set to 1) after the RF pulse
adc = timeax.*0;
adc(timeax>Tp) = 1; % Set it to 0 when off (during RF) and 1 after

% Loop over all time and simulate magnetization
Mt = zeros(Nt,3); % Track a 3D vector over all Nt timepoints
Mt(1,:) = M0; % Initial conditions
fprintf('Starting Bloch simulation...');
for tdx=2:Nt
    % Calculate magnetization at the next timepoint using bloch simulation
    Mt(tdx,:) = blochRK4(Mt(tdx-1,:), B1x(tdx), B1y(tdx), offset, 1/T1, 1/T2, dT);
    %Mt(tdx,:) = blochODE45(Mt(tdx-1,:), B1x(tdx), B1y(tdx), offset, 1/T1, 1/T2, dT);
    %Mt(tdx,:) = blochRotMatrix(Mt(tdx-1,:), B1x(tdx), B1y(tdx), offset, 1/T1, 1/T2, dT);
end
fprintf(' done.\n');
%
% Animate M
figure(1)
% The full trace takes a while, even with a tiny delay
animate_Mvector(Mt, .001);







