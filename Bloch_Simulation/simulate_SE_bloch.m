% This code builds on the inhomogeneous FID simulation to create a spin
% echo using the full bloch equations.

% Specify the RF pulses
Tp = 0.005; % in s, a short pulse
flipAngleExc = 90; % Specified flip angle, degrees
flipAngleRefoc = 180; % Specified flip angle, degrees

% Configure the the time axis to simulate
dT = 0.00001; % Time resolution. Needs to be finer than the RF pulsewidth
timeax = 0:dT:0.5; % Time axis for simulation
Nt = size(timeax,2); % number of timepoints
TE = 0.200; % Specify the echo time

% Add T1 and T2
T1 = 1000; % s, Making it very long, essentially neglecting T1 relaxation
T2 = 0.300; % 20 ms

% Define the B1 function.
% To keep the refocused magnetization along Mx, I'll use a -90y then a 180x
B1maxExc = flipAngleExc/180 * 500 * (0.001/Tp) * 2 * pi;
B1maxRefoc = flipAngleRefoc/180 * 500 * (0.001/Tp) * 2 * pi;
B1x = zeros(Nt,1);
B1y = zeros(Nt,1);

% I'll do this in a loop to be clearer (but slower than indexing)
for tdx=1:Nt
    timept = timeax(tdx);
    
    % The 90 degree pulse goes from t=[0,Tp]
    if timept>0 && timept<=Tp
       B1y(tdx) = -B1maxExc; 
    end

    % The 180 goes from t=[TE/2, Tp+(TE/2)]
    if timept>TE/2 && timept<=(Tp+(TE/2))
        B1x(tdx) = B1maxRefoc;
    end 
end

% Convert to complex 
B1 = B1x + 1j.*B1y;

% Consider a set of Ni isochromats, each with a discrete off-resonance
% value. The offsets will be normally distributed around zero, with a stdev
% of X Hz
Ni = 8;
spreadHz = 4; 
offsets = rand(1,Ni) .* 2 * pi * spreadHz;

% Uniformly distributed
offsets = linspace(-spreadHz/2, spreadHz/2, Ni) .* (2*pi);

% All isochromats have the same starting magnetization, Mz=M0=1
M0 = [0 0 1];

% Loop over all time and simulate magnetization
% We will track the magnetization for each isochromat, each timepoint, and
% each of the three cartesian components
Mti = zeros(Nt,Ni,3); % Track a 3D vector over all Nt timepoints
Mti(1,:,1) = M0(1); % Initial conditions for all isochromats
Mti(1,:,2) = M0(2); 
Mti(1,:,3) = M0(3); 

fprintf('Starting Bloch simulation...');
for tdx=2:Nt
    % Calculate each isochromat separately
    for idx=1:Ni
        % Calculate magnetization at the next timepoint using bloch simulation
        % Note: this fn expects a 1x3 and chokes on a 3x1, need to select
        % the Mtmp from this isochromat and this timepoint 
        Mtmp = squeeze(Mti(tdx-1,idx,:)).'; 
        Mti(tdx,idx,:) = blochRK4(Mtmp, B1x(tdx), B1y(tdx), offsets(idx), 1/T1, 1/T2, dT);
    end
end
fprintf(' done.\n');

% Let's observe the entire signal throughout the experiment
Mxyi = Mti(:,:,1) + 1j*Mti(:,:,2);
Mzi = Mti(:,:,3);

% The total signal is the vector average of all the individual signals
Mxyt = sum(Mxyi,2)./Ni;
Mzt = sum(Mzi,2) ./ Ni;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display the results

% So here's the pulse sequence
figure(1)
clf
subplot(3,1,1)
plot(timeax, abs(B1), 'b');
xlabel('time (s)');
ylabel('gamma-B1 (Hz)')

% Transferse Magnetization, magnitude
subplot(3,1,2)
hold on
% Plot the signal from each of the isochromats in blue
for idx=1:Ni
    plot(timeax, abs(Mxyi(:,idx)), 'b');
end

% ... and the sum in black
plot(timeax, abs(Mxyt), 'k');
hold off
xlabel('time (s)');
ylabel('abs{Mxy(t)}');

% Transferse Magnetization, phase
subplot(3,1,3)
hold on
% Plot the signal from each of the isochromats in blue
for idx=1:Ni
    plot(timeax, (angle(Mxyi(:,idx))), 'b');
end

% ... and the sum in black
plot(timeax, (angle(Mxyt)), 'k');
hold off
%set(gca, 'ylim', [-pi pi].*1.2);
xlabel('time (s)');
ylabel('angle{Mxy(t)}');



