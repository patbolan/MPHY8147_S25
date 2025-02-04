% Simulation of the fid experiment, in rotating frame, using a numerical
% solver for the bloch equations

% Specify the RF pulse
Tp = 0.0001; % a 100 us pulse. 
% Note if you use a longer pulse you can see off-resonance and relaxation
% effects
flipAngle = 90; % Specified flip angle, degrees

% Configure the the time axis to simulate
dT = 0.000001; % Time resolution. Needs to be finer than the RF pulsewidth
timeax = 0:dT:0.5; % Time axis for simulation
Nt = size(timeax,2); % number of timepoints

% Add T1 and T2 (set them to ~1000s to disable relaxation)
T1 = 1.0; % typ ~1s
T2 = 0.020; % typ is 20 ms

% Consider a single isochromat, starting at equilibrium (Mz=1)
offset = 100 * 2 * pi; % Off-resonance term
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

% Animate M
figure(1)
% The full trace takes a while, even with a tiny delay
%animate_Mvector(Mt, .001);

% Slow version to see the excitation
%animate_Mvector(Mt(1:200,:), .01, true);
%%

% Plot results
figure(1)
ax1=subplot(4,1,1);
plot(timeax.*1000, -B1y ./(2*pi));
ylabel('-B1y (Hz)');
%maxY = max(-B1y(:))/(2*pi);
%set(gca, 'ylim', [-0.1 1.1]* maxY);
xlabel('time (ms)');

ax2=subplot(4,1,2);
plot(timeax.*1000, Mt(:,1));

set(gca, 'ylim', [-1 1])
ylabel('Mx');
xlabel('time (ms)');

ax3=subplot(4,1,3);
plot(timeax.*1000, Mt(:,2));
set(gca, 'ylim', [-1 1])
ylabel('My');
xlabel('time (ms)');

ax4=subplot(4,1,4);
plot(timeax.*1000, Mt(:,3));
set(gca, 'ylim', [-1 1])
ylabel('Mz');
xlabel('time (ms)');

linkaxes([ax1, ax2, ax3, ax4],'x')

%% View Spectrum
% Take the Fourier Transform of the transverse magnetization
% complex Mxy = Mx * i My;
Mxy = Mt(:,1) + 1j*Mt(:,2); % Complex, transverse magnetization

% Only take the points after the RF pulse. Experimentally, we cannot
% measure during the RF pulse; we wait a little after the RF is turned off
% and then turn on our receiver. We'll call this signal the FID. 

% Take only time points while the ADC is turned on (after RF pulse)
pointsSampled = find(adc==1);
fid = Mxy(pointsSampled);
timeOfFid = timeax(pointsSampled);

% Before a FT we scale the first point of an FID by 1/2. This is a subtle 
% effect that is often ignored. (G Otting et al, JMR 66 p187 1986) 
tmpfid = fid;
tmpfid(1) = tmpfid(1).*0.5; 
spec = fftshift(fft(tmpfid));

% Calculate frequency axis
SW = 1/dT; % Hz
Npoints = size(timeOfFid,2);
dFreq = SW/(Npoints-1); % Hz
freqax = (-SW/2:dFreq:SW/2)-dFreq/2;

figure(2)
plot(freqax, real(spec), 'k', ...
    freqax, imag(spec), 'b')
set(gca, 'xlim', [-200 200])
legend('real', 'imag')
xlabel('frequency (Hz)')
ylabel('S(f)');








