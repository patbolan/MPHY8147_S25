% Simulate a Lorentzian lineshape 
% This is done by simulating a exponentially decaying complex sinusoid in the time domain and 
% performing an FT to look at the Lorentzian in the frequency domain. 
% For speed, this code uses an analytical expression for the free induction decay. 
% This could also be done using a numerical solution of the Bloch eqns

% Sample in rotating frame (don't need to sample very fast)
dt = 0.0001;
timeax = 0:dt:0.5;

% Now simulate an exponentially decaying sinusoid
T2 = 0.020; % 100 ms 
R2 = 1/T2;
S0 = 1; % Initial signal at t=0
offset = 0; % Hz. Let it be a little off-resonance

% Consider a group of Ni isochromats, each with a discrete off-resonance
% value. The offsets will be normally distributed around zero, with a stdev
% of X Hz
Ni = 4;
spreadHz = 1; 
offsetsHz = randn(1,Ni) .* spreadHz;

% Calculate the time-domain signal in the rotating frame
% Create an array of FID signals
fid_array = zeros(size(timeax,2), Ni);
for idx=1:Ni
   fid_array(:,idx) = S0 .* exp(-timeax .* R2) .* exp(-1j * 2 * pi * offsetsHz(idx) * timeax);     
end

% Sum them and normalize
fid_sum = sum(fid_array,2) / Ni;

% Convert to spectral domain
% Before a FT we scale the first point by 1/2. Long story, will explain
% later.....
tmp = fid_array;
tmp(1,:) = tmp(1,:) * 0.5;
spec_array = fftshift(fft(tmp,[],1)); % Notice the shifting
spec_sum = sum(spec_array,2) / Ni;

% Calculate the frequency domain axis
SW = 1/dt; % Hz
Npoints = size(timeax,2);
dFreq = SW/(Npoints-1); % Hz
freqax = -SW/2:dFreq:SW/2;

% Display results
figure(1)
clf
subplot(2,2,1)
plot(timeax.*1000, abs(fid_sum), 'k', ...
    timeax.*1000, real(fid_sum), 'b');
set(gca, 'ylim', [-.2, 1]);
set(gca, 'xlim', [0, 100]);

legend('abs', 'real');
xlabel('time (ms)')
ylabel('s(t)');
title('Sum, time domain')

subplot(2,2,2)
plot(freqax, real(spec_sum), 'k', ...
    freqax, imag(spec_sum), 'b');
set(gca, 'xlim', [-100 100]); % zoom in a bit
legend('real', 'imag')
xlabel('frequency (Hz)')
ylabel('S(f)');
title('Sum, frequency domain')

subplot(2,2,3)
hold on
for idx=1:Ni
    plot(timeax.*1000, real(fid_array(:,idx)), ':b');
end
plot(timeax.*1000, real(fid_sum), 'k');
set(gca, 'ylim', [-.2, 1]);
set(gca, 'xlim', [0, 100]);
xlabel('time (ms)')
ylabel('s(t)');
hold off

subplot(2,2,4)
hold on
for idx=1:Ni
    tmpfid = fid_array(:,idx);
    tmpfid(1) = tmpfid(1).*0.5; 
    spec = fftshift(fft(tmpfid));    
    plot(freqax, real(spec), ':b');
end
plot(freqax, real(spec_sum), 'k');
set(gca, 'xlim', [-100 100]); % zoom in a bit
xlabel('frequency (Hz)')
ylabel('S(f)');
hold off

zoom on % makes measurement easier

