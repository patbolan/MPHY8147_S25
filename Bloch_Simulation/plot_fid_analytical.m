% Plot the analytical expression for a simple FID

% Consider a 0.5T MR system
f0 = 42.58E6 * 0.5; % 1H Larmor Frequency at 0.5T, Hz
dT = 0.0001; % Sampling period, in seconds
BW = 1/dT; % Hz
time_s = 0:dT:0.1; % time = 0 --> 100ms 

% Now simulate an exponentially decaying sinusoid
T2 = 0.020; % 20 ms 
S0 = 1; % Initial signal at t=0
offset = 100; % Hz Let it be a little off-resonance

% Signal in rotating frame
s = S0 .* exp(-time_s ./ T2) .* exp(-1j * 2 * pi * offset * time_s);

% Common to plot real and a
figure(1)
plot(time_s, real(s), 'k', ...
    time_s, imag(s), 'b', ...
    time_s, abs(s), '-r')
xlabel('time (s)')
ylabel('signal')
legend('real', 'imaginary', 'abs')

figure(2)
plot3(time_s, real(s), imag(s))
xlabel('time (s)')
ylabel('real(s)')
zlabel('imaginary(s)')

