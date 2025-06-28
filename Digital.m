clear all
clc
close all
pkg load signal
%%% 						Part I							%%%
%% Shihab Mohamed

numOfBits =64;
Data = randi([0 1], 1, numOfBits);  % Random binary data

bitRate = 1;
T= length(Data)/bitRate;
n=200;
N = n*length(Data);
dt = T/N;
t = (0:(N-1)) * dt;
x = zeros(1, N);

%AMI is Alternative mark inversion (bipolar)
last_polarity = -1;  % start with −1 so that the first ‘1’ becomes +1
for i = 0:(length(Data)-1)
  if Data(i+1) == 1
    % flip polarity
    last_polarity = -last_polarity;
    x(i*n + (1:n)) = last_polarity;
  else
    x(i*n + (1:n)) = 0;
  end
end

%%%%%%%%%%%%%%%%%% AMI in Time Domain %%%%%%%%%%%
figure(1)
xlabel('Time');
ylabel('Volt');
grid on;
plot(t,x);

%%%%%%%%%%%%%%%%%% AMI in Freq Domain %%%%%%%%%%%
fs= 1/dt;
df = fs/length(x);
X = fftshift(fft(x))/N;
if(rem(N,2)==0)
  f = (-0.5*fs):df:((0.5*fs)-df);
else
  f = (-0.5*fs+0.5*df):df:(0.5*fs-0.5*df);
end
figure(2)
xlabel('freq');
ylabel('X');
grid on;
plot(f,abs(X));


% Unipolar NRZ reference waveform
%-------------------
z = zeros(1, N);
for i = 0:(length(Data)-1)
  if Data(i+1) == 1
    z(i*n + (1:n)) = +1;
  else
    z(i*n + (1:n)) =  0;
  end
end
figure (5)
xlabel('Time'); ylabel('V');
plot(t,z);

Z = abs(fftshift(fft(z))/N);
figure(6)
xlabel("freq"); ylabel("Z(f)");
plot(f,Z);



%%% 						Part II							%%%
%% 			Mostafa El Hassan
clear;


N = 100;                    % Number of bits
Rb = 1000;                  % Bit rate (bps)
Tb = 1 / Rb;                % Bit duration
Fs = 100000;                % Sampling frequency
Fc = 5000;                  % Carrier frequency (Hz)
t = 0:1/Fs:N*Tb-(1/Fs);     % Time vector
data = randi([0 1], 1, N);  % Random binary data

% Non-return to zero signal because to modulate the data there will be a size
% mismatch between the data and the no. samples of the carrier and ensures
% that the data is held const for its duration Tp
nrz = repelem(data, Fs*Tb);
% Carrier
carrier = cos(2*pi*Fc*t);

% ASK modulation
ask_mod = nrz .* carrier;

% Plot Time-Domain Signal (Temporal Plot)
figure;
plot(t(1:1000), ask_mod(1:1000)); title('ASK Modulated Signal (Time Domain)');
xlabel('Time (s)'); ylabel('Amplitude');

% Plot Spectrum
figure;
Nt = length(ask_mod); % number of time samples of the signal
df = Fs/Nt ;

if rem(Nt, 2) == 0    % even
    f = (-Nt/2 : Nt/2 - 1)* df ;
else                   % odd
    f = (-(Nr- 1)/2 : (Nt- 1)/2)*df ;
end


ASK_FFT = abs(fftshift(fft(ask_mod,Nt)));
plot(f, ASK_FFT); title('ASK Modulated Signal (Frequency Domain)');
xlabel('Frequency (Hz)'); ylabel('|X(f)|');



% Coherent Demodulation
phases = [30, 60, 90];
for i = 1:length(phases)
    phase_rad = deg2rad(phases(i));
    local_osc = cos(2*pi*Fc*t + phase_rad);
    ask_demod = ask_mod .* local_osc;

    % Low-pass filter using moving average
    lpf_output = filter(ones(1, Fs*Tb), 1, ask_demod);
    sample_points = Fs*Tb*(0:N-1) + round(Fs*Tb/2);     
    demod_bits = lpf_output(sample_points) > 0.5;

    % Plot demodulated signal (optional)
    figure;
    plot(t(1:1000), ask_demod(1:1000)); title(['ASK Demodulated Signal with Phase = ', num2str(phases(i)), '°']);
    xlabel('Time (s)'); ylabel('Amplitude');

    % Bit Error Rate
    ber = sum(demod_bits ~= data) / N;
    fprintf('ASK with phase %d°: Bit Error Rate = %.2f%%\n', phases(i), ber*100);
end


