clear;
clc;

%%%%%%%%   definitions  %%%%%%%%%% Eiad Hassan
%%
fs = 100;           % sampling frequency in Hz
df = 0.01;          % desired frequency resolution
N =ceil (fs/df);    % number of points needed
T = N/fs;           % total time
ts = 1/fs;          % time step
t=-50.495:ts:((N-1)*ts-50.495);
if(rem(N,2)==0) %% Even
  f = - (0.5*fs) : df : (0.5*fs-df) ; %% Frequency vector if x/f is even
else %% Odd
  f = - (0.5*fs-0.5*df) : df : (0.5*fs-0.5*df) ; %% Frequency vector if x/f is odd
end
x1= zeros(size(t));
x2= zeros(size(t));
x3= zeros(size(t));
x2(t >= -4 & t <= 0) = t(t >= -4 & t <= 0) + 5;
x3(t > 0 & t <= 4) = -t(t > 0 & t <= 4) + 5;
x=x1+x2+x3; %%%%%%%%  x(t) %%%%%%%%%

%%%%%%%%   Q1  %%%%%%%%%%

figure;
plot(t,x);
xlim([-5,5]);
xlabel('t');
ylabel('x(t)');
grid on;
legend('x(t)');

%%%%%%%%     Q2&Q3  %%%%%%%%%%    analytical: i wrote it in a paper with steps

X_num = fftshift(fft(x)) * ts; %% * ts -> non-periodic
figure
X = (sin(8*pi*f) ./ (pi*f)) + (sin(4*pi*f) ./ (pi*f)).* (sin(4*pi*f) ./ (pi*f));
plot(f,abs(X_num),'b');
hold on
plot(f,abs(X),'r--','LineWidth', 2);
hold off
xlabel('Frequency (Hz)');
ylabel('|X(f)|');
grid on;
legend('FT computed','FT analytical');

%%%%%%%%     Q4  %%%%%%%%%%
X_value=abs(X_num);
X_squ=X_value .^2;
max_amplitude=max(X_squ);
max_amplitude_index=find(X_squ==max_amplitude);
last_index=length(X_squ);
for i=max_amplitude_index:last_index
  if(abs(X_squ(i))<=0.05*max_amplitude)
    BW=f(i);
    break;
  end
end
fprintf('5%% power bandwidth of x(t) ≈ %.4f Hz\n', BW);




%------------------------------------------------------------------------------------------------------
%%%%%%%%    Q5  %%%%%%%%%% Ali Hany
%LPF of 1Hz
BW_lpf_5 = 1;
H_lpf_5 = zeros(size(f));
H_lpf_5 = abs(f)<1;

% Apply the filter in the frequency domain
Y_5_freq = X_num .* H_lpf_5;
%plot(f,Y_5_freq);
%xlabel('frequency');
%ylabel('Amplitude');
%title('x_f_bandlimited');
% Convert back to the time domain
y_5_time = ifft(fftshift(Y_5_freq)) * fs;

% Plot the input and output signals
figure;
plot(t, x, 'b', t, real(y_5_time), 'r');
xlabel('Time (s)');
ylabel('Amplitude');
title(['Input Signal and Output of LPF with BW = ', num2str(BW_lpf_5), ' Hz']);
legend('Input Signal', 'Output Signal');
grid on;

%%%%%%%%    Q6  %%%%%%%%%%

BW_lpf_6 = 0.3; % Bandwidth of the LPF
% Design the ideal LPF in the frequency domain
H_lpf_6 = zeros(size(f));
indices_passband_6 = find(abs(f) <= BW_lpf_6/2);
H_lpf_6(indices_passband_6) = 1;

% Apply the filter in the frequency domain
Y_6_freq = X_num .* H_lpf_6;

% Convert back to the time domain
y_6_time = ifft(fftshift(Y_6_freq)) * fs;

% Plot the input and output signals
figure;
plot(t, x, 'b', t, real(y_6_time), 'r');
xlabel('Time (s)');
ylabel('Amplitude');
title(['Input Signal and Output of LPF with BW = ', num2str(BW_lpf_6), ' Hz']);
legend('Input Signal', 'Output Signal');
grid on;

%%%%%%%%    Q7  %%%%%%%%%%

% Define m(t)
t_m =0:ts:((N-1)*ts); % Adjust time vector for m(t)
m1= zeros(size(t_m));
m2= zeros(size(t_m));
m1(t_m >= 0 & t_m <= 4) = cos(2 * pi * 0.5 * t_m(t_m >= 0 & t_m <= 4));
m=m1+m2;

%%%%%%%%    Q1 repeated for m(t) %%%%%%%%%%

figure;
plot(t_m,m);
xlim([0,5]);
xlabel('t');
ylabel('m(t)');
title('m(t) in Time Domain');
grid on;

%%%%%%%%    Q2&3 repeated for m(t) %%%%%%%%%%

M_num = fftshift(fft(m)) *ts; % Using N for consistent frequency vector

f_m = f; % Use the same frequency vector as before

figure;
plot(f_m,abs(M_num));
xlabel('Frequency (Hz)');
ylabel('|M(f)|');
title('FT of m(t) (Computed)');
grid on;

% Analytical FT of m(t):
% m(t) = cos(pi*t) for 0 < t < 4
% M(f) = integral from 0 to 4 of cos(pi*t) * exp(-j2*pi*f*t) dt
% Using Euler's formula: cos(pi*t) = 0.5 * (exp(j*pi*t) + exp(-j*pi*t))
% M(f) = 0.5 * integral from 0 to 4 of [exp(j*pi*t(1-2f)) + exp(-j*pi*t(1+2f))] dt
% M(f) = 0.5 * [ sin(pi*(1-2f)*4) / (pi*(1-2f)) + sin(pi*(1+2f)*4) / (pi*(1+2f)) ]
M_analytical = zeros(size(f_m));
M_analytical = 2 * (sinc(4*(f - 0.5)) .* exp(-1j*4*pi*(f - 0.5)) +sinc(4*(f + 0.5)) .* exp(-1j*4*pi*(f + 0.5)));

hold on;
plot(f_m, abs(M_analytical), 'r--','LineWidth', 2);
legend('FT computed', 'FT analytical');
hold off;

%%%%%%%%    Q4 repeated for m(t) %%%%%%%%%%
% Compute power spectrum of m(t)

M_pos=M_num(f>=0);
M_value=abs(M_pos);
M_squ=M_value .^2;
max_amplitude=max(M_squ);
max_amplitude_index = find(abs(M_squ - max_amplitude) < 1e-10, 1);
last_index=length(M_squ);
for j=max_amplitude_index:last_index
  if(abs(M_squ(j))<=0.05*max_amplitude)
    BW=f(j+5000);
    break;
  end
end
fprintf('5%% power bandwidth of m(t) ≈ %.4f Hz\n', BW);






%--------------------------------------------------------------------------------------------------
%%%%%%%%%                Q8                 %%%%%%%%%% Karim Walid
%% FDM:   x(t) -> LPF -> DSB-SC -> BPF
%%        m(t) -> LPF -> SSB -> BPF

fc1 = 20;
fc2 = 23;       % 20 + 1 (bw) + 2 (guard band) + 1 (bw)
c1_t = cos(2*pi*fc1*t);
c2_t = cos(2*pi*fc2*t);

%%%%%%%%    Q10      %%%%%%%
%%          TX          %%
%%        s(t) -> LPF -> DSB-SC -> BPF
%%        s(t) -> LPF -> SSB -> BPF
%% We already have x(t)_lpf from Q5: y_5_time
x_t_lpf=y_5_time;
s1_t = x_t_lpf .* c1_t;   %%%% DSB-SC modulation %%%%%

%%%%%%%%%%%%%          LPF for m(t)          %%%%%%%%%%%
% We modulate the signal by DSB-SC then apply a BPF for the Upper side band.

% Design the ideal LPF in the frequency domain


% The Expression for Single Side Band with Upper Side Band (USB)  %% Q9 %%:
m_t_dsbsc = m .* c2_t;
H = zeros(size(f));
H(f>fc2 & f<(fc2+1))=1;%% +ve frequency
H(f<-fc2 & f>-(fc2+1))=1;%% -ve frequency
M_f = fftshift(fft(m_t_dsbsc)) * ts;
%figure;
%plot(f,M_f);
%xlabel("Frequency");
%ylabel("M(f)_dsb");
%title('M(f)_dsb')

Musb = M_f .* H; %% SSB


%%                          FDM                         %%
%%%%%%%%%             BPF for s1_t             %%%%%%%%%%


Y5 = fftshift(fft(s1_t)) * ts;  % FT of s1_t

% BPF parameters
fc1 = 20;
BW_bpf = 2;

% Ideal BPF
H1_bpf = zeros(size(f));
H1_bpf((f >= fc1 - 1) & (f <= fc1 + 1)) = 1;
H1_bpf((f <= -fc1 + 1) & (f >= -fc1 - 1)) = 1;

% Apply BPF
S1_f = Y5 .* H1_bpf;

% FT of S1(f) after BPF
s1_t_bpf = ifft(ifftshift(S1_f)) * fs;



%%%%%%%%%%%%%%             BPF for s2_t                    %%%%%%%%%%%%
%S2f_temp = fftshift(fft(Musb)) * ts;  % FT of s2_t
% Band Pass Filter Parameters
%fc2 = 24;
% Ideal BPF
%H2_bpf = zeros(size(f));
%H2_bpf((f >= fc2 - 1) & (f <= fc2 + 1)) = 1;
%H2_bpf((f <= -fc2 + 1) & (f >= -fc2 - 1)) = 1;
% Apply BPF
%S2_f = S2f_temp .* H2_bpf;

% FT of S1(f) after BPF
s2_t_bpf = ifft(ifftshift(Musb)) * fs;

s_t = s1_t_bpf + s2_t_bpf;

figure;
plot(t,s_t);
xlabel("t(time)");
ylabel("s(t))");
title('s(t) in Time Domain');
legend('s(t)');
grid on;



%%%%%%%%%%%%%%%%  Q12 RX  %%%%%%%%%%%
%% RX same steps as TX but in reverse
%%      s(t) -> BPF -> DSB-SC (*cos(2pi*f1)) -> LPF
%%      s(t) -> BPF -> SSB  (*cos(2pi*f2)) -> LPF
S_f = fftshift(fft(s_t)) * ts;
%bpf
H2_bpf = zeros(size(f));
H2_bpf((f >= fc2 - 1) & (f <= fc2 + 1)) = 1;
H2_bpf((f <= -fc2 + 1) & (f >= -fc2 - 1)) = 1;

s11_f = S_f .* H1_bpf;
s22_f = S_f .* H2_bpf;

%demodulator
s11_t = ifft(ifftshift(s11_f)) * fs;
s22_t = ifft(ifftshift(s22_f)) * fs;

s11_t = s11_t .* cos(2*pi*fc1*t);
s22_t = s22_t .* cos(2*pi*fc2*t);

%lpf
s11_f = fftshift(fft(s11_t)) * ts;
s22_f = fftshift(fft(s22_t)) * ts;

H = zeros(size(f));
H = abs(f)<1;
s11_f = s11_f .* H;

H = zeros(size(f));
H = abs(f)<2;
s22_f = s22_f .* H;

s11_t = ifft(ifftshift(s11_f)) * fs;
s22_t = ifft(ifftshift(s22_f)) * fs;



%%Plotting
figure;
plot(t,x_t_lpf);
hold on
plot(t,s11_t,'r--');
xlabel("time");
legend('x(t)_(sent)','x(t)_(recieved)');
title('x(t)_(sent) & x(t)_(recieved)')


figure;
plot(t_m,m);
xlim([-2,6]);
hold on
plot(t_m,s22_t,'r--');
xlabel("time");
legend('m(t)_(sent)','m(t)_(recieved)');
title('m(t)_(sent) & m(t)_(recieved)');
grid on;
