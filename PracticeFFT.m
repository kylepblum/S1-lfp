Fs = 1000; % sampling frequency
Ts= 1/Fs; %sampling period or time step
dt = 0:Ts:5-Ts; %signal duration
f1=10;
f2=30;
f3=70;

% y=Asin(2pift+theta);

y1=10*sin(2*pi*f1*dt);
y2=5*sin(2*pi*f2*dt);
y3=1*sin(2*pi*f3*dt);
y4=y1+y2+y3;

% subplot(4,1,1)
% plot(dt,y1,'r')
% subplot(4,1,2)
% plot(dt,y2,'r')
% subplot(4,1,3)
% plot(dt,y3,'r')
% subplot(4,1,4)
% plot(dt,y4,'r')

nfft=length(y4); % length of time domain signal
nfft2=2^nextpow2(nfft); % length of signal in power of 2
ff=fft(y4,nfft2);
fff=ff(1:nfft2/2);
xfft = Fs*(0:nfft2/2-1)/nfft2;

subplot(2,1,1)
plot(dt,y4)
xlabel('Time (s)')
ylabel('Amplitude (V)')
title('Time Domain Signal')

subplot(2,1,2)
plot(xfft,abs(fff));
xlabel('Frequency (Hz)')
ylabel('Normalized Amplitude')
title('Frequency Domain Signal')

