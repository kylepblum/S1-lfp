%% Load and process cds File
 
pathname = '~/LimbLab/Projects/S1-LFP/data';
% filename = 'Duncan_20191106_CObumpmove_cds.mat';
filename = 'Han_20191101_CObumpmove_cds.mat';

load([pathname filesep filename])

% Split into trials with lfp data
trials=table2struct(cds.trials);
lfp = table2array(cds.lfp);

%Get rid of non reward trials 
result = [trials.result]=='R';
trials = trials(result);
result = ~isnan([trials.goCueTime]);
trials = trials(result);

%Align trials to start/stop time
for n = 1:numel(trials)
    start = find(lfp(:,1) == trials(n).startTime);
    fin = find(lfp(:,1) == trials(n).endTime);
    trials(n).lfp = lfp(start:fin,2:numel(lfp(1,:))); 
end

%% Test filterLFP script
params = [];
trials = filterLFP(trials,params);
%% create movement onset from spike td data
%use movement onset from td with spike data
for n =1:numel(trials)
trials(n).idx_movement_on = (td(n).idx_movement_on*50)/2;
end
result = ~isnan([trials.idx_movement_on]);
trials = trials(result);

% align by movement onset
for n = 1:numel(trials)
    if trials(n).idx_movement_on < numel(trials(n).lfp(:,1))
    if trials(n).idx_movement_on + 1000 < numel(trials(n).lfp(:,1))
    start = trials(n).idx_movement_on - 450;
    fin = trials(n).idx_movement_on + 1000;
    trials(n).lfp = trials(n).lfp(start:fin,:);
    else
        trials(n) = [];
    end
    end
end
 
%% Apply 1 Hz High pass filter
[a,b] = butter(1, 1/1000, 'high');
for n = 1:numel(trials)
     trials(n).lfp = filtfilt(a,b,trials(n).lfp);
end

%% Split into act, comb, and pas
% Get bump-move trials
trials_temp = [];
for i = 1:numel(trials)
    if ~isnan(trials(i).idx_movement_on) && ~isnan(trials(i).bumpDir) && trials(i).bumpTime > trials(i).goCueTime
        if strcmp(trials(i).result,'R') || strcmp(trials(i).result,'F')
            if isempty(trials_temp)
                trials_temp = trials(i);
            else
                trials_temp = [trials_temp; trials(i)];
            end
        end
    end
end

trials_comb = trials_temp;

% Get bump-only trials
trials_temp = [];
for i = 1:numel(trials)
    if ~isnan(trials(i).idx_movement_on) && ~isnan(trials(i).bumpDir) && trials(i).bumpTime < trials(i).goCueTime
        if strcmp(trials(i).result,'R') || strcmp(trials(i).result,'F')
            if isempty(trials_temp)
                trials_temp = trials(i);
            else
                trials_temp = [trials_temp; trials(i)];
            end
        end
    end
end

trials_bump = trials_temp; %td bump

% Get non-bump trials
trials_temp = [];
for i = 1:numel(trials)
    if ~isnan(trials(i).idx_movement_on) && isnan(trials(i).bumpDir)
        if strcmp(trials(i).result,'R') && ~isnan(trials(i).goCueTime)
            if isempty(trials_temp)
                trials_temp = trials(i);
            else
                trials_temp = [trials_temp; trials(i)];
            end
        end
    end
end
trials_act = trials_temp; %td non-bump



td_bl = trials;
%% separate into dirs

% all donditions
dir = [trials.tgtDir]==0;
trials0 = trials(dir);

dir = [trials.tgtDir]==45;
trials45 = trials(dir);

dir = [trials.tgtDir]==90;
trials90 = trials(dir);

dir = [trials.tgtDir]==135;
trials135 = trials(dir);

dir = [trials.tgtDir]==180;
trials180 = trials(dir);

dir = [trials.tgtDir]==225;
trials225 = trials(dir);

dir = [trials.tgtDir]==270;
trials270 = trials(dir);

dir = [trials.tgtDir]==315;
trials315 = trials(dir);

%Active dirs
dir = [trials_act.tgtDir]==0;
atrials0 = trials_act(dir);

dir = [trials_act.tgtDir]==45;
atrials45 = trials_act(dir);

dir = [trials_act.tgtDir]==90;
atrials90 = trials_act(dir);

dir = [trials_act.tgtDir]==135;
atrials135 = trials_act(dir);

dir = [trials_act.tgtDir]==180;
atrials180 = trials_act(dir);

dir = [trials_act.tgtDir]==225;
atrials225 = trials_act(dir);

dir = [trials_act.tgtDir]==270;
atrials270 = trials_act(dir);

dir = [trials_act.tgtDir]==315;
atrials315 = trials_act(dir);

% center bump dirs
for n=1:numel(trials_bump)
    if trials_bump(n).bumpDir > 315
        trials_bump(n).bumpDir = trials_bump(n).bumpDir - 360;
    end
end

dir = [trials_bump.bumpDir]==0;
btrials0 = trials_bump(dir);

dir = [trials_bump.bumpDir]==45;
btrials45 = trials_bump(dir);

dir = [trials_bump.bumpDir]==90;
btrials90 = trials_bump(dir);

dir = [trials_bump.bumpDir]==135;
btrials135 = trials_bump(dir);

dir = [trials_bump.bumpDir]==180;
btrials180 = trials_bump(dir);

dir = [trials_bump.bumpDir]==225;
btrials225 = trials_bump(dir);

dir = [trials_bump.bumpDir]==270;
btrials270 = trials_bump(dir);

dir = [trials_bump.bumpDir]==315;
btrials315 = trials_bump(dir);
%% plot LFPs
figure;
for n= 1:96
%     subplot(12,8,n)
    plot(trials(78).lfp(:,n), 'b')
hold on
end
xline((114*10)/2, 'LineWidth', 2)
xlabel('time (0.5 ms bins)')
ylabel('µV') 
title('Plot of  LFP signals for all electrodes across one trial (bump during reach)')
%% Plot multiple trials for LFP in each direction with an example electrode
%Pick an example electrode
e=44;
randt = [14, 25, 36];
for n=1:numel(trials)
    tgtDirs(n) = trials(n).tgtDir;
end
tgtDirs = unique(tgtDirs);
A=find(isnan(tgtDirs));
tgtDirs(A) = [];

idx_subplot = [6,3, 2, 1, 4,7,8,9];

for t = 1:numel(randt)
    for s = 1:numel(idx_subplot)
        subplot(3,3,idx_subplot(s))
        plot((eval(strcat('trials',string(tgtDirs(s)),'(',string(randt(t)),')','.lfp(:,',string(e),')'))))
        hold on
        axis([0 1450 -250 250])
        xlabel('time (ms)')
        ylabel(strcat('µV in ', string(tgtDirs(s)), ' dir'))
        xline(450);
    end
end

%% Calculate mean LFP of one electrode across trials in each direction
sm = [];
meanplot = [];
e = 22;
for n=1:numel(trials)
    tgtDirs(n) = td(n).tgtDir;
end
tgtDirs = unique(tgtDirs);
A=find(isnan(tgtDirs));
tgtDirs(A) = [];


for s = 1:numel(tgtDirs)
    for t = 1:numel(eval(strcat('trials', string(tgtDirs(s)))))
        if ~isempty(eval(strcat('trials', string(tgtDirs(s)),'(', string(t),')', '.lfp')))
        sm(t,:,s) = eval(strcat('trials', string(tgtDirs(s)),'(', string(t),')', '.lfp(:,',string(e),')'));    
        end
    end
    meanplot(1,:,s) = sum(sm(:,:,s))/numel(eval(strcat('trials', string(tgtDirs(s)))));
end

idx_subplot = [6, 2, 4, 8,];
for s = 1:numel(idx_subplot)
        subplot(3,3,idx_subplot(s))
        plot(meanplot(:,:,s))
        axis([0 1500 -100 100])
        xlabel('time (0.5 ms bins)')
        ylabel(strcat( 'in ', string(tgtDirs(s)), ' dir'))
        xline(450)
end

%% FFT Practice
t = 200; % trial number
Fs = 2000; % sampling frequency
Ts= 1/Fs; %sampling period or time step

 y4 = trials(t).lfp(:,22); %for a single trial

% % for calculating mean signal over one trial
% clear sm
% 
% for n = 1:numel(trials(t).lfp(:,1))
%     sm(n,:) = sum(trials(t).lfp(n,:));
% end
%  y4=sm/numel(trials(t).lfp(1,:));
% 
nfft=length(y4); % length of time domain signal
nfft2=2^nextpow2(nfft); % length of signal in power of 2
ff=fft(y4,nfft2);
fff=ff(1:nfft2/2);
xfft = Fs*(0:nfft2/2-1)/nfft2;

subplot(2,1,1);
plot(y4);
xlabel('Time (0.5 ms bins)');
ylabel('Amplitude (µV)');
title('Time Domain Signal');
xline(450);


subplot(2,1,2);
plot(xfft,abs(fff));
xlabel('Frequency (Hz)');
ylabel('Normalized Amplitude');
axis([0 100 0 30000]);
title('Frequency Domain Signal');


%% FFT on different movement conditions each trial
%act case 
 e = 46; % example electrode
 Freq = [0 100]; % band of interest
 Pow = [0 50000];% power of interest
 lfp_act = [];
 for n= 1:numel(trials_act)
     lfp_act(:,n) = trials_act(n).lfp(:,e);
 end
x = lfp_act; 
 
Fs = 2000; % sampling frequency

for n = 1:numel(trials_act)
    x = trials_act(n).lfp(:,e);
    nfft=length(x); % length of time domain signal
    nfft2=2^nextpow2(nfft); % length of signal in power of 2
    ff=fft(x,nfft2);
    fff=ff(1:nfft2/2);
    xfft = Fs*(0:nfft2/2-1)/nfft2;
    
    subplot(3,1,1);
    plot(xfft,abs(fff));
    hold on
end
    xlabel('Frequency (Hz)');
    ylabel('Normalized Amplitude');
    title('Active Frequency Bands');
    axis([Freq(1) Freq(2) Pow(1) Pow(2)])


    
%pas case 
 lfp_bump = [];
 
 for n= 1:numel(trials_bump)
     lfp_bump(:,n) = trials_bump(n).lfp(:,e);
 end

x = lfp_bump; 
Fs = 2000; % sampling frequency 

for n = 1:numel(trials_bump)
    x = trials_bump(n).lfp(:,e);
    nfft=length(x); % length of time domain signal
    nfft2=2^nextpow2(nfft); % length of signal in power of 2
    ff=fft(x,nfft2);
    fff=ff(1:nfft2/2);
    xfft = Fs*(0:nfft2/2-1)/nfft2;
    
    subplot(3,1,2);
    plot(xfft,abs(fff));
    hold on
end

xlabel('Frequency (Hz)');
ylabel('Normalized Amplitude');
title('Center Hold Bump Frequency Bands');
axis([Freq(1) Freq(2) Pow(1) Pow(2)])



%combined movebump case 
lfp_comb = [];
 for n= 1:numel(trials_comb)
     lfp_comb(:,n) = trials_comb(n).lfp(:,e);
 end
 
x=lfp_comb;
Fs = 2000; % sampling frequency 

for n = 1:numel(trials_comb)
    x = trials_comb(n).lfp(:,e);
    nfft=length(x); % length of time domain signal
    nfft2=2^nextpow2(nfft); % length of signal in power of 2
    ff=fft(x,nfft2);
    fff=ff(1:nfft2/2);
    xfft = Fs*(0:nfft2/2-1)/nfft2;
    
    subplot(3,1,3);
    plot(xfft,abs(fff));
    hold on
end
xlabel('Frequency (Hz)');
ylabel('Normalized Amplitude');
title('Move-bump Frequency Bands');
    axis([Freq(1) Freq(2) Pow(1) Pow(2)])

%% mean FFT in different directions for an example electrode, act and pas conditions
close all;
e = 22; %example electrode
Freq = [0 100]; % band of interest
Pow = [0 3500];% power of interest
sm = []; meanplot = [];
for n=1:numel(trials)
    tgtDirs(n) = td(n).tgtDir;
end

tgtDirs = unique(tgtDirs);
A=find(isnan(tgtDirs));
tgtDirs(A) = [];

%act
for s = 1:numel(tgtDirs)
    for t = 1:numel(eval(strcat('atrials', string(tgtDirs(s)))))
        if ~isempty(eval(strcat('atrials', string(tgtDirs(s)),'(', string(t),')', '.lfp')))
            if isnan(eval(strcat('atrials', string(tgtDirs(s)),'(', string(t),').bumpTime')))
                sm(t,:,s) = eval(strcat('atrials', string(tgtDirs(s)),'(', string(t),')', '.lfp(:,',string(e),')'));
            end
        end
    end
     meanplot(1,:,s) = sum(sm(:,:,s))/numel(sm(:,1,s));
end

idx_subplot = [6, 2,  4,8];

for s = 1:numel(idx_subplot)
    x = meanplot(1,:,s);
    nfft=length(x); % length of time domain signal
    nfft2=2^nextpow2(nfft); % length of signal in power of 2
    ff=fft(x,nfft2);
    fff=ff(1:nfft2/2);
    xfft = Fs*(0:nfft2/2-1)/nfft2;
    
    
    subplot(3,3,idx_subplot(s))
    plot(xfft,abs(fff));
    hold on
    xlabel('Frequency (Hz)');
    ylabel('Norm Amp.')
    title(strcat( 'Frequency bands in ', string(tgtDirs(s)), ' dir'));
     axis([Freq(1) Freq(2) Pow(1) Pow(2)])
end
% Pas
sm = []; meanplot = [];
for s = 1:numel(tgtDirs)
    for t = 1:numel(eval(strcat('btrials', string(tgtDirs(s)))))
        if ~isempty(eval(strcat('btrials', string(tgtDirs(s)),'(', string(t),')', '.lfp')))
            if eval(strcat('btrials', string(tgtDirs(s)),'(', string(t),').ctrHoldBump')) == 1
                sm(t,:,s) = eval(strcat('btrials', string(tgtDirs(s)),'(', string(t),')', '.lfp(:,',string(e),')'));
            end
        end
    end
     meanplot(1,:,s) = sum(sm(:,:,s))/numel(sm(:,1,s));
end

idx_subplot = [6, 2,  4,8];

for s = 1:numel(idx_subplot)
    x = meanplot(1,:,s);
    nfft=length(x); % length of time domain signal
    nfft2=2^nextpow2(nfft); % length of signal in power of 2
    ff=fft(x,nfft2);
    fff=ff(1:nfft2/2);
    xfft = Fs*(0:nfft2/2-1)/nfft2;
    
    
    subplot(3,3,idx_subplot(s))
    plot(xfft,abs(fff),'r');
    hold on
    xlabel('Frequency (Hz)');
    ylabel('Norm Amp.')
    title(strcat( 'Frequency bands in ', string(tgtDirs(s)), ' dir'));
     axis([Freq(1) Freq(2) Pow(1) Pow(2)])
end
legend('act','pas')
%% FFT after mean signal subtracted from trial for electrode
close all;
e = 22; %example electrode
Freq = [0 1000]; % band of interest
Pow = [0 7000];% power of interest
sm = []; meanplot = [];
for n=1:numel(trials)
    tgtDirs(n) = td(n).tgtDir;
end

tgtDirs = unique(tgtDirs);
A=find(isnan(tgtDirs));
tgtDirs(A) = [];

% act mean calculation
for s = 1:numel(tgtDirs)
    for t = 1:numel(eval(strcat('atrials', string(tgtDirs(s)))))
        if ~isempty(eval(strcat('atrials', string(tgtDirs(s)),'(', string(t),')', '.lfp')))
            if isnan(eval(strcat('atrials', string(tgtDirs(s)),'(', string(t),').bumpTime')))
                sm(t,:,s) = eval(strcat('atrials', string(tgtDirs(s)),'(', string(t),')', '.lfp(:,',string(e),')'));
            end
        end
    end
     meanplot(1,:,s) = sum(sm(:,:,s))/numel(sm(:,1,s));
end

% act mean subtraction
sig = [];
for s = 1:numel(tgtDirs)
    for t = 1:numel(eval(strcat('atrials', string(tgtDirs(s)))))
        if ~isempty(eval(strcat('atrials', string(tgtDirs(s)),'(', string(t),')', '.lfp')))
            if isnan(eval(strcat('atrials', string(tgtDirs(s)),'(', string(t),').bumpTime')))
            sig(t,:,s) = sm(t,:,s) - meanplot(1,:,s);
            end
        end
    end
end

    

idx_subplot = [6,3, 2, 1, 4,7,8,9];

for s = 1:numel(idx_subplot)
    for t = 1:numel(sig(:,1,1))
    x = sig(t,:,s);
    nfft=length(x); % length of time domain signal
    nfft2=2^nextpow2(nfft); % length of signal in power of 2
    ff=fft(x,nfft2);
    fff=ff(1:nfft2/2);
    xfft = Fs*(0:nfft2/2-1)/nfft2;
    
    
    subplot(3,3,idx_subplot(s))
    plot(xfft,abs(fff), 'r');
    hold on
    xlabel('Frequency (Hz)');
    ylabel('Norm Amp.')
    title(strcat( 'Frequency bands in ', string(tgtDirs(s)), ' dir'));
     axis([Freq(1) Freq(2) Pow(1) Pow(2)])
    end
end
% Pas
sm = []; meanplot = [];
for s = 1:numel(tgtDirs)
    for t = 1:numel(eval(strcat('btrials', string(tgtDirs(s)))))
        if ~isempty(eval(strcat('btrials', string(tgtDirs(s)),'(', string(t),')', '.lfp')))
            if eval(strcat('btrials', string(tgtDirs(s)),'(', string(t),').ctrHoldBump')) == 1
                sm(t,:,s) = eval(strcat('btrials', string(tgtDirs(s)),'(', string(t),')', '.lfp(:,',string(e),')'));
            end
        end
    end
     meanplot(1,:,s) = sum(sm(:,:,s))/numel(sm(:,1,s));
end

% pas mean subtraction
sig = [];
for s = 1:numel(tgtDirs)
    for t = 1:numel(eval(strcat('btrials', string(tgtDirs(s)))))
        if ~isempty(eval(strcat('btrials', string(tgtDirs(s)),'(', string(t),')', '.lfp')))
            if eval(strcat('btrials', string(tgtDirs(s)),'(', string(t),').ctrHoldBump')) == 1
            sig(t,:,s) = sm(t,:,s) - meanplot(1,:,s);
            end
        end
    end
end


idx_subplot = [6,3, 2, 1, 4,7,8,9];

for s = 1:numel(idx_subplot)
    for t = 1:numel(sig(:,1,1))
    x = sig(t,:,s);
    nfft=length(x); % length of time domain signal
    nfft2=2^nextpow2(nfft); % length of signal in power of 2
    ff=fft(x,nfft2);
    fff=ff(1:nfft2/2);
    xfft = Fs*(0:nfft2/2-1)/nfft2;
    
    
    subplot(3,3,idx_subplot(s))
    plot(xfft,abs(fff),'b');
    hold on
    xlabel('Frequency (Hz)');
    ylabel('Norm Amp.')
    title(strcat( 'Frequency bands in ', string(tgtDirs(s)), ' dir'));
    axis([Freq(1) Freq(2) Pow(1) Pow(2)])
    end
end
legend('act','pas')

%% Subtract mean LFP signal for a trial from an electrode

%insert trial number you wish to look at
t=82;
%modulated electrodes
e=[1, 3, 8, 9, 10, 11, 13, 15, 16];

%Calculate mean signal for that trial
clear sm
for n = 1:numel(trials(t).lfp(:,1))
    sm(n,:) = sum(trials(t).lfp(n,:));
end
sm=sm/numel(trials(t).lfp(1,:));

%plot electrode signal minus mean signal for that trial
for n = 1:numel(e)
    subplot(3,3,n)
    plot(trials(t).lfp(:,e(n))-sm)
    axis([0 numel(sm) -200 200])
    xlabel('time (0.5 ms)')
    ylabel('LFP µV') 
end


%% Calculate mean LFP signal

for n = 1:numel(trials)
    sig(:,:,n) = trials(n).lfp;
end
for t = 1:numel(sig(1,:,1))
    for n = 1:numel(sig(:,1,1))
        sm(n,t) = sum(sig(n,t,:));
        err(n,t)=std(sig(n,t,:));
    end
end
err= err/sqrt(numel(trials));
sm = sm/numel(trials);
%% plot mean LFP
plot(sm(:,9),'b')
xlabel('time (0.5 ms)')
ylabel('LFP µV') 
title('Plot of mean LFP signal for an electrode')
figure;
errorbar(sm(:,9), err(:,9), 'b')
xlabel('time (0.5 ms)')
ylabel('LFP µV') 
title('Plot of Error Bars LFP signal for an electrode (std/sqrt(N))')



