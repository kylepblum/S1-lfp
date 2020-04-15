%% Load and process cds File
 
pathname = '~/LimbLab/Projects/S1-LFP/data';
filename = 'Duncan_20191106_CObumpmove_cds.mat';

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
    trials(n).lfp = lfp(start:fin,:); 
end
%% create movement onset from spike td data
%use movement onset from td with spike data
for n =1:numel(trials)
trials(n).idx_movement_on = (td(n).idx_movement_on*10)/2;
end
result = ~isnan([trials.idx_movement_on]);
trials = trials(result);

% align by movement onset
for n = 1:numel(trials)
    start = find(trials(n).lfp(:,1) == trials(n).lfp((trials(n).idx_movement_on - 450),1));
    fin = find(trials(n).lfp(:,1) == trials(n).lfp((trials(n).idx_movement_on + 1000),1));
    trials(n).lfp = lfp(start:fin,2:numel(lfp(1,:))); 
end
%% Apply 1000 Hz High pass filter
[a,b] = butter(1, 1/15000, 'high');
for n = 1:numel(trials)
     trials(n).lfp = filtfilt(a,b,trials(n).lfp);
end

%% plot LFPs
figure;
for n= 1:96
%     subplot(12,8,n)
    plot(trials(78).lfp(:,n), 'b')
hold on
end
xline((114*10)/2, 'LineWidth', 2)
xlabel('time (0.5 ms bins)')
ylabel('�V') 
title('Plot of  LFP signals for all electrodes across one trial (bump during reach)')

%% separate into dirs
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
%% Plot multiple trials for LFP in each direction with an example electrode
%Pick an example electrode
e=1;
randt = [1, 2, 3];
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
        plot(eval(strcat('trials',string(tgtDirs(s)),'(',string(randt(t)),')','.lfp(:,',string(e),')')))
        hold on
        axis([0 1500 -250 250])
        xlabel('time (ms)')
        ylabel(strcat('�V in ', string(tgtDirs(s)), ' dir'))
        xline(225);
    end
end


    


%% Subtract mean signal for a trial from an electrode

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
    ylabel('LFP �V') 
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
ylabel('LFP �V') 
title('Plot of mean LFP signal for an electrode')
figure;
errorbar(sm(:,9), err(:,9), 'b')
xlabel('time (0.5 ms)')
ylabel('LFP �V') 
title('Plot of Error Bars LFP signal for an electrode (std/sqrt(N))')



