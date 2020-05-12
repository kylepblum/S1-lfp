%% load & process TD file

pathname = '~/LimbLab/Projects/S1-gamma/data';
% filename = 'Duncan_20190524_CObumpmove_10ms.mat';
% filename = 'Duncan_20190710_CObumpmove_10ms.mat';
% filename = 'Duncan_20191016_CObumpmove_10ms.mat';
% filename = 'Duncan_20191106_CObumpmove_10ms.mat';
% filename = 'Han_20191010_CObumpmove_10ms.mat';
filename = 'Han_20191101_CObumpmove_10ms.mat';
% filename = 'Han_20191031_CObumpmove_10ms.mat';
 

load([pathname filesep filename]);


% Smooth spikes, get split TD, get movement onsets, get rid of unsorted units

%Smooth spikes
% smoothParams.signals = {'S1_spikes'};
% smoothParams.width = 0.03;
% smoothParams.calc_rate = true;
 trial_data.S1_spikes_bins = trial_data.S1_spikes;
% td = smoothSignals(trial_data,smoothParams);
td=trial_data;

%Get speed
td.speed = sqrt(td.vel(:,1).^2 + td.vel(:,2).^2);

%Get rectified velocity
td.vel_rect = [td.vel(:,1) td.vel(:,2) -td.vel(:,1) -td.vel(:,2)];
td.vel_rect(td.vel_rect < 0) = 0;

%Get accel
td.acc = diff(td.vel)./td.bin_size;
td.acc(end+1,:) = td.acc(end,:);

%Remove offset
td.pos(:,1) = td.pos(:,1)+0;
td.pos(:,2) = td.pos(:,2)+32;

% Smooth kinematic variables
smoothParams.signals = {'pos','vel','acc','force'};
smoothParams.width = 0.03;
smoothParams.calc_rate = false;
td = smoothSignals(td,smoothParams);

%Get rid of unsorted units
sorted_idx = find(td.S1_unit_guide(:,2)~=0);
td.S1_spikes = td.S1_spikes(:,sorted_idx);
td.S1_spikes_bins = td.S1_spikes_bins(:,sorted_idx);
td.S1_unit_guide = td.S1_unit_guide(sorted_idx,:);

%Split TD
splitParams.split_idx_name = 'idx_startTime';
splitParams.linked_fields = {'result','bumpDir','tgtDir'};
td = splitTD(td,splitParams);

%Get movement onset
td(isnan([td.idx_goCueTime])) = [];

moveOnsetParams.start_idx = 'idx_goCueTime';
moveOnsetParams.end_idx = 'idx_endTime';
td = getMoveOnsetAndPeak(td,moveOnsetParams);

% Separate TD into bump, act, and comb structures

% Get bump-move trials
td_temp = [];
for i = 1:numel(td)
    if ~isnan(td(i).idx_movement_on) && ~isnan(td(i).bumpDir) && td(i).idx_bumpTime > td(i).idx_goCueTime
        if strcmp(td(i).result,'R') || strcmp(td(i).result,'F')
            if isempty(td_temp)
                td_temp = td(i);
            else
                td_temp = [td_temp; td(i)];
            end
        end
    end
end

td_comb = td_temp;

% Get bump-only trials
td_temp = [];
for i = 1:numel(td)
    if ~isnan(td(i).idx_movement_on) && ~isnan(td(i).bumpDir) && td(i).idx_bumpTime < td(i).idx_goCueTime
        if strcmp(td(i).result,'R') || strcmp(td(i).result,'F')
            if isempty(td_temp)
                td_temp = td(i);
            else
                td_temp = [td_temp; td(i)];
            end
        end
    end
end

td_bump = td_temp; %td bump

% Get non-bump trials
td_temp = [];
for i = 1:numel(td)
    if ~isnan(td(i).idx_movement_on) && isnan(td(i).bumpDir)
        if strcmp(td(i).result,'R') && ~isnan(td(i).idx_goCueTime)
            if isempty(td_temp)
                td_temp = td(i);
            else
                td_temp = [td_temp; td(i)];
            end
        end
    end
end
td_act = td_temp; %td non-bump



td_bl = td;

% rebin to 50 ms size
td = binTD(td, 5);

% Get rid of non reward trials 
result = [td.result]=='R';
td = td(result);

%% Align movement onset
%use movement onset from td with spike data

result = ~isnan([td.idx_movement_on]);
td = td(result);

% align by movement onset
for n = 1:numel(td)
    if td(n).idx_movement_on < numel(td(n).S1_spikes(:,1))
    if td(n).idx_movement_on + 10 < numel(td(n).S1_spikes(:,1))
    start = td(n).idx_movement_on - 4;
    fin = td(n).idx_movement_on + 10;
    td(n).S1_spikes = td(n).S1_spikes(start:fin,:);
    else
        td(n).S1_spikes = [];
    end
    end
end

%% plot LFPs
figure;
for n= 1:96
%     subplot(12,8,n)
    plot(td(78).lfp(:,n), 'b')
hold on
end
xline((114*10)/2, 'LineWidth', 2)
xlabel('time (0.5 ms bins)')
ylabel('µV') 
title('Plot of  LFP signals for all electrodes across one trial (bump during reach)')
%% separate into dirs
dir = [td.tgtDir]==0;
td0 = td(dir);

dir = [td.tgtDir]==45;
td45 = td(dir);

dir = [td.tgtDir]==90;
td90 = td(dir);

dir = [td.tgtDir]==135;
td135 = td(dir);

dir = [td.tgtDir]==180;
td180 = td(dir);

dir = [td.tgtDir]==225;
td225 = td(dir);

dir = [td.tgtDir]==270;
td270 = td(dir);

dir = [td.tgtDir]==315;
td315 = td(dir);
%% Plot multiple trials for LFP in each direction with an example electrode
%Pick an example electrode
e=35;
randt = [14, 25, 36];
for n=1:numel(td)
    tgtDirs(n) = td(n).tgtDir;
end
tgtDirs = unique(tgtDirs);
A=find(isnan(tgtDirs));
tgtDirs(A) = [];

idx_subplot = [6,3, 2, 1, 4,7,8,9];

for t = 1:numel(randt)
    for s = 1:numel(idx_subplot)
        subplot(3,3,idx_subplot(s))
        plot((eval(strcat('trials',string(tgtDirs(s)),'(',string(randt(t)),')','.S1_spikes(:,',string(e),')'))))
        hold on
%         axis([0 15 0 2])
        xlabel('time (50 ms bins)')
        ylabel(strcat( 'in ', string(tgtDirs(s)), ' dir'))
    end
end
%% Calculate mean firing rate of one electrode across trials in each direction
sm = [];
meanplot = [];
e = 6;
for n=1:numel(td)
    tgtDirs(n) = td(n).tgtDir;
end
tgtDirs = unique(tgtDirs);
A=find(isnan(tgtDirs));
tgtDirs(A) = [];


for s = 1:numel(tgtDirs)
    for t = 1:numel(eval(strcat('td', string(tgtDirs(s)))))
        if ~isempty(eval(strcat('td', string(tgtDirs(s)),'(', string(t),')', '.S1_spikes')))
        sm(t,:,s) = eval(strcat('td', string(tgtDirs(s)),'(', string(t),')', '.S1_spikes(:,',string(e),')'));    
        end
    end
    meanplot(1,:,s) = sum(sm(:,:,s))/numel(eval(strcat('td', string(tgtDirs(s)))));
end

idx_subplot = [6,3, 2, 1, 4,7,8,9];
for s = 1:numel(idx_subplot)
        subplot(3,3,idx_subplot(s))
        plot(meanplot(:,:,s))
        axis([0 15 0 10])
        xlabel('time (50 ms bins)')
        ylabel(strcat( 'in ', string(tgtDirs(s)), ' dir'))
        xline(4)
end


%% Subtract mean signal for a trial from an electrode

%insert trial number you wish to look at
t=82;
%modulated electrodes
e=[1, 3, 8, 9, 10, 11, 13, 15, 16];

%Calculate mean signal for that trial
clear sm
for n = 1:numel(td(t).lfp(:,1))
    sm(n,:) = sum(td(t).lfp(n,:));
end
sm=sm/numel(td(t).lfp(1,:));

%plot electrode signal minus mean signal for that trial
for n = 1:numel(e)
    subplot(3,3,n)
    plot(td(t).lfp(:,e(n))-sm)
    axis([0 numel(sm) -200 200])
    xlabel('time (0.5 ms)')
    ylabel('LFP µV') 
end

%% Calculate mean LFP signal

for n = 1:numel(td)
    sig(:,:,n) = td(n).lfp;
end
for t = 1:numel(sig(1,:,1))
    for n = 1:numel(sig(:,1,1))
        sm(n,t) = sum(sig(n,t,:));
        err(n,t)=std(sig(n,t,:));
    end
end
err= err/sqrt(numel(td));
sm = sm/numel(td);
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



