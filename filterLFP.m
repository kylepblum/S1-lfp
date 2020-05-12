%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function td_s = processLFP(trial_data,params)
%
%   Processes LFP by high pass filtering, rectifying, then low pass
%   filtering. Uses a Butterworth filter for both high and low pass
%   filtering. 
%
% INPUTS:
%   trial_data : the struct
%   params     : parameter structf
%       .highPass           :       frequency you want to high pass filter your LFP
%                                   (default is 1 Hz)
%       .highPassOrder      :       order of Butterworth filter you want to use (1)
%       .lowPass            :       frequency you want to low pass filter your LFP
%                                   (default is 500 Hz)
%       .lowPassOrder       :       order of Butterworth filter you want to use (1)
%       .idx_Signal          :       (string input) field where signal is in
%                                   struct (default is 'lfp')
%       .Fs                 :       sampling frequency (default is 2000 Hz)
%
% OUTPUTS:
%   td_s : the struct with processed LFP signal
%
% Written by Nathan Schimpf. Updated May 2020.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function td_s = filterLFP(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETER DEFAULTS
highPass        =  1;
highPassOrder   =  1;
lowPass         =  500;
lowPassOrder    =  1;
idxSignal      = 'lfp';
Fs              = 2000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(params), assignParams(who,params); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

td_s = trial_data;

% High pass filter
[a,b] = butter(highPassOrder, highPass/(Fs/2), 'high');
for trials = 1:numel(td_s)
    td_s(trials).(idxSignal) = filtfilt(a,b,td_s(trials).(idxSignal));
end

% Rectify
for trials = 1:numel(td_s)
    td_s(trials).(idxSignal) = abs(td_s(trials).(idxSignal));
end

% Low pass filter
[a,b] = butter(lowPassOrder, lowPass/(Fs/2), 'low');
for trials = 1:numel(td_s)
    td_s(trials).(idxSignal) = filtfilt(a,b,td_s(trials).(idxSignal));
end


     
