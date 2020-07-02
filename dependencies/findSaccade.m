function [sacTime, sacIx]= findSaccade(time,trace)

%% parameters to tweak
pos_filter_length=5; % short filter to smooth for noisy eye traces
trace=filtfilt(ones(pos_filter_length,1)/pos_filter_length, 1, trace')';
filter_length=20; %number of amples to average for baseline velocity
detect_thresh=300; % 300 leo, 800 nancy, 1000 nancy at 2khz threshold in deg/s to detect a saccade
start_thresh=5; % threshold in deg/s to determine the start and end of a saccade
min_isi=50; % 50 for 1000Hz SR....minimum number of samples between any two saccades
minDur=5; %minimum duration of a saccade
blink_isi=50; % ignore saccades this many samples around a blink
reproject_onto_saccade_direction=true; %work on signed velocity instead of speed
%minBin =11; %minimum/earliest possible bin allowed to have a sacccade
%% calculate veocities
tps=median(diff(time));
vel=[[0;0],(trace(:,3:end)-trace(:,1:end-2))/(tps*2),[0;0]];
speed=sqrt(sum(vel.^2));

a = 1;
b = ones(1,filter_length/2)/filter_length/2;
speedf = filter(b,a,speed);


speedspeedf=speed-speedf;

potential_saccades=diff([0 (speedspeedf)>detect_thresh])==1;%check data start condition, consider loader more data before...
potential_saccades=find(potential_saccades);
if numel(potential_saccades)>1
    %potential_saccades(potential_saccades < minBin) = [];
    potential_saccades([min_isi diff(potential_saccades)]<min_isi)=[];
    
    [~, maxix]     = max(speedspeedf); 
    potential_saccades(potential_saccades > maxix) = [];
    [~, closestIx] = min(abs(potential_saccades - maxix));
    sacIx = potential_saccades(closestIx);
else
    sacIx   = potential_saccades;
end

sacTime = time(sacIx);

assert(~isempty(sacIx), 'failed to identify a saccade!')