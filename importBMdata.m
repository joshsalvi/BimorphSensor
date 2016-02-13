function importBMdata(datapath)
% This script imports data, without spacefiles, LabVIEW's
% variableloadclamp.vi to MATLAB. You will need to input the full directory
% name containing all files with your data, including the logfile.
%
% importVLCdata(datapath)
%
% datapath : string of the path containing your data
% MSVLC : Mech Stim or VLC? (1=VLC, 2=Mech Stim)
%
% Example:
% importVLCdata('/Users/username/Data/2014-09-30.01/Ear 1/Cell 2/')
%
% In some cases, the logfile will have incorrect indices. This becomes
% obvious for cases where "Fs" (the sampling rate) and "pulse" are not
% their expected values. If this occurs, modify the script and run it
% again.
%
% jsalvi@rockefeller.edu
%


% Find all appropriate files
file = dir(sprintf('%s%s',datapath,'*.txt'));   % find all data files
logfile = dir(sprintf('%s%s',datapath,'*.log'));    % find the logfile
a = length(file);       % number of sessions in the directory
%{
% Find all appropriate files
file = dir(sprintf('%s%s',datapath,'*VLC*.txt'));   % find all data files
logfile = dir(sprintf('%s%s',datapath,'*.csv'));    % find the logfile
a = length(file);       % number of sessions in the directory
%}
% Import logdata
logdata = importdata(sprintf('%s%s',datapath,logfile.name));  % logdata.data contains log data of interest
if isstruct(logdata) == 0
%    logdata.data = logdata;
end
if isempty('logdata.textdata(isnan(logdata.data(:,8))==0,3)')==0
    comments = logdata.textdata(isnan(logdata.data(:,8))==0,2); % import comments
end

Fs = logdata.data(1,1);


% Note that some logfiles will have formatting issues and the indices
% listed above may be off by ±10 in certain cases. If this is true, open
% logdata.data and modify the script accordingly.



%Number of traces
ntrace=logdata.data(isnan(logdata.data(:,8))==0,8);
numavg=logdata.data(isnan(logdata.data(:,8))==0,10);
% Find all raw files
for j = 1:a
    rawfiles(j) = isempty(findstr(file(j).name,'raw'));
end
nonraw = find(rawfiles==1);raw = find(rawfiles==0);
ntraceraw = ntrace(ntrace~=0);
numavgraw = numavg(numavg~=0);

% Import the data, some of this may be redundant
for i = 1:a
    data2{i}=importdata(sprintf('%s%s',datapath,file(i).name));   % initial import
end

for i =1:a
    Xd{i}=data2{i};     % extract data from its structure format (not necessary, but easier)
end
clear data2 data 
%}
% Import non-raw time traces



dt = 1/Fs(1);
sizeX = length(Xd);
for j = 1:sizeX
        tvec{j} = 0:dt:length(Xd{j})*dt-dt;
end

clear i j 

% Save extracted data
disp('Saving...');
savefile = sprintf('%s%s',datapath,'Extracted Data.mat');
save(savefile);
disp('Finished.');

end
