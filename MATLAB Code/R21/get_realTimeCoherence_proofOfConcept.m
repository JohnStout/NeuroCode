%% get_realTimeCoherence
clear; clc

%% initialize

% define LFPs to use
LFP1name = 'CSC1';
LFP2name = 'CSC9';

% setup
[srate,timing] = realTimeDetect_setup(LFP1name,LFP2name,0.25);

%% coherence detection

% for multitapers
params.tapers = [3 5];
params.Fs     = srate;
params.fpass  = [4 12];

% define a looping time
loop_time = .5; % minutes - note that this isn't perfect, but its a few seconds behind dependending on the length you set. The lag time changes incrementally because there is a 10-20ms processing time that adds up

% define amount of data to collect
amountOfData = .25; % seconds

% define number of samples that correspond to the amount of data in time
numSamples2use = amountOfData*srate;

% define for loop
looper = (loop_time*60)/amountOfData; % N minutes * 60sec/1min * (1 loop is about .250 ms of data)

%% initial
% start with this variable set to zero, it tells the code whether to clear
% the stream
clearIt = 0;

%% now add/remove time and calculate coherence
    
% for loop start
coh_temp = [];
coh = [];
for i = 1:looper
    
    tic;
    if i == 1
        % clear the stream
        clearStream(LFP1name,LFP2name);

        % get data
        pause(2); % grab one sec of data if this is the first loop
        try
            [succeeded, dataArray, timeStampArray, channelNumberArray, samplingFreqArray, ...
            numValidSamplesArray, numRecordsReturned, numRecordsDropped , funDur.getData ] = NlxGetNewCSCData_2signals(LFP1name, LFP2name);  
        catch
        end  
        
        % calculate og coherence
        coh_temp = [];
        coh_temp = coherencyc(dataArray(1,:),dataArray(2,:),params);
        coh(i)   = nanmean(coh_temp);
    else
        
    end

    % sometimes when collecting data, it will fail. But eventually, it will
    % collect the data successfully. It's not clear if the instance of
    % failing to collect the data results in lost data, or if its collected
    % later
    %{
    pause(0.25);
    next = 0;
    while next == 0
        try
            [~, dataArray_new, timeStampArray, ~, ~, ...
            numValidSamplesArray, numRecordsReturned, numRecordsDropped , funDur.getData ] = NlxGetNewCSCData_2signals(LFP1name, LFP2name);  
            next = 1;
        catch
            [~, dataArray_new, timeStampArray, ~, ~, ...
            numValidSamplesArray, numRecordsReturned, numRecordsDropped , funDur.getData ] = NlxGetNewCSCData_2signals(LFP1name, LFP2name);  
            next = 1;
        end  
    end
    %}
    
    % running this a ton of times until we get failures           
    [~, dataArray_new, timeStampArray, ~, ~, ...
    numValidSamplesArray, numRecordsReturned, numRecordsDropped , funDur.getData ] = NlxGetNewCSCData_2signals(LFP1name, LFP2name);  

    % define number of samples for continuous updating
    numSamples2use = [];
    numSamples2use = size(dataArray_new,2);
    
    % remove number of samples
    dataArray(:,1:numSamples2use)=[];
    
    % add data to end of dataArray
    dataArray = horzcat(dataArray,dataArray_new);

    % calculate coherence - chronux toolbox is way faster. Like sub 0.01
    % seconds sometimes, while wcoherence is around 0.05 sec.
    %[coh(i+1),phase,~,~,~,freq] = coherencyc(dataArray(1,:),dataArray(2,:),params);
    coh_temp = [];
    coh_temp = coherencyc(dataArray(1,:),dataArray(2,:),params);
    coh(i+1) = nanmean(coh_temp);
    
    % amount of data in consideration
    timings(i) = length(dataArray)/srate;
    
    % take averages
    outtoc(i) = toc;
    
    disp(['loop # ',num2str(i),'/',num2str(looper)])
       
end 

figure(); stem(coh)

% remove the very first 2 events, they will be excluded anyway. They take
% too long bc of initializing stuff
outtoc(1:2)=[];
outtoc_ms = outtoc*1000; % to ms
figure('color','w');
histogram(outtoc_ms);
ylabel('Number of streaming events')
xlabel('Time (ms) to stream and calculate coherence')
box off
title(['Events streamed at a rate of ', num2str(amountOfData),' seconds'])

% remove path of github download directory
rmpath(github_download_directory);

