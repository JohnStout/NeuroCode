load('data_122minutes')

figure('color','w');
stem(coh)
axis tight;
xlim([100 200])

% load in data_20min to keep working

%% coherence magnitudes 
coh(1) = [];

% plot the histogram of coherence magnitudes
figure('color','w')
h_og = histogram(coh,'FaceColor',[.6 .6 .6]); hold on;
h_og.FaceColor = [.6 .6 .6];
% zscore coherence
zscoredCoh = zscore(coh);
% remove first calculation
% get 1 std
cohHighStd = coh(zscoredCoh > 1);
cohHighVal = coh(dsearchn(zscoredCoh',1));
% plot
h_highS = histogram(cohHighStd,'FaceColor',[0 .6 0]);
h_highS.BinWidth = h_og.BinWidth;
% get -1std
cohLowStd = coh(zscoredCoh < -1);
cohLowVal = coh(dsearchn(zscoredCoh',-1));
% plot
h_lowS = histogram(cohLowStd,'FaceColor','r');
h_lowS.BinWidth = h_og.BinWidth;
box off;
legend(['All data'],['High Thresh. (1std = ',num2str(cohHighVal),')'],['Low Thresh. (-1std = ',num2str(cohLowVal),')'])
ylabel('Iteration')
xlabel('Coherence magnitude')

% chance of observing high/low coh
figure('color','w'); hold on;
prob_high = numel(cohHighStd)/numel(coh);
prob_low  = numel(cohLowStd)/numel(coh);
bar(1,prob_high,'FaceColor',[0 .6 0]);
bar(2,prob_low,'FaceColor','r');
ylabel('p(coherence)')
axes = gca;
axes.XTick = [ 1 2 ];
axes.XTickLabel = [{'High Coherence'} {'Low Coherence'}]
axes.XTickLabelRotation = 45;

%% coherence durations

% each sample is sampled at .25 seconds, therefore, if the size is 1, it
% was sampled at .25, if size is 2, then 0.5, etc..

% get sizes of timeStamps - to use cellfun2, run Startup in pipeline
timeSizes = cell2mat(cellfun2(timeStamps,'size',{'2'}));

% max of 
[maxTime, maxLoc] = max(timeSizes);

% convert to seconds
timeConv = zeros(size(timeSizes));
prob_observe = [];
for i = 1:maxTime
    timeConv(timeSizes == i) = i*amountOfData;
    
    % get p(i)
    prob_observe(i) = numel(find(timeSizes == i))/numel(timeConv);
end

% histogram of sampled times
figure('color','w')
plot(prob_observe,'k','LineWidth',2)
box off
ylabel('Coherence Magnitude')
xlabel('Duration Possibilities (ms)')
axes = gca;
axes.XTick = [1:maxTime];
for i = 1:maxTime
    axes.XTickLabel{i} = i*amountOfData;
end
axes.XTickLabelRotation = 45;

    
%% script function

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

figure('color','w')