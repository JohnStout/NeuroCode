%% number of VTEs before/after muscimol
clear; clc;
place2store = 'X:\01.Experiments\John n Andrew\Data and Exploratory\VTEs and Re suppression';
addpath(place2store)

% main Datafolder
Datafolders   = 'X:\01.Experiments\RERh Inactivation Recording';

% int name and vt name
int_name     = 'Int_JS_fixed';
vt_name      = 'VT1.mat';
missing_data = 'interp'; % interpolate missing data

% preSmooth smooths the position data by 1 second
preSmooth = 1; % 1 for yes, 0 for no

% define middle of stem
middleStemPosition = 500;

% what orientation is the stem?
stemOrientation = 'X'; % going forward in the x direction or y direction

% rat names - this was confirmed using 'Behavior' excel sheet in henrys
% folder. The 'Behavior' excel sheet was used to make fig8B in JNeuro paper
% LEFT OFF: need to make 'Int_VTE_JS.mat' files in Morty folders
rat{1} = 'Usher';
rat{2} = 'Rick';
rat{3} = 'Ratticus';
rat{4} = 'Ratdle';
rat{5} = 'Eric2';
rat{6} = '14-22';
rat{7} = 'Morty';

% condition
Day_cond{1} = 'Saline';
Day_cond{2} = 'Muscimol';
Day_cond{3} = 'Baseline';

% infusion treatment
infusion{1} = '\Baseline';
infusion{2} = '\Saline';
infusion{3} = '\Muscimol'; % within day infusion

% loop across rats
for rat_i = 1:length(rat)
    
    % datafolder for the rat
    datafolder = [];
    datafolder_rat = [Datafolders,'\',rat{rat_i}];
    
    for day_cond_i = 1:length(Day_cond)
        
        % datafolder that includes condition
        datafolder_cond = [datafolder_rat,'\',Day_cond{day_cond_i}];
        
        % get content in directory
        dir_content = dir(datafolder_cond);
        dir_content = extractfield(dir_content,'name');
 
        % find non folders
        remIdx = contains(dir_content,'.mat') | contains(dir_content,'.');
        
        % remove .mat files
        dir_content(remIdx)=[];
        
        % loop across infusion conditions
        for inf_i = 1:length(dir_content)
        
            % define datafolder
            datafolder = [datafolder_cond,'\',dir_content{inf_i}];

            % the struct array below that organizes the data cannot create
            % fields with space characters, therefore rename them if they
            % exist.
            if contains(dir_content{inf_i},'Baseline 1')
                dir_content{inf_i} = 'Baseline1';
            elseif contains(dir_content{inf_i},'Baseline 2')
                dir_content{inf_i} = 'Baseline2';
            end
            
            % if rat has a '-' character in it, remove
            if contains(rat{rat_i},'14-22')
                rat{rat_i} = 'rat1422';
            end
            
            % get VTE data in a structure array that is organized based on
            % rat, day condition, and infusion type
            [IdPhi.(rat{rat_i}).(Day_cond{day_cond_i}).(dir_content{inf_i}),...
             xPos.(rat{rat_i}).(Day_cond{day_cond_i}).(dir_content{inf_i}),...
             yPos.(rat{rat_i}).(Day_cond{day_cond_i}).(dir_content{inf_i}),...
             tsPos.(rat{rat_i}).(Day_cond{day_cond_i}).(dir_content{inf_i}),...
             xPosOG.(rat{rat_i}).(Day_cond{day_cond_i}).(dir_content{inf_i}),...
             yPosOG.(rat{rat_i}).(Day_cond{day_cond_i}).(dir_content{inf_i}),...
             tsPosOG.(rat{rat_i}).(Day_cond{day_cond_i}).(dir_content{inf_i}),...
             timeSpent.(rat{rat_i}).(Day_cond{day_cond_i}).(dir_content{inf_i})] ...
             = get_session_IdPhi(datafolder,int_name,vt_name,missing_data,middleStemPosition,stemOrientation,preSmooth);
                
            % restore appropriate naming for 14-22
            if contains(rat{rat_i},'rat1422')
                rat{rat_i} = '14-22';
            end
                    
        end       
    end
end
   
% -- z score data across sessions, but within each rat -- %

% rat names
rat{1} = 'Usher';
rat{2} = 'Rick';
rat{3} = 'Ratticus';
rat{4} = 'Ratdle';
rat{5} = 'Eric2';
rat{6} = 'rat1422';
rat{7} = 'Morty';

for rati = 1:length(rat)
    
    % col 1 is baseline, col 2 is infusion
    IdPhi_bl_saline{rati} = IdPhi.(rat{rati}).Saline.Baseline;
    IdPhi_sa_saline{rati} = IdPhi.(rat{rati}).Saline.Saline;
    
    % baseline
    IdPhi_bl_baseline1{rati} = IdPhi.(rat{rati}).Baseline.Baseline1;
    IdPhi_bl_baseline2{rati} = IdPhi.(rat{rati}).Baseline.Baseline2;
  
    % muscimol
    IdPhi_mu_baseline{rati} = IdPhi.(rat{rati}).Muscimol.Baseline;
    IdPhi_mu_muscimol{rati} = IdPhi.(rat{rati}).Muscimol.Muscimol;
    
    % concatenate data
    IdPhi_data{rati} = horzcat(IdPhi_bl_saline{rati},IdPhi_sa_saline{rati},IdPhi_bl_baseline1{rati},...
        IdPhi_bl_baseline2{rati},IdPhi_mu_baseline{rati},IdPhi_mu_muscimol{rati});
    
    % array required to return 
    return_sizes{rati} = horzcat(numel(IdPhi_bl_saline{rati}),numel(IdPhi_sa_saline{rati}),...
        numel(IdPhi_bl_baseline1{rati}),numel(IdPhi_bl_baseline2{rati}),...
        numel(IdPhi_mu_baseline{rati}),numel(IdPhi_mu_muscimol{rati}));
    
end

% zscore across sessions, then return data
return_sizes2 = cellfun(@numel,IdPhi_data,'UniformOutput',false);
IdPhi_cat  = horzcat(IdPhi_data{:});
%zIdPhi_cat = zscore(IdPhi_cat);
zIdPhi_cat = zscore(log(IdPhi_cat));
zIdPhi_dist = zIdPhi_cat; % save the distribution
for rati = 1:length(rat)
    zIdPhi_all{rati} = zIdPhi_cat(1:return_sizes2{rati}(1));
    zIdPhi_cat(1:return_sizes2{rati}(1))=[]; % erase this data so that we can continue extracting
    return_sizes2{rati}(1) = [];
end

for rati = 1:length(rat)
    
    % reorganize data, delete zIdPhi_all data and return sizes after
    % finishing with it
    zIdPhi.(rat{rati}).Saline.Baseline = zIdPhi_all{rati}(1:return_sizes{rati}(1)); % extract
    zIdPhi_all{rati}(1:return_sizes{rati}(1))=[]; % erase this data so that we can continue extracting
    return_sizes{rati}(1) = [];
    
    % now do saline condition within saline
    zIdPhi.(rat{rati}).Saline.Saline = zIdPhi_all{rati}(1:return_sizes{rati}(1));
    zIdPhi_all{rati}(1:return_sizes{rati}(1))=[]; % erase this data so that we can continue extracting
    return_sizes{rati}(1) = [];
    
    % baseline - 
    zIdPhi.(rat{rati}).Baseline.Baseline1 = zIdPhi_all{rati}(1:return_sizes{rati}(1)); % extract
    zIdPhi_all{rati}(1:return_sizes{rati}(1))=[]; % erase this data so that we can continue extracting
    return_sizes{rati}(1) = [];
    
    % now do saline condition within saline
    zIdPhi.(rat{rati}).Baseline.Baseline2 = zIdPhi_all{rati}(1:return_sizes{rati}(1));
    zIdPhi_all{rati}(1:return_sizes{rati}(1))=[]; % erase this data so that we can continue extracting
    return_sizes{rati}(1) = [];    
  
    % muscimol
    zIdPhi.(rat{rati}).Muscimol.Baseline = zIdPhi_all{rati}(1:return_sizes{rati}(1)); % extract
    zIdPhi_all{rati}(1:return_sizes{rati}(1))=[]; % erase this data so that we can continue extracting
    return_sizes{rati}(1) = [];
    
    % now do saline condition within saline
    zIdPhi.(rat{rati}).Muscimol.Muscimol = zIdPhi_all{rati}(1:return_sizes{rati}(1));
    zIdPhi_all{rati}(1:return_sizes{rati}(1))=[]; % erase this data so that we can continue extracting
    return_sizes{rati}(1) = [];      
end

% house keeping
%clearvars -except zIdPhi zIdPhi_dist IdPhi xPos xPosOG yPos yPosOG tsPos tsPosOG zIdPhi_dist rat Day_cond zThresh

%% determine a threshold using redish method with gaussian
%{
% what z-score is in the 75th percentile, per rat?
figure('color','w')
h1 = histogram(zIdPhi_dist); hold on;
h1.FaceColor = 'r';
ylimits = ylim;
title('Distribution of zIdPhi estimates across sessions')
box off
ylabel('Number of Trials')
xlabel('zIdPhi')
h1.BinWidth = 0.1;
bins = h1.BinEdges;
%}

% try fitting a complex 
fig = figure('color','w');
h1 = histogram(zIdPhi_dist); hold on;
%hist(zIdPhi_dist); hold on;
%h1 = findobj(gca,'Type','patch');
h1.FaceColor = 'r';
ylimits = ylim;
title('Distribution of z[ln(IdPhi)] estimates across sessions')
box off
ylabel('Number of Trials')
xlabel('z[ln(IdPhi)]')
h1.BinWidth = 0.1;
zIdPhi_bins = h1.BinEdges;
axis tight;
% fit kernel distirbution 
%pd = fitdist(zIdPhi_dist','Kernel','Kernel','epanechnikov');
pd = fitdist(zIdPhi_dist','Kernel','Kernel','normal');
stepVal = 1/length(zIdPhi_dist);
x_curve = h1.BinEdges(1):stepVal:h1.BinEdges(end);
y_curve = pdf(pd,x_curve);
yyaxis right; hold on;
plot(x_curve,y_curve,'k','LineWidth',2); 

% -- identify the first largely positive deflection, this will be our
% threshold -- %
% find first point of deflection after second half of peak in curve
y_axis = h1.BinCounts;
x_axis = h1.BinEdges;
[maxVal,maxIdx] = max(y_axis); % use the index for the max val to extract x-axis
maxBins = x_axis(maxIdx:maxIdx+1); % +1 because we're working with bins
% now using maxBins(2) to identify the distribution after the most common
% value
y_curve_remainder = y_curve(x_curve > maxBins(2));
x_curve_remainder = x_curve(x_curve > maxBins(2));
% take derivative to find first point of positive deflection
dy_curve_remainder = gradient(y_curve_remainder);
positiveDeflections = find(dy_curve_remainder > 0);
firstPosDeflection = positiveDeflections(1);
binFirstPosDef = round(x_curve_remainder(firstPosDeflection),1); % round to the first decimal to the right

% ---- UPDATE ---- %
% 9-17-2020, after visual inspection of behavior, time spent, and head
% velocity, a threshold of 1 seems better and corresponds to the second
% bump in the distribution. This should ensure minimal false positive VTE
% classifications
% ----  END   ---- %

% find the second maximima in the VTE distribution to ensure minimal false
% positives (this is overly conservative). Get the bin of first maximima
% after the bin of first positive deflection
%peaks()
x_after1posDef = x_curve_remainder(x_curve_remainder > binFirstPosDef);
y_after1posDef = y_curve_remainder(x_curve_remainder > binFirstPosDef);

% we'll consider the VTE distribution the data after the first positive
% deflection point, indicating some change in the normal distribution to
% the left
yyaxis right;
p2_right = plot(x_after1posDef,y_after1posDef,'b','LineWidth',2);
p2_right.LineStyle = '-';

% use peaks to find local maxima
idxOfy_maxima = findpeaks(y_after1posDef');
xOFy_maxima   = x_after1posDef(idxOfy_maxima.loc);

% extract the second local maxima to be well within the distribution of
% VTEs to avoid false positives - false positives are worse than false
% negatives. The second local maxima defines our VTE threshold.
if idxOfy_maxima.loc(1) == 1 % for some reason, I got a weird thing where the local maxima was detected as the first value of the distribution
    threshold = xOFy_maxima(2);
else    
    threshold = xOFy_maxima(1); % this corresponds to the z-score
end
% round threshold
threshold = round(threshold,2);
% threshold will be the bin that corresponds to the first positive
% deflection in the data. this was chosen because it should signify some
% change in the distribution. I used the normal curve to elucidate this
% deflection to detect the general trend.
% plot a line
yyaxis left;
ylimits = ylim;
l1 = line([threshold threshold],[ylimits(1) ylimits(2)]);
l1.Color = 'b';
l1.LineStyle = '--';
l1.LineWidth = 1;
text([threshold+0.1],[ylimits(2)/2],['VTE threshold, the 1st local maxima' ...
     newline 'in the VTE distribution (blue), is ',num2str(threshold)])
% set color of axes
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'b';
% remove right axis - not necessary
ax.YAxis(2).Visible = 'off';
% set painters so you can edit it in illustrator, save file
cd(place2store)
print -painters -depsc fig_zlnIdPhi_distribution.eps

%% second definition
% define based on timespent distribution and vte behavior
for rati = 1:length(rat)
    
    % col 1 is baseline, col 2 is infusion
    tsPos_bl_saline{rati} = tsPos.(rat{rati}).Saline.Baseline;
    tsPos_sa_saline{rati} = tsPos.(rat{rati}).Saline.Saline;
    
    % baseline
    tsPos_bl_baseline1{rati} = tsPos.(rat{rati}).Baseline.Baseline1;
    tsPos_bl_baseline2{rati} = tsPos.(rat{rati}).Baseline.Baseline2;
  
    % muscimol
    tsPos_mu_baseline{rati} = tsPos.(rat{rati}).Muscimol.Baseline;
    tsPos_mu_muscimol{rati} = tsPos.(rat{rati}).Muscimol.Muscimol;
    
    % concatenate data
    tsPos_data{rati} = horzcat(tsPos_bl_saline{rati},tsPos_sa_saline{rati},tsPos_bl_baseline1{rati},...
        tsPos_bl_baseline2{rati},tsPos_mu_baseline{rati},tsPos_mu_muscimol{rati});
    
end
 
% cat across all observations and create a distribution of timespents
tsPos_data = horzcat(tsPos_data{:});
for i = 1:length(tsPos_data)
    timeSpent_array(i) = (tsPos_data{i}(end)-tsPos_data{i}(1))/1e6;
end

% consider the time if its 1 std above the mean
percentile = 75; % qualitatively seems to match up for the most part
%std_aboveTS = 1;
fig_timeSpent = figure('color','w');
h1 = histogram(timeSpent_array); hold on;
h1.FaceColor = [0 0.5 0.5];
ylimits = ylim;
box off
ylabel('Number of Trials')
xlabel('Timespent at CP (seconds)')
h1.BinWidth = 1.5;
ts_percentile = prctile(timeSpent_array,percentile);
%ts_std = zscore(timeSpent_array);
%ts_1std = timeSpent_array(dsearchn(ts_std',1));
axis tight
ylimits = ylim;
l_ts = line([ts_percentile ts_percentile],[ylimits(1) ylimits(2)]);
l_ts.LineStyle = '--';
l_ts.Color = 'k';
l_ts.LineWidth = 1.5;
text([ts_percentile+1],[ylimits(2)/2],['Time-spent at ',num2str(percentile),'th',' percentile = ',num2str(ts_percentile)])

% set painters so you can edit it in illustrator, save file
cd(place2store)
print -painters -depsc fig_timeSpent_distribution.eps

% second threshold
threshold2 = ts_percentile;

%% plot examples
% after looking through rat1422 with timespent and angular velocity, I feel 
% like a score of 1 safely gets VTE events. Is there a way to draw this
% from the distribution? I need to look at more rats. Remember when
% interpreting velocity colors, if you see a solid color, it just means
% there was no change in velocity, not that they were going slow/fast. What
% would this metrics units be? rad/sec?

%threshold = 0.5;
% next, lets plot some examples
rat_get = 'Usher';    % rat
day_get = 'Muscimol'; % day (baseline, saline, muscimol day)
con_get = 'Muscimol'; % condition (saline, baseline, or muscimol infusion)

% get session IdPhi
session_zIdPhi = zIdPhi.(rat_get).(day_get).(con_get);

% get session position information
session_x = xPosOG.(rat_get).(day_get).(con_get);
session_y = yPosOG.(rat_get).(day_get).(con_get);
session_ts = tsPosOG.(rat_get).(day_get).(con_get);

% plot location specific data
loc_x = xPos.(rat_get).(day_get).(con_get);
loc_y = yPos.(rat_get).(day_get).(con_get);
loc_ts = tsPos.(rat_get).(day_get).(con_get);
%angVel = dphi.(rat_get).(day_get).(con_get); % not sure if this is actually angular velocity

% plot based off distribution - 0 seems a good spot
for i = 1:length(zIdPhi_bins)-1
    idx2plot = find(session_zIdPhi >= zIdPhi_bins(i) & session_zIdPhi <= zIdPhi_bins(i+1));
    if isempty(idx2plot) == 1
        disp(['No events between zIdPhi of ',num2str(zIdPhi_bins(i)),' and ',num2str(zIdPhi_bins(i+1))]);
       continue
    else
        for ii = 1:numel(idx2plot)
            fig_vtes = figure('color','w');
            p1 = plot(session_x,session_y,'Color',[.8 .8 .8]); %box off;
            ylim([-100 500]);
            hold on;             
            %plot(loc_x{idx2plot(ii)},loc_y{idx2plot(ii)},'b','LineWidth',1.4)
            
            % plot position with color heat indicating angular velocity
            z = zeros(size(loc_x{idx2plot(ii)}));
            %col = angVel{idx2plot(ii)}';
            [~,head_velocity] = instant_speed2(loc_x{idx2plot(ii)},loc_y{idx2plot(ii)},loc_ts{idx2plot(ii)},'y'); % 'y' to convert to seconds
            head_velocity = normalize(head_velocity,'range');
            
            %figure()            
            %s = surf([loc_x{idx2plot(ii)};loc_x{idx2plot(ii)}],[loc_y{idx2plot(ii)};loc_y{idx2plot(ii)}],[z;z],[head_velocity;head_velocity],...
                %'EdgeColor','interp','LineWidth',2)
            
            %surf([loc_x{idx2plot(ii)};loc_x{idx2plot(ii)}],[loc_y{idx2plot(ii)};loc_y{idx2plot(ii)}],[head_velocity;head_velocity],...
                %'EdgeColor','interp','LineWidth',2)    
            s = scatter(loc_x{idx2plot(ii)},loc_y{idx2plot(ii)},[],head_velocity);
            s.Marker = '.';
            s.SizeData = 75;
            
            
            %{
            surface([loc_x{idx2plot(ii)};loc_x{idx2plot(ii)}],[loc_y{idx2plot(ii)};loc_y{idx2plot(ii)}],[z;z],[head_velocity;head_velocity],...
                        'facecol','no',...
                        'edgecol','interp',...
                        'linew',2); 
            %}
            % trying another way to plot color since surface gives so many
            % issues with illustrator
            
            %{
            colorMap = parula(numel(head_velocity));
            sortedVel = sort(head_velocity); % sort velocity to make an index
            for iii = 1:length(head_velocity)-1
                colorDot_idx = ([find(head_velocity(iii) == sortedVel),find(head_velocity(iii+1) == sortedVel)]);
                colorDot(iii,:) = mean(colorMap(colorDot_idx,:),1);
                plot([loc_x{idx2plot(ii)}(iii),loc_x{idx2plot(ii)}(iii+1)],[loc_y{idx2plot(ii)}(iii),loc_y{idx2plot(ii)}(iii+1)] ,'-','Color',colorDot(iii,:),'LineWidth',2)
                shading interp
            end
            %}
      
            if session_zIdPhi(idx2plot(ii)) > threshold
                text_thresh_indicator = 'VTE';
            else
                text_thresh_indicator = 'Non VTE';
            end
            title(['Bins: ',num2str(zIdPhi_bins(i)),':',num2str(zIdPhi_bins(i+1)),' | ','z[ln(IdPhi)] = ',num2str(session_zIdPhi(idx2plot(ii))),...
                newline ' Classified as a ',text_thresh_indicator])
             
            % orient axis
            xlim([400 700]);
            ylim([120 330]);
            
            % plot time spent
            timeSpent_run = (loc_ts{idx2plot(ii)}(end)-loc_ts{idx2plot(ii)}(1))/1e6;
            xlimits = xlim;
            ylimits = ylim;
            text([xlimits(2)/1.5],[ylimits(2)/1.2],['TimeSpent at CP = ',num2str(timeSpent_run),' sec.'])
                    
            % remove axes - they're distracting
            set(gca,'xtick',[]);
            set(gca,'ytick',[]);
            
            % size of fig
            set(gcf,'Position',[300 250 275 260])
            
            % saving data
            prompt = 'Do you want to save this fig? [Y/N] ';
            answer = input(prompt,'s');
            if contains(answer,'y') | contains(answer,'Y')
                savePrompt = 'Enter a name to save the file ';
                saveName = input(savePrompt,'s')
                %print('-painters','-dsvg',[saveName,'.svg'],'-r0')
                print('-painters',[saveName,'.eps'],'-depsc','-r0')
                %print2eps(saveName, fig_vtes)
            end
            disp('Press any key to continue')
            pause; 
            close;
        end
    end
end

%% ______________ ANALYSIS MODE _______________ %%
%% QUESTION 1

% -- Do rats deliberate more without Re? -- %

% define a threshold
%threshold = 1; % zscore
VTE=[]; Non=[]; trial_VTE = []; propVTE = [];
for rati = 1:length(rat)
    
    VTE(rati,1) = numel(find(zIdPhi.(rat{rati}).Baseline.Baseline1 > threshold));
    VTE(rati,2) = numel(find(zIdPhi.(rat{rati}).Baseline.Baseline2 > threshold));    
    VTE(rati,3) = numel(find(zIdPhi.(rat{rati}).Saline.Baseline > threshold));
    VTE(rati,4) = numel(find(zIdPhi.(rat{rati}).Saline.Saline > threshold));
    VTE(rati,5) = numel(find(zIdPhi.(rat{rati}).Muscimol.Baseline > threshold));
    VTE(rati,6) = numel(find(zIdPhi.(rat{rati}).Muscimol.Muscimol > threshold));

    Non(rati,1) = numel(find(zIdPhi.(rat{rati}).Baseline.Baseline1 < threshold));
    Non(rati,2) = numel(find(zIdPhi.(rat{rati}).Baseline.Baseline2 < threshold));    
    Non(rati,3) = numel(find(zIdPhi.(rat{rati}).Saline.Baseline < threshold));
    Non(rati,4) = numel(find(zIdPhi.(rat{rati}).Saline.Saline < threshold));
    Non(rati,5) = numel(find(zIdPhi.(rat{rati}).Muscimol.Baseline < threshold));
    Non(rati,6) = numel(find(zIdPhi.(rat{rati}).Muscimol.Muscimol < threshold));   
    
    % get a percentage
    propVTE(rati,1) = (numel(find(zIdPhi.(rat{rati}).Baseline.Baseline1 > threshold)))/(numel(zIdPhi.(rat{rati}).Baseline.Baseline1));
    propVTE(rati,2) = (numel(find(zIdPhi.(rat{rati}).Baseline.Baseline2 > threshold)))/ (numel(zIdPhi.(rat{rati}).Baseline.Baseline2));    
    propVTE(rati,3) = (numel(find(zIdPhi.(rat{rati}).Saline.Baseline > threshold)))/ (numel(zIdPhi.(rat{rati}).Saline.Baseline));
    propVTE(rati,4) = (numel(find(zIdPhi.(rat{rati}).Saline.Saline > threshold)))/(numel(zIdPhi.(rat{rati}).Saline.Saline));
    propVTE(rati,5) = (numel(find(zIdPhi.(rat{rati}).Muscimol.Baseline > threshold)))/(numel(zIdPhi.(rat{rati}).Muscimol.Baseline));
    propVTE(rati,6) = (numel(find(zIdPhi.(rat{rati}).Muscimol.Muscimol > threshold)))/(numel(zIdPhi.(rat{rati}).Muscimol.Muscimol)); 
    
    % trial-based distribution
    trial_VTE.baseline.baseline1{rati} = zIdPhi.(rat{rati}).Baseline.Baseline1;
    trial_VTE.baseline.baseline2{rati} = zIdPhi.(rat{rati}).Baseline.Baseline2;    
    trial_VTE.saline.baseline{rati} = zIdPhi.(rat{rati}).Saline.Baseline;
    trial_VTE.saline.saline{rati} = zIdPhi.(rat{rati}).Saline.Saline;
    trial_VTE.muscimol.baseline{rati} = zIdPhi.(rat{rati}).Muscimol.Baseline;
    trial_VTE.muscimol.muscimol{rati} = zIdPhi.(rat{rati}).Muscimol.Muscimol;    
end

% LEFT OFF - check Mortys data, he has the same proportion of VTEs across
% multiple conditions, kinda weird

% notice that we're not seeing an effect on the single rat scale - probably
% because there are so few trials and so few rats, therefore lets look
% across sessions
figure()
colorLine = [];% colorLine.decrease = 'm'; colorLine.increase = 'b'; colorLine.noChange = [.8 .8 .8]
colorLine = 1; % default colors
jitterIdx = 1:7;
colorIdx{1} = 'r'; colorIdx{2} = 'b'; colorIdx{3} = 'm'; colorIdx{4} = 'y'; colorIdx{5} = 'k'; colorIdx{6} = 'g'; colorIdx{7} = 'c';
BarPlotsJitteredData(propVTE,1,0,1,1,[],colorLine,jitterIdx,colorIdx)
set(gca,'FontSize',12)
ylabel('Proportion of Deliberative-like events');
ax = gca;
ax.XTick = [1:6];
ax.XTickLabel = [{'Baseline 1'},{'Baseline 2'},{'Baseline'},{'Saline'},{'Baseline'},{'Muscimol'}];
ax.XTickLabelRotation = 45;
% plot legend with each rat
f=get(gca,'Children');
% had to find which f data corresponded to which rat. I compared this to
% the rat cell array to ensure the colors were correct
legend(f(22:22+7),flipud(rat')','Location','NorthWest')
title('QUESTION: Does Re inactivation induce deliberation?')

% like henry, do a change score
prop_baselineDay = (propVTE(:,2)-propVTE(:,1))./(propVTE(:,2)+propVTE(:,1));
prop_salineDay   = (propVTE(:,4)-propVTE(:,3))./(propVTE(:,4)+propVTE(:,3));
prop_muscimolDay = (propVTE(:,6)-propVTE(:,5))./(propVTE(:,6)+propVTE(:,5));
diffScores = horzcat(prop_baselineDay,prop_salineDay,prop_muscimolDay);

figure()
jitterIdx = 1:7;
colorIdx{1} = 'r'; colorIdx{2} = 'b'; colorIdx{3} = 'm'; colorIdx{4} = 'y'; colorIdx{5} = 'k'; colorIdx{6} = 'g'; colorIdx{7} = 'c';
fig = BarPlotsJitteredData(diffScores,1,0,1,0,[],[],jitterIdx,colorIdx)
set(gca,'FontSize',12)
ylabel('Deliberative Difference (testing - baseline)');
ax = gca;
ax.XTick = [1:3];
ax.XTickLabel = [{'Baseline'},{'Saline'},{'Muscimol'}];
ax.XTickLabelRotation = 45;
[h,p]=ttest(diffScores,0)
% below isnt correct yet
f=get(gca,'Children');
% had to find which f data corresponded to which rat. I compared this to
% the rat cell array to ensure the colors were correct
legend(f(1:1+7),flipud(rat')','Location','NorthWest')

%{
% -- using trials, find prop of vte events -- %
trialFormat = [];
trialFormat.baseline_data{1} = horzcat(trial_VTE.baseline.baseline1{:});
trialFormat.baseline_ratTrials{1} = cellfun(@numel,trial_VTE.baseline.baseline1);
trialFormat.baseline_data{2} = horzcat(trial_VTE.baseline.baseline2{:});
trialFormat.baseline_ratTrials{2} = cellfun(@numel,trial_VTE.baseline.baseline2);
%}

% ANSWER 1: 

%% QUESTION
% Ratdle does worse no matter what on the second condition, so what happens
% without his data?
% notice that we're not seeing an effect on the single rat scale - probably
% because there are so few trials and so few rats, therefore lets look
% across sessions
propVTE_noRatdle = propVTE;
propVTE_noRatdle(find(contains(rat,'Ratdle')==1),:)=[];
figure()
colorLine = [];% colorLine.decrease = 'm'; colorLine.increase = 'b'; colorLine.noChange = [.8 .8 .8]
colorLine = 1; % default colors
jitterIdx = 1:6;
colorIdx{1} = 'r'; colorIdx{2} = 'b'; colorIdx{3} = 'm'; colorIdx{4} = 'k'; colorIdx{5} = 'g'; colorIdx{6} = 'c';
BarPlotsJitteredData(propVTE_noRatdle,1,0,1,1,[],colorLine,jitterIdx,colorIdx)
set(gca,'FontSize',12)
ylabel('Proportion of Deliberative-like events');
ax = gca;
ax.XTick = [1:6];
ax.XTickLabel = [{'Baseline 1'},{'Baseline 2'},{'Baseline'},{'Saline'},{'Baseline'},{'Muscimol'}];
ax.XTickLabelRotation = 45;
% plot legend with each rat
f=get(gca,'Children');
% had to find which f data corresponded to which rat. I compared this to
% the rat cell array to ensure the colors were correct
rat_noRatdle = rat;
rat_noRatdle(find(contains(rat,'Ratdle')==1))=[];
legend(f(19:19+6),flipud(rat_noRatdle')','Location','NorthWest')
title('QUESTION: Does Re inactivation induce deliberation? Exclude Ratdle')

% like henry, do a change score
prop_baselineDay_noRatdle = (propVTE_noRatdle(:,2)-propVTE_noRatdle(:,1))./(propVTE_noRatdle(:,2)+propVTE_noRatdle(:,1));
prop_salineDay_noRatdle   = (propVTE_noRatdle(:,4)-propVTE_noRatdle(:,3))./(propVTE_noRatdle(:,4)+propVTE_noRatdle(:,3));
prop_muscimolDay_noRatdle = (propVTE_noRatdle(:,6)-propVTE_noRatdle(:,5))./(propVTE_noRatdle(:,6)+propVTE_noRatdle(:,5));
diffScores_noRatdle = horzcat(prop_baselineDay_noRatdle,prop_salineDay_noRatdle,prop_muscimolDay_noRatdle);

figure()
jitterIdx = 1:6;
colorIdx{1} = 'r'; colorIdx{2} = 'b'; colorIdx{3} = 'm'; colorIdx{4} = 'k'; colorIdx{5} = 'g'; colorIdx{6} = 'c';
fig = BarPlotsJitteredData(diffScores_noRatdle,1,0,1,0,[],[],jitterIdx,colorIdx)
set(gca,'FontSize',12)
ylabel('Deliberative Difference (testing - baseline)');
ax = gca;
ax.XTick = [1:3];
ax.XTickLabel = [{'Baseline'},{'Saline'},{'Muscimol'}];
ax.XTickLabelRotation = 45;
f=get(gca,'Children');
% had to find which f data corresponded to which rat. I compared this to
% the rat cell array to ensure the colors were correct
legend(f(1:1+5),flipud(rat_noRatdle')','Location','NorthWest')

[h,p]=ttest(diffScores_noRatdle,0)

%% QUESTION 2
% what about vte events following or corresponding to an incorrect
% decision?

%% QUESTION 2
% -- in terms of deliberation, within a rat, do certain rats deliberate
% more significantly without Re? This would have to be compared across
% multiple conditions -- %

%% QUESTION 3
% -- Do rats revert to perseveration strategies without Re? -- %

% have to redefine this
rat{6} = '14-22';

persevIdx = []; alternIdx = [];
for rat_i = 1:length(rat)
    
    % datafolder for the rat
    datafolder = [];
    datafolder_rat = [Datafolders,'\',rat{rat_i}];
    
    for day_cond_i = 1:length(Day_cond)
        
        % datafolder that includes condition
        datafolder_cond = [datafolder_rat,'\',Day_cond{day_cond_i}];
        
        % get content in directory
        dir_content = dir(datafolder_cond);
        dir_content = extractfield(dir_content,'name');
 
        % find non folders
        remIdx = contains(dir_content,'.mat') | contains(dir_content,'.');
        
        % remove .mat files
        dir_content(remIdx)=[];
        
        % loop across infusion conditions
        for inf_i = 1:length(dir_content)
        
            % define datafolder
            datafolder = [datafolder_cond,'\',dir_content{inf_i}];

            % the struct array below that organizes the data cannot create
            % fields with space characters, therefore rename them if they
            % exist.
            if contains(dir_content{inf_i},'Baseline 1')
                dir_content{inf_i} = 'Baseline1';
            elseif contains(dir_content{inf_i},'Baseline 2')
                dir_content{inf_i} = 'Baseline2';
            end
            
            % if rat has a '-' character in it, remove
            if contains(rat{rat_i},'14-22')
                rat{rat_i} = 'rat1422';
            end
            
            % get alternation and perseveration indices
            [alternIdx.(rat{rat_i}).(Day_cond{day_cond_i}).(dir_content{inf_i}),...
            persevIdx.(rat{rat_i}).(Day_cond{day_cond_i}).(dir_content{inf_i}),...
            traj_change.(rat{rat_i}).(Day_cond{day_cond_i}).(dir_content{inf_i})] = ...
            alternPersevIdx(datafolder,int_name);
        
            % performance
            performance_accuracy.(rat{rat_i}).(Day_cond{day_cond_i}).(dir_content{inf_i}) = ...
            performanceAccuracy(datafolder,int_name);            
            
            % restore appropriate naming for 14-22
            if contains(rat{rat_i},'rat1422')
                rat{rat_i} = '14-22';
            end
                    
        end       
    end
end

% have to redefine this rat
rat{6} = 'rat1422';

%% QUESTION: is performance affected without Re?
percent_accurate=[];
for rati = 1:length(rat)    
    % get performance
    percent_accurate(rati,1) = performance_accuracy.(rat{rati}).Baseline.Baseline1;
    percent_accurate(rati,2) = performance_accuracy.(rat{rati}).Baseline.Baseline2;
    percent_accurate(rati,3) = performance_accuracy.(rat{rati}).Saline.Baseline;
    percent_accurate(rati,4) = performance_accuracy.(rat{rati}).Saline.Saline;
    percent_accurate(rati,5) = performance_accuracy.(rat{rati}).Muscimol.Baseline;
    percent_accurate(rati,6) = performance_accuracy.(rat{rati}).Muscimol.Muscimol;    
end

colorLine = [];% colorLine.decrease = 'm'; colorLine.increase = 'b'; colorLine.noChange = [.8 .8 .8]
colorLine = 1; % default colors
jitterIdx = 1:7;
colorIdx{1} = 'r'; colorIdx{2} = 'b'; colorIdx{3} = 'm'; colorIdx{4} = 'y'; colorIdx{5} = 'k'; colorIdx{6} = 'g'; colorIdx{7} = 'c';
BarPlotsJitteredData(percent_accurate,1,0,1,1,[],colorLine,jitterIdx,colorIdx)
set(gca,'FontSize',12)
ylabel('percent accuracy (%)');
ax = gca;
ax.XTick = [1:6];
ax.XTickLabel = [{'Baseline 1'},{'Baseline 2'},{'Baseline'},{'Saline'},{'Baseline'},{'Muscimol'}];
ax.XTickLabelRotation = 45;
% plot legend with each rat
f=get(gca,'Children');
% had to find which f data corresponded to which rat. I compared this to
% the rat cell array to ensure the colors were correct
legend(f(22:22+7),flipud(rat')','Location','NorthWest')
title('Question: Does Re inactivation affect behavior? - replication')


%% QUESTION: do rats perseverate more without Re?
persev=[]; percent_accurate=[];
for rati = 1:length(rat)    
    % get perveration data
    persev(rati,1) = persevIdx.(rat{rati}).Baseline.Baseline1;
    persev(rati,2) = persevIdx.(rat{rati}).Baseline.Baseline2;
    persev(rati,3) = persevIdx.(rat{rati}).Saline.Baseline;
    persev(rati,4) = persevIdx.(rat{rati}).Saline.Saline;
    persev(rati,5) = persevIdx.(rat{rati}).Muscimol.Baseline;
    persev(rati,6) = persevIdx.(rat{rati}).Muscimol.Muscimol;    
end

colorLine = [];% colorLine.decrease = 'm'; colorLine.increase = 'b'; colorLine.noChange = [.8 .8 .8]
colorLine = 1; % default colors
jitterIdx = 1:7;
colorIdx{1} = 'r'; colorIdx{2} = 'b'; colorIdx{3} = 'm'; colorIdx{4} = 'y'; colorIdx{5} = 'k'; colorIdx{6} = 'g'; colorIdx{7} = 'c';
BarPlotsJitteredData(persev,1,0,1,1,[],colorLine,jitterIdx,colorIdx)
set(gca,'FontSize',12)
ylabel('Proportion of perseveration behaviors');
ax = gca;
ax.XTick = [1:6];
ax.XTickLabel = [{'Baseline 1'},{'Baseline 2'},{'Baseline'},{'Saline'},{'Baseline'},{'Muscimol'}];
ax.XTickLabelRotation = 45;
% plot legend with each rat
f=get(gca,'Children');
% had to find which f data corresponded to which rat. I compared this to
% the rat cell array to ensure the colors were correct
legend(f(22:22+7),flipud(rat')','Location','NorthWest')
title('Question: Does Re inactivation induce perseveration?')

%% QUESTION
% whats the relationship between perseveration behaviors and VTEs?

% globally, what is the relationship? This considers all data
% use propVTE, make into vector, compare to vectorized persev
persev_vector = persev(:)';
propVTE_vector = propVTE(:)';
[r,p] = corrcoef(persev_vector,propVTE_vector);

figure('color','w');
scat1 = scatter(persev_vector,propVTE_vector);
scat1.MarkerEdgeColor = 'k';
scat1.MarkerFaceColor = [.8 .8 .8];
l_scat1 = lsline;
l_scat1.Color = 'r';
axis tight
ylabel('Proportion of VTE events per session')
xlabel('Proportion of Perseverations per session')
title('QUESTION: Whats the relationship between VTEs and perseverations?')

% answer: there is no significant relationship between vtes and
% perseveration

%% QUESTION
% what is the relationship between vtes and alternation? Does Re
% inactivation result in random alternations?
altern=[]; 
for rati = 1:length(rat)    
    % get perveration data
    altern(rati,1) = alternIdx.(rat{rati}).Baseline.Baseline1;
    altern(rati,2) = alternIdx.(rat{rati}).Baseline.Baseline2;
    altern(rati,3) = alternIdx.(rat{rati}).Saline.Baseline;
    altern(rati,4) = alternIdx.(rat{rati}).Saline.Saline;
    altern(rati,5) = alternIdx.(rat{rati}).Muscimol.Baseline;
    altern(rati,6) = alternIdx.(rat{rati}).Muscimol.Muscimol;
end

colorLine = [];% colorLine.decrease = 'm'; colorLine.increase = 'b'; colorLine.noChange = [.8 .8 .8]
colorLine = 1; % default colors
jitterIdx = 1:7;
colorIdx{1} = 'r'; colorIdx{2} = 'b'; colorIdx{3} = 'm'; colorIdx{4} = 'y'; colorIdx{5} = 'k'; colorIdx{6} = 'g'; colorIdx{7} = 'c';
BarPlotsJitteredData(altern,1,0,1,1,[],colorLine,jitterIdx,colorIdx)
set(gca,'FontSize',12)
ylabel('Proportion of alternation behaviors');
ax = gca;
ax.XTick = [1:6];
ax.XTickLabel = [{'Baseline 1'},{'Baseline 2'},{'Baseline'},{'Saline'},{'Baseline'},{'Muscimol'}];
ax.XTickLabelRotation = 45;
% plot legend with each rat
f=get(gca,'Children');
% had to find which f data corresponded to which rat. I compared this to
% the rat cell array to ensure the colors were correct
legend(f(22:22+7),flipud(rat')','Location','NorthWest')
title('Question: Does Re inactivation induce random alternation?')

% answer: Re inactivation, as is seen in the perseveration plot, does not
% results in random alternation (50%). Instead, rats revert to a
% perseverative strategy. Is this a ballistic strat?

altern_vector = altern(:)';
[r,p] = corrcoef(altern_vector,propVTE_vector);

figure('color','w');
scat2 = scatter(altern_vector,propVTE_vector);
scat2.MarkerEdgeColor = 'k';
scat2.MarkerFaceColor = [.8 .8 .8];
l_scat2 = lsline;
l_scat2.Color = 'r';
axis tight
ylabel('Proportion of VTE events per session')
xlabel('Proportion of alternations per session')
title('QUESTION: Whats the relationship between VTEs and alternations?')

%% QUESTION
% With vs without Re, what is the relationship between VTEs and
% perseverations

% get data from Re inactivation day - maybe I should examine all trials
% somehow?
withRe_persev  = persev(:,5);
woutRe_persev  = persev(:,6);
withRe_propVTE = propVTE(:,5);
woutRe_propVTE = propVTE(:,6);

[r,p] = corrcoef(withRe_persev,withRe_propVTE);
[r,p] = corrcoef(woutRe_persev,woutRe_propVTE);

figure('color','w');
subplot 211
    scat3 = scatter(woutRe_persev,woutRe_propVTE);
    scat3.MarkerEdgeColor = 'k';
    scat3.MarkerFaceColor = [.8 .8 .8];
    l_scat3 = lsline;
    l_scat3.Color = 'r';
    axis tight
    ylabel('Proportion of VTE events without Re')
    xlabel('Proportion of perseverations without Re')
    title('QUESTION: Whats the relationship between VTEs and perseverations with/without Re?')
subplot 212
    scat4 = scatter(withRe_persev,withRe_propVTE);
    scat4.MarkerEdgeColor = 'k';
    scat4.MarkerFaceColor = [.8 .8 .8];
    l_scat4 = lsline;
    l_scat4.Color = 'r';
    axis tight
    ylabel('Proportion of VTE events with Re')
    xlabel('Proportion of perseverations with Re')
    
%% QUESTION
% what behavioral feature are vtes linked to in control conditions? is the
% same true ewithout re?

%% QUESTION
% after an incorrect decision

%% QUESTION
% for each rat, are they performing worse after Re inactivation? Maybe
% deliberation is a product of Re inactivation, but Re isn't directly
% supporting VTEs

corAlpha = 0.05/3;

for rati = 1:length(rat)   
    [h_musc(rati),p_musc(rati)] = ttest2(zIdPhi.(rat{rati}).Muscimol.Baseline,zIdPhi.(rat{rati}).Muscimol.Muscimol,'Alpha',corAlpha)
    [h_sali(rati),p_sali(rati)] = ttest2(zIdPhi.(rat{rati}).Saline.Baseline,zIdPhi.(rat{rati}).Saline.Saline,'Alpha',corAlpha)
    [h_base(rati),p_base(rati)] = ttest2(zIdPhi.(rat{rati}).Baseline.Baseline1,zIdPhi.(rat{rati}).Baseline.Baseline2,'Alpha',corAlpha)
end

figure('color','w'); hold on;
bar(1,numel(find(h_base == 1))); 
bar(2,numel(find(h_sali == 1))); 
bar(3,numel(find(h_musc == 1))); 
ylim([0 numel(rat)])
ax = gca;
ax.XTick = [1:3];
ax.XTickLabel = [{'Baseline'},{'Saline'},{'Muscimol'}];
ax.XTickLabelRotation = 45;
ylabel(['# rats with significant change ' newline 'after experimental condition'])
xlabel('Condition')

% matrix of effects

%numRatsVTE = table([;1],[6;7],'VariableNames',{'VTE','NonVTE'},'RowNames',{'Baseline','Experiment'})

%% QUESTION 
% whats the relationship between VTEs and choice accuracy per session?
accuracy_vector = percent_accurate(:)';
propVTE_vector = propVTE(:)';
[r,p] = corrcoef(accuracy_vector,propVTE_vector);

figure('color','w');
scat1 = scatter(accuracy_vector,propVTE_vector);
scat1.MarkerEdgeColor = 'k';
scat1.MarkerFaceColor = [.8 .8 .8];
l_scat1 = lsline;
l_scat1.Color = 'r';
axis tight
ylabel('Proportion of VTE events per session')
xlabel('Percent accuracy per session')
title('QUESTION: Whats the relationship between VTEs and choice accuracy?')

%% anova code setup
propVTE_switch = propVTE';
propVTE_setup  = propVTE_switch(:);

accuracy_switch = percent_accurate';
accuracy_setup  = accuracy_switch(:);

persev_switch = persev';
persev_setup  = persev_switch(:);

% now take this, copy and paste into the excel sheet with the appropriate
% name

%% within condition ttests
corAlpha = 0.05/3;
cond_idx = 2:2:size(propVTE,2);
h_vte = []; p_vte = []; h_per = []; p_per = []; h_acc = []; p_acc = [];
for condi = 1:length(cond_idx)
    [h_vte(condi),p_vte(condi),ci,stat_vte{condi}] = ttest(propVTE(:,cond_idx(condi)-1),propVTE(:,cond_idx(condi)),'Alpha',corAlpha);
    [h_per(condi),p_per(condi),ci,stat_per{condi}] = ttest(persev(:,cond_idx(condi)-1),persev(:,cond_idx(condi)),'Alpha',corAlpha);
    [h_acc(condi),p_acc(condi),ci,stat_acc{condi}] = ttest(percent_accurate(:,cond_idx(condi)-1),percent_accurate(:,cond_idx(condi)),'Alpha',corAlpha);
end

% corrected pvalue
p_vte = p_vte*3;
p_per = p_per*3;
p_acc = p_acc*3;

%% QUESTION: do VTEs significantly differ on alternation trials that follow an error or perseveration (synonymous)
% use traj_change and vte variables. However, make sure you only look at
% the vte variables from trial 2 and on. This is because traj_change is the
% quantification of left/rights, so you inherintly lose a trial
% (specifically the first one). Think of it this way, how can you have an
% alternation value on the first trial? Alternation estimates like
% perseveration depend on the previous trajectory.
% IMPORTANT: a "1" indicates an alternation in the traj_change variable. a
% "0" indicates a perseveration

VTE_Persev = [];
for rati = 1:length(rat)    
    % get perveration data
    VTE_Persev{rati,1} = [zIdPhi.(rat{rati}).Baseline.Baseline1(2:end)',traj_change.(rat{rati}).Baseline.Baseline1];
    VTE_Persev{rati,2} = [zIdPhi.(rat{rati}).Baseline.Baseline2(2:end)',traj_change.(rat{rati}).Baseline.Baseline2];
    VTE_Persev{rati,3} = [zIdPhi.(rat{rati}).Saline.Baseline(2:end)',traj_change.(rat{rati}).Saline.Baseline];
    VTE_Persev{rati,4} = [zIdPhi.(rat{rati}).Saline.Saline(2:end)',traj_change.(rat{rati}).Saline.Saline];
    VTE_Persev{rati,5} = [zIdPhi.(rat{rati}).Muscimol.Baseline(2:end)',traj_change.(rat{rati}).Muscimol.Baseline];
    VTE_Persev{rati,6} = [zIdPhi.(rat{rati}).Muscimol.Muscimol(2:end)',traj_change.(rat{rati}).Muscimol.Muscimol];
end

% per rat, are vtes more likely to occur after a perseveration?
numRats = size(VTE_Persev,1);
numCond = size(VTE_Persev,2); % number of conditions on rows

% find zIdPhis for the multiple conditions
zIdPhi_errors = []; zIdPhi_correct = []; zIdPhi_correctAfterError = []; zIdPhi_errorAfterError = [];
for rati = 1:numRats
    
    for condi = 1:numCond
    
        % find alternations
        alt_idx = find(VTE_Persev{rati,condi}(:,2) == 1);

        % find alternations that follow an error
        alt_after_per = alt_idx; % temporary
        alt_after_per(alt_after_per == 1) = []; % if the first trial included is an alternation, we need to remove it (for extraction purposes - not actually removing), so that we can define alt after error
        alt_after_error = alt_after_per((VTE_Persev{rati,condi}(alt_after_per-1,2) == 0));
        
        % get zIdPhi values when errors follow another error
        per_idx = find(VTE_Persev{rati,condi}(:,2) == 0); % find perseverations (errors)    
        per_after_per = per_idx; % temporary
        per_after_per(per_after_per == 1) = []; % need to remove the first trial bc it cant possibly be defined by another perseveration
        error_after_error = per_after_per((VTE_Persev{rati,condi}(per_after_per-1,2) == 0));
        
        % find zIdPhis for the multiple conditions
        zIdPhi_errors{rati,condi}            = VTE_Persev{rati,condi}(per_idx,1);
        zIdPhi_correct{rati,condi}           = VTE_Persev{rati,condi}(alt_idx,1);
        zIdPhi_correctAfterError{rati,condi} = VTE_Persev{rati,condi}(alt_after_error,1);
        zIdPhi_errorAfterError{rati,condi}   = VTE_Persev{rati,condi}(error_after_error,1);

        % find consecutive perseverations?
    
    end
    
end
    

% or, do vte probabilities increase after a couple errors? Do they result
% in perseverations?







