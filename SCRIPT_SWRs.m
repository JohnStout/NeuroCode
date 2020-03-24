%% get_SWRs
%
% this function is used to extract sharp-wave ripple counts
%
% INPUTS: 
% datafolder: a string variable used to define directory
% input: a structure array containing directions for which if else
%        statements to run
% Int: a matrix containing timestamp, accuracy, and trajectory information
% lfp: a 1xN vector
% Timestamps: 1xN vector; same length of lfp 
% phase_bandpass: a vector that tells the function which phases to extract
%                 ripples between
% srate: sampling rate of LFP
%
% OUTPUTS:
% swr_count: a variable containing the number of sharp-wave ripples
% swr_events: contains lfp data from the ripple events
% swr_event_idx: the index of events
%
% Future iterations of this function should include controlling for speed.
% instantaneous_speed_DNMP was written as a way to do so, but it hasn't
% been controlled for here.
%
% written by John Stout and Andrew Garcia
    
%% identify and isolate SWRs



%% preperatory things
    
    datafolder = 'X:\01.Experiments\FASD\FASD_LFP_1913 111519 DA SWP';
    load(strcat(datafolder,'\Events.mat'))

    % get lfp data
    addpath('X:\03. Lab Procedures and Protocols\MATLABToolbox\chronux\spectral_analysis\continuous');
    hpc_lfp = load(strcat(datafolder,'\CSC3'));
    [Timestamps, lfp] = interp_TS_to_CSC_length_non_linspaced(hpc_lfp.Timestamps, hpc_lfp.Samples);     

    % define epochs of interest
    pre_sleep_start  = TimeStamps_EV(2);
    pre_sleep_end    = TimeStamps_EV(3);
    post_sleep_start = TimeStamps_EV(7);
    post_sleep_end   = TimeStamps_EV(8);  
    
    % define the bandpass filter for swrs
    phase_bandpass = [150 250];  

    % get timestamps from between pre_sleep_start to pre_sleep_end. Do the same
    % for post_sleep_start and post_sleep_end. Save as separate variables.
    lfp_pre  = lfp(Timestamps>pre_sleep_start & Timestamps<pre_sleep_end);
    %? %lfp_post = lfp(Timestamps>post_sleep_start & Timestamps<post_sleep_end);
    
    %figure; 
    %plot(lfp_pre)
    
    %figure; 
    %plot(lfp_pre)
    
    % clean lfp?
    input.clean = 1; 
    
    if input.clean == 1 
        addpath('X:\03. Lab Procedures and Protocols\MATLABToolbox\chronux\spectral_analysis\continuous')
        addpath('X:\03. Lab Procedures and Protocols\MATLABToolbox\Basic Functions')
        % define parameters for chronux toolbox functions
        params.tapers   = [5 9]; %[NW k] [timebandwidth product leading tapers] but where do you find these parameters {s}
        params.pad      = 0;
        params.Fs       = hpc_lfp.SampleFrequencies(1,1);
        params.fpass    = [0 300]; %[0 Fs/2] why300? {s}
        params.err      = [2 .05];
        params.trialave = 0;
    %lfp_data_pre    = DetrendDenoise(lfp_pre,params.Fs); % detrend
    %because im skipping detrend denoise 
        lfp_data_pre = lfp_pre
        lfp_data_pre    = lfp_data_pre'; %why do you need the data to be horiz?
    else
        lfp_data_pre = lfp_pre; 
    end

    %define sampling rate of lfp (is this synonymous to time (2) would be 2s?)
    srate = hpc_lfp.SampleFrequencies(1);
    
    figure('color','white'); hold on;
    plot(lfp_data_pre(1:2000),'black'); % 1 second if 1:2000
    xlabel('Samples')
    ylabel('LFP Raw')
    box off
    
%% filtering data
        
    % bandpass filter (3rd deg butterworth filter)
    [lfp_filtered] = skaggs_filter_var(lfp_data_pre, phase_bandpass(1),...
        phase_bandpass(2), srate);  
    
    figure('color','w'); hold on;
    subplot 311
        plot(lfp_data_pre(1:2000),'k');
        ylabel('lfp raw data')
    subplot 312
        plot(lfp_filtered(1:2000),'k');
        ylabel('lfp filtered')
    subplot 313
        plot(lfp_hilbert(1:2000),'r');
        hold on;
        plot(lfp_filtered(1:2000),'b');
        ylabel('LFP hilbert')
    % Hilbert Transform, abs is absolute because hilbert is complex vals
    lfp_hilbert = abs(hilbert(lfp_filtered));
    
    %figure('color','w'); hold on;
    %subplot 313
        %plot(lfp_hilbert(1:2000),'r')
        %ylabel('LFP hilbert')

    % smooth with gaussian
    input.smooth_w_gauss = 1;
    if input.smooth_w_gauss == 1

        % generate gaussian (will be a 64 pnt gaussian window)
        cd 'X:\03. Lab Procedures and Protocols\MATLABToolbox\Toolboxes\signal\signal'
        gauss_width = (srate/1000)*32; % 32ms * samples per ms = 64 samples = 32 ms of data
        alpha       = 4; % std = (N)/(alpha*2) -> https://www.mathworks.com/help/signal/ref/gausswin.html
        w           = gausswin(gauss_width,alpha);   

        % convolve the gaussian kernel with lfp data
        lfp_smoothed = conv(lfp_hilbert,w,'same');

        % save as common name variable for smooth running
        lfp_data_pre = lfp_smoothed;

        % display
        disp('Smoothed with Gaussian Filter')

    elseif input.smooth_w_gauss == 0

       lfp_data_pre = lfp_hilbert;

       % display
       disp('Did not smooth with gaussian filter')

    end

%
figure('color','w'); hold on;
subplot 311;
plot(lfp_filtered(1:4000),'r');
ylabel('LFP BANDPASS');
subplot 312;
plot(lfp_hilbert(1:4000),'b');
ylabel('LFP HILBERT')
subplot 313;
plot(lfp_smoothed(1:4000),'b');
ylabel('LFP GAUSSIAN')
%
    % how many standard deviations above the mean does Jadhav use?
    std_above_mean = 4; %JADHAV AND OTHERS USE 4-6 FOR MORE THAN 15MS
    
    % zscore
    input.zscore = 1;
    if input.zscore == 1

        lfp_z    = zscore(lfp_data);
        lfp_data = [];
        lfp_data = lfp_z;

        % find above threshold
        aboveThreshold = (lfp_data >= std_above_mean);  
        aboveThreshold = double(aboveThreshold); 

        % display
        disp('zscored')

    else
        % find above threshold
        % http://www.jneurosci.org/content/jneuro/37/49/11789.full.pdf
        aboveThreshold = (lfp_data >= ((mean(lfp_data)+...
            ((std(lfp_data))*std_above_mean))));  
        aboveThreshold = double(aboveThreshold); 

        % display
        disp('not zscored')
    end  
    
    figure(); plot(lfp_z(1:4000))

    clear edges rising fallilng spanWidth wideEnough startPos endPos...
        swr_event_string idx swr_events  zerotimes ...
        swr_onset_temp swr_offset_temp

    % gives the loc of the "0" to 1 transition, 0 not the 1, or vice versa, 
    % will be 1 or -1
    edges       = diff(aboveThreshold);
    rising      = find(edges == 1);
    falling     = find(edges == -1);
    spanWidth   = falling - rising;

    % find when the swr meets the criteria of being >= 15ms
    wideEnough  = find(spanWidth >= 30); 

    % index the start and end of swrs that meet the wideEnough criteria
    startPos = []; endPos = [];
    startPos    = rising(wideEnough);                       
    endPos      = falling(wideEnough) - 1;    

    % concatenate indexed swr times
    swr_event_string = cell2mat(arrayfun(@(x,y) x:1:y, ...
        startPos', endPos', 'uni', false));  

    % find lfp data for each ripple
    for j = 1:length(startPos)
        swr_events{j} = lfp_data(startPos(j):endPos(j));
        swr_event_index{j} = startPos(j):endPos(j);
    end

    % find how many swrs
    swr_count = size(swr_event_index,2); 
    
    
    
    
    
    % are we sure that each event is a swr? Is it possible that we're also
    % extracting noise?
    
    % plot out SWR events

