%% hilbert_phase
%
% INPUTS: signal_filtered - this is the filtered signal using some bandpass
% filtering method (I used 3rd degree butterworth aka skaggs filter)
%
% OUTPUTS: inst_phase and inst_power - both instantaneous measures of phase
%          and power. inst_phase is in radians but centered around zero.
%
%          PhaseRadians is also in radians, however, it's not zero
%          centered. Thus when multiplying by 180/pi you should get a
%          circle that ranges from 0 to 360
% 
% Written by John Stout

function [phase_corrected,inst_phase,inst_power]=hilbert_metrics(signal_filtered)

% Hilbert Transform 
lfp_hilbert = (hilbert(signal_filtered)); 

% hilbert angle (phase)
inst_phase = angle(lfp_hilbert);

% hilbert power (length squared)
inst_power = (abs(lfp_hilbert)).^2; 

%% convert radians
% since radians will be centered around zero, you'll need to convert them
% or you'll have a max phase degree of 180 and a min of -180.
    % find the smallest inst_phase radian
        min_phase = min(inst_phase);
    % add the magnitude of the smallest value to the entire vector. This
    % should scale everything so that the smallest value is now zero. Now
    % when you convert to degrees the largest will be roughly 360, smallest
    % is 0.
    phase_corrected = inst_phase+abs(min_phase);

%% plot to check
%{
figure();
subplot 411
plot(inst_power(1:2000),'k');
subplot 412
plot(inst_phase(1:2000),'m');
subplot 413
plot(phase_corrected(1:2000),'r');
subplot 414
plot(phase_corrected(1:2000)*180/pi,'b');
%}

end