%% Shannons entropy
% prob for bin 1 = (bin count/sum of all bin count)
% log_val for bin 1 = prob*logbase2
% probability-log-probability for bin 1 = prob*log_val
% sum the probability-log-probability values and *-1
%
% Variable formatting: 'binned_spks' is a cell array containing doubles for
%                           each element. The doubles represent a spk count
%                           across the entire session for a given bin. For
%                           example if there are 7 bins and binned_spks{1}
%                           = [21 1 17 18 9 8 1] then 21 represents 21
%                           total spks across the session for bin 1.
%
%                      'clusters' is a Nx1 struct containing the names of
%                      all clusters in the session                     
%
% written by John Stout
% last edit 12/16/18

function [shannons_entropy] = shannon_entropy(binned_spks,clusters)

    for ci = 1:length(clusters)

        % convert your data to probability values
        prob{ci} = binned_spks{ci}./sum(binned_spks{ci});   

        % take the log2 for bits. add eps to prevent taking log2 of zero
        log_var{ci} = log2(prob{ci}+eps);

        % multiply the log-prob by the prob
        log_prob{ci} = prob{ci}.*log_var{ci};

        % sum the log-probability values
        shannons_entropy{ci} = (sum(log_prob{ci}))*-1;

        clear hdat xval
    end

    % convert to double
    shannons_entropy = cell2mat(shannons_entropy);

    % get rid of cells that have NaNs (didn't fire)
    shannons_entropy(isnan(shannons_entropy))=0;

end