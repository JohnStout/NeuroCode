%% Freedman-Dianconis
%
% This function utilizes a statistical estimate for the number of bins to
% use for calculations of entropy
%
% Written by John Stout
% last edit on 12/11/18

%% fun
function [nbins] = freedman_diaconis(Int,ExtractedY,TimeStamps_VT,mazei)
% index timestamps during stem to get Y values    
    for i = 1:size(Int,1);
        time_idx{i} = find(TimeStamps_VT>Int(i,mazei(1))...
            & TimeStamps_VT<Int(i,mazei(2)));
    end
    space_idx = horzcat(time_idx{:});

% index y values
  x = ExtractedY(space_idx);
            
% Freedman-Diaconis rule to statistically estimate number of bins
  nbins = ceil((max(x)-min(x))/(2*(iqr(x))*(length(x)^(-1/3))));

end