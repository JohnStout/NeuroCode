%% 
% this function takes a variable, normalizes it between 0 and 1, then sorts
% it based on a separate input. The output is a graph that is a heatplot,
% but sorted so that the max rate of each neuron is organized based on the
% column, or x-axis
%
% INPUTS
% x: a matrix that you want sorted. Should be formated with rows being cell
% number and columns being observation
% y: a matrix that you want to sort by. Same formatting as above
%
% ~~~ VERY IMPORTANT ~~~
% This analysis is extremely misleading unless you average population vectors across trials
% like the buzsaki paper on time cells did. If you normalize each population vector to its
% peak rate within the epoch you're interested in examining, you will observe a beautiful diagonal
% band. But this diagonal band is more likely an artifact of probability, that what the brain is doing.
% This code forces a vector to be between 0 and 1. Therefore, future iterations of this code will require
% a matrix for each individual neuron, then do some sort of averaging to make this better reflect true
% neuronal activity and not artificially determined relevant activity.
%
% OUTPUTS:
% figure
%
% writen by John Stout

function [] = SortedRateMap(x,y)

% normalize each row between 0 and 1
numcells = size(x,1);
numbins  = size(x,2);

for celli = 1:numcells
    for bini = 1:numbins
        x_norm(celli,bini) = (x(celli,bini)-min(x(celli,:)))/...
            (max(x(celli,:))-min(x(celli,:)));   
        y_norm(celli,bini) = (y(celli,bini)-min(y(celli,:)))/...
            (max(y(celli,:))-min(y(celli,:)));           
    end
end

% get indices to sort by
[maxval, idxmax]  = max(y_norm');
[matsort,idxsort] = sort(idxmax);

% use idxsort to sort the x_norm variable
x_sort = x_norm(idxsort,:);

% make figure
figure('color','w');
imagesc(x_sort);
colorbar
ylabel('Neuron ID')
xlabel('Bin Number')
set(gca,'FontSize',13);
colormap('jet')

end
