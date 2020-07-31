%% Multi-class classification (leave N out)
% this code uses the svmLIB toolbox to run a multi-class linear SVM
% classifier on your dataset. It does so using a leave 1 out approach.
%
% This classifier leaves N (1 per class) trials out of training and those
% trials are used for testing. This saves considerable amount of
% training/testing time. For example if you have 100 labels, with 10
% classes, a leave 1 out approach would do 100 different testings, 1000
% different times on 1000 different organizations of neurons (numIter*100 iterations). 
% However, this trains on each of the 10 labels in a successive fashion such
% that the model is trained 10 different times, 1000 times (numIter*10 iterations)
%
%
% INPUTS
% data: should be in the following format cell array format:
%           data{1} = bin (so if you were using 7 stem bins, you would have
%           a length(data) = 7. data{n} (where n is a bin) will be of size
%           (row = classes of interest, column = cluster). If you were to
%           select data{1}{1,1} you will have a vector whereby firing rates
%           are in each element, and the firing rates are all in the first
%           column. In other words data{1}{1,1}(:,1) will be your vector of
%           firing rates whereby rows indicate trial number and elements
%           are firing rates. if you were to do data{1}{1,1}(1,:) you
%           should get an error because the data should not be formatted
%           this way.
%
% numClass: a scalar indicating the number of classes. So size(data{1},1)
%           tells you the number of classes if you've correctly formatted
%
% classNames: a cell array with the names of each class
%
% numIter: number of iterations whereby you are drawing random neurons and
%     mixing up the trial ordering. This controls for any spurious effects,
%     however removes any temporal component of the data.
%
% numFeats: number of features (neurons). This will be apparent in your
%           data variable as data{1} columns are the different clusters
%
% numObs: number of observations (trials). Note this should be the smallest
%          number of trials observed. So if one of your neurons only had 12
%          trials, but the rest had 18, and you want to include all data,
%          set numObs = 12
%
% numbins: the number of bins you're analyzing. This is good for analyzing
%           data from multiple stem bins or something. This the primary
%           organization factor for your data cell array, so data{1} is
%           data from say stem bin 1, with the rows of data{1} being the
%           multiple classes of interest (say various things like sample
%           left choice left or something).
%
% ~~~ OUTPUTS ~~~
% probMat: probability matrix
%
% NOTE - This does not work yet, it is in progress to becoming a function
% ALSO NOTE: THIS IS ONLY FOR TJUNCTION RIGHT NOW - NEED TO FIX

%% define these
% the script format implies you have the data variable in your workspace
%clear; clc
%load('data_mazeLocations_multiclass_Format')

numIter    = 1000; % number of iterations
numFeats   = 187;  % number of cells
numObs     = 14;   % number of trials
numbins    = 1;    % number of spatial or time bins
numClasses = 6;    % number of classes to classify

%% Initialization
disp('Checking that inputs are organized appropriately...')
pause(1);

% do some brief checks
checkForm1 = length(data)    == numbins;
checkForm2 = size(data{1},1) == numClasses;
checkForm3 = size(data{1},2) == numFeats;

if checkForm1 ~= 1 || checkForm2 ~= 1 || checkForm3 ~= 1
    disp('Error in data formatting, please open this function and read instructions')
    return
else
    disp('Data organized correctly')
    pause(1);
end

% make an array containing all data (organized by class on rows and
% features on columns)
dataCell = data;

% make labels vector
clear labels
for feati = 1:numClasses
    labels(:,feati) = ones(numObs,1)+(feati-1);
end
labels = labels(:);

% extract data, controlling for sample size
for bini = 1:numbins
    for nIt = 1:numIter 
        clear dataNew dataForm dataPerm

        disp(['shuffle trials ', num2str(nIt)])

        % this randomizes the ordering of trials for each
        % cluster. This, therefore, breaks any temporal
        % component to the population coding
        for classi = 1:numClasses   % loop across classes (rows)
            for clusti = 1:numFeats % loop across features (col)

                % random permutate and draw n observations (trials)
                clear randPull1
                randPull1 = randsample(randperm(size(dataCell{bini}{classi,clusti},1)),numObs);

                % extract data
                dataNew{classi,clusti} = dataCell{bini}{classi,clusti}(randPull1,:);

            end
        end
        % reformat - don't have to loop separately. This new
        % format will be an array where each cell element is an
        % observation, and within the first shell is the number
        % of neurons concatenated horizontally
        for m = 1:numbins       % loop across number of bins
            for classi = 1:numClasses
                for mm = 1:numFeats % loop across number of feats
                    dataForm{classi,m}(:,mm) = dataNew{classi,mm}(:,m);
                end
            end
        end

        % concatenate the data vertically so that one class is
        % on top, one class is on bottom                    
        dataPerm = [];
        % generalized formatting for classifier
        for binii = 1:numbins % loop across bins
             clear temp
             % isolate one bin
             temp = dataForm(:,binii);

             % concatenate vertically
             dataPerm{binii} = vertcat(temp{:});
        end

        % replace any NaNs with 0 - note that NaNs will drive
        % the accuracy of your model to 50%
        for binii = 1:numbins
            dataPerm{binii}(find(isnan(dataPerm{binii})==1)) = 0;
        end

        addpath('X:\03. Lab Procedures and Protocols\MATLABToolbox\chronux\spectral_analysis\continuous\libsvm-3.20\matlab')           

        % need to make this a loop across bins
        clear performance predict_label accuracy p

        for nLab = 1:numel(labels)/numClasses % this needs to change, it needs to be a permutative test, so prob 1000

            % need to test one of each class
            clear rand_label
            for classi = 1:numClasses
                rand_label(classi) = randperm(numObs,1); % get a random number between 1 and 6
            end

            % use labels to grab training data
            clear idx_temp idx_testing
            for classi = 1:numClasses
                idx_tmp{classi}     = find(labels(:)==classi);
                idx_testing{classi} = idx_tmp{classi}(rand_label(classi));
            end
            idx_testing = cell2mat(idx_testing)';            

            % check that the testing index is formatted correctly
            for classi = 1:numClasses
                if labels(idx_testing(classi)) ~= classi
                    disp('warning - youve erroneously selected training/testing data')
                    break
                end
            end  
            
            % clear important variables
            clear trainData testData trainLabel testLabel

            % ~~~~~~~~~~~~~~~~~~~~~~~~~ CHANGE
            % training data
            trainData          = dataPerm{bini}; % CHANGE ME TO DYNAMIC
            trainLabel         = labels;
            trainData(idx_testing,:)  = [];
            trainLabel(idx_testing,:) = [];

            % train model
            clear model
            model = svmtrain(trainLabel, trainData, '-c 1 -t 0 -b 1');

            % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CHANGE
            % testing data - need one observation per class
            testData  = dataPerm{bini}(idx_testing,:); % IMPORTANT - CHANGE TO DYNAMIC
            testLabel = labels(idx_testing,:);

            % test classifier
            [predict_label(nLab,:), accuracy, p{nLab}] = svmpredict(testLabel, testData, model,'-b 1');

            % store accuracy
            performance(nLab,:) = accuracy(1);

        end

        % averaged accuracy
        meanPerf{nIt} = mean(performance);

        % probability of predicting classes
        classMat = [1:numClasses]'; % variable indicating different classes

        % make a matrix that contains probabilities of hits and misses, where
        % mat(1,1) = 1 predicted as 1, mat(1,2) = 1 predicted as 2, mat(2,1) =
        % 2 predicted as 1, mat(2,2) = 2 predicted as 2
        for i = 1:numClasses
            clear idxClass trueClass predClass idxHit numTest idx
            % define true and predicted classes
            idxClass  = find(labels == classMat(i)); % an index of the first group of classes
            trueClass = labels(idxClass);            % true classes        
            predClass = predict_label(idxClass);     % predicted classes

            % total number of labels per class
            numTest = length(trueClass); % total number of cases where 1 should have been

            % fill in matrix
            for ii = 1:numClasses
                % predClass is defined once above and limited to one class
                idx = find(predClass == classMat(ii));
                probMat{bini}{nIt}(i,ii) = length(idx)/numTest;
            end
        end
    end   
    % get probability matrix
    prob3d{bini}  = cat(3,probMat{bini}{:});
    probAvg{bini} = mean(prob3d{bini},3);
    probStd{bini} = std(prob3d{bini},[],3);

    % shuffling labels to create a chance level distribution
    for nShuff = 1:numIter
        disp(['shuffle labels # ',num2str(nShuff)]);

        % shuffle labels
        shuffIdx = randperm(length(labels));
        labelShuff = labels(shuffIdx);

        for nLab = 1:numel(labelShuff)/numClasses % this needs to change, it needs to be a permutative test, so prob 1000

            % need to test one of each class
            clear rand_label
            for classi = 1:numClasses
                rand_label(classi) = randperm(numObs,1); % get a random number between 1 and 6
            end

            % use labels to grab training data
            clear idx_temp idx_testing
            for classi = 1:numClasses
                idx_tmp{classi}     = find(labels(:)==classi);
                idx_testing{classi} = idx_tmp{classi}(rand_label(classi));
            end
            idx_testing = cell2mat(idx_testing)';            

            % clear important variables
            clear trainData testData trainLabel testLabel

            % training data
            trainData                 = dataPerm{bini};
            trainLabel                = labelShuff;
            trainData(idx_testing,:)  = [];
            trainLabel(idx_testing,:) = [];

            % train model
            clear model
            model = svmtrain(trainLabel, trainData, '-c 1 -t 0 -b 1');

            % testing data - need one observation per class
            testData  = dataPerm{bini}(idx_testing,:); 
            testLabel = labels(idx_testing,:);

            % test classifier
            clear accuracy
            [predict_label(nLab,:), accuracy, p{nLab}] = svmpredict(testLabel, testData, model,'-b 1');

            % store accuracy
            performanceShuff(nLab,:) = accuracy(1);

        end

        % averaged accuracy
        meanShuff{nShuff} = mean(performanceShuff);

        % probability of predicting classes
        classMat = [1:numClasses]'; % variable indicating different classes

        % make a matrix that contains probabilities of hits and misses, where
        % mat(1,1) = 1 predicted as 1, mat(1,2) = 1 predicted as 2, mat(2,1) =
        % 2 predicted as 1, mat(2,2) = 2 predicted as 2
        for i = 1:numClasses
            clear idxClass trueClass predClass idxHit numTest idx
            % define true and predicted classes
            idxClass  = find(labelShuff == classMat(i)); % an index of the first group of classes
            trueClass = labelShuff(idxClass);            % true classes        
            predClass = predict_label(idxClass);     % predicted classes

            % total number of labels per class
            numTest = length(trueClass); % total number of cases where 1 should have been

            % fill in matrix
            for ii = 1:numClasses
                % predClass is defined once above and limited to one class
                idx = find(predClass == classMat(ii));
                probShuff{nShuff}(i,ii) = length(idx)/numTest;
            end
        end    
    end
 
    prob3dShuff{bini}  = cat(3,probShuff{:});
    probAvgShuff{bini} = mean(prob3dShuff{bini},3);
    probStdShuff{bini} = std(prob3dShuff{bini},[],3);

    % statistics
    for i = 1:numClasses
        for ii = 1:numClasses
            [hZ{bini}(i,ii),pZ{bini}(i,ii),ci,statZ{bini}{i,ii}] = ztest(probAvg{bini}(i,ii),probAvgShuff{bini}(i,ii),probStdShuff{bini}(i,ii));
        end
    end
end

% figure
binSelect = 1;
figure('color','w')
subplot 211
    imagesc(probAvg{binSelect})
    c = colorbar;
    shading 'interp'
    title('Multi-class performance')
    ylabel(c,'Prob. of prediction')
    ax = gca;
    ax.XTick      = 1:numClasses;
    ax.YTick      = 1:numClasses;    
    ax.XTickLabel = classNames;
    ax.YTickLabel = classNames;
    ax.XTickLabelRotation = 45;
subplot 212
    imagesc(probAvgShuff{binSelect})
    c = colorbar;
    shading 'interp'
    title('Multi-class shuffle')
    ylabel(c,'Prob. of prediction')
    ax = gca;
    ax.XTick = 1:numClasses;
    ax.YTick = 1:numClasses;    
    ax.XTickLabel = classNames;
    ax.YTickLabel = classNames; 
    ax.XTickLabelRotation = 45;

% plot as a line graph
probAvgMazeLocs   = diag(probAvg{binSelect});
probStdMazeLocs   = diag(probStd{binSelect});
probAvgShuffDiag  = diag(probAvgShuff{binSelect});
probAvgShuffStdDi = diag(probStdShuff{binSelect});

figure('color','w')
errorbar(probAvgMazeLocs,probStdMazeLocs,'k','LineWidth',2);
hold on;
errorbar(probAvgShuffDiag,probAvgShuffStdDi,'Color',[0.8 0.8 0.8],'LineWidth',2)
ylabel('Probability of accurate prediction')
box off
ax = gca;
ax.XTick = 1:numClasses;
ax.XTickLabel = classNames;
ax.XTickLabelRotation = 45;

    
