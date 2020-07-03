%% this will be a function for multi class classification
%
% INPUTS
% data1 and data2: should be cell arrays containing your data. However, note
%     that each row should contain separate data. Each column will be a
%     separate feature.
%
%
% NOTE - This does not work yet

%% things to define as inputs
% predefine in function inputs*******
numClasses = 4;

% define number of observations
numObs = 6; % smallest number of trials (observations)

% define number of features
numFeats = 187; % will be the number of features (i.e. dimensions, cells, etc..)

% define number of iterations
numIter = 1000;

% will be defined in funtion
data1 = vertcat(sampleLvR.lefts, sampleLvR.rights);
data2 = vertcat(choiceLvR.lefts, choiceLvR.rights);

%% stuff inside the function

% make an array containing all data (organized by class on rows and
% features on columns)
dataCell = vertcat(data1,data2);

% make labels vector
clear labels
for feati = 1:numClasses
    labels(:,feati) = ones(numObs,1)+(feati-1);
end
labels = labels(:);

% extract data, controlling for sample size
for nIt = 1:numIter 
    disp(['shuffle trials ', num2str(nIt)])

    % this randomizes the ordering of trials for each
    % cluster. This, therefore, breaks any temporal
    % component to the population coding
    for classi = 1:numClasses   % loop across classes (rows)
        for clusti = 1:numFeats % loop across features (col)

            % random permutate and draw n observations (trials)
            clear randPull
            randPull1 = randsample(randperm(size(dataCell{classi,clusti},1)),numObs);

            % extract data
            dataNew{classi,clusti} = dataCell{classi,clusti}(randPull1,:);

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
    for bini = 1:numbins % loop across bins
         clear temp
         % isolate one bin
         temp = dataForm(:,bini);
         
         % concatenate vertically
         dataPerm{bini} = vertcat(temp{:});
    end

    % replace any NaNs with 0 - note that NaNs will drive
    % the accuracy of your model to 50%
    for bini = 1:numbins
        dataPerm{bini}(find(isnan(dataPerm{bini})==1)) = 0;
    end

    addpath('X:\03. Lab Procedures and Protocols\MATLABToolbox\chronux\spectral_analysis\continuous\libsvm-3.20\matlab')           

    % need to make this a loop across bins
    clear predict_label accuracy p
    
    % ~~~~~~~~~ THIS DOES NOT WORK RIGHT YET - ALSO, NOT DYNAMIC - ALSO,
    % DOESNT WORK WITH large scale iteration
    
    for nIt2 = 1:(numIter/2) % this needs to change, it needs to be a permutative test, so prob 1000

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
        
        % check
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
        trainData              = dataPerm{7}; % CHANGE ME TO DYNAMIC
        trainLabel             = labels;
        trainData(idx_testing,:)  = [];
        trainLabel(idx_testing,:) = [];
        
        % final check for correct indexing
        clear check
        for classi = 1:numClasses
             check(classi) = length(find(trainLabel == classi));
        end
        check = unique(check);
        
        % if there are more than one lengths of training labels, that means
        % there was incorrect indexing of labels. I pulled 1 of each label
        % for testing
        if length(check) > 1
            disp('error in indexing labels')
            return
        end
        
        % train model
        model{nIt}{nIt2} = svmtrain(trainLabel, trainData, '-c 1 -t 0 -b 1');

        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CHANGE
        % testing data - need one observation per class
        testData  = dataPerm{7}(idx_testing,:); % IMPORTANT - CHANGE TO DYNAMIC
        testLabel = labels(idx_testing,:);

        % test classifier
        [predict_label{nIt}{nIt2}, accuracy{nIt}{nIt2}, p{nIt}{nIt2}] = svmpredict(testLabel, testData, model{nIt}{nIt2},'-b 1');
    
        % get probabilities
        %prob(:,k) = p(:,model{nIt}{nIt2}.Label==1);    %# probability of class==k

    end
    
    % ~~~ now we need to average across classification selections ~~~ %
    
    % make 3D array
    probMat{nIt} = cat(3,p{nIt}{:});
    
    % take average across classification selections
    probAvg{nIt} = mean(probMat{nIt},3);
    
    % after this, we'll need to average across trial selections   
end   
    % make 3D array
    probItMat = cat(3,probAvg{:});
    
    % take average across classification selections
    probItAvg = mean(probItMat,3);
    
    figure('color','w')
    imagesc(probItAvg)
    colorbar
    


%% other ideas
% get accuracies
accuracy = horzcat(accuracy{:});
ClassifierPerf = accuracy(1,:);

% find match/ mismatch in labels
predictedLabels = cell2mat(predict_label)';
totalAcc = (length(find(predictedLabels == labels))/(numel(labels)))*100;

% which labels were often mixed up?
combos = combntns(1:numClasses,2);
numCombos = size(combos,1);

for comboi = 1:numCombos
    predictionErrors{comboi} = find(labels == combos(comboi,1) & predictedLabels == combos(comboi,2));
end
    
    
    % loop iteratively
    for nn = 1:1000
        disp(['shuffle ratio of training/testing ',num2str(nn)]);
        for bini = 1:numbins
            % randperm
            perm_cond1 = randperm(numObs); % use the number of observations per condition to pull randomly from
            perm_cond2 = randperm(numObs); 
            % get labels
            labels = vertcat(ones(numObs,1),-ones(numObs,1)); 
            % training and testing labels
            floor_75 = floor(numObs*.75);
            ceil_end = ceil(numObs*.75);
            train_cond1 = perm_cond1(1:floor_75);
            train_cond2 = perm_cond2(1:floor_75); % take 1 through the floor of 75% of the data
            test_cond1  = perm_cond1(ceil_end:end);
            test_cond2  = perm_cond2(ceil_end:end); % take ceiling of 75% of data till the end
            % data
            data_cond1 = data1{bini}(perm_cond1,:); % get random data (not really random, but random trials)
            data_cond2 = data2{bini}(perm_cond2,:);
            % training data
            data_training(1:floor_75,:) = data_cond1(train_cond1,:);
            data_training(ceil_end:floor_75*2,:) = data_cond2(train_cond2,:);
            % testing data
            data_testing(1:length(ceil_end:numObs),:) = data_cond1(test_cond1,:);
            data_testing(length(ceil_end:numObs)+1:length(ceil_end:numObs)*2,:) = data_cond2(test_cond2,:);
            % labels
            labels_training = vertcat(repmat(1,floor_75,1),repmat(-1,floor_75,1));
            labels_testing  = vertcat(repmat(1,length(ceil_end:numObs),1),repmat(-1,length(ceil_end:numObs),1));
            % svm training
            model = svmtrain(labels_training, data_training, '-c 1 -t 0');
            [predict_label, accuracy, dec_value] = svmpredict(labels_testing, data_testing, model);
            % store accuracy - row is shuffled trials, column
            % is shuffled ratio of training/testing
            Accuracy_iteration{bini}(n,nn) = accuracy(1);
            % clear data
            clear labels_training labels_testing data_testing data_training data_sample ...
                data_choice train_sample test_sample train_choice test_choice perm_sample ...
                perm_choice labels_sample labels_choice model predict_label accuracy dec_value
        end
    end
end
 

