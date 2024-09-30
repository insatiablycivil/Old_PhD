%% Main

function AnomalyFeat01
    global diagMode
    global verbose
    global testID
    global placeholder
   
    %%% Run Time Settings %%%
    diagMode = false; % True to plot stuff
    verbose = false; % True to plot more stuff
%     testID = 65;
%     placeholder(2) = placeholder(2) + placeholder(3);
%     for i = placeholder(2) + placeholder(3) - 1 : - 1: placeholder(2)
%placeholder = placeholder + 1;
%  PlotTrace(7,1)
%     end
%     PlotTrace(969)
%                1  2  3  4  5  6  7  8  9 10 
%       steps = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
        steps = [0, 0, 0, 0, 0, 0, 0, 1, 0, 0];
%       steps = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
%     step0 = false;
%     step1 = false;
%     step2 = false;
%     step3 = false;
%     step4 = false;
%     step5 = false;
%     step6 = false;
%     step7 = false;
%     step8 = false;
%     step9 = false;
    %%% Run Time Settings %%%
    
    if steps(1) == true
        Import % Import Data
    end
    if steps(2) == true
        PreProcess % Not Used
    end
    if steps(3) == true
        Process % Find Start and End Times
    end    
    if steps(4) == true
        FeatureExtraction; % Feature Extraction
    end
    if steps(5) == true
        FeatureSelection % Feature Selection
    end
    if steps(6) == true
        ModelTraining
    end
    if steps(7) == true
        ModelApplying
    end
    if steps(8) == true
        ViewResults
    end
    if steps(9) == true
        Diagnose
    end
    if steps(10) == true
        DiagnoseTesting
    end
end

%% Import

function Import
    % Imported Data
    global data
    %global codes 
    global numRecords 
    global maxLength
    
    % Imports Data and Codes
    ImportNPG;    
    
    % Sets Maximum Lengths and Number of Records
    numRecords = length(data);
    maxLength = 0;
    for i = 1 : numRecords
       maxLength = max([maxLength length(data{i})]); 
    end
end

%% Import DataNPG

function ImportNPG
    % Imported Data
    global model
    global sourceCurrent
    global sourceVoltage
    global data
    global codes
    global labels
    
    % Taken from WorkSpace
    model = evalin('base','A7T');
    sourceCurrent = evalin('base','allNPG');
    sourceVoltage = evalin('base','allNPGVoltage');
    kelmanFeatures = evalin('base','kelman_features');
    labels = evalin('base','A7TFeatLabels');
    % For now, just use first column
    labels = labels(:,1);
    
    % Initialise
    data = cell(length(model),1);
    codes = zeros(length(data), 9);
    
    % Reads and Asigns Values
    for i = 1 : length(model)
        ID = i;
        recordID = model(ID);
        
        for j = 1 : 2
            if j == 1
                record = sourceCurrent{1,recordID};
            else
                record = sourceVoltage{1,recordID};
            end
            record(end + 1 : 13000) = 0;
            record(isnan(record)) = 0;

            recordMax = record(1);
            recordScale = record(2);
            recordFreq = record(3);
            recordDelay = record(4);
            
            if j == 1
                recordStart = (recordDelay + 4)*1000/recordFreq;
                recordLatch = kelmanFeatures(recordID,1);
                recordBuffer = kelmanFeatures(recordID,2);
                recordMCon = kelmanFeatures(recordID,3);
                recordACon = kelmanFeatures(recordID,4);
                recordEndTime = kelmanFeatures(recordID,5);

                recordpk1 = (kelmanFeatures(recordID,6));
                recordplt = (kelmanFeatures(recordID,7));
                
                codes(i,1)  = recordStart;
                codes(i,2)  = recordLatch;
                codes(i,3)  = recordBuffer;
                codes(i,4)  = recordACon;
                codes(i,5)  = recordEndTime;
                codes(i,6)  = recordMCon;
                codes(i,7)  = recordpk1;
                codes(i,8)  = recordplt;
                codes(i,9)  = recordMax;
                codes(i,10) = recordScale;
                codes(i,11) = recordFreq;
                codes(i,12) = recordDelay;
                codes(i,13) = 0;%allCodes(i);
            end

            record = (record/recordScale)*recordMax;
            record(1) = recordMax;
            record(2) = recordScale;
            record(3) = recordFreq;
            record(4) = recordDelay;

            data{i,j} = record;
        end
    end
end

%% PreProcess

function PreProcess
    % For Now, No PreProcessing    
    % Imported Data
    global data
    global recordList
    global numRecords
    global maxLength
    global codes
    global processedCodes
    
    recordList = zeros(numRecords,maxLength,2);
    processedCodes = codes;
    
    for i = 1 : numRecords
        if i == 217
            disp('hi')
        end
        recordList(i,:,1) = data{i,1};
        recordList(i,:,2) = data{i,2};
        % Check Current Max, Scale, and Frequency: Use Mode as Backup
        for j = 9 : 11
            if isnan(codes(i,j)) || isinf(codes(i,j)) || (codes(i,j) < 1)
                 processedCodes(i,j) = mode(recordList(:,j-8,1)); 
            end
        end
        % Set Min Time (equiv. 5)
        minTime = 5 * 1000 / processedCodes(i,11);
        % Check Delay: Use Min. Time as Backup
        if isnan(codes(i,12)) || isinf(codes(i,12)) || (codes(i,12) < minTime)
            processedCodes(i,12) = minTime;
        else
            processedCodes(i,12) = processedCodes(i,12) + minTime; % HERE
        end
        % Set Start: Use Delay (Min. Time) as Backup
        if isnan(processedCodes(i,1)) || isinf(processedCodes(i,1)) || ...
                (processedCodes(i,1) < minTime) || ...
                (processedCodes(i,1) > ((maxLength/processedCodes(i,11))*1000)-1)
            processedCodes(i,1) = processedCodes(i,12);
        end  
        % Set Max Time (Length2Time(MaxLength) - Start)
        maxTime = ((maxLength / processedCodes(i,11)) * 1000) - processedCodes(i,1);
        % Set End: Use Max Length as Backup
        if isnan(processedCodes(i,5)) || isinf(processedCodes(i,5)) || ...
                (processedCodes(i,5) < minTime + 1) || ...
                (processedCodes(i,5) > maxTime)
            processedCodes(i,5) = maxTime + processedCodes(i,1);
        else
            processedCodes(i,5) = processedCodes(i,5) + processedCodes(i,1);
        end
        % Set Latch: Use 0 as Flag
        if isnan(processedCodes(i,2)) || isinf(processedCodes(i,2)) || ...
                (processedCodes(i,2) < 1) || ...
                (processedCodes(i,2) >= processedCodes(i,5) - 1)
            processedCodes(i,2) = 0;
        end
        % Set ACon: Use 0 as Flag
        if isnan(processedCodes(i,4)) || isinf(processedCodes(i,4)) || ...
                (processedCodes(i,4) < 1) || ...
                (processedCodes(i,4) < processedCodes(i,2)) ||...
                (processedCodes(i,4) >= processedCodes(i,5) - 1)
            processedCodes(i,4) = 0;
        end
        % Set Buffer: Use 0 as Flag
        if isnan(processedCodes(i,3)) || isinf(processedCodes(i,3)) || ...
                (processedCodes(i,3) < 1) || ...
                (processedCodes(i,3) < processedCodes(i,2)) ||...
                (processedCodes(i,3) >= processedCodes(i,4))
            processedCodes(i,3) = 0;
        end
        % Spread Evenly for Invalid Values
        if processedCodes(i,2) == 0
        % Blank Latch
            if processedCodes(i,3) == 0
            % Blank Latch + Buffer
                if processedCodes(i,4) == 0
                % Blank Latch + Buffer + ACon
                    % Use End: ((End - Start) / 4) + Start
                    for j = 1 : 4
                        processedCodes(i,j+1) = (j * (processedCodes(i,5) - processedCodes(i,1)) / 4) + processedCodes(i,1);
                    end
                else
                % Working Acon
                    % Use ACon: (ACon / 3) + Start
                    for j = 1 : 3
                        processedCodes(i,j+1) = (j * processedCodes(i,4) / 3) + processedCodes(i,1);
                    end
                end
            else
            % Working Buffer
                % Use Buffer: (Buffer / 2) + Start
                for j = 1 : 2
                    processedCodes(i,j+1) = (j * processedCodes(i,3) / 2) + processedCodes(i,1);
                end
                if processedCodes(i,4) == 0
                % Blank ACon
                    % Use End: ((End - Buffer) / 2) + Buffer + Start
                    for j = 1 : 2
                        processedCodes(i,j+3) = (j * (processedCodes(i,5) - processedCodes(i,3)) / 2) + processedCodes(i,3) + processedCodes(i,1);
                    end
                else
                % Working ACon
                    % Use ACon: ACon + Start
                	processedCodes(i,4) = processedCodes(i,4) + processedCodes(i,1);
                end
            end
        else
        % Working Latch
            % Use Latch: Latch + Start
            processedCodes(i,2) = processedCodes(i,2) + processedCodes(i,1);
            if processedCodes(i,3) == 0
            % Blank Buffer
                if processedCodes(i,4) == 0
                % Blank ACon
                    % Use End: ((End - Latch) / 3) + Latch
                    for j = 1 : 3
                        processedCodes(i,j+2) = (j * (processedCodes(i,5) - processedCodes(i,2)) / 3) + processedCodes(i,2);
                    end
                else
                % Working ACon
                    % Use ACon: ((Acon - (Latch - Start)) / 2) + Latch
                    for j = 1 : 2
                        processedCodes(i,j+2) = (j * (processedCodes(i,4) - (processedCodes(i,2) - processedCodes(i,1))) / 3) + processedCodes(i,2);
                    end
                end
            else
            % Working Buffer
                % Use Buffer: Buffer + Start
                processedCodes(i,3) = processedCodes(i,3) + processedCodes(i,1);
                if processedCodes(i,4) == 0
                % Blank ACon
                    % Use End: ((End - Buffer) / 2) + Start
                    for j = 1 : 2
                        processedCodes(i,j+3) = (j * (processedCodes(i,5) - processedCodes(i,3)) / 2) + processedCodes(i,1);
                    end
                else
                % Working ACon
                    % Use ACon: ACon + Start
                    processedCodes(i,4) = processedCodes(i,4) + processedCodes(i,1);
                end
            end
        end
     end
end       

%% Process

function Process
    % Imported Data
    global diagMode
    %global verbose
    global testID
    global recordList
    global startsEnds
    global numRecords
    global maxLength

    startsEnds = ones(numRecords, 2);
    noiseToleranceEnd = 500;
    noiseToleranceStart = noiseToleranceEnd / 20; % Make sure divisible
    
    for i = 1 : numRecords

        threshhold = max(recordList(i,5:end,1)) / 50;
        flagStart = true;
        flagEnd = true;
        minDuration = 50;
        noiseBufferEnd = noiseToleranceEnd;
        noiseBufferStart = noiseToleranceStart;
        if diagMode == true && i == testID
            % 1 = timeStep;       2 = current;     3 = threshhold; 
            % 4 = flagStart;      5 = bufferStart;
            % 6 = currentStart;   7 = potentialStart;
            % 8 = flagEnd;        9 = bufferEnd;
            %10 = currentEnd;    11 = potentialEnd;
            diag = zeros(maxLength,11);
        end

        for j = 1 : maxLength
            if flagStart == true && recordList(i,j,1) > threshhold  
                noiseBufferStart = noiseBufferStart - 1;
                if noiseBufferStart < 1   
                    if j <= 11 + noiseToleranceStart
                        startsEnds(i,1) = 1;
                    else
                        startsEnds(i,1) = j - (10 + noiseToleranceStart);
                    end
                    flagStart = false;
                end
            else
                noiseBufferStart = noiseToleranceStart;
            end
            if flagStart == false
                minDuration = minDuration - 1;
                if minDuration < 0 && flagEnd == true
                    if j == maxLength
                        startsEnds(i,2) = j - (noiseToleranceEnd - noiseBufferEnd);
                        flagEnd = false;
                        if diagMode == true && i == testID
                            disp("maxLength")
                        end
                    elseif recordList(i,j,1) < threshhold
                        noiseBufferEnd = noiseBufferEnd - 1;
                        if threshhold < max(recordList(i,:,1))
                            threshhold = threshhold * ...
                                (1 + ( 1 * ((noiseToleranceEnd - noiseBufferEnd) / noiseToleranceEnd)));
                        end
                        if noiseBufferEnd < 1
                            if (j - noiseToleranceEnd) + 5 < length(recordList(i,:,1))
                                startsEnds(i,2) = (j - noiseToleranceEnd) + 5;
                                if diagMode == true && i == testID
                                    disp("foundEnd")
                                end
                            else
                                startsEnds(i,2) = length(recordList(i,:,1));
                                if diagMode == true && i == testID
                                    disp("reachedEnd")
                                end
                            end
                            flagEnd = false;
                        end
                    else
                        noiseBufferEnd = noiseToleranceEnd;
                    end
                end
            end
            if diagMode == true && i == testID
                %threshhold
                diag(j,1)  = j;
                diag(j,2)  = recordList(i,j,1);
                diag(j,3)  = threshhold;
                diag(j,4)  = flagStart;
                diag(j,5)  = noiseBufferStart;
                diag(j,6)  = startsEnds(i,1);
                diag(j,7)  = j - (noiseToleranceStart - noiseBufferStart);
                diag(j,8)  = flagEnd;
                diag(j,9)  = noiseBufferEnd;
                diag(j,10) = startsEnds(i,2);
                diag(j,11)  = j - (noiseToleranceEnd - noiseBufferEnd);
            end
        end
        
        % Wipe Excess
        %recordList(i,startsEnds(i,2):end,1) = 0;
        
        if diagMode == true && i == testID
            disp(startsEnds(i,:))
            disp(startsEnds(i,2) - startsEnds(i,1))
        end        
    end

    if diagMode == true
        figure
        record = recordList(testID,startsEnds(testID,1):startsEnds(testID,2),1);
        plot(record)
        threshhold = max(recordList(i,5:end,1)) / 50;
        line([0 length(record)] ,[threshhold threshhold], 'Color', 'g', 'LineWidth', 1);
    end   
end

%% Feature Extraction

function FeatureExtraction
    % Imported Data
    global numRecords
    global maxLength
    global recordList
    global startsEnds
    global calcMode
    global chunks
    global processedCodes
    global codes
    global allFeatures
    global CGEARFeatLabels
    
    allFeatures = cell(3,1);
    %for now stick to 1
    
    for z = 1 : size(allFeatures,1) 
        chunkSize = 3; % Number of Chunks
        calcMode  = z;
        % Mode 1  = Fixed Window
        % Mode 2  = Variable Window using Proportions
        % Mode 3  = Variable Window using Kelman Features
        % Mode 4  = Variable Window using RMS Method
        % Mode 5  = Variable Window using Change Point Method

        switch calcMode
            case 1
                ChunkingFixed(chunkSize);
            case 2
                ChunkingProportions(chunkSize);
            case 3
                ChunkingKelman;
            case 4
                ChunkingRMS;
            case 5
                ChunkingChangePoint;
        end
        
        featureSet = zeros(numRecords,1);
        for k = 1 : 1 % for now stick to 1
        % Current, then Voltage
            for i = 1 : numRecords
                % featureNumber = Tracking Index for Feature
                fN = 1;
                i
                if i == 473
                    disp('ha')
                end
                % Whole-Trace
                record = recordList(i,chunks(i,1):chunks(i,end));
                features = GeneralStats(record);
                featureSet(i,fN:fN+length(features)-1) = features; fN=fN+length(features);
                features = DetailedStats(record);
                featureSet(i,fN:fN+length(features)-1) = features; fN=fN+length(features);
                % Need Contextual Information
                record = recordList(i,:);
                features = KnowledgeFeatures(record,chunks(i,1),chunks(i,end),codes(i,:),startsEnds(i,:));
                featureSet(i,fN:fN+length(features)-1) = features; fN=fN+length(features);
                
                if k == 1
                    % Chunked-Trace Features
                    for j = 1 : size(chunks,2)-1
                        record = recordList(i,chunks(i,j):chunks(i,j + 1));
                        features = GeneralStats(record);
                        featureSet(i,fN:fN+length(features)-1) = features; fN=fN+length(features);
                        features = StdDevFilter(record);
                        featureSet(i,fN:fN+length(features)-1) = features; fN=fN+length(features);
                        features = GradientFilter(record);
                        featureSet(i,fN:fN+length(features)-1) = features; fN=fN+length(features);
                    end

                    % Labelling Features
                    features = codes(i,1:5);
                    featureSet(i,fN:fN+length(features)-1) = features; fN=fN+length(features);
                    %         features = processedCodes(i,1:5);
                    %         allFeatures(i,fN:fN+length(features)-1) = features; fN=fN+length(features);
                end
            end

        % Save Feature Set
        allFeatures{z,k} = featureSet;
        
        end
    end

end

%% Chunking Record into Fixed Windows

function ChunkingFixed(chunkSize)
    % Imported Data
    global numRecords
    global recordList
    global processedCodes
    global startsEnds
    global chunks
    
    kelman = false; % False to use startsEnds
    if ~exist('chunkSize','var')
        chunkSize = 5;
    end
    chunks = zeros(size(recordList,1),chunkSize + 1);
    if kelman == true
        % Picking Max Range Is Too Much Variance
        % 98% of Range in ms
        chunkLength = quantile(processedCodes(:,5) - processedCodes(:,1),0.98) / chunkSize;  
    else
        % Picking Max Range Is Too Much Variance
        % 98% of Range in ms
        chunkLength = quantile(((startsEnds(:,2) * 1000 ./ processedCodes(:,11)) - ...
            (startsEnds(:,1) * 1000 ./ processedCodes(:,11))),0.98) / chunkSize;       
    end
    
    for i = 1 : numRecords
        if kelman == true
            % in Time Steps
            chunks(i,1) = ceil(processedCodes(i,1) * 1000 / processedCodes(i,11)); 
        else
            chunks(i,1) = startsEnds(i,1);
        end
        for j = 1 : chunkSize
            % in Time Steps
            chunks(i,j + 1) = chunks(i,1) + floor((chunkLength * processedCodes(i,11) * j / 1000));
        end
    end
end

%% Chunking Record into Proportioned Windows

function ChunkingProportions(chunkSize)
    % Imported Data
    global numRecords
    global maxLength
    global recordList
    global processedCodes
    global startsEnds
    global chunks
    
    kelman = false; % False to use startsEnds
    if ~exist('chunkSize','var')
        chunkSize = 5;
    end
    chunks = zeros(size(recordList,1),chunkSize + 1);
    
    for i = 1 : numRecords
        if kelman == true
            chunkLength = (processedCodes(i,5) - processedCodes(i,1)) / chunkSize;   
            chunks(i,1) = ceil(processedCodes(i,1) * processedCodes(i,11) / 1000); 
        else
            chunkLength = ((startsEnds(i,2) - startsEnds(i,1)) * 1000 ./ processedCodes(i,11)) / chunkSize;
            chunks(i,1) = startsEnds(i,1);
        end
        for j = 1 : chunkSize
            % in Time Steps
            chunks(i,j + 1) = chunks(i,1) + floor(chunkLength * processedCodes(i,11) * j / 1000);
        end
    end
end

%% Chunking Record using Kelman Features

function ChunkingKelman
    % Imported Data
    global numRecords
    global maxLength
    global recordList
    global processedCodes
    global chunks
    
    incBuffer = true;
    % There are two variants:
    % 1: incBuffer == False:  Start, Latch, ACon, End
    % 2: incBuffer == True:   Start, Latch, Buffer, ACon, End
    
    if incBuffer == false
        chunkSize = 4;
    else
        chunkSize = 3;
    end
    
    chunks = zeros(size(recordList,1),chunkSize + 1);
    
    for i = 1 : numRecords

        if incBuffer == true
            chunks(i,1:5) = floor(processedCodes(i,1:5) * processedCodes(i,11) / 1000);
        else
            chunks(i,1:4) = floor(processedCodes(i,[1 2 4 5]) * processedCodes(i,11) / 1000);
        end
    end    
end

%% Chunking Using RMS [Not Implemented for Now]
%  Not convinced it works across models

function ChunkingRMS
end

%% Chunking Using ChangePoint [Not Implemented for Now]
%  I didn't fully understand implementation

function ChunkingChangePoint
end

%% Calculating General Statistical Spread Features

function features = GeneralStats(record)
    features(1) = sum(record);
    features(2) = max(record);
    features(3) = mean(record);
    features(4) = median(record);
    features(5) = iqr(record);
    features(6) = length(record);
    features(7) = features(2) / features(3);
    if isfinite(features(7)) ~= true
        features(7) = 9999;
    end
    if isnan(features(7)) == true
        features(7) = 9999;
    end
end

%% Calculating Detailed Statistical Spread Features

function features = DetailedStats(record)
    features(1) = min(record);
    features(2) = quantile(record, 0.05);
    features(3) = quantile(record, 0.95);
    features(4) = mode(record);
    features(5) = std(record);   
end

%% Calculating the Knowledge-Based Features (Timing Based Ones Etc.)
function features = KnowledgeFeatures(record,starting,ending,codes,startsEnds)
    
    curFreq = codes(11);
%     kelman = false; % False to use startsEnds
%     if kelman == true
%         curDelay = codes(12);
%         recordEndTime = codes(5);
%     else
%         recordStart = startsEnds(1);
%         recordEndTime = startsEnds(2);
%         curDelay = recordStart * curFreq / 1000;
%     end
    
    record = record(starting:ending);
    if length(record) < 10 | sum(isnan(record)) > 0 | sum(isinf(record)) > 0% Invalid Trace
        features(1:10) = -100;
    else
        % Start and End Times
        features(1) = starting * 1000 / curFreq;
        features(7) = ending * 1000 / curFreq;    

        %%% Need to Save These for Later!!! %%% ALSO GO THROUGH DIFFERENT
        %%% CHANGE POINTSSSSSSSSSS
        for i = 3 : 5
            [ipt, residual] = findchangepts(record,'Statistic','linear','MaxNumChanges',i);
            features(5 + i) = residual;
        end

        length(ipt)
        if length(ipt) < 2 || length(record) < 50 % Something is too wrong to continue
            features(2:6) = -99;
        else
        %%% Need to Save These for Later!!! %%%     

            % Looking for Buffer, First Flip, Then Pick First Peak
            recordFlipped = quantile(record,0.99) - record();
            % Could try with and without smoothing?
            truncStart = ceil(length(recordFlipped)*0.05);
            truncEnd = ceil(length(recordFlipped)*0.95);
            [pks,locs] = findpeaks(smooth(recordFlipped(truncStart:truncEnd)),'NPeaks',2, ...
                'Annotate','extents','MinPeakProminence',0.4, ...
                'MinPeakDistance',50);
            % Buffer
            if isempty(locs) % If Peak Method doesn't work, ChangePoint backup
                [ipt, residual] = findchangepts((recordFlipped(truncStart:truncEnd)), ...
                    'Statistic','linear','MaxNumChanges',4);
                if length(ipt) < 3 % Something is too wrong to continue
                    features(2:6) = -98;
                else
                    buff = truncStart + ipt(2);
                end  
            else
                buff = (truncStart + locs(1));
            end
            if features(2) ~= -98
                features(4) = buff * 1000 / curFreq;

            %     % Looking
            %     figure
            %     findpeaks(smooth(recordFlipped(truncStart:truncEnd)),'NPeaks',2, ...
            %         'Annotate','extents','MinPeakProminence',0.4,'MinPeakDistance',50);  
            %     figure
            %     plot(smooth(record(1:buff)))

                % Looking for Ipk, Widest Peak %
                peakBounds = 0.53;
                peakDists = 56;
                pks2 = [];
                try
                    while isempty(pks2) && peakBounds > 0.01 % Reduce Threshold until a Peak is found
                        peakBounds = peakBounds * 0.75;
                        peakDists = min([peakDists * 0.9, length(record(1:buff))-5]);
                        [pks2, locs2, widths2] = findpeaks(smooth(record(1:buff)), 'NPeaks',2, ...
                            'Annotate','extents','MinPeakProminence',peakBounds, ...
                            'MinPeakDistance',peakDists,'MinPeakWidth',peakDists);
                    end
                    if peakBounds <= 0.01
                        features([2:3 5:6]) = -97;
                    end
                catch % Something is too wrong to continue
                    features([2:3 5:6]) = -97;
                end

                if features(2) ~= -97
                    [M2,I2] = max(widths2);
                    peak = locs2(I2);
                    % Peak
                    features(2) = locs2(I2) * 1000 / curFreq;

                    % Looking for Latch, Flip Then Check Change-Point, Looking for Inflection
                    % Between Peak and Buffer
                    % CHECK WHEN BUMP IS LARGE
                    recordPeakToBufferFlipped = quantile(record(peak: buff),0.99) ...
                        - record(peak:buff);
                    ipt3 = findchangepts(recordPeakToBufferFlipped,'Statistic','linear','MaxNumChanges',3);
                    if length(ipt3) < 2
                        features([2 4:5]) = -96;
                    else
                        % Latch
                        try % Something is too wrong to continue
                            features(3) = (peak + ipt3(2)) * 1000 / curFreq; % Latch
                        catch
                            features([3 5:6]) = -95;
                        end

                        if features(3) ~= -95 
                            % Looking for ACon, Elbow Point and Nearest Change-Point
                            % Distance Perpendicular to Line Between Buffer and End Point
                            % THIS IS SO JANKY, ALTERNATIVE?
                            recordBufferToEnd = record(buff :end);

                            nPoints = length(recordBufferToEnd);
                            allCoord = [1:length(recordBufferToEnd);recordBufferToEnd]';
                            firstPoint = allCoord(1,:);
                            lineVec = allCoord(end,:) - firstPoint;
                            lineVecN = lineVec / sqrt(sum(lineVec.^2));
                            vecFromFirst = bsxfun(@minus, allCoord, firstPoint);
                            scalarProduct = dot(vecFromFirst, repmat(lineVecN,nPoints,1),2);
                            vecFromFirstParallel = scalarProduct * lineVecN;
                            vecToLine = vecFromFirst - vecFromFirstParallel;
                            fun = @(a,b) not(or(a,b)) ;
                            testing99 = bsxfun(fun,  sign(vecToLine(:,1))-1,sign(vecToLine(:,2))-1);
                            distToLine = sqrt(sum(vecToLine.^2,2)) .* testing99;
                            [maxDist,idxOfBestPoint] = max(distToLine);
                            % ACon 
                            features(5) = (buff + allCoord(idxOfBestPoint,1)) * 1000 / curFreq;

                            % Acon, Nearest Change Point

                            [ipt] = findchangepts(record,'Statistic','linear','MaxNumChanges',5);
                            ipt = ipt .* (1000 / curFreq);
                            [k i] = min(abs(ipt - features(5)));
                            features(6) = ipt(i);
                        end
                    end
                end
            end
        end
    end
end

%% Calculating Standard Deviation Filtered Features

function features = StdDevFilter(record)
    record = stdfilt(record);
    features = GeneralStats(record);
    temp = DetailedStats(record);
    features = [features, temp];
end

%% Calculating Gradient Filtered Features

function features = GradientFilter(record)
    record = gradient(record);
    features = GeneralStats(record);
    temp = DetailedStats(record);
    features = [features, temp];
end

%% Feature Selection

function FeatureSelection()
    % Imported Data
    global allFeatures
    global labels
    global labelledSet
    global unlabelledSet
    
    labelledSet = cell(size(allFeatures,1),size(allFeatures,2));
    unlabelledSet = cell(size(allFeatures,1),size(allFeatures,2));
    
    for i = 1 : size(allFeatures,1)
        for j = 1 : size(allFeatures,2)
            temp = allFeatures{i,j};
            unlabelledSet{i,j} = temp(size(labels,1)+1:end,:);
            temp = temp(1:size(labels,1),:);
            temp(:,end+1) = labels(:,j);
            labelledSet{i,j} = temp;
        end
    end
    
end

function ModelTraining
    % Imported Data
    global labelledSet
    global classifiers
    global finalClassifier
    global labels
    
    
    classifiers     = cell(size(labelledSet,1),3,4,size(labels,2));
    finalClassifier = cell(1,4,size(labels,2));

    for m = 1 : size(labels,2)
        trees           = cell(size(labelledSet,1),3,5,2);
        svms            = cell(size(labelledSet,1),4,5,2);
        ensembles       = cell(size(labelledSet,1),3,5,2);
        
        for i = 1 : size(labelledSet,1)
            % Convert to Table
            inputTable = array2table(labelledSet{i,m});
            predictors = inputTable(:, 1 : (end-size(labels,2)));
            response = inputTable(:, (end-size(labels,2))+m);
            isCategoricalPredictor = false(width(predictors),1);

            %%% Trees

            for j = 1 : 3
                for k = 1 : 5
                    switch j
                        case 1 % Fine Tree Classifier
                            trainedClassifier = ClassifierTreeFine(predictors, response);
                        case 2 % Medium Tree Classifier
                            trainedClassifier = ClassifierTreeMedium(predictors, response);
                        case 3 % Coarse Tree Classifier
                            trainedClassifier = ClassifierTreeCoarse(predictors, response);
                    end
                    trees = AsignResults(trees, i, j, k, trainedClassifier);
                end
            end
            classifiers = SelectBest(classifiers,trees,i,1,m);
            %Wipe trees (Memory?)
            if i == size(labelledSet,1)
                trees = {};
            end

            %%% SVMs
            for j = 1 : 4
                for k = 1 : 5
                    switch j
                        case 1 % Linear SVM
                            trainedClassifier = ClassifierSVMLinear(predictors, response);
                        case 2 % Cubic SVM
                            trainedClassifier = ClassifierSVMCubic(predictors, response);
                        case 3 % Fine SVM
                            trainedClassifier = ClassifierSVMFine(predictors, response);
                        case 4 % Coarse SVM
                            trainedClassifier = ClassifierSVMCoarse(predictors, response);
                    end
                    svms = AsignResults(svms, i, j, k, trainedClassifier); 
                end
            end
            classifiers = SelectBest(classifiers,svms,i,2,m);
            %Wipe SVMs (Memory?)
            if i == size(labelledSet,1)
                svms = {};
            end

            %%% Ensembles
            for j = 1 : 3
                for k = 1 : 5
                    switch j
                        case 1 % Ensemble Boosted
                            trainedClassifier = ClassifierEnsembleBoosted(predictors, response);
                        case 2 % Ensemble Bagged
                            trainedClassifier = ClassifierEnsembleBagged(predictors, response);
                        case 3 % Ensemble RUS Boosted
                            trainedClassifier = ClassifierEnsembleRUSBoosted(predictors, response);
                    end
                    ensembles = AsignResults(ensembles, i, j, k, trainedClassifier);
                end
            end
            classifiers = SelectBest(classifiers,ensembles,i,3,m);
            %Wipe Ensembles (Memory?)
            if i == size(labelledSet,1)
                ensembles = {};
            end
        end

        combinedResults = horzcat(classifiers{:,:,4,m});
        for i = 1 : size(combinedResults,2)/2
            if combinedResults(1,(2*(i-1))+1) + combinedResults(1,(2*(i-1))+2) == 1
            % 0-1 Scale Already, just use first value for 0-1 Scale
                temp(:,i) = combinedResults(:,(2*(i-1))+1);
            else
                if combinedResults(1,(2*(i-1))+1) + combinedResults(1,(2*(i-1))+2) == 0
                % -1-1 Scale, (Add 1, Divide by 2) for 0-1 Scale
                    temp(:,i) = (combinedResults(:,(2*(i-1))+1) + 1) / 2;
                else
                % Different Scale, (#1 / (#1 + #2)) for 0-1 Scale
                    temp(:,i) = combinedResults(:,(2*(i-1))+1) ./ (combinedResults(:,(2*(i-1))+1) + combinedResults(:,(2*(i-1))+2));
                end
            end
        end
        combinedResults = array2table(temp);
        predictors = combinedResults;
        %combinedResults(:,end) = table2array(response);
        trainedClassifier = ClassifierSubspaceDiscriminant(predictors, response);
        partitionedModel = crossval(trainedClassifier.classifier, 'KFold', 5);
        [validationPredictions, validationScores] = kfoldPredict(partitionedModel);
        validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');
        finalClassifier{1,1,m} = trainedClassifier;
        finalClassifier{1,2,m} = validationAccuracy;
        finalClassifier{1,3,m} = validationPredictions;
        finalClassifier{1,4,m} = validationScores;
    end
end

function classifierGroup = AsignResults(classifierGroup, i, j, k, trainedClassifier)
    partitionedModel = crossval(trainedClassifier.classifier, 'KFold', 5);
    [validationPredictions, validationScores] = kfoldPredict(partitionedModel);
    validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');
    classifierGroup{i,j,k,1} = trainedClassifier;
    classifierGroup{i,j,k,2} = validationAccuracy;
    classifierGroup{i,j,k,3} = validationPredictions;
    classifierGroup{i,j,k,4} = validationScores;
end

function classifiers = SelectBest(classifiers, classifierGroup, i, j,m)
        % Compare Means of Class and Pick Max Then Pick Max Value of Class
        [maxMean, idMean] = max(mean(transpose(squeeze(cell2mat(classifierGroup(i,:,:,2))))));
        [maxVal, idVal] = max(squeeze(cell2mat(classifierGroup(i,idMean,:,2))));
        minVal = min(squeeze(cell2mat(classifierGroup(i,idMean,:,2))));
        classifiers{i,j,1,m} = classifierGroup{i,idMean,idVal,1};
        classifiers{i,j,2,m} = [minVal, maxMean, maxVal];
        classifiers{i,j,3,m} = squeeze(cell2mat(classifierGroup(i,idMean,idVal,3)));
        classifiers{i,j,4,m} = squeeze(cell2mat(classifierGroup(i,idMean,idVal,4)));
end

function trainedClassifier = ClassifierTreeFine(predictors, response)
    predictorNames = predictors.Properties.VariableNames;
    isCategoricalPredictor = false(width(predictors),1);
    
    classifier = fitctree(...
    predictors, ...
    response, ...
    'SplitCriterion', 'gdi', ...
    'MaxNumSplits', 100, ...
    'Surrogate', 'off', ...
    'ClassNames', [0; 1]);
    predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
    PredictFcn = @(x) predict(classifier, x);
    trainedClassifier.predictFcn = @(x) PredictFcn(predictorExtractionFcn(x));
    trainedClassifier.classifier = classifier;
end

function trainedClassifier = ClassifierTreeMedium(predictors, response)
    predictorNames = predictors.Properties.VariableNames;
    isCategoricalPredictor = false(width(predictors),1);
    classifier = fitctree(...
        predictors, ...
        response, ...
        'SplitCriterion', 'gdi', ...
        'MaxNumSplits', 20, ...
        'Surrogate', 'off', ...
        'ClassNames', [0; 1]);
    predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
    PredictFcn = @(x) predict(classifier, x);
    trainedClassifier.predictFcn = @(x) PredictFcn(predictorExtractionFcn(x));
    trainedClassifier.classifier = classifier;
end

function trainedClassifier = ClassifierTreeCoarse(predictors, response)
predictorNames = predictors.Properties.VariableNames;    
isCategoricalPredictor = false(width(predictors),1);
    classifier = fitctree(...
        predictors, ...
        response, ...
        'SplitCriterion', 'gdi', ...
        'MaxNumSplits', 20, ...
        'Surrogate', 'off', ...
        'ClassNames', [0; 1]);
    predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
    PredictFcn = @(x) predict(classifier, x);
    trainedClassifier.predictFcn = @(x) PredictFcn(predictorExtractionFcn(x));
    trainedClassifier.classifier = classifier;
end

function trainedClassifier = ClassifierSVMLinear(predictors, response)
predictorNames = predictors.Properties.VariableNames;    
isCategoricalPredictor = false(width(predictors),1);
    classifier = fitcsvm(...
        predictors, ...
        response, ...
        'KernelFunction', 'linear', ...
        'PolynomialOrder', [], ...
        'KernelScale', 'auto', ...
        'BoxConstraint', 1, ...
        'Standardize', true, ...
        'ClassNames', [0; 1]);
    predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
    PredictFcn = @(x) predict(classifier, x);
    trainedClassifier.predictFcn = @(x) PredictFcn(predictorExtractionFcn(x));
    trainedClassifier.classifier = classifier;
end

function trainedClassifier = ClassifierSVMCubic(predictors, response)
predictorNames = predictors.Properties.VariableNames;    
isCategoricalPredictor = false(width(predictors),1);
    classifier = fitcsvm(...
        predictors, ...
        response, ...
        'KernelFunction', 'polynomial', ...
        'PolynomialOrder', 3, ...
        'KernelScale', 'auto', ...
        'BoxConstraint', 1, ...
        'Standardize', true, ...
        'ClassNames', [0; 1]);
    predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
    PredictFcn = @(x) predict(classifier, x);
    trainedClassifier.predictFcn = @(x) PredictFcn(predictorExtractionFcn(x));
    trainedClassifier.classifier = classifier;
end

function trainedClassifier = ClassifierSVMFine(predictors, response)
predictorNames = predictors.Properties.VariableNames;    
isCategoricalPredictor = false(width(predictors),1);
    classifier = fitcsvm(...
        predictors, ...
        response, ...
        'KernelFunction', 'gaussian', ...
        'PolynomialOrder', [], ...
        'KernelScale', 1.6, ...
        'BoxConstraint', 1, ...
        'Standardize', true, ...
        'ClassNames', [0; 1]);
    predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
    PredictFcn = @(x) predict(classifier, x);
    trainedClassifier.predictFcn = @(x) PredictFcn(predictorExtractionFcn(x));
    trainedClassifier.classifier = classifier;
end

function trainedClassifier = ClassifierSVMCoarse(predictors, response)
predictorNames = predictors.Properties.VariableNames;    
isCategoricalPredictor = false(width(predictors),1);
    classifier = fitcsvm(...
        predictors, ...
        response, ...
        'KernelFunction', 'gaussian', ...
        'PolynomialOrder', [], ...
        'KernelScale', 25, ...
        'BoxConstraint', 1, ...
        'Standardize', true, ...
        'ClassNames', [0; 1]);
    predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
    PredictFcn = @(x) predict(classifier, x);
    trainedClassifier.predictFcn = @(x) PredictFcn(predictorExtractionFcn(x));
    trainedClassifier.classifier = classifier;
end

function trainedClassifier = ClassifierEnsembleBoosted(predictors, response)
predictorNames = predictors.Properties.VariableNames;    
isCategoricalPredictor = false(width(predictors),1);
    template = templateTree(...
        'MaxNumSplits', 20);
    classifier = fitcensemble(...
        predictors, ...
        response, ...
        'Method', 'AdaBoostM1', ...
        'NumLearningCycles', 30, ...
        'Learners', template, ...
        'LearnRate', 0.1, ...
        'ClassNames', [0; 1]);
    predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
    PredictFcn = @(x) predict(classifier, x);
    trainedClassifier.predictFcn = @(x) PredictFcn(predictorExtractionFcn(x));
    trainedClassifier.classifier = classifier;
end

function trainedClassifier = ClassifierEnsembleBagged(predictors, response)
predictorNames = predictors.Properties.VariableNames;    
isCategoricalPredictor = false(width(predictors),1);
    template = templateTree(...
        'MaxNumSplits', 99);
    classifier = fitcensemble(...
        predictors, ...
        response, ...
        'Method', 'Bag', ...
        'NumLearningCycles', 30, ...
        'Learners', template, ...
        'ClassNames', [0; 1]);
    predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
    PredictFcn = @(x) predict(classifier, x);
    trainedClassifier.predictFcn = @(x) PredictFcn(predictorExtractionFcn(x));
    trainedClassifier.classifier = classifier;
end

function trainedClassifier = ClassifierEnsembleRUSBoosted(predictors, response)
    predictorNames = predictors.Properties.VariableNames;    
    isCategoricalPredictor = false(width(predictors),1);
    template = templateTree(...
        'MaxNumSplits', 20);
    classifier = fitcensemble(...
        predictors, ...
        response, ...
        'Method', 'RUSBoost', ...
        'NumLearningCycles', 30, ...
        'Learners', template, ...
        'LearnRate', 0.1, ...
        'ClassNames', [0; 1]);
    predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
    PredictFcn = @(x) predict(classifier, x);
    trainedClassifier.predictFcn = @(x) PredictFcn(predictorExtractionFcn(x));
    trainedClassifier.classifier = classifier;
end

function trainedClassifier = ClassifierSubspaceDiscriminant(predictors, response)
    predictorNames = predictors.Properties.VariableNames;    
    isCategoricalPredictor = false(width(predictors),1);
    subspaceDimension = max(1, min(4, width(predictors) - 1));
    classifier = fitcensemble(...
        predictors, ...
        response, ...
        'Method', 'Subspace', ...
        'NumLearningCycles', 30, ...
        'Learners', 'discriminant', ...
        'NPredToSample', subspaceDimension, ...
        'ClassNames', [0; 1]);
    predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
    PredictFcn = @(x) predict(classifier, x);
    trainedClassifier.predictFcn = @(x) PredictFcn(predictorExtractionFcn(x));
    trainedClassifier.classifier = classifier;
end

function ModelApplying
    % Imported Data
    global classifiers
    global finalClassifier
    global unlabelledSet
    global testResults
    global finalResults
    global labels
    
    finalResults = []; % Clear previous stuff
    
    testResults = {}; %% SET UP
    for m = 1 : size(labels,2)
        for i = 1 : size(unlabelledSet,1)
            % Convert to Table
            inputTable = array2table(unlabelledSet{i});
            predictors = inputTable;
            isCategoricalPredictor = false(width(predictors),1);
            
            for j = 1 : 3
                % Do I need to run each multiple times??
                trainedClassifier = classifiers{i,j,1,m};
                partitionedModel = crossval(trainedClassifier.classifier, 'KFold', 5);
                [validationPredictions, validationScores] = kfoldPredict(partitionedModel);
                validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');
                [predictions scores] = trainedClassifier.predictFcn(table2array(predictors));
                testResults{i,j,1,m} = predictions;
                testResults{i,j,2,m} = scores;
            end
        end 
        
        combinedResults = horzcat(testResults{:,:,2,m});
        for i = 1 : size(combinedResults,2)/2
            if combinedResults(1,(2*(i-1))+1) + combinedResults(1,(2*(i-1))+2) == 1
                % 0-1 Scale Already, just use first value for 0-1 Scale
                temp(:,i) = combinedResults(:,(2*(i-1))+1);
            else
                if combinedResults(1,(2*(i-1))+1) + combinedResults(1,(2*(i-1))+2) == 0
                    % -1-1 Scale, (Add 1, Divide by 2) for 0-1 Scale
                    temp(:,i) = (combinedResults(:,(2*(i-1))+1) + 1) / 2;
                else
                    % Different Scale, (#1 / (#1 + #2)) for 0-1 Scale
                    temp(:,i) = combinedResults(:,(2*(i-1))+1) ./ (combinedResults(:,(2*(i-1))+1) + combinedResults(:,(2*(i-1))+2));
                end
            end
        end
        % Previous code chunk like this uses cells?
        combinedResults = array2table(temp);
        predictors = combinedResults;
        trainedClassifier = finalClassifier{1,1,m};
        [predictions, scores] = trainedClassifier.predictFcn(table2array(predictors));
        finalResults(:,:,m) = [predictions, abs((scores(:,1)-0.5) .* 2)];
    end
end

function ViewResults
%%% Creates a window for each of the time-based features.
%   Wordy so following abbreviations used:
%   N  = Population (Unlabelled Set), n  = Sample (Labelled Set)
%   g  = Good in Sample,              G  = Good in Population
%   g' = Predicted Good in Sample     G' = Predicted Good in Population
%   b  = Bad in Sample,               B  = Bad in Population
%   b' = Predicted Bad in Sample      B' = Predicted Bad in Population
%   Note: G and B are unknowns.
%   
%   Each window will have 5 figures (colours are red-green-red)
%   Figure 1 (All   Row 1): g' vs g       = How good is filter
%   Figure 2 (All   Row 2): n  vs N       = How good is sample
%   Figure 3 (Third Row 3): G' vs g vs g' = Another view
%   Figure 4 (Third Row 3): B' vs b vs b' = Another view
%   Figure 5 (Third Row 3): STATS: NEED TO WORK ON THIS!!!!
%                           precision, recall, similarity of dists, etc...

%%% To Do:
%           Get Stats Section Filled Out
%           Think how to make it work with multiple labels, e.g. Voltage

    % Imported Data
    global finalResults
    global finalClassifier
    global codes
    global labels
    global thresholds
    
    thresholds = cell(4,4);
    
    labelledResults = finalClassifier{3};
    
    for i = 1 : 4
        switch i    % It's hard to fit an overall title in...
            case 1  % ...so using window name instead
                message = 'Latch Times - Start Times';
                subtractFlag = 0; % Whether to subtract from prev. feat.
            case 2
                message = 'Buffer Times - Latch Times';
                subtractFlag = 1;
            case 3
                message = 'Auxiliary Contact Times - Buffer Times';
                subtractFlag = 1;
            case 4
                message = 'End Times - Auxiliary Contact Times';
                subtractFlag = 1;
        end

        figure('name',message)
        labelLen = length(labelledResults);

        % codes           =  N  
        % finalResults    = (N-n) + G' + B' (w/o g' + b')
        % labelledResults =  n    + g' + b'
        % labels          =  n    + g  + b
        
        
        % Figure 1 (All Row 1): g' vs g = How good is filter
        subplot(3,1,1)
        % Select g'
        selectedData = codes(1:labelLen,i+1) - ...
            codes(1:labelLen,i) * subtractFlag;
        selectedData = selectedData(labelledResults > 0); % CHANGE FOR MULTI-LABEL
        % Select g
        selectedData2 = codes(1:labelLen,i+1) - ...
        codes(1:labelLen,i) * subtractFlag;
        selectedData2 = selectedData2(labels(:) > 0); % CHANGE FOR MULTI-LABEL
        % Contrast g' vs g
        selectedData = {selectedData, selectedData2}; 
        PlotResults(selectedData)
        title(strcat({'Comparing Distribution of Labelled Known ('}, ...
            {num2str(length(selectedData{2}))}, ...
            {' / '}, {num2str(labelLen)}, ...
            {') and Labelled Predicted ('}, ...
            {num2str(length(selectedData{1}))}, ...
            {' / '}, {num2str(labelLen)}, ...
            {') Valid Results'}),'FontSize',9)
        
        % Save Distribution of g' for later use
        [f,xi] = ksdensity(selectedData2);
        thresholds{i,1} = f;
        thresholds{i,2} = xi;
        thresholds{i,3} = ksdensity(selectedData2,'Function','cdf');
        where = find(thresholds{i,3} < 0.05);
        temp  = find(thresholds{i,3} > 0.95);
        thresholds{i,4} = [xi(where(end)) xi(temp(1))];
        

        % Figure 2 (All Row 2): N-n vs N vs n = How good is sample
        subplot(3,1,2)
        % Select N - n
        selectedData = codes(labelLen + 1 : end,i+1) - ...
            codes(labelLen + 1 : end,i) * subtractFlag;
        % Select n
        selectedData2 = codes(1 : labelLen,i+1) - ...
            codes(1 : labelLen,i) * subtractFlag;
        % Contrast (N-n) vs n+(N-n) vs n
        selectedData = {selectedData, [selectedData2;selectedData], selectedData2};
        PlotResults(selectedData)
        title(strcat({'Comparing Distribution of Population ('}, ...
            {num2str(length(selectedData{2}))}, ...
            {') and Labelled Sample ('}, ...
            {num2str(length(selectedData{3}))}, ...
            {') and Population - Sample ('}, ...
            {num2str(length(selectedData{1}))}, ...
            {') Valid Results '}),'FontSize',9)
        
        % Figure 3 (Third Row 3): G' vs g vs g' = Another View
        subplot(3,3,7)
        % Select G' (w/o g')
        selectedData = codes(labelLen + 1 : end,i+1) - ...
            codes(labelLen + 1 : end,i) * subtractFlag;
        selectedData = selectedData(finalResults(:,1) > 0);
        % Select g'
        selectedData2 = codes(1 : labelLen,i+1) - ...
            codes(1 : labelLen,i) * subtractFlag;
        selectedData2 = selectedData2(labelledResults(:) > 0); % CHANGE FOR MULTI-LABELS
        % Contrast (G') vs TEMP vs g'
        selectedData = {[selectedData; selectedData2], selectedData, selectedData2};
        % Select g
        selectedData2 = codes(1 : labelLen,i+1) - ...
            codes(1 : labelLen,i) * subtractFlag;
        selectedData2 = selectedData2(labels(:) > 0); % CHANGE FOR MULTI-LABELS        
        % Contrast (G') vs g vs g'
        selectedData{2} = selectedData2;
        PlotResults(selectedData)
        title(strcat({'Comparing Distribution of Labelled Known ('}, ...
            {num2str(length(selectedData{2}))}, ...
            {' / '}, {num2str(labelLen)}, ...    
            {') and Labelled Predicted ('}, ...
            {num2str(length(selectedData{3}))}, ...
            {' / '}, {num2str(labelLen)}, ...    
            {') and Population Predicted ('}, ...
            {num2str(length(selectedData{1}))}, ...
            {' / '}, {num2str(size(codes,1))}, ...
            {') Valid Results '}),'FontSize',9)
              
        % Figure 4 (Third Row 3): B' vs b vs b' = Another View
        subplot(3,3,8)
        % Select B' (w/o b')
        selectedData = codes(labelLen + 1 : end,i+1) - ...
            codes(labelLen + 1 : end,i) * subtractFlag;
        selectedData = selectedData(finalResults(:,1) == 0);
        % Select b'
        selectedData2 = codes(1 : labelLen,i+1) - ...
            codes(1 : labelLen,i) * subtractFlag;
        selectedData2 = selectedData2(labelledResults(:) == 0); % CHANGE FOR MULTI-LABELS
        % Contrast (B') vs TEMP vs b'
        selectedData = {[selectedData; selectedData2], selectedData, selectedData2};
        % Select b
        selectedData2 = codes(1 : labelLen,i+1) - ...
            codes(1 : labelLen,i) * subtractFlag;
        selectedData2 = selectedData2(labels(:) == 0); % CHANGE FOR MULTI-LABELS        
        % Contrast (B') vs b vs b'
        selectedData{2} = selectedData2;
        PlotResults(selectedData)
        title(strcat({'Comparing Distribution of Labelled Known ('}, ...
            {num2str(length(selectedData{2}))}, ...
            {' / '}, {num2str(labelLen)}, ...    
            {') and Labelled Predicted ('}, ...
            {num2str(length(selectedData{3}))}, ...
            {' / '}, {num2str(labelLen)}, ...    
            {') and Population Predicted ('}, ...
            {num2str(length(selectedData{1}))}, ...
            {' / '}, {num2str(size(codes,1))}, ...
            {') Invalid Results '}),'FontSize',9)
        
        % Figure 5 (Third Row 3): Stats = NEED TO DECIDE!
        
        
%         subplot(3,1,1)
%         % All Unlabelled Valid
%         where = finalResults(:,1) > 0;
%         selectedData = codes(length(labelledResults) + 1 : end,i+1);
%         if subtractFlag == 1
%             selectedData = selectedData - codes(length(labelledResults) + 1 : end,i);
%         end
%         selectedData = selectedData(where);
%         % All Labelled Valid
%         where = labelledResults(:) > 0;
%         temp  = codes(where,i+1);
%         if subtractFlag == 1
%             temp = temp - codes(where,i);
%         end
%         % Combined
%         selectedData = [selectedData; temp];
%         PlotResults(selectedData)
%         title('Distribution of All Predicted Valid Results','FontSize',9) 
%         title(strcat({'Distribution of All ('}, {num2str(length(selectedData))}, ...
%             {' / '}, {num2str(size(codes,1))},{') Predicted Valid Results'}),'FontSize',9)
%         % Save Distribution
%         [f,xi] = ksdensity(selectedData);
%         thresholds{i,1} = f;
%         thresholds{i,2} = xi;
%         thresholds{i,3} = ksdensity(selectedData,'Function','cdf');
%         where = find(thresholds{i,3} < 0.05);
%         temp  = find(thresholds{i,3} > 0.95);
%         thresholds{i,4} = [xi(where(end)) xi(temp(1))];
%         
%         % Predicted Good of Unlabelled vs Predicted Labelled Good of Labelled
%         subplot(3,2,3)
%         % All Unlabelled Valid
%         where = finalResults(:,1) > 0;
%         selectedData = codes(length(labelledResults) + 1: end,i+1);
%         if subtractFlag == 1
%             selectedData = selectedData - codes(length(labelledResults) + 1: end,i);
%         end
%         selectedData = selectedData(where);
%         % All Labelled Valid
%         where = labelledResults(:) > 0;
%         temp  = codes(where,i+1);
%         if subtractFlag == 1
%             temp = temp - codes(where,i);
%         end
%         % Contrasted
%         selectedData = {selectedData, temp};
%         PlotResults(selectedData)
%         title(strcat({'Comparing Distribution of Labelled ('}, {num2str(length(selectedData{2}))}, ...
%             {') and Unlabelled ('}, {num2str(length(selectedData{1}))},{') Predicted Valid Results'}),'FontSize',9)
% 
% 
%         % Predicted Bad of Unlabelled vs Predicted Bad of Labelled
%         subplot(3,2,4)
%         % All Unlabelled Invalid
%         where = finalResults(:,1) == 0;
%         selectedData = codes(length(labelledResults) + 1: end,i+1);
%         if subtractFlag == 1
%             selectedData = selectedData - codes(length(labelledResults) + 1: end,i);
%         end
%         selectedData = selectedData(where);
%         % All Labelled Invalid
%         where = labelledResults(:) == 0;
%         temp = codes(where,i+1);
%         if subtractFlag == 1
%             temp = temp - codes(where,i);
%         end
%         % Contrasted
%         selectedData = {selectedData, temp};
%         PlotResults(selectedData)
%         title(strcat({'Comparing Distribution of Labelled ('}, {num2str(length(selectedData{2}))}, ...
%             {') and Unlabelled ('}, {num2str(length(selectedData{1}))},{') Predicted Invalid Results'}),'FontSize',9)
% 
%         % Known Combined of Unlabelled vs Predicted Combined of Labelled
%         subplot(3,2,5)
%         % All Unlabelled
%         selectedData = codes(length(labelledResults) + 1: end,i+1);
%         if subtractFlag == 1
%             selectedData = selectedData - codes(length(labelledResults) + 1: end,i);
%         end
% 
%         % All Labelled
%         temp = codes(1 : length(labelledResults),i+1);
%         if subtractFlag == 1
%             temp = temp - codes(1 : length(labelledResults),i);
%         end
%         % Contrasted
%         selectedData = {selectedData, temp};
%         PlotResults(selectedData)
%         title(strcat({'Comparing Distribution of Labelled ('}, {num2str(length(selectedData{2}))}, ...
%             {') and Unlabelled ('}, {num2str(length(selectedData{1}))},{') Results'}),'FontSize',9)
% 
% 
%         % Predicted Good of Labelled vs Known Good of Labelled
%         subplot(3,2,6)
%         % All Labelled Valid Predicted
%         where = labels(:) > 0;
%         selectedData = codes(1 : length(labelledResults),i+1);
%         if subtractFlag == 1
%             selectedData = selectedData - codes(1 : length(labelledResults),i);
%         end
%         selectedData = selectedData(where);
%         % All Labelled Valid Known
%         where = labelledResults(:) > 0;
%         temp = codes(1 : length(labelledResults),i+1);
%         if subtractFlag == 1
%             temp = temp - codes(1 : length(labelledResults),i);
%         end
%         % Contrasted
%         selectedData = {selectedData, temp(where)};
%         PlotResults(selectedData)
%         title(strcat({'Comparing Distribution of Labelled ('}, {num2str(length(selectedData{1}))}, ...
%             {') and Predicted ('}, {num2str(length(selectedData{2}))},{') Valid Results'}),'FontSize',9)

    %     % Construct a Legend with the data from the sub-plots
    %     hL = legend([1,2,3,4],{'Data Axes 1','Data Axes 2','Data Axes 3','Data Axes 4'});
    %     % Programatically move the Legend
    %     newPosition = [0.4 0.4 0.2 0.2];
    %     newUnits = 'normalized';
    %     set(hL,'Position', newPosition,'Units', newUnits);
    end
end

function PlotResults(selectedData)
    
    if iscell(selectedData) == true
    % Multiple Stacks for Same Plot
        for i = 1 : size(selectedData,2)
            if size(selectedData{i}) == 1
                temp = selectedData{i};
                temp(2) = temp(1) + 100;
                temp(3) = temp(1) - 100;
                selectedData{i} = transpose(temp);
            end
        end
        if size(selectedData,2) == 2
            if isempty(selectedData{1}) == 1
                temp = selectedData;
                selectedData{1} = temp{3};
                selectedData{3} = [];
                selectedData = selectedData(~cellfun('isempty',selectedData));
            end

            if isempty(selectedData{2}) == 1
                disp('oh no...');
                selectedData{2} = [999,999,999];
            end
        end
        if size(selectedData,2) == 3
            if isempty(selectedData{3}) == 0
                selectedData = selectedData(~cellfun('isempty',selectedData));
            end
        end
        

        for i = 1 : size(selectedData,2)
            comparisons(i) = length(selectedData{i});
        end
        numIntervals = ceil(1 + 3.3*log(max(comparisons)));
        numIntervals = max([numIntervals 20*size(selectedData,2)]);
        for i = 1 : size(selectedData,2)
            comparisons(i)  = min(selectedData{i});
            comparisons2(i) = quantile(selectedData{i},0.95)*1.1;
        end
        intervalWidth = (max(comparisons2) ...
            - min(comparisons)) / numIntervals;
        intervalWidth = min([intervalWidth 3]);
        intervalWidth = max([intervalWidth 1]);
        x = min(comparisons) : intervalWidth : ...
            max(comparisons2);
        
        nCount = histc(selectedData{1},x);
        for i = 2 : size(selectedData,2)
            nCount  = [nCount, histc(selectedData{i},x)];
        end
        temp = selectedData;
        temp2 = nCount;
        for i = 1 : size(selectedData,2)
            relativeFreq(:,i) = nCount(:,i) / length(selectedData{i});
        end

    else
        numIntervals = 1 + 3.3*log(length(selectedData));
        numIntervals = max([numIntervals 40]);
        intervalWidth = (quantile(selectedData,0.95)*1.1 - min(selectedData)) / numIntervals;
        intervalWidth = min([intervalWidth 3]);
        intervalWidth = max([intervalWidth 1]);
        x = min(selectedData) : intervalWidth : quantile(selectedData,0.95)*1.1;
        nCount = histc(selectedData,x);
        relativeFreq = nCount / length(selectedData);
    end
    
    

    
    graph = bar(x , relativeFreq, 0.8,'FaceColor',[0.2 .4 .4], ...
        'FaceAlpha', 0.8, 'EdgeColor',[0.2 .4 .4]);
    ylim([0 1]);
    ylim manual
    xlim manual
    
    
    if iscell(selectedData) == true
        for i = 1 : size(selectedData,2)
            comparison = quantile(selectedData{i},0.95)*1.1;
        end        
    else
        comparison = quantile(selectedData,0.95)*1.1;
    end
    xlim([0 max(comparison)]);
    if iscell(selectedData) == true
        colour = [[0.55,0.3,0.55];[0.3 .55 .55]; [0.55,0.3,0.55]];
        set(graph, 'BarWidth', 1)
        for i = 1 : size(selectedData,2)
            set(graph(i),'FaceColor',colour(i,:))
        end
    end
    hold on;
    yyaxis right
    ylabel('Relative Frequency')
    ax = gca;
    ax.YColor = [0.2 .4 .4];
    
    
    yyaxis left
    ylabel('Probability Density Estimate')
    ax = gca;
    ax.YColor = [0.3 .65 .65];
    ylim([0 1]);
    
    if iscell(selectedData) == true
    % Multiple Stacks for Same Plot
        flag = size(selectedData,2) - 1;
    else
        flag = 0;
    end
    
    ylim([0 0.5]);
    for i = 1 + (flag ~= 0)  : 1 + (flag ~= 0) + flag
        if i ~= 1
            selectedData = temp{i - 1};
            nCount = temp2(:,i - 1);
            colour = [[0.55,0.3,0.55];[0.3 .55 .55]; [0.55,0.3,0.55]];
            alpha  = 0.6;
        else
            colour = [0.4 .75 .75];
            alpha  = 0.6;
        end
        [f,xi] = ksdensity(selectedData);
        plot(xi,f, 'Color', colour(i - (flag ~= 0),:),'Marker', 'none',  'LineWidth', 2, 'LineStyle', '-');
        ar = area(xi,f, 'FaceColor', colour(i - (flag ~= 0),:), 'FaceAlpha', alpha);
        tally = ksdensity(selectedData,'Function','cdf');
        plot(xi,tally, 'Color', colour(i - (flag ~= 0),:),'Marker', 'none',  'LineWidth', 1, 'LineStyle', '-')
        
        where = tally < (0.05*tally(end));
        where = find(where);
        plot([xi(where(end)) xi(where(end))], [0 max(nCount)], 'Color', colour(i - (flag ~= 0),:),'Marker', 'none',  'LineWidth', 2, 'LineStyle', '- -')
        plot(xlim, [0.95 0.95], 'Color', colour(i - (flag ~= 0),:),'Marker', 'none',  'LineWidth', 0.5, 'LineStyle', '- -')
        message = strcat({'5th Percentile'},{'\n'},{num2str(xi(where(end)))});
        message = sprintf(message{1});
        %text(xi(where(end)),max(nCount)/2,message, 'BackgroundColor', 'w', 'HorizontalAlignment', 'right');

        where = tally > (0.95*tally(end));
        where = find(where);
        plot([xi(where(1)) xi(where(1))], [0 max(nCount)], 'Color', colour(i - (flag ~= 0),:),'Marker', 'none',  'LineWidth', 2, 'LineStyle', '- -')
        plot(xlim, [0.05 0.05], 'Color', colour(i - (flag ~= 0),:),'Marker', 'none',  'LineWidth', 0.5, 'LineStyle', '- -')
        message = strcat({'95th Percentile'},{'\n'},{num2str(xi(where(1)))});
        message = sprintf(message{1});
        %text(xi(where(1)),max(nCount)/2,message, 'BackgroundColor', 'w', 'HorizontalAlignment', 'left');
        if max(f) > 0.1
            ylim([0 1]);
        end
    end
    %legend();
    hold off;
end

function Diagnose
    % Imported Data
    global thresholds
    global codes
    global finalClassifier
    global finalResults
    global diagnoses
    global recordList
    
    labelledResults = finalClassifier{3};
    diagnoses = zeros(size(codes,1),4,2);
    
    for i = 1 : 4
        if i == 1
            flag  = 0;
        else
            flag  = 1;
            
        end
        % Diagnose All
        selectedData = codes(:,i+1);
        if flag == 1
            selectedData = selectedData - codes(:,i);
        end
        % Record Percentile for All
        for j = 1 : length(selectedData)
            where = find(selectedData(j) > thresholds{i,2});
            if isempty(where) == true
                where = 100;
            end
            where = where(end);
            temp  = thresholds{i,3};
            temp  = temp(where);
            diagnoses(j,i,2) = temp;
        end
        
        boundary = thresholds{i,4};
        boundary = boundary(1);
        where = selectedData(:) <= boundary;
        diagnoses(where,i,1) = -1;

        boundary = thresholds{i,4};
        boundary = boundary(2);
        where = selectedData(:) >= boundary;
        diagnoses(where,i,1) = 1;
        
        % Wipe Predicted Bad of All
        % All Unlabelled Invalid
        where = find(finalResults(:,1) == 0);
        where = where + 100;
        where = [find(labelledResults(:,1) == 0) ; where];
        diagnoses(where,:,1) = -9;
        
        
        figure
        subplot(3,1,1)
        % Diagnose All
        
        selectedData = codes(:,i+1);
        if flag == 1
            selectedData = selectedData - codes(:,i);
        end
        % Record Percentile for All
        where = diagnoses(:,i,1) == -1;
        selectedData = selectedData(where);
        
        temp = selectedData;
        selectedData = codes(:,i+1);
        if flag == 1
            selectedData = selectedData - codes(:,i);
        end
        % Record Percentile for All
        where = diagnoses(:,i,1) == 0;
        selectedData = selectedData(where);
        selectedData = {temp, selectedData};
       
        temp = selectedData;
        selectedData = codes(:,i+1);
        if flag == 1
            selectedData = selectedData - codes(:,i);
        end
        % Record Percentile for All
        where = diagnoses(:,i,1) == 1;
        selectedData = selectedData(where);
        
        temp2 = selectedData;
        selectedData = temp;
        selectedData{1,3} = temp2;
        
        PlotResults(selectedData)
        
        subplot(3,1,2)
        % Diagnose All
        selectedData = codes(:,i+1);
        if flag == 1
            selectedData = selectedData - codes(:,i);
        end
        % Record Percentile for All
        where = diagnoses(:,i,1) == 0;
        selectedData = selectedData(where);
        PlotResults(selectedData)
        
        subplot(3,1,3)
        % Diagnose All
        selectedData = codes(:,i+1);
        if flag == 1
            selectedData = selectedData - codes(:,i);
        end
        % Record Percentile for All
        where = diagnoses(:,i,1) > 0;
        selectedData = selectedData(where);
        PlotResults(selectedData)
        
        
        
        
        
        figure;
        temp = diagnoses(:,i,1);
        temp2 = find(temp == -1);
        subplot(4,1,1)
        for j = 1 : min([50 length(temp2)])
            hold on
            plot(0 : (1000/codes(temp2(j),10)) : ((1000/codes(temp2(j),10))*(3000 - codes(temp2(j),12))),(recordList(temp2(j),codes(temp2(j),12):3000)));
            hold off
        end
        
        temp2 = find(temp == 0);
        subplot(4,1,2)
        for j = 1 : min([50 length(temp2)])
            hold on
            plot(0 : (1000/codes(temp2(j),10)) : ((1000/codes(temp2(j),10))*(3000 - codes(temp2(j),12))),(recordList(temp2(j),codes(temp2(j),12):3000)));
        end
        
        temp2 = find(temp == 1);
        subplot(4,1,3)
        for j = 1 : min([50 length(temp2)])
            hold on
            plot(0 : (1000/codes(temp2(j),10)) : ((1000/codes(temp2(j),10))*(3000 - codes(temp2(j),12))),(recordList(temp2(j),codes(temp2(j),12):3000)));
        end
        
        temp2 = find(temp == -9);
        subplot(4,1,4)
        for j = 1 : min([50 length(temp2)])
            hold on
            plot(0 : (1000/codes(temp2(j),10)) : ((1000/codes(temp2(j),10))*(3000 - codes(temp2(j),12))),(recordList(temp2(j),codes(temp2(j),12):3000)));
        end
    end
end

function DiagnoseTesting
     % Imported Data
    global finalResults
    global finalClassifier
    global codes
    global labels
    global allFeatures
    global recordList
%     global thresholds
%     
%     thresholds = cell(4,4);
    
    labelledResults = finalClassifier{3};
    relevantFeatures = allFeatures{2,1};
    relevantFeatures = relevantFeatures(:,18) ./ relevantFeatures(:,30);
    relevantFeatures(isnan(relevantFeatures)) = 0;
    relevantFeatures(~isfinite(relevantFeatures)) = 100;
    relevantFeatures = relevantFeatures * 10;
    
    for i = 1 : 1
        switch i
            case 1
                message = 'Latch Times - Start Times';
                flag  = 0;
            case 2
                message = 'Buffer Times - Latch Times';
                flag  = 1;
            case 3
                message = 'Auxiliary Contact Times - Buffer Times';
                flag  = 1;
            case 4
                message = 'End Times - Auxiliary Contact Times';
                flag  = 1;
        end

        figure('name',message)

        % Predicted Good of All
        subplot(3,1,1)
        % All Unlabelled Valid
        where = finalResults(:,1) > 0;
        selectedData = relevantFeatures(length(labelledResults) + 1 : end);
        if flag == 1
            selectedData = selectedData - codes(length(labelledResults) + 1 : end,i);
        end
        selectedData = selectedData(where);
        % All Labelled Valid
        where = labelledResults(:) > 0;
        temp = relevantFeatures(where);
        if flag == 1
            temp = temp - codes(where,i);
        end
        % Combined
        selectedData = [selectedData; temp];
        PlotResults(selectedData)
        title('Distribution of All Predicted Valid Results','FontSize',9) 
        title(strcat({'Distribution of All ('}, {num2str(length(selectedData))}, ...
            {' / '}, {num2str(size(codes,1))},{') Predicted Valid Results'}),'FontSize',9)
        % Save Distribution
        [f,xi] = ksdensity(selectedData);
        thresholds{i,1} = f;
        thresholds{i,2} = xi;
        thresholds{i,3} = ksdensity(selectedData,'Function','cdf');
        where = find(thresholds{i,3} < 0.05);
        temp  = find(thresholds{i,3} > 0.95);
        thresholds{i,4} = [xi(where(end)) xi(temp(1))];
        
        % Predicted Good of Unlabelled vs Predicted Labelled Good of Labelled
        subplot(3,2,3)
        % All Unlabelled Valid
        where = finalResults(:,1) > 0;
        selectedData = relevantFeatures(length(labelledResults) + 1 : end);
        if flag == 1
            selectedData = selectedData - codes(length(labelledResults) + 1: end,i);
        end
        selectedData = selectedData(where);
        % All Labelled Valid
        where = labelledResults(:) > 0;
        temp = relevantFeatures(where);
        if flag == 1
            temp = temp - codes(where,i);
        end
        % Contrasted
        selectedData = {selectedData, temp};
        PlotResults(selectedData)
        title(strcat({'Comparing Distribution of Labelled ('}, {num2str(length(selectedData{2}))}, ...
            {') and Unlabelled ('}, {num2str(length(selectedData{1}))},{') Predicted Valid Results'}),'FontSize',9)


        % Predicted Bad of Unlabelled vs Predicted Bad of Labelled
        subplot(3,2,4)
        % All Unlabelled Invalid
        where = finalResults(:,1) == 0;
        selectedData = relevantFeatures(length(labelledResults) + 1 : end);
        if flag == 1
            selectedData = selectedData - codes(length(labelledResults) + 1: end,i);
        end
        selectedData = selectedData(where);
        % All Labelled Invalid
        where = labelledResults(:) == 0;
        temp = relevantFeatures(where)
        if flag == 1
            temp = temp - codes(where,i);
        end
        % Contrasted
        selectedData = {selectedData, temp};
        PlotResults(selectedData)
        title(strcat({'Comparing Distribution of Labelled ('}, {num2str(length(selectedData{2}))}, ...
            {') and Unlabelled ('}, {num2str(length(selectedData{1}))},{') Predicted Invalid Results'}),'FontSize',9)

        % Known Combined of Unlabelled vs Predicted Combined of Labelled
        subplot(3,2,5)
        % All Unlabelled
        selectedData = relevantFeatures(length(labelledResults) + 1 : end);
        if flag == 1
            selectedData = selectedData - codes(length(labelledResults) + 1: end,i);
        end
        selectedData = selectedData(where);
        % All Labelled
        where = labelledResults(:) == 0;
        temp = relevantFeatures(1 : length(labelledResults));
        if flag == 1
            temp = temp - codes(1 : length(labelledResults),i);
        end
        % Contrasted
        selectedData = {selectedData, temp};
        PlotResults(selectedData)
        title(strcat({'Comparing Distribution of Labelled ('}, {num2str(length(selectedData{2}))}, ...
            {') and Unlabelled ('}, {num2str(length(selectedData{1}))},{') Results'}),'FontSize',9)


        % Predicted Good of Labelled vs Known Good of Labelled
        subplot(3,2,6)
        % All Labelled Valid Predicted
        where = labels(:) > 0;
        selectedData = relevantFeatures(1 : length(labelledResults));
        if flag == 1
            selectedData = selectedData - codes(1 : length(labelledResults),i);
        end
        selectedData = selectedData(where);
        % All Labelled Valid Known
        where = labelledResults(:) > 0;
        temp = relevantFeatures(1 : length(labelledResults));
        if flag == 1
            temp = temp - codes(1 : length(labelledResults),i);
        end
        % Contrasted
        selectedData = {selectedData, temp(where)};
        PlotResults(selectedData)
        title(strcat({'Comparing Distribution of Labelled ('}, {num2str(length(selectedData{1}))}, ...
            {') and Predicted ('}, {num2str(length(selectedData{2}))},{') Valid Results'}),'FontSize',9)
    end
    
    
    labelledResults = finalClassifier{3};
    diagnoses = zeros(size(codes,1),4,2);
    for i = 1 : 1
        if i == 1
            flag  = 0;
        else
            flag  = 1;
            
        end
        % Diagnose All
        selectedData = relevantFeatures(:);
        if flag == 1
            selectedData = selectedData - codes(:,i);
        end
        % Record Percentile for All
        for j = 1 : length(selectedData)
            where = find(selectedData(j) > thresholds{i,2});
            if isempty(where) == true
                where = 100;
            end
            where = where(end);
            temp  = thresholds{i,3};
            temp  = temp(where);
            diagnoses(j,i,2) = temp;
        end
        
        boundary = thresholds{i,4};
        boundary = boundary(1);
        where = selectedData(:) < boundary;
        diagnoses(where,i,1) = -1;

        boundary = thresholds{i,4};
        boundary = boundary(2);
        where = selectedData(:) > boundary;
        diagnoses(where,i,1) = 1;
        
        % Wipe Predicted Bad of All
        % All Unlabelled Invalid
        where = finalResults(:,1) < 0;
        selectedData = relevantFeatures(length(labelledResults) + 1 : end);
        if flag == 1
            selectedData = selectedData - codes(length(labelledResults) + 1 : end,i);
        end
        selectedData = selectedData(where);
        % All Labelled Invalid
        where = labelledResults(:) < 0;
        temp  = relevantFeatures(where);
        if flag == 1
            temp = temp - codes(where,i);
        end
        % Combined
        selectedData = [selectedData; temp];
        where = selectedData < 99999;
        diagnoses(where,i,1) = 0;
        
        
        figure
        subplot(3,1,1)
        % Diagnose All
        
        % Diagnose All

        
        selectedData = relevantFeatures(:);
        if flag == 1
            selectedData = selectedData - codes(:,i);
        end
        % Record Percentile for All
        where = diagnoses(:,i,1) == -1;
        selectedData = selectedData(where);
        
        temp = selectedData;
        selectedData = relevantFeatures(:);
        if flag == 1
            selectedData = selectedData - codes(:,i);
        end
        % Record Percentile for All
        where = diagnoses(:,i,1) == 0;
        selectedData = selectedData(where);
        selectedData = {temp, selectedData};
       
        temp = selectedData;
        selectedData = relevantFeatures(:);
        if flag == 1
            selectedData = selectedData - codes(:,i);
        end
        % Record Percentile for All
        where = diagnoses(:,i,1) == 1;
        selectedData = selectedData(where);
        
        temp2 = selectedData;
        selectedData = temp;
        selectedData{1,3} = temp2;
        
        PlotResults(selectedData)
        
        subplot(3,1,2)
        % Diagnose All
        selectedData = relevantFeatures(:);
        if flag == 1
            selectedData = selectedData - codes(:,i);
        end
        % Record Percentile for All
        where = diagnoses(:,i,1) == 0;
        selectedData = selectedData(where);
        PlotResults(selectedData)
        
        subplot(3,1,3)
        % Diagnose All
        selectedData = relevantFeatures(:);
        if flag == 1
            selectedData = selectedData - codes(:,i);
        end
        % Record Percentile for All
        where = diagnoses(:,i,1) == 1;
        selectedData = selectedData(where);
        PlotResults(selectedData)
    end
    
    
temp = diagnoses(:,1,1);    
temp2 = find(temp == 0);
figure;
title('')
for j = 1 : min([50 length(temp2)])
    hold on
            plot(0 : (1000/codes(temp2(j),10)) : ((1000/codes(temp2(j),10))*(3000 - codes(temp2(j),12))),(recordList(temp2(j),codes(temp2(j),12):3000)));
    hold off
end

temp = diagnoses(:,1,1);    
temp2 = find(temp == 1);
figure;
for j = 1 : min([50 length(temp2)])
    hold on
            plot(0 : (1000/codes(temp2(j),10)) : ((1000/codes(temp2(j),10))*(3000 - codes(temp2(j),12))),(recordList(temp2(j),codes(temp2(j),12):3000)));
    hold off
end

temp = diagnoses(:,1,1);    
temp2 = find(temp == -1);
figure;
for j = 1 : min([50 length(temp2)])
    hold on
            plot(0 : (1000/codes(temp2(j),10)) : ((1000/codes(temp2(j),10))*(3000 - codes(temp2(j),12))),(recordList(temp2(j),codes(temp2(j),12):3000)));
    hold off
end


end

function PlotTrace(i,newFig)
        % Imported Data
    global model
    global source
    global allData
    global startsEnds
	kelmanFeatures = evalin('base','kelman_features');
    
    % New: 2 = New each time, 1 = New and add to same, 0 = add to same
    if newFig == 1
        figure
        hold on
    else
        hold on
    end

    for j = 1 : length(i)
        if newFig == 2
            figure
        end
        % Reads and Asigns Values
        subplot(4,1,1)
        ID = i(j);
        recordID = model(ID);
        record = source{1,recordID};
        record(end + 1 : 13000) = 0;
        record(isnan(record)) = 0;

        curMax = record(1);
        curScale = record(2);
        curFreq = record(3);
        curDelay = record(4);

        recordStart = (curDelay + 4)*1000/curFreq;
        recordLatch = kelmanFeatures(recordID,1);
        recordBuffer = kelmanFeatures(recordID,2);
        recordMCon = kelmanFeatures(recordID,3);
        recordACon = kelmanFeatures(recordID,4);
        recordEndTime = kelmanFeatures(recordID,5);

        recordpk1 = (kelmanFeatures(recordID,6));
        recordplt = (kelmanFeatures(recordID,7));

        record = (record/curScale)*curMax;
        
        
        recordV = allData{2,recordID};
        
        volMax = recordV(1);
        volScale = recordV(2);
        volFreq = recordV(3);
        volDelay = recordV(4);
        
        recordV = (recordV/volScale)*volMax;

        hold on;
        plot(0 : (1000/curFreq) : ((1000/curFreq)*(length(record)))-(1000/curFreq), record);
%         plot(0 : (1000/curFreq) : ((1000/curFreq)*(length(record)))-(1000/curFreq), stdfilt(record));
%         plot(0 : (1000/curFreq) : ((1000/curFreq)*(length(record)))-(1000/curFreq), abs(gradient(record)));
        ylim([0 10])        
        hold on;
        
        plot([recordLatch + recordStart recordLatch + recordStart], [0 10])
        plot([recordBuffer + recordStart recordBuffer + recordStart], [0 10])
        plot([recordACon + recordStart recordACon + recordStart], [0 10])
        plot([recordEndTime + recordStart recordEndTime + recordStart], [0 10])
        
        yyaxis right
        
        %plot(( (1000/volFreq)) : (1000/volFreq) : ((1000/volFreq)*(length(recordV)))-(5 * 1000/volFreq), recordV(6:end));
        plot((5 * (1000/volFreq)) : (1000/volFreq) : ((1000/volFreq)*(length(recordV)))-(1000/volFreq), movmean(recordV(6:end),5));
        plot((5 * (1000/volFreq)) : (1000/volFreq) : ((1000/volFreq)*(length(recordV)))-(1000/volFreq), recordV(6:end));
        %         plot([recordStart recordStart], [0 10])

        plot([recordMCon + recordStart recordMCon + recordStart], [0 50])
        ylim([30 40])
        hold off;
        
        xlim([(recordStart-5)*0.9 (recordEndTime+recordStart)*1.1])
        
        title(i)
        global setA7T1
        global setA7T2
        global setA7T3
        
        
        subplot(4,1,2)
        
        plot(0 : (1000/curFreq) : ((1000/curFreq)*(length(record)))-(1000/curFreq), record);
        xlim([recordStart recordEndTime+recordStart])
        hold on;
        plot([setA7T1(i,1) setA7T1(i,1)], [0 10]) % Start
        plot([setA7T1(i,2)+setA7T1(i,1) setA7T1(i,2)+setA7T1(i,1)], [0 10]) % Peak
        plot([setA7T1(i,3)+setA7T1(i,1) setA7T1(i,3)+setA7T1(i,1)], [0 10]) % Latch
        plot([setA7T1(i,4)+setA7T1(i,1) setA7T1(i,4)+setA7T1(i,1)], [0 10]) % Buffer
        plot([setA7T1(i,5)+setA7T1(i,1) setA7T1(i,5)+setA7T1(i,1)], [0 10]) % ACon
        plot([setA7T1(i,6)+setA7T1(i,1) setA7T1(i,6)+setA7T1(i,1)], [0 10]) % ACon
        plot([setA7T1(i,7) setA7T1(i,7)], [0 10]) % No Offset for End Time (?)
        hold off;

        subplot(4,1,3)
        
        plot(0 : (1000/curFreq) : ((1000/curFreq)*(length(record)))-(1000/curFreq), record);
        xlim([recordStart recordEndTime+recordStart])
        hold on;
        plot([setA7T2(i,1) setA7T2(i,1)], [0 10]) % Start
        plot([setA7T2(i,2)+setA7T2(i,1) setA7T2(i,2)+setA7T2(i,1)], [0 10]) % Peak
        plot([setA7T2(i,3)+setA7T2(i,1) setA7T2(i,3)+setA7T2(i,1)], [0 10]) % Latch
        plot([setA7T2(i,4)+setA7T2(i,1) setA7T2(i,4)+setA7T2(i,1)], [0 10]) % Buffer
        plot([setA7T2(i,5)+setA7T2(i,1) setA7T2(i,5)+setA7T2(i,1)], [0 10]) % ACon
        plot([setA7T2(i,6)+setA7T2(i,1) setA7T2(i,6)+setA7T2(i,1)], [0 10]) % ACon
        plot([setA7T2(i,7) setA7T2(i,7)], [0 10]) % No Offset for End Time (?)
        hold off;        
        
        subplot(4,1,4)
        
        plot(0 : (1000/curFreq) : ((1000/curFreq)*(length(record)))-(1000/curFreq), record);
        xlim([recordStart recordEndTime+recordStart])
        hold on;
        plot([setA7T3(i,1) setA7T3(i,1)], [0 10]) % Start
        plot([setA7T3(i,2)+setA7T3(i,1) setA7T3(i,2)+setA7T3(i,1)], [0 10]) % Peak
        plot([setA7T3(i,3)+setA7T3(i,1) setA7T3(i,3)+setA7T3(i,1)], [0 10]) % Latch
        plot([setA7T3(i,4)+setA7T3(i,1) setA7T3(i,4)+setA7T3(i,1)], [0 10]) % Buffer
        plot([setA7T3(i,5)+setA7T3(i,1) setA7T3(i,5)+setA7T3(i,1)], [0 10]) % ACon
        plot([setA7T3(i,6)+setA7T3(i,1) setA7T3(i,6)+setA7T3(i,1)], [0 10]) % ACon
        plot([setA7T3(i,7) setA7T3(i,7)], [0 10]) % No Offset for End Time (?)
        hold off;
         
%          try
%             testing = (quantile(record((curDelay+4+1):(curDelay+4 + (recordEndTime * curFreq / 1000))),0.99) -record((curDelay+4+1):(curDelay+4 + (recordEndTime * curFreq / 1000))));
%          catch
%             testing = (quantile(record((curDelay+4+1):end),0.99) -record((curDelay+4):end)); 
%          end
%          figure
%          findchangepts(testing,'Statistic','linear','MaxNumChanges',5)
%          [ipt] = findchangepts(testing,'Statistic','linear','MaxNumChanges',5);
%          hold on;
% %          findpeaks(record(5:700),'MinPeakProminence',1,'Annotate','extents')
%          findpeaks(smooth(testing),'NPeaks',2,'Annotate','extents','MinPeakProminence',0.4,'MinPeakDistance',50)
%          [pks,locs] = findpeaks(smooth(testing),'NPeaks',2,'Annotate','extents','MinPeakProminence',0.4,'MinPeakDistance',50)
%          figure
%          if isempty(locs)
%              Buff = ipt(2)+(curDelay+4);  
%          else
%              Buff = locs(1)+(curDelay+4);
%          end
%          Buff * 1000 / curFreq % Buffer
%          plot(smooth(record((curDelay+4): Buff)))
%          hold  on;
%          Bounds = 0.533;
%          pks2 = [];
%          while isempty(pks2) % Loosen Rules
%              Bounds = Bounds * 0.75;
%              [pks2, locs2, widths2] = findpeaks(smooth(record((curDelay+4): Buff)),'NPeaks',2,'Annotate','extents','MinPeakProminence',Bounds,'MinPeakDistance',50);         
%          end
%          findpeaks(smooth(record((curDelay+4): Buff)),'NPeaks',2,'Annotate','extents','MinPeakProminence',Bounds,'MinPeakDistance',50)
%          [M2,I2] = max(widths2);
%          ((locs2(I2)+(curDelay+4)) * 1000 / curFreq) % Peak
%          figure
%          
%          testing2 = (quantile((record((curDelay+4) + locs2(I2): Buff)),0.99) -(record((curDelay+4) + locs2(I2): Buff)));
%          
%          %findpeaks(smooth(testing2),'Annotate','extents')
%          %hold on
%          %plot(gradient(gradient(smooth(testing2))))
%          %findpeaks(gradient(gradient(smooth(testing2))),'Annotate','extents')
%          findchangepts(testing2,'Statistic','linear','MaxNumChanges',3)
%          ipt3 = findchangepts(testing2,'Statistic','linear','MaxNumChanges',3);
%          %[pks3, locs3, widths3] = findpeaks(gradient(gradient(smooth(testing2))),'Annotate','extents')
%          %[M3,I3] = max(widths3)
%          %(((curDelay+4) + locs2(I2) + locs3(I3)) * 1000 / curFreq) % Latch
%          (((curDelay+4) + locs2(I2) + ipt3(2)) * 1000 / curFreq) % Latch
%          %figure
%          %plot(smooth(record((curDelay+4) + locs(1): (curDelay+4 + (recordEndTime * curFreq / 1000)))))
%          
%          testing3 = record(((Buff+1) :(curDelay+4 + (recordEndTime * curFreq / 1000))));
%          
%          
%          nPoints = length(testing3);
%          allCoord = [1:length(testing3);testing3]';
%          firstPoint = allCoord(1,:);
%          lineVec = allCoord(end,:) - firstPoint;
%          lineVecN = lineVec / sqrt(sum(lineVec.^2));
%          vecFromFirst = bsxfun(@minus, allCoord, firstPoint);
%          
%          scalarProduct = dot(vecFromFirst, repmat(lineVecN,nPoints,1),2);
%          vecFromFirstParallel = scalarProduct * lineVecN;
%          vecToLine = vecFromFirst - vecFromFirstParallel;
%          fun = @(a,b) not(or(a,b)) ;
%          testing99 = bsxfun(fun,  sign(vecToLine(:,1))-1,sign(vecToLine(:,2))-1);
%          distToLine = sqrt(sum(vecToLine.^2,2)) .* testing99;
%          figure('Name','distance from curve to line'), plot(distToLine)
%         [maxDist,idxOfBestPoint] = max(distToLine);
%          figure, plot(testing3)
%         hold on
%         plot(allCoord(idxOfBestPoint,1), allCoord(idxOfBestPoint,2), 'or')
%          %refLine = linspace(testing3(1),testing3(end),length(testing3));
%         allCoord(idxOfBestPoint,1) 
%         acon1 = ((Buff) + allCoord(idxOfBestPoint,1)) * 1000 / curFreq % Peak
% 
%         [ipt] = findchangepts(testing,'Statistic','linear','MaxNumChanges',5);
%         ipt = ipt .* (1000 / curFreq);
%         [k i] = min(abs(ipt-acon1))
%         ipt(i)
        
%          testing4 = record(((locs(1)+(curDelay+4+1))+idxOfBestPoint :(curDelay+4 + (recordEndTime * curFreq / 1000))));
%          
%          
%          nPoints = length(testing4);
%          allCoord2 = [1:length(testing4);testing4]';
%          firstPoint = allCoord2(1,:);
%          lineVec = allCoord2(end,:) - firstPoint;
%          lineVecN = lineVec / sqrt(sum(lineVec.^2));
%          vecFromFirst = bsxfun(@minus, allCoord2, firstPoint);
%          
%          scalarProduct = dot(vecFromFirst, repmat(lineVecN,nPoints,1),2);
%          vecFromFirstParallel = scalarProduct * lineVecN;
%          vecToLine = vecFromFirst - vecFromFirstParallel;
%          fun = @(a,b) not(or(a,b)) ;
%          testing99 = bsxfun(fun,  sign(vecToLine(:,1))-1,sign(vecToLine(:,2))-1);
%          distToLine = sqrt(sum(vecToLine.^2,2)) .* testing99;
%          figure('Name','distance from curve to line'), plot(distToLine)
%         [maxDist2,idxOfBestPoint2] = max(distToLine);
%          figure, plot(testing4)
%         hold on
%         plot(allCoord2(idxOfBestPoint2,1), allCoord2(idxOfBestPoint2,2), 'or')
%          %refLine = linspace(testing3(1),testing3(end),length(testing3));
%         allCoord2(idxOfBestPoint2,1) 
%         ((locs(1)+(curDelay+4)) + allCoord(idxOfBestPoint,1) + allCoord2(idxOfBestPoint2,1)) * 1000 / curFreq % Peak
    end
    

end
