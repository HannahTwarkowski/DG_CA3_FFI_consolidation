%% Extract final spatial and temporal traces
% This script processes the entire pipeline of post-CNMF-E processing. It
% creates a copy of the CNMF-E file. Neurons that did not qualify to be
% kept, are being replaced by 0s to maintain the original index. They should
% be removed for further analysis (e.g. calcium event rate, correlation).

% It is used in "A dentate  gyrus-CA3  inhibitory  circuit  promotes  evolution 
% of  hippocampal-cortical ensembles during memory consolidation" 
% https://doi.org/10.1101/2021.05.21.445117
% script by Victor Steininger & edits by Hannah Twarkowski



%Input
%1) exel file names of files to load & safe
% column A: name for new file postpipeline [mouse_session_postp] e.g.:(ms3_16A_postp)
  % this file will mimic CNMF-E's structure (neuron.A, neuron.C,
  % neuronC_raw)
% column B: name of CNMF-E file to be used e.g. 18-Apr_20_35_23;
% column C: name of new file with information regardign calcium traces e.g.
  % spike timing (e.g. ms3_16A_postp_peakProps)
  
%2 training data set (TsSet.mat)
%3 polygeom.m matlab function from H.J. Sommer;
% https://www.mathworks.com/matlabcentral/fileexchange/319-polygeom-m;

%output: 
% 1)'postp'-file: filtered file of CNMF-E 'postp', mimics CNMF-E's structure (neuron.A, neuron.C,
  % neuronC_raw)
% 2) 'postp_peakProps' file: includes peak properties such as sparse spike-timing.
% Used in future analysis (e.g. calcium event rate, correlation);
% 

% Read the names for new file postpipeline the file to extract
[~, txt] = xlsread('CFC_BL_Pipeline_Names_Edited_example','A:B');
oname = txt(:,1); % file names given at the end by the pipeline
iname = txt(:,2); % CNMF-E file names
nfiles = size(oname,1);

for n = 1:nfiles

    %% Part 1 - Data extraction
    % CNMF-E data extraction
    %insert pathway for CNMF-E files under fname
    fname = "C:\Users\HTwar\Desktop\CA1_adult_onFinalData\Github\" + iname{n} + ".mat"; %put name of the CNMF-E data file to process
    CNMFE_Data = load(fname);
    neuron = CNMFE_Data.neuron;
    Coor = [];
    neuron.Coor = [];
    Coor = neuron.get_contours;
    neuron.Coor = Coor;
    [numb, siz] = size(neuron.C);

    %% Part 2 - Shape-based neuron cleaning
    % In this part CNMF-E neurons are sorted based on their contours. If their
    % shape is too elongated, they are more likely to be a dendrite or an
    % artifact. The threshold of inertia axes ratio was determined from the
    % dataset of an ACC video and can be adjusted.
    % N.B. you need the file 'polygeom' to run the code. The package can be
    % also found on internet.

    InertRatioThresh = 5; %Set the selection threshold for non-elongated neurons
    InertRatio = zeros(1,numb);
    for i = 1:numb
        cont = Coor(i);
        contD = cont{1}; %extract neuron contour
        xValMax = max(contD(1,:));
        xValMin = min(contD(1,:));
        yValMax = max(contD(2,:));
        yValMin = min(contD(2,:));
        if (xValMax-xValMin <= 0)||(yValMax-yValMin <=0) % if the contour is flat, polygeom is crashing
            InertRatio(i) = InertRatioThresh + 999; % so the neuron is eliminated before hand by setting a value above threshold
        else
            [geom, iner, cpmo] = polygeom(contD(1,:), contD(2,:));
            InertRatio(i) = max(cpmo(1), cpmo(3))/min(cpmo(1), cpmo(3)); %calculates the ratio of long and short axis of the neuron's shape
        end
    end

    % Index neurons with an elongated shape (wrt the threshold)
    IndexShape = [];
    for j = 1:numb
        if InertRatio(j) > InertRatioThresh
            IndexShape = [IndexShape j];
        end
    end


    %% Part 3 - Peaks detection and cleaning
    clear lesLocs lesLocsC lesPeaks lesPeaksC lesWHR lesWHR_C lesProm lesPromC lesPNR lesPNRC;

    baselineSTD = zeros(1,numb);
    thresholdFluo = 3; %this low PNR threshold captures true calcium events and noise

    for j = 1:numb
        TempNeuronRaw = neuron.C_raw(j,:);
        TempNeuronC = neuron.C(j,:);
        %Calculation of the noise variance
        baseline = TempNeuronRaw - neuron.C(j,:); %C subtracted to C_raw gives an estimate of the baseline
        STDTemp = std(baseline);
        baselineSTD(j) = STDTemp;
        threshCompa = STDTemp*thresholdFluo; %The threshold is expressed as a multiple of the noise std
        %Broad extraction of peaks from C and C_raw traces.
        [pks, locs, wid, prom] = findpeaks(TempNeuronRaw,'MinPeakHeight',threshCompa, 'MinPeakDistance',15,'MinPeakWidth',3);
        [pksC, locsC, widC, promC] = findpeaks(TempNeuronC, 'MinPeakHeight',threshCompa, 'MinPeakDistance',15,'MinPeakWidth',3);
        
        lesLocs(j,1:size(locs,2)) = locs; %Temporal location of the peaks
        lesLocsC(j,1:size(locsC,2)) = locsC;
        lesPeaks(j,1:size(pks,2)) = pks;%Height of the peaks (fluorescence amplitude)
        lesPeaksC(j,1:size(pksC,2)) = pksC;
        lesWHR(j,1:size(pks,2)) = pks./wid;%Width-to-height ratio of every peak
        lesWHR_C(j,1:size(pksC,2)) = pksC./widC;
        lesProm(j,1:size(prom,2)) = prom; %Prominence of the peaks
        lesPromC(j,1:size(promC,2)) = promC;
        lesPNR(j,1:size(pks,2)) = pks./STDTemp; %Peak-to-noise ratio of the peaks
        lesPNRC(j,1:size(pksC,2)) = pksC./STDTemp;
    end

    % Activate the following line if you want to store all parameters in one
    % variable
    PeakProp = struct('locs',lesLocs,'locsC',lesLocsC,'peaks',lesPeaks,'peaksC',lesPeaksC,'WHR',lesWHR,'WHRC',lesWHR_C,'prom',lesProm,'promC',lesPromC,'PNR',lesPNR,'PNRC',lesPNRC);

    % The following section will pair detected peaks between the C and C_raw
    % traces. For each peak detected in the C trace, it calculates the temporal 
    % distance to every peak from the raw trace (within a neuron). It finds
    % the peak in the raw trace the closest to C peak. If this distance is still
    % too large, the peaks are not considered perks and discarded. Otherwise, the
    % peaks are paired and share the same 'index number'. The indices of paired
    % peaks are stored in indexRaw and indexC variables.
    indexRaw = zeros(size(lesLocs));
    indexC = zeros(size(lesLocsC));
    peakDistThresh = 12; %maximum acceptable distance between 2 peaks in C_raw and C to be classified as a pair
    indTmp = 0;
    dsTmp = 999;
    for i = 1:numb
        indTmp = -2;
        lenC = size(lesLocsC,2)-size(find(lesLocsC(i,:)==0),2);
        lenRaw = size(lesLocs,2)-size(find(lesLocs(i,:)==0),2);
        if lenRaw >= lenC
            len = size(lesLocsC,2)-size(find(lesLocsC(i,:)==0),2); %number of peaks in the C trace of the ith neuron
            counter = 1;
            for j = 1:len
                dist = abs(lesLocs(i,:) - lesLocsC(i,j)); %calculates the distance between one peak in C and all the raw peaks
                [mini, ind] = min(dist); %find the raw peak closest to the C peak
                if(mini > peakDistThresh)
                    indexC(i,j) = -1; % -1 is used as a code for 'unpaired'
                else
                    if(ind == indTmp) % the new peak in C is associated to a previously paired peak of C_raw
                        if(mini < dstTmp) % the new peak in C is closer to the C_raw than the former partner
                            jjj = find(indexC(i,:)==indexRaw(i,indTmp)); % finding the location of the previous paired peak in C
                            indexC(i,j) = indexRaw(i,indTmp); % the index number is given to the new peak
                            indexC(i,jjj) = -1; % the former partner is discarded
                            dstTmp = dist; % the distance of the pair is updated
                        else
                            indexC(i,j) = -1; % if the new peak is further away it is assigned as unpaired
                        end
                    else % the new peak in C is not yet associated to a peak from C_raw
                        indexC(i,j) = counter;
                        indexRaw(i,ind) = counter;
                        counter = counter + 1;
                        indTmp = ind;
                        dstTmp = dist;
                    end
                end
            end
            len2 = size(lesLocs,2)-size(find(lesLocs(i,:)==0),2);%Iteration in the C_raw trace. Peaks that have not been paired in the previous loop are marked as unpaired
            for k = 1:len2
                if(indexRaw(i,k) == 0)
                    indexRaw(i,k) = -1;
                end
            end
        else
            len = size(lesLocs,2)-size(find(lesLocs(i,:)==0),2); %number of peaks in the C trace of the ith neuron
            counter = 1;
            for j = 1:len
                dist = abs(lesLocsC(i,:) - lesLocs(i,j)); %calculates the distance between one peak in C and all the raw peaks
                [mini, ind] = min(dist); %find the raw peak closest to the C peak
                if(mini > peakDistThresh)
                    indexRaw(i,j) = -1; % -1 is used as a code for 'unpaired'
                else
                    if(ind == indTmp) % the new peak in C is associated to a previously paired peak of C_raw
                        if(mini < dstTmp) % the new peak in C is closer to the C_raw than the former partner
                            jjj = find(indexRaw(i,:)==indexC(i,indTmp)); % finding the location of the previous paired peak in C
                            indexRaw(i,j) = indexC(i,indTmp); % the index number is given to the new peak
                            indexRaw(i,jjj) = -1; % the former partner is discarded
                            dstTmp = dist; % the distance of the pair is updated
                        else
                            indexRaw(i,j) = -1; % if the new peak is further away it is assigned as unpaired
                        end
                    else % the new peak in C is not yet associated to a peak from C_raw
                        indexRaw(i,j) = counter;
                        indexC(i,ind) = counter;
                        counter = counter + 1;
                        indTmp = ind;
                        dstTmp = dist;
                    end
                end
            end
            len2 = size(lesLocsC,2)-size(find(lesLocsC(i,:)==0),2);%Iteration in the C_raw trace. Peaks that have not been paired in the previous loop are marked as unpaired
            for k = 1:len2
                if(indexC(i,k) == 0)
                    indexC(i,k) = -1;
                end
            end
        end
    end

    C = neuron.C;
    C_raw = neuron.C_raw;
    %Calculating correlation between paired peaks in C and C_raw and the local
    %PNR of every peak. These are two other parameters used for the
    %peak classification.
    PairCorr = zeros(size(indexC));
    LocalPNR = zeros(size(indexC));
    for i = 1:numb
        len = size(lesLocsC,2)-size(find(lesLocsC(i,:)==0),2);
        for j = 1:len
            if (indexC(i,j) > 0)
                %Definition of the trace in the local window.
                %These multiple 'if' checks the location of the peaks because the 
                %size of the correlation window has to be adapted if the peak it 
                %too close to the beginning or the end.
                if (lesLocsC(i,j) <= 5) 
                    ind1 = 1;
                    ind2 = lesLocsC(i,j) + 19;
                    ind1b = 1;
                    ind2b = lesLocsC(i,j) + 34;
                elseif ((lesLocsC(i,j)+20) > size(C,2))
                    ind1 = lesLocsC(i,j) - 5;
                    ind2 = size(C,2);
                    ind1b = lesLocsC(i,j) - 15;
                    ind2b = size(C,2);
                else
                    ind1 = lesLocsC(i,j)-5;
                    ind2 = lesLocsC(i,j)+19;
                    if(lesLocsC(i,j) <= 15)
                        ind1b = 1;
                    else
                        ind1b = lesLocsC(i,j) - 15;
                    end
                    if((lesLocsC(i,j)+34) > size(C,2))
                        ind2b = size(C,2);
                    else
                        ind2b = lesLocsC(i,j) + 34;
                    end
                end
                tempCorr = corr(C(i,ind1:ind2)',C_raw(i,ind1:ind2)'); %calculate the correlation between the C and C_raw peaks
                PairCorr(i,j) = tempCorr;
                tempBL = C_raw(i,ind1b:ind2b) - C(i,ind1b:ind2b); %calculate the BL in a local window
                stdt = std(tempBL);
                LocalPNR(i,j) = lesPeaksC(i,j)./stdt; %local PNR is defined as the PNR only taking the BL locally to the peak, and not the overall trace BL
            end
        end
    end

    % Activate the following line to add LocalPNR and PairCorr to
    % your saved variable
    PeakProp.LocalPNR = LocalPNR;
    PeakProp.PairCorr = PairCorr;

    % Buidling the data matrix for subsequent PCA transformation.
    % This simply rearranges the data of all the peaks from all neurons in a 
    % single matrix. 
    nobs = size(find(indexC > 0),1);
    Data = zeros(nobs,8);
    countP = 1;
    group = zeros(1,nobs);
    for i = 1:numb
        len = size(find(indexC(i,:)>0),2);
        indyC = find(indexC(i,:)>0);
        indyRaw = find(indexRaw(i,:)>0);
        for j = 1:len
            Data(countP,1) = PairCorr(i,indyC(j));
            Data(countP,2) = lesPNR(i,indyRaw(j));
            Data(countP,3) = lesPNRC(i,indyC(j));
            Data(countP,4) = lesWHR(i,indyRaw(j));
            Data(countP,5) = lesWHR_C(i,indyC(j));
            Data(countP,6) = lesProm(i,indyRaw(j));
            Data(countP,7) = lesPromC(i,indyC(j));
            Data(countP,8) = LocalPNR(i,indyC(j));
            countP = countP + 1;
        end
    end

    % Here labels for the input CNMF-E peaks are estimated using one or two
    % models. The PCA model is estimated from the first two parameters of a PCA
    % analysis of the training set. The PNR uses the parameters 'PNR' and
    % 'LocalPNR' as the classificaiton parameters. Prediction on the new set
    % are saved in 'labelPCA' and 'labelPNR'.
    load('TrSet'); %Load the training set that has been used to fix the classification parameters
    PCAtf = TrSet.PCAres;
    MdlPCA = TrSet.MdlPCA;
    coeff = PCAtf.coeff;
    mu = PCAtf.mu;
    scoreTest = (inv(coeff)*(Data-ones(size(Data,1),1)*mu)')'; %This calculates the PCA of the dataset using the parameters from the training set.
    PC1Test = scoreTest(:,1)';
    PC2Test = scoreTest(:,2)'; %The 2 first PC are sufficient to classify the peaks. We only perform a two-dimensional classification.
    [labelPCA, ~] = predict(MdlPCA, [PC1Test', PC2Test']); %classify the peaks as 'good' or 'bad' according to the model parameters.

    % Activate the line below to save the label of the classification
    PeakProp.labelPred = labelPCA;

    % The following part clean the unpaired and low quality peaks from the
    % dataset. It thens register all 'good' peaks and their properties in
    % a new vector. 
    peakNum = 1;
    deleteNum = 0;
    locsFilt = zeros(size(lesLocs));
    locsCFilt = zeros(size(lesLocsC));
    peaksFilt = zeros(size(lesPeaks));
    peaksCFilt = zeros(size(lesPeaksC));
    PNRFilt = zeros(size(lesPNR));
    PNRCFilt = zeros(size(lesPNRC));
    promFilt = zeros(size(lesProm));
    promCFilt = zeros(size(lesPromC));
    WHRFilt = zeros(size(lesWHR));
    WHRCFilt = zeros(size(lesWHR_C));
    LocalPNRFilt = zeros(size(LocalPNR));
    for i = 1:numb
        %Initialization of the arguments of the loop
        clear toDelete locsTmp peaksTmp PNRTmp promTmp WHRTmp LocalPNRTmp;
        peakTmp = peakNum;
        len = size(lesLocsC,2)-size(find(lesLocsC(i,:)==0),2);
        toDelete = zeros(1,len);
        locsTmp = lesLocsC(i,:);
        peaksTmp = lesPeaksC(i,:);
        PNRTmp = lesPNRC(i,:);
        promTmp = lesPromC(i,:);
        WHRTmp = lesWHR_C(i,:);
        LocalPNRTmp = LocalPNR(i,:);
        %Removal of 'bad' and 'unpaired' neurons
        for j = 1:len
            if(indexC(i,j) == -1) %if peak is marked as 'unpaired' --> delete
                toDelete(j) = 1;
            elseif(indexC(i,j) > 0)
                if(labelPCA(peakNum) == 1) %if peak is marked as 'low quality' --> delete
                    toDelete(j) = 1;
                    peakNum = peakNum + 1;
                elseif(labelPCA(peakNum) == 3)
                    peakNum = peakNum + 1;
                end
            end
        end
        idDel = find(toDelete == 1);
        numDel = size(idDel,2);
        filler = zeros(1,numDel);
        % Registration of the good peaks' parameters
        locsTmp(idDel) = []; locsTmp = [locsTmp, filler]; locsCFilt(i,:) = locsTmp;
        peaksTmp(idDel) = []; peaksTmp = [peaksTmp, filler]; peaksCFilt(i,:) = peaksTmp;
        PNRTmp(idDel) = []; PNRTmp = [PNRTmp, filler]; PNRCFilt(i,:) = PNRTmp;
        promTmp(idDel) = []; promTmp = [promTmp, filler]; promCFilt(i,:) = promTmp;
        WHRTmp(idDel) = []; WHRTmp = [WHRTmp, filler]; WHRCFilt(i,:) = WHRTmp;
        LocalPNRTmp(idDel) = []; LocalPNRTmp = [LocalPNRTmp, filler]; LocalPNRFilt(i,:) = LocalPNRTmp;

        %The following lines repeat the previous operations but in the C_raw trace
        clear toDelete locsTmp peaksTmp PNRTmp promTmp WHRTmp LocalPNRTmp
        peakNum = peakTmp;
        len2 = size(lesLocs,2)-size(find(lesLocs(i,:)==0),2);
        toDelete = zeros(1,len2);
        locsTmp = lesLocs(i,:);
        peaksTmp = lesPeaks(i,:);
        PNRTmp = lesPNR(i,:);
        promTmp = lesProm(i,:);
        WHRTmp = lesWHR(i,:);
        for k = 1:len2
            if(indexRaw(i,k) == -1)
                toDelete(k) = 1;
            elseif(indexRaw(i,k) > 0)
                if(labelPCA(peakNum) == 1)
                    toDelete(k) = 1;
                    peakNum = peakNum + 1;
                elseif(labelPCA(peakNum) == 3)
                    peakNum = peakNum + 1;
                end
            end
        end
        idDel = find(toDelete == 1);
        numDel = size(idDel,2);
        filler = zeros(1,numDel);
        locsTmp(idDel) = []; locsTmp = [locsTmp, filler]; locsFilt(i,:) = locsTmp;
        peaksTmp(idDel) = []; peaksTmp = [peaksTmp, filler]; peaksFilt(i,:) = peaksTmp;
        PNRTmp(idDel) = []; PNRTmp = [PNRTmp, filler]; PNRFilt(i,:) = PNRTmp;
        promTmp(idDel) = []; promTmp = [promTmp, filler]; promFilt(i,:) = promTmp;
        WHRTmp(idDel) = []; WHRTmp = [WHRTmp, filler]; WHRFilt(i,:) = WHRTmp;
    end

    % Activate the lines below to save the filtered peaks and their properties
    PeakProp.locsFilt = locsFilt;
    PeakProp.peaksFilt = peaksFilt;
    PeakProp.peaksCFilt = peaksCFilt;
    PeakProp.PNRFilt = PNRFilt;
    PeakProp.PNRCFilt = PNRCFilt;
    PeakProp.WHRFilt = WHRFilt;
    PeakProp.WHRCFilt = WHRCFilt;
    PeakProp.promFilt = promFilt;
    PeakProp.promCFilt = promCFilt;
    PeakProp.LocalPNRFilt = LocalPNRFilt;

    %Register all the final peaks' parameters in one structure
    PeaksFinal = struct('locs',locsFilt,'peaks', peaksFilt,'PNR',PNRFilt,'WHR',WHRFilt,'prom',promFilt,'localPNR',LocalPNRFilt);
    %This saves the peak properties
    peakName = oname{n} + "_peakProps";
    save(peakName,'PeaksFinal');

    

    %% Part 4 - Activity-based neuron cleaning 
    % In this part neurons are further sorted based on their activity. 
    % It uses the cleared peaks from Part 3 to determine what neurons to keep
    % and what neurons to discard. Neurons will be classified according to 2
    % parameters:
    % 1) the number of peaks. Neurons with less than 2 peaks will be thrown
    % away.
    % 2) the PNR value of the highest peak. If the highest peak has a low PNR
    % value, the neuron is probably of low quality and is discarded.
    % The selection is made following the following equation, a neuron is
    % discarded if:
    % Number of peaks + 2*PNR(maxPeak) - 3 < 0

    % Calculation of the PNR of the maximum peak for every neuron
    ss = size(find(PeaksFinal.peaks > 0),1);
    meanPNR = sum(sum(PeaksFinal.PNR))/ss;
    neuData = zeros(numb,1);
    for i = 1:numb
        ssTmp = size(find(PeaksFinal.peaks(i,:)>0),2);
        %The loop only runs if there is at least one peak that is not equal to
        %0
        if ssTmp > 0
            meanPNRTmp = sum(PeaksFinal.PNR(i,:))/ssTmp;
            maxPNRTmp = max(PeaksFinal.PNR(i,:));
            neuData(i,1) = maxPNRTmp/meanPNR;
        end
    end

    % This registers all neurons that does not qualify to the activity equation
    % to be deleted.
    activities = zeros(numb, 1);
    IndexDynamic = [];
    numPks = zeros(1, numb);
    for i = 1:numb
        numPks(i) = size(find(PeaksFinal.peaks(i,:)>0),2);
        if(numPks(i)+2*neuData(i,:)-3 < 0)
            IndexDynamic = [IndexDynamic i];
        end
    end

    %% Part 5 - Cleaning overlapping neurons
    % This part resolves case of overlapping neurons. First, overlapping
    % candidates are selected on the basis of the proximity of their center of
    % mass. The overlap of their area is then used to determine what neurons
    % overlap too much and should be cleaned. When two neurons are found to
    % overlap, the average value of the PNR of their 3 highest peak is used to
    % determine the better quality neuron that should be retained and the lower
    % quality neuron that should be discarded. 

    %load('spatCorrCA1')%This matrix calculates the overlap % between every pair of neurons
    %spatCorrMat = spatCorr
    d_min = 36; %Maximum distance between the centers of mass of two neurons to
                % be considered a potential overlapping pair
    threshOverlap = 0.5; %If 2 neurons overlap more than the threshold, one of them is discarded 
    IndexOverlap = [];

    % The average PNR value of the 3 highest peaks is estimated for every
    % neuron even though not all them are going to be tested. It is however not
    % computationally heavy.
    max3PNR = zeros(1,numb);
    ss = size(find(PeaksFinal.peaks > 0),1);
    meanPNR = sum(sum(PeaksFinal.PNR))/ss;
    for i = 1:numb
        ssTmp = size(find(PeaksFinal.peaks(i,:)>0),2);
        if ssTmp > 2 
            meanPNRTmp = sum(PeaksFinal.PNR(i,:))/ssTmp;
            maxPNRTmp = max(PeaksFinal.PNR(i,:));
            [valMx indMx] = sort(PeaksFinal.PNR(i,:),'descend');
            pk1 = valMx(1)/meanPNR;
            pk2 = valMx(2)/meanPNR;
            pk3 = valMx(3)/meanPNR;
            mean3PNR = mean([pk1, pk2, pk3]);
        % If the neuron does not have 3 peaks, only the 2 peaks are being used
        % to determine the averaged PNR
        elseif ssTmp == 2 
            meanPNRTmp = sum(PeaksFinal.PNR(i,:))/ssTmp;
            maxPNRTmp = max(PeaksFinal.PNR(i,:));
            [valMx indMx] = sort(PeaksFinal.PNR(i,:),'descend');
            pk1 = valMx(1)/meanPNR;
            pk2 = valMx(2)/meanPNR;
            mean3PNR = mean([pk1, pk2]);
        end
        max3PNR(i) = mean3PNR;
    end

    % The distance between the neuron centers is calculated across the entire
    % dataset
    ctr = neuron.estCenter();
    yy = ctr(:,1); 
    xx = ctr(:,2);
    dist_v = sqrt(bsxfun(@minus, xx, xx').^2 + bsxfun(@minus, yy, yy').^2); %calculates the distance between every neuron
    dist_v = tril(dist_v); %only takes the lower triangular matrix (because of redundant information)
    [ind2, ind1] = find((dist_v <= d_min)&(dist_v > 0)); %detects close neurons but ignores the zero in the matrix
    closePairs = zeros(size(dist_v));

    % Close candidates are being further tested. If their area overlaps more
    % than the authorized threshold, the one with the weaker average PNR is
    % discarded
    Coor = neuron.get_contours;
    IndexOverlap = [];
    for i = 1:size(ind1,1)
        cont1 = Coor{ind1(i)};
        cont2 = Coor{ind2(i)};
        p1 = polyshape(cont1(1,:), cont1(2,:)); %polyshape is a type of object that represents a shape (polygon)
        p2 = polyshape(cont2(1,:), cont2(2,:));
        pInt = intersect(p1,p2); %this function claculates the polygon that is contained in both polygons
        ratio1 = area(pInt)/area(p1);
        ratio2 = area(pInt)/area(p2);
        % The mean of the overlapping ratio is taken as the total overlapping
        % percentage
        if mean([ratio1,ratio2]) >= threshOverlap
            PNRdiff = max3PNR(ind1(i)) - max3PNR(ind2(i));
            % tests what neuron has the higher average PNR
            if PNRdiff >= 0
                IndexOverlap = [IndexOverlap, ind2(i)];
            else
                IndexOverlap = [IndexOverlap, ind1(i)];
            end
        end
    end

    %% Part 6 - Wrinkly neurons
    % The final cleaning operation is to detect neurons with a wrinkly contour.
    % Although the roughness of their contour is not a reason to consider as
    % bad, this feature is often found in neurons with a bad signal.

    %Reset contours of the neuron
    Coor = [];
    neuron.Coor = [];
    Coor = neuron.get_contours;
    neuron.Coor = Coor;

    % Extraction of the perimeter and area to calculate the compacity
    peri = zeros(1,numb);
    areaa = zeros(1,numb);
    ratioz = zeros(1,numb);
    for i = 1:numb
        if sum(neuron.C(i,:) > 0.01)
            contemp = neuron.Coor{i};
            x = contemp(1,:)';
            y = contemp(2,:)';
            shptemp = alphaShape(x,y,2);
            pertemp = perimeter(shptemp);
            areatemp = area(shptemp);
            peri(i) = pertemp;
            areaa(i) = areatemp;
            ratiotemp = (pertemp*pertemp)/(4*pi*areatemp);
            ratioz(i) = ratiotemp;
        end
    end

    threshCompa = 35; %Set the treshold of the compacity score for neuron selection
    IndexCompa = [];
    for i = 1:numb
        if ratioz(i) <= threshCompa
            IndexCompa = [IndexCompa i];
        end
    end

    % Assembles the 4 stages of neuron classification
    FinInd = zeros(1, numb);
    indouille = [IndexShape, IndexDynamic, IndexCompa, IndexOverlap];
    for i = 1:numb
        if isempty(find(indouille==i))
             FinInd(i) = 3;% good neuron
         else
             FinInd(i) = 1; % bad neurons
         end
    end

    % Deletes the neurons' activity traces and location. 
    % Note that removing the neurons from the set is not preferred because it
    % changes the code of every neuron in the dataset. 
    for i = 1:numb
        if FinInd(i) == 1;
            neuron.C(i,:) = 0;
            neuron.C_raw(i,:) = 0;
            neuron.A(:,i) = 0;
        end
    end
    C = neuron.C;
    C_raw = neuron.C_raw;
    A = neuron.A;
    
    neuronClean = {'C', C, 'C_raw', C_raw, 'A', A};

    % Save C, C_raw, and C with deleted neurons replaced by zeros
    save(oname{n},'neuronClean');
    clear neuron numb siz Coor CNMFE_Data Data IndexOverlap IndexShape IndexCompa IndexDynamic FinInd neuronClean;
end