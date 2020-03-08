clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Generate dFoverF curve from raw data, Extract and output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reference:In vivo two-photon imaging of sensory-evoked dendritic calcium
%signals in cortical neuron, Nature Protocols
%Input: Directory containing Trials and Stims xlsx files.
%Output:
% dFoF_TrialSorted;
% dFoF_FreqSorted;
% dFoF_ROISorted;
% dFoF_ROISorted_Norm;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add path
Main_dir = fileparts(which('dFoF_Cal_Sort.m'));
addpath(fullfile(Main_dir, 'Plot_functions'));

maindir = uigetdir;
totaloutdir = strcat(maindir, '\ExtractOutput');
mkdir(totaloutdir);


%% User defined parameters

% Output files: 

% Filter Option: 1.Select Sound reponsive curve to plot; 0. Plot all curve
% without filtering
FilterOn = 1;


% Default parameters of input files;
framerate = 5;  %Hz;
AVGWindow = [-0.4 10]; % unit: s;          AVG(1) should be <= 0;
baseFreq = 3;
Interval = AVGWindow(2) - AVGWindow(1); % unit: s ; should be smaller than stimulation interval;

% Parameters related to matrix dimensions;
nROIs = 3;
FreqNum = 25;   % 3K~48K totally 25 frequencies: freq = baseFreq*2^((i-1)/6); i = 1 + 6*log2(freq/baseFreq);
nFrames  = 600;
minrows = Interval*framerate;
nStimTrial = 25;
nStimAll = 20;   % Max stim num of a single frequency among all trials; 

% Time Constants
dt = 1/framerate;    %Frame interval/0.2s;
t0 = 2*dt;   %Time constant for noise filtering using exponentially weighted moving averaged.
t1 = 3*dt;   %Time window for smoothing;
t2 = 50;     %Time window for finding the baseline F0

BaseTWin = round(t2/dt)-1;

savepath0 = strcat(totaloutdir,'\Stims_TrialSorted.mat');
savepath1 = strcat(totaloutdir,'\dFoF_TrialSorted.mat');
savepath2 = strcat(totaloutdir,'\dFoF_FreqSorted.mat');
savepath3 = strcat(totaloutdir,'\dFoF_ROISorted.mat');
savepath4 = strcat(totaloutdir,'\dFoF_ROISorted_Norm.mat');



subdirpath = fullfile( maindir, '*.xlsx' );
dat = dir( subdirpath );               
nTrials = length( dat )/2;


%% Key matrix 

% input matrix of fluorescence data: nFrames * nROIs * nTrials;
inFtRaw = nan(nFrames,nROIs,nTrials);
% input matrix of stimulation sequence: 
%  nStims * 1 * nTrials -- Stim frame,  nStims * 2 * nTrials -- Stim frequency;
inStimSeq = nan(nStimTrial,2,nTrials);

% Output matrix
dFoF_TrialSorted = nan(nFrames,nROIs,nTrials);   %2nd dim: nROI, 3rd dim: nFrames;

dFoF_FreqSorted = nan(minrows,nROIs,nStimAll,FreqNum); 

dFoF_ROISorted = nan(minrows,nStimAll,FreqNum,nROIs);  

dFoF_ROISorted_Norm = nan(minrows,nStimAll,FreqNum,nROIs); 

tseries = transpose((AVGWindow(1):dt:AVGWindow(2)-dt));

%% Section 0: read files into matrix;
for i = 1 : nTrials
    TrialFile = strcat(maindir,'\Trial',num2str(i),'.xlsx');
    StimFile = strcat(maindir, '\Stim',num2str(i),'.xlsx');
    
    %open excel file;
    %tmpFtRaw = xlsread(TrialFile); % raw data; first column is time sequence; the others are F of differenct ROIs;
    %system('taskkill /F /IM EXCEL.EXE');
    tmpFtRaw = readmatrix(TrialFile); % raw data; first column is time sequence; the others are F of differenct ROIs;

    inFtRaw(:,:, i) = tmpFtRaw(:,2:end);    
    
    tmpStimSeq = readmatrix(StimFile); % stim data; first row is time sequence of stimulation; the other is audio frequency;
    %system('taskkill /F /IM EXCEL.EXE');
    
    num = size(tmpStimSeq,1);
    inStimSeq(1:num,1,i) = round((tmpStimSeq(1:num,1) + AVGWindow(1)).*framerate + 1);
    inStimSeq(1:num,2,i) = tmpStimSeq(1:num,2);     
end


%% Section 1: Process each data set and output "dFoF_TrialSorted" and "dFoF_FrequencySorted";
for N = 1 : nTrials
    
    [n, m, trow] = size(inFtRaw(:,:,N)); % n = nFrames; m = nROIs;
    
    Ft_Sub(:,:) = transpose(inFtRaw(:,:,N));    % [m,n] dimension, array containing F of different ROIs;
    
    %% smoothed version of Ft by averaging within a time window; time window equals 3dt=0.6s in this case;
    Ft_Sm = zeros(m,n);   % smoothed version of Ft during as time window before frame j;
    Ft_Sm(:,1) = Ft_Sub(:,1);
    Ft_Sm(:,n) = Ft_Sub(:,n);
    
    
    for j = 2: (n - 1)
        Ft_Sm(:,j) = ( Ft_Sub(:,j-1) + Ft_Sub(:,j)+  Ft_Sub(:,j+1))/3;
    end
    
    %% Baseline F0(t) is the average among the 40% smallest value of the smooth Ft_Sm within a time window before t; time window equals 20dt= 4s
    F0_t = zeros(m,n);
    
    for j = 1:BaseTWin
        lwsort = sort(Ft_Sm(:,1:j),2);                            % sorted by ascending;
        F0_t(:,j) = mean(lwsort(:,1:1+round(0.4*j)),2);
    end
    for j = BaseTWin + 1:n
        lwsort = sort(Ft_Sm(:,j-BaseTWin:j),2);
        F0_t(:,j) = mean(lwsort(:,1:1+round(0.4*BaseTWin)),2);    % the average of 40% smallest values;
    end
    
    %% Calculate relative change of fluorescence signal Rt;
    Rt = (Ft_Sub - F0_t)./F0_t;
    
    %% Apply noise filtering(exponentially weighted moving average: EWMA) to get final result dFoF;
    Wt = zeros(n);
    for p = 1:n
        Wt(p) = exp(-p*dt/t0);
    end
    
    dFoF  = zeros(m,n);
    SumWt = zeros(n);
    for t = 1:n
        for tau = 1:t
            dFoF(:,t) = dFoF(:,t) + Rt(:,t-tau+1).*Wt(tau);
            SumWt(t) = SumWt(t) + Wt(tau);
        end
        
        dFoF(:,t) = dFoF(:,t)./SumWt(t);
    end
    
    %% Save dFoF_t Curve
    dFoF_TrialSorted(:,:,N) =  transpose(dFoF);
    
    %% Extract data for different stimulation
    tmp = inStimSeq(:,1,N);
    tlength1 = length(tmp(~isnan(tmp))); %length of non-nan rows;
    
    for i = 1:tlength1
        % Extract data with different frequency;
        dFoF_Seq = dFoF(:,inStimSeq(i,1,N):(inStimSeq(i,1,N) + minrows - 1));
        
        tStimFreq = inStimSeq(i,2,N);
        tFreqNum =  round(1 + 6*log2(tStimFreq/baseFreq));
        % assign the dFoF_Seq to first j of Non-NaN matrix;
        for j = 1: size(dFoF_FreqSorted,3)
            B = isnan(dFoF_FreqSorted(:,:,j,tFreqNum));
            if all(B == 1)
                dFoF_FreqSorted(:,:,j,tFreqNum) = transpose(dFoF_Seq);
                break;
            end
        end
    end
    
end



%% Section 2: Generate and output "dFoF_ROISorted" based on "dFoF_FreqSorted";   
for i = 1:FreqNum
    for m = 1:nROIs
        
        % find the number of nStims that have values, the following are
        % NaN;
        tmp1 = dFoF_FreqSorted(1,m,:,i);
        tlength1 = length(tmp1(~isnan(tmp1)));
        if tlength1 > 0                     % non all rows are NaN for tStim;
            for n = 1:tlength1
                %% Selected the sound responsive curve for plotting
                if FilterOn == 1
                    %determine whether the neurons is tone responsive using
                    %peak detection; Criteria: the maxmum peak should be
                    %within 2.0s after the stimulation;
                    tmp2 = dFoF_FreqSorted(:,m,n,i);
                    tmp2 = smooth(tmp2,5);
                    
                    [pks, locs] = findpeaks(tmp2,'SortStr','descend');
                    
                    % peak time should be [0.5, 3s]; peak value should equal to maximum;
                    if ~isempty(locs) && locs(1) > (0.5-AVGWindow(1))*framerate && locs(1) < (3-AVGWindow(1))*framerate && pks(1)== max(tmp2) %&& pks(1) < threshold 
                        tmp3 = dFoF_ROISorted(1,:,i,m);
                        tcol = length(tmp3(~isnan(tmp3)));
                        dFoF_ROISorted(:,tcol + 1,i,m) = dFoF_FreqSorted(:,m,n,i);
                        
                    end
                    
                    %% Plot all curve without filtering
                elseif FilterOn == 0
                    tmp3 = dFoF_ROISorted(1,:,i,m);
                    tcol = length(tmp3(~isnan(tmp3)));
                    dFoF_ROISorted(:,tcol + 1,i,m) = dFoF_FreqSorted(:,m,n,i);
                end
            end
        end
    end
end



%% Normalize the tuning curve to the largest amplitude of each ROI;
for m = 1:nROIs
    % Normalized to max vaule of the mean reponse curve of different frequency;
    AVGdFoF =  mean(dFoF_ROISorted(:,:,:,m),2,'omitnan');
    dFoF_ROISorted_Norm(:,:,:,m) = dFoF_ROISorted(:,:,:,m)./max(max(AVGdFoF));    
end

save(savepath0, 'inStimSeq','-v7.3','-nocompression');
save(savepath1, 'dFoF_TrialSorted','-v7.3','-nocompression');
save(savepath2, 'dFoF_FreqSorted','-v7.3','-nocompression');
save(savepath3, 'dFoF_ROISorted','-v7.3','-nocompression');
save(savepath4, 'dFoF_ROISorted_Norm','-v7.3','-nocompression');

%% Section 3: Save the data into excel files (Optional); Note: this will take a lot of time;  

% savepath5 = strcat(totaloutdir,'\dFoF_TrialSorted.xls');
% for i = 1: nTrials
%     strTrial = strcat('Trial',num2str(i,'%02d'));
%     writematrix(dFoF_TrialSorted(:,:,i),savepath5,'Sheet',i);
% end