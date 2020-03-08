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
ResponseWindow = [-2 4]; % unit: s;          0 corresponds to the time of peak amplitude.
baseFreq = 3;
Interval = ResponseWindow(2) - ResponseWindow(1); % unit: s ; should be smaller than stimulation interval;
PeakThreshold = 0.4;                              %dFoF threshold for peak detection;


% Parameters related to matrix dimensions;
nROIs = 3;
%FreqNum = 25;   % 3K~48K totally 25 frequencies: freq = baseFreq*2^((i-1)/6); i = 1 + 6*log2(freq/baseFreq);
nFrames  = 600;
minrows = Interval*framerate;
nStimTrial = 25;
nResonpseAll = 20;   % Max stim num of a single frequency among all trials; 




% Time Constants
dt = 1/framerate;    %Frame interval/0.2s;
t0 = 2*dt;   %Time constant for noise filtering using exponentially weighted moving averaged.
t1 = 3*dt;   %Time window for smoothing;
t2 = 40;     %Time window for finding the baseline F0

BaseTWin = round(t2/dt)-1;

%savepath0 = strcat(totaloutdir,'\Stims_TrialSorted.mat');
savepath1 = strcat(totaloutdir,'\dFoF_TrialSorted.mat');
savepath2 = strcat(totaloutdir,'\RespSeg_TrialSorted.mat');
savepath3 = strcat(totaloutdir,'\Amplitude_TrialSorted.mat');
savepath4 = strcat(totaloutdir,'\RespSeg_ROISorted.mat');
savepath5 = strcat(totaloutdir,'\Struct_Amplitude_TrialSorted.mat');
savepath6 = strcat(totaloutdir,'\Struct_RespSeg_ROISorted.mat');

subdirpath = fullfile( maindir, '*.xlsx' );
dat = dir( subdirpath );               
nTrials = length( dat );


%% Key matrix 

% input matrix of fluorescence data: nFrames * nROIs * nTrials;
inFtRaw = nan(nFrames,nROIs,nTrials);
% input matrix of stimulation sequence: 
%  nStims * 1 * nTrials -- Stim frame,  nStims * 2 * nTrials -- Stim frequency;
inStimSeq = nan(nStimTrial,2,nTrials);

% Output matrix
dFoF_TrialSorted = nan(nFrames,nROIs,nTrials);   %2nd dim: nROI, 3rd dim: nFrames;

RespSeg_TrialSorted = nan(minrows,nROIs,nResonpseAll,nTrials); 

Amplitude_TrialSorted = nan(nROIs,nResonpseAll,nTrials); 

RespSeg_ROISorted = nan(minrows,nTrials*nResonpseAll,nROIs);  

tseries = transpose((ResponseWindow(1):dt:ResponseWindow(2)-dt));

%% Section 0: read files into matrix;
for i = 1 : nTrials
    
    TrialFile = strcat(maindir,'\Trial',num2str(i),'.xlsx');
    tmpFtRaw = readmatrix(TrialFile); % raw data; first column is time sequence; the others are F of differenct ROIs;
    inFtRaw(:,:, i) = tmpFtRaw(:,2:end);    
    
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
    
    %% Extract data for different response segment [peak-2s, peak+4s]  
    tmp2 = dFoF_TrialSorted(:,1,N);
    tmp2 = smooth(tmp2,5);
    [pks, locs] = findpeaks(tmp2);
      
    % find the response seg based on the peak of the first ROI;
    k=0;
    for i = 1:length(locs)       
        if (locs(i)>10)&&(locs(i)<nFrames-20)&&(pks(i)>PeakThreshold)    %dFoF > 1;
            k=k+1;
            dFoF_Seq = dFoF(:,(locs(i)-10):(locs(i)+19));   %[peaktime - 2s, peaktime+4s]
            RespSeg_TrialSorted(:,:,k,N) = transpose(dFoF_Seq);
            dFoF_Seq_temp = dFoF(:,(locs(i)-5):(locs(i)+5));
            Amplitude_TrialSorted(:,k,N) = max(dFoF_Seq_temp,[],2);
            
        end
    end
    
end

nResponse = size(RespSeg_TrialSorted,3);

for a=1: nROIs   
    x = 1;
    for b=1:nTrials
        for c=1:nResponse            
            if ~isnan(RespSeg_TrialSorted(:,a,c,b))
                RespSeg_ROISorted(:,x,a) =  RespSeg_TrialSorted(:,a,c,b);
                x = x+1;
            end
        end
    end
end


for num1 = 1:nTrials    
    Amplitude.('Trial_'+string(num1)) = Amplitude_TrialSorted(:,:,num1);
end

for num2 = 1:nROIs    
    RespSeg.('ROI_'+string(num2)) = RespSeg_ROISorted(:,:,num2);
end


%save(savepath0, 'inStimSeq','-v7.3','-nocompression');
save(savepath1, 'dFoF_TrialSorted','-v7.3','-nocompression');
save(savepath2, 'RespSeg_TrialSorted','-v7.3','-nocompression');
save(savepath3, 'Amplitude_TrialSorted','-v7.3','-nocompression');
save(savepath4, 'RespSeg_ROISorted','-v7.3','-nocompression');
save(savepath5,'Amplitude','-v7.3','-nocompression');
save(savepath6,'RespSeg','-v7.3','-nocompression');
%% Section 3: Save the data into excel files (Optional); Note: this will take a lot of time;  

% savepath5 = strcat(totaloutdir,'\dFoF_TrialSorted.xls');
% for i = 1: nTrials
%     strTrial = strcat('Trial',num2str(i,'%02d'));
%     writematrix(dFoF_TrialSorted(:,:,i),savepath5,'Sheet',i);
% end