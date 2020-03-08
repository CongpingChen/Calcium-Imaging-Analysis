clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Calculate dFoF amplitude for all Stims of different Trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input: Directory containing "dFoF_StimSorted.mat" file.
%Output: "Amplitude_ROISorted.mat",amplitude for all stims; and
%"PeakTime_ROISorted", peaktime for all stims;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maindir = uigetdir;


%% Default parameters of input files;
framerate = 5;  %Hz;
AVGWindow = [-1 7];  % unit: s;
Interval = AVGWindow(2) - AVGWindow(1);        % unit: s ; smaller than stimulation interval;
minrows = Interval*framerate;

%% Parameters related to matrix dimensions;
nTrials = 25;
nROIs = 20;
FreqNum = 25;   % 3K~48K totally 25 frequencies: freq = baseFreq*2^((i-1)/6); i = 1 + 6*log2(freq/baseFreq);
nFrames  = 600;

dt = 1/framerate;

nStimTrial = 25;
nStimAll = 20;   % Max stim num of a single frequency among all trials; 

%% Calculate Option: 1.Response value; Version 2: Local maximum; Version 3: Maximum during the responsive window.
MaxOption = 2;

subdirpath = fullfile( maindir, '*.xlsx' );
dat = dir( subdirpath );   

savepath1 = strcat(maindir,'\Amplitude_ROISorted.mat');
savepath2 = strcat(maindir,'\PeakTime_ROISorted.mat');

%% Key matrix

% input matrix
% dFoF_StimSorted = nan(minrows,nROIs,nStimAll,nTrials); 

% Output matrix
Amplitude_ROISorted = nan(nTrials,nStimTrial,nROIs);   
PeakTime_ROISorted = nan(nTrials,nStimTrial,nROIs); 
%dFoF_ROISorted = nan(minrows,nStimAll,nTrials, nROIs);


%tseries = transpose((0:1:minrows-1)*dt);
tseries = transpose((AVGWindow(1):dt:AVGWindow(2)-dt));


%% Start processing;
filepath = strcat(maindir,'\dFoF_ROISorted.mat');

for i = 1:nROIs
    for j = 1: nTrials
        
        dFoF_ROISorted = importdata(filepath);
        tmp0 = dFoF_ROISorted(1,:,j,i);
        tlength = length(tmp0(~isnan(tmp0))); %length of non-nan rows;
        
        for k = 1:tlength
            
            tmp1 = Amplitude_ROISorted(j,:,i);
            tcol1 = length(tmp1(~isnan(tmp1)));
            
            tmp2 = PeakTime_ROISorted(j,:,i);
            tcol2 = length(tmp2(~isnan(tmp2)));
            
            %% calculate max amplitude;
            switch (MaxOption)
                case 1
                    %% Version 1: Response value
                    % amplitude equals the peak value (1:3) s minus average value (0:1) s;
                    Res_Value = mean(dFoF_ROISorted((1:3)*framerate,k,j,i), 1) - mean(dFoF_ROISorted(1:framerate,k,j,i), 1);
                    
                    % assign the dFoF_Seq to first j of Non-NaN matrix;
                    if Res_Value > 0
                        Amplitude_ROISorted(j,tcol1+1,i) = max(dFoF_ROISorted((1:3)*framerate,k,j,i),[],1);
                    else
                        Amplitude_ROISorted(j,tcol1+1,i) = nan;
                        
                    end
                case 2
                    %% Version 2: Local maximum
                    % amplitude is local maximum of stimulation-reponsive
                    % curves
                    tempb = transpose(dFoF_ROISorted(:,k,j,i));
                    tempb = smooth(tempb,5);
                    [pks, locs] = findpeaks(tempb,'SortStr','descend');
                    
                    if isempty(locs)
                        Amplitude_ROISorted(j,tcol1+1,i) = nan;                                                     % not ES responsive, amplitude is 0;
                    elseif locs(1) > 2 && locs(1) < 2*framerate && pks(1)== max(tempb)                             % peak time should be [1s, 3s]; peak value should equal to maximum;
                        Amplitude_ROISorted(j,tcol1+1,i) = max(dFoF_ROISorted(1:4*framerate,k,j,i),[],1);            % ES responsive;
                    else
                        Amplitude_ROISorted(j,tcol1+1,i)= nan;                                                      % not tone responsive neuron, set as 0;
                    end
                case 3
                    %% Version 3: Maximum during the responsive window
                    % amplitude equals the peak value (1-21 frame) minus average value (1-3 frame);
                    Amplitude_ROISorted(j,tcol1+1,i) = max(dFoF_ROISorted(1:4*framerate,k,j,i),[],1);
            end
            
            %% calculate peak time;
            tempb = transpose(dFoF_ROISorted(:,k,j,i));
            tempb = smooth(tempb,5);
            [pks, locs] = findpeaks(tempb,'SortStr','descend');
            if isempty(locs)
                PeakTime_ROISorted(j,tcol2+1,i) = nan;                                                     % not ES responsive, amplitude is 0;
            elseif locs(1) > 2 && locs(1) < 4*framerate && pks(1)== max(tempb)                           % peak time should be [1s, 4s]; peak value should equal to maximum;
                PeakTime_ROISorted(j,tcol2+1,i) = tseries(locs(1));                                          % peak time;
            else
                PeakTime_ROISorted(j,tcol2+1,i) = nan;                                                     % not tone responsive neuron, set as 0;
            end
        end
    end
end


save(savepath1, 'Amplitude_ROISorted','-v7.3','-nocompression');
save(savepath2, 'PeakTime_ROISorted','-v7.3','-nocompression');

savepath3 = strcat(maindir,'\Amplitude_ROISorted_Mo',num2str(MaxOption),'.xlsx');
savepath4 = strcat(maindir,'\PeakTime_ROISorted_Mo',num2str(MaxOption),'.xlsx');

%% Output to excel files;
for i = 1: nROIs
    strROI = strcat('ROI_',num2str(i));    
    xlswrite(savepath3,Amplitude_ROISorted(:,:,i),strROI);
    system('taskkill /F /IM EXCEL.EXE');
    xlswrite(savepath4,PeakTime_ROISorted(:,:,i),strROI);
    system('taskkill /F /IM EXCEL.EXE');  
end