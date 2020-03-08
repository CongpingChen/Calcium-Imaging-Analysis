clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Generate dFoverF curve from raw data, Extract and output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reference:In vivo two-photon imaging of sensory-evoked dendritic calcium
%signals in cortical neuron, Nature Protocols
%Input: Directory containing Trials and Stims xml files.
%Output: Sorted files based on trials and frequency;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maindir = uigetdir;

%% Default parameters of input files;
framerate = 5;  %Hz;
AVGWindow = [-1 9]; % unit: s;          AVG(1) should be <= 0;
Interval = AVGWindow(2) - AVGWindow(1); % unit: s ; should be smaller than stimulation interval;
minrows = Interval*framerate;

%% Parameters related to matrix dimensions;
nROIs = 15;
nFrames  = 600;
nStimTrial = 25;    
nStimAll = 50;      % Max stim num of a single frequency among all trials;

%% Baseline Calculation: 1. moving smallest; 2. moving averge of lowest; 3. 25% smallest among all frames; 4. Constant value;
BaselineOption = 1;

%% Filter Option: 1.Select Sound reponsive curve to plot; 0. Plot all curve without filtering
FilterOn = 1;

threshold = 1.5;  % threshold to filter out the spontaneous firing;

%% Calculate Option: 1.Response value; Version 2: Local maximum; Version 3: Maximum during the responsive window.
MaxOption = 2;


%% Time Constants (default)
dt = 1/framerate;    %Frame interval/0.2s;
t0 = 2*dt;   %Time constant for noise filtering using exponentially weighted moving averaged.
t1 = 3*dt;   %Time window for smoothing;
t2 = 15;     %Time window for finding the baseline F0

BaseTWin = round(t2/dt)-1;

%% Fileanme and path
totaloutdir = strcat(maindir, '\ExtractOutput_Bv', num2str(BaselineOption),'_Fo',num2str(FilterOn),'_Mv', num2str(MaxOption));
mkdir(totaloutdir);

savepath0 = strcat(totaloutdir,'\Stims_TrialSorted.mat');
savepath1 = strcat(totaloutdir,'\dFoF_TrialSorted.mat');
savepath2 = strcat(totaloutdir,'\dFoF_StimSorted.mat');
savepath3 = strcat(totaloutdir,'\dFoF_ROISorted.mat');

savepath7 = strcat(totaloutdir,'\Amplitude_ROISorted.mat');
savepath8 = strcat(totaloutdir,'\PeakTime_ROISorted.mat');

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
dFoF_StimSorted = nan(minrows,nROIs,nStimAll,nTrials); 
dFoF_ROISorted = nan(minrows,nStimAll,nTrials, nROIs); 
Amplitude_ROISorted = nan(nTrials,nStimTrial,nROIs);   
PeakTime_ROISorted = nan(nTrials,nStimTrial,nROIs); 

tseries = transpose((AVGWindow(1):dt:AVGWindow(2)-dt));

%% Section 0: read files into matrix;
for i = 1 : nTrials
    TrialFile = strcat(maindir,'\Trial',num2str(i),'.xlsx');
    StimFile = strcat(maindir,'\Stim',num2str(i),'.xlsx');
    
    %open excel file;
    tmpFtRaw = readmatrix(TrialFile); % raw data; first column is time sequence; the others are F of differenct ROIs;
    %system('taskkill /F /IM EXCEL.EXE');
    inFtRaw(:,:, i) = tmpFtRaw(:,2:end);    
    
    tmpStimSeq = readmatrix(StimFile); % stim data; first row is time sequence of stimulation; the other is audio frequency;
    %system('taskkill /F /IM EXCEL.EXE');
    
    num = size(tmpStimSeq,1);
    inStimSeq(1:num,1,i) = round((tmpStimSeq(1:num,1) + AVGWindow(1)).*framerate + 1);
    inStimSeq(1:num,2,i) = tmpStimSeq(1:num,2);     
end


%% Section 1: Output files 'dFoF_TrialSorted','dFoF_StimSorted';
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
    
    %% Calculate fluorescence baseline F0;
    switch (BaselineOption)
        case 1 %% Version 1: moving smallest as baseline:
            %Baseline F0(t) is the minimum value of smooth Ft_Sm within a time window before t; time window equals 20dt= 4s
            F0_t = zeros(m,n);
            for j = 1:BaseTWin
                F0_t(:,j) = min(Ft_Sm(:,1:j),[],2);
            end
            for j = BaseTWin + 1:n
                F0_t(:,j) = min(Ft_Sm(:,j-BaseTWin:j),[],2);
            end
        case 2  %% Version 2: moving smallest as baseline:
            %Baseline F0(t) is the average among the 40% smallest value of the smooth Ft_Sm within a time window before t;
            F0_t = zeros(m,n);
            for j = 1:BaseTWin
                lwsort = sort(Ft_Sm(:,1:j),2);                                  % sorted by ascending;
                F0_t(:,j) = mean(lwsort(:,1:1+round(0.4*j)),2);
            end
            for j = BaseTWin + 1:n
                lwsort = sort(Ft_Sm(:,j-BaseTWin:j),2);
                F0_t(:,j) = mean(lwsort(:,1:1+round(0.4*BaseTWin)),2);          % the average of 40% smallest values;
            end
            
        case 3 %% Version 3: 25% smallest among all frames
            lwsort = sort(Ft_Sm(:,:),2);
            F0_t = zeros(m,n);
            
            for i = 1:n
                F0_t(:,i) = mean(lwsort(:,1:1+round(0.15*size(Ft_Sm,2))),2);      % the average of 25% smallest values among all frames
            end
            
        case 4 %% Version 4: Constant value
            F0_t = 4;
            
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
    
    for i = 2:tlength1
        % Extract data with different stimulation;
        dFoF_Seq = dFoF(:,inStimSeq(i,1,N):(inStimSeq(i,1,N) + minrows - 1));        
        % assign the dFoF_Seq to first j of Non-NaN matrix;               
        tmp = dFoF_StimSorted(1,1,:,N);
        tcol = length(tmp(~isnan(tmp)));
        dFoF_StimSorted(:,:,tcol + 1,N) = transpose(dFoF_Seq);           
    end
end


%% Section 2: Generate and output "dFoF_ROISorted" based on "dFoF_FreqSorted";   
for i = 1:nTrials
    for m = 1:nROIs        
        % find the number of nStims that have values, the following are
        % NaN;
        tmp1 = dFoF_StimSorted(1,m,:,i);
        tlength1 = length(tmp1(~isnan(tmp1)));
        
        if tlength1 > 0                     % non all rows are NaN for tStim;
            for n = 1:tlength1
                %% Selected the sound responsive curve for plotting
                if FilterOn == 1
                    %determine whether the neurons is tone responsive using
                    %peak detection; Criteria: the maxmum peak should be
                    %within 2.0s after the stimulation;
                    tmp2 = dFoF_StimSorted(:,m,n,i);
                    tmp2 = smooth(tmp2,5);
                    
                    [pks, locs] = findpeaks(tmp2,'SortStr','descend');
                    
                    % % peak time should be [0.5, 3s]; peak value should equal to maximum;
                    if ~isempty(locs) && locs(1) > (0.5-AVGWindow(1))*framerate && locs(1) < (3-AVGWindow(1))*framerate && pks(1)== max(tmp2) %&& pks(1) < threshold                    
                        tmp3 = dFoF_ROISorted(1,:,i,m);
                        tcol = length(tmp3(~isnan(tmp3)));
                        dFoF_ROISorted(:,tcol + 1,i,m) = dFoF_StimSorted(:,m,n,i);                        
                    end
                    
                    %% Plot all curve without filtering
                elseif FilterOn == 0
                    tmp3 = dFoF_ROISorted(1,:,i,m);
                    tcol = length(tmp3(~isnan(tmp3)));
                    dFoF_ROISorted(:,tcol + 1,i,m) = dFoF_StimSorted(:,m,n,i);
                end
            end
        end
    end
end



save(savepath0, 'inStimSeq','-v7.3','-nocompression');
save(savepath1, 'dFoF_TrialSorted','-v7.3','-nocompression');
save(savepath2, 'dFoF_StimSorted','-v7.3','-nocompression');
save(savepath3, 'dFoF_ROISorted','-v7.3','-nocompression');


%% Section 3: Save the data into excel files (Optional); Note: this will take a lot of time;  

% savepath4 = strcat(totaloutdir,'\dFoF_TrialSorted_Bv',num2str(BaselineOption),'.xlsx');
% savepath5 = strcat(totaloutdir,'\dFoF_StimSorted_Bv',num2str(BaselineOption),'.xlsx');
% savepath6 = strcat(totaloutdir,'\dFoF_ROISorted_Bv',num2str(BaselineOption),'_Fv',num2str(FilterOn),'.xlsx');
% 
% xlswrite(savepath4,[0 0]);
% system('taskkill /F /IM EXCEL.EXE');
% xlswrite(savepath5,[0 0]);
% system('taskkill /F /IM EXCEL.EXE');
% xlswrite(savepath6,[0 0]);
% system('taskkill /F /IM EXCEL.EXE');
% 
% %% Output to excel files;
% for i = 1: nTrials
%     strTrial = strcat('Trial',num2str(i,'%02d'));
%     xlswrite(savepath4,dFoF_TrialSorted(:,:,i),strTrial);
%     system('taskkill /F /IM EXCEL.EXE');
%     
%     tmp = dFoF_StimSorted(1,1,:,i);
%     tcol = length(tmp(~isnan(tmp)));
%     for j=1:tcol        
%         if j==1
%             xlswrite(savepath5,dFoF_StimSorted(:,:,1,i),strTrial,'B1');
%             system('taskkill /F /IM EXCEL.EXE');
%         else
%             b = xlsread(savepath5,strTrial);
%             system('taskkill /F /IM EXCEL.EXE');
%             nRows = size(b,1) + 2;
%             pos = strcat('B',num2str(nRows));
%             xlswrite(savepath5,dFoF_StimSorted(:,:,j,i),strTrial,pos);
%             system('taskkill /F /IM EXCEL.EXE');
%         end
%         
%     end
% end

%% Section 4: 
for i = 1:nROIs
    for j = 1: nTrials
        
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
                    % response value equals the peak value (0.5:3) s minus average value (-0.5:0) s;
                    Res_Value = mean(dFoF_ROISorted((0.5-AVGWindow(1))*framerate:(3-AVGWindow(1))*framerate,k,j,i), 1) - mean(dFoF_ROISorted(-0.5-AVGWindow(1)*framerate:-AVGWindow(1)*framerate,k,j,i), 1);
                    
                    % assign the dFoF_Seq to first j of Non-NaN matrix;
                    if Res_Value > 0
                        Amplitude_ROISorted(j,tcol1+1,i) = max(dFoF_ROISorted((0.5-AVGWindow(1))*framerate:(3-AVGWindow(1))*framerate,k,j,i),[],1);
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
                    elseif locs(1) > (0.5-AVGWindow(1))*framerate && locs(1) < (3-AVGWindow(1))*framerate && pks(1)== max(tempb)                             % peak time should be [0.5s, 3s]; peak value should equal to maximum;
                        Amplitude_ROISorted(j,tcol1+1,i) = max(dFoF_ROISorted((0.5-AVGWindow(1))*framerate:(3-AVGWindow(1))*framerate,k,j,i),[],1);          % ES responsive;
                    else
                        Amplitude_ROISorted(j,tcol1+1,i)= nan;                                                      % not tone responsive neuron, set as 0;
                    end
                case 3
                    %% Version 3: Maximum during the responsive window
                    % amplitude equals the peak value in (0.5:3)s 
                    Amplitude_ROISorted(j,tcol1+1,i) = max(dFoF_ROISorted((0.5-AVGWindow(1))*framerate:(3-AVGWindow(1))*framerate,k,j,i),[],1);
            end
            
            %% calculate peak time;
            tempb = transpose(dFoF_ROISorted(:,k,j,i));
            tempb = smooth(tempb,5);
            [pks, locs] = findpeaks(tempb,'SortStr','descend');
            if isempty(locs)
                PeakTime_ROISorted(j,tcol2+1,i) = nan;                                                     % not ES responsive, amplitude is 0;
            elseif locs(1) > (0.5-AVGWindow(1))*framerate && locs(1) < (3-AVGWindow(1))*framerate && pks(1)== max(tempb)                           % peak time should be [1s, 3s]; peak value should equal to maximum;
                PeakTime_ROISorted(j,tcol2+1,i) = tseries(locs(1));                                        % peak time;
            else
                PeakTime_ROISorted(j,tcol2+1,i) = nan;                                                     % not tone responsive neuron, set as 0;
            end
        end
    end
end


save(savepath7, 'Amplitude_ROISorted','-v7.3','-nocompression');
save(savepath8, 'PeakTime_ROISorted','-v7.3','-nocompression');

%% Output to excel files;
% savepath9 = strcat(totaloutdir,'\Amplitude_ROISorted_Mo',num2str(MaxOption),'.xlsx');
% savepath10 = strcat(totaloutdir,'\PeakTime_ROISorted_Mo',num2str(MaxOption),'.xlsx');
% 
% for i = 1: nROIs
%     strROI = strcat('ROI_',num2str(i));    
%     xlswrite(savepath9,Amplitude_ROISorted(:,:,i),strROI);
%     system('taskkill /F /IM EXCEL.EXE');
%     xlswrite(savepath10,PeakTime_ROISorted(:,:,i),strROI);
%     system('taskkill /F /IM EXCEL.EXE');  
% end