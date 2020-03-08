clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Plot curves of different trials among all ROIs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input: Directory containing "Amplitude_ROISorted.mat".
%Output: "dFoF_TuningCurves.svg";
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Add path;
Main_dir = fileparts(which('PlotMultiCurves.m'));
addpath(fullfile(Main_dir, 'Plot_functions'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% User defined parameters;
maindir = uigetdir;

% Select ROI to plot;
ROI_Array = [1 2 3 4 5 6];
%ROI_Array = [11 12 13 14 15 16 17 18 19 20];
%ROI_Array = [21 22 23 24 25 26 27 28 18 19];
%ROI_Array = [5 8 9 11 12 13 14 15 16 20];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default parameters of input files;
% framerate = 5;  %Hz;
% baseFreq = 3;
% Interval = 4;   % unit: s; smaller than stimulation interval;
nCol = 3;

% Parameters related to matrix dimensions;
%nROIs = 10;
nTrials = 27; 
nBaseTrials = 7;
TrialInt = 3;   % unit: min; Trial interval;
%nFrames  = 600;

% dt = 1/framerate;
% minrows = Interval*framerate;
% nStimTrial = 25;
% nStimAll = 20;   % Max stim num of a single frequency among all trials; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plotting parameters;
ylmin = -2;
% ylmax = 5;     % Max dFoF should be normalized to 1;
% xyratio = 1.2;  
ylmax = 2;     % Max dFoF should be normalized to 1;
xyratio = 0.4;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Amplitude_ROISorted = nan(nTrials,nStimTrial,nROIs);   
% PeakTime_ROISorted = nan(nTrials,nStimTrial,nROIs); 


%% Start processing;
tseries = transpose([(-nBaseTrials:1:-1) (1:1:nTrials-nBaseTrials)]*TrialInt);

xlmin = tseries(1)-1;
xlmax = tseries(end)+1;

filepath = strcat(maindir,'\Amplitude_ROISorted.mat');
Amplitude_ROISorted = importdata(filepath);

nrows = length(tseries);
%Time_N = minrows;
%Amp_ROISelected = Amplitude_ROISorted(1:nTrials,:,ROI_Array);

for i = 1:length(ROI_Array)    
    
    h = subaxis(length(ROI_Array)/nCol,nCol,i,'SpacingVert',0.04,'SpacingHoriz',0.04);
    
    tmpAmp(:,:) = Amplitude_ROISorted(:,:,ROI_Array(i));
    %AVG_Amp = mean(tmpAmp,2,'omitnan');
        
    % Offset the baseline as 0 to better see the trend after LTP induction;   
    AVG_Amp = mean(tmpAmp,2,'omitnan') - mean(mean(tmpAmp(1:nBaseTrials,:),2,'omitnan'));
        
    % remove nan from AVG_Amp;
    nan_k = find(isnan(AVG_Amp));        
    for k = 1:length(nan_k)
        if nan_k(k) == 1
            AVG_Amp(1) = AVG_Amp(2);
        elseif nan_k(k) == length(AVG_Amp)
            AVG_Amp(end) = AVG_Amp(end-1);
        else
            AVG_Amp(nan_k(k)) =  (AVG_Amp(nan_k(k)-1) + AVG_Amp(nan_k(k)+1))/2;
        end
    end
        
    tmp0 = Amplitude_ROISorted(:,:,ROI_Array(i));
    tlength = sum(~isnan(tmp0),2);                    %length of non-nan rows;
    
    STD_Amp = std(tmpAmp,0,2,'omitnan');    
    SEM_Amp = STD_Amp./sqrt(tlength);
    
    STD_Amp(isnan(STD_Amp)) = 0;
    SEM_Amp(isnan(SEM_Amp)) = 0;
       
    tdata = [tseries AVG_Amp SEM_Amp];                              %       first column is time series, the others are dFoF of different stims;
    plotshaded(tdata,xlmin,xlmax,ylmin,ylmax,xyratio);
        
    hPos = get(h,'Position');
    text(xlmin-18,(ylmax+ylmin)/2,num2str(ROI_Array(i)),'FontName','AvantGarde','FontSize',20,'FontWeight','bold');

end

%% Plot dF/F scale bar
%plot([xlmax; xlmax], [ylmax/2-0.05; ylmax/2+0.05], '-k', 'LineWidth', 4);
%h = text(8.5,ylmax/2, 'dF/F', 'HorizontalAlignment','center','FontName','AvantGarde','FontSize',12);
%set(h,'Rotation',90);

%% Plot Stim annotation;

% for k = 1:nTrials    
%     freq = baseFreq*2^((k-1)/6);
%     subaxis(length(ROI_Array)+1,nTrials,length(ROI_Array)*nTrials + k,'Spacing',0);    
%     plot([0 0],[0 0.05],'Color',[0.6 0.6 0.6],'LineWidth',2);
%     hold on;
%     plot([xlmin-2 xlmax+2],[0 0],'Color',[0.6 0.6 0.6],'LineWidth',2);  % Plot baseline
%     xlim([xlmin xlmax])
%     ylim([ylmin ylmax/2])
%     pbaspect([1 xyratio xyratio])         %Display ratio for X and Y axis;
%     
%     axis off
%     set(gca, 'xtick',[])      %remove x tick;
%     set(gca, 'XTickLabel',[]) %remove x tick label;
%     set(gca, 'ytick',[])      %remove y tick;
%     set(gca, 'YTickLabel',[]) %remove y tick label;
%     %plotshaded([0 0],[0 ylmax/2],'k',xlmin,xlmax,ylmin,ylmax,xyratio);
%     %line([0 0],[0 48],'Color',[0.8 0.8 0.8],'LineWidth',2); % draw vertical line
% %     if k==1
% %         plot([xlmin; xlmin],[3; 48],'-k',[xlmin;xlmin+0.8],[3;3],'-k',[xlmin;xlmin+0.8],[48;48],'-k','LineWidth',4); %Freq annotation bar;
% %         %text(-10,-5,'3KHz');
% %         %text(-10,48,'48');
% %     end
% end

%% Plot time scale bar;
%plot([xlmax-2; xlmax], [-4; -4], '-k', 'LineWidth', 4);
%text(6,-18, '5s', 'HorizontalAlignment','right','FontName','AvantGarde','FontSize',12);

%% Save figure;
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
savepath = strcat(maindir,'\TrialAmp.svg');
saveas(gcf,savepath);
close;