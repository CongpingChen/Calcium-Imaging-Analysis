clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Plot Tuning curves of all ROIs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input: Directory containing "RespSeg_ROISorted.mat" file.
%Output: "dFoF_SegAVG.svg";
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Add path;
Main_dir = fileparts(which('PlotMultiCurves.m'));
addpath(fullfile(Main_dir, 'Plot_functions'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% User defined parameters;
maindir = uigetdir;

% Select ROI to plot;
ROI_Array = [1 2 3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Default parameters of input files;
framerate = 5;  %Hz;
baseFreq = 3;
AVGWindow = [0 5];  % unit: s;  
Interval = AVGWindow(2) - AVGWindow(1);        % unit: s 

% Parameters related to matrix dimensions;
nROIs = 3;
FreqNum = 25;   % 3K~48K totally 25 frequencies: freq = baseFreq*2^((i-1)/6); i = 1 + 6*log2(freq/baseFreq);
nFrames  = 600;

dt = 1/framerate;
minrows = Interval*framerate;
nStimTrial = 25;
nStimAll = 20;   % Max stim num of a single frequency among all trials; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plotting parameters;
xlmin = AVGWindow(1);
xlmax = AVGWindow(2);
ylmin = -0.2;
% ylmax = 5;     % Max dFoF should be normalized to 1;
% xyratio = 1.2;  
ylmax = 10;     % Max dFoF should be normalized to 1;
xyratio = 2.0;  

%tseries = transpose((0:1:minrows-1)*dt);
tseries = transpose((AVGWindow(1)+dt:dt:AVGWindow(2)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Start processing;
filepath = strcat(maindir,'\RespSeg_ROISorted.mat');

for i = 1:length(ROI_Array)
    
    RespSeg_ROISorted = importdata(filepath);
    
    Time_N = minrows;
    
    %xlmax = Interval;
      
    h = subaxis(length(ROI_Array),1,i,'Spacing',0);
    
    tmp = RespSeg_ROISorted(1,:,ROI_Array(i));
    tlength = length(tmp(~isnan(tmp)));         %length of non-nan rows;
    
    if tlength == 0
        dFoF = zeros(minrows,1);
    else
        dFoF = zeros(minrows,tlength);
        dFoF(:,:) = RespSeg_ROISorted(1:minrows,1:tlength,ROI_Array(i));
    end
    
    tdFoF = [tseries dFoF];                 % first column is time series, the others are dFoF of different stims;
    
    if i==1
        plotcurves(tdFoF,minrows,'k',xlmin,xlmax,ylmin,ylmax,xyratio);
    elseif i==2
        plotcurves(tdFoF,minrows,'b',xlmin,xlmax,ylmin,ylmax,xyratio);
    else
        plotcurves(tdFoF,minrows,'r',xlmin,xlmax,ylmin,ylmax,xyratio);
    end
       
    hPos = get(h,'Position');
    text(-5,ylmax/2,num2str(ROI_Array(i)),'FontName','AvantGarde','FontSize',20,'FontWeight','bold');
    
    
end

%% Plot dF/F scale bar
%plot([xlmax; xlmax], [ylmax/2-0.25; ylmax/2+0.25], '-k', 'LineWidth', 4);
plot([xlmax-0.5; xlmax-0.5], [3; 4], '-k', 'LineWidth', 4);
%h = text(8.5,ylmax/2, 'dF/F', 'HorizontalAlignment','center','FontName','AvantGarde','FontSize',12);
%set(h,'Rotation',90);

%% Plot time scale bar;
plot([xlmax-1.5; xlmax-0.4], [3; 3], '-k', 'LineWidth', 4);
%text(6,-18, '5s', 'HorizontalAlignment','right','FontName','AvantGarde','FontSize',12);

%% Save figure;
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
savepath = strcat(maindir,'\dFoF_SegAVG.svg');
saveas(gcf,savepath);
close;