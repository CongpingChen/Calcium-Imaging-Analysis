clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Plot dFoF curves of all trials for all ROIs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input: Directory containing "AllTrials.xlsx" and "AllStims.xlsx" files.
%Output: "dFoF_MultiCurves.svg";
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Add path
Main_dir = fileparts(which('PlotMultiCurves.m'));
addpath(fullfile(Main_dir, 'Plot_functions'));

maindir = uigetdir;
Trial_N = 3;                          
ROI_Array = [1 2 3];


% Default parameters of input files;
framerate = 5;  %Hz;
dt = 1/framerate;
baseFreq = 3;
Interval = 4;   % unit: s; smaller than stimulation interval;

% Parameters related to matrix dimensions;
nROIs = 12;
FreqNum = 25;   % 3K~48K totally 25 frequencies: freq = baseFreq*2^((i-1)/6); i = 1 + 6*log2(freq/baseFreq);
PlotFrames  = [100,600];
minrows = Interval*framerate;
nStimTrial = 25;
nStimAll = 20;   % Max stim num of a single frequency among all trials; 

%% plotting parameters;
ROI_N = length(ROI_Array);
xlmin = PlotFrames(1)/framerate-0.5;
xlmax = PlotFrames(2)/framerate;
ylmin = -0.2;
ylmax = 8;
xyratio = 0.5;


%filepath1 = strcat(maindir,'\Stims_TrialSorted.mat');
filepath = strcat(maindir,'\dFoF_TrialSorted.mat');

%Stims_All = importdata(filepath1);
dFoF_All  = importdata(filepath);

tseries = transpose((PlotFrames(1):1:PlotFrames(2))*dt);
%% Plot Trial X ROI subplot;

for j=1:Trial_N
    for i = 1:ROI_N
             
      
        dFoF(:,1) = dFoF_All(PlotFrames(1):PlotFrames(2),ROI_Array(i),j);
        
        datax = tseries;
        datay = [dFoF zeros(length(dFoF),1)];
        
        %h = subaxis(ROI_N,Trial_N,(i-1)*Trial_N + j,'Spacing',0);
        %h = subaxis(ROI_N,1,i,'Spacing',0);

        if (i==1)
            plotshaded(datax,datay,'k',xlmin,xlmax,ylmin,ylmax,xyratio); %5/60
        elseif (i==2)
            plotshaded(datax,datay,'b',xlmin,xlmax,ylmin,ylmax,xyratio);
        else
            plotshaded(datax,datay,'r',xlmin,xlmax,ylmin,ylmax,xyratio);
        end
        %hPos = get(h,'Position');
        %text(-10,1,num2str(ROI_Array(i)),'FontName','AvantGarde','FontSize',20,'FontWeight','bold');

    end
    
    %% Plot dF/F scale bar
    plot([xlmax-1; xlmax-1], [4; 5], '-k', 'LineWidth', 4);
    %h = text(xlmax+2,ylmax/2, 'dF/F', 'HorizontalAlignment','center','FontName','AvantGarde','FontSize',20);
    %set(h,'Rotation',90);
    
    
    %% Plot time scale bar
    hold on
    plot([xlmax-6; xlmax-0.8], [4; 4], '-k', 'LineWidth', 4);
    %text(xlmax,-3, '5s', 'HorizontalAlignment','right','FontName','AvantGarde','FontSize',15);
    
    %% Save figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    savepath = strcat(maindir,'\dFoF_MultiCurves_T',num2str(j),'.svg');
    saveas(gcf,savepath);
    close;
    
end


