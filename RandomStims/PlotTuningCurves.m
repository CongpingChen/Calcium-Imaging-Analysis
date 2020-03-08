clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Plot Tuning curves of all ROIs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input: Directory containing "dFoF_ROISorted.mat" or "dFoF_ROISorted_Norm.mat" file.
%Output: "dFoF_TuningCurves.svg";
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Add path;
Main_dir = fileparts(which('PlotMultiCurves.m'));
addpath(fullfile(Main_dir, 'Plot_functions'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% User defined parameters;
maindir = uigetdir;

% Select ROI to plot;
ROI_Array = [1 2 3 4 5 6 7 8 9 10];
%ROI_Array = [11 12 13 14 15 16 17 18 19 20];
%ROI_Array = [25 26 27 28 29 30 31 32];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Default parameters of input files;
framerate = 5;  %Hz;
baseFreq = 3;
AVGWindow = [-0.6 4];  % unit: s;  
Interval = AVGWindow(2) - AVGWindow(1);        % unit: s ; smaller than stimulation interval;

% Parameters related to matrix dimensions;
nROIs = 30;
FreqNum = 25;   % 3K~48K totally 25 frequencies: freq = baseFreq*2^((i-1)/6); i = 1 + 6*log2(freq/baseFreq);
nFrames  = 600;

dt = 1/framerate;
minrows = Interval*framerate;
nStimTrial = 25;
nStimAll = 20;   % Max stim num of a single frequency among all trials; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plotting parameters;
xlmin = -1.2;
xlmax = AVGWindow(2);
ylmin = -0.2;
% ylmax = 5;     % Max dFoF should be normalized to 1;
% xyratio = 1.2;  
ylmax = 1.0;     % Max dFoF should be normalized to 1;
xyratio = 1.5;  

%tseries = transpose((0:1:minrows-1)*dt);
tseries = transpose((AVGWindow(1)+dt:dt:AVGWindow(2)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Start processing;
filepath = strcat(maindir,'\dFoF_ROISorted_Norm.mat');

for i = 1:length(ROI_Array)
    
        dFoF_ROISorted = importdata(filepath);
        
        Time_N = minrows;
        
        %xlmax = Interval;    
        
        for j = 1:FreqNum
            
            h = subaxis(length(ROI_Array)+1,FreqNum,(i-1)*FreqNum + j,'Spacing',0);           
            
            tmp = dFoF_ROISorted(1,:,j,ROI_Array(i));
            tlength = length(tmp(~isnan(tmp)));         %length of non-nan rows;
                        
            if tlength == 0               
                dFoF = zeros(minrows,1);
            else
                dFoF = zeros(minrows,tlength);
                dFoF(:,:) = dFoF_ROISorted(:,1:tlength,j,ROI_Array(i));
            end
            
                tdFoF = [tseries dFoF];                 % first column is time series, the others are dFoF of different stims;                
                plotcurves(tdFoF,minrows,xlmin,xlmax,ylmin,ylmax,xyratio);
                

                if j==1
                    hPos = get(h,'Position');
                    text(-5,ylmax/2,num2str(ROI_Array(i)),'FontName','AvantGarde','FontSize',20,'FontWeight','bold');
                end
        end
end

%% Plot dF/F scale bar
%plot([xlmax; xlmax], [ylmax/2-0.25; ylmax/2+0.25], '-k', 'LineWidth', 4);
plot([xlmax; xlmax], [ylmax/2-0.5; ylmax/2+0.5], '-k', 'LineWidth', 4);
%h = text(8.5,ylmax/2, 'dF/F', 'HorizontalAlignment','center','FontName','AvantGarde','FontSize',12);
%set(h,'Rotation',90);

%% Plot Stim annotation;

for k = 1:FreqNum    
    freq = baseFreq*2^((k-1)/6);
    subaxis(length(ROI_Array)+1,FreqNum,length(ROI_Array)*FreqNum + k,'Spacing',0);
    plotshaded([0 0.5],[freq freq],'k',xlmin,xlmax,ylmin,ylmax,xyratio);
    line([0 0],[0 48],'Color',[0.8 0.8 0.8],'LineWidth',2); % draw vertical line
    if k==1
        plot([xlmin; xlmin],[3; 48],'-k',[xlmin;xlmin+0.8],[3;3],'-k',[xlmin;xlmin+0.8],[48;48],'-k','LineWidth',4); %Freq annotation bar;
        %text(-10,-5,'3KHz');
        %text(-10,48,'48');
    end
end

%% Plot time scale bar;
plot([xlmax-5; xlmax], [-4; -4], '-k', 'LineWidth', 4);
%text(6,-18, '5s', 'HorizontalAlignment','right','FontName','AvantGarde','FontSize',12);

%% Save figure;
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
savepath = strcat(maindir,'\dFoF_TuningCurves_Fv1_Norm.svg');
saveas(gcf,savepath);
close;