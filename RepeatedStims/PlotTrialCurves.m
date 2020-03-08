clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Plot amplitude of different trials among all ROIs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input: Directory containing "dFoF_ROISorted.mat" or "dFoF_ROISorted_Norm.mat" file.
%Output: "dFoF_TuningCurves.svg";
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Add path;
Main_dir = fileparts(which('PlotTrialCurves.m'));
addpath(fullfile(Main_dir, 'Plot_functions'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% User defined parameters;
maindir = uigetdir;

% Select ROI to plot;
ROI_Array = [1 2 3 4 5 6 7];
%ROI_Array = [11 12 13 14 15 16 17 18 19 20];
%ROI_Array = [16 17 18 19 20 21 22 23 24 25];
%ROI_Array = [25 26 27 28 29 30 31 32];
%ROI_Array = [5 8 9 11 12 13 14 15 16 20];
%ROI_Array = [3 4 6 7 10 17];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Default parameters of input files;
framerate = 5;  %Hz;
AVGWindow = [-1 4];  % unit: s;  
Interval = AVGWindow(2) - AVGWindow(1);        % unit: s ; smaller than stimulation interval;
minrows = Interval*framerate;

% Parameters related to matrix dimensions;
%nROIs = 10;
nTrials = 2;   
nFrames  = 600;

dt = 1/framerate;
nStimTrial = 25;
nStimAll = 20;   % Max stim num of a single frequency among all trials; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plotting parameters;
xlmin = AVGWindow(1)-1.2;
ylmin = -0.2;
% ylmax = 5;     % Max dFoF should be normalized to 1;
% xyratio = 1.2;  
ylmax = 2;     % Max dFoF should be normalized to 1;
%xyratio = 1.3;   % 1.3 for 10 ROI * 27 Trials(TBS experiment);
xyratio = 1.4;    % 0.6 for 10 ROI* 11 Trials (ES test);

%tseries = transpose((0:1:minrows-1)*dt);
tseries = transpose((AVGWindow(1)+dt:dt:AVGWindow(2)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%dFoF_ROISorted = nan(minrows,nStimAll,nTrials, nROIs); 


%% Start processing;
filepath = strcat(maindir,'\dFoF_ROISorted.mat');
dFoF_ROISorted = importdata(filepath);

xlmax = AVGWindow(2);

for i = 1:length(ROI_Array)
    
        for j = 1:nTrials
            
            h = subaxis(length(ROI_Array)+1,nTrials,(i-1)*nTrials + j,'Spacing',0);
            
            tmp = dFoF_ROISorted(1,:,j,ROI_Array(i));
            
            tlength = length(tmp(~isnan(tmp)));          % length of non-nan rows;
            
            if tlength == 0
                dFoF = zeros(minrows,1);
            else
                dFoF = zeros(minrows,tlength);
                dFoF(:,:) = dFoF_ROISorted(1:minrows,1:tlength,j,ROI_Array(i));
            end
            
            tdFoF = [tseries dFoF];                      % first column is time series, the others are dFoF of different stims;                        
            plotcurves(tdFoF,minrows,xlmin,xlmax,ylmin,ylmax,xyratio);
            
            if j==1
                hPos = get(h,'Position');
                text(-8,ylmax/2,num2str(ROI_Array(i)),'FontName','AvantGarde','FontSize',28,'FontWeight','bold');
            end
        end
end

%% Plot dF/F scale bar
%plot([xlmax; xlmax], [ylmax/2-0.25; ylmax/2+0.25], '-k', 'LineWidth', 4);
plot([xlmax; xlmax], [ylmax/2-0.05; ylmax/2+0.05], '-k', 'LineWidth', 4);
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
%hold on
plot([xlmax-2; xlmax], [ylmin; ylmin], '-k', 'LineWidth', 4);
%text(6,-18, '5s', 'HorizontalAlignment','right','FontName','AvantGarde','FontSize',12);

%% Save figure;
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
savepath = strcat(maindir,'\dFoF_TrialCurves_Fv1.svg');
saveas(gcf,savepath);
close;