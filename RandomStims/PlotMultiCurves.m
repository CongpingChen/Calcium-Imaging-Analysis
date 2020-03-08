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
Trial_N = 1;                          
ROI_Array = [1 2 3];


% Default parameters of input files;
framerate = 5;  %Hz;
dt = 1/framerate;
baseFreq = 3;
Interval = 4;   % unit: s; smaller than stimulation interval;

% Parameters related to matrix dimensions;
nROIs = 12;
FreqNum = 25;   % 3K~48K totally 25 frequencies: freq = baseFreq*2^((i-1)/6); i = 1 + 6*log2(freq/baseFreq);
nFrames  = 600;
minrows = Interval*framerate;
nStimTrial = 25;
nStimAll = 20;   % Max stim num of a single frequency among all trials; 

%% plotting parameters;
ROI_N = length(ROI_Array);
xlmin = -1;
xlmax = 120;
ylmin = -0.2;
ylmax = 5.5;
xyratio = 0.3;


filepath1 = strcat(maindir,'\Stims_TrialSorted.mat');
filepath2 = strcat(maindir,'\dFoF_TrialSorted.mat');

Stims_All = importdata(filepath1);
dFoF_All  = importdata(filepath2);

tseries = transpose((1:1:nFrames)*dt);
%% Plot Trial X ROI subplot;
for i = 1:ROI_N
    
    for j = 1:Trial_N
                
            stim(:,:) = Stims_All(:,:,j);
            dFoF(:,1) = dFoF_All(:,ROI_Array(i),j);
            
            datax = tseries;            
            datay = [dFoF zeros(length(dFoF),1)];
           
            stimx = transpose(stim(:,1));      %time of stims
            %stimy = transpose(stim(:,2));      %freq of stims         
            %h = subaxis(ROI_N+1,Trial_N,(i-1)*Trial_N + j,'Spacing',0);
            h = subaxis(ROI_N,Trial_N,(i-1)*Trial_N + j,'Spacing',0);
%             for k = 1:length(stimx)
%                 line([stimx(k) stimx(k)],[ylmin-0.5 ylmax+0.5],'Color',[0.8 0.8 0.8],'LineWidth',2); % draw Stim Start Line
%             end
%           hold on;
            if (i==1)
                plotshaded(datax,datay,'k',xlmin,xlmax,ylmin,ylmax,xyratio); %5/60
            elseif (i==2)
                plotshaded(datax,datay,'b',xlmin,xlmax,ylmin,ylmax,xyratio);
            else
                plotshaded(datax,datay,'r',xlmin,xlmax,ylmin,ylmax,xyratio);
            end
            %line([0 0],[ylmin-0.5 ylmax+0.5],'Color',[0.8 0.8 0.8],'LineWidth',2); % draw Trial Start line
            if j==1
                hPos = get(h,'Position');
                text(-10,1,num2str(ROI_Array(i)),'FontName','AvantGarde','FontSize',20,'FontWeight','bold');
            end

    end
end

%% Plot dF/F scale bar
plot([xlmax-1; xlmax-1], [1; 2], '-k', 'LineWidth', 4);
%h = text(xlmax+2,ylmax/2, 'dF/F', 'HorizontalAlignment','center','FontName','AvantGarde','FontSize',20);
%set(h,'Rotation',90);

%% Plot Stim annotation;
% for l = 1:Trial_N
%     
%     stim(:,:) = Stims_All(:,:,l);
%     stimx = transpose(stim(:,1));      %time of stims
%     stimy = transpose(stim(:,2));      %freq of stims
%     
%     subaxis(ROI_N+1,Trial_N,Trial_N*ROI_N + l,'Spacing',0);
%     
%     for m = 1:length(stimx)
%         plotshaded([stimx(m) stimx(m)+0.5],[stimy(m) stimy(m)],'k',xlmin,xlmax,ylmin,ylmax,xyratio);
%         line([stimx(m) stimx(m)],[0 48],'Color',[0.8 0.8 0.8],'LineWidth',2); % draw vertical line
%     end
%     
%     if l==1
%         plot([xlmin; xlmin],[3; 48],'-k',[xlmin;xlmin+4],[3;3],'-k',[xlmin;xlmin+4],[48;48],'-k','LineWidth',3); %Freq annotation bar;
%         %text(xlmin-10,-5,'3KHz');
%         %text(xlmin-10,48,'48');
%     end
% end

%% Plot time scale bar
hold on
plot([xlmax-6; xlmax-0.63], [1; 1], '-k', 'LineWidth', 4);
%text(xlmax,-3, '5s', 'HorizontalAlignment','right','FontName','AvantGarde','FontSize',15);

%% Save figure;
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
savepath = strcat(maindir,'\dFoF_MultiCurves.svg');
saveas(gcf,savepath);
close;


