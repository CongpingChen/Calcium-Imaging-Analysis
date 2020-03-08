clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Plot Receptive Field of all ROIs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input: Directory containing "ROISorted.xlsx" of different Sound Pressure Level (SPL) or treatment;
%Output: "dFoF_ReceptiveField.svg" of each neurons.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Add path;
Main_dir = fileparts(which('PlotMultiCurves.m'));
addpath(fullfile(Main_dir, 'Plot_functions'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% User defined parameters;
maindir = uigetdir;
ROI_Array = [2]; % Select ROI to plot;
%ROI_Array = [2 3 4 5 7 10 14 17 22]; % Select ROI to plot;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Default parameters of input files;
framerate = 5;  %Hz;
baseFreq = 3;
AVGWindow = [-0.6 4];  % unit: s;  
Interval = AVGWindow(2) - AVGWindow(1);        % unit: s ; smaller than stimulation interval;
dt = 1/framerate;  
%% Parameters related to matrix dimensions;
nROIs = 12;
FreqNum = 25;   % 3K~48K totally 25 frequencies: freq = baseFreq*2^((i-1)/6); i = 1 + 6*log2(freq/baseFreq);
nFrames  = 600;
minrows = Interval*framerate;
nStimTrial = 25;
nStimAll = 20;   % Max stim num of a single frequency among all trials; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plotting parameters;
xlmin = -1.2;
xlmax = AVGWindow(2);
ylmin = -0.2;
ylmax = 1;  % Max dFoF should be normalized to 1;
xyratio = 3.0;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subdirpath = fullfile( maindir, '*.mat' );
dat = dir( subdirpath );

%% Find max amplitude among different files for normalization;
MaxAmplitude = zeros(1,length(ROI_Array));
for i = 1:length(ROI_Array)
    
    for n = 1 : length(dat)
        
        datapath = fullfile( maindir, dat(n).name);
        
        dFoF_ROISorted = importdata(datapath);
        
        dFoF =  dFoF_ROISorted(:,:,:,ROI_Array(i));
        
        TempAmp = max(max(max(dFoF)));
        
        if TempAmp > MaxAmplitude(i)
            MaxAmplitude(i) = TempAmp;
        end
    end
end

%% Normalization and plot receptive field for individual ROIs
%tseries = transpose((1:1:minrows)*dt);
tseries = transpose((AVGWindow(1)+dt:dt:AVGWindow(2)));

for i = 1:length(ROI_Array)
    
    for n = 1 : length(dat)
        
        datapath = fullfile( maindir, dat(n).name);
        dFoF_ROISorted = importdata(datapath);
        
        str_dB = dat(n).name;
        str_dB = str_dB(1:end-4);
                
        Time_N = minrows;       
        %xlmax = Interval;
          
        for j = 1:FreqNum
            
            h = subaxis(length(dat)+1,FreqNum,(n-1)*FreqNum + j,'Spacing',0);
            
            tmp = dFoF_ROISorted(1,:,j,ROI_Array(i));
            tlength = length(tmp(~isnan(tmp)));         %length of non-nan rows;
                        
            if tlength == 0               
                dFoF = zeros(minrows,1);
            else
                dFoF = zeros(minrows,tlength);
                dFoF(:,:) = dFoF_ROISorted(:,1:tlength,j,ROI_Array(i))./MaxAmplitude(i);% Normalization;
            end
            
                Norm_tdFoF = [tseries dFoF];                 % first column is time series, the others are dFoF of different stims; 
             
            plotcurves(Norm_tdFoF,minrows,xlmin,xlmax,ylmin,ylmax,xyratio);
            line([0 0],[ylmin-0.5 ylmax+0.5],'Color',[0.8 0.8 0.8],'LineWidth',2); % draw vertical line           
            if j==1
                hPos = get(h,'Position');
                text(-10,ylmax/2,str_dB,'FontName','AvantGarde','FontSize',20,'FontWeight','bold');
            end
        end
    end
    
    %% Plot dF/F scale bar
    plot([xlmax; xlmax], [ylmax/2-0.25; ylmax/2+0.25], '-k', 'LineWidth', 4);
    %h = text(8.5,ylmax/2, '0.5dF/F', 'HorizontalAlignment','center','FontName','AvantGarde','FontSize',18);
    %set(h,'Rotation',90);
    
    %% Plot Stim annotation;
    for k = 1:FreqNum
        freq = baseFreq*2^((k-1)/6);
        subaxis(length(dat)+1,FreqNum,length(dat)*FreqNum + k,'Spacing',0);
        plotshaded([0 0.5],[freq freq],'k',xlmin,xlmax,ylmin,ylmax,xyratio);
        line([0 0],[0 48],'Color',[0.8 0.8 0.8],'LineWidth',2); % draw vertical line
        if k==1
            plot([xlmin; xlmin],[3; 48],'-k',[xlmin;xlmin+0.8],[3;3],'-k',[xlmin;xlmin+0.8],[48;48],'-k','LineWidth',4); %Freq annotation bar;
%             text(-10,-5,'3KHz','FontSize',20);
%             text(-10,48,'48','FontSize',20);
        end
    end
    
    %% Plot time scale bar;
    plot([xlmax-2; xlmax], [-4; -4], '-k', 'LineWidth', 4);
    %text(4,-18, '5s', 'HorizontalAlignment','right','FontName','AvantGarde','FontSize',18);
       
    %% Save figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    savepath = strcat(maindir,'\ReceptiveField_ROI_',num2str(ROI_Array(i)),'.svg');
    saveas(gcf,savepath);
    close;
    
end