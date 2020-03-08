function varargout = plotcurves(data,minrows,xlmin,xlmax,ylmin,ylmax,xyratio)
% data: minrows*(1+nStims) matrix; first dimension is time, others are dFoF of different stims;
%
% example
% x=[-10:.1:10];plotshaded(x,[sin(x.*1.1)+1;sin(x*.9)-1],'r');

n_Stims = size(data,2)-1;

plot([xlmin xlmax],[0 0],'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',2);  % Plot baseline
hold on;

plot([0 0],[ylmin-0.5 ylmax+0.5],'Color',[0.8 0.8 0.8],'LineWidth',2); % draw vertical line
hold on;

% if n_Stims == 0
%     x = transpose(data(1:minrows,1));
%     y = transpose(zeros(minrows,1));
%     plot(x,y,'k','LineWidth',4);
%     xlim([xlmin xlmax])
%     ylim([ylmin ylmax])
%     pbaspect([1 xyratio xyratio])         %Display ratio for X and Y axis;
% else
    for i = 1: n_Stims        
        % plot curves of different trials;
        x = transpose(data(1:minrows,1));
        y = transpose(data(1:minrows,1+i));

        plot(x,y,'Color',[0.6 0.6 0.6],'LineWidth',2);
        hold on
        xlim([xlmin xlmax])
        ylim([ylmin ylmax])
        pbaspect([1 xyratio xyratio])         %Display ratio for X and Y axis;
    end
    
    % plot average curves of different trials;
    x = transpose(data(1:minrows,1));
    y = transpose(mean(data(1:minrows,2:1+n_Stims),2));
    plot(x,y,'k','LineWidth',4);
    %hold on;
    xlim([xlmin xlmax])
    ylim([ylmin ylmax])
    pbaspect([1 xyratio xyratio])         %Display ratio for X and Y axis;
    

%     plot(x,y,fstr,'LineWidth',2);
%     hold on;
%     plot([xlmin xlmax],[0 0],'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',2);  % Plot baseline
%     xlim([xlmin xlmax])
%     ylim([ylmin ylmax])
%     pbaspect([1 xyratio xyratio])         %Display ratio for X and Y axis;

axis off
set(gca, 'xtick',[])      %remove x tick;
set(gca, 'XTickLabel',[]) %remove x tick label;
set(gca, 'ytick',[])      %remove y tick;
set(gca, 'YTickLabel',[]) %remove y tick label;