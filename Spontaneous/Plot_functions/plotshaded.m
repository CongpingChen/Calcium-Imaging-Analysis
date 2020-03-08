function varargout = plotshaded(x,y,fstr,xlmin,xlmax,ylmin,ylmax,xyratio)
% x: x coordinates
% y: either just one y vector, or 2xN  y-data
% for 2xN matrix of y-data, y(1,:) should be the mean while y(2,:) should be the boundary, e.g. mean+/-std;
% fstr: format ('r' or 'b--' etc)
%
% example
% x=[-10:.1:10];plotshaded(x,[sin(x.*1.1)+1;sin(x*.9)-1],'r');

if size(y,1)>size(y,2)
    y=y';
end

if size(y,1)==1 % just plot one line
    
    plot(x,y,fstr,'LineWidth',3);
    hold on;
    plot([xlmin xlmax],[0 0],'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',2);  % Plot baseline
    xlim([xlmin xlmax])
    ylim([-4 48])
    pbaspect([1 xyratio xyratio])         %Display ratio for X and Y axis;
    
    %     plot(x,y(1,:),fstr,'LineWidth',2);
    %     xlim([-0.2 5.2])
    %     ylim([-0.2 3])
    %     pbaspect([1 1.8 1.8])         %Display ratio for X and Y axis;
    
    %     axis off
    %     set(gca, 'xtick',[])      %remove x tick;
    %     set(gca, 'XTickLabel',[]) %remove x tick label;
    %     set(gca, 'ytick',[])      %remove y tick;
    %     set(gca, 'YTickLabel',[]) %remove y tick label;
    
end

if size(y,1)==2 %plot shaded area
    
    % if y(:,2) = 0; plot just one line;
    if y(2,:) == 0
        plot([xlmin xlmax],[0 0],'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',2);  % Plot baseline
        hold on;
        plot(x,y(1,:),fstr,'LineWidth',2);
        xlim([xlmin xlmax])
        ylim([ylmin ylmax])
        pbaspect([1 xyratio xyratio])         %Display ratio for X and Y axis;
        
        % else plot mean and std patch
    else
        yt =[ y(1,:) - y(2,:);y(1,:) + y(2,:)];
        px=[x,fliplr(x)]; % make closed patch
        py=[yt(1,:), fliplr(yt(2,:))];
        patch(px,py,1,'FaceColor',[0.2 0.2 0.2],'EdgeColor','none');
        hold on;
        plot([xlmin xlmax],[0 0],'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',2);  % Plot baseline
        hold on;
        plot(x,y(1,:),fstr,'LineWidth',2);
        
        xlim([xlmin xlmax])
        ylim([ylmin ylmax])
        pbaspect([1 xyratio xyratio])         %Display ratio for X and Y axis;
    end
    
    
end


alpha(.2); % make patch transparent

axis off
set(gca, 'xtick',[])      %remove x tick;
set(gca, 'XTickLabel',[]) %remove x tick label;
set(gca, 'ytick',[])      %remove y tick;
set(gca, 'YTickLabel',[]) %remove y tick label;