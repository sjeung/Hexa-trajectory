function VN_animate_baseline(travelPoints)
% plot the trajectories
%
% calls function
%       buildgridmap
%
% Inputs
%       travelPoints
%                   : 2X(number of points), x,y coordinates of positions to
%                   be visited
%       trialSpecification
%                   : trialSpecification(1,:) stores angular offsets
%                   : trialSpecification(2,:) stores trial lengths
%       nRun
%                   : number of the corresponding run
%       targetIndex
%                   : 1XnTrial, indices of targets 
%
% Outputs
%       none
%
% Last edit         : 2018.12.11
% Author            : Sein Jeung
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%                           Parameters
%--------------------------------------------------------------------------

% trial number
nTrial                  = 1; 

% arena ceter coordinates
xCenter                 = 0;
yCenter                 = 0;

% length of the sqaure arena axes in meters
xArenaAxis              = 13;
yArenaAxis              = 8;

% distal landmark coordinates
xLandmark               = [];
yLandmark               = [];

%--------------------------------------------------------------------------
%                              Animation
%--------------------------------------------------------------------------

% iterate over number of trials in one run
for n = 1:nTrial
    
    figure
    hold on
    
    % set figure limits
    xlim([-xArenaAxis, xArenaAxis])
    ylim([-yArenaAxis, yArenaAxis])
    
    % place global landmarks
    scatter(xLandmark,yLandmark, 100,'filled','MarkerFaceColor',[.7 .7 .4])
    
    % move along the trajectory
    %----------------------------------------------------------------------
    h = animatedline('Color',[.3 .3 .8],'LineWidth',2);
    
    if n == 1
        startingSegmentIndex  = 1;
    else
        startingSegmentIndex  = 37;
    end
    
    endingSegmentIndex  = 36;
    
    % mark the starting point of each trial 
    scatter(travelPoints(1,startingSegmentIndex),travelPoints(2,startingSegmentIndex),50,[.6 .9 .6],'filled')
    
    for p = startingSegmentIndex : endingSegmentIndex
        a       = [travelPoints(1,p),travelPoints(2,p)];
        b       = [travelPoints(3,p),travelPoints(4,p)];
        
        % straight line function from a to b
        func = @(x)a(2) + (a(2)-b(2))/(a(1)-b(1))*(x-a(1));
        
        % determine the x values
        x = linspace(a(1),b(1),50);
        
        % determine the y values
        y = func(x);
        
        % y value when func is not defined
        if  a(1) == b(1)
            y = linspace(a(2),b(2),50);
        end
        
        % get a handle to a plot graphics object
        hPlot = plot(NaN,NaN,'o','Color','k');
        
 
        for f = 1:length(x)
            
            addpoints(h,x(f),y(f));
            
            % update the plot graphics object with the next position
            set(hPlot,'XData',x(f),'YData',y(f));
            drawnow
        end
        
        % label each point visited
        text(x(length(x)),y(length(x)),num2str(p),'FontSize', 15)
        
    end
    
    % mark the ending point of each trial 
    scatter(travelPoints(3,endingSegmentIndex),travelPoints(4,endingSegmentIndex),50,[.9 .6 .6],'filled')
    
end

end