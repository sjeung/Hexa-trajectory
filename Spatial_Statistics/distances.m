function  distances(admissiblePaths, X, Y, xyBounds)


% four corners 
cornersX    = [xyBounds(1,1) xyBounds(1,1) xyBounds(1,2) xyBounds(1,2)];
cornersY    = [xyBounds(2,2) xyBounds(2,1) xyBounds(2,1) xyBounds(2,2)];

% initialize a vector to store mean distances
meanCornerDs        = NaN(1,36); 
meanBoundDs         = NaN(1,36); 
meanXBoundDs        = NaN(1,36); 
meanYBoundDs        = NaN(1,36); 

for i = 1:36
    
    % linear indices for segments associated with each heading
    headingIndices = admissiblePaths(:) == i;
    
    % x, y positions associated with each heading
    targetsX        = X(headingIndices);
    targetsY        = Y(headingIndices);
    
    % distances to the closest corner
    cornersD        = min(sqrt((repmat(targetsX,1,4)-cornersX).^2 + (repmat(targetsY,1,4)-cornersY).^2),[],2);

    % distances to the closest bound
    boundsD         = min([abs(repmat(targetsX,1,2)-xyBounds(1,:)) abs(repmat(targetsY,1,2)-xyBounds(2,:))],[],2);
    
    % distances to the closer x bound 
    xBoundsD        = min(abs(repmat(targetsX,1,2)-xyBounds(1,:)),[],2);
    
    % distances to the closer y bound
    yBoundsD        = min(abs(repmat(targetsY,1,2)-xyBounds(2,:)),[],2);
    
    
    meanCornerDs(i)     = mean(cornersD);
    meanBoundDs(i)      = mean(boundsD);
    meanXBoundDs(i)     = mean(xBoundsD);
    meanYBoundDs(i)     = mean(yBoundsD);
    
end

figure
%--------------------------------------------------------------------------
subplot(2,2,1)
polarhistogram('BinEdges',deg2rad(0:10:360),'BinCounts',meanCornerDs,'EdgeColor','none')
% hold on 
% 
% [offset,beta1,beta2,phi,p] = periodicGLM(meanCornerDs, 1);
% 
% theta = 0:0.01:2*pi;
% rho = beta1*cos(4*theta) + beta2*sin(4*theta) + offset;
% polarplot(theta,rho,'LineWidth',1)
title(['Distance to a corner'])

%--------------------------------------------------------------------------
subplot(2,2,2)
polarhistogram('BinEdges',deg2rad(0:10:360),'BinCounts',meanBoundDs,'EdgeColor','none')
% hold on 
% 
% [offset,beta1,beta2,phi,p] = periodicGLM(meanBoundDs, 1);
% 
% theta = 0:0.01:2*pi;
% rho = beta1*cos(4*theta) + beta2*sin(4*theta) + offset;
% polarplot(theta,rho,'LineWidth',1)
title(['Distance to a bound'])

%--------------------------------------------------------------------------
subplot(2,2,3)
polarhistogram('BinEdges',deg2rad(0:10:360),'BinCounts',meanXBoundDs,'EdgeColor','none')
% hold on 
% [offset,beta1,beta2,phi,p] = periodicGLM(meanXBoundDs, 2);
% 
% theta = 0:0.01:2*pi;
% rho = beta1*cos(4*theta) + beta2*sin(4*theta) + offset;
% polarplot(theta,rho,'LineWidth',1)
title(['Distance to a vertical bound'])

%--------------------------------------------------------------------------
subplot(2,2,4)
polarhistogram('BinEdges',deg2rad(0:10:360),'BinCounts',meanYBoundDs,'EdgeColor','none')
% hold on 
% [offset,beta1,beta2,phi,p] = periodicGLM(meanYBoundDs, 4);
% 
% theta = 0:0.01:2*pi;
% rho = beta1*cos(4*theta) + beta2*sin(4*theta) + offset;
% polarplot(theta,rho,'LineWidth',1)
title(['Distance to a horizontal bound'])

end
