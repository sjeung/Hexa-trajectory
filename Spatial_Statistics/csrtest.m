function [R, pval] = csrtest(points, environmentBounds)

% Assessement of spatial randomenss using nearest neighbor statistics
%--------------------------------------------------------------------------
%
% Inputs
%                ID : integer, participant ID 
%           session : integer, session index
%
% Outputs
%              none
%--------------------------------------------------------------------------
%points              = rand(10,2); 
%points            = [[.1 .1 .1 .2 .2 .2 .3 .3 .3]' [.1 .2 .3 .1 .2 .3 .1 .2 .3]'];
%environmentBounds  = [0 1; 0 1]; 

% area
A           = (environmentBounds(1,2)-environmentBounds(1,1))*(environmentBounds(2,2)-environmentBounds(2,1));

% number of points in the area
N           = size(points,1);


% expected nearest neighbor distance under randomness
%--------------------------------------------------------------------------
nullR           = .5*1/sqrt(N/A);

scatter(points(:,1),points(:,2))

% % create a distance matrix : due to memory limit, doesn't work 
% %--------------------------------------------------------------------------
% distanceMatrix  = NaN(N,N); 
% 
% for r = 1:N
%     
%     for c = 1:N 
%         
%         if r == c
%         
%             % the distance to self is given as infinity
%             distanceMatrix(r,c)     = inf; 
%             
%         else
%             
%             distanceMatrix(r,c)     = sqrt((points(r,1)-points(c,1))^2 + (points(r,2)-points(c,2))^2);
%         
%         end
%     end
% end


% empirical nearest neighbor distance 
%--------------------------------------------------------------------------

% initialize a vector of nearst neighbor distances for all points 
nnDistances         = NaN(1,N); 

% iterate over some points
for p = 1:N

    distances       = sqrt((points(p,1)-points(:,1)).^2 + (points(p,2)-points(:,2)).^2);
    distances(p)    = inf; 
    
    % distance to all other points 
    R               = min(distances(:)); 
     
    nnDistances(p)  = R; 

end 

empiricalR = mean(nnDistances); 
R = empiricalR/nullR;

variance = (4-pi)/(4*pi*N*N/A);%std(nnDistances);
z = (empiricalR-nullR)/sqrt(variance); 
pval = normcdf(z); 

end 