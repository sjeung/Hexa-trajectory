function [stopperCount,countMatrix,segments] = samplemoves(stopperCount,arenaAxes,startCoords,headings,countMatrix,nSegments,maxRecursion,varargin)
% sample movements from the given starting position.
% avoid consecutive movements of same or opposite direction.
% caution... inputs with impossible values will result in endless recursion.
%
% Inputs
%       stopperCount 
%                   : integer, counting the number of recursion 
%                     (how many times the dead end has been reached)
%       arenaAxes 
%                   : 1 x 2 array,  
%                     lengths of the horizontal and the vertical axes 
%                     of the rectangular arena, in that order
%       startCoords 
%                   : 1 x 2 array, 
%                     Cartesian coordinates of starting postion 
%                     in the order of x and y
%       headings  
%                   : 1 x nHeadings array, 
%                     possible movement directions in Degrees
%       countMatrix     
%                   : 1 x nHeadings array of integers, 
%                     remaining number of available samples per heading  
%       nSegments   
%                   : integer, 
%                     number of segments in the trajectory 
%
%   
%
% Outputs
%       stopperCount
%                   : integer, counting the number of recursion 
%                     (how many times the dead end has been reached)
%       countMatrix     
%                   : 1 x nHeadings array of integers, 
%                     remaining number of available samples per heading                              
%       segments
%                   : struct with fields
%                       directions  1 x nSegments array 
%                       distances   1 x nSegments array 
%                       starts      2 x nSegments array 
%                       ends        2 x nSegments array
%
%
% Last edit         : 2020.09.04
% Author            : Sein Jeung
%--------------------------------------------------------------------------

% if a heading is given as an input argument, 
% control so that the first segment of the trajectory is not parallel to
% whatever came before it and turn angles should be under 120 degrees. 
if isempty(varargin{1})
    turnAngleControl        = 0; 
else
    turnAngleControl        = 1; 
    encodingPathHeading     = varargin{1}; 
end

% recursivity control
stopper             = 0;

% upper and lower bound of the interval of the uniform distribution
% from which the segment lengths will be sampled
moveLength          = [5 7];

% arena dimensions 
arenaH              = arenaAxes(1); 
arenaV              = arenaAxes(2); 

% starting point
x                   = startCoords(1);
y                   = startCoords(2);

% center
centerX             = 0; 
centerY             = 0; 

% initialize the count matrix to be updated
countMatrixNew      = countMatrix;

% add the angular offset corresponding to the trial
moveDirection           = headings;

% construct a movement pool
%--------------------------------------------------------------------------
% list all combinations of movement direction and distances
movePool            = [moveDirection; ...
                       repmat(max(moveLength),1,numel(moveDirection))];

% add updates in x, y coordinates according to each heading and maximum
% distance 
movePool            = [movePool; cosd(movePool(1,:)).*movePool(2,:); sind(movePool(1,:)).*movePool(2,:)];

% sampled directions within a trial
sampledDirections   = [];

% sampled distances within a trial
sampledDistances    = [];

% dot coordinates of the starting points of segments
segmentStartPoints   = [];

% dot coordinates of the ending points of segments
segmentEndPoints     = [];

% sample from the movement pool
%--------------------------------------------------------------------------
for s = 1:nSegments
    
    % define a vector storing the remaining number of samples for each direction
    sampleCount         = countMatrixNew(:);
    
    % updated coordinates defining position bounds (farthest they can go) 
    positionBounds     = [movePool(3,:) + x; movePool(4,:) + y];
    
    % probability weights for available moves
    probabilityWeights  = sampleCount;
    
    % assign 0 probability to moves that go beyond the area
    probabilityWeights((positionBounds(1,:)-centerX).^2 > ((arenaH-centerX)/2)^2) = 0;
    probabilityWeights((positionBounds(2,:)-centerY).^2 > ((arenaV-centerY)/2)^2) = 0;
    
    % control turn angles (only in the navigation phase)
    if turnAngleControl == 1 
        
        if s == 1
            % the first segment of a navigation path
            % can't have the same orientation as the encoding path
            for i = 1:numel(probabilityWeights)
                if movePool(1,i) == encodingPathHeading || ...
                    movePool(1,i) == encodingPathHeading + 180 || ...
                    movePool(1,i) == encodingPathHeading - 180
                
                    % assign 0 probability to segments parallel to 
                    % encoding path 
                    probabilityWeights(i) = 0;
                end
            end
            
        else
                
            for i = 1:numel(probabilityWeights)
                if movePool(1,i) == sampledDirections(s-1)|| ...
                    min([abs(movePool(1,i) - sampledDirections(s-1)), 360-abs(movePool(1,i) - sampledDirections(s-1))]) > 120
                    
                    % assign 0 probability to segments with turn angles 
                    % that are too large (above 120 degrees)
                    probabilityWeights(i) = 0;
                end
            end
            
        end
        
    end
    
    % check if the sum of probabilities is 0
    if sum(countMatrixNew(:)) > 0 && sum(probabilityWeights) <= 0
        stopper = 1;
        break
    end
    
    % starting point of the segment
    segmentStartPoints     = [segmentStartPoints, [x;y]];
    
    % normalize probabilities
    probabilityWeights  = probabilityWeights/sum(probabilityWeights);

    % sample one of the admissible moves according to the given weights
    index               = randsample(numel(probabilityWeights),1,true,probabilityWeights);
    
    % sampled heading 
    sampledHeading      = moveDirection(index); 
    
    % sample length from a uniform distibution 
    sampledLength       = rand(1,1)*(max(moveLength) - min(moveLength)) + min(moveLength); 
    
    % updated the coordinates given the heading and the length 
    updatedPosition     = [x + cosd(sampledHeading)*sampledLength; y + sind(sampledHeading)*sampledLength];
    
    % update x, y coordinates and store in the output matrix
    x                   = updatedPosition(1);
    y                   = updatedPosition(2);
    segmentEndPoints    = [segmentEndPoints, [x;y]];
    sampledDirections   = [sampledDirections movePool(1,index)];
    sampledDistances    = [sampledDistances sampledLength];
    
    % count down the remaining sample number for the corresponding direction
    countMatrixNew(index) = countMatrixNew(index)-1;
    
end

% output struct 
segments.directions     = sampledDirections;
segments.distances      = sampledDistances;
segments.starts         = segmentStartPoints;
segments.ends           = segmentEndPoints;


% Recursion control 
%--------------------------------------------------------------------------

% if the loop terminated without a break, update the count matrix
if stopper == 0
    countMatrix = countMatrixNew;
    
    % if the loop was broken in the middle, count the number of breaks
else
    stopperCount = stopperCount + 1;
    
    % repeat the process if maximum number of recursion is not reached
    if stopperCount < maxRecursion
        [stopperCount,countMatrix,segments] = samplemoves(stopperCount,arenaAxes,startCoords,headings,countMatrix,nSegments,maxRecursion, varargin{1});
        
        % otherwise terminate the process
    else
        stopperCount = maxRecursion;
    end
end

end