function [stopperCount,countMatrix,segments] = VN_sample_baseline(stopperCount,xStart,yStart,maxRecursion)
% sample movements from the given starting position.
% avoid consecutive movements of same or opposite direction.
% caution... inputs with impossible values will result in endless recursion.
%
% Inputs
%       xStart
%                   : x coordinate of the starting position
%       yStart
%                   : y coordinate of the starting position
%       moves      
%                   :  struct with fields
%         direction : array of run directions
%            length : array of segment lengths
%           repeats : number of repeats per run direction
%            center : 1X2,[center x coordinate, center y coordinate]
%              axes : 1X2,[length of horizontal axis of elliptical arena,
%                          length of vertical axis of elliptical arena]
%
% Outputs
%       stopper
%                   : recursion control
%                     (repeat if dead end is reached during sampling)
%       segment
%                   : 2X(number of segments), array of points to be visited
%                     on the trajectory
%
% Author            : Sein Jeung
%--------------------------------------------------------------------------

% recursivity control
stopper             = 0;

% segment length
movelength          = 5;

% number of trials
nTrial              = 1;

% starting point
x                   = xStart;
y                   = yStart;

% center
centerX             = 0; 
centerY             = 0; 

% arena dimensions
% Dark Hall physical dimensions             : 15 x 9 
% Tracking space                            : 13 x 8
% Tracking space in the desktop version     : 26 x 16
arenaH              = 24; 
arenaV              = 12; 

% initialize an array of starting points
segmentStartRun     = [];

% initialize an array of ending points
segmentEndRun       = [];

% initialize the count matrix to be updated
countMatrixNew      = ones(2,36);

% sampled directions in the run
sampledDirectionsRun   = [];

% sampled distancess in the run
sampledDistancesRun   = [];

% specification of variable size
%--------------------------------------------------------------------------


for t = 1:nTrial
    
    % add the angular offset corresponding to the trial
    moveDirection       = 10:10:360;

    % construct a movement pool
    %----------------------------------------------------------------------
    % list all combinations of movement direction and distances
    movePool            = [moveDirection; ...
                           repmat(movelength,1,numel(moveDirection))];
    
    % add updates in x, y coordinates according to the movement specified
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
    %----------------------------------------------------------------------
    for s = 1:36

        % define a vector storing the remaining number of samples for each direction
        sampleCount         = countMatrixNew(t,:);
        
        % updated coordinates after available moves
        updatedPosition     = [movePool(3,:)+x; movePool(4,:)+y];
        
        % probability weights for available moves
        probabilityWeights  = sampleCount;
               
        % assign 0 probability to moves that go beyond the area
        probabilityWeights((updatedPosition(1,:)-centerX).^2 > ((arenaH-centerX)/2)^2) = 0;
        probabilityWeights((updatedPosition(2,:)-centerY).^2 > ((arenaV-centerY)/2)^2) = 0;
        
        
        % assign 0 probability to moves with repeating motion directions or following a move to the opposite direction
        if s >= 2
            for i = 1:numel(probabilityWeights)
                if movePool(1,i) == sampledDirections(s-1)|| ...
                        movePool(1,i) == sampledDirections(s-1) + 180 || ...
                        movePool(1,i) == sampledDirections(s-1) - 180
                    probabilityWeights(i) = 0;
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
        index = randsample(numel(probabilityWeights),1,true,probabilityWeights);
        
        % update x, y coordinates and store in the output matrix
        x                   = updatedPosition(1,index);
        y                   = updatedPosition(2,index);
        segmentEndPoints    = [segmentEndPoints, [x;y]];
        sampledDirections   = [sampledDirections movePool(1,index)];
        sampledDistances    = [sampledDistances movePool(2,index)];
        
        % count down the remaining sample number for the corresponding direction
        countMatrixNew(t,movePool(1,index)/10) = countMatrixNew(t,movePool(1,index)/10)-1;
        
    end
   
    if stopper == 1
        break
    end
    
    % store sampled directions
    sampledDirectionsRun    = [sampledDirectionsRun sampledDirections];
    
    % store sampled distances
    sampledDistancesRun     = [sampledDistancesRun sampledDistances];
    
    % store segment starting points
    segmentStartRun         = [segmentStartRun segmentStartPoints];
    
    % store segment starting points
    segmentEndRun           = [segmentEndRun segmentEndPoints];
    
end

segments.directions     = sampledDirectionsRun;
segments.distances      = sampledDistancesRun;
segments.starts         = segmentStartRun;
segments.ends           = segmentEndRun;

% if the loop terminated without a break, update the count matrix
if stopper == 0
    countMatrix = countMatrixNew;
    
% if the loop was broken in the middle, count the number of breaks
else
    stopperCount = stopperCount + 1;

    % repeat the process if maximum number of recursion is not reached
    if stopperCount < maxRecursion 
       [stopperCount,countMatrix,segments] = VN_sample_baseline(stopperCount,xStart,yStart,maxRecursion);
    
    % otherwise terminate the process
    else
       stopperCount = maxRecursion;
    end
end

end
