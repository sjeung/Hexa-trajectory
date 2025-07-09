function [admissiblePaths,admissibleHistoryX, admissibleHistoryY] = sampling_bounded(headings, length, xyBounds, noise)


if isempty(xyBounds)
    xBound = inf;
    yBound = inf;
else
    % here it is assumed that x and y bounds are symmetric wrt the origin
    xBound  = xyBounds(1,2);
    yBound  = xyBounds(2,2);
end

nSegments = numel(headings);

%--------------------------------------------------------------------------
% Simulation  
%--------------------------------------------------------------------------

% initialization
%--------------------------------------------------------------------------
admissiblePaths         = [1:nSegments]';
%admissiblePaths         = [1:9]';

for i = 1:nSegments
    
    randomLength = length + rand(1,1)*noise; 
    
    admissiblePoints(i,:)    = [cosd(headings(i))*randomLength sind(headings(i))*randomLength];
end

admissibleHistoryX      = admissiblePoints(:,1);
admissibleHistoryY      = admissiblePoints(:,2);


%--------------------------------------------------------------------------
for i = 2:nSegments
    
    disp(['Simulation step ' num2str(i) ', assessing ' num2str(size(admissiblePaths,1)) ' paths.'])
    
    newPaths        = [];
    newPoints       = [];
    newHistoryX      = [];
    newHistoryY      = [];
    
    for j = 1:nSegments        
        
        numCandidates       = nnz(~any(admissiblePaths == j,2));
        
        candidatePaths      = [admissiblePaths(~any(admissiblePaths == j,2),:) repelem(j,numCandidates,1)];
        
        randomLengths       = repelem(length,numCandidates,1) + rand(numCandidates,1)*noise;
        
        candidatePoints     = admissiblePoints(~any(admissiblePaths == j,2),:) + [cosd(headings(j)).*randomLengths sind(headings(j)).*randomLengths];
        
        candidateHistoryX   = [admissibleHistoryX(~any(admissiblePaths == j,2),:) candidatePoints(:,1)]; 
        
        candidateHistoryY  = [admissibleHistoryY(~any(admissiblePaths == j,2),:) candidatePoints(:,2)]; 
        
        % identify points that go beyond the area 
        rejectedIndices     = find(abs(candidatePoints(:,1)) > xBound | abs(candidatePoints(:,2))> yBound);
        
        % remove the paths 
        candidatePaths(rejectedIndices,:)       = [];
        candidatePoints(rejectedIndices,:)      = []; 
        candidateHistoryX(rejectedIndices,:)    = [];
        candidateHistoryY(rejectedIndices,:)    = [];

        newPaths        = [newPaths;  candidatePaths];
        newPoints       = [newPoints;  candidatePoints];
        newHistoryX     = [newHistoryX; candidateHistoryX];
        newHistoryY     = [newHistoryY; candidateHistoryY];
    
    end
%     
%     % find unique rows only
%     %----------------------------------------------------------------------
%     [newPaths uniqueIndices] = unique(sort(newPaths,2),'rows');
%     
%     newPoints           = newPoints(uniqueIndices,:);
%     newHistoryX         = newHistoryX(uniqueIndices,:);
%     newHistoryY         = newHistoryY(uniqueIndices,:);
%     %----------------------------------------------------------------------
%     
    % This part is a patchy solution
    %----------------------------------------------------------------------
    if size(newPaths,1) > 200000
        
        pruningIndices      = randsample(size(newPaths,1),200000);
        
        newPaths      = newPaths(pruningIndices,:);
        newPoints     = newPoints(pruningIndices,:);
        newHistoryX   = newHistoryX(pruningIndices,:);
        newHistoryY   = newHistoryY(pruningIndices,:);
        
    end

    admissiblePaths     = newPaths;
    admissiblePoints    = newPoints; 
    admissibleHistoryX  = newHistoryX; 
    admissibleHistoryY  = newHistoryY;

    
%     %----------------------------------------------------------------------
%     % Plotting 
%     %----------------------------------------------------------------------
%     
%     % 1. Positions 
%     %----------------------------------------------------------------------
%     if i == 2 ||  i == 5 || i == 10
%     
%     nPaths      = size(newPoints,1);
%     sampling_heatmap(i,nPaths,newPoints,xyBounds)
%     
%     %figure
%     %scatter(newPoints(:,1),newPoints(:,2))
%     
%     end
%     
%     % 2. Headings at step i
%     %----------------------------------------------------------------------
%     %if i == 10
%         
%         figure     
%         hist(admissiblePaths(:,i),nSegments)
%         title(['Histogram of adimissible headings at step ' num2str(i)])
%     %end
    
end

end