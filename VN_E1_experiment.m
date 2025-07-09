function VN_E1_experiment(ID)
% construct a trajectory consisting of straight segments 
% with given headings and lengths 
% 
% 
% calls functions 
%                   : buildgridmap
%                     samplemoves
%
% Inputs
%                ID : integer, participant ID 
%
% Outputs
%              none
%
%
% Last edit         : 2019.10.30
% Author            : Sein Jeung
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Saved data
%
%       filename:   'traj_subXX.m'
%       variable:   trajectoryData
%       column 1:   participant ID
%       column 2:   block index
%       column 3:   trial index
%       column 4:   trial length
%       column 5:   angular offset
%       column 6:   segment index
%       column 7:   basis run direction 
%       column 8:   distance
%       column 9:   update in x 
%      column 10:   update in y
%      column 11:   start position x
%      column 12:   start position y
%      column 13:   end position x
%      column 14:   end position y
%      column 15:   segment index (0 for base, 1 for encoding, 2 for end) 
%--------------------------------------------------------------------------

%%
%--------------------------------------------------------------------------
%                           Directory Management  
%--------------------------------------------------------------------------
% participant ID (function input)
participantID           = num2str(ID,'%02d');

% find current directory  
currentDirectory        = pwd;

% directory of the trajectory data matrix
dataDirectory           = fullfile(currentDirectory, 'traj_pilot');

% subject-specific file name
subjectFileName         = fullfile(dataDirectory,['traj_sub' participantID ]);

% directory of the trajectory data matrix
if  ~exist(dataDirectory)
    mkdir(dataDirectory)
end


%--------------------------------------------------------------------------
%                            Parameters 
%--------------------------------------------------------------------------
% headings for encoding path 
encodingHeadings        = 0:5:355;

% headings for pseudorandom navigation 
navigationHeadings      = 0:5:355;

% number of runs 
nRuns                   = 8;

% number of sessions
nSessions               = 1; 

% number of runs in a session 
nRunsInSession          = nRuns/nSessions; 

% number of trials in a run 
nTrials                 = 9; 

% average number of segments in a trial 
nMeanSegments           = 5;

% number of conditions 
% in this case, half of the experiment covers all conditions 
nConditions             = numel(encodingHeadings); 

% number of segments in a session
nSegmentsInSession     = nRunsInSession*nTrials*nMeanSegments; 

% x,y coordinates of the starting position 
startCoords             = [0 0];

% length of the rectangular arena horizonal, vertical axes in meters
arenaAxes               = [26 16]; 

% maximal number of recursion per trial 
maxRecursion            = 500; 


% Initialization 
%--------------------------------------------------------------------------
% initialize a matrix to store all data 
trajectoryData          = NaN(nSegmentsInSession*nSessions,15);

% initialize a run-by-run data array 
directionArray            = {};
distanceArray             = {};
positionArray             = {};
segmentFlagArray          = {};
targetIndexArray          = {};

% store subject ID in the data matrix
trajectoryData(:,1)     = ID; 


%%
%--------------------------------------------------------------------------
%        Randomizing: Sample target position and encoding paths
%--------------------------------------------------------------------------
% write a script that creates a pseudorandom trajectory consisting of
% N segments given the lengths, headings and probability weights
% iterate over half the runs

%%
%--------------------------------------------------------------------------
%             Randomizing: Random Navigation Paths
%--------------------------------------------------------------------------

% iterate over sessions
for session = 1:nSessions
    
    % initilaize repetition control
    repeat                  = 1;
    
    while repeat == 1
        
        % repetition control
        repeat = 0;
        
        % initialize counts for encoding paths
        encodingPathCounter = ones(1,nConditions);
        
        % sample movements to construct the encoding trajectory 
        [stopperCount,countMatrix,encodingSegments] = samplemoves(0,arenaAxes,startCoords,encodingHeadings,encodingPathCounter,numel(encodingHeadings),maxRecursion,[]);
        
        % initialize counts for navigating paths
        navigationPathCounter = (nMeanSegments-1)*ones(1,numel(navigationHeadings));
        
        % randomize segment numbers for random navigation phase in each trial
        nSegmentsArray = repmat([nMeanSegments-2,nMeanSegments-1,nMeanSegments],1,nConditions/3);
        nSegmentsArray = nSegmentsArray(randperm(nConditions));

        % iterate over runs in each half
        for r = 1:nRunsInSession
            
            runIndex        = r + (session-1)*nRunsInSession;   
            
            % find the index of first segment in each run
            if r == 1
                runStartSegment = (session-1)*nSegmentsInSession + 1; 
            else
                runStartSegment = (session-1)*nSegmentsInSession + sum(nSegmentsArray(1:(r-1)*nTrials)) + (r-1)*nTrials + 1; 
            end
            
            % find the index of last segment in each run 
            runEndSegment = (session-1)*nSegmentsInSession + sum(nSegmentsArray(1:r*nTrials)) + r*nTrials;  
            
            % add run number to the data matrix
            trajectoryData(runStartSegment:runEndSegment,2) = r + nRunsInSession*(session-1);
            
            % number of segments per run
            nSegmentsPerRun         = sum(nSegmentsArray((r-1)*nTrials+1:r*nTrials)) + nTrials;
            
            % initialize empty arrays to store segment data in run  
            directionsInRun             = [];
            distancesInRun              = [];
            startPositionsInRun         = [];
            endPositionsInRun           = [];
            segmentFlagsInRun           = [];
            
            % iterate over trials in a run 
            % to construct the random navigation trajectories 
            for t = 1:nTrials
                
                % index of a trial in one whole session 
                trialIndex  = (r-1)*nTrials + t;
                
                % number of segments in a trial 
                nSegments   = nSegmentsArray((r-1)*nTrials + t);
                
                % starting point of the navigation phase
                xStart      = encodingSegments.ends(1,trialIndex);
                yStart      = encodingSegments.ends(2,trialIndex);
                
                % sample movements to construct the navigation trajectory
                [stopperCount,navigationPathCounter,navigatingSegments] = samplemoves(0,arenaAxes,[xStart, yStart],navigationHeadings,navigationPathCounter,nSegments,maxRecursion, encodingSegments.directions(trialIndex));
                
                % store the sampled data in a run data arrays
                directionsInRun         = [directionsInRun, encodingSegments.directions(trialIndex), navigatingSegments.directions];
                distancesInRun          = [distancesInRun, encodingSegments.distances(trialIndex), navigatingSegments.distances];
                startPositionsInRun     = [startPositionsInRun encodingSegments.starts(:,trialIndex) navigatingSegments.starts];
                endPositionsInRun       = [endPositionsInRun encodingSegments.ends(:,trialIndex) navigatingSegments.ends];
                segmentFlagsInRun       = [segmentFlagsInRun, 1, zeros(1,nSegments-1), 2];
                
                % if recursion count reaches the maximum,
                % it is likely that the trial is left with impossible
                % combination of remaining directions and position 
                % restart sampling for the whole session  
                if stopperCount == maxRecursion
                    
                    % restart sampling from the first run
                    repeat = 1;
                    break;
                    
                end
                
            end
            
            % break if impossible 
            if repeat == 1
                break;
            end
              
            % store the sampled run data in cell arrays 
            directionArray(runIndex)            = {directionsInRun};
            distanceArray(runIndex)             = {distancesInRun};
            positionArray(1,runIndex)           = {startPositionsInRun};
            positionArray(2,runIndex)           = {endPositionsInRun};
            segmentFlagArray(runIndex)          = {segmentFlagsInRun };
            targetIndexArray(runIndex)          = {1};
            
            % add one (encoding path) to the segment lengths 
            trialLengths                 = nSegmentsArray((r-1)*nTrials+1:r*nTrials) + 1;
            
            % store trial index
            trajectoryData(runStartSegment:runEndSegment,3)      = repelem([1:nTrials],trialLengths);
            
            % store trial length
            trajectoryData(runStartSegment:runEndSegment,4)      = repelem(trialLengths,trialLengths);
            
            % store angular offset
            trajectoryData(runStartSegment:runEndSegment,5)      = mod(directionsInRun,60);
            
            % store segment index
            trajectoryData(runStartSegment:runEndSegment,6)      = 1:nSegmentsPerRun;
            
        end
    end
end

% store basis movement direction
trajectoryData(:,7)                  = horzcat(directionArray{1,:})'-trajectoryData(:,5);

% store move distances
trajectoryData(:,8)                  = horzcat(distanceArray{1,:})';

% update in x
trajectoryData(:,9)                  = cosd(trajectoryData(:,5)+trajectoryData(:,7)).*trajectoryData(:,8);

% update in y
trajectoryData(:,10)                 = sind(trajectoryData(:,5)+trajectoryData(:,7)).*trajectoryData(:,8);

% convert position array to matrix
positionMatrix                       = cell2mat(positionArray);

% start and end positon coordinates
trajectoryData(:,11:14)              = positionMatrix(1:4,:)';

% store segment flags
trajectoryData(:,15)                 = horzcat(segmentFlagArray{1,:})';



%%
%--------------------------------------------------------------------------
%                               Save Data
%--------------------------------------------------------------------------
save(subjectFileName,'trajectoryData')

% Write csv data to be added to the resource folder in the Unity project 
%--------------------------------------------------------------------------
% waypoints consisting of end points of all segments + one at the
% beginning, guiding the player to the starting point of the first segment
% of each run 
csvdata             = [];

for r = 1:nRuns
    csvdata = [csvdata; trajectoryData(find(trajectoryData(:,2) == r,1,'first') ,11:12) 0; trajectoryData(trajectoryData(:,2) == r ,13:14) trajectoryData(trajectoryData(:,2) == r,15)];
end 

% name the csv file
csvname         = ['segments_sub' num2str(participantID) '.txt' ];
size(csvdata')

% write csv file
dlmwrite(csvname,csvdata');


%% 
%--------------------------------------------------------------------------
%                              Animation
%--------------------------------------------------------------------------
%iterate over runs
for r = 1
    
    % find numer of segments in each trial
    trialLengthsInRun   = trajectoryData(trajectoryData(:,2) == r & trajectoryData(:,15) == 1 , 4);
    
    % angular offsets of the grid template (not relevant in this version)
    angularOffsets      = zeros(1,numel(trialLengthsInRun));
   
    % target index in each trial (all ones in this version, where the first segment is the encoding segment)
    targetIndices       = trajectoryData(trajectoryData(:,2) == r & trajectoryData(:,15) == 1 , 6); 
    
    % animate movement along the trajectory 
    %trajectory_animation(cell2mat(positionArray(:,r)),angularOffsets,trialLengthsInRun',targetIndices,r);
    
end

end