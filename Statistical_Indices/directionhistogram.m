function directionhistogram
% Last edit         : 2018.12.12
% Author            : Sein Jeung
%-------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Saved data
%
%       filename:   'traj_subXX.m'
%       variable:   trajectoryData
%       column 1:   participant ID
%       column 2:   run index
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
%      column 15:   turn angle 
%
%--------------------------------------------------------------------------
%%
%--------------------------------------------------------------------------
%                           Specifications
%--------------------------------------------------------------------------

% total number of subjects 
nSubjectTotal   = 2; 

% number of sessions per subject
%nSession        = 2; 

% vector of ID of subjects to be excluded 
excludedSubject = [1]; 

% initialize vector of subject ID 
subjects        = 1:nSubjectTotal;

% iterate over excluded subjects
for i = 1:numel(excludedSubject)
    
    % remove the ID from subject ID vector
    subjects(find(subjects == excludedSubject(i))) = [];
   
end

% number of subjects included in the analysis
nSubject        = nSubjectTotal - numel(excludedSubject); 

% angular offsets
angularOffset      = [0 10 20 30 40 50];

% bsae directions 
baseDirection      = [0 60 120 180 240 300];

%%
%--------------------------------------------------------------------------
%                  Load and Concatenate Data Matrix
%--------------------------------------------------------------------------

% initialize an empty matrix 
allTrajectories = [];

for n = subjects
    %for s = 1:nSession
        load(fullfile(pwd,'traj_pilot',['traj_sub' num2str(n,'%02d')]))
        
        % append data matrix
        allTrajectories = [allTrajectories; trajectoryData];
    %end
end


%%
%--------------------------------------------------------------------------
%                       Base run directions 
%--------------------------------------------------------------------------

% specify bin edges for polar histograms
polarBinEdges = deg2rad([0 60 120 180 240 300 360]);

figure
polarhistogram(deg2rad(allTrajectories(:,7)),polarBinEdges,'FaceColor','blue','FaceAlpha',.3)
title('Histogram Base Directions')

% % run direction histogram run-by-run
% %--------------------------------------------------------------------------
% figure
% % iterate over runs 
% for r = 1:6
%     subplot(2,3,r)
%     polarhistogram(deg2rad(allTrajectories(mod(allTrajectories(:,2),6) == mod(r,6),7)),polarBinEdges,'FaceColor','red','FaceAlpha',.3)
%     rlim([0 300])
%     title(['Run ' num2str(r)])
% end
% 
% % run direction histogram trial-by-trial
% %--------------------------------------------------------------------------
% figure
% 
% % iterate over trials 
% for t = 1:6
%     subplot(2,3,t)
%     polarhistogram(deg2rad(allTrajectories(mod(allTrajectories(:,3),6) == mod(t,6),7)),polarBinEdges,'FaceColor','yellow','FaceAlpha',.3)
%     rlim([0 300])
%     title(['Trial ' num2str(t)])
% end

% 1-way ANOVA
%--------------------------------------------------------------------------

% % initialize data matrix
% y       = [];
% 
% for n = 1:nSubject 
% y       = [y; sum(allTrajectories(:,1) == subjects(n) & allTrajectories(:,7) == 0),...
%             sum(allTrajectories(:,1) == subjects(n) & allTrajectories(:,7) == 60),...
%             sum(allTrajectories(:,1) == subjects(n) & allTrajectories(:,7) == 120),...
%             sum(allTrajectories(:,1) == subjects(n) & allTrajectories(:,7) == 180),...
%             sum(allTrajectories(:,1) == subjects(n) & allTrajectories(:,7) == 240),...
%             sum(allTrajectories(:,1) == subjects(n) & allTrajectories(:,7) == 300)];
% end
% 
%[p,tbl,stats] = anova1(y)


%%
%--------------------------------------------------------------------------
%                        Base + Angular Offset
%--------------------------------------------------------------------------

% specify bin edges for polar histograms
polarBinEdges = deg2rad(0:10:360);

figure
h = polarhistogram(deg2rad(allTrajectories(:,7)+allTrajectories(:,5)),polarBinEdges,'FaceColor','blue','FaceAlpha',.3);
rlim([0 size(allTrajectories,1)/36 + 2]);
title('Histogram Run Directions')
% 
% % run direction histogram run-by-run
% %--------------------------------------------------------------------------
% figure
% % iterate over runs 
% for r = 1:6
%     subplot(2,3,r)
%     polarhistogram(deg2rad(allTrajectories(mod(allTrajectories(:,2),6) == mod(r,6),7)...
%         + allTrajectories(mod(allTrajectories(:,2),6) == mod(r,6),5)),...
%         polarBinEdges,'FaceColor','red','FaceAlpha',.3)
%     rlim([0 70])
%     title(['Run ' num2str(r)])
% end
% 
% % run direction histogram trial-by-trial
% %--------------------------------------------------------------------------
% figure
% 
% % iterate over trials 
% for t = 1:6
%     subplot(2,3,t)
%     polarhistogram(deg2rad(allTrajectories(mod(allTrajectories(:,3),6) == mod(t,6),7)...
%         + allTrajectories(mod(allTrajectories(:,3),6) == mod(t,6),5)),...
%         polarBinEdges,'FaceColor','yellow','FaceAlpha',.3)
%     rlim([0 70])
%     title(['Trial ' num2str(t)])
% end
% 

%%
%--------------------------------------------------------------------------
%                        Segment lengths
%--------------------------------------------------------------------------

meanLengths   = zeros(numel(angularOffset)*numel(baseDirection),0);
for a = 1:numel(angularOffset)
    for b = 1:numel(baseDirection)
        lengths     = allTrajectories(allTrajectories(:,5) == angularOffset(a) &...
                       allTrajectories(:,7) == baseDirection(b),8);
        meanLengths(a + (b-1)*numel(angularOffset)) = mean(lengths);
    end
end 

Options = {'labels',{'0','180','90','270'},...
    'cMap','summer','freqlabelangle',45,'lablegend','segment lengths',...
    'legendvariable','length'}; 

WindRose(allTrajectories(:,5)+allTrajectories(:,7),allTrajectories(:,8),Options)
title('')
%title('Polar histogram of segment lengths sampled','FontSize',15)


% % run direction histogram run-by-run
% %--------------------------------------------------------------------------
% figure
% % iterate over runs 
% for r = 1:6
%     subplot(2,3,r)
%     polarhistogram(deg2rad(allTrajectories(mod(allTrajectories(:,2),6) == mod(r,6),7)...
%         + allTrajectories(mod(allTrajectories(:,2),6) == mod(r,6),5)),...
%         polarBinEdges,'FaceColor','red','FaceAlpha',.3)
%     rlim([0 60])
%     title(['Run ' num2str(r)])
% end
% 
% % run direction histogram trial-by-trial
% %--------------------------------------------------------------------------
% figure
% 
% % iterate over trials 
% for t = 1:6
%     subplot(2,3,t)
%     polarhistogram(deg2rad(allTrajectories(mod(allTrajectories(:,3),6) == mod(t,6),7)...
%         + allTrajectories(mod(allTrajectories(:,3),6) == mod(t,6),5)),...
%         polarBinEdges,'FaceColor','yellow','FaceAlpha',.3)
%     rlim([0 60])
%     title(['Trial ' num2str(t)])
%end

%%
%--------------------------------------------------------------------------
%                        Direction*Length
%--------------------------------------------------------------------------

% direction bin edges
directionBinEdges = deg2rad(0:10:360);

figure
polarhistogram('BinEdges',directionBinEdges,'BinCounts',meanLengths)
title('Lengths')

figure
polarhistogram('BinEdges',directionBinEdges,'BinCounts',(h.BinCounts/mean(h.BinCounts)).*meanLengths)
title('Lengths weighted by direction count')
end