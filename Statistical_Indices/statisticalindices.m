function statisticalindices
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
    subjects(subjects == excludedSubject(i)) = [];
   
end

% number of subjects included in the analysis
nSubject        = nSubjectTotal - numel(excludedSubject); 

% length of the rectangular arena axes in meters
xArenaAxis              = 26;
yArenaAxis              = 16;

% segment lengths
segmentLengthArray      = [3 3.5 4 4.5 5 5.5 6]*2;

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
%                        Position Heatmap
%--------------------------------------------------------------------------
% bin numbers
xBins       = 40;
yBins       = 20;

% bin edges along x, y axis
xPoints = linspace(-xArenaAxis/2, xArenaAxis/2, xBins);
yPoints = linspace(-yArenaAxis/2, yArenaAxis/2, yBins);

% sample points
x = allTrajectories(:,13);
y = allTrajectories(:,14);

%figure
%scatter(x,y,'filled','MarkerFaceAlpha',.3)

% bin the data
N = histcounts2(y(:), x(:), yPoints, xPoints);

% plot heatmap
figure
imagesc(xPoints, yPoints, N);
set(gca, 'XLim', xPoints([1 end]), 'YLim', yPoints([1 end]), 'YDir', 'normal');
title('Visited points aggregated over all trials','FontSize',15)
colorbar

%% 
%--------------------------------------------------------------------------
%                             Moran's I 
%--------------------------------------------------------------------------

% % mean number of data points per bin 
% meanHist            = mean(N(:));
% 
% % number of spatial units 
% nUnits              = numel(N);
% 
% % weight matrix initialization
% weightMatrix        = zeros(nUnits,nUnits);
% 
% points = [reshape(N,[nUnits,1]), repelem([1:size(N,2)]',size(N,1)), repmat([1:size(N,1)]',size(N,2),1)];
% csvwrite('sample.csv',points)
% 
% for i = 1:nUnits 
%     for j = 1:nUnits
%         
%         [a,b] = ind2sub(size(N),i);
%         [c,d] = ind2sub(size(N),j);
%         
%         if (a-c)^2+(b-d)^2 == 1
%             weightMatrix(i,j) = 1;
%         elseif abs(a-c) == 1 && abs(b-d) == 1
%             weightMatrix(i,j) = 1;
%         end
%         
%     end
% end
% 
% % row standardization of weight matrix
% for i = 1:nUnits
%    if sum(weightMatrix(i,:)) > 0
%     weightMatrix(i,:) = weightMatrix(i,:)/sum(weightMatrix(i,:));
%    end    
% end
% 
% csvwrite('weight.csv',weightMatrix,1)
% 
% sumSq = 0;
% ijSumProd = 0;
% 
% for i = 1:nUnits
%     
%     jSumProd = 0; 
%     
%     for j = 1:nUnits
%             jSumProd = jSumProd + weightMatrix(i,j)*(N(i)- meanHist)*(N(j)- meanHist);     
%     end
%     
%     ijSumProd = ijSumProd + jSumProd;
%     sumSq = sumSq + (N(i) - meanHist)^2;
% 
% end
% weightMatrixSum     = sum(weightMatrix(:));
% moransI             = nUnits/weightMatrixSum*ijSumProd/sumSq
% 
% % calculate the z-score
% %--------------------------------------------------------------------------
% 
% % normalize N 
% % expected value under null hypothesis
% nullEV              = -1/(nUnits-1);
% 
% % variance 
% sumWeights = 0;
% 
% for i = 1:nUnits
%    for j = 1:nUnits
%         sumWeights = sumWeights + (weightMatrix(i,j)+ weightMatrix(j,i))^2;
%    end
% end
% 
% S1 = sumWeights/2; 
% S2 = 0;
% 
% for i = 1:nUnits
%    ijSum = 0;
%    jiSum = 0;
%    for j = 1:nUnits
%         ijSum = ijSum + weightMatrix(i,j);
%         jiSum = jiSum + weightMatrix(j,i);
%    end
%    S2 = S2 + (ijSum + jiSum)^2; 
% end
% 
% quartSum    = 0;
% quadSum     = 0; 
% 
% for i = 1:nUnits
%     diff        = N(i)-meanHist;
%     quartSum    = quartSum + diff^4;
%     quadSum     = quadSum + diff^2;
% end
% 
% S3 = nUnits^(-1)*quartSum/(nUnits^(-1)*quadSum)^2;
% 
% 
% S4 = (nUnits^2-3*nUnits+3)*S1 - nUnits*S2 + 3*weightMatrixSum^2;
% 
% S5 = (nUnits^2-nUnits)*S1 - 2*nUnits*S2 + 6*weightMatrixSum^2;
% 
% % variance of Moran's I 
% varI = (nUnits*S4 - S3*S5)/((nUnits-1)*(nUnits-2)*(nUnits-3)*weightMatrixSum^2)-nullEV^2;
% zScore = (moransI - nullEV)/sqrt(varI);
% p = erfc(abs(moransI - nullEV)/sqrt(2*varI))


%%
% histogram heatmap run-by-run
%--------------------------------------------------------------------------

figure
% iterate over runs
for r = 1:6
    subplot(2,3,r)
    x = allTrajectories(mod(allTrajectories(:,2),6) == mod(r,6),13);
    y = allTrajectories(mod(allTrajectories(:,2),6) == mod(r,6),14);
    
    % bin the data
    N = histcounts2(y(:), x(:), yPoints, xPoints);
    
    % plot heatmap
    imagesc(xPoints, yPoints, N);
    set(gca, 'XLim', xPoints([1 end]), 'YLim', yPoints([1 end]), 'YDir', 'normal');
    title(['Run ' num2str(r)])
end

% histogram heatmap trial-by-trial
%--------------------------------------------------------------------------

figure
% iterate over trials
for t = 1:6
    subplot(2,3,t)
    x = allTrajectories(mod(allTrajectories(:,3),6) == mod(t,6),13);
    y = allTrajectories(mod(allTrajectories(:,3),6) == mod(t,6),14);
    
    % bin the data
    N = histcounts2(y(:), x(:), yPoints, xPoints);
    
    % plot heatmap
    imagesc(xPoints, yPoints, N);
    set(gca, 'XLim', xPoints([1 end]), 'YLim', yPoints([1 end]), 'YDir', 'normal');
    title(['Trial ' num2str(t)])
end

%%
%--------------------------------------------------------------------------
%                         Segment Lengths
%--------------------------------------------------------------------------
figure
histogram(allTrajectories(:,8),[segmentLengthArray max(segmentLengthArray)+.5])
title('Histogram Segment Lengths')

% segment length histogram run-by-run
%--------------------------------------------------------------------------
figure

% iterate over runs 
for r = 1:6
    subplot(2,3,r)
    histogram(allTrajectories(mod(allTrajectories(:,2),6) == mod(r,6),8),[segmentLengthArray max(segmentLengthArray)+.5])
    title(['Run ' num2str(r)])
end

% segment length histogram trial-by-trial
%--------------------------------------------------------------------------
figure

% iterate over trials 
for t = 1:6
    subplot(2,3,t)
    histogram(allTrajectories(mod(allTrajectories(:,3),6) == mod(t,6),8),[segmentLengthArray max(segmentLengthArray)+.5])
    title(['Trial ' num2str(t)])
end


%%
%--------------------------------------------------------------------------
%                           Run directions
%--------------------------------------------------------------------------

% specify bin edges for polar histograms
polarBinEdges = deg2rad([0 60 120 180 240 300 360]);

figure
polarhistogram(deg2rad(allTrajectories(:,7)),polarBinEdges,'FaceColor','blue','FaceAlpha',.3)
title('Histogram Run Directions')

% run direction histogram run-by-run
%--------------------------------------------------------------------------
figure
% iterate over runs 
for r = 1:6
    subplot(2,3,r)
    polarhistogram(deg2rad(allTrajectories(mod(allTrajectories(:,2),6) == mod(r,6),7)),polarBinEdges,'FaceColor','red','FaceAlpha',.3)
    rlim([0 30])
    title(['Run ' num2str(r)])
end

% run direction histogram trial-by-trial
%--------------------------------------------------------------------------
figure

% iterate over trials 
for t = 1:6
    subplot(2,3,t)
    polarhistogram(deg2rad(allTrajectories(mod(allTrajectories(:,3),6) == mod(t,6),7)),polarBinEdges,'FaceColor','yellow','FaceAlpha',.3)
    rlim([0 30])
    title(['Trial ' num2str(t)])
end

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
% [p,tbl,stats] = anova1(y)
% xlabel('Directions in Degrees') 
% ylabel('Mean Sample Number per Participant')
% xticklabels({'0','60','120','180','240','300'})
% 
% title('Sample Number per Direction')

end