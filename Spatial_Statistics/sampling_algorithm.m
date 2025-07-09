% script for testing algorithms

% default figure size 
set(groot, 'defaultFigurePosition',[500 500 800 500])

%headings        = 0:60:300; 
%headings        = 0:45:315; 
headings        = 0:10:350; 
%headings        = 0:10:350; 

length          = 5; 
xyBounds        = [-13 13; -8 8]; 
noiseRange      = 2; 

% 1. unbounded, noiseless
%--------------------------------------------------------------------------
%sampling_unbounded(headings, length, [], 0);

% 2. unbounded, noisy
%--------------------------------------------------------------------------
%sampling_bounded(headings, length, [], noiseRange);

% 3. bounded, noiseless
%--------------------------------------------------------------------------
% [admissiblePaths, X, Y] = sampling_bounded(headings, length, xyBounds, 0);

% for i = 1:numel(headings)
%     
%     headingIndices  = find(admissiblePaths == i);
%    
%     % where the targets are 
%     targetPositionsX  = X(headingIndices);
%     targetPositionsY  = Y(headingIndices);
% %     
% %     figure    
% %     rectangle('Position',[xyBounds(1,1) xyBounds(2,1) xyBounds(1,2)-xyBounds(1,1) xyBounds(2,2)-xyBounds(2,1)])
% %     hold on
% % 
% %     
% %     scatter(targetPositionsX,targetPositionsY)
% %    
%     figure
%     binSize  = .5;
%     xEdges   = xyBounds(1,1):binSize:xyBounds(1,2);
%     yEdges   = xyBounds(2,1):binSize:xyBounds(2,2);
%     Edges    = {xEdges yEdges};
%     hist3([targetPositionsX targetPositionsY],'CdataMode','auto','EdgeColor','none','Edges',Edges)
%     
%     xlabel('X')
%     ylabel('Y')
%     colorbar
%     caxis([0, 300])
%     view(2)
% 
% end
% 
% % Make a gif 
% %--------------------------------------------------------------------------
% h = figure;
% axis tight manual % this ensures that getframe() returns a consistent size
% filename = 'Overheadings_nn.gif';
% 
% for i = 1:numel(headings)
%     
%     headingIndices  = find(admissiblePaths == i);
%    
%     % where the targets are 
%     targetPositionsX  = X(headingIndices);
%     targetPositionsY  = Y(headingIndices);
% %     
% %     figure    
% %     rectangle('Position',[xyBounds(1,1) xyBounds(2,1) xyBounds(1,2)-xyBounds(1,1) xyBounds(2,2)-xyBounds(2,1)])
% %     hold on
% % 
% %     
% %     scatter(targetPositionsX,targetPositionsY)
% %    
%     
%     binSize  = .5;
%     xEdges   = xyBounds(1,1):binSize:xyBounds(1,2);
%     yEdges   = xyBounds(2,1):binSize:xyBounds(2,2);
%     Edges    = {xEdges yEdges};
%     hist3([targetPositionsX targetPositionsY],'CdataMode','auto','EdgeColor','none','Edges',Edges)
%     xlabel('X')
%     ylabel('Y')
%     colorbar
%     caxis([0, 300])
%     view(2)
%     title(['Heading ' num2str(headings(i)) ', all steps, no noise'])
%     
%     drawnow
%     % Capture the plot as an image
%     frame = getframe(h);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     % Write to the GIF File
%     if i == 1
%         imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
%     else
%         imwrite(imind,cm,filename,'gif','WriteMode','append');
%     end
% end
% 
% 
% % Make another gif 
% %--------------------------------------------------------------------------
% h = figure;
% axis tight manual % this ensures that getframe() returns a consistent size
% filename = 'Overtime_nn.gif';
% 
% for i = 1:numel(headings)
%     
%     % where the targets are 
%     targetPositionsX  = X(:,i);
%     targetPositionsY  = Y(:,i);
% 
%     %     figure    
% %     rectangle('Position',[xyBounds(1,1) xyBounds(2,1) xyBounds(1,2)-xyBounds(1,1) xyBounds(2,2)-xyBounds(2,1)])
% %     hold on
% % 
% %     
% %     scatter(targetPositionsX,targetPositionsY)
% 
%     % sampled points over time 
%     % subplot(1,2,1)
%     binSize  = .5;
%     xEdges   = xyBounds(1,1):binSize:xyBounds(1,2);
%     yEdges   = xyBounds(2,1):binSize:xyBounds(2,2);
%     Edges    = {xEdges yEdges};
%     hist3([targetPositionsX targetPositionsY],'CdataMode','auto','EdgeColor','none','Edges',Edges)
%     xlabel('X')
%     ylabel('Y')
%     colorbar
%     caxis([0, 300])
%     view(2)    
%     title(['Step ' num2str(i) ', no noise'])
%      
% %     subplot(1,2,2)
% %     % sampled headings over time 
% %     headingsRad         = deg2rad(headings(admissiblePaths(:,i))); 
% %     polarhistogram(headingsRad,numel(headings))
% %     
%     drawnow
%     % Capture the plot as an image
%     frame = getframe(h);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     
%     % Write to the GIF File
%     if i == 1
%         imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
%     else
%         imwrite(imind,cm,filename,'gif','WriteMode','append');
%     end
%     
% end

% 
% % Make gif representing heading  
% %--------------------------------------------------------------------------
% h = figure;
% axis tight manual % this ensures that getframe() returns a consistent size
% filename = 'Overtime_lines_nn.gif';
% 
% for i = 1:numel(headings)
%     
%     % where the targets are 
%     targetPositionsX  = X(:,i);
%     targetPositionsY  = Y(:,i);
% 
%     if i == 1 
%         pathStartX  = zeros(numel(targetPositionsX),1);
%         pathStartY  = zeros(numel(targetPositionsX),1); 
%     else
%         pathStartX  = X(:,i-1);
%         pathStartY  = Y(:,i-1);
%     end 
%     
%     hold off
%     for j = randsample(numel(targetPositionsX),1000)'
%         plot([pathStartX(j) targetPositionsX(j)], [pathStartY(j) targetPositionsY(j)])
%         hold on
%         
%     end
%     
%     title(['Step ' num2str(i) ', no noise'])
%     
%     drawnow
%     
%     % Capture the plot as an image
%     frame = getframe(h);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     
%     % Write to the GIF File
%     if i == 1
%         imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
%     else
%         imwrite(imind,cm,filename,'gif','WriteMode','append');
%     end
%     
% 
% end

% 
% 
% frequencyMatrix = zeros(36,36);
% 
% % Plot when each heading occurs the most frequently
% for i = 1:numel(headings)
%     
%     for j = 1:numel(headings)
%         
%         frequencyMatrix(i,j) = nnz(admissiblePaths(:,j) == i);
%         
%     end
%     
% end
% 
% figure 
% 
% meanfrequency = mean(mean(frequencyMatrix,2));
% 
% heatmap(frequencyMatrix, 'Colormap',parula)
% caxis([meanfrequency*0.9 meanfrequency*1.1])
% set(gcf,'position',[500,500,700,500])
% 
% %--------------------------------------------------------------------------
% % 4. bounded, noisy
% %--------------------------------------------------------------------------
 [admissiblePaths, X, Y] = sampling_bounded(headings, length, xyBounds, noiseRange);


% for i = 1:numel(headings)
%     
%     headingIndices  = find(admissiblePaths == i);
%    
%     % where the targets are 
%     targetPositionsX  = X(headingIndices);
%     targetPositionsY  = Y(headingIndices);
% %     
% %     figure    
% %     rectangle('Position',[xyBounds(1,1) xyBounds(2,1) xyBounds(1,2)-xyBounds(1,1) xyBounds(2,2)-xyBounds(2,1)])
% %     hold on
% % 
% %     
% %     scatter(targetPositionsX,targetPositionsY)
% %    
%     figure
%     binSize  = .5;
%     xEdges   = xyBounds(1,1):binSize:xyBounds(1,2);
%     yEdges   = xyBounds(2,1):binSize:xyBounds(2,2);
%     Edges    = {xEdges yEdges};
%     hist3([targetPositionsX targetPositionsY],'CdataMode','auto','EdgeColor','none','Edges',Edges)
%     
%     xlabel('X')
%     ylabel('Y')
%     colorbar
%     caxis([0, 300])
%     view(2)
% 
% end

% Make a gif 
%--------------------------------------------------------------------------
h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'Overheadings_yn.gif';

for i = 1:numel(headings)
    
    headingIndices  = find(admissiblePaths == i);
   
    % where the targets are 
    targetPositionsX  = X(headingIndices);
    targetPositionsY  = Y(headingIndices);
%     
%     figure    
%     rectangle('Position',[xyBounds(1,1) xyBounds(2,1) xyBounds(1,2)-xyBounds(1,1) xyBounds(2,2)-xyBounds(2,1)])
%     hold on
% 
%     
%     scatter(targetPositionsX,targetPositionsY)
%    
    
    binSize  = .5;
    xEdges   = xyBounds(1,1):binSize:xyBounds(1,2);
    yEdges   = xyBounds(2,1):binSize:xyBounds(2,2);
    Edges    = {xEdges yEdges};
    hist3([targetPositionsX targetPositionsY],'CdataMode','auto','EdgeColor','none','Edges',Edges)
    xlabel('X')
    ylabel('Y')
    colorbar
    caxis([0, 300])
    view(2)
    title(['Heading ' num2str(headings(i)) ', all steps, noisy'])
   
    drawnow
    
    % Capture the plot as an image
    frame = getframe(h);
    im = frame2im(frame);
    
    imwrite(im,[ 'byheading\byheading' num2str(i) '.png']);
    [imind,cm] = rgb2ind(im,256);
    
    % Write to the GIF File
    if i == 1
        imwrite(imind,cm,filename,'gif','Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
    
end

Rs = [];
Pvals = [];
% Make another gif 
%--------------------------------------------------------------------------
h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'Overtime_yn.gif';

for i = 1:numel(headings)
    
    % where the targets are 
    targetPositionsX    = X(:,i);
    targetPositionsY    = Y(:,i);
    
    if numel(targetPositionsX) > 1000
    
        % calculate the uniformity index
        thinnedX            = targetPositionsX(randsample(numel(targetPositionsX),1000));
        thinnedY            = targetPositionsY(randsample(numel(targetPositionsX),1000));
        [R, Pval] = csrtest([thinnedX thinnedY],xyBounds);
    
    else 
        [R, Pval] = csrtest([targetPositionsX targetPositionsY],xyBounds);
    end
    
    binSize  = .5;
    xEdges   = xyBounds(1,1):binSize:xyBounds(1,2);
    yEdges   = xyBounds(2,1):binSize:xyBounds(2,2);
    Edges    = {xEdges yEdges};
    hist3([targetPositionsX targetPositionsY],'CdataMode','auto','EdgeColor','none','Edges',Edges)
    xlabel('X')
    ylabel('Y')
    colorbar
    caxis([0, 300])
    view(2)    
    title(['Step ' num2str(i) ', noisy, R = ' num2str(R) ', p = ' num2str(Pval)])
    
    drawnow
    
    % Capture the plot as an image
    frame = getframe(h);
    
    im = frame2im(frame);
        
    imwrite(im,[ 'bystep\bystep' num2str(i) '.png']);
    
    [imind,cm] = rgb2ind(im,256);

    
    % Write to the GIF File
    if i == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
    
Rs  = [Rs R];
Pvals = [Pvals Pval]; 

end

%--------------------------------------------------------------------------
figure
bar(Rs,'EdgeColor','none','FaceColor',[.3 .5 .5])
title(['Index of spatial randomness, 1000 points at each step'])
xlabel('Time step')
ylabel('Clark-Evans Index')

frequencyMatrix = zeros(36,36);

% Plot when each heading occurs the most frequently
for i = 1:numel(headings)
    
    for j = 1:numel(headings)
        
        frequencyMatrix(i,j) = nnz(admissiblePaths(:,j) == i);
        
    end
    
end

figure 

heatmap(frequencyMatrix, 'Colormap',parula)

title('Occurence of each heading at each time step')
xlabel('Time step')
ylabel('Heading')

% Visualize and test distances to corners and boundaries 
distances(admissiblePaths, X, Y, xyBounds);


% hist3([targetPositionsX targetPositionsY],'CdataMode','auto')
% xlabel('X')
% ylabel('Y')
% colorbar
% %caxis([0, 300])
% view(2)
% title(['Step ' num2str(i) ', noisy'])


% figure
% for i = 50   
%     plot([0 X(i,:)],[0 Y(i,:)])
%     hold on
% 
% end

%
% N = hist3([targetPositionsX targetPositionsY],'CdataMode','auto','Edges',Edges)
% figure
% plot(unique(N))