function sampling_heatmap(step,nPaths,points,xyBounds)


    binSize  = 0.2; 
    nBinsX   = 50; 
    nBinsY   = 50; 
    
    figure
    
    % unbounded case
    if  isempty(xyBounds)
        hist3(points,[nBinsX nBinsY],'CdataMode','auto','EdgeColor','none')
    % bounds are specified
    else
        xEdges   = xyBounds(1,1):binSize:xyBounds(1,2);
        yEdges   = xyBounds(2,1):binSize:xyBounds(2,2);
        Edges    = {xEdges yEdges};
        hist3(points,'CdataMode','auto','EdgeColor','none','Edges',Edges)
    end
    
    xlabel('X')
    ylabel('Y')
    colorbar
    view(2)
    
    %scatter(xUpdate,yUpdate,'MarkerEdgeColor','none','MarkerFaceColor',[0 i*1/nSegments 0])
    hold on
    
    title([num2str(nPaths) ' paths at step ' num2str(step)])

   
end 