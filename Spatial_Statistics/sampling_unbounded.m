% demo for N = 10, unbounded 

function sampling_unbounded(headings, length, noise)

nSegments   = numel(headings); 
xyBounds    = []; 

for i = 1:nSegments %:ceil(nSegments/2) 
    
 
    sampledIndices  = nchoosek(1:nSegments,i);
    
    xUpdate = [];
    yUpdate = [];
    
    for j = 1:size(sampledIndices,1)
        
        randomLength = length + rand(1,1)*noise;
        
        xUpdate = [xUpdate sum(cosd(headings(sampledIndices(j,:))).*randomLength)];
        yUpdate = [yUpdate sum(sind(headings(sampledIndices(j,:))).*randomLength)];
   
    end
    
    
    points  = [xUpdate',yUpdate'];
    nLocs   = size(unique(points,'rows'),1); 

    step        = i; 
    nPaths      = prod(nSegments-i+1:nSegments);
    
   sampling_heatmap(step,nPaths,points,xyBounds)
     
%     figure
%     
%     hist3(points,[50,50],'CdataMode','auto','EdgeColor','none')
%     xlabel('X')
%     ylabel('Y')
%     colorbar
%     view(2)
% 
%     %scatter(xUpdate,yUpdate,'MarkerEdgeColor','none','MarkerFaceColor',[0 i*1/nSegments 0])
%     hold on
%     
%     title([num2str(nLocs) ' locations, ' num2str(prod(nSegments-i+1:nSegments)) ' paths at step ' num2str(i)])
%     
end



end 




