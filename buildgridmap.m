function  gridMap = buildgridmap(xCenter,yCenter,xGridMapSize,yGridMapSize,angularOffset)
% creates a grid map with given center and size 
%
% Inputs
%       xCenter
%                   : x coordinate of the center
%       yCenter
%                   : y coordinate of the center
%       xGridMapSize
%                   : length of the horizontal side of rectangular grid map
%       yGridMapSize
%                   : length of the vertical side of rectangular grid map
%       angularOffset
%                   : rotates the map about the center
%
% Outputs
%       gridMap
%                   : vertex coordinates of the grid map
%
% author: Sein Jeung
%--------------------------------------------------------------------------


% number of grids in one row of the map, defines the resolution
nGridsHorizontal        = 16;

% grid size given as the length of one side of the triangle
gridSide                = xGridMapSize/nGridsHorizontal;

% height of the triangle
gridHeight              = sqrt(3)/2*gridSide;

% number of grids in one column of the map
nGridsVertical          = ceil(yGridMapSize/gridHeight);

% for symmetric centering, half of nGridsVertical is odd if nGridsHorizontal is odd, otherwise even
if  mod(nGridsHorizontal,2)     == 0
    nGridsVertical      = nGridsVertical + 4 - mod(nGridsVertical,4);
else
    nGridsVertical      = nGridsVertical + 6 - mod(nGridsVertical,4);
end

% x coordinate templates of the grid vertices
xVertices1              = [-.5*xGridMapSize + xCenter:gridSide:.5*xGridMapSize + xCenter];
xVertices2              = linspace(-.5*xGridMapSize + xCenter + .5*gridSide,.5*xGridMapSize + xCenter -.5*gridSide ,numel(xVertices1)-1);

% y coordinates templates of the grid vertices
yVertices1              = -.5*nGridsVertical*gridHeight*ones(1,numel(xVertices1)) + yCenter;
yVertices2              = -.5*nGridsVertical*gridHeight*ones(1,numel(xVertices2)) + yCenter + gridHeight;

% array of x coordinates for the whole map
xVertices               = repmat([xVertices1 xVertices2],1,nGridsVertical/2);
xVertices               = [xVertices xVertices1];

% array of y coordinates for the whole map
%--------------------------------------------------------------------------
% initialize an empty array of y coordinates
yVertices               = [];

% iterate over half the number of grids in one column 
for i = 1:nGridsVertical/2
    
    % append y coordinates row by row in steps of two
    yVertices           = [yVertices [yVertices1 yVertices2] + (i-1)*gridHeight*2];
end
yVertices               = [yVertices yVertices1 + nGridsVertical*gridHeight];


% rotation matrix
R = [cosd(angularOffset) -sind(angularOffset); sind(angularOffset) cosd(angularOffset)];

% rotate around the center
center                  = repmat([xCenter;yCenter],1,length(xVertices));
v                       = [xVertices; yVertices];
vo                      = R*(v-center) + center;
vo                      = vo';
TRI                     = delaunayTriangulation(vo(:,1),vo(:,2));

% output assignment
gridMap         = TRI;

end