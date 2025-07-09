
function [offset,beta1,beta2,phi,p] = periodicGLM(magnitude, symmetry)

% Fitting GLM
%--------------------------------------------------------------------------
symmetry = 4; 

% extract the heading
theta = deg2rad(0:10:350);

% sine and cosine of movement direction
cosTheta = cos(symmetry*theta);
sinTheta = sin(symmetry*theta);

% Arrange the data in table format
dataTable = array2table([cosTheta',sinTheta',magnitude'], 'VariableNames',{'cos_theta','sin_theta','magnitude'});

% model specification
modelspec = 'linear';

% Fit the general linear model to the data
% (by default, the last variable is the response variable)
model                       = fitglm(dataTable, modelspec);
offset                      = model.Coefficients{1,1};
beta1                       = model.Coefficients{2,1};
beta2                       = model.Coefficients{3,1};
phi                         = atan2(beta2,beta1)/symmetry;
p                           = coefTest(model);

end