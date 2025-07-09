function moran_test

N = randi(100,20);
N(1:10,:) = 50;

%% 
%--------------------------------------------------------------------------
%                             Moran's I 
%--------------------------------------------------------------------------

% mean number of data points per bin 
meanHist            = mean(N(:));

% number of spatial units 
nUnits              = numel(N);

% weight matrix initialization
weightMatrix        = zeros(nUnits,nUnits);

sumSq = 0;
ijSumProd = 0;

for i = 1:nUnits
    
    jSumProd = 0; 
    
    for j = 1:nUnits
        
        [a,b] = ind2sub(size(N),i);
        [c,d] = ind2sub(size(N),j);
        
        if (a-c)^2+(b-d)^2 == 1
            weightMatrix(i,j) = 1;
            jSumProd = jSumProd + (N(i)- meanHist)*(N(j)- meanHist);
        end
        
    end
    
    ijSumProd = ijSumProd + jSumProd;
    sumSq = sumSq + (N(i) - meanHist)^2;

end
weightMatrixSum     = sum(weightMatrix(:));
moransI             = nUnits/weightMatrixSum*ijSumProd/sumSq

% calculate the z-score
%--------------------------------------------------------------------------

% expected value under null hypothesis
nullEV              = -1/(nUnits-1)

% variance 
sumWeights = 0;
for i = 1:nUnits
   for j = 1:nUnits
        sumWeights = sumWeights + (weightMatrix(i,j)+ weightMatrix(j,i))^2;
   end
end
S1 = sumWeights/2; 

S2 = 0;
for i = 1:nUnits
   ijSum = 0;
   jiSum = 0;
   for j = 1:nUnits
        ijSum = ijSum + weightMatrix(i,j);
        jiSum = jiSum + weightMatrix(j,i);
   end
   S2 = S2 + (ijSum + jiSum)^2; 
end

quartSum    = 0;
quadSum     = 0; 
for i = 1:nUnits
    diff        = N(i)-meanHist;
    quartSum    = quartSum + diff^4;
    quadSum     = quadSum + diff^2;
end
S3 = nUnits^(-1)*quartSum/(nUnits^(-1)*quadSum)^2;

S4 = (nUnits^2-3*nUnits+3)*S1 - nUnits*S2 + 3*weightMatrixSum^2;

S5 = (nUnits^2-nUnits)*S1 - 2*nUnits*S2 + 6*weightMatrixSum^2;

% variance of Moran's I 
varI = (nUnits*S4 - S3*S5)/((nUnits-1)*(nUnits-2)*(nUnits-3)*weightMatrixSum^2) 

%- nullEV^2
%varI2 = (nUnits*S4-S3*S1*(1-2*nUnits))/((nUnits-1)*(nUnits-2)*(nUnits-3)*weightMatrixSum^2)

zScore = (moransI - nullEV)/sqrt(varI)

end