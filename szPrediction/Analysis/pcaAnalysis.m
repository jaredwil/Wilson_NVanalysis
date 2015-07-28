%Jared D. Wilson
%7/28/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PCA analysis for feature selection (suppliment script)
%trainFeats and testFeats must already be computed before running this

%loop through all features to see effect it has on pca analysis
Latent = zeros(size(trainFeats,2)-1,size(trainFeats,2)-1);
for p = 1:size(trainFeats,2)
    tmp = trainFeats;
    tmp(:,p) = [];
    [COEFF, SCORE, Latent(:,p)] = pca(tmp); 
    
end


%%
% manual PCA analysis
%trainFeats already normalized


S = cov(trainFeats);
[P,EigVals] = eig(S);
EigVals = diag(EigVals);
%flip Eigvals & P matrix
pcaCoef = fliplr(P);
Lat = flipud(EigVals);
Explained = cumsum(Lat)./sum(Lat);
precExp = (Lat./sum(Lat))*100;

coef90 = pcaCoef(:,Explained <= .90);

pcaFeats = trainFeats*coef90;

pcaTestFeats = testFeats*coef90;






