szHorizon = 2;

train   = data.train;
test    = data.test;
% avgFeats = data.mean;
% stdFeats = data.std;

%isolate labels and feats in training data
trainLabels = train(:,1);
trainFeats = train(:,2:end);
testLabels = test(:,1);
testFeats = test(:,2:end);


%remove the interictal data so model is only trained on szhorizon for
%reression
% trainFeats(trainLabels > szHorizon*60*60,:) = [];
% trainLabels(trainLabels > szHorizon*60*60,:) = [];


%normalize features
avgFeats = mean(trainFeats,1);
stdFeats = std(trainFeats,[],1);
trainFeats = bsxfun(@rdivide, bsxfun(@minus,trainFeats,avgFeats), stdFeats);
% trainFeats(trainLabels > szHorizon*60*60,:) = [];
% trainLabels(trainLabels > szHorizon*60*60,:) = [];


tInt = lassoRes.int;
fLasso = lassoRes.coef;
% 
% testFeats(testLabels > szHorizon*60*60,:) = [];
% testLabels(testLabels > szHorizon*60*60,:) = [];
testFeats = bsxfun(@rdivide, bsxfun(@minus,testFeats,avgFeats), stdFeats);

utLasso = testFeats*fLasso + tInt;

lasso_coor = corr(testLabels,utLasso);
timeMin = linspace(0,(size(testLabels,1)-1)*30,size(testLabels,1))/60;

figure(60)
plot(timeMin,utLasso/60)
hold on;
plot(timeMin,testLabels/60)
hline(30,'r:');
ylabel('Time to Sz (min)')
xlabel('Time (min)')
xlim([0 max(timeMin)])
ylim([min(testLabels/60)-std(testLabels/60) , max(testLabels/60)+std(testLabels/60)])