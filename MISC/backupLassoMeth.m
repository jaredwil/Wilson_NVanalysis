%split up for a simple 2 fold CV
% tr = length(retrainLabels)*0.7;
% evalFeats = retrainFeats(tr+1:end,:);
% evalLabels = retrainLabels(tr+1:end);
% retrainFeats = retrainFeats(1:tr,:);
% retrainLabels = retrainLabels(1:tr); 
% 
% tic;
% disp(['Training Lasso Model on Patient: ' pt{i} ' (NO INTERICTAL)'])
% opts = statset('UseParallel',true);  %Do dis shit in parallel!!!  use number of availble workers for CV
% [fLasso, tINFO] = lasso(retrainFeats, retrainLabels,'CV',5,'Options',opts);
% lassoTime(i) = toc;
% %BOOM!!!
% disp(['DONE Training Lasso Model on Patient: ' pt{i}])
% dzT = tINFO.DF; 
% tInt = tINFO.Intercept;
% % xAlpha = xINFO.Alpha;
% tInt2 = repmat(tInt, size(retrainFeats,1),1);
% utLasso = retrainFeats*fLasso + tInt2;
% tTrain = repmat(retrainLabels,1,100);
% tLasso_coor = corr(tTrain,utLasso);
% tLasso_coor = tLasso_coor(1,:);
% 
% figure(11)
% plot(dzT,tLasso_coor*100,'r.');
% xlabel('Number of Non-Zero Feature Weights in Lasso Model')
% ylabel('Resulting Test Correlation (%)')
% title('TRAIN COULD BE MISLEADING')
% grid on;
% hold on;
% 
% tInt2 = repmat(tInt, size(evalFeats,1),1);
% utLasso = evalFeats*fLasso + tInt2;
% tTrain = repmat(evalLabels,1,100);
% tLasso_coor = corr(tTrain,utLasso);
% tLasso_coor = tLasso_coor(1,:);
% 
% minIdx = tINFO.IndexMinMSE;
% seIdx  = tINFO.Index1SE;
% corrIdx = tLasso_coor == max(tLasso_coor);
% 
% figure(12)
% plot(dzT,tLasso_coor*100,'r.');
% xlabel('Number of Non-Zero Feature Weights in Lasso Model')
% ylabel('Resulting Test Correlation (%)')
% title('EVALUATION')
% grid on;
% hold on;
% h = plot(dzT(minIdx),tLasso_coor(minIdx)*100,'gs','LineWidth',4);
% h = vline(dzT(minIdx),'g:');
% h = plot(dzT(seIdx),tLasso_coor(seIdx)*100,'bs','LineWidth',4);
% h = vline(dzT(seIdx),'b:');
% 
% figure(14)
% h2 = lassoPlot(fLasso,tINFO,'PlotType','CV');
% % Use a log scale for MSE to see small MSE values better
% % set(gca,'YScale','log');
% % plotName = ['H:\jaredwil\Lasso Results\matLasso_history\' pt{i} '_mseRes_CCS'];
% 
% bestLasso_corr = fLasso(:,corrIdx);
% bestLasso_1SE = fLasso(:,seIdx);
% bestLasso_Min = fLasso(:,minIdx);
% 
% bestInt_corr = tInt(corrIdx); 
% bestInt_1SE = tInt(seIdx); 
% bestInt_Min = tInt(minIdx); 
% 
% numFeats_corr = dzT(corrIdx); 
% numFeats_1SE = dzT(seIdx); 
% numFeats_Min = dzT(minIdx);