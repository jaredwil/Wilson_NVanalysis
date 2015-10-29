# Wilson_NVanalysis
**Pipeline and Analysis Tools used for Master's Thesis Project**<br />
**Jared Wilson**<br />
**2015**<br />


##TIME 2 SEIZURE PREDICTION MODEL
This project used data collected by Neuro Vista as part of the first implanted seziure likelihood study conducted 2010-2011. The dataset was stored on the cloud (IEEG.org) for this project, and was accessed through the used of a developed MATLAB tookbox avialable on the site. Many supporting functions and libraries were used in this project for ML applications, interacting with the cloud, and parallel proecessing. Broad overviews of these areas are provided in my Thesis (link provided). The majority of the paper describes the method and algorithm used. This readme will explain the DEMO files provided for a new user that has access to this dataste. This project was published in 2015 as part of my Master's Thesis in Electrical Engineering at the University of Pennsylvania. 

###Model Description -- LASSO
Thesis link is below to provide an overview of the methods and model used in this project.

LINK TO PAPER: 

###DEMO
two demo files are provided and will work given you have the IEEG Matlab Toolbox and access to the dataset. <br />

**DEMOv1** <br />
This is a very simple demo that loads in the lasso model for a selected patient the makes predictions for all two hour seizure horizons (only for Clinically Confirmed Seizures). The data are in .mat file contained on FOUIER with a path included in the function. Two plots are generated from the results. One is a plot that shows the correlation between true and predicted lables for each seizures throughout the recording, each line corrosponding to a seizure. The X-axis is time of recording (days) and the Y-axis is correlation for that seizure. The second generated plot is a visulization of the model prediction compared to the time to seizure label. The seizures plotted in this figure and test seizures that have a correlation greater than a value set by the user (default threshold = 55%).

**DEMO_interictal** <br />
This function is used to test models on an extended seizure horizon. The function interacts with the IEEG portal and calculates windows manually unlike the DEMOv1. This makes the script take MUCH MUCH longer to run as it must download and process hours of data depending on the extended seizure horizon selected. The plots generated in this function can also look a little weird if trying to look at too many seizures at once. Once again a corr threshold is used to determine the seizures you would like to analyze.   <br />

This demo file is a little more abstract and messy, read through it if you desire.

###Want to get more feature vectors to test model???
A pipeline is established to do this that uses the code found int szPrediction\Pipeline. A brief overiew of each function is provided below (additional comments and documentation are provided in functions): <br />
<br />
**FeatExt.m**  <br />
function [ feats ] = FeatExt( y , fs)  <br />
Input to this function is a time window y and sampling freq (fs). Output is array of features extracted in function.  <br />

**getSzFeats.m** <br />
[trainFeats, predIdx] = getSzFeats(ptSession, szStartT, szEndT, szHorizon, winLen, winDisp) <br />
This function takes in the current patient iEEG session and extracts windowed features based on winLen and winDisp entered by the user. Features are only extracted from times that are within the determined szHorizon (entered by user) and the seizure start times (szStartT). <br />

**getSzTime.m** <br />
[ startT, endT ] = getSzTime(ptSession)  <br />
Get ALL seizure start times from the pt of interest. The seizure start times are extracted from the annotation layers. All layers containing "seizure" are extracted.

**szPred_train.m**/**szPred_trainNOint.m**  -- variation that doesn't retrieve interictal data <br />
[trainFeats, trainLabels ] = szPred_test(pt, usernm, pswdBin, trPct, winLen, winDisp, szHorizon  )  <br />
szPred_train takes the current pt and retrieves the trainFeats and trainLabels based on the inputs provided by the user. The szHorizon defines the 'preictal' phase where features will be extracted from all lables that are larger than the defined szHorizon are interictal feature vectors. (NOint version does not get interictal data)  <br />

**szPred_test.m**/**szPred_testNOint.m**  -- variation that doesn't retrieve interictal data <br />
[ testFeats, testLabels ] = szPred_test(pt, usernm, pswdBin, trPct, winLen, winDisp, szHorizon  )
szPred_test takes the current pt and retrieves the testFeats and testLabels based on the inputs provided by the user. The szHorizon defines the 'preictal' phase where features will be extracted from all lables that are larger than the defined szHorizon are interictal feature vectors. (NOint version does not get interictal data)  <br />

**test_interPar.m**/**train_interPar.m** --does the same thing just gets interictal data from either side of train/test divide <br />
[ interFeats ] = test_interPar( ptSession, szpredIdx, winLen, winDisp, Nlim ) <br />
This function is used to extract interictal data from the NV patient data to be trained on for interictal classification. The times that included sz and szpred Horizon are inputs to be skipped when searching for valid interictal data.



