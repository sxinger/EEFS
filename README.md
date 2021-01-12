# Ensemble Embedded Feature Selection Framework (EEFS)

This is an open repository sharing codes of the stability-predicatbility feature selection framework, which we initially developed for DKD risk factor discovery. However, this framework can be adoptted for other feature selection tasks. 

The figure above also shows the 4 steps of the feature selection framework: 
1) feature ranking – where we apply selective methods of embedded feature selection on bootstrap samples (provides better accuracy than traditional feature selection) for ensemble classifiers and rank the features in accordance with their “importance” or their discriminatory power in association with the outcomes. Shapely Additive Explanation (SHAP) values will be used to rank features based on their average association with outcomes.  This not only provides better interpretability but better consistency over other ranking methods.
- a) regularized regression which assumes an additive structure and benefits linear relationships (LASSO); 
- b) tree-based gradient boosting ML which assumes hierarchical structure and benefits nonlinear relationships (GBM),
- c) neural network which incorporates a structure handling more complex correlations (DNN). All three methods have demonstrated success in a variety of biomedical studies.


2) feature ensemble – where we aggregate the rankings with different ensemble techniques to mediate feature ranking variations against sampling perturbations. 6 ensemble strategies are currently supported: 
- a) Mean Aggregation
- b) Stability Aggregation
- c) Exponential Aggregation
- d) Weighted Mean Aggregation
- e) Weighted Stability Aggregation
- f) Weighted Exponential Aggregation

3) feature sizing – for each combination of feature selection method and ensemble technique, or ensemble-feature-selection model, we will conduct a golden-section search to estimate the minimal feature size sufficient for that particular classifier to achieve a close-to-optimal accuracy based on 5-fold cross-validation. We adopted an iterative golden-section search procedure for approximating a minimal feature size to achieve accuracy significantly close to the optimum using the Delong test. 

4) evaluation – we will identify the optimal ensemble-feature-selection model based on feature predictability and stability. The feature selection results will identify a minimal set of interventional features, simple or amalgamated or bundling, that are strongly associated with pain improvement. The associations, usually measured by average absolute changes in the likelihood of improving outcomes, are used to rank the features with a sign (positive or negative) to identify the direction of the association. Below we describe each step in more detail. 
- a) Predictability: will be evaluated by areas under receiver operating and precision-recall curves, and as by calibration on validation sets. 
- b) Stability: will be evaluated by the Kuncheva Index and a Weighted Consistency Index (WCI) for fixed feature sizes and a Relative WCI for flexible feature sizes adapted to different feature selection algorithms.


Please cite our paper if you are using our package: 

•	X. Song, L.R. Waitman, Y. Hu, ASL. Yu, D. Robbins, M. Liu. Robust Clinical Market Identification for Diabetic Kidney Disease with Ensemble Feature Selection. Journal of the American Medical Informatics Association. 2019. 26(3):242-253.
