### Ensemble methods for survival analysis

### The functions are now also packed into "ensemblesurv" [package](https://github.com/dianashams/ensemblesurv).

#### !!! Warning: Note the project is at the development stage, so the code can change without any notification. The functions are not user-friendly to the extent they are intended to be once the project is completed.  Debugging and further nested cross-validation functionality is being done. 
#### Feel free to upload and use the code, however some bits may not be consistent (yet), use it responsibly in your project 

The project aims to provide researchers in the field of predictive modelling for health outcomes (e.g. predicting a disease incidence) with user-friendly functions to embed tree-based machine learning algorithms into the classical Cox proportionate hazards model (Cox PH). 

The algorithms can provide better predictive performance compared to the baseline CoxPH and quantify contribution of the non-linear relationships into the predictive preformance. The situation, when CoxPH have not been outperformed by the ensemble methods, can be used as a sign that strong non-linearities are unlikely to be present in the data, and classical CoxPH is as suited for analyzing such data.

Project page https://dianashams.github.io/ensemble-methods-for-survival-analysis/

#### If you use this project's code, please cite:

Shamsutdinova, D., Stamate, D., Roberts, A., & Stahl, D. (2022). Combining Cox Model and Tree-Based Algorithms to Boost Performance and Preserve Interpretability for Health Outcomes. In IFIP International Conference on Artificial Intelligence Applications and Innovations (pp. 170-181). Springer, Cham.

### How to use this code ###
This is not a package (yet), so please download the files, then source "R_simulating_survival_data.R", "Ensemble_Methods.R". You can then run "Examples_simulated_data.R" or "Example_GBSG2.R" which illustrates application of the methods, i.e. training and computing its apparent and internally-validated performance metrics.
**All functions assume the data is a data frame,  "time" and "event" columns define survival outcome; most functions require a list of column names that correspond to model's predictors**
