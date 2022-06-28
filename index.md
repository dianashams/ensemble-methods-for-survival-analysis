## Ensemble methods for survival analysis project

Survival analysis is a brunch of methods that analyse time-to-event data. That is, it estimates the time it takes for an event to happen in a given population, and how  it dependis on various parameters. The applications could be estimating what proportion of the population will fail by certain point, or how the rate of event changes with time. Population maybe a population of bacteria, people, or group of financial companies, or a company client base; an event can be death, or disease onset or recovery from a disease, bankruptcy, or moving to a competitor. 

This project contains several functions written in R language that mix classical Cox regression for survival data and decision tree-based machine learning algorithms (survival decision tree and survival random forest). We mixed classical regression and ML methods while thinking about the health outcomes (e.g. a disease onset), and to develop a -friendly way of incorportating machine learning for health researchers. In this, we aimed to preserve algorithm interpretability and get an insight into the data structure, as much as we targeted a better predictive performance compared to a baseline and widely used Cox Proportionate Hazards model.

The link to the code is [here](https://github.com/dianashams/ensemble-methods-for-survival-analysis)

#### Methods  
The methods' description  and initial results were presented at the [talk](https://youtu.be/1Z8C0pAi_Cs) at the [NIHR Maudsley Biomedical Research Center Prediction Modelling Group](https://www.maudsleybrc.nihr.ac.uk/facilities/prediction-modelling-presentations/) meeting (September, 2021). 

More details are in the [AIAI 2022 conference paper](https://link.springer.com/chapter/10.1007/978-3-031-08337-2_15)[1].

#### Links and references
[1] Shamsutdinova, D., Stamate, D., Roberts, A., & Stahl, D. (2022). Combining Cox Model and Tree-Based Algorithms to Boost Performance and Preserve Interpretability for Health Outcomes. In IFIP International Conference on Artificial Intelligence Applications and Innovations (pp. 170-181). Springer, Cham.

[2] Amunategui, M.: Data Exploration & Machine Learning, Hands-on. https://amunategui.github.io/survival-ensembles/index.html

[3] Marmerola, G.D.: Calibration of probabilities for tree-based models | Guilherme’s Blog. https://gdmarmerola.github.io/probability-calibration/

[4] Su, X., Tsai, C.-L.: Tree-augmented Cox proportional hazards models. Biostatistics 6, 486–499 (2005)

[5] Breiman, L., Friedman, J.H., Olshen, R.A., Stone, C.J.: Classification and Regression Trees. Wadsworth & Brooks/Cole Advanced Books & Software, Monterey, CA (1984)

[6] Shimokawa, A., Kawasaki, Y., Miyaoka, E.: Comparison of splitting methods on survival tree. Int. J. Biostat. 11, 175–188 (2015)

[7] Ishwaran, H., Lauer, M.S., Blackstone, E.H., Lu, M.: randomForestSRC: Random Survival Forests Vignette (2021)

[8] Royston, P., Altman, D.G.: Regression using fractional polynomials of continuous covariates: parsimonious parametric modelling. Appl. Stat. 43, 429 (1994)

[9] Therneau, T., Atkinson, E.: An introduction to recursive partitioning using the RPART routines (2019)

[10] Ishwaran, H.: Variable importance in binary regression trees and forests. Electron. J. Stat. 1, 519–537 (2007)

#### Support or Contact
diana.shamsutdinova@kcl.ac.uk

