# RRP
Rainfall-Runoff-Phosphorus Model

This rainfall-runoff-phosphorus model is a python implementation of the RRP model by Hahn et al. (2013). The following modifications are made:

Forest HRU:
In the original model, no discharge was generated in forest HRUs. However, although fast flow discharge is unlikely to be generated in forested areas, we assume that slow flow discharge is likely to occur. Accordingly, in this model version, slow flow discharge is generated in forest HRUs and the necessary parameters are added.

Subsurface Drainage:
Several authors have highlighted the importance of subsurface drainage systems for phosphorus losses (for example Stamm et al. (1998)). Accordingly, we implemented an option to reflect this important pathway. If a map of drainage probability (the probability that a pixel is artificially drained) is provided, this information is used to adjust the likelihood of fast flow occurrence by modifying the topographic wetness index:

TWI_adj=TWI_ori*(P(D)+0.5)                       for P(D)>0.5

With

TWI_adj=Adjusted topographic wetness index
TWI_ori=Original topographic wetness index
P(D)=Drainage probabilty








Hahn, C.; Prasuhn, V.; Stamm, C.; Lazzarotto, P.; Evangelou, M. W. H.; Schulin, R. (2013): Prediction of dissolved reactive phosphorus losses from small agricultural catchments: calibration and validation of a parsimonious model, Hydrology and Earth System Sciences, 17(10), 3679-3693, doi:10.5194/hess-17-3679-2013.

Stamm, C.; Flühler, H.; Gächter, R.; Leuenberger, J.; Wunderli, H. (1998): Preferential transport of phosphorus in drained grassland soils, Journal of Environmental Quality, 27(3), 515-522, doi:10.2134/jeq1998.00472425002700030006x. 

