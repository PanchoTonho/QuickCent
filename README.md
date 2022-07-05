# QuickCent

Support code for "QuickCent: a fast and frugal heuristic for centrality estimation on networks", Plana & Pérez (2018), research article published on 2018 IEEE / WIC / ACM International Conference on Web Intelligence, which is part of the (in progress) PhD dissertation "Computational models for network societies" by Francisco Plana. The paper may be revised in the following link.

https://ieeexplore.ieee.org/abstract/document/8609599

All required code is on QCent_script.R. Next, a brief explanation of the main methods:

MAE_dist_train_sizes() is the method that produced data required to produce Fig 1 (Robustness of estimates). This dataset is recorded on resultados_simulaciones8.csv, which may be plotted by plot_MAE_dist_train_sizes().

produce_comparison_methods() is the method that produced data required to produce Fig 2 and 3 (Comparison with other methods).  These datasets are recorded on resultados_metodos8_0.1.csv (10 % as training set) and resultados_metodos8_1.csv (100 % as training set), and can be plotted by plot_comparison_methods(). These files also record time measurements on Tables III and IV, which may be reproduced with get_statistics_Table_III_IV(). Original paper has some minor round errors in the elapsed time measurements, which do not alter at all the conclusions.

assumption_validation() is the method that produced the dataframe storing the data used to validate the assumptions of the heuristic (Setting and assumptions verification). This dataset is stored on valid_assump.csv, and the statistics from Tables I and II may be reproduced by get_statistics_Table_I_II().

All the code was written on R, where the dependencies are the following packages:

RWeka
igraph
Matrix
ggplot2
poweRlaw.


