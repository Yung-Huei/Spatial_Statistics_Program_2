# Nonstationary Spatial Modeling and Estimation Using a Divide-and-Conquer Approach
Collaborate with Hsin-Cheng Huang and Chun-Shu Chen
## Abstract
Regional heterogeneity in spatial covariance structure motivates nonstationary models that allow covariance behavior to vary across subregions. However, treating region-specific processes as independent can produce artifacts near subregion boundaries. We propose a globally valid nonstationary covariance framework that embeds region-specific stationary Gaussian components within a multivariate Mat ́ern process. The regional components retain their own variance, range, and smoothness parameters, while cross-component dependence and compactly supported transition weights allow them to blend coherently across space. A transition-width parameter controls the extent of this mixing, accommodating both relatively abrupt changes and gradual transitions in spatial dependence. The construction includes the stationary Mat ́ern model as a special case and converges to it as the regional covariance parameters become common, regardless of the modeling partition or transition width. We develop a scalable estimation and prediction strategy based on approximate composite likelihood and nested aggregation of local kriging predictors. Simulation studies examine covariance estimation, predictive performance, and robustness to the modeling partition. An application to a large land surface temperature dataset demonstrates the practical utility of the approach and its competitive predictive performance relative to established spatial methods.

## Code
1. simulation:  
   - **code1e03new.R**: The functions for our proposed method are used in simulation and real data analysis, and are stored as **code1e03new.RData**, which will be needed and loaded into R if the simulation code is executed.
   - **code1e03new_20260917.R**: Add two functions to generate identical Vehhia approximate sampling points for parallel and non-parallel computations in **code1e03new.R**.
   - **spatdiv.R**: This code is from Tzeng et al. (2024), who propose a decomposition method that partitions a spatial domain into distinct stationary components. For detailed information, please refer to the article:  
     [Tzeng, ShengLi, Bo-Yu Chen, and Hsin-Cheng Huang. "Assessing Spatial Stationarity and Segmenting Spatial Processes into Stationary Components." *Journal of Agricultural, Biological and Environmental Statistics* 29.2 (2024): 301–319.](https://doi.org/10.1007/s13253-023-00588-5)  Similarly, it is stored as **spatdiv.RData**, which will be needed and loaded into R if the simulation code is executed.
   - **Scenario I**: The simulation for scenario I includes six different estimation methods.
   - **Scenario II**: The simulation for scenario II includes six different estimation methods.
   - **Scenario III**: The simulation for scenario III includes six different estimation methods.
   - **Figure 6**: The simulation for LDK with differrent $M$ (Figure 6).
      
3. real data:
   - **AllSatelliteTemps.RData**: The dataset we use in Section 5.2: An Application of Daytime Land Surface Temperature Data.
   - **code1e03new.RData**: If the code is executed, it will be needed and loaded into R.
   - **spatdiv.RData**: If the code is executed, it will be needed and loaded into R.
   - **K5**: The data analysis for our proposed method is conducted under $K=5$.
   - **K6**: The data analysis for our proposed method is conducted under $K=6$.
   - **K7**: The data analysis for our proposed method is conducted under $K=7$.
   - **K8**: The data analysis for our proposed method is conducted under $K=8$.
   - **K9**: The data analysis for our proposed method is conducted under $K=9$.
   - **K10**: The data analysis for our proposed method is conducted under $K=10$.
   - **DeepKriging**: Data analysis of DeepKriging using five different random seeds.
