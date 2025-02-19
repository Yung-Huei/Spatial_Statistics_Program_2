# Fast Spatial Prediction for Nonstationary Processes Using a Divide-and-Conquer Approach 
Collaborate with Hsin-Cheng Huang
## Abstract
Analyzing spatial data over extensive domains is crucial for understanding complex environmental and geographical phenomena, yet it often reveals nonstationary features that challenge traditional statistical methods. These challenges include specifying appropriate covariance functions and implementing efficient kriging for large datasets. This paper introduces a novel methodology that leverages a linear combination of stationary processes, each with spatially varying weights, to model nonstationary processes. While traditional methods often treat subregional stationary processes as independent units, our approach integrates them into a comprehensive multivariate spatial process. This innovative integration reduces boundary discontinuities and enables a seamless transition between subregions without compromising computational efficiency. Our model allows the spatial covariance function to adapt smoothly or abruptly across different subregions through a tunable parameter. Notably, the model simplifies to a global stationary process when the covariance structures of all components are identical. For parameter estimation, we develop an efficient maximum composite likelihood approach that exploits the structured nature of local weight functions. Furthermore, we propose a novel double-kriging technique that ensures rapid spatial prediction with minimal efficiency loss through a divide-and-conquer strategy. Numerical experiments validate the effectiveness of our methodology in accurately estimating nonstationary spatial covariance functions and enhancing prediction accuracy.
## Code
1. simulation:  
   - **code1e03new.R**: The functions for our proposed method are used in simulation and real data analysis, and are stored as **code1e03new.RData**, which will be needed and loaded into R if the simulation code is executed.
   - **spatdiv.R**: This code is from Tzeng et al. (2024), who propose a decomposition method that partitions a spatial domain into distinct stationary components. For detailed information, please refer to the article:  
     [Tzeng, ShengLi, Bo-Yu Chen, and Hsin-Cheng Huang. "Assessing Spatial Stationarity and Segmenting Spatial Processes into Stationary Components." *Journal of Agricultural, Biological and Environmental Statistics* 29.2 (2024): 301–319.](https://doi.org/10.1007/s13253-023-00588-5)  Similarly, it is stored as **spatdiv.RData**, which will be needed and loaded into R if the simulation code is executed.
   - **Scenario I**: The simulation for scenario I includes six different estimation methods.
   - **Scenario II**: The simulation for scenario II includes six different estimation methods.
   - **Scenario III**: The simulation for scenario III includes six different estimation methods.
   - **Figure10**: The simulation for LDK with differrent $M$ (Figure 10).
      
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
