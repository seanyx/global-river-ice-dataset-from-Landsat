
This GitHub repository includes Google Earth Engine (GEE) python/JavaScript scripts for calculating river ice and projecting future ice condition, R scripts for importing, cleaning, analyzing, and visualizing the data.

Files contained:

1. Calculate river ice extent from Landsat images
Script that calculate one river ice extent value per Landsat image for >800,000 Landsat images (1984–2018). The python script automates the calculation, exporting the data one month per GEE task.  

__GEE scripts/Download_landsat_river_ice_data_with_era_temp.py__: Calculate global river ice extent dataset.  
__GEE scripts/functions.py__: contains functions used in the "GEE scripts/Download_landsat_river_ice_data_with_era_temp.py" script.

2. Clean, analyze, and model river ice extent
The global river ice extent dataset was analyzed in R using the script "river ice analysis.Rmd".
Script "Calculae_river_ice_for_validation.js" contains code to calculate river ice extent near the gauge where in situ ice observations were available. The result was used to evaluate the river ice extent data.

3. Estimate future changes in river ice extent
The river ice extent model was implemented on GEE to calculate future river ice condition using CMIP5 climate model simulations. The predicted river ice extent data were then downloaded to local computer and analyzed using R scripts "river ice analysis.Rmd".  

__GEE scripts/Future annual river ice extent.js__: Calculate the annual river ice extent from 2008 to 2100.  
__GEE scripts/Future river ice extent comparing two periods.js__: Calculate the ice duration and ice extent difference between the periods 2009–2029 and 2080–2100.  
__GEE scripts/climate model SAT biases.js__: Calculate SAT biases for all models in the CMIP5 simulations, used as justification for choosing the three models used in the paper.
