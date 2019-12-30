This folder contains JavaScript code written to calculate river ice extent for gauge locations to evaluate river ice extent against in situ observations of river ice.

1. Scripts that read in and reformat the in situ data, export gauge locations then used to calculate ice extent from Landsat images on GEE
__01_import_NWS_ice_data.Rmd__  
__01_import_WSC_ice_data.Rmd__  
__01_merge_site_metadata.Rmd__  

2. Script calculate ice extent based on Landsat images
__calculate_river_ice_for_validation.js__: Calculate river ice extent for given gauge locations.  
__02_import_landsat_ice_merged.Rmd__: Import the Landsat-derive ice extent.  

3. Compare the gauge based river ice observations to the Landsat-derived river ice extent.
__03_validation_merged_sites.Rmd__: Carry out the comparison.  
