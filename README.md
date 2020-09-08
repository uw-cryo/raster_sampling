# raster_sampling
Compilation of resources and inter-comparison of raster sampling strategies

# Introduction and Motivation
Sampling a raster or ndarray at points is a common, basic need, but there is surprisingly little discussion or intercomparison of available options.

This is such an important topic, and I know that the literature is full of errors related to a lack of expertise when sampling (ie using nearest neighbor for large grid cells).  It might also be a good opportunity to demonstrate how these decisions can affect results.  I think most new users just go with rasterio sample function, without considering consequences.

I think there is real value in centralizing these in a library, potentially with different wrappers depending on input object (DataFrame, numpy arrays, xarray DataSet).  I’m envisioning a separate repo for this, so people can easily install and use for different applications.

# Existing approaches
1. `rasterio` sample for 2D raster
1. `rasterstats` zonal statistics
1. `scipy.ndimage.map_coordinates`
1. `scipy` interpolation options (for a local window of raster pixels)
1. `demquery` https://github.com/kylebarron/demquery (thanks to @scottyhq!)
1. `xarray` sampling
1. Custom implementations in existing tools from @SmithB and @tsutterley
  * https://github.com/SmithB/pointCollection
  * https://github.com/tsutterley/spatial-interpolators
  * captoolkit?

I started compiling newer code in the IS2 Hackweek tutorial: Appendix A, https://github.com/ICESAT-2HackWeek/2020_ICESat-2_Hackweek_Tutorials/blob/master/05.Geospatial_Analysis/shean_ICESat-2_hackweek_tutorial_GeospatialAnalysis_rendered.ipynb.  Realizing that the scipy.ndimage.map_coordinates approach has an issue in the rendered notebook…

# Issues
* Handling of nodata gaps
* Performance
* Input raster cell dimensions vs point sampling density
* Handling of different raster data types (e.g., only use nearest for classified rasters)

# Goals
* Start compiling existing resources
* Prepare test data
* Design intercomparison 
