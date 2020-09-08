# raster_sampling
Compilation of resources and inter-comparison of raster sampling strategies

# Introduction and Motivation
Sampling a raster or ndarray at points is a common, basic need, but there is surprisingly little discussion or intercomparison of available options.

This is such an important topic, and I know that the literature is full of errors related to a lack of expertise when sampling (ie using nearest neighbor for large grid cells).  It might also be a good opportunity to demonstrate how these decisions can affect results.  I think most new users just go with rasterio sample function, without considering consequences.

I think there is real value in centralizing these in a library, potentially with different wrappers depending on input object (DataFrame, numpy arrays, xarray DataSet).  I’m envisioning a separate repo for this, so people can easily install and use for different applications.

# Goals for 9/8/20 sprint
* Continue compiling notes/resources on available approaches
* Prepare test data
* Design potential intercomparison function/notebook
* Discuss potential integration with existing efforts (demquery, icepyx, rioxarray)

# Existing approaches
1. `rasterio` sample for 2D raster
1. `rasterstats` zonal statistics
1. `scipy.ndimage.map_coordinates`
1. `scipy` interpolation options (for a local window of raster pixels)
1. `demquery` https://github.com/kylebarron/demquery (thanks to @scottyhq!)
1. `xarray` interpolation and sampling (http://xarray.pydata.org/en/stable/interpolation.html#advanced-interpolation)
1. Custom implementations in existing tools from @SmithB and @tsutterley
  * https://github.com/SmithB/pointCollection
  * https://github.com/tsutterley/spatial-interpolators
  * captoolkit?
1. regionmask: https://regionmask.readthedocs.io/en/stable/
 * Emilio https://github.com/waterhackweek/waterdata
 * https://github.com/waterhackweek/waterdata/blob/master/mashup_waterbudget.ipynb
1. https://geoviews.org/user_guide/Resampling_Grids.html
1. https://github.com/OGGM/oggm

I started compiling newer code in the IS2 Hackweek tutorial: Appendix A, https://github.com/ICESAT-2HackWeek/2020_ICESat-2_Hackweek_Tutorials/blob/master/05.Geospatial_Analysis/shean_ICESat-2_hackweek_tutorial_GeospatialAnalysis_rendered.ipynb.  Realizing that the scipy.ndimage.map_coordinates approach has an issue in the rendered notebook…

@scottyhq provided the following for xarray interpolation:
xarray interpolation example here https://github.com/ICESAT-2HackWeek/pangeo-demo. Note that xarray is using numpy behind the scenes, and in this case scipy.interpolate http://xarray.pydata.org/en/stable/generated/xarray.DataArray.interp.html. Note that you can get rendered versions of hvplot if the repo is public and you use an nbviewer link (https://nbviewer.jupyter.org/github/ICESAT-2HackWeek/pangeo-demo/blob/master/atl06-sample-rendered.ipynb)
http://xarray.pydata.org/en/stable/interpolation.html#advanced-interpolation
Thread from pangeo and IS2 Hackweek Slack:
* you can change x and y from lists to DataArrays with a new dimension (‘z’) e.g.  x = xr.DataArray([x1, x2, x3], dims='z')

# Scaling
1. scipy ndinterp - create interpolator once, use multiple times: https://github.com/SmithB/pointCollection
1. https://github.com/SmithB/LSsurf use lsq solver to do this, build equations
1. MPI tile-based loading/sampling, using scipy splines: https://github.com/tsutterley/read-ICESat-2/blob/master/scripts/MPI_DEM_ICESat2_ATL06.py
 * challenges around parallel HDF5, https://docs.h5py.org/en/stable/mpi.html

# Issues
* Handling of nodata gaps
* Performance for large point arrays, large rasters
* Input raster cell dimensions vs point sampling density
* Handling of different raster data types (e.g., only use nearest for classified rasters)
* Best solution for raster properties (e.g., large gradients, spatial frequency content)
* xarray - great for ndarrays from GCMs, most efficient ways to chunk large point coordinate arrays, perform operations

## ATL06
* Already has DEM attributes in ATL06, ArcticDEM/REMA (reverts to GDEM) - Jeff Lee implemented

# Candidate use cases
* Simple use case: starting with point coordinates in NumPy arrays, raster DEM as GeoTiff
* More advanced use case: GCM reanalysis data, xyzt point coordinates

# Next steps
* Review resources listed above, compile notes
* Start creating notebooks with snippet for approach or specific problem, links to documentation
