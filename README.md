# raster_sampling
Compilation of resources and inter-comparison of raster sampling strategies

## Introduction and Motivation
Sampling a raster or ndarray at points is a common, basic need, but there is surprisingly little discussion or intercomparison of available options.

This is such an important topic, and the scientific literature is full of errors related to a lack of expertise when sampling (i.e. using nearest neighbor for large grid cells).  Most new users will go with rasterio sample function, without considering consequences.

There is value in general implementation of different approaches in a centralized library, potentially with different wrappers depending on input object (DataFrame, numpy arrays, xarray DataSet).  This could be installed as a dependency and used for different applications.

# Goals for sprint
* Continue compiling notes/resources on available approaches
* Discuss and prepare test data
* Design intercomparison function/notebook to document spread of values/results when using different methods
* Discuss potential integration with existing efforts (demquery, icepyx, rioxarray)

# Existing approaches
1. `rasterio` sample for 2D raster
    * @dshean started compiling raster sampling code in the IS2 Hackweek tutorial 
    * Appendix A of mega-notebook: https://github.com/ICESAT-2HackWeek/2020_ICESat-2_Hackweek_Tutorials/blob/master/05.Geospatial_Analysis/shean_ICESat-2_hackweek_tutorial_GeospatialAnalysis_rendered.ipynb.
1. `rasterstats` zonal statistics
    * See @dshean tutorial above for some references
    * `regionmask`: https://regionmask.readthedocs.io/en/stable/ (suggested by Emilio)
        * https://github.com/waterhackweek/waterdata
        * https://github.com/waterhackweek/waterdata/blob/master/mashup_waterbudget.ipynb
1. `scipy.ndimage.map_coordinates`
1. `scipy` interpolation options (for a local window of raster pixels)
1. `demquery` https://github.com/kylebarron/demquery (thanks to @scottyhq!)
1. `xarray` interpolation and sampling (http://xarray.pydata.org/en/stable/interpolation.html#advanced-interpolation)
    * Great discussion on these topics (and cross-referenced issues/PRs): https://github.com/pydata/xarray/issues/475
    * @scottyhq provided the following for xarray interpolation:
    * https://github.com/ICESAT-2HackWeek/pangeo-demo
        * Note that xarray is using numpy behind the scenes, and in this case scipy.interpolate http://xarray.pydata.org/en/stable/generated/xarray.DataArray.interp.html. 
        * Note that you can get rendered versions of hvplot if the repo is public and you use an nbviewer link (https://nbviewer.jupyter.org/github/ICESAT-2HackWeek/pangeo-demo/blob/master/atl06-sample-rendered.ipynb)
    * Official doc (could use some improvement - another potential task): http://xarray.pydata.org/en/stable/interpolation.html#advanced-interpolation
    * Thread from pangeo and IS2 Hackweek Slack: "you can change x and y from lists to DataArrays with a new dimension (‘z’) e.g.  x = xr.DataArray([x1, x2, x3], dims='z')"
1. `pangeo-pyinterp`: https://github.com/CNES/pangeo-pyinterp
1. Custom implementations in existing tools from @SmithB and @tsutterley
    * https://github.com/SmithB/pointCollection
    * https://github.com/tsutterley/spatial-interpolators
    * captoolkit?
1. Build bilinear equations in lsq framework, then solve: https://github.com/SmithB/LSsurf

1. https://geoviews.org/user_guide/Resampling_Grids.html (suggested by Friedrich)
1. https://github.com/OGGM/oggm

# Scaling Considerations
1. scipy ndinterp - create interpolator once, use multiple times
1. @tsutterley has implemention for MPI tile-based loading/sampling, using scipy splines:
    * https://github.com/tsutterley/read-ICESat-2/blob/master/scripts/MPI_DEM_ICESat2_ATL06.py
    * challenges around parallel HDF5, https://docs.h5py.org/en/stable/mpi.html

# Issues
* Handling of raster nodata/nan - many algorithms will propagate nans differently
* Performance for large point arrays, large rasters
* Input raster cell dimensions vs point sampling density
* Handling of different raster data types (e.g., only use nearest for classified rasters)
* xarray - great for ndarrays from GCMs, most efficient ways to chunk large point coordinate arrays, perform operations

# How to choose?
* Functions to assess input raster properties (e.g., classified vs. continuous, large gradients, spatial frequency content) and point properties
* Relative point density for points and raster
* Potentially add some logic to select best method based on these metrics

# Method Comparison
* For same input points and raster/ndarray, extract values from all samples
* No "truth", but can look at various statistics 
* ATL06 Already has DEM attributes in ATL06, ArcticDEM/REMA (reverts to GDEM) - details on how sampling is performed unknown (Jeff Lee)

# Candidate use cases
* Simple use case: starting with point coordinates in NumPy arrays, raster DEM as GeoTiff
* More advanced use case: GCM reanalysis data, xyzt point coordinates

# User Guide and Documentation
* Provide general recommendations based on typical scenarios
* Document case where different methods yield very different results

# Next steps
* Review resources listed above, compile notes
* Start creating notebooks with snippet for approach or specific problem, links to documentation
