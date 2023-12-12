# latentMAR
Package accompanying LMAR outcome paper

The purpose of this package is to share the code that was used to obtain the results for the data example (Experience Corps study) in the paper, specifically the results in Tables 4 (sample description) and 5 (effect estimates under different assumptions).

The package includes the conventional R folder that is the main part of the R-package structure, and an "analysis" folder that is specific to the application.

## The R folder

The R folder contains mostly functions that can be useful generally, including functions that compute mixture weights based on specific missingness assumptions, stratum-specific outcome means under control based on those mixture weights and principal identification assumptions, and effects based on those outcome means, etc. 

There is a single function in this folder that is specific to the application, which computes the nuisance functions based on models specified for the application. This function could alternatively be put in the "analysis" folder; putting it in the R folder seemed to reduce friction, but maybe it's a lazy choice on our part. Anyway, this function is there for transparency, but it is unlikely anyone will be interested in looking at it.

## The analysis folder

This folder shows a single subfolder "EC". (The reason for this additional level of hierarchy is that we have other examples in other subfolders that are not shared on github.)

The code for the example resides in the folder "script" within "EC", with five files for data cleaning, sample description (ie Table 4), nuisance estimation, point estimates and bootstrapping, which are all self-obvious. We'd like to highlight only two pieces:

1. To see how the point estimates are computed, you can look in the file EC.3-point_estimates.R. It shows three steps after estimating the nuisance functions: 

- computing the mixture weights based on the specific missingness assumption of choice -- using the latentMAR::mixture_weights() function
- computing the stratum-specific conditional outcome means under control based on the mixture weights and the principal identification assumption of choice -- using the latentMAR::mus_under_control() function
- computing the effect estimates (via the plugin estimator) based on the above outcome means -- using the latentMAR::plugin_estimate() function

2. The whole procedure above is put in the latentMAR::effects_all_assumptions() function, and this function is used in the bootstrap. You can see how the bootstrap is done in file EC.4-bootstrap.R. This file obtains both the point estimates and confidence intervals in Table 5 in the paper.
