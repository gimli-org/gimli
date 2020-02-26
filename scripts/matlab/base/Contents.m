% GIMLab base functions
% arms           - absolute root mean squared misfit
% b2r            - Blue-white-red colormap for difference plots
% cbar           - draw colorbar
% chi2           - calculate weighted L2 residual
% curvature      - curvature of a parametric function
% destrip        - strips string from comments (#...)
% epsprint       - export figure into eps (and pdf) file
% exportfig      - general export function for figures
% exportpng      - export figure as png image
% fcmcluster     - Fuzzy c-means cluster analysis
% gaulag         - Gauss-Laguerre quadreature points/weights
% gauleg         - Gauss-Legendre quadreature points/weights
% getcindex      - get color indices for graphics patching
% getgrfile      - UI function for exporting graphics file (eps/png)
% gps2xyz        - converts GPS coordinates (HHMMSS) to xyz values
% harmfit        - least squares fit of time series by harmonic functions
% int2hum        - converts number in human readable (12.3M) string
% interperc      - computes interpercentiles from vector histogram
% inter2dpoint   - interpolate 2d to points by delaunay interpolation
% minmax         - yields minimum/maximum value [min(x) max(x)]
% num2strcell    - converts numerical array into string cell (for labels)
% readconfigfile - reads confic file of (token=value) into structure
% rrms           - relative root mean squared misfit
% rndig          - round value(s) to number of counting digits
% savematbin     - save matrix to binary file
% snaptoline     - snap 2d points onto straight line
% turnmatrix     - compute turn matrix out of three angles
% writeconfigfile - writes data structure into config file (token=value)
% writepng        - write figure into png image
% writesurfergrid - write matrix into surfer (GRD) file
% zeil2vec	- converts line into vector of numerical values

%deprecated functions for backward compatibility
% antiindex - overall index out of I,J,K position in a 3d grid
% cbar2.m   - alternative version of cbar
% finex1,finex2,finex3 - specialized export functions
% finexbw.m - export function in grayscale
% kruemm    - alternative version of curvature
% tomkurv   - alternative version of curvature
% message3  - alternative version of message 
% rms2.m    - alternative version of rms
