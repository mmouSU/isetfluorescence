Image Systems Engineering Toolbox for Fluorescence (ISETfluorescence)

The ISETfluorescence toolbox provides a set of functions to simulate scenes with fluorescent components. Here, we extend the open-source and freely available tools in Image Systems Engineering Toolbox for cameras (ISETCam).  The  toolbox also implements several reflectance and fluorescence separation and estimation algorithms.

If you use this code in your work, please cite

Henryk Blasinski, Joyce Farrell and Brian Wandell - "Simultaneous reflectance
and fluorescence spectra estimation" 2016 
(https://arxiv.org/abs/1605.04243)


0. License
----------

The code is provided as is. You are free to use and modify the code in 
non-commercial and research applications. 
If you are interested in commercial applications, please contact us as the 
method and apparatus is a subject of a US Patent 20,160,116,410.


1. Dependencies (required)
--------------------------

a. Image Systems Engineering Toolbox (ISETCam) (https://github.com/ISET/ISETCam).
b. cvx Convex Optimization toolbox for Matlab (www.cvxr.com).


2. Dependencies (optional)
--------------------------

The following code is used by some data plotting scripts and are not required
to simulate or analyze reflectance and fluorescence.

a. Barwitherr function from Matlab File exchange. 
   (http://www.mathworks.com/matlabcentral/fileexchange/30639-barwitherr-errors-varargin).
b. Code accompanying the Computational Colour Science using Matlab from
   Matlab File exchange.
   (http://www.mathworks.com/matlabcentral/fileexchange/40640-computational-colour-science-using-matlab-2e)


3. Sample data
--------------

The sample inputs to our algorithms as well as the results of the analyses 
can be downloaded from Stanford's Digital Repository: 
http://purl.stanford.edu/xc528jd5098
The zip file downloaded from the above side contains two directories: 'data'
with raw input data from real captures and simulations and 'results' storing
algorithm outputs. Please place these two folders directly in the fiToolbox
root folder.




