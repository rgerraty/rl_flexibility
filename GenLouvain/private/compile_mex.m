%	Script to compile the mex-files used by GENLOUVAIN and GENLOUVAINRAND

mex -largeArrayDims modchange_y.c

mex -largeArrayDims tidyconfig_c.c

mex -largeArrayDims modf_bipord.c