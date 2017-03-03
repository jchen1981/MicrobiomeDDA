# MicrobiomeDDA
1. mbzinb_0.2.tar.gz is a package implements a method of differential abundance tests between two groups.
To install the package, download the package, open a terminal, change path to the directory that contains the downloaded .gz file, and then type 'R CMD INSTALL mbzinb_0.2.tar.gz'. The usage of the package is documented as a standard R package.
2. zeroinfl.plus.github.R
It implements a fully generalized regression framework allowing the prevalance, abundance, and dispersion to depend on covariates. Existing packages does not allowing covariate-dependent dispersion, which could lead to either reduced power or inflated type I error if the heterogenety is not taken into account. It has similar syntax as 'zeroinf' from 'pscl' package.
3. zeroinfl.plus.example.R
Artifical simulations to illustrate the proposed method. 
4. zeroinfl.plus.realdata.R
Using a real microbiome data to illustrate the proposed method (GMPR + winsorization + omnibus test).
