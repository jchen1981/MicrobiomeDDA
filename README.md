# MicrobiomeDDA
1. mbzinb_0.2.tar.gz - a package implements differential distribution analysis between two groups based on zeroinflated negative binomial model. To install the package, download the package, open a terminal, change path to the directory that contains the downloaded .gz file, and then type 'R CMD INSTALL mbzinb_0.2.tar.gz'. The usage of the package is documented as a standard R package.
2. zeroinfl.plus.github.R - implements a fully generalized regression framework allowing the prevalance, abundance, and dispersion to depend on covariates. Existing packages do not allow covariate-dependent dispersion, which could lead to either reduced power or inflated type I error if the heterogenety is not taken into account. It has similar syntax as 'zeroinfl' from 'pscl' package.
3. zeroinfl.plus.daa.R - a wrapper function (ZISeq) to perform differential distribution analysis based on OTU table and meta dat.
4. zeroinfl.plus.example.R - Artificial simulations to illustrate the proposed method. 
5. zeroinfl.plus.realdata.R - Using a real microbiome data to illustrate the proposed method (GMPR + winsorization + omnibus test).
