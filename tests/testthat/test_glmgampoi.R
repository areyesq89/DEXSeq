context("glmGamPoi")
test_that( "Funcion to create DEXSeqDataSet has the correct inputs", {

    data("pasillaDEXSeqDataSet", package="pasilla")

    dxd <- estimateSizeFactors( dxd )

    library(glmGamPoi)
    expect_message(
        dxd_glmgampoi <- estimateDispersions( dxd, fitType="glmGamPoi", BPPARAM=MulticoreParam(2), niter=2 ),
        "Parallelization has not been implemented for estimation" )
    
    dxd_glmgampoi <- testForDEU( dxd_glmgampoi, fitType="glmGamPoi" )
    expect_s4_class(dxd_glmgampoi, "DEXSeqDataSet")

    dxd_default <- estimateDispersions( dxd )
    dxd_default <- testForDEU( dxd_default )
    expect_s4_class(dxd_default, "DEXSeqDataSet")

    dxd_default_p2 <- estimateDispersions( dxd_default, BPPARAM=MulticoreParam(2) )
    dxd_default_p2 <- testForDEU( dxd_default_p2, BPPARAM=MulticoreParam(2) )
    expect_s4_class(dxd_default, "DEXSeqDataSet")

    res_default <- DEXSeqResults( dxd_default )
    res_default_p2 <- DEXSeqResults( dxd_default_p2 )

    expect_identical(res_default$pvalue, res_default_p2$pvalue)
    ## res_glmgampoi <- DEXSeqResults(dxd_glmgampoi, independentFiltering=FALSE)
    ## res_default <- DEXSeqResults(dxd_default, independentFiltering=FALSE)
    ## table( glmgampoi=res_glmgampoi$padj < 0.1, def=res_default$padj < 0.1 )

    dxd_glmgampoi_p2 <- testForDEU( dxd_glmgampoi, fitType="glmGamPoi", BPPARAM=MulticoreParam(3) )
    res_glmgampoi <- DEXSeqResults( dxd_glmgampoi )
    
    res_glmgampoi_p2 <- DEXSeqResults( dxd_glmgampoi_p2 )
    expect_identical( res_glmgampoi$pvalue, res_glmgampoi_p2$pvalue )

} )
