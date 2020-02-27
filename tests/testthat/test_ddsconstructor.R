context( "Create DEXSeqDataSet" )
test_that( "Funcion to create DEXSeqDataSet has the correct inputs", {
  data(pasillaDEXSeqDataSet, package="pasilla")
  mat <- featureCounts( dxd )
  pDat <- as.data.frame( colData(dxd) )
  pDat <- droplevels(pDat[pDat$exon == "this",])
  pDat$exon <- NULL
  dxdSub <- DEXSeqDataSet( mat, pDat, featureID = featureIDs( dxd ), groupID=geneIDs(dxd) )
  expect_equal( counts(dxd), counts(dxdSub) )
  expect_s4_class( dxdSub, "DEXSeqDataSet" )
  expect_error(
      DEXSeqDataSet( as.vector(mat), pDat, featureID = featureIDs( dxd ), groupID=geneIDs(dxd) ),
      "Unexpected input: the parameter 'countData' must be either a matrix or a data.frame" )
  expect_error(
      DEXSeqDataSet( mat, pDat, featureID = seq_along( featureIDs( dxd ) ), groupID=geneIDs(dxd) ),
      "Unexpected input: the parameter 'featureID' must be either a character or a factor" )
  expect_error(
      DEXSeqDataSet( mat, pDat, featureID = featureIDs( dxd ), groupID=seq_along( geneIDs(dxd) ) ),
      "Unexpected input: the parameter 'groupID' must be either a character or a factor" )
  expect_error(
      DEXSeqDataSet( mat, NULL, featureID = featureIDs( dxd ), groupID=geneIDs(dxd) ),
      "Unexpected input: the parameter 'sampleData' must be a data.frame" )
  expect_error(
      DEXSeqDataSet( mat, pDat, featureID = featureIDs( dxd ), groupID=geneIDs(dxd)[seq_len(nrow(mat)-1)] ),
      "Unexpected length of 'groupID' parameter, it must be the same as the number of rows of countData" )
  expect_error(
    DEXSeqDataSet( mat, pDat, featureID = featureIDs( dxd )[seq_len(nrow(mat)-1)], groupID=geneIDs(dxd) ),
    "Unexpected length of 'featureID' parameter, it must be the same as the number of rows of countData" )
  expect_error(
    DEXSeqDataSet( mat, pDat[seq_len(ncol(mat)-1),], featureID = featureIDs( dxd ), groupID=geneIDs(dxd) ),
    "Unexpected number of rows of the 'sampleData' parameter, it must be the same as the number of columns of countData" )
  expect_error(
    DEXSeqDataSet( mat, design=~1, pDat, featureID = featureIDs( dxd ), groupID=geneIDs(dxd) ),
    "The design formula does not specify an interaction contrast with the variable 'exon'" )
  expect_error(
    DEXSeqDataSet( mat, pDat[,1,drop=FALSE], featureID = featureIDs( dxd ), groupID=geneIDs(dxd) ),
    "The variables 'condition', present in the design formula must be columns of 'sampleData'" )
} )
