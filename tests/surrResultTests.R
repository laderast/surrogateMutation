require(testthat)

setwd("../data/")
load("surrogateResult.RData")
load("intome.RData")

testthat("node properties work for getSubnetworksFromSample",{
  test <- getSubnetworkForSamples("BRCA1(672)", intome, surrResult = surrogateResult)
  expect_equal(length(nodes(test)), 7)
  nodeCounts <- nodeData(test, attr="counts")
  expect_equal(nodeCounts$MYC, 3)
  #expect_equal()
}) 

testthat("buildSurrogateTable is Correct", {
  test <- buildSurrogateTable(surrogateResult)
  expect_equal(class(test$Gene), "factor")
  expect_equal(class(test$ID), "factor")
  expect_equal(class(test$degree), "numeric")
})

