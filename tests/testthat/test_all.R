require(testthat)
test_that("InferMorphologyTree_solution",{	
	results <- InferMorphologyTree_solution(in.place=TRUE)
	expect_is(results$parsimony.tree, "phylo")
	expect_is(results$ml.tree, "phylo")
	expect_true(is.null(results$parsimony.tree$edge.length))
	expect_false(is.null(results$ml.tree$edge.length))
})

test_that("InferDNATreeWithBootstrappingAndPartitions_solution",{
	results <- InferDNATreeWithBootstrappingAndPartitions_solution(in.place=TRUE)
	expect_is(results$ml.tree, "phylo")
	expect_false(is.null(results$ml.with.bs.tree$node.label))
	expect_equal(length(results[[3]]), 100)
})

test_that("InferMorphologyTree_exercise",{	
	results <- InferMorphologyTree_exercise(in.place=TRUE, output.path="~/Desktop")
	expect_is(results$parsimony.tree, "phylo")
	expect_is(results$ml.tree, "phylo")
	expect_true(is.null(results$parsimony.tree$edge.length))
	expect_false(is.null(results$ml.tree$edge.length))
})

test_that("InferDNATreeWithBootstrappingAndPartitions_exercise",{
	results <- InferDNATreeWithBootstrappingAndPartitions_exercise(in.place=TRUE, output.path="~/Desktop")
	expect_is(results$ml.tree, "phylo")
	expect_false(is.null(results$ml.with.bs.tree$node.label))
	expect_equal(length(results[[3]]), 100)
})

