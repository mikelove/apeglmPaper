generateZINB <- function(deprob, groupprob){
	params <- newSplatParams()
	# LFC centered at 2
	params <- setParam(params, "de.facLoc", log(2) * 2)
	params <- setParam(params, "de.facScale", log(2) * 1)
	params <- setParam(params, "de.prob", deprob)
	# include drop-out 
	params <- setParam(params, "dropout.type", "experiment")
	params <- setParam(params, "dropout.mid", 2)
	# three groups, we will show apeglm improvements on the two small groups
	sim <- splatSimulate(params, group.prob=c(groupprob, groupprob, (1-2*groupprob)), method="groups", seed=1)
	table(sim$Group)
	# define the true LFC
	rowData(sim)$log2FC <- with(rowData(sim), log2(DEFacGroup2/DEFacGroup1))
	sim
}