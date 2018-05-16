#' A function to visualy inspect module robustness
#'
#' Testing robustness of modules by dropping half of samples
#'
#' @param metabolites_data metabolite dataframe
#' @param module_colors list of module colors generated
#' @param power1 power of scale free topology
#' @param main1 name of the plots
#' @param filename filename of the output plot
#'
#' @return plot
#' @export
#'
#' @examples set.seed(1)
#' pdf(file = "Robustness_dropping_samples.pdf", width = 12, height = 9)
#' par(mfrow=c(8,2), mar=c(2,2,2,1))
#' for (i in 1:8) {
#'  subsetindivs=sample(1:no.indivs,344,replace=F)
#'  robustness_drop_samples(metabolites_data = allsamples_norm_cl_trimmed[subsetindivs,],
#'                       module_colors = moduleColors,power1=6 ,main1="SCD Metabolites Network",
#'                       filename=paste("metabolites_clust_",i,".txt",sep=""))
#' }
#' dev.off()
robustness_drop_samples <- function(metabolites_data,module_colors, power1, main1,filename){
  if (dim(metabolites_data)[[2]] != length(module_colors)) {
    print("Dimension Error")
    }  else {
    adjacency = WGCNA::adjacency(metabolites_data, power = power1)
    TOM =  WGCNA::TOMsimilarity(adjacency)
    dissTOM = 1-TOM
    WGCNA::collectGarbage();
    geneTree = stats::hclust(stats::as.dist(dissTOM), method = "average")
    plot(geneTree,labels=F, xlab="", sub="", main=main1)
    plotColorUnderTree(geneTree,module_colors, main="SCD Metabolites Network Module Colors",rowLabels="")
    write.table(file=filename,data.frame(cbind(colnames(metabolites_data),module_colors)), quote = FALSE, row.names = FALSE, col.names = TRUE,sep="\t")
  }
}
