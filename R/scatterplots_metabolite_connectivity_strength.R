#' A function plot module connection strength
#'
#' Plotting Module membership (x-axis) and metabolite significance (y-axis)
#'
#' @param module_color module color
#' @param trait name of trait
#' @param dataset dataframe with module membership and metabolite significance
#' @param alt_color module color different from the original module name
#'
#' @return plot
#' @export
#'
#' @examples scatterplots_metabo(module_color="darkorange",trait="RBC.beta",
#' dataset=metabo_bloodtrait_cor_final,alt_color="")
#' scatterplots_metabo(module_color="lightcyan",trait="eGFR.beta",
#' dataset="metabo_bloodtrait_cor_final",module_membership_cor="module_membership_cor",alt_color="blue")
scatterplots_metabo = function (module_color, trait,metabolite_significance,module_membership_cor,alt_color) {
  module = module_color
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  if (alt_color == ""){
    WGCNA::verboseScatterplot(abs(module_membership_cor[moduleGenes, column]),
                       abs(data.frame(metabolite_significance[,which(names(metabolite_significance)==trait)])[moduleGenes, 1]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = paste("Metabolite significance for", gsub(".beta","",trait), sep=" "),
                       main = paste("Module membership vs. metabo significance\n"),
                       cex.main =.9, cex.lab = .9, cex.axis = 1.2, col = module_color,pch=19)
  } else {
    WGCNA::verboseScatterplot(abs(module_membership_cor[moduleGenes, column]),
                       abs(data.frame(metabolite_significance[,which(names(metabolite_significance)==trait)])[moduleGenes, 1]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = paste("Metabolite significance for", gsub(".beta","",trait), sep=" "),
                       main = paste("Module membership vs. metabo significance\n"),
                       cex.main =.9, cex.lab = .9, cex.axis = 1.2, col = alt_color,pch=19)

  }
}
