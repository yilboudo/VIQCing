#' A function to apply combat to your dataframe konwing the batch
#'
#' Apply batch effect correction using function ComBat
#'
#' @param df dataframe with individuals in rows and metabolites in columns
#' @param batch_df dataframe with information on batches (which indivudual belong to which batch)
#' @param batch_colname colname for actual batch
#' @return dataframe to which ComBat was applied
#' @export
#'
#' @examples
#' removing batch effect from dataset using combat
#' new_df <- batch_combat(allsamples, batch_data, individual_batches)
run_combat <- function(df, batch_df,batch_colnumb) {
  edata = df
  mod = model.matrix(~batch_df[,batch_colnumb], data=batch_df)
  mod0 = model.matrix(~1,data=batch_df)
  combat_edata = sva::ComBat(dat=base::t(edata), batch=batch_df[,batch_colnumb], mod=mod0, par.prior=TRUE, prior.plots=F)
  return(as.data.frame(t(combat_edata)))
}
