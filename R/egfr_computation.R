#' A function to compute estimated glomular filtration rate eGFR
#'
#' Calculate eGFR
#'
#' @param df dataframe with individuals in rows and phenotypes (blood measures, age, sex, ect...) in columns
#' @param sex_colname column name for gender value
#' @param creatine_colname column name for creatine value
#' @param age_colname age column name for age value
#' @return dataframe with eGFR values as scores
#' @export
#'
#' @examples
#' calculate eGFR from dataset with phenotypes
#' new_df <- egfr_calculation(df,sex_colname="Sex",
#'           creatine_colname="Creatinine",
#'           age_colname="Age")
egfr_calculation <- function(df,sex_colname,creatine_colname,age_colname) {

  for (i in 1:dim(df)[1])
  {
    if (is.na(df[which(colnames(df)==creatine_colname)][i,]) == FALSE &&
        df[which(colnames(df)==sex_colname)][i,] == 2  &&
        df[which(colnames(df)==creatine_colname)][i,] <= 0.7) {

      df$eGFR[i] <-  166 * (df[which(colnames(df)==creatine_colname)][i,] / 0.7)^(-0.329) * (0.993)^(df[which(colnames(df)==age_colname)][i,])  }

    else if (is.na(df[which(colnames(df)==creatine_colname)][i,]) == FALSE &&
             df[which(colnames(df)==sex_colname)][i,] == 2  &&
             df[which(colnames(df)==creatine_colname)][i,] > 0.7) {

      df$eGFR[i] <- 166 * (df[which(colnames(df)==creatine_colname)][i,] / 0.7)^(-1.209) * (0.993)^(df[which(colnames(df)==age_colname)][i,])  }

    else if (is.na(df[which(colnames(df)==creatine_colname)][i,]) == FALSE &&
             df[which(colnames(df)==sex_colname)][i,] == 1  &&
             df[which(colnames(df)==creatine_colname)][i,] <= 0.9) {

      df$eGFR[i] <- 163 * (df[which(colnames(df)==creatine_colname)][i,] / 0.9)^(-0.411) * (0.993)^(df[which(colnames(df)==age_colname)][i,])}

    else if (is.na(df[which(colnames(df)==creatine_colname)][i,]) == FALSE &&
             df[which(colnames(df)==sex_colname)][i,] == 1 &&
             df[which(colnames(df)==creatine_colname)][i,] > 0.9) {

      df$eGFR[i] <- 163 * (df[which(colnames(df)==creatine_colname)][i,] / 0.9)^(-1.209) * (0.993)^(df[which(colnames(df)==age_colname)][i,])}

    else if (is.na(df[which(colnames(df)==creatine_colname)][i,]) == TRUE) {

      df$eGFR[i] <- NA}
  }
  return(df)
}
