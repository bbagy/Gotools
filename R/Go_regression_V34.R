
#' Perform Regression Analysis on Given Data
#'
#' @param data The data frame containing the variables for analysis.
#' @param project The project name for output file naming.
#' @param cate.outs Categorical outcome variables (optional).
#' @param con.outs Continuous outcome variables (optional).
#' @param cate.vars Categorical predictor variables (optional).
#' @param con.vars Continuous predictor variables (optional).
#' @param mul.vars Flag for multivariate analysis (TRUE/FALSE).
#' @param interaction Interaction terms for the model (optional).
#' @param randomEff Random effects for mixed models (optional).
#' @param orders Order of factor levels for categorical variables.
#' @param pvalue P-value threshold for significance.
#' @param name Additional name identifier for output files (optional).
#'
#' @details
#' This function performs various types of regression analysis (linear, logistic, mixed models) based on the specified parameters. It supports univariate, multivariate, and mixed-effects modeling. The function adjusts the p-values using the False Discovery Rate (FDR) method.
#'
#' @return
#' The function generates a CSV file for each outcome variable containing the model coefficients, standard errors, p-values, adjusted p-values, and confidence intervals. It also saves RDS files for multivariate models.
#'
#' @examples
#' Go_regression(data = my_data,
#'               project = "MyProject",
#'               cate.outs = c("Outcome1", "Outcome2"),
#'               con.outs = "Outcome3",
#'               cate.vars = c("Predictor1", "Predictor2"),
#'               con.vars = "Predictor3",
#'               mul.vars = "Predictor1 + Predictor2 + (1 | Subject)",
#'               randomEff = "Subject",
#'               orders = c("Level1", "Level2"),
#'               pvalue = 0.05,
#'               name = "Analysis1")
#'
#' @export

Go_regression <- function(data, project,
                          cate.outs=NULL,
                          con.outs=NULL,
                          cate.vars=NULL,
                          con.vars=NULL,
                          mul.vars=FALSE,
                          randomEff=NULL,
                          orders, pvalue=0.05, name=NULL){
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_table <- file.path(sprintf("%s_%s/table/regression",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out_table)) dir.create(out_table)

  # data control
  # fix outcome column types
  if(!is.null(cate.outs)){
    for (cate.out in cate.outs) {
      data[,cate.out] <- factor(data[,cate.out], levels = intersect(orders, data[,cate.out]))
    }
    outcomes <- c(cate.outs)
  }

  if(!is.null(con.outs)){
    for (con.out in con.outs) {
      data[,con.out] <- as.character(data[[con.out]])
      data[,con.out][data[,con.out]==""] <- NA  # "NA" 대신 NA로 처리
      data[,con.out] <- as.numeric(data[[con.out]])
    }
    outcomes <- c(con.outs)
  }

  if(!is.null(cate.outs) & !is.null(con.outs)){
    outcomes <- unique(c(cate.outs, con.outs))
  }

  # fix variation column types
  if(!is.null(cate.vars)){
    for (cate in cate.vars) {
      data[,cate] <- factor(data[,cate], levels = intersect(orders, data[,cate]))
    }
    varis <- c(cate.vars)
  }

  if(!is.null(con.vars)){
    for (con in con.vars) {
      data[,con] <- as.character(data[[con]])
      data[,con][data[,con]==""] <- NA  # "NA" 대신 NA로 처리
      data[,con] <- as.numeric(data[[con]])
    }
    varis <- c(con.vars)
  }

  if(!is.null(cate.vars) & !is.null(con.vars)){
    varis <- unique(c(cate.vars, con.vars))
  }

  #----------------------------------------------------#
  #--------------    regression model     -------------#
  #----------------------------------------------------#
  set.seed(1)
  for (outcome in outcomes) {
    if (!is.null(mul.vars)) {
      # Multivariate analysis

      if (class(data[,outcome]) == "character" && length(unique(data[,outcome])) > 3) {
        print("Not able to analyze due to more than 3 unique levels")
        next
      }

      # Prepare formula
      form <- as.formula(sprintf("%s ~ %s", outcome, mul.vars))

      print("Multivariate analysis")
      type <- "multi"

      print(form)

      # Fit model
      if(!is.null(randomEff)){
        m <- "Regression (LMEM)"
        mod <- lmerTest::lmer(form, data=data)
      } else {
        if (class(data[,outcome]) == "numeric"){
          m <- "Regression (glm-poisson)"
          mod <- glm(form, data=data,  family = poisson(link='log'))
        } else if (length(unique(data[,outcome])) == 2){
          m <- "Logistic regression (glm-binomial)"
          mod <- glm(form, data=data,  family=binomial("logit"))
        }
      }

      print(m)

      # Extract results
      coef <- as.data.frame(summary(mod)$coefficients)
      coef <- coef[setdiff(rownames(coef), "(Intercept)"),,drop=F]

      if(!is.null(randomEff)){
        colnames(coef) <- c("Estimate", "SE", "df","t", "pval")
      }else{
        colnames(coef) <- c("Estimate", "SE", "t", "pval")
      }

      if (dim(coef)[1] == 0){
        next
      }

      conf <- data.frame(confint(mod))
      conf <- conf[setdiff(rownames(conf), "(Intercept)"),,drop=F]
      conf.na <- na.omit(conf)
      if(dim(conf.na)[1] == 0){
        conf.na <- conf
      }

      colnames(conf.na) <- c("2.5 %", "97.5 %")

      if(!is.null(randomEff)){
        coef <- coef
      }else{
        coef$`2.5 %` <- conf.na$`2.5 %`
        coef$`97.5 %` <- conf.na$`97.5 %`
      }

      coef$outcome <- outcome
      coef$model <- m
      coef$formula <- rep(paste(deparse(form), collapse = ""), nrow(coef))

      if(!is.null(randomEff)){
        coef <- coef
      }else{
        coef$deviance <- pchisq(q=mod$null.deviance-mod$deviance, df=mod$df.null-mod$df.residual, lower.tail = FALSE)
      }

      res <- coef
      res$padj <- p.adjust(res$pval, method="fdr")
      res$comp <- factor(rownames(res), levels=rownames(res))
      res$dir <- ifelse(res$pval < pvalue, ifelse(sign(res$Estimate)==1, "up", "down"), "NS")

      write.csv(res, quote = FALSE, col.names = T, file=sprintf("%s/regression_%s.%s.%s.%s%s.csv",
                                                                 out_table, project, outcome, type,
                                                                 ifelse(is.null(name), "", paste(name, ".", sep = "")),
                                                                 format(Sys.Date(), "%y%m%d"), sep="/"))

      saveRDS(mod, sprintf("%s/regression_%s.%s.%s.%s%s.rds",
                           out_table, project, outcome, type,
                           ifelse(is.null(name), "", paste(name, ".", sep = "")),
                           format(Sys.Date(), "%y%m%d"), sep="/"))
    }else {
      # Univariate analysis

      if (class(data[,outcome]) == "character") {
        if (length(unique(data[,outcome])) > 3){
          print("Not able to analyze due to more than 3 unique levels")
          next
        }
        data[,outcome] <- factor(data[,outcome], levels = intersect(orders, data[,outcome]))

        if(length(unique(data[,outcome])) == 2){
          data[,outcome] <- factor(ifelse(data[,outcome]== levels(data[,outcome])[1], 0, 1), levels=c(0,1),
                                   labels = levels(data[,outcome]))
        }
      } else if (class(data[,outcome])  == "numeric") {
        data[,outcome] <- as.character(data[[outcome]])
        data[,outcome][data[,outcome]==""] <- "NA"
        data[,outcome] <- as.numeric(as.character(data[[outcome]]))
      }

      res <- {}

      for (mvar in varis) {
        if (outcome == mvar) next
        if (length(unique(data[,mvar])) == 1) next

        form <- as.formula(sprintf("%s ~ %s%s", outcome,
                                   ifelse(is.null(randomEff), "", sprintf("(1 | %s) +", randomEff)),
                                   mvar))
        print("Univariate analysis")
        type <- "uni"

        print(form)

        if(!is.null(randomEff)){
          m <- "Regression (LMEM)"
          mod <- lmerTest::lmer(form, data=data)
        } else {
          if (class(data[,outcome]) == "numeric"){
            m <- "Regression (glm-poisson)"
            mod <- glm(form, data=data,  family = poisson(link='log'))
          } else if (length(unique(data[,outcome])) == 2){
            m <- "Logistic regression (glm-binomial)"
            mod <- glm(form, data=data,  family=binomial("logit"))
          }
        }

        print(m)

        # Extract results
        coef <- as.data.frame(summary(mod)$coefficients)
        coef <- coef[setdiff(rownames(coef), "(Intercept)"),,drop=F]

        if(!is.null(randomEff)){
          colnames(coef) <- c("Estimate", "SE", "df","t", "pval")
        }else{
          colnames(coef) <- c("Estimate", "SE", "t", "pval")
        }

        if (dim(coef)[

          1] == 0){
          next
        }

        conf <- data.frame(confint(mod))
        conf <- conf[setdiff(rownames(conf), "(Intercept)"),,drop=F]
        conf.na <- na.omit(conf)
        if(dim(conf.na)[1] == 0){
          conf.na <- conf
        }

        colnames(conf.na) <- c("2.5 %", "97.5 %")

        if(!is.null(randomEff)){
          coef <- coef
        }else{
          coef$`2.5 %` <- conf.na$`2.5 %`
          coef$`97.5 %` <- conf.na$`97.5 %`
        }

        coef$outcome <- outcome
        coef$mvar <- mvar
        coef$model <- m

        if(!is.null(randomEff)){
          coef <- coef
        }else{
          coef$deviance <- pchisq(q=mod$null.deviance-mod$deviance, df=mod$df.null-mod$df.residual, lower.tail = FALSE)
        }

        res <- rbind(res, coef)


      }
      res$padj <- p.adjust(res$pval, method="fdr")
      res$comp <- factor(rownames(res), levels=rownames(res))
      res$dir <- ifelse(res$pval < pvalue, ifelse(sign(res$Estimate)==1, "up", "down"), "NS")

      write.csv(res, quote = FALSE, col.names = T, file=sprintf("%s/regression_%s.%s.%s.%s%s.csv",
                                                                 out_table, project, outcome, type,
                                                                 ifelse(is.null(name), "", paste(name, ".", sep = "")),
                                                                 format(Sys.Date(), "%y%m%d"), sep="/"))
    }
  }
}
