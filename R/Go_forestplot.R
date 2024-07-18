#' Go_forestplot
#'
#' This function generates forest plots for linear mixed models (LMM) based on the given outcomes and formula. It saves the plots as PDF files in an output directory.
#'
#' @param df A data frame containing the data to be used for the models.
#' @param project A character string specifying the project name. Used for naming output directories and files.
#' @param outcomes A character vector specifying the names of the outcome variables to be modeled.
#' @param formular A character string representing the formula for the model. It should contain a placeholder ("%s") for the outcome variable.
#' @param name An optional character string to include in the output file names. Default is NULL.
#' @param height A numeric value specifying the height of the output plots in inches.
#' @param width A numeric value specifying the width of the output plots in inches.
#' 
#' @import lme4
#' @import ggplot2
#' 
#' @return This function does not return a value. It creates and saves forest plot PDFs in the specified directory.
#' @examples
#' \dontrun{
#' df <- your_data_frame
#' outcomes <- c("outcome1", "outcome2")
#' formular <- "%s ~ predictor1 + predictor2 + (1|random_effect)"
#' Go_forestplot(df, project = "MyProject", outcomes = outcomes, formular = formular, height = 6, width = 8)
#' }

Go_forestplot <- function(df, 
                          project, 
                          outcomes, 
                          formular,
                          name=NULL,
                          height,
                          width){
  
  #===== out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  for (outcome in outcomes){
    form <-  as.formula(sprintf(formular, outcome))
    mod <- lmer(form, data = df)
    
    # summary(mod)
    
    # ===== detection method from model
    patterns <- list(
      glm = "\\bglm\\(",  # Ensures "glmer(" is treated as a separate word
      lm = "\\blm\\(", 
      glmer = "\\bglmer\\(",  # Ensures "glmer(" is treated as a separate word
      lmer = "\\blmer\\("     # Ensures "lmer(" is treated as a separate word
    )
    
    # Function to detect model type
    detect_model_type <- function(model_call_string, patterns) {
      detected_models <- sapply(patterns, function(pat) grepl(pat, model_call_string))
      names(detected_models[detected_models == TRUE])
    }
    
    model_call <- summary(mod)$call
    model_call_string <- deparse(model_call)
    # Use the function to detect model type again
    model_call_string <- model_call_string[1]
    model_type <- detect_model_type(model_call_string, patterns)
    
    # Print the detected model types
    print(model_type)
    
    # ===== extract table
    coef <- as.data.frame(summary(mod)$coefficients)
    coef <- coef[setdiff(rownames(coef), "(Intercept)"),,drop=F]
    
    if(model_type == "lm" | model_type == "glm"){
      colnames(coef) <- c("Estimate", "SE", "df","stat", "pval")
    }else if(model_type =="lmer") {
      colnames(coef) <- c("Estimate", "SE", "df","stat", "pval")
    }
    
    coef$padj <- p.adjust(coef$pval, method="fdr")
    print(coef)
    
    # ===== calculation for the confidence interval 
    coef$lower_ci <- coef$Estimate - qt(0.975, df=Inf) * coef$SE
    coef$upper_ci <- coef$Estimate + qt(0.975, df=Inf) * coef$SE
    
    # colnames(conf.na) <- c("2.5 %", "97.5 %")
    
    coef$outcome <- outcome
    coef$model <- model_type
    
    print(coef)
    
    coef$dir <- ifelse(coef$pval < 0.05, ifelse(sign(coef$Estimate)== 1, "up", "down"), "NS")
    coef$dirPadj <- ifelse(coef$padj < 0.05, TRUE, FALSE)
    
    padj_shape <- c(19,1); names(padj_shape) <- c(TRUE, FALSE)
    
    mycols <- c("#2366C0FF", "#FF6435FF")
    
    if(!is.null(mycols)){
      dircolors <- c("down" = mycols[1], "NS" = "grey", "up" = mycols[2])
      
    }else{
      dircolors <- c("down" = "#f8766d", "NS" = "grey", "up" = "#7cae00")
    }
    
    coef$fix_effect <-rownames(coef) ; head(coef)
    
    # Forest plot을 위한 색상 설정
    lims <- max(abs(coef$upper_ci - coef$lower_ci)) * 1.1
    
    # Forest plot 그리기
    forest_plot <- ggplot(coef, aes(x=reorder(fix_effect, Estimate), y=Estimate, ymin=lower_ci, ymax=upper_ci)) + #reorder(fix_effect, Estimate)
      geom_pointrange(aes(color=dir, shape=dirPadj), size=0.5) +
      geom_hline(yintercept=0, linetype="dashed", color = "black") +
      coord_flip() + 
      scale_color_manual(values=dircolors) +
      scale_shape_manual(values=c("TRUE" = 19, "FALSE" = 1)) +  # Closed circle for significant, open for not
      labs(title=sprintf("Forest Plot of %s LMM Results", outcome), y="Effect Size", x="Variable") +
      #theme_minimal()
      theme_classic()
    
    # Plot 출력
    print(forest_plot)
    
    ggsave(sprintf("%s/forest.%s.%s.%s%s.pdf", out_path, project, 
                   outcome,
                   ifelse(is.null(name), "", paste(name, ".", sep = "")), 
                   format(Sys.Date(), "%y%m%d")), plot = forest_plot, device = "pdf", width = width, height = height)
  }
}


