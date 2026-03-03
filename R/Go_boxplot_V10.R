
#' Generate Boxplots for Multiple Variables
#'
#' @param df Data frame containing the data to be plotted.
#' @param cate.vars Categorical variables to be used for the boxplot's x-axis.
#' @param project Project name used for output file naming.
#' @param outcomes Numeric variables to be plotted on the y-axis.
#' @param orders Vector of ordered factor levels for the categorical variables.
#' @param mycols Custom color palette for the plots.
#' @param combination Number of group combinations to display in the boxplot.
#' @param ylim Y-axis limits for the plot.
#' @param title Title of the plot.
#' @param facet Optional variable for creating facetted plots.
#' @param paired Indicates if the data points are paired.
#' @param name Optional name for saving the plot.
#' @param ncol Number of columns for facet wrapping in the plot.
#' @param addnumber Boolean to add the number of samples in each group.
#' @param standardsize Boolean to control the size of the plot based on group size.
#' @param statistics Whether to perform statistical tests.
#' @param parametric Whether to use parametric tests.
#' @param model Statistical engine: "nonparametric", "parametric", or "lmm". If NULL, inferred from `parametric`.
#' @param covariates Optional covariate column names for adjusted models (ANCOVA/LMM).
#' @param p_adjust P-value adjustment method for engine-based pairwise tests (e.g., "BH", "bonferroni").
#' @param xangle Angle of x-axis labels.
#' @param cutoff Significance level for statistical tests.
#' @param height Height of the plot.
#' @param width Width of the plot.
#' @param plotCols Number of columns in the multiplot layout.
#' @param plotRows Number of rows in the multiplot layout.
#'
#' @details
#' This function creates boxplots for multiple categorical variables against one or more numeric outcomes. It supports faceting, custom color palettes, and statistical testing with options for parametric or non-parametric methods. The function can handle paired data and allows customization of plot size and layout.
#'
#' @return
#' A PDF file containing the boxplots.
#'
#' @examples
#' Go_boxplot(df = my_data, cate.vars = c("Group", "Condition"), project = "MyProject",
#'            outcomes = c("Variable1", "Variable2"), orders = c("Control", "Treatment"),
#'            height = 10, width = 8, plotCols = 2, plotRows = 1)
#'
#' @export

Go_boxplot <- function(df, cate.vars, project, outcomes,
                       orders=NULL,
                       mycols=NULL,
                       combination=NULL,
                       ylim =NULL,
                       title= NULL,
                       facet= NULL,
                       paired=NULL,
                       name= NULL,
                       ncol=NULL,
                       addnumber=TRUE,
                       standardsize=TRUE,
                       statistics = TRUE,
                       parametric= FALSE,
                       model = NULL,
                       covariates = NULL,
                       p_adjust = "BH",
                       xangle=90,
                       cutoff = 0.1,
                       height, width, plotCols, plotRows){

  has_facet <- function(x) !is.null(x) && length(x) >= 1
  build_comparisons <- function(cbn_mat) {
    comparisons <- vector("list", ncol(cbn_mat))
    for (idx in seq_len(ncol(cbn_mat))) {
      comparisons[[idx]] <- cbn_mat[, idx]
    }
    comparisons
  }
  build_plot_title <- function(base_title, test_name, pval) {
    sprintf("%s%s%s%s",
            base_title,
            ifelse(is.null(test_name), "", paste("\n", test_name, " ", sep = "")),
            ifelse(is.null(pval), "", paste("p=", " ", sep = "")),
            ifelse(is.null(pval), "", paste(pval, " ", sep = "")))
  }
  build_plot_subtitle <- function(stat_res, p_adjust, use_covariates, covariates_label) {
    subtitle_parts <- character(0)
    method_label <- if (!is.null(stat_res$test.name) && stat_res$test.name %in% c("ANCOVA", "LMM")) stat_res$test.name else NULL
    if (!is.null(method_label)) {
      subtitle_parts <- c(subtitle_parts, paste0("method=", method_label))
    }
    if (!is.null(stat_res$annotation) && !is.null(method_label) && method_label %in% c("ANCOVA", "LMM")) {
      subtitle_parts <- c(subtitle_parts, sprintf("pairwise (adjust=%s)", p_adjust))
    }
    if (isTRUE(use_covariates)) {
      subtitle_parts <- c(subtitle_parts, sprintf("covariates=%s", covariates_label))
    }
    if (length(subtitle_parts) == 0) return(NULL)
    paste(subtitle_parts, collapse = "\n")
  }
  sanitize_tag <- function(x) gsub("[^A-Za-z0-9._-]+", "-", x)

  if(!is.null(dev.list())) dev.off()

  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out_path)) dir.create(out_path)
  set.seed(151)


  # out file
  # "name" definition
  if (class(name) == "function"){
    name <- NULL
  }

  tt <- try(mycols,T)
  if(class(tt) == "try-error"){
    print("mycols is not defined.")
    mycols <- NULL
  }

  tt <- try(orders,T)
  if(class(tt) == "try-error"){
    print("orders is not defined.")
    orders <- NULL
  }

  print("box.tickness1")

  has_covariates <- !is.null(covariates) && length(covariates) > 0
  model_is_auto <- is.null(model)
  resolved_model <- if (model_is_auto) {
    if (has_covariates) {
      "parametric"
    } else if (isTRUE(parametric)) {
      "parametric"
    } else {
      "nonparametric"
    }
  } else {
    tolower(as.character(model))
  }
  if (!resolved_model %in% c("nonparametric", "parametric", "lmm")) {
    stop("`model` must be one of: 'nonparametric', 'parametric', 'lmm'.")
  }
  covariates_label <- if (has_covariates) paste(covariates, collapse = "+") else NULL
  use_covariates <- !is.null(covariates_label) && resolved_model %in% c("parametric", "lmm")
  if (!is.null(covariates_label) && resolved_model == "nonparametric" && !model_is_auto) {
    warning("`covariates` are ignored when model='nonparametric'.")
  }
  if (resolved_model == "lmm" && is.null(paired)) {
    warning("model='lmm' requested but `paired` is NULL. LMM cannot be fitted.")
  }
  method_file_tag <- if (resolved_model == "lmm" && !is.null(paired)) {
    "(method=LMM)."
  } else if (isTRUE(use_covariates) && resolved_model == "parametric") {
    "(method=ANCOVA)."
  } else {
    ""
  }
  covariates_file_tag <- if (isTRUE(use_covariates)) paste0("(cov=", sanitize_tag(covariates_label), ").") else ""

  # plot design
  if (height*width <= 6){
    dot.size = 0.7
    box.tickness = 0.3
  }else if (height*width > 6 & height*width < 10){
    dot.size = 1
    box.tickness = 0.4
  }else{
    dot.size = 1.5
    box.tickness = 0.5
  }

  file_name <- paste0(
    "box.", project, ".",
    ifelse(is.null(facet), "", paste(facet, ".", sep = "")),
    ifelse(is.null(paired), "", paste("(paired=",paired, ").", sep = "")),
    ifelse(is.null(combination), "", paste("(cbn=",combination, ").", sep = "")),
    method_file_tag,
    covariates_file_tag,
    ifelse(is.null(name), "", paste(name, ".", sep = "")),
    format(Sys.Date(), "%y%m%d"),
    ".pdf"
  )
  pdf(file.path(out_path, file_name), height = height, width = width)

  # plot
  plotlist <- list()
  for (mvar in cate.vars) {
    if (length(unique(df[,mvar])) < 2){
      next
    }

    if (has_facet(facet)){
      if (facet == mvar){
        next
      }
    }

    # remove Na
    print("Control NA")
    df <- data.frame(df)
    df[,mvar] <- as.character(df[,mvar]);df[,mvar]
    df[,mvar][df[,mvar]==""] <- "NA";df[,mvar]
    df.na <- subset(df, df[,mvar] != "NA");df.na[,mvar]  # subset 를 사용한 NA 삭제
    df.na[,mvar] <- as.factor(df.na[,mvar]);df.na[,mvar]

    df.na[,mvar] <- factor(df.na[,mvar], levels = intersect(orders, df.na[,mvar]))


    # Add number of samples in the group
    print("Add sample number informations")
    if(addnumber==TRUE){
    renamed_levels <- as.character(levels(df.na[,mvar]));renamed_levels
    oldNames <- unique(df.na[,mvar]);oldNames
    if (length(renamed_levels) == 0) {
      renamed_levels <- oldNames
    }
    for (grp_name in oldNames) {
      total <- length(which(df.na[,mvar] == grp_name));total
      new_n <- paste(grp_name, " (n=", total, ")", sep="");new_n
      levels(df.na[[mvar]])[levels(df.na[[mvar]])== grp_name] <- new_n
      renamed_levels <- replace(renamed_levels, renamed_levels == grp_name, new_n);renamed_levels
    }
    }



    print(sprintf("##-- %s (total without NA: %s/%s) --##",
                  mvar, dim(df.na)[1], dim(df)[1]))

    if (length(unique(df.na[,mvar])) ==1) {
      next
    }

    summary.df.na <- summary(df.na[,mvar])

    #------------------------------#
    # for group combination or not #
    #------------------------------#


    if (!is.null(combination)){
      print(sprintf("Combination n=%s", combination))
      group.cbn <- combn(x = levels(df.na[,mvar]), m = combination)

      #print(count(group.cbn))

      group_comparisons <- build_comparisons(group.cbn);group_comparisons

      print(1)
      for(i in 1:length(group_comparisons)){
        print(group_comparisons[i])
        group.combination <- unlist(group_comparisons[i]);group.combination

        if (combination >= 2 && combination <= 10){
          df.cbn <- subset(df.na, df.na[,mvar] %in% group.combination[1:combination])
        }  else{
          print("combination should be 2~10 only.")
          break
        }

        unique(df.cbn[,mvar])


        # make a comnination for stat
        df.cbn[,mvar] <- factor(df.cbn[,mvar])
        cbn <- combn(x = levels(df.cbn[,mvar]), m = 2)


        my_comparisons <- build_comparisons(cbn);my_comparisons

        if(combination != 2){
          combination.N <- combination - 1
          my_comparisons <- my_comparisons[1:combination.N]
        }




        for(oc in outcomes){
          # remove NA for facet
          if (has_facet(facet)) {
            for (fc in facet){
              df.cbn[,fc] <- as.character(df.cbn[,fc]);df.cbn[,fc]
              df.cbn[,fc][df.cbn[,fc] == ""] <- "NA"
              df.cbn.sel <- df.cbn[!is.na(df.cbn[,fc]), ]
              df.cbn <- df.cbn.sel
              # facet or not
              df.cbn[,fc] <- factor(df.cbn[,fc], levels = orders)
            }
          }

          print(oc)

          if (statistics){
            stat_res <- Go_boxplot_stats_engine(
              df = df.cbn,
              mvar = mvar,
              oc = oc,
              comparisons = my_comparisons,
              model = model,
              parametric = parametric,
              covariates = covariates,
              paired = paired,
              p_adjust = p_adjust
            )
            test.name <- stat_res$test.name
            pval <- stat_res$pval
          }else{
            stat_res <- list(test.name = NULL, pval = NULL, testmethod = NULL, annotation = NULL)
            test.name <- NULL
            pval <- NULL
          }


          p1 <- ggplot(df.cbn, aes_string(x = mvar, y = oc))  + labs(y=oc, x=NULL) + #theme_bw() +
            theme(strip.background = element_blank()) +
            theme(text=element_text(size=8), axis.text.x=element_text(angle=xangle,hjust=1,vjust=0.5,size=8),
                  plot.title=element_text(size=8),
                  plot.subtitle=element_text(size=6, lineheight=0.9)) # ,face="bold"



          title_text <- build_plot_title(ifelse(is.null(title), mvar, title), test.name, pval)
          subtitle_text <- if (statistics) build_plot_subtitle(stat_res, p_adjust, use_covariates, covariates_label) else NULL
          p1 <- p1 + labs(title = title_text, subtitle = subtitle_text)

          if (statistics){
            p1 <- Go_boxplot_add_stats_layer(
              p1 = p1,
              stat_res = stat_res,
              my_comparisons = my_comparisons,
              paired = paired,
              cutoff = cutoff
            )
          }

          if(!is.null(ylim)){
            if(oc == "Chao1"){
              p1 = p1
            }else{
              p1 = p1 + ylim(ylim[1] , ylim[2])
            }
          }


          # paired plot type
          if (!is.null(paired)) {

            if(!is.null(mycols)){
              p1 <- p1 + scale_color_manual(values = mycols)
            }else{
              p1 <- p1
            }

            p1 = p1 + geom_boxplot(aes_string(colour = mvar), outlier.shape = NA, lwd = box.tickness)  + theme(legend.position="none")
            p1 = p1 + geom_point(aes_string(colour = mvar, group = paired), alpha = 0.8, size = dot.size, position = position_dodge(0.3), show.legend = F)  #scale_shape_manual(values = c(1, 16, 8, 0,15, 2,17,11, 10,12,3,4,5,6,7,8,9,13,14))
            p1 = p1 + geom_line(aes_string(group = paired), color = "grey50", linewidth = 0.3, position = position_dodge(0.3))
            p1 = p1 + theme(legend.title = element_blank(), legend.position="bottom", legend.justification="left",legend.box.margin = ggplot2::margin(0,0,0,-1,"cm"))
          }  else{

            # count or table for number of variable
            if (max(table(df.cbn[,mvar])) > 150){

              if(!is.null(mycols)){
                p1 <- p1 + scale_fill_manual(values = mycols)
              }else{
                p1 <- p1
              }
              p1 = p1 + geom_boxplot(aes_string(fill = mvar), outlier.shape = NULL, lwd = box.tickness)   + theme(legend.position="none")
              # outlier.shape = NA
            } else {

              if(!is.null(mycols)){
                p1 <- p1 + scale_color_manual(values = mycols)
              }else{
                p1 <- p1
              }

              p1 = p1 + geom_boxplot(aes_string(colour = mvar), outlier.shape = NA, lwd = box.tickness) + theme(legend.position="none")
              p1 = p1 + geom_jitter(aes_string(colour = mvar), shape = 16, alpha = 0.8, size = dot.size, position = position_jitter(0.2))

            }
          }

          # facet
          if (has_facet(facet)) {
            if(is.null(ncol)){
              ncol <- length(unique(df[,facet]))
            }


            p1 = p1 + facet_wrap(as.formula(sprintf("~ %s" , paste(setdiff(facet, "SocpleType"), collapse="+"))),
                                 scales="free_x", ncol = ncol)
            p1 = p1 + guides(color = "none", size = "none", shape= "none")
          } else {
            p1 = p1 + guides(color = "none", size = "none", shape= "none")
          }


          # plot size ratio
          if (length(unique(df.cbn[,mvar])) < 5){
            if(standardsize==TRUE){
              num.subgroup <- length(unique(df.cbn[,mvar]))*0.1
            }else{
              num.subgroup <- 0.9
            }
          }else{
            num.subgroup <- length(unique(df.cbn[,mvar]))*0.1
          }

          p1  <- p1  + theme(panel.grid = element_blank(),
                             panel.background = element_rect(fill = "white", colour = "Black", linewidth = 0.5, linetype = "solid"),
                             aspect.ratio = 1/num.subgroup)

          plotlist[[length(plotlist)+1]] <- p1
        }
      }
    }else{
      # make a comnination for stat

      print("Check combination for statistics")
      cbn <- combn(x = levels(df.na[,mvar]), m = 2)

      my_comparisons <- build_comparisons(cbn);my_comparisons
      # check statistics method
      for(oc in outcomes){
        # remove NA for facet
        if (!is.null(facet)) {
          for (fc in facet){
            df.na[,fc] <- as.character(df.na[,fc]);df.na[,fc]
            df.na[,fc][df.na[,fc] == ""] <- "NA"
            df.na.sel <- df.na[!is.na(df.na[,fc]), ]
            df.na <- df.na.sel
            # facet or not
            df.na[,fc] <- factor(df.na[,fc], levels = orders)
          }
        }

        print(oc)

        if (statistics){
          stat_res <- Go_boxplot_stats_engine(
            df = df.na,
            mvar = mvar,
            oc = oc,
            comparisons = my_comparisons,
            model = model,
            parametric = parametric,
            covariates = covariates,
            paired = paired,
            p_adjust = p_adjust
          )
          test.name <- stat_res$test.name
          pval <- stat_res$pval
        }else{
          stat_res <- list(test.name = NULL, pval = NULL, testmethod = NULL, annotation = NULL)
          test.name <- NULL
          pval <- NULL
        }


        p1 <- ggplot(df.na, aes_string(x = mvar, y = oc))  + labs(y=oc, x=NULL) +
          theme(strip.background = element_blank()) +
          theme(text=element_text(size=8), axis.text.x=element_text(angle=xangle,hjust=1,vjust=0.5,size=8),
                plot.title=element_text(size=8),
                plot.subtitle=element_text(size=6, lineheight=0.9))





        # paired plot type
        if (!is.null(paired)) {

          if(!is.null(mycols)){
            p1 <- p1 + scale_color_manual(values = mycols)
          }else{
            p1 <- p1
          }

          p1 = p1 + geom_boxplot(aes_string(colour = mvar), outlier.shape = NA, lwd = box.tickness)  + theme(legend.position="none")
          p1 = p1 + geom_point(aes_string(colour = mvar, group = paired), alpha = 0.8, size = dot.size, position = position_dodge(0.3), show.legend = F)   #scale_shape_manual(values = c(1, 16, 8, 0,15, 2,17,11, 10,12,3,4,5,6,7,8,9,13,14))
          p1 = p1 + geom_line(aes_string(group = paired), color = "grey50", linewidth = 0.3, position = position_dodge(0.3))
          p1 = p1 + theme(legend.title = element_blank(), legend.position="bottom", legend.justification="left",legend.box.margin = ggplot2::margin(0,0,0,-1,"cm"))

        } else{
          # count or table for number of variable
          if (max(table(df.na[,mvar])) > 150){

            if(!is.null(mycols)){
              p1 <- p1 + scale_fill_manual(values = mycols)
            }else{
              p1 <- p1
            }
            p1 = p1 + geom_boxplot(aes_string(fill = mvar), outlier.shape = NULL, lwd = box.tickness)   + theme(legend.position="none")

          } else {

            if(!is.null(mycols)){
              p1 <- p1 + scale_color_manual(values = mycols)
            }else{
              p1 <- p1
            }

            p1 = p1 + geom_boxplot(aes_string(colour = mvar), outlier.shape = NA, lwd = box.tickness) + theme(legend.position="none")
            p1 = p1 + geom_jitter(aes_string(colour = mvar), shape = 16, alpha = 0.8, size = dot.size, position = position_jitter(0.2))

          }
        }



        if (statistics){
          p1 <- Go_boxplot_add_stats_layer(
            p1 = p1,
            stat_res = stat_res,
            my_comparisons = my_comparisons,
            paired = paired,
            cutoff = cutoff
          )
        }



        # Close an image
        title_text <- build_plot_title(ifelse(is.null(title), mvar, title), test.name, pval)
        subtitle_text <- if (statistics) build_plot_subtitle(stat_res, p_adjust, use_covariates, covariates_label) else NULL
        p1 <- p1 + labs(title = title_text, subtitle = subtitle_text)

        # y axis limit
        if(!is.null(ylim)){
          if(oc == "Chao1"){
            p1 = p1
          }else{
            p1 = p1 + ylim(ylim[1] , ylim[2])
          }
        }
        # facet
        if (has_facet(facet)) {
          if(is.null(ncol)){
            ncol <- length(unique(df[,facet]))
          }


          p1 = p1 + facet_wrap(as.formula(sprintf("~ %s" , paste(setdiff(facet, "SocpleType"), collapse="+"))), scales="free_x", ncol = ncol)
          p1 = p1 + guides(color = "none", size = "none", shape= "none")
        } else {
          p1 = p1 + guides(color = "none", size = "none", shape= "none")
        }

        # plot size ratio
        if (length(unique(df.na[,mvar])) < 5){
          if(standardsize==TRUE){
            num.subgroup <- length(unique(df.na[,mvar]))*0.1
          }else{
            num.subgroup <- 0.9
          }
        }else{
          num.subgroup <- length(unique(df.na[,mvar]))*0.1
        }

        p1  <- p1  + theme(panel.grid = element_blank(),
                           panel.background = element_rect(fill = "white", colour = "Black", linewidth = 0.5, linetype = "solid"),
                           aspect.ratio = 1/num.subgroup)

        plotlist[[length(plotlist)+1]] <- p1

      }
    }
  }
  multiplot <- function(..., plotlist=NULL, file, cols=1, rows=1) {
    require(grid)
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    numPlots = length(plots)

    i = 1
    while (i < numPlots) {
      numToPlot <- min(numPlots-i+1, cols*rows)
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(i, i+cols*rows-1), ncol = cols, nrow = rows, byrow=T)
      if (numToPlot==1) {
        print(plots[[i]])
      } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        # Make each plot, in the correct location
        for (j in i:(i+numToPlot-1)) {
          # Get the i,j matrix positions of the regions that contain this subplot
          matchidx <- as.data.frame(which(layout == j, arr.ind = TRUE))
          print(plots[[j]], vp = viewport(layout.pos.row = matchidx$row,
                                          layout.pos.col = matchidx$col))
        }
      }
      i <- i+numToPlot
    }
  }
  multiplot(plotlist=plotlist, cols=plotCols, rows=plotRows)
  dev.off()
}
