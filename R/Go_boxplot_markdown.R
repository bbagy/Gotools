#' Create Customized Boxplots for Given Categories and Outcomes
#'
#' This function generates boxplots for specified categorical variables and outcomes from a data frame.
#' It includes options for ordering, color customization, combination of groups, statistical annotations,
#' and layout adjustments.
#'
#' @param df The data frame containing the variables and outcomes for the boxplots.
#' @param cate.vars A vector of column names in 'df' that represent categorical variables for boxplot grouping.
#' @param project A character string representing the project name or identifier.
#' @param outcomes A vector of column names in 'df' that represent outcomes to be plotted.
#' @param orders Optional ordering of the factor levels for the categorical variables.
#' @param mycols Optional custom color vector for the boxplots.
#' @param combination Optional number of group combinations to compare in the boxplots.
#' @param ylim Optional y-axis limits for the boxplots.
#' @param title Optional title for the boxplots.
#' @param facet Optional column name in 'df' for faceting the boxplots.
#' @param paired Optional column name in 'df' for creating paired boxplots.
#' @param name Optional name for saving the output plots.
#' @param addnumber Logical indicating whether to add sample size to group labels.
#' @param standardsize Logical indicating whether to use standard size ratios for the plots.
#' @param statistics A character string ("yes" or "no") indicating whether to perform statistical tests.
#' @param parametric A character string ("yes" or "no") indicating whether to use parametric tests.
#' @param star A character string ("yes" or "no") indicating whether to use star annotations for significance.
#' @param xangle Numeric value for the angle of x-axis labels.
#' @param font Numeric value specifying the font size for text elements in the plot.
#' @param cutoff Numeric value specifying the significance cutoff for statistical annotations.
#' @param plotCols Numeric value specifying the number of columns in the multiplot layout.
#' @param plotRows Numeric value specifying the number of rows in the multiplot layout.
#'
#' @return The function generates boxplots and optionally saves them but does not return a value.
#'
#' @examples
#' # Assuming 'df' is a data frame with appropriate data
#' Go_boxplot_markdown(df, cate.vars = c("Group"), project = "MyProject", outcomes = c("Outcome1"),
#'                     orders = NULL, mycols = c("#00BFC4", "#F8766D"), combination = NULL,
#'                     ylim = c(0, 10), title = "Boxplot Example", facet = NULL, paired = NULL,
#'                     name = "Boxplot", addnumber = TRUE, standardsize = TRUE, statistics = "yes",
#'                     parametric = "no", star = "no", xangle = 90, font = 10, cutoff = 0.1,
#'                     plotCols = 2, plotRows = 1)
#'
#' @export
#' @importFrom ggplot2 ggplot aes_string xlab ylab geom_boxplot geom_jitter facet_wrap theme
#' @importFrom ggsignif stat_compare_means


Go_boxplot_markdown <- function(df, cate.vars, project, outcomes,
                                orders=NULL,
                                mycols=NULL,
                                combination=NULL,
                                ylim =NULL,
                                title= NULL,
                                facet= NULL,
                                paired=NULL,
                                name= NULL,
                                addnumber=TRUE,
                                standardsize=TRUE,
                                statistics = "yes",
                                parametric= "no",
                                star="no",
                                xangle=90,
                                font=10,
                                cutoff = 0.1,  plotCols, plotRows){

  # out file
  # "name" definition
  if (class(name) == "function"){
    name <- NULL
  }

  tt <- try(mycols,T)
  if(class(tt) == "try-error"){
    mycols <- NULL
  }

  tt <- try(orders,T)
  if(class(tt) == "try-error"){
    orders <- NULL
  }

  # plot design
  "if (height*width <= 6){
   dot.size = 0.7
    box.tickness = 0.3
  }else if (height*width > 6 & height*width < 10){
    dot.size = 1
    box.tickness = 0.4
  }else{
    dot.size = 1.5
    box.tickness = 0.5
  }"


  dot.size = 2
  box.tickness = 1


  # plot
  plotlist <- list()
  for (mvar in cate.vars) {
    if (length(unique(df[,mvar])) < 2){
      next
    }

    if (length(facet) >= 1){
      if (facet == mvar){
        next
      }
    } else {}

    # remove Na
    df <- data.frame(df)
    df[,mvar] <- as.character(df[,mvar]);df[,mvar]
    df[,mvar][df[,mvar]==""] <- "NA";df[,mvar]
    df.na <- subset(df, df[,mvar] != "NA");df.na[,mvar]  # subset 를 사용한 NA 삭제
    df.na[,mvar] <- as.factor(df.na[,mvar]);df.na[,mvar]

    df.na[,mvar] <- factor(df.na[,mvar], levels = intersect(orders, df.na[,mvar]))


    # Add number of samples in the group
    if(is.null(facet) && addnumber == TRUE) {
      renamed_levels <- as.character(levels(df.na[,mvar]));renamed_levels
      oldNames <- unique(df.na[,mvar]);oldNames
      if (length(renamed_levels) == 0) {
        renamed_levels <- oldNames
      }
      for (name in oldNames) {
        total <- length(which(df.na[,mvar] == name));total
        new_n <- paste(name, " (n=", total, ")", sep="");new_n
        levels(df.na[[mvar]])[levels(df.na[[mvar]])== name] <- new_n
        renamed_levels <- replace(renamed_levels, renamed_levels == name, new_n);renamed_levels
      }
    }else{
      df.na <- df.na
    }


    if (length(unique(df.na[,mvar])) ==1) {
      next
    }

    summary.df.na <- summary(df.na[,mvar])

    #------------------------------#
    # for group combination or not #
    #------------------------------#


    if (!is.null(combination)){
      group.cbn <- combn(x = levels(df.na[,mvar]), m = combination)

      #print(count(group.cbn))

      group_comparisons <- {}
      for(i in 1:ncol(group.cbn)){
        x <- group.cbn[,i]
        group_comparisons[[i]] <- x
      };group_comparisons

      for(i in 1:length(group_comparisons)){
        #print(group_comparisons[i])
        group.combination <- unlist(group_comparisons[i]);group.combination

        if(combination ==2){
          basline <- group.combination[1]
          smvar <- group.combination[2]
          df.cbn <- subset(df.na, df.na[,mvar] %in% c(basline,smvar))
        } else if(combination ==3){
          basline <- group.combination[1]
          smvar1 <- group.combination[2]
          smvar2 <- group.combination[3]
          df.cbn <- subset(df.na, df.na[,mvar] %in% c(basline,smvar1, smvar2))
        }else if(combination ==4){
          basline <- group.combination[1]
          smvar1 <- group.combination[2]
          smvar2 <- group.combination[3]
          smvar3 <- group.combination[4]
          df.cbn <- subset(df.na, df.na[,mvar] %in% c(basline,smvar1, smvar2,smvar3))
        }else if(combination ==5){
          basline <- group.combination[1]
          smvar1 <- group.combination[2]
          smvar2 <- group.combination[3]
          smvar3 <- group.combination[4]
          smvar4 <- group.combination[5]
          df.cbn <- subset(df.na, df.na[,mvar] %in% c(basline,smvar1, smvar2,smvar3,smvar4))
        }else if(combination ==6){
          basline <- group.combination[1]
          smvar1 <- group.combination[2]
          smvar2 <- group.combination[3]
          smvar3 <- group.combination[4]
          smvar4 <- group.combination[5]
          smvar5 <- group.combination[6]
          df.cbn <- subset(df.na, df.na[,mvar] %in% c(basline,smvar1, smvar2,smvar3,smvar4,smvar5))
        }else if(combination ==7){
          basline <- group.combination[1]
          smvar1 <- group.combination[2]
          smvar2 <- group.combination[3]
          smvar3 <- group.combination[4]
          smvar4 <- group.combination[5]
          smvar5 <- group.combination[6]
          smvar6 <- group.combination[7]
          df.cbn <- subset(df.na, df.na[,mvar] %in% c(basline,smvar1, smvar2,smvar3,smvar4,smvar5,smvar6))
        }else if(combination ==8){
          basline <- group.combination[1]
          smvar1 <- group.combination[2]
          smvar2 <- group.combination[3]
          smvar3 <- group.combination[4]
          smvar4 <- group.combination[5]
          smvar5 <- group.combination[6]
          smvar6 <- group.combination[7]
          smvar7 <- group.combination[8]
          df.cbn <- subset(df.na, df.na[,mvar] %in% c(basline,smvar1, smvar2,smvar3,smvar4,smvar5,smvar6,smvar7))
        }else if(combination ==9){
          basline <- group.combination[1]
          smvar1 <- group.combination[2]
          smvar2 <- group.combination[3]
          smvar3 <- group.combination[4]
          smvar4 <- group.combination[5]
          smvar5 <- group.combination[6]
          smvar6 <- group.combination[7]
          smvar7 <- group.combination[8]
          smvar8 <- group.combination[9]
          df.cbn <- subset(df.na, df.na[,mvar] %in% c(basline,smvar1, smvar2,smvar3,smvar4,smvar5,smvar6,smvar7,smvar8))
        }else if(combination ==10){
          basline <- group.combination[1]
          smvar1 <- group.combination[2]
          smvar2 <- group.combination[3]
          smvar3 <- group.combination[4]
          smvar4 <- group.combination[5]
          smvar5 <- group.combination[6]
          smvar6 <- group.combination[7]
          smvar7 <- group.combination[8]
          smvar8 <- group.combination[9]
          smvar9 <- group.combination[10]
          df.cbn <- subset(df.na, df.na[,mvar] %in% c(basline,smvar1, smvar2,smvar3,smvar4,smvar5,smvar6,smvar7,smvar8,smvar9))
        }  else{
          break
        }

        unique(df.cbn[,mvar])


        # make a comnination for stat
        df.cbn[,mvar] <- factor(df.cbn[,mvar])
        cbn <- combn(x = levels(df.cbn[,mvar]), m = 2)


        my_comparisons <- {}
        for(i in 1:ncol(cbn)){
          x <- cbn[,i]
          my_comparisons[[i]] <- x
        };my_comparisons

        if(combination != 2){
          combination.N <- combination - 1
          my_comparisons <- my_comparisons[1:combination.N]
        }




        for(oc in outcomes){
          # remove NA for facet
          if (length(facet) >= 1) {
            for (fc in facet){
              df.cbn[,fc] <- as.character(df.cbn[,fc]);df.cbn[,fc]
              df.cbn[,fc][df.cbn[,fc] == ""] <- "NA"
              df.cbn.sel <- df.cbn[!is.na(df.cbn[,fc]), ]
              df.cbn <- df.cbn.sel
              # facet or not
              df.cbn[,fc] <- factor(df.cbn[,fc], levels = orders)
            }
          }

          # check statistics method
          if (statistics == "yes"| statistics == "YES"|statistics == "Yes"){
            if (parametric == "yes"| parametric == "YES"|parametric == "Yes"){
              if (nlevels(factor(df.cbn[,mvar])) > 2) {
                test <- aov(as.formula(sprintf("%s ~ %s", oc, mvar)), df.cbn)
                pval <- round(summary(test)[[1]][["Pr(>F)"]][1],4)
                test.name <- "ANOVA"
                testmethod <-  "t.test"
              } else {
                testmethod <-  "t.test"
                pval <- NULL
                test.name <- "Pairwise T-Test"
              }
            }else{
              if (nlevels(factor(df.cbn[,mvar])) > 2) {
                test <- kruskal.test(as.formula(sprintf("%s ~ %s", oc, mvar)), df.cbn)
                pval <- round(test$p.value, 4)
                test.name <- "KW"
                testmethod <- "wilcox.test"
              } else {
                testmethod <- "wilcox.test"
                pval <- NULL
                test.name <- "Pairwise Wilcoxon"
              }
            }
          }else{
            test.name<-NULL
            pval <- NULL
          }


          p1 <- ggplot(df.cbn, aes_string(x=mvar, y=oc))  + labs(y=oc, x=NULL) + #theme_bw() +
            theme(strip.background = element_blank()) +
            theme(text=element_text(size=font), axis.text.x=element_text(angle=xangle,hjust=1,vjust=0.5,size=font),
                  plot.title=element_text(size=font)) # ,face="bold"



          if (!is.null(title)) {
            p1 <- p1 + ggtitle(sprintf("%s%s%s%s", title,
                                       ifelse(is.null(test.name), "", paste("\n",test.name, " ", sep = "")),
                                       ifelse(is.null(pval), "", paste("p=", " ", sep = "")),
                                       ifelse(is.null(pval), "", paste(pval, " ", sep = "")), sep=""))
          } else{
            p1 <- p1 + ggtitle(sprintf("%s%s%s%s", mvar,
                                       ifelse(is.null(test.name), "", paste("\n",test.name, " ", sep = "")),
                                       ifelse(is.null(pval), "", paste("p=", " ", sep = "")),
                                       ifelse(is.null(pval), "", paste(pval, " ", sep = "")), sep=""))
          }

          # control statistic on the plot

          if(is.null(test.name)){
            p1 <- p1
          } else if(test.name == "KW" | test.name == "ANOVA"){
            if(pval < cutoff){
              if (statistics == "yes"| statistics == "YES"|statistics == "Yes"){
                if (star == "no") {
                  p1 <- p1 + stat_compare_means(method= testmethod, label = "p.format", comparisons = my_comparisons, size = 2)
                }  else if (star == "yes") {
                  p1 <- p1 + stat_compare_means(method= testmethod, label = "p.signif", comparisons = my_comparisons, hide.ns = TRUE, size = 3)
                }
              }else if(statistics == "no"| statistics == "NO"|statistics == "No"){
                p1 <- p1
              }
            }else {
              p1 <- p1
            }
          }else if(testmethod == "wilcox.test" | testmethod == "t.test"){
            if (statistics == "yes"| statistics == "YES"|statistics == "Yes"){
              if (star == "no") {
                p1 <- p1 + stat_compare_means(method= testmethod, label = "p.format", comparisons = my_comparisons, size = 2)
              }  else if (star == "yes") {
                p1 <- p1 + stat_compare_means(method= testmethod, label = "p.signif", comparisons = my_comparisons, hide.ns = TRUE, size = 3)
              }
            }else if(statistics == "no"| statistics == "NO"|statistics == "No"){
              p1 <- p1
            }
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

            p1 = p1 + geom_boxplot(aes_string(colour=mvar),outlier.shape = NA,lwd=box.tickness)  + theme(legend.position="none")
            p1 = p1 + geom_point(aes_string(fill=mvar,group=paired),alpha = 0.8, size = dot.size, position = position_dodge(0.3),show.legend = F)  #scale_shape_manual(values = c(1, 16, 8, 0,15, 2,17,11, 10,12,3,4,5,6,7,8,9,13,14))
            p1 = p1 + geom_line(aes_string(group=paired), color="grey50", size=0.3,position = position_dodge(0.3))
            p1 = p1 + theme(legend.title = element_blank(), legend.position="bottom", legend.justification="left",legend.box.margin = ggplot2::margin(0,0,0,-1,"cm"))
          }  else{

            # count or table for number of variable
            if (max(table(df.cbn[,mvar])) > 100){

              if(!is.null(mycols)){
                p1 <- p1 + scale_fill_manual(values = mycols)
              }else{
                p1 <- p1
              }
              p1 = p1 + geom_boxplot(aes_string(fill=mvar),outlier.shape = NULL,lwd=box.tickness)   + theme(legend.position="none")
              # outlier.shape = NA
            } else {

              if(!is.null(mycols)){
                p1 <- p1 + scale_color_manual(values = mycols)
              }else{
                p1 <- p1
              }

              p1 = p1 + geom_boxplot(aes_string(colour=mvar),outlier.shape = NA,lwd=box.tickness) + theme(legend.position="none")
              p1 = p1 + geom_jitter(aes_string(colour=mvar),shape=16, alpha = 0.8, size = dot.size, position=position_jitter(0.2))

            }
          }

          # facet
          if (length(facet) >= 1) {
            facetCol <- length(unique(df[,facet]))
            p1 = p1 + facet_wrap(as.formula(sprintf("~ %s" , paste(setdiff(facet, "SocpleType"), collapse="+"))), scales="free_x", ncol = facetCol)
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
                             panel.background = element_rect(fill = "white", colour = "Black",size = 0.5, linetype = "solid"),
                             aspect.ratio = 1/num.subgroup)

          plotlist[[length(plotlist)+1]] <- p1
        }
      }
    }else{
      # make a comnination for stat

      cbn <- combn(x = levels(df.na[,mvar]), m = 2)


      my_comparisons <- {}
      for(i in 1:ncol(cbn)){
        x <- cbn[,i]
        my_comparisons[[i]] <- x
      };my_comparisons
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


        if (statistics == "yes"| statistics == "YES"|statistics == "Yes"){
          if (parametric == "yes"| parametric == "YES"|parametric == "Yes"){
            if (nlevels(factor(df.na[,mvar])) > 2) {
              test <- aov(as.formula(sprintf("%s ~ %s", oc, mvar)), df.na)
              pval <- round(summary(test)[[1]][["Pr(>F)"]][1],4)
              test.name <- "ANOVA"
              testmethod <-  "t.test"
            } else {
              testmethod <-  "t.test"
              pval <- NULL
              test.name <- "Pairwise T-Test"
            }
          }else{
            if (nlevels(factor(df.na[,mvar])) > 2) {
              test <- kruskal.test(as.formula(sprintf("%s ~ %s", oc, mvar)), df.na)
              pval <- round(test$p.value, 4)
              test.name <- "KW"
              testmethod <- "wilcox.test"
            } else {
              testmethod <- "wilcox.test"
              pval <- NULL
              test.name <- "Pairwise Wilcoxon"
            }
          }
        }else{
          test.name<-NULL
          pval <- NULL
        }


        p1 <- ggplot(df.na, aes_string(x=mvar, y=oc))  + labs(y=oc, x=NULL) +
          theme(strip.background = element_blank()) +
          theme(text=element_text(size=8), axis.text.x=element_text(angle=xangle,hjust=1,vjust=0.5,size=8),
                plot.title=element_text(size=8))

        # paired plot type
        if (!is.null(paired)) {

          if(!is.null(mycols)){
            p1 <- p1 + scale_color_manual(values = mycols)
          }else{
            p1 <- p1
          }

          p1 = p1 + geom_boxplot(aes_string(colour=mvar),outlier.shape = NA,lwd=box.tickness)  + theme(legend.position="none")
          p1 = p1 + geom_point(aes_string(fill=mvar,group=paired),alpha = 0.8, size = dot.size, position = position_dodge(0.3), show.legend = F)   #scale_shape_manual(values = c(1, 16, 8, 0,15, 2,17,11, 10,12,3,4,5,6,7,8,9,13,14))
          p1 = p1 + geom_line(aes_string(group=paired), color="grey50", size=0.3,position = position_dodge(0.3))
          p1 = p1 + theme(legend.title = element_blank(), legend.position="bottom", legend.justification="left",legend.box.margin = ggplot2::margin(0,0,0,-1,"cm"))

        } else{
          # count or table for number of variable
          if (max(table(df.na[,mvar])) > 100){

            if(!is.null(mycols)){
              p1 <- p1 + scale_fill_manual(values = mycols)
            }else{
              p1 <- p1
            }
            p1 = p1 + geom_boxplot(aes_string(fill=mvar),outlier.shape = NULL,lwd=box.tickness)   + theme(legend.position="none")

          } else {

            if(!is.null(mycols)){
              p1 <- p1 + scale_color_manual(values = mycols)
            }else{
              p1 <- p1
            }

            p1 = p1 + geom_boxplot(aes_string(colour=mvar),outlier.shape = NA,lwd=box.tickness) + theme(legend.position="none")
            p1 = p1 + geom_jitter(aes_string(colour=mvar),shape=16, alpha = 0.8, size = dot.size, position=position_jitter(0.2))

          }
        }



        # control statistic on the plot
        if(is.null(paired)){
          if(is.null(test.name)){
            p1 <- p1
          } else if(test.name == "KW" | test.name == "ANOVA"){
            if(pval < cutoff){
              if (statistics == "yes"| statistics == "YES"|statistics == "Yes"){
                if (star == "no") {
                  p1 <- p1 + stat_compare_means(method= testmethod, label = "p.format", comparisons = my_comparisons, size = 2)
                } else if (star == "yes") {
                  p1 <- p1 + stat_compare_means(method= testmethod, label = "p.signif", comparisons = my_comparisons, hide.ns = TRUE, size = 3)
                }
              }else if(statistics == "no"| statistics == "NO"|statistics == "No"){
                p1 <- p1
              }
            }else {
              p1 <- p1
            }
          }else if(testmethod == "wilcox.test" | testmethod == "t.test"){
            if (statistics == "yes"| statistics == "YES"|statistics == "Yes"){
              if (star == "no") {
                p1 <- p1 + stat_compare_means(method= testmethod, label = "p.format", comparisons = my_comparisons, size = 2)
              }  else if (star == "yes") {
                p1 <- p1 + stat_compare_means(method= testmethod, label = "p.signif", comparisons = my_comparisons, hide.ns = TRUE, size = 3)
              }
            }else if(statistics == "no"| statistics == "NO"|statistics == "No"){
              p1 <- p1
            }
          }
        }else{

          if(is.null(test.name)){
            p1 <- p1
          } else if(test.name == "KW" | test.name == "ANOVA"){
            if(pval < cutoff){
              if (statistics == "yes"| statistics == "YES"|statistics == "Yes"){
                if (star == "no") {
                  p1 <- p1 + stat_compare_means(method= testmethod, label = "p.format", comparisons = my_comparisons, size = 2, paired = TRUE)
                } else if (star == "yes") {
                  p1 <- p1 + stat_compare_means(method= testmethod, label = "p.signif", comparisons = my_comparisons, hide.ns = TRUE,
                                                size = 3,paired = TRUE)
                }
              }else if(statistics == "no"| statistics == "NO"|statistics == "No"){
                p1 <- p1
              }
            }else {
              p1 <- p1
            }
          }else if(testmethod == "wilcox.test" | testmethod == "t.test"){
            if (statistics == "yes"| statistics == "YES"|statistics == "Yes"){
              if (star == "no") {

                if (data.frame(table(df.na[,mvar]))$Freq[1] ==  data.frame(table(df.na[,mvar]))$Freq[2]){
                  p1 <- p1 + stat_compare_means(method= testmethod, label = "p.format", comparisons = my_comparisons, size = 2,paired = TRUE)
                }else{
                  test.name <- paste(test.name, "\n","(not fully paired) ", sep = "")
                  p1 <- p1 + stat_compare_means(method= testmethod, label = "p.format", comparisons = my_comparisons, size = 2)
                }
              }  else if (star == "yes") {

                if (data.frame(table(df.na[,mvar]))$Freq[1] ==  data.frame(table(df.na[,mvar]))$Freq[2]){
                  p1 <- p1 + stat_compare_means(method= testmethod, label = "p.signif", comparisons = my_comparisons, size = 2,paired = TRUE)
                }else{
                  p1 <- p1 + stat_compare_means(method= testmethod, label = "p.signif", comparisons = my_comparisons, size = 2)
                }

              }
            }else if(statistics == "no"| statistics == "NO"|statistics == "No"){
              p1 <- p1
            }
          }
        }



        # Close an image
        if (!is.null(title)) {
          p1 <- p1 + ggtitle(sprintf("%s%s%s%s", title,
                                     ifelse(is.null(test.name), "", paste("\n",test.name, " ", sep = "")),
                                     ifelse(is.null(pval), "", paste("p=", " ", sep = "")),
                                     ifelse(is.null(pval), "", paste(pval, " ", sep = "")), sep=""))
        } else{
          p1 <- p1 + ggtitle(sprintf("%s%s%s%s", mvar,
                                     ifelse(is.null(test.name), "", paste("\n",test.name, " ", sep = "")),
                                     ifelse(is.null(pval), "", paste("p=", " ", sep = "")),
                                     ifelse(is.null(pval), "", paste(pval, " ", sep = "")), sep=""))
        }


        # y axis limit
        if(!is.null(ylim)){
          if(oc == "Chao1"){
            p1 = p1
          }else{
            p1 = p1 + ylim(ylim[1] , ylim[2])
          }
        }
        # facet
        if (length(facet) >= 1) {
          facetCol <- length(unique(df[,facet]))
          p1 = p1 + facet_wrap(as.formula(sprintf("~ %s" , paste(setdiff(facet, "SocpleType"), collapse="+"))), scales="free_x", ncol = facetCol)
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
                           panel.background = element_rect(fill = "white", colour = "Black",size = 0.5, linetype = "solid"),
                           aspect.ratio = 1/num.subgroup)

        plotlist[[length(plotlist)+1]] <- p1

      }
    }
  }
  multiplot(plotlist=plotlist, cols=plotCols, rows=plotRows)
}

