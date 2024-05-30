
#' Generate Biodiversity Plots for Microbiome Data
#'
#' @param psIN Phyloseq object containing microbiome data.
#' @param cate.vars Categorical variables to be used for grouping in the plot.
#' @param project Name of the project for labeling.
#' @param orders Ordering of factors in the plot.
#' @param distance_metrics Vector of distance metrics to be used.
#' @param cate.conf Additional confounder variables, if any.
#' @param plot Type of plot to generate, e.g., "PCoA".
#' @param ellipse Logical indicating whether to include ellipses.
#' @param statistics Logical indicating whether to include statistical tests.
#' @param mycols Custom color palette.
#' @param paired Variable for paired analysis, if applicable.
#' @param combination Number of groups to combine for comparison.
#' @param shapes Shape of points in the plot.
#' @param ID Identifier for labeling points.
#' @param facet Faceting variable for the plot.
#' @param name Optional name for the output plot.
#' @param addnumber Logical indicating whether to add sample numbers to groups.
#' @param height Height of the plot.
#' @param width Width of the plot.
#'
#' @details
#' This function creates biodiversity plots such as PCoA for microbiome data, facilitating comparisons across different conditions or factors. The function supports various distance metrics and provides options for customization.
#'
#' @return
#' A PDF file containing the biodiversity plot(s).
#'
#' @examples
#' Go_bdiv(psIN = ps_object,
#'         cate.vars = c("Condition", "Treatment"),
#'         project = "MyMicrobiomeStudy",
#'         orders = c("Condition1", "Condition2"),
#'         distance_metrics = c("bray", "unifrac"),
#'         cate.conf = "AgeGroup",
#'         plot = "PCoA",
#'         ellipse = TRUE, or group names
#'         statistics = TRUE,
#'         mycols = c("blue", "red"),
#'         paired = "PatientID",
#'         combination = 2,
#'         shapes = 16,
#'         ID = "SampleID",
#'         facet = "Group",
#'         name = "BiodivPlot",
#'         addnumber = TRUE,
#'         height = 6,
#'         width = 8)
#'
#' @export

Go_bdiv <- function(psIN, cate.vars, project, orders, distance_metrics,
                    cate.conf=NULL,
                    plot="PCoA",
                    ellipse=TRUE,
                    statistics = TRUE,
                    mycols=NULL,
                    paired = NULL,
                    combination=NULL,
                    shapes = NULL,
                    ID = NULL,
                    facet=NULL,
                    name=NULL,
                    addnumber=TRUE,
                    height, width){

  if(!is.null(dev.list())) dev.off()

  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out_path)) dir.create(out_path)

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


  pdf(sprintf("%s/ordi.%s.%s%s%s%s%s%s%s.pdf", out_path,
              project,
              ifelse(is.null(facet), "", paste(facet, ".", sep = "")),
              ifelse(is.null(combination), "", paste("(cbn=",combination, ").", sep = "")),
              ifelse(is.null(cate.conf), "", paste("with_confounder", ".", sep = "")),
              ifelse(is.null(paired), "", paste("(paired=",paired, ").", sep = "")),
              ifelse(is.null(name), "", paste(name, ".", sep = "")),
              ifelse(ellipse == FALSE, "ellipse_FALSE.",
                     ifelse(ellipse == TRUE, "", paste("ellipse_", ellipse, ".", sep = ""))),
              format(Sys.Date(), "%y%m%d")), height = height, width = width)

  plotlist <- list()
  for (mvar in cate.vars) {
    mapping <- data.frame(sample_data(psIN))
    mapping[,mvar] <- factor(mapping[,mvar])

    sample_data(psIN) <- mapping


    if (length(facet) >= 1){
      if (facet == mvar){
        next
      }
    } else {}

    if (length(shapes) >= 1){
      if (shapes == mvar){
        next
      }
    } else {}

    #------------------------------#
    # for group combination or not #
    #------------------------------#

    if (!is.null(combination)){
      mapping[,mvar] <- factor(mapping[,mvar], levels = intersect(orders, mapping[,mvar]))
      group.cbn <- combn(x = levels(mapping[,mvar]), m = combination)

      #print(count(group.cbn))

      group_comparisons <- {}
      for(i in 1:ncol(group.cbn)){
        x <- group.cbn[,i]
        group_comparisons[[i]] <- x
      };group_comparisons

      print(1)
      for(i in 1:length(group_comparisons)){
        print(group_comparisons[i])
        group.combination <- unlist(group_comparisons[i]);group.combination

        if(combination ==2){
          basline <- group.combination[1]
          smvar <- group.combination[2]
          mapping.cbn <- subset(mapping, mapping[,mvar] %in% c(basline,smvar))
        } else if(combination ==3){
          basline <- group.combination[1]
          smvar1 <- group.combination[2]
          smvar2 <- group.combination[3]
          mapping.cbn <- subset(mapping, mapping[,mvar] %in% c(basline,smvar1, smvar2))
        }else if(combination ==4){
          basline <- group.combination[1]
          smvar1 <- group.combination[2]
          smvar2 <- group.combination[3]
          smvar3 <- group.combination[4]
          mapping.cbn <- subset(mapping, mapping[,mvar] %in% c(basline,smvar1, smvar2,smvar3))
        }else{
          print("combination should be 2, 3, and 4 only.")
          break
        }

        psIN.cbn <- psIN
        sample_data(psIN.cbn) <- mapping.cbn
        for(distance_metric in distance_metrics){
          # remove na
          mapping.sel <- data.frame(sample_data(psIN.cbn))
          mapping.sel[mapping.sel==""] <- "NA"
          mapping.sel.na <- mapping.sel[!is.na(mapping.sel[,mvar]), ]
          na.count <- length(mapping.sel.na)
          psIN.cbn.na <- prune_samples(rownames(mapping.sel[!is.na(mapping.sel[,mvar]), ]), psIN)
          mapping.sel.na.rem <- data.frame(sample_data(psIN.cbn.na ))


          if (!is.null(facet)) {
            print(sprintf("##-- %s-%s (total without NA: %s/%s) --##",
                          facet,mvar, dim(mapping.sel.na.rem)[1], dim(mapping.sel)[1]))
          } else{
            print(sprintf("##-- %s (total without NA: %s/%s) --##",
                          mvar, dim(mapping.sel.na.rem)[1], dim(mapping.sel)[1]))
          }



          ## fix factor  and  numeric
          mapping.sel.na.rem[,mvar] <- factor(mapping.sel.na.rem[,mvar])


          ord_meths = plot # c("DCA", "CCA", "RDA", "DPCoA", "NMDS","PCoA")
          # Execute ordination and collect necessary outputs
          plist = plyr::llply(as.list(ord_meths), function(i, psIN.na, distance_metric){
            ordi = ordinate(psIN.na, method=i, distance=distance_metric)
            df = as.data.frame(ordi$vectors[, 1:2])  # Assuming 'vectors' has the axes coordinates
            colnames(df) = c("Axis_1", "Axis_2")

            # Calculate percentage of variance explained for the first two axes
            if ("Eigenvalues" %in% names(ordi$values)) {
              var_explained = ordi$values$Eigenvalues / sum(ordi$values$Eigenvalues) * 100
              df$Axis1_Percent = var_explained[1]  # Percentage for Axis 1
              df$Axis2_Percent = var_explained[2]  # Percentage for Axis 2
            }

            # Merge with sample metadata
            metadata = as.data.frame(sample_data(psIN.na))
            return(cbind(df, metadata))  # Combine ordination data with sample metadata
          }, psIN.cbn.na, distance_metric)

          # Name the list elements according to the ordination methods
          names(plist) <- ord_meths

          # Convert the list to a dataframe
          pdataframe = plyr::ldply(plist, identity)


          names(pdataframe)[1] = "method"

          pdataframe[,facet] <- factor(pdataframe[,facet], levels = orders)

          pdataframe[,mvar] <- factor(pdataframe[,mvar], levels = orders)


           # Add number of samples in the group
          if(addnumber==TRUE){
              renamed_levels <- as.character(levels(pdataframe[,mvar]));renamed_levels
              oldNames <- unique(pdataframe[,mvar]);oldNames
          if (length(renamed_levels) == 0) {
            renamed_levels <- oldNames
          }
          for (name in oldNames) {
            total <- length(which(pdataframe[,mvar] == name));total
            new_n <- paste(name, " (n=", total, ")", sep="");new_n
            levels(pdataframe[[mvar]])[levels(pdataframe[[mvar]])== name] <- new_n
            renamed_levels <- replace(renamed_levels, renamed_levels == name, new_n);renamed_levels
          }
          }else{
            pdataframe <- pdataframe
          }

          # Plots
          axis1_percent_avg <- mean(pdataframe$Axis1_Percent, na.rm = TRUE)
          axis2_percent_avg <- mean(pdataframe$Axis2_Percent, na.rm = TRUE)

          if (!is.null(shapes)) {
            pdataframe[,shapes] <- factor(pdataframe[,shapes], levels = orders)
            p = ggplot(pdataframe, aes_string(x = "Axis_1", y = "Axis_2", color = mvar))
            p = p + geom_point(aes_string(shape=shapes), size=0.9, alpha = 1) +  # Add points to the plot
              scale_shape_manual(values = c(1, 16, 8, 0,15, 2,17,11, 10,12,3,4,5,6,7,8,9,13,14))

          }else{
            p = ggplot(pdataframe, aes_string(x = "Axis_1", y = "Axis_2", color = mvar))
            p = p + geom_point(size=0.9, alpha = 1)
          }


          p = p + labs(
              x = paste("Axis 1 (", sprintf("%.2f", axis1_percent_avg), "%)", sep = ""),
              y = paste("Axis 2 (", sprintf("%.2f", axis2_percent_avg), "%)", sep = "")
            )# Add points to the plot
          p = p + ggtitle(sprintf("%s (%s)",mvar,distance_metric))
          p = p + facet_wrap(~ method, scales="free") + theme_bw() + theme(strip.background = element_blank())# open(1), cross(10), closed(2)
          p = p + theme(legend.position = "bottom",
                        legend.title = element_blank(),
                        legend.justification="left",
                        legend.box = "vertical",
                        legend.box.margin = ggplot2::margin(0,0,0,-1,"cm"))


          if(!is.null(mycols)){
            p <- p + scale_color_manual(values = mycols)
          }else{
            p <- plist[[1]]
          }

          # ID variation
          if (!is.null(ID)) {
            p = p + geom_text_repel(aes_string(label = ID), size = 2)
          } else {
            p = p
          }

          # ellipse variation
          if (!is.null(ellipse) && ellipse != TRUE) {
            p <- p + stat_ellipse(aes_string(group = ellipse, color = ellipse), type = "norm", linetype = 2)
          } else if (ellipse == TRUE) {
            p <- p + stat_ellipse(type = "norm", linetype = 2)
          }




          if (!is.null(facet)) {
            ncol <- length(unique(mapping.sel.na.rem[,facet]))
            p = p + facet_wrap(as.formula(sprintf("~ %s", facet)), scales="free_x", ncol = ncol)
          }
          else {
            p = p
          }

          if(!is.null(paired)){
            p <- p + geom_line(aes_string(group = paired),color="grey", size=0.2,  arrow = arrow(type = "closed",
                                                                                      length=unit(0.025, "inches")))
          }
          #===================================#
          # Add permanova for two combination #
          #===================================#
          if (statistics){
            set.seed(1)
            distance <- Go_dist(psIN = psIN.cbn.na, project = project, cate.vars = mvar, distance_metrics = distance_metric)

            x <- as.dist(distance[[distance_metric]])
            factors <-  mapping.sel.na.rem[,mvar]

            R2 <- c()
            p.value <- c()
            F.Model <- c()
            pairs <- c()
            SumsOfSqs <- c()
            Df <- c()

            x1=as.matrix(x)[factors %in% unique(factors), factors %in% unique(factors)]

            # run
            map.pair <- subset(mapping.sel.na.rem, mapping.sel.na.rem[,mvar] %in% unique(factors))

            # count to table

            if (!is.null(cate.conf)) {
              for(conf in cate.conf){
                map.pair[,conf] <- factor(map.pair[,conf])
              }
              form <- as.formula(sprintf("x1 ~ %s + %s", mvar, paste(setdiff(cate.conf, "SampleType"), collapse="+")))
              print(form)
            }else{
              form <- as.formula(sprintf("x1 ~ %s", mvar))
              print(form)
            }

            ad <- adonis2(form, data = map.pair, permutations=999, by="terms")# "terms"  "margin" NULL

            Df <- c(Df,ad[1,1])
            SumsOfSqs <- c(SumsOfSqs, ad[1,2])
            R2 <- round(c(R2,ad[1,3]), digits=3)
            F.Model <- c(F.Model,ad[1,4]);
            p.value <- c(p.value,ad[1,5])

            pairw.res <- data.frame(Df,SumsOfSqs,R2,F.Model,p.value)

            class(pairw.res) <- c("pwadonis", "data.frame")
            # end adonis end
            tmp <- as.data.frame(pairw.res)
            tmp$distance_metric <- distance_metric
            tmp$padj <- p.adjust(tmp$p.value, method="bonferroni")

            grob <- grobTree(textGrob(paste(distance_metric, "\nR2=",R2,"\nPERMANOVA p=",tmp$padj,sep=""), x=0.01,  y=0.15, hjust=0,
                                      gp=gpar(fontsize=8))) #, fontface="italic"


            if (table(map.pair[,mvar])[1] <=2 | table(map.pair[,mvar])[2] <=2){
              p=p
            }else{
              p = p + annotation_custom(grob)
            }
          }else{
            p=p
          }

          #plotlist[[length(plotlist)+1]] <- p
          print(p)
        }
      }
    }  else{
      for(distance_metric in distance_metrics){
        # remove na
        mapping.sel <- data.frame(sample_data(psIN))
        #mapping.sel[mapping.sel==""] <- "NA"
        mapping.sel.na <- mapping.sel[!is.na(mapping.sel[,mvar]), ]
        na.count <- length(mapping.sel.na)
        psIN.na <- prune_samples(rownames(mapping.sel[!is.na(mapping.sel[,mvar]), ]), psIN)
        mapping.sel.na.rem <- data.frame(sample_data(psIN.na ))


        if (!is.null(facet)) {
          print(sprintf("##-- %s-%s (total without NA: %s/%s) --##",
                        facet,mvar, dim(mapping.sel.na.rem)[1], dim(mapping.sel)[1]))
        } else{
          print(sprintf("##-- %s (total without NA: %s/%s) --##",
                        mvar, dim(mapping.sel.na.rem)[1], dim(mapping.sel)[1]))
        }



        ## fix factor  and  numeric
        mapping.sel.na.rem[,mvar] <- factor(mapping.sel.na.rem[,mvar])



        ord_meths= plot # c("DCA", "CCA", "RDA", "DPCoA", "NMDS","PCoA")
        # Execute ordination and collect necessary outputs
        plist = plyr::llply(as.list(ord_meths), function(i, psIN.na, distance_metric){
          ordi = ordinate(psIN.na, method=i, distance=distance_metric)
          df = as.data.frame(ordi$vectors[, 1:2])  # Assuming 'vectors' has the axes coordinates
          colnames(df) = c("Axis_1", "Axis_2")

          # Calculate percentage of variance explained for the first two axes
          if ("Eigenvalues" %in% names(ordi$values)) {
            var_explained = ordi$values$Eigenvalues / sum(ordi$values$Eigenvalues) * 100
            df$Axis1_Percent = var_explained[1]  # Percentage for Axis 1
            df$Axis2_Percent = var_explained[2]  # Percentage for Axis 2
          }

          # Merge with sample metadata
          metadata = as.data.frame(sample_data(psIN.na))
          return(cbind(df, metadata))  # Combine ordination data with sample metadata
        }, psIN.na, distance_metric)

        # Name the list elements according to the ordination methods
        names(plist) <- ord_meths

        # Convert the list to a dataframe
        pdataframe = plyr::ldply(plist, identity)



        names(pdataframe)[1] = "method"

        pdataframe[,facet] <- factor(pdataframe[,facet], levels = orders)

        pdataframe[,mvar] <- factor(pdataframe[,mvar], levels = orders)

           # Add number of samples in the group
          if(addnumber==TRUE){
              renamed_levels <- as.character(levels(pdataframe[,mvar]));renamed_levels
              oldNames <- unique(pdataframe[,mvar]);oldNames
          if (length(renamed_levels) == 0) {
            renamed_levels <- oldNames
          }
          for (name in oldNames) {
            total <- length(which(pdataframe[,mvar] == name));total
            new_n <- paste(name, " (n=", total, ")", sep="");new_n
            levels(pdataframe[[mvar]])[levels(pdataframe[[mvar]])== name] <- new_n
            renamed_levels <- replace(renamed_levels, renamed_levels == name, new_n);renamed_levels
          }
          }else{
            pdataframe <- pdataframe
          }


        # Plots
        axis1_percent_avg <- mean(pdataframe$Axis1_Percent, na.rm = TRUE)
        axis2_percent_avg <- mean(pdataframe$Axis2_Percent, na.rm = TRUE)



        if (!is.null(shapes)) {
          pdataframe[,shapes] <- factor(pdataframe[,shapes], levels = orders)
          p = ggplot(pdataframe, aes_string(x = "Axis_1", y = "Axis_2", color = mvar))
          p = p + geom_point(aes_string(shape=shapes), size=0.9, alpha = 1) +  # Add points to the plot
            scale_shape_manual(values = c(1, 16, 8, 0,15, 2,17,11, 10,12,3,4,5,6,7,8,9,13,14))

        }else{
          p = ggplot(pdataframe, aes_string(x = "Axis_1", y = "Axis_2", color = mvar))
          p = p + geom_point(size=0.9, alpha = 1)
        }


        p = p + labs(
          x = paste("Axis 1 (", sprintf("%.2f", axis1_percent_avg), "%)", sep = ""),
          y = paste("Axis 2 (", sprintf("%.2f", axis2_percent_avg), "%)", sep = "")
        )# Add points to the plot

        p = p + ggtitle(sprintf("%s (%s)",mvar,distance_metric))
        p = p + facet_wrap(~ method, scales="free") + theme_bw() + theme(strip.background = element_blank())# open(1), cross(10), closed(2)
        p = p + theme(legend.position = "bottom",
                      legend.title = element_blank(),
                      legend.justification="left",
                      legend.box = "vertical",
                      legend.box.margin = ggplot2::margin(0,0,0,-1,"cm"),
                      plot.title=element_text(size=8,face="bold"))

        if(!is.null(mycols)){
          p <- p + scale_color_manual(values = mycols)
        }else{
          p <- p
        }

        # ID variation
        if (!is.null(ID)) {
          p <- p + geom_text_repel(aes_string(label = ID), size = 2)
        } else {
          p <- p
        }

        # ellipse variation
        if (!is.null(ellipse) && ellipse != TRUE) {
          p <- p + stat_ellipse(aes_string(group = ellipse, color = ellipse), type = "norm", linetype = 2)
        } else if (ellipse == TRUE) {
          p <- p + stat_ellipse(type = "norm", linetype = 2)
        }

        if (!is.null(facet)) {
          ncol <- length(unique(mapping.sel.na.rem[,facet]))
          p <- p + facet_wrap(as.formula(sprintf("~ %s", facet)), scales="free_x", ncol = ncol)
        }
        else {
          p <- p
        }


        if(!is.null(paired)){
          p <- p + geom_line(aes_string(group = paired),color="grey", size = 0.2,arrow = arrow(type = "closed",
                                                                                    length=unit(0.025, "inches")))
        }

        #===================================#
        # Add permanova                     #
        #===================================#

        if (statistics){
          set.seed(1)
          distance <- Go_dist(psIN = psIN.na, project = project, name = NULL, cate.vars = mvar,distance_metrics = distance_metric)

          x <- as.dist(distance[[distance_metric]])
          factors <-  mapping.sel.na.rem[,mvar]

          R2 <- c()
          p.value <- c()
          F.Model <- c()
          pairs <- c()
          SumsOfSqs <- c()
          Df <- c()

          x1=as.matrix(x)[factors %in% unique(factors), factors %in% unique(factors)]

          # run
          map.pair <- subset(mapping.sel.na.rem, mapping.sel.na.rem[,mvar] %in% unique(factors))

          # count to table

          if (!is.null(cate.conf)) {
            for(conf in cate.conf){
              map.pair[,conf] <- factor(map.pair[,conf])
            }
            form <- as.formula(sprintf("x1 ~ %s + %s", mvar, paste(setdiff(cate.conf, "SampleType"), collapse="+")))
            print(form)
          }else{
            form <- as.formula(sprintf("x1 ~ %s", mvar))
            print(form)
          }

          ad <- adonis2(form, data = map.pair, permutations=999, by="terms")# "terms"  "margin" NULL

          Df <- c(Df,ad[1,1])
          SumsOfSqs <- c(SumsOfSqs, ad[1,2])
          R2 <- round(c(R2,ad[1,3]), digits=3)
          F.Model <- c(F.Model,ad[1,4]);
          p.value <- c(p.value,ad[1,5])

          pairw.res <- data.frame(Df,SumsOfSqs,R2,F.Model,p.value)

          class(pairw.res) <- c("pwadonis", "data.frame")
          # end adonis end
          tmp <- as.data.frame(pairw.res)
          tmp$distance_metric <- distance_metric
          tmp$padj <- p.adjust(tmp$p.value, method="bonferroni")

          grob <- grobTree(textGrob(paste(distance_metric, "\nR2=",R2,"\nPERMANOVA p=",tmp$padj,sep=""), x=0.01,  y=0.15, hjust=0,
                                    gp=gpar(fontsize=8))) #, fontface="italic"

          if (table(map.pair[,mvar])[1] <=2 | table(map.pair[,mvar])[2] <=2){
            p=p
          }else{
            p = p + annotation_custom(grob)
          }
        }else{
          p=p
        }
        #plotlist[[length(plotlist)+1]] <- p

        p  <- p  + theme(panel.grid = element_blank(),
                         legend.key = element_blank(), # remove square
                         panel.background = element_rect(fill = "white", colour = "Black",size = 0.7, linetype = "solid"),
                         aspect.ratio = 1) +
          geom_vline(xintercept = 0, size = 0.1) + geom_hline(yintercept = 0, size = 0.1)
        print(p)
      }
    }
  }
  dev.off()
}

