#' Beta Diversity Analysis
#'
#' This function performs beta diversity analysis on a given phyloseq object.
#' It supports various diversity metrics, group combinations, and visualization options such as PCoA plots with ellipses.
#'
#' @param psIN Phyloseq object containing the OTU/ASV counts and associated sample data.
#' @param cate.vars Vector of column names in the sample data of 'psIN' representing categorical variables.
#' @param project Name or identifier for the project or dataset being analyzed.
#' @param orders Optional ordering for the levels of the categorical variables.
#' @param distance_metrics Vector of distance metrics to be used for beta diversity analysis.
#' @param cate.conf Optional vector of column names in the sample data of 'psIN' representing categorical confounding variables.
#' @param plot Type of plot to generate, default is "PCoA".
#' @param ellipse Option to add ellipses to the plot ("yes" or "no").
#' @param statistics Option to perform statistical analysis ("yes" or "no").
#' @param mycols Color palette for the plot.
#' @param paired Optional parameter for paired analysis.
#' @param combination Numeric value for group combination analysis.
#' @param shapes Optional vector for shape aesthetics in the plot.
#' @param ID Optional parameter for identifying points on the plot.
#' @param facet Optional parameter for facetting the plot.
#' @param name Optional name for the analysis output.
#' @param addnumber Option to add sample size numbers to the plot categories.
#'
#' @return A list of ggplot objects, each representing a beta diversity plot for the specified distance metric and categorical variable.
#'
#' @examples
#' # Assuming 'ps' is a phyloseq object with appropriate data
#' results <- Go_bdiv_markdown(ps, cate.vars = c("Group"),
#'                             project = "my_project",
#'                             orders = NULL,
#'                             distance_metrics = c("bray", "unifrac"),
#'                             plot = "PCoA",
#'                             ellipse = "yes")
#'
#' @export
#' @importFrom phyloseq phyloseq sample_data
#' @importFrom ggplot2 ggplot aes_string geom_point geom_text_repel stat_ellipse theme_bw
#' @importFrom gridExtra grobTree textGrob
#' @importFrom grid gpar

Go_bdiv_markdown <- function(psIN, cate.vars, project, orders, distance_metrics,
                    cate.conf=NULL,
                    plot="PCoA",
                    ellipse="yes",
                    statistics = "yes",
                    mycols=NULL,
                    paired = NULL,
                    combination=NULL,
                    shapes = NULL,
                    ID = NULL,
                    facet=NULL,
                    name=NULL,
                    addnumber=TRUE){

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

      for(i in 1:length(group_comparisons)){
        #print(group_comparisons[i])
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


          ## fix factor  and  numeric
          mapping.sel.na.rem[,mvar] <- factor(mapping.sel.na.rem[,mvar])


          ord_meths = plot # c("DCA", "CCA", "RDA", "DPCoA", "NMDS","PCoA")
          plist = llply(as.list(ord_meths), function(i, psIN.cbn.na, distance_metric){
            ordi = ordinate(psIN.cbn.na, method=i, distance=distance_metric)
            plot_ordination(psIN.cbn.na, ordi, type = "samples", color= mvar)
          }, psIN.cbn.na, distance_metric)



          names(plist) <- ord_meths

          pdataframe = ldply(plist, function(x){
            df = x$data[, 1:2]
            colnames(df) = c("Axis_1", "Axis_2")
            return(cbind(df, x$data))
          })
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
          p = ggplot(pdataframe, aes_string("Axis_1", "Axis_2", color=mvar))


          if (!is.null(shapes)) {

            pdataframe[,shapes] <- factor(pdataframe[,shapes], levels = orders)
            p = p +  geom_point(aes_string(shape=shapes), size=0.8, alpha = 1) + scale_shape_manual(values = c(1, 16, 8, 0,15, 2,17,11, 10,12,3,4,5,6,7,8,9,13,14))

          }else{
            p = p + geom_point(size=0.8, alpha = 1)+ ggtitle(sprintf("%s (%s)",mvar,distance_metric))
          }

          p = p + ggtitle(sprintf("%s (%s)",mvar,distance_metric))
          p = p + facet_wrap(~ method, scales="free") + theme_bw() + theme(strip.background = element_blank())# open(1), cross(10), closed(2)
          p = p + theme(legend.position = "bottom",
                        legend.title = element_blank(),
                        legend.text=element_text(size=11),
                        legend.justification="left",
                        legend.box = "vertical",
                        legend.box.margin = ggplot2::margin(0,0,0,-1,"cm"))


          if(!is.null(mycols)){
            p <- p + scale_color_manual(values = mycols)
          }else{
            p <- p
          }

          # ID variation
          if (!is.null(ID)) {
            p = p + geom_text_repel(aes_string(label = ID), size = 2)
          } else {
            p = p
          }

          # ellipse variation
          if (ellipse == "yes" | ellipse == "Yes" ) {
            p = p + stat_ellipse(type = "norm", linetype = 2)
          } else if (ellipse == "no" | ellipse == "No" ){
            p = p
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
          if (statistics == "yes"| statistics == "YES"|statistics == "Yes"){
            set.seed(1)
            distance <- Go_dist_markdown(psIN = psIN.cbn.na, project = project, cate.vars = mvar, distance_metrics = distance_metric)

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
              #print(form)
            }else{
              form <- as.formula(sprintf("x1 ~ %s", mvar))
              #print(form)
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



        ## fix factor  and  numeric
        mapping.sel.na.rem[,mvar] <- factor(mapping.sel.na.rem[,mvar])



        ord_meths= plot # c("DCA", "CCA", "RDA", "DPCoA", "NMDS","PCoA")
        plist = llply(as.list(ord_meths), function(i, psIN.na, distance_metric){
          ordi = ordinate(psIN.na, method=i, distance=distance_metric)
          plot_ordination(psIN.na, ordi, type = "samples", color= mvar)
        }, psIN.na, distance_metric)

        names(plist) <- ord_meths

        pdataframe = ldply(plist, function(x){
          df = x$data[, 1:2]
          colnames(df) = c("Axis_1", "Axis_2")
          return(cbind(df, x$data))
        })
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
        p = ggplot(pdataframe, aes_string("Axis_1", "Axis_2", color=mvar))


        if (!is.null(shapes)) {
          pdataframe[,shapes] <- factor(pdataframe[,shapes], levels = orders)
          p = p +  geom_point(aes_string(shape=shapes), size=0.8, alpha = 1) + scale_shape_manual(values = c(1, 16, 8, 0,15, 2,17,11, 10,12,3,4,5,6,7,8,9,13,14))

        }else{
          p = p + geom_point(size=0.8, alpha = 1)+ ggtitle(sprintf("%s (%s)",mvar,distance_metric))
        }

        p = p + ggtitle(sprintf("%s (%s)",mvar,distance_metric))
        p = p + facet_wrap(~ method, scales="free") + theme_bw() + theme(strip.background = element_blank())# open(1), cross(10), closed(2)
        p = p + theme(legend.position = "bottom",
                      legend.title = element_blank(),
                      legend.text=element_text(size=11),
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
        if (ellipse == "yes" | ellipse == "Yes" ) {
          p <- p + stat_ellipse(type = "norm", linetype = 2)
        } else if (ellipse == "no" | ellipse == "No" ){
          p <- p
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

        if (statistics == "yes"| statistics == "YES"|statistics == "Yes"){
          set.seed(1)
          distance <- Go_dist_markdown(psIN = psIN.na, project = project, name = NULL, cate.vars = mvar,distance_metrics = distance_metric)

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
            #print(form)
          }else{
            form <- as.formula(sprintf("x1 ~ %s", mvar))
            #print(form)
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
}

