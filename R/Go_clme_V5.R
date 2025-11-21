
#' Conduct Conditional Linear Mixed Effects (CLME) Analysis on Phyloseq Data
#'
#' This function performs a conditional linear mixed effects analysis on data from a Phyloseq object.
#' It supports various customization options including grouping, node specification, and plotting parameters.
#'
#' @param psIN Phyloseq object containing the data for analysis.
#' @param cate.vars Categorical variables used in the analysis.
#' @param project Name of the project or analysis.
#' @param paired Logical value indicating if the data is paired.
#' @param mycols Color specifications for the plot.
#' @param addnumber Logical value to add the number of samples in each group.
#' @param standardsize Logical value to standardize plot size.
#' @param node Node parameter for umbrella pattern in CLME.
#' @param decreasing Logical value indicating the decreasing trend in CLME.
#' @param timepoint Time point variable in the data.
#' @param ID Identifier for individual samples.
#' @param orders Order of levels in the categorical variables.
#' @param xangle Angle of x-axis labels in the plot.
#' @param name Optional name for the analysis.
#' @param height Height of the output plot.
#' @param width Width of the output plot.
#' @param plotCols Number of columns for multiple plots layout.
#' @param plotRows Number of rows for multiple plots layout.
#'
#' @return Generates CLME plots as PDF files and saves them in a specified directory.
#'
#' @details
#' The function applies a CLME model to the Phyloseq data, considering various parameters and
#' constraints. It can handle paired data, apply node constraints, and consider decreasing trends.
#' The function generates plots for each categorical variable and subgroup within the data.
#'
#' @examples
#' # psIN is a Phyloseq object
#' # Example usage:
#' Go_clme(psIN = psIN,
#'         cate.vars = c("Treatment", "Condition"),
#'         project = "MyProject",
#'         paired = TRUE,
#'         node = 10,
#'         decreasing = TRUE,
#'         timepoint = "Day",
#'         ID = "PatientID",
#'         orders = c("Day1", "Day2", "Day3"),
#'         xangle = 45,
#'         name = "Analysis1",
#'         height = 10,
#'         width = 10,
#'         plotCols = 2,
#'         plotRows = 2)
#'
#' @export



Go_clme <- function(psIN, cate.vars, project, paired,
                    mycols=NULL,
                    addnumber=TRUE,
                    standardsize=TRUE,
                    node,
                    decreasing,
                    timepoint,
                    ID,
                    orders,
                    xangle,
                    name,
                    height,
                    width,
                    plotCols,
                    plotRows){

  if(!is.null(dev.list())) dev.off()

  alpha_metrics = c("Chao1","Shannon")

  # Descriptions 분석 하고자 하는 variation에 subgroup
  # paired 환자나 같은 사람 ID
  # node 전반적인 패턴을 보고 가장 높은 time point 에 node를 설정
  # decreasing 패턴에 증가 하는지 감소 하는지 판단 하고 decreazing = true and false 를 판단, mean and median 으로 판단

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

  pdf(sprintf("%s/clme.%s.%s.%s.%s.%s.pdf", out_path,
              project,
              ifelse(is.null(name), "", paste(name, ".", sep = "")),
              node,
              decreasing,
              format(Sys.Date(), "%y%m%d")), height = height, width = width)

  # adiv
  adiv <- estimate_richness(psIN, measures=alpha_metrics)
  mapping <-data.frame(sample_data(psIN))
  rownames(adiv) <- gsub("^X", "", rownames(adiv))
  adiv$SampleID <- rownames(adiv)
  rownames(adiv) <- rownames(mapping)
  adiv <- merge(adiv, mapping, by="row.names"); rownames(adiv) <- adiv$SampleID

  if (length(orders) >= 1) {
    adiv[,timepoint] <- factor(adiv[,timepoint], levels = orders)
  }


  # clme
  cons <- list(order = "umbrella" , node=node, decreasing = decreasing)
  # 전반적인 패턴을 보고 가장 높은 time point 에 node를 설정
  # 패턴에 증가 하는지 감소 하는지 판단 하고 decreazing = true and false 를 판단, mean and median 으로 판단

  print(cons)

  plotlist <- list()
  for (mvar in cate.vars) {
    print(mvar)

    if (length(unique(adiv[,mvar])) < 2){
      next
    }



    # Na 제거

    adiv[,mvar] <- data.frame(adiv[,mvar]);adiv[,mvar]
    adiv[,mvar][adiv[,mvar]==""] <- "NA";adiv[,mvar]
    adiv[,mvar]<- as.factor(adiv[,mvar]);adiv[,mvar]


    # adiv.na <- adiv[!(is.na(adiv[,mvar])), ];adiv.na[,mvar] 틀린건 없는 거 같은데 지워지지 않는다.
    adiv.na <- subset(adiv, adiv[,mvar] != "NA");adiv.na[,mvar]  # subset 를 사용한 NA 삭제
    adiv <- adiv.na



    # Add number of samples in the group
    if(addnumber==TRUE){
    renamed_levels <- as.character(levels(adiv[,mvar]));renamed_levels
    oldNames <- unique(adiv[,mvar]);oldNames
    if (length(renamed_levels) == 0) {
      renamed_levels <- oldNames
    }
    for (name in oldNames) {
      total <- length(which(adiv[,mvar] == name));total
      new_n <- paste(name, " (n=", total, ")", sep="");new_n
      levels(adiv[[mvar]])[levels(adiv[[mvar]])== name] <- new_n
      renamed_levels <- replace(renamed_levels, renamed_levels == name, new_n);renamed_levels
    }
    }else{
      adiv <- adiv
    }




    if (mvar == timepoint){
      for (am in alpha_metrics){
        form <-as.formula(sprintf("%s ~ %s + (1|%s)" , am, timepoint, paired))

        clme.mod <- clme(form, data = adiv, constraints = cons, seed = 2, nsim = 200) #
        clme.sum <- summary(clme.mod, seed=2)
        clme.globalp <- function(model) { label <- substitute(
          italic(p) == globalp,
          list(globalp <- model$p.value) )
        as.character(as.expression(format(globalp, nsmall=3)))
        }

        clme.globalp <- paste("CLME P=",clme.globalp(clme.sum))


        # plot design
        if (height*width <= 6){
          dot.size = 0.6
          line.tickness = 0.3
        }else if (height*width > 6 & height*width < 10){
          dot.size = 0.9
          line.tickness = 0.4
        }else{
          dot.size = 1.3
          line.tickness = 0.5
        }

        # plot
        p <- ggplot(adiv, mapping = aes_string(x=timepoint, y=am, color=paired, group=paired)) +
          geom_line(color="grey",size=line.tickness,position = position_dodge(0.3)) +
          geom_point(alpha = 0.8, size = dot.size,position = position_dodge(0.3)) + ylab(sprintf("%s Index\n", am)) +
          ggtitle(sprintf("%s \n (%s) ", mvar, clme.globalp))  +
          theme_bw() + theme(strip.background = element_blank()) +
          theme(title=element_text(size=8), axis.text.x=element_text(angle=xangle,hjust=1,vjust=0.5)) + theme(legend.position= "NONE" )+ theme(aspect.ratio = 1)+
          theme(panel.background = element_rect(fill = "white", colour = "grey50"))


          if(!is.null(mycols)){
           p <- p + scale_color_manual(values = mycols)
           }else{
           p <- p
           }



        if (length(ID) == 1) {
          p= p + geom_text_repel(aes_string(label = ID), size = 2)
        }
        plotlist[[length(plotlist)+1]] <- p
      }
    }else{
      for (des in unique(adiv[,mvar])){
        if(dim(subset(adiv, adiv[,mvar] == des))[1] < 3){
          next
        }
        if(timepoint == mvar){
          next
        }
        print(des)
        for (am in alpha_metrics){
          form <-as.formula(sprintf("%s ~ %s + (1|%s)" , am, timepoint, paired))

          clme.mod <- clme(form, data = adiv[adiv[,mvar] == des,], constraints = cons, seed = 2, nsim = 1000)
          clme.sum <- summary(clme.mod, seed=2)
          clme.globalp <- function(model) { label <- substitute(
            italic(p) == globalp,
            list(globalp <- model$p.value) )
          as.character(as.expression(format(globalp, nsmall=3)))
          }

          clme.globalp <- paste("CLME P=",clme.globalp(clme.sum))


                  # plot design
        if (height*width <= 6){
          dot.size = 0.6
          line.tickness = 0.3
        }else if (height*width > 6 & height*width < 10){
          dot.size = 0.9
          line.tickness = 0.4
        }else{
          dot.size = 1.3
          line.tickness = 0.5
        }

          # plot
          p <- ggplot(adiv[adiv[,mvar]==des,], mapping = aes_string(x=timepoint, y=am, color=paired, group=paired)) +
            geom_line(color="grey",size=line.tickness,position = position_dodge(0.3)) +
            geom_point(alpha = 0.8, size = dot.size,position = position_dodge(0.3)) +
            xlab(timepoint) + ylab(sprintf("%s Index\n", am)) +
            ggtitle(sprintf("%s-%s \n (%s) ", mvar, des, clme.globalp))   + #theme_bw() +
            theme(title=element_text(size=8), axis.text.x=element_text(angle=xangle,hjust=1,vjust=0.5)) + theme(legend.position= "NONE" )


           if(!is.null(mycols)){
           p <- p + scale_color_manual(values = mycols)
           }else{
           p <- p
           }



          # plot size ratio
          if (length(unique(adiv[,mvar])) < 5){
            if(standardsize==TRUE){
              num.subgroup <- length(unique(adiv[,mvar]))*0.2
            }else{
              num.subgroup <- 1
            }
          }else{
            num.subgroup <- length(unique(adiv[,mvar]))*0.2
          }


          p  <- p  + theme(panel.grid = element_blank(),
                           panel.background = element_rect(fill = "white", colour = "Black",size = 0.5, linetype = "solid"),
                           aspect.ratio = 1/num.subgroup)
          plotlist[[length(plotlist)+1]] <- p
        }
      }
    }


  }
  multiplot(plotlist=plotlist, cols=plotCols, rows=plotRows)
  dev.off()
}


