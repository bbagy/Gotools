
#' Create a Dual Y-Axis Plot
#'
#' @param df Data frame containing the data to be plotted.
#' @param TaxTab Path to the taxonomy table file.
#' @param cate.vars Vector of categorical variables for grouping.
#' @param project Name of the project associated with the data.
#' @param Box Column name in the data frame for boxplot values.
#' @param Line1 Column name in the data frame for the first line plot values.
#' @param Line2 Optional column name for the second line plot values.
#' @param title Optional title for the plot.
#' @param name Optional name for output file.
#' @param mycols Optional color palette for the plot.
#' @param orders Optional order of factor levels for the categorical variables.
#' @param xangle Rotation angle for x-axis labels.
#' @param height Height of the output plot.
#' @param width Width of the output plot.
#'
#' @details
#' This function creates a dual y-axis plot, typically combining a boxplot and one or two line plots.
#' It is useful for visualizing relationships between different types of data (e.g., abundance and environmental variables) across categorical groups.
#'
#' @return
#' A PDF file containing the generated plot.
#'
#' @examples
#' Go_dualYplot(df = data_frame,
#'              TaxTab = "tax_table.csv",
#'              cate.vars = c("Group1", "Group2"),
#'              project = "MyProject",
#'              Box = "Abundance",
#'              Line1 = "Temperature",
#'              Line2 = "pH",
#'              title = "Abundance vs Environmental Factors",
#'              name = "Abundance_EnvFactors",
#'              height = 6,
#'              width = 8)
#'
#' @export

Go_dualYplot <- function(df, TaxTab, cate.vars, project,  Box, Line1, Line2=NULL,
                       title= NULL, 
                       name= NULL, 
                       mycols=NULL, 
                       orders=NULL,
                       xanlgle=90,  height, width){
  
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
  # out file
  pdf(sprintf("%s/dualYplot.%s.%s%s%s.pdf", out_path, 
              project, 
              ifelse(is.null(Line1), "", paste(Line1, ".", sep = "")), 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              format(Sys.Date(), "%y%m%d")), height = height, width = width)


  ## fix factor  and  numeric
  df$etc <- NULL
  df2 <- read.csv(sprintf("%s",TaxTab),header=T,row.names=1,check.names=FALSE);head(df2)
  
  rownames(df2)<- df2$Species
  df2$Species <-NULL
  rownames(df2) <- gsub(" ","_",rownames(df2));rownames(df2)
  
  
  df2 <- as.data.frame(t(df2))

  for (var in cate.vars) {
    print(var)
      df[,var] <- factor(df[,var])
  }
  

  
  # plot
  for (mvar in cate.vars) {
    if (length(unique(df[,mvar])) < 2){
      next
    }
    
    
    # merge merged.df and taxa table
    merged.df <- merge(df, df2, by="row.names");head(merged.df)
    
    # NA remove
    merged.df[,mvar] <- as.character(merged.df[,mvar]);merged.df[,mvar]
    merged.df[,mvar][merged.df[,mvar]==""] <- "NA";merged.df[,mvar]
    merged.df.na <- subset(merged.df, merged.df[,mvar] != "NA");merged.df.na[,mvar]  # subset 를 사용한 NA 삭제
    merged.df.na[,mvar] <- as.factor(merged.df.na[,mvar]);merged.df.na[,mvar]  
    
    # re-order
    if (length(orders) >= 1) {
      merged.df.na[,mvar] <- factor(merged.df.na[,mvar], levels = orders)
    } else {
      merged.df.na[,mvar] <- factor(merged.df.na[,mvar])
    }
    
    # Add number of samples in the group
    renamed_levels <- as.character(levels(merged.df.na[,mvar]));renamed_levels
    oldNames <- unique(merged.df.na[,mvar]);oldNames
    if (length(renamed_levels) == 0) {
      renamed_levels <- oldNames
    }
    for (name in oldNames) {
      total <- length(which(merged.df.na[,mvar] == name));total
      new_n <- paste(name, " (n=", total, ")", sep="");new_n
      levels(merged.df.na[[mvar]])[levels(merged.df.na[[mvar]])== name] <- new_n
      renamed_levels <- replace(renamed_levels, renamed_levels == name, new_n);renamed_levels
    }
    
    
    
    print(sprintf("##-- %s (total without NA: %s/%s) --##", 
                  mvar, dim(merged.df.na)[1], dim(merged.df)[1]))
    
    if (length(unique(merged.df.na[,mvar])) ==1) {
      next
    }
    
    summary.merged.df.na <- summary(merged.df.na[,mvar])
    

    
    #===============================#
    # Visualization for Dual Y axis #
    #===============================#
  
    # for Line1
    mean.line1 <- aggregate(merged.df.na[,Line1], list(merged.df.na[,mvar]), FUN=mean)
    colnames(mean.line1) <- c(mvar, Line1);mean.line1
    mean.line1[,Line1] <- mean.line1[,Line1]*10

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
    
    p <- ggplot() + theme_bw() + theme(strip.background = element_blank()) + #theme_ipsum() +
      geom_boxplot(data=merged.df.na, mapping=aes(x=!!sym(mvar), y=!!sym(Box), colour=!!sym(mvar)), outlier.shape = NA, show.legend = FALSE) +
      theme(text=element_text(size=9), axis.text.x=element_text(angle=xanlgle,hjust=1,vjust=0.5),
            plot.title=element_text(size=9,face="bold"))
      # theme(legend.position="none") +

     if(!is.null(mycols)){
      p <- p + scale_color_manual(values = mycols)
     }else{
       p <- p
     }
    
    
    # count or table for number of variable
    if (max(table(merged.df.na[,mvar])) > 250 & max(table(merged.df.na[,mvar])) < 500){
      dot.size <- dot.size/2
      p = p + geom_jitter(data=merged.df.na, mapping=aes(x=!!sym(mvar), y=!!sym(Box), colour=!!sym(mvar)), 
                          shape=16, alpha = 0.8, size = dot.size, position=position_jitter(0.2), show.legend = FALSE) 
    } else  if (max(table(merged.df.na[,mvar])) < 250 ){
      p = p + geom_jitter(data=merged.df.na, mapping=aes(x=!!sym(mvar), y=!!sym(Box), colour=!!sym(mvar)), 
                          shape=16, alpha = 0.8, size = dot.size, position=position_jitter(0.2), show.legend = FALSE)  
    }else if(max(table(merged.df.na[,mvar])) > 500) {
      dot.size <- dot.size/3
      p = p + geom_jitter(data=merged.df.na, mapping=aes(x=!!sym(mvar), y=!!sym(Box), colour=!!sym(mvar)), 
                          shape=16, alpha = 0.8, size = dot.size, position=position_jitter(0.2), show.legend = FALSE) 
    }
    
    mean.line1.melt <- melt(mean.line1)
    # get label for geom_text_repel
    label.line1 <- subset(mean.line1.melt, variable == Line1);
    n <- dim(label.line1)[1]
    label.line1.sel <- label.line1[n,]
    label.line1.sel$variable <- gsub("_"," ", label.line1.sel$variable)
    
    
    
    p1 <- p + geom_line(data = mean.line1.melt, 
                        mapping = aes(x = !!sym(mvar), y = value, group=variable, color= variable), 
                        inherit.aes = FALSE, size=1)  + guides(color = "none") +
      scale_linetype_manual(values=c("solid", "solid")) + theme(legend.position="top")+
      geom_text_repel(data = label.line1.sel, aes(x = !!sym(mvar), y = value,label = variable),
                      size=3, fontface="italic")
    

   
      
    # for Line2
    if (!is.null(Line2)){
      mean.line1 <- aggregate(merged.df.na[,Line1], list(merged.df.na[,mvar]), FUN=mean)
      colnames(mean.line1) <- c(mvar, Line1);mean.line1
      mean.line1[,Line1] <- mean.line1[,Line1]*10
      
      
      mean.line2 <- aggregate(merged.df.na[,Line2], list(merged.df.na[,mvar]), FUN=mean)
      colnames(mean.line2) <- c(mvar, Line2);mean.line1
      mean.line2[,Line2] <- mean.line2[,Line2]*10
      
      mean.line <- merge(mean.line1, mean.line2, by=mvar);head(mean.line)
      
      mean.line.melt <- melt(mean.line)

      # get label for geom_text_repel
      label.line1 <- subset(mean.line.melt, variable == Line1);
      n <- dim(label.line1)[1]
      label.line1.sel <- label.line1[n,]
      label.line1.sel$variable <- gsub("_"," ", label.line1.sel$variable)
      
      label.line2 <- subset(mean.line.melt, variable == Line2);
      n <- dim(label.line2)[1]
      label.line2.sel <- label.line2[n,]
      label.line2.sel$variable <- gsub("_"," ", label.line2.sel$variable)
      
      
      p1 <- p + geom_line(data = mean.line.melt, 
                          mapping = aes(x = !!sym(mvar), y = value, group=variable, color= variable), 
                          inherit.aes = FALSE, size=1)  + guides(color = "none") +
        scale_linetype_manual(values=c("solid", "solid")) + theme(legend.position="top") +
        geom_text_repel(data = label.line1.sel, aes(x = !!sym(mvar), y = value,label = variable),
                        size=3, fontface="italic")  +
        geom_text_repel(data = label.line2.sel, aes(x = !!sym(mvar), y = value,label = variable),
                        size=3, fontface="italic") 
          
      
      

    }
    

    
    p1 <- p1 + scale_y_continuous(sec.axis = sec_axis(~.*10, name="Relative abundance (%)")) 
    
    if (!is.null(title)) {
      p1 <- p1 + ggtitle(title)
    } else{
      p1 <- p1 + ggtitle(sprintf("%s", mvar))
    }
    print(p1)
  }
  dev.off()
}

