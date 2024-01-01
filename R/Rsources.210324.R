
#' Remove specified rows from a data frame or matrix
#'
#' @param rowNum A vector of row numbers to be removed.
#' @param data The data frame or matrix from which rows will be removed.
#'
#' @return A data frame or matrix with specified rows removed.
#'
#' @examples
#' data <- data.frame(x = 1:5, y = 6:10)
#' removeRows(c(2, 4), data)
#'
#' @export
removeRows <- function(rowNum, data) {
    newData <- data[-rowNum, , drop = FALSE]
    rownames(newData) <- NULL
    newData
}


#' Retrieve Taxonomy Information for Given OTUs
#'
#' @param otus Vector of OTUs for which taxonomy information is required.
#' @param tax_tab Taxonomy table containing taxonomic information for each OTU.
#' @param level Taxonomic level at which information is sought.
#' @param taxRanks Vector of all possible taxonomic ranks.
#' @param na_str Character string to represent missing data in taxonomy.
#'
#' @return A data frame with taxonomy information for the specified OTUs.
#'
#' @examples
#' # Example usage pending specific data structure
#'
#' @export

getTaxonomy <- function(otus, tax_tab, level, taxRanks,na_str = c( "NA")) {
    ranks <- taxRanks
    sel <- ranks[1:match(level, ranks)]
    inds <- apply(tax_tab[otus,sel], 1, function(x) max(which(!(x %in% na_str))))
    retval <- as.data.frame(tax_tab)[cbind(otus, ranks[inds])]
    retval[inds!=match(level, ranks)] <- paste(na_str[1], retval[inds!=match(level, ranks)], sep="_")
    return(retval)
}

#' Retrieve Comprehensive Taxonomy Information
#'
#' @param otus Vector of OTUs for taxonomy information.
#' @param tax_tab Taxonomy table with taxonomic details.
#' @param level Desired taxonomic level.
#' @param na_str String to represent missing data.
#'
#' @return A data frame with comprehensive taxonomy information.
#'
#' @examples
#' # Example usage pending specific data structure
#'
#' @export

getTaxonomyAll <- function(otus, tax_tab, level, na_str = c( "NA")) {
    rank <- rank
    sel <- rank[1:match(level, rank)]
    inds <- apply(tax_tab[otus,sel], 1, function(x) max(which(!(x %in% na_str))))
    retval <- as.data.frame(tax_tab)[cbind(otus, rank[inds])]
    retval[inds!=match(level, rank)] <- paste(na_str[1], retval[inds!=match(level, rank)], sep="_")
    return(retval)
}



#' Calculate Matthews Correlation Coefficient
#'
#' @param preds Predicted values (optional if `x` is provided).
#' @param actuals Actual values (optional if `y` is provided).
#' @param x Indicator matrix for predictions (if `preds` is not provided).
#' @param y Indicator matrix for truth (if `actuals` is not provided).
#'
#' @return Matthews Correlation Coefficient value.
#'
#' @examples
#' # Example usage pending specific data structure
#'
#' @export
mcc <- function(preds=NULL, actuals=NULL, x=NULL, y=NULL) {
    # if preds and actuals are provided, x and y will be ignored
    if (!is.null(preds)) {
        nclasses <- length(union(preds, actuals))
        x <- matrix(0, nrow=length(preds), ncol=nclasses)
        y <- matrix(0, nrow=length(actuals), ncol=nclasses)
        x[cbind(1:nrow(x), preds+1)] <- 1
        y[cbind(1:nrow(y), actuals+1)] <- 1
    }
    if (!all(dim(x) == dim(y))) {
        stop("X and Y must have the same dimensions")
    }
    
    cov_biased <- function(x, y) {
        sum(sapply(1:ncol(x), function(k) {
            cov(x[,k], y[,k]) # unbiased estimate with (n-1) denominator as opposed to (n), but cancels out anyways so identical result
        }))
    }
    numerator <- cov_biased(x,y)
    denominator <- sqrt(cov_biased(x,x) * cov_biased(y,y))
    numerator / denominator
}


#' Plot Multiple ggplot Objects in a Single Layout
#'
#' @param ... ggplot objects to be plotted.
#' @param plotlist List of ggplot objects (alternative to ...).
#' @param file Filename if the plot is to be saved.
#' @param cols Number of columns in the plot layout.
#' @param rows Number of rows in the plot layout.
#'
#' @return A grid layout combining multiple ggplot objects.
#'
#' @examples
#' # Example usage with ggplot2 plots
#'
#' @export
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

#' Normalize Data by Rows or Columns
#'
#' @description
#' `normalizeByRows` normalizes the data by rows.
#' `normalizeByCols` normalizes the data by columns.
#'
#' @param df Data frame or matrix to be normalized.
#' @param rsum Target sum for row normalization.
#' @param csum Target sum for column normalization.
#' @param level Taxonomic level for normalization (for `normalizeByCols`).
#' @param delim Delimiter for splitting taxonomy strings (for `normalizeByCols`).
#'
#' @return Normalized data frame or matrix.
#'
#' @examples
#' data <- matrix(rnorm(20), nrow = 4)
#' normalizeByRows(data)
#' normalizeByCols(data)
#'
#' @export

normalizeByRows <- function (df, rsum=1)
{
    while (any(abs((rowSums(df)-rsum))>1e-13)) {
        df <- rsum*(df / rowSums(df))
    }
    return(df)
}



normalizeByCols <- function (df, csum=1, level=NULL, delim="\\|")
{
    if (is.null(level)) {
        while (any(abs((colSums(df)-csum))>1e-13 & colSums(df)!=0, na.rm=T)) {
            missing <- which(colSums(df)==0)
            df <- sweep(df, 2, colSums(df)/csum, "/")
            df[,missing] <- 0
        }
    } else {
        tmp <- df
        tmp$taxa <- rownames(tmp)
        tmp$splitter <- factor(unlist(lapply(rownames(tmp), function(x) unlist(strsplit(x, delim))[level])))
        names <- rownames(tmp)[order(tmp$splitter)]
        tmp <- ddply(tmp, .(splitter), function(x) {
            x <- x[, setdiff(colnames(x), c("taxa", "splitter"))]
            while (any(abs((colSums(x)-csum))>1e-13 & colSums(df)!=0, na.rm=T)) {
                x <- sweep(x, 2, colSums(x)/csum, "/")
            }
            x
        })
        rownames(tmp) <- names
        df <- tmp[, setdiff(colnames(tmp), "splitter")]
    }
    return(df)
}


#' Rename Levels with Counts
#'
#' @param fvec Factor vector with levels to be renamed.
#' @param originalLevelsAsNames Logical; if TRUE, original levels are used as names.
#'
#' @return A factor vector with levels renamed to include counts.
#'
#' @examples
#' vec <- factor(c("A", "B", "A", "C"))
#' renameLevelsWithCounts(vec)
#'
#' @export

renameLevelsWithCounts <- function(fvec, originalLevelsAsNames=FALSE) {
    tab <- table(fvec)
    retval <- sprintf("%s (n=%d)", fvec, tab[unlist(lapply(fvec, function(x) match(x, names(tab))))])
    #    newlevels <- sprintf("%s (n=%d)", levels(fvec), tab[levels(fvec)])
    newlevels <- sprintf("%s (n=%d)", levels(fvec), tab[unlist(lapply(names(tab), function(x) which(levels(fvec)==x)))])
    retval <- factor(retval, levels=newlevels)
    if (originalLevelsAsNames) {
        names(retval) <- fvec
    }
    return(retval)
}

#' Summarize Data by Group
#'
#' This function provides a summary (mean and standard deviation) of a specified variable within a data frame, grouped by one or more factors.
#'
#' @param data A data frame containing the data to be summarized.
#' @param varname The name of the variable to summarize.
#' @param groupnames A vector of names of grouping variables.
#'
#' @return A data frame with the mean and standard deviation of the specified variable for each group.
#'
#' @examples
#' data <- data.frame(group = rep(c("A", "B"), each = 5), value = rnorm(10))
#' data_summary(data, varname = "value", groupnames = "group")
#'
#' @export
#' @import plyr
data_summary <- function(data, varname, groupnames){
    require(plyr)
    summary_func <- function(x, col){
        c(mean = mean(x[[col]], na.rm=TRUE),
        sd = sd(x[[col]], na.rm=TRUE))
    }
    data_sum<-ddply(data, groupnames, .fun=summary_func,
    varname)
    data_sum <- rename(data_sum, c("mean" = varname))
    return(data_sum)
}


#' Calculate Geometric Mean of Pairwise Ratios (GMPR) Size Factors
#'
#' This function computes the GMPR size factors for microbiome sequencing data or zero-inflated sequencing data. The size factors can be used as offsets in count-based regression models or as divisors for normalized data.
#'
#' @param comm A matrix of counts (rows - features like OTUs or genes, columns - samples).
#' @param intersect.no Minimum number of shared features between sample pairs for ratio calculation.
#' @param ct.min Minimum count threshold for inclusion in ratio calculation.
#' @param trace Logical; if TRUE, progress messages will be printed.
#'
#' @return A vector of GMPR size factors. Samples with distinct sets of features will be output as NA.
#'         Includes an attribute 'NSS' indicating the number of samples with significant sharing.
#'
#' @examples
#' # Example usage pending specific data structure
#'
#' @export

GMPR <- function (comm, intersect.no = 10, ct.min = 1, trace = TRUE) {
  comm[comm < ct.min] <- 0
  
  if (is.null(colnames(comm))) {
    colnames(comm) <- paste0('S', 1:ncol(comm))
  }
  
  if (trace) cat('Begin GMPR size factor calculation ...\n')
  
  comm.no <- numeric(ncol(comm))
  gmpr <- sapply(1:ncol(comm),  function(i) {
    if (i %% 50 == 0) {
      cat(i, '\n')
    }
    x <- comm[, i]
    # Compute the pairwise ratio
    pr <- x / comm
    # Handling of the NA, NaN, Inf
    pr[is.nan(pr) | !is.finite(pr) | pr == 0] <- NA
    # Counting the number of non-NA, NaN, Inf
    incl.no <- colSums(!is.na(pr))
    # Calculate the median of PR
    pr.median <- colMedians(pr, na.rm=TRUE)
    # Record the number of samples used for calculating the GMPR
    comm.no[i] <<- sum(incl.no >= intersect.no)
    # Geometric mean of PR median
    if (comm.no[i] > 1) {
      return(exp(mean(log(pr.median[incl.no >= intersect.no]))))
    } else {
      return(NA)
    }
  }
  )
  
  if (sum(is.na(gmpr))) {
    warning(paste0('The following samples\n ', paste(colnames(comm)[is.na(gmpr)], collapse='\n'),
                   '\ndo not share at least ', intersect.no, ' common taxa with the rest samples! ',
                   'For these samples, their size factors are set to be NA! \n',
                   'You may consider removing these samples since they are potentially outliers or negative controls!\n',
                   'You may also consider decreasing the minimum number of intersecting taxa and rerun the procedure!\n'))
  }
  
  if (trace) cat('Completed!\n')
  if (trace) cat('Please watch for the samples with limited sharing with other samples based on NSS! They may be outliers! \n')
  names(gmpr) <- names(comm.no) <- colnames(comm)
  
  attr(gmpr, 'NSS') <- comm.no
  
  return(gmpr)
}
