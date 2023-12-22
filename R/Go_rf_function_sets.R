
#' Random Forest Grid Search
#'
#' This function performs a grid search to find the best hyperparameters for a Random Forest model.
#'
#' @param trainSetDf DataFrame containing the training dataset.
#' @param testSetDf DataFrame containing the test dataset.
#' @param outcome Name of the outcome variable in the dataset.
#' @param orders Vector of levels for the outcome variable.
#' @param formula The formula for the Random Forest model.
#'
#' @return A list containing the best mtry, ntree, maxnodes, nodesize, and best_accuracy found during the grid search.
#'
#' @details
#' The function explores combinations of hyperparameters (mtry, ntree, maxnodes, and nodesize) and evaluates the model's performance on the test dataset to find the optimal configuration.
#'
#' @examples
#' # Example usage:
#' Go_rfGridSearch(trainSetDf = trainData, testSetDf = testData, outcome = "Species", orders = c("setosa", "versicolor", "virginica"), formula = Species ~ .)
#'
#' @export

Go_rfGridSearch <- function(trainSetDf,
                            testSetDf,
                            outcome,
                            orders,
                            formula){
  
  
  trainSetDf[,outcome] <- factor(trainSetDf[,outcome], levels = intersect(orders, trainSetDf[,outcome]))
  
  # Define the range of hyperparameters to search over
  mtry_values <- seq(1, floor(sqrt(ncol(trainSetDf) - 1)), 1)
  ntree_values <- c(100, 200, 500, 1000, 1500)
  maxnodes_values <- c(5, 10, 15, 20, 25, 30, 35)
  nodesize_values <- c(1, 3, 5, 7, 10, 15)
  
  # Initialize an array to store the accuracy scores
  accuracy_scores <- array(0, dim = c(length(mtry_values), length(ntree_values), length(maxnodes_values), length(nodesize_values)),
                           dimnames = list(mtry = mtry_values, ntree = ntree_values, maxnodes = maxnodes_values, nodesize = nodesize_values))
  
  # Perform the grid search using loops for mtry, ntree, maxnodes, and nodesize
  for (i in seq_along(mtry_values)) {
    for (j in seq_along(ntree_values)) {
      for (k in seq_along(maxnodes_values)) {
        for (l in seq_along(nodesize_values)) {
          set.seed(123)  # Set a seed for reproducibility
          rf_model <- randomForest(formula, data = trainSetDf, importance = TRUE, mtry = mtry_values[i], ntree = ntree_values[j], maxnodes = maxnodes_values[k], nodesize = nodesize_values[l])
          predicted_classes <- predict(rf_model, newdata = testSetDf, type = "class")
          accuracy_scores[i, j, k, l] <- mean(testSetDf[,outcome] == predicted_classes)
        }
      }
    }
  }
  
  # Find the values of mtry, ntree, maxnodes, and nodesize with the highest accuracy
  best_indices <- which(accuracy_scores == max(accuracy_scores), arr.ind = TRUE)
  best_mtry <- mtry_values[best_indices[1, 1]]
  best_ntree <- ntree_values[best_indices[1, 2]]
  best_maxnodes <- maxnodes_values[best_indices[1, 3]]
  best_nodesize <- nodesize_values[best_indices[1, 4]]
  best_accuracy <- accuracy_scores[best_indices[1, 1], best_indices[1, 2], best_indices[1, 3], best_indices[1, 4]]
  
  # Print the results
  print(paste("Best mtry value:", best_mtry))
  print(paste("Best ntree value:", best_ntree))
  print(paste("Best maxnodes value:", best_maxnodes))
  print(paste("Best nodesize value:", best_nodesize))
  print(paste("Best accuracy:", round(best_accuracy, 3)))
  
  # Return the results
  return(list(best_mtry = best_mtry, best_ntree = best_ntree, best_maxnodes = best_maxnodes, best_nodesize = best_nodesize, best_accuracy = best_accuracy))
}



#' Random Forest Model Training Using Caret Package
#'
#' This function uses the caret package to perform Random Forest model training with hyperparameter tuning.
#'
#' @param trainSetDf DataFrame containing the training dataset.
#' @param testSetDf DataFrame containing the test dataset.
#' @param outcome Name of the outcome variable in the dataset.
#' @param orders Vector of levels for the outcome variable.
#' @param formula The formula for the Random Forest model.
#'
#' @return A dataframe with the best ntree, nodesize, and Accuracy.
#'
#' @details
#' The function applies a cross-validation approach to tune hyperparameters like ntree and nodesize to optimize model performance.
#'
#' @examples
#' # Example usage:
#' Go_rfcaret(trainSetDf = trainData, testSetDf = testData, outcome = "Species", orders = c("setosa", "versicolor", "virginica"), formula = Species ~ .)
#'
#' @export


Go_rfcaret <- function(trainSetDf,
                            testSetDf,
                            outcome,
                            orders,
                            formula){
  
  
  trainSetDf[,outcome] <- factor(trainSetDf[,outcome], levels = intersect(orders, trainSetDf[,outcome]))
  
  # Define the hyperparameter search grid
  control <- trainControl(method = "cv", number = 5, savePredictions = "all", summaryFunction = multiClassSummary)
  hyper_grid <- expand.grid(mtry = seq(1, ncol(train_data), 1)) 
  
  # List of potential values for ntree and nodesize
  ntree_values <- c(100, 200, 300, 400, 500)
  nodesize_values <- c(1, 3, 5, 7, 9)
  
  # Empty data frame to store results
  results <- data.frame()
  
  print(paste("ntree =", ntree_values))
  print(paste("nodesize =", nodesize_values))
  
  # Nested loop to train model with each combination of hyperparameters
  for (ntree in ntree_values) {
    for (nodesize in nodesize_values) {
      set.seed(123)
      
      # Print loop status
      print(paste("Loop status: ntree =", ntree, ", nodesize =", nodesize))
      
      model <-  train(formula, 
                      data = train_data, 
                      method = "rf", 
                      trControl = control,
                      metric = "Accuracy",
                      tuneGrid = hyper_grid,
                      importance = TRUE,
                      proximity = TRUE,
                      ntree = ntree,
                      nodesize = nodesize)
      
      result <- data.frame(ntree = ntree, 
                           nodesize = nodesize, 
                           Accuracy = max(model$results$Accuracy))
      results <- rbind(results, result)
    }
  }
  
  # Print the results
  print(results)
  
  # Identify the row with the maximum accuracy
  best_row <- results[which.max(results$Accuracy), ]
  
  best_row$mtry <- model$bestTune
  #print(best_row)
  return(best_row)
}


#' Random Forest Confusion Matrix Plot
#'
#' This function generates a confusion matrix plot for a Random Forest model.
#'
#' @param rf Random Forest model object.
#' @param project Name of the project or analysis.
#' @param trainSetDf DataFrame containing the training dataset.
#' @param testSetDf DataFrame containing the test dataset.
#' @param outcome Name of the outcome variable.
#' @param orders Vector of levels for the outcome variable.
#' @param name Optional name for the analysis.
#' @param height Height of the plot.
#' @param width Width of the plot.
#'
#' @details
#' The function creates a visual representation of the model's performance on the test dataset using a confusion matrix.
#'
#' @examples
#' # Example usage:
#' Go_rfCMplot(rf = model, project = "MyProject", trainSetDf = trainData, testSetDf = testData, outcome = "Species", orders = c("setosa", "versicolor", "virginica"), name = "Analysis1")
#'
#' @export

Go_rfCMplot <- function(rf, 
                           project, 
                           trainSetDf,
                           testSetDf, 
                           outcome, 
                           orders,
                           name = NULL,
                           height = 2, width=3){
  
  if(!is.null(dev.list())) dev.off()
  
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_rf <- file.path(sprintf("%s_%s/pdf/rf_plot",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_rf)) dir.create(out_rf)
  
  trainSetDf[,outcome] <- factor(trainSetDf[,outcome], levels = intersect(orders, trainSetDf[,outcome]))
  testSetDf[,outcome] <- factor(testSetDf[,outcome], levels = intersect(orders, testSetDf[,outcome]))
  
  # confusion matrix plot for test set
  
  if(class(rf)=="list"){
    print("Ensemble RF model")
    # Create a matrix to store individual predictions for the test data
    predictions_label_test <- matrix(NA, nrow = nrow(testSetDf), ncol = ntree)
    
    for (i in 1:ntree) {
      # Obtain predictions from each individual model for the test data
      predictions_label_test[, i] <- predict(ensemble[[i]], newdata = testSetDf, type = "prob")[, 2]
    }
    
    # Obtain the predicted classes from the ensemble model for the test data
    predicted_classes <- ifelse(rowMeans(predictions_label_test) >= 0.5,orders[1], orders[2])
    

    
    
    # Create the confusion matrix
    confusionTab <- data.frame(table(testSetDf[,outcome], predicted_classes, dnn = c("Actual", "Predicted")))
    
    
    
    
    
    Accuracy <- round(mean(testSetDf[,outcome] == predicted_classes, na.rm=TRUE) * 100, 3);Accuracy
    print(paste("Prediction Accuracy:", Accuracy))
    

    
    
    p.confusion <- ggplot(data =  confusionTab, mapping = aes(x = Actual, y = Predicted)) +
      geom_tile(aes(fill = Freq), colour = "white") +geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1) +
      scale_fill_gradient(low = "#2F86FFFF" , high = "#F14809FF") +
      theme(plot.title = element_text(size = 10), legend.position = "none")+ ggtitle(paste("Ensemble RF Prediction Accuracy (test set) \n", Accuracy,"%",sep=""))+
      theme(panel.grid = element_blank(),  panel.background = element_rect(fill = "white", colour = "Black",size = 0.7, linetype = "solid"), 
            aspect.ratio = 0.5/1) 
    
    
    
  }else{
    
    tt <- try(predicted_classes <- predict(rf, newdata = testSetDf, type = "class"),T)
    
    if(class(tt)== "try-error"){
      cm <- confusionMatrix(rf) 
      
      confusionTab <- as.data.frame(cm$table)
      colnames(confusionTab) <- c("Actual", "Predicted","Freq")
      # Manually calculate the accuracy
      Accuracy <- round(sum(diag(cm$table)) / sum(cm$table)*100,3)
      oob_error <- round(rf$finalModel$err.rate[nrow(rf$finalModel$err.rate), "OOB"]*100,3)
      
      cat(sprintf("Cross-Validated (5 fold) Confusion Matrix\n Prediction Accuracy (test set): %s%% \n Out-of-Bag (OOB) error rate: %s%% \n\n", Accuracy, oob_error))
      
      
      p.confusion <- ggplot(data =  confusionTab, mapping = aes(x = Actual, y = Predicted)) +
        geom_tile(aes(fill = Freq), colour = "white") +geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1) +
        scale_fill_gradient(low = "#2F86FFFF" , high = "#F14809FF") +
        theme(plot.title = element_text(size = 10), legend.position = "none")+ ggtitle(paste("Ensemble RF Prediction Accuracy (test set) \n", Accuracy,"%",sep=""))+
        theme(panel.grid = element_blank(),  panel.background = element_rect(fill = "white", colour = "Black",size = 0.7, linetype = "solid"), 
              aspect.ratio = 0.5/1) 
      
    }else{
      predicted_classes <- predict(rf, newdata = testSetDf, type = "class")
      confusionTab <- data.frame(table(testSetDf[,outcome], predicted_classes, dnn = c("Actual", "Predicted"))) # 대각선은 예측에 성공/ 비대각선은 실폐
      Accuracy <- round(mean(testSetDf[,outcome] == predicted_classes, na.rm=TRUE) * 100, 3);Accuracy
      print(paste("Prediction Accuracy:", Accuracy))
      
      p.confusion <- ggplot(data =  confusionTab, mapping = aes(x = Actual, y = Predicted)) +
        geom_tile(aes(fill = Freq), colour = "white") +geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1) +
        scale_fill_gradient(low = "#2F86FFFF" , high = "#F14809FF") +
        theme(plot.title = element_text(size = 10), legend.position = "none")+ ggtitle(paste("Prediction Accuracy (test set) \n", Accuracy,"%",sep=""))+
        theme(panel.grid = element_blank(),  panel.background = element_rect(fill = "white", colour = "Black",size = 0.7, linetype = "solid"), 
              aspect.ratio = 0.5/1) 
      
      # Compute the OOB confusion matrix
      oob_error_rate <- round(rf$err.rate[nrow(rf$err.rate), "OOB"] *100,3)
      oob_confusion_matrix <- data.frame(confusionMatrix(rf$predicted, trainSetDf[,outcome])$table)
      
      colnames(oob_confusion_matrix) <- c("Predicted", "Actual", "Freq")
      
      print(paste("OOB (out of bag) estimate of error rate:", oob_error_rate))
      p.odd.cm <- ggplot(data =  oob_confusion_matrix, mapping = aes(x = Actual, y = Predicted)) +theme_bw() + 
        geom_tile(aes(fill = Freq), colour = "white") +geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1) +
        scale_fill_gradient(low = "#2F86FFFF" , high = "#F14809FF") +
        theme(plot.title = element_text(size = 10), legend.position = "none")+ ggtitle(paste("OOB estimate of error rate \n", oob_error_rate,"%",sep=""))
      
    }
  }
  
  
  

  
  
  pdf(sprintf("%s/rf.cm.plot.%s.%s.%s%s.pdf", out_rf, 
              outcome,
              project, 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              format(Sys.Date(), "%y%m%d")), height = height, width = width)

  if(class(rf)=="list"){
    print(p.confusion)
  }else{
    if(class(tt)== "try-error"){
      print(p.confusion)
    }else{
      print(p.confusion)
      print(p.odd.cm)
    }
  }
 
  
  
  dev.off()
}


#' Out-of-Bag Error Rate Plot for Random Forest
#'
#' This function generates a plot of the Out-of-Bag (OOB) error rate across different numbers of trees in a Random Forest model.
#'
#' @param rf Random Forest model object.
#' @param project Name of the project or analysis.
#' @param trainSetDf DataFrame containing the training dataset.
#' @param orders Vector of levels for the outcome variable.
#' @param outcome Name of the outcome variable.
#' @param mycols Colors to use in the plot.
#' @param name Optional name for the analysis.
#' @param height Height of the plot.
#' @param width Width of the plot.
#'
#' @details
#' The function visualizes how the OOB error rate changes as more trees are added to the model, providing insight into model convergence.
#'
#' @examples
#' # Example usage:
#' Go_rfoobError(rf = model, project = "MyProject", trainSetDf = trainData, orders = c("setosa", "versicolor", "virginica"), outcome = "Species", name = "Analysis1")
#'
#' @export


Go_rfoobError <- function(rf, 
                        project, 
                        trainSetDf,
                        orders, 
                        outcome, 
                        mycols=NULL,
                        name= NULL,
                        height = 2, width=3){
  if(!is.null(dev.list())) dev.off()
  
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_rf <- file.path(sprintf("%s_%s/pdf/rf_plot",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_rf)) dir.create(out_rf)
  
  tt <- try(mycols,T)
  if(class(tt) == "try-error"){
    print("mycols is not defined.")
    mycols <- NULL
  }
  
  set.seed(123) 

  trainSetDf[,outcome] <- factor(trainSetDf[,outcome], levels = intersect(orders, trainSetDf[,outcome]))
  label1 <- as.vector(sort(unique(trainSetDf[,outcome]), levels = intersect(orders, trainSetDf[,outcome])))[1]
  label2 <- as.vector(sort(unique(trainSetDf[,outcome]), levels = intersect(orders, trainSetDf[,outcome])))[2]
  
  


  oob_error_rates <- rf$err.rate
  
  # Create a data frame
  oob_data <- data.frame(Trees = 1:nrow(oob_error_rates),
                         OOB_Error_Rate_label1 = oob_error_rates[, label1],
                         OOB_Error_Rate_label2 = oob_error_rates[, label2],
                         OOB_Error_Rate_Total = oob_error_rates[, "OOB"])
  
  
  
  # Calculate the OOB error rates at the end of the trees
  oob_error_rate_label1 <- round(oob_data$OOB_Error_Rate_label1[nrow(oob_data)] * 100, 1)
  oob_error_rate_label2 <- round(oob_data$OOB_Error_Rate_label2[nrow(oob_data)] * 100, 1)
  oob_error_rate_total <- round(oob_data$OOB_Error_Rate_Total[nrow(oob_data)] * 100, 1)
  
  # Update the labels for the legend
  group1_label <- paste0(label1, "(", oob_error_rate_label1, "%)")
  group2_label <- paste0(label2, "(", oob_error_rate_label2, "%)")
  total_label <- paste0("Total OOB (", oob_error_rate_total, "%)")
  
  # Convert to long format
  oob_data_long <- oob_data %>%
    pivot_longer(cols = -Trees,
                 names_to = "Class",
                 values_to = "Error_Rate")
  
  
  # Plot with ggplot2
  p.oob<- ggplot(data = oob_data_long, aes(x = Trees, y = Error_Rate, color = Class)) + geom_line() +
    labs(title = "OOB Error Rate", x = "Number of Trees",y = "OOB Error Rate") +
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", colour = "Black", size = 0.7, linetype = "solid"),
          aspect.ratio = 1,
          legend.position = c(0.8, 0.2))
  
  if(!is.null(mycols)){
    p.oob <- p.oob +  scale_color_manual(values = mycols, name = "",
                                    labels = c(group1_label, group2_label, total_label))
  }else{
    p.oob <- p.oob +scale_color_discrete(name = "", labels = c(group1_label, group2_label, total_label))
  }
  
  pdf(sprintf("%s/rf.ooberrorrate.plot.%s.%s.%s%s.pdf", out_rf, 
              outcome,
              project, 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              format(Sys.Date(), "%y%m%d")), height = height, width = width)
  print(p.oob)
  dev.off()
}



#===== Go_rfrocplot

Go_rfrocplot <- function(rf, 
                          project, 
                          trainSetDf,
                          testSetDf,
                          orders, 
                          outcome, 
                          mycols=NULL,
                          name= NULL,
                          height = 2, width=3){
  if(!is.null(dev.list())) dev.off()
  set.seed(123) 
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_rf <- file.path(sprintf("%s_%s/pdf/rf_plot",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_rf)) dir.create(out_rf)
  
  tt <- try(mycols,T)
  if(class(tt) == "try-error"){
    print("mycols is not defined.")
    mycols <- NULL
  }
  
  #===== ROC plot ========#
  df <- rbind(data.frame(predictor = predict(rf, newdata = trainSetDf, type = "prob")[,2],
                         known.truth = trainSetDf[,outcome],
                         model = "Train_set"),
              data.frame(predictor = predict(rf, newdata = testSetDf, type = "prob")[,2],
                         known.truth = testSetDf[,outcome],
                         model = "Test_set"))
  
  library(plotROC)
  
  
  ggroc <- ggplot(df, aes(d = known.truth, m = predictor, color = model, group = model)) + 
    geom_roc(n.cuts = 0,size=0.7) +
    theme(panel.grid = element_blank(),  panel.background = element_rect(fill = "white", colour = "Black",size = 0.7, linetype = "solid"), 
          aspect.ratio = 1)
  
  
  if(!is.null(mycols)){
    ggroc <- ggroc +  scale_color_manual(values = mycols, name="",
                                         limits = c("Train_set", "Test_set"))
  }else{
    ggroc <- ggroc +scale_color_discrete(name="",
                                         limits = c("Train_set", "Test_set"))
  }
  
  
  
  
  
  # AUC (area under ROC )
  aucTab <- calc_auc(ggroc)
  
  # Get the AUC values
  train_auc <- aucTab$AUC[2]
  test_auc <- aucTab$AUC[1]
  
  # Print the AUC values
  cat("The AUC value for the training set is: ", train_auc, "\n")
  cat("The AUC value for the test set is: ", test_auc, "\n")
  
  
  
  annotate <- ggroc + theme(plot.title = element_text(hjust = 0,size = 10), legend.position = c(0.75, 0.17)) +
    ggtitle("ROC plots") +xlab("False positive") + ylab("True positive") +
    annotate("text", x = 0.75, y = 0.4, label = paste("Train_set AUC =", round(aucTab$AUC[2], 3), "\n","Test_set AUC =", round(aucTab$AUC[1], 3) ),size = 2.5) +geom_abline(intercept = 0, slope = 1, color = "darkgrey", linetype = "dashed")
  
  pdf(sprintf("%s/rf.roc.plot.%s.%s.%s%s.pdf", out_rf, 
              outcome,
              project, 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              format(Sys.Date(), "%y%m%d")), height = height, width = width)
  print(annotate)
  dev.off()
}




#===== Go_rfproximitiesplot
Go_rfproximityplot <- function(rf, 
                         project, 
                         trainSetDf,
                         testSetDf,
                         orders, 
                         outcome, 
                         mycols=NULL,
                         name= NULL,
                         height = 2, width=3){
  if(!is.null(dev.list())) dev.off()
  
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_rf <- file.path(sprintf("%s_%s/pdf/rf_plot",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_rf)) dir.create(out_rf)
  
  tt <- try(mycols,T)
  if(class(tt) == "try-error"){
    print("mycols is not defined.")
    mycols <- NULL
  }
  
  set.seed(123) 
  trainSetDf[,outcome] <- factor(trainSetDf[,outcome], levels = intersect(orders, trainSetDf[,outcome]))
  label1 <- unique(trainSetDf[,outcome])[1]
  label2 <- unique(trainSetDf[,outcome])[2]
  
  
  # Calculate the proximity matrix
  
  proximity_matrix <- rf$proximity
  
  
  if(class(rf) == "list"){
    proximity_matrices <- lapply(rf, function(model) {
      return(model$finalModel$proximity)
    })
    # Check that all matrices are the same size
    sizes <- sapply(proximity_matrices, dim)
    if (!all(sizes == sizes[1,])) {
      stop("All proximity matrices must be the same size")
    }
    
    # Assuming proximity_matrices is a list of matrices
    proximity_matrix <- Reduce("+", proximity_matrices) / length(proximity_matrices)
    
    
    
  }else if (is.null(proximity_matrix)){
    proximity_matrix <- rf$finalModel$proximity
  }
  
  
  if(class(rf) == "list"){
    title <- "Average Proximity Plot (t-SNE)"
  }else{
    title <- "Proximity Plot (t-SNE)"
  }
  
  # Perform multidimensional scaling (MDS) to reduce the dimensionality of the proximity matrix
  # mds_coords <- cmdscale(prox_matrix, k = 2)
  
  # Create a color vector based on the classes in your dataset (assuming you have a binary outcome)
  colors <- ifelse(trainSetDf[,outcome] == label1, "blue", "red")
  
  
  # Obtain the proximity matrix for the test set


  
  # Load the required package
  library(Rtsne)
  
  # Perform t-SNE
  tsne_result <- Rtsne(proximity_matrix, is_distance = TRUE, perplexity = 30, check_duplicates = FALSE)
  
  
  # Create a data frame with the t-SNE results and labels
  tsne_data <- data.frame(tsne_x = tsne_result$Y[, 1], tsne_y = tsne_result$Y[, 2], label = trainSetDf[,outcome])
  
  # Plot the t-SNE results using ggplot2
  proximity <- ggplot(tsne_data, aes(x = tsne_x, y = tsne_y, color = label)) + 
    geom_point(size = 1, alpha = 0.8) + ggtitle(title) + labs(x = "t-SNE 1", y = "t-SNE 2")+ 
    theme(panel.grid = element_blank(),  panel.background = element_rect(fill = "white", colour = "Black",size = 0.7, linetype = "solid"), 
          aspect.ratio = 1) +
    geom_vline(xintercept = 0, size = 0.1) + geom_hline(yintercept = 0, size = 0.1) 
  
  
  proximity <-ggplot(tsne_data, aes(x = tsne_x, y = tsne_y, color = label)) + 
    geom_point(size = 1, alpha = 0.8) + ggtitle("Proximity Plot (t-SNE)") + 
    labs(x = "t-SNE 1", y = "t-SNE 2") + 
    theme(panel.grid = element_blank(),  
          panel.background = element_rect(fill = "white", colour = "Black", size = 0.7, linetype = "solid"), 
          aspect.ratio = 1, 
          plot.title = element_text(hjust = 0,size = 10), legend.position = c(0.75, 0.17)) +
    geom_vline(xintercept = 0, size = 0.1) + geom_hline(yintercept = 0, size = 0.1)
  
  
  
  
  if(!is.null(mycols)){
    proximity <- proximity +  scale_color_manual(values = mycols, name="")
  }else{
    proximity <- proximity +scale_color_discrete(name="")
  }
  
  
  
  pdf(sprintf("%s/rf.proximityPlot.%s.%s.%s%s.pdf", out_rf, 
              outcome,
              project, 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              format(Sys.Date(), "%y%m%d")), height = height, width = width)
  print(proximity)
  dev.off()
}



Go_rfimportanceplot <- function(rf, 
                                 project, 
                                 name= NULL,
                                 height = 2, width=3){
  if(!is.null(dev.list())) dev.off()
  
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_rf <- file.path(sprintf("%s_%s/pdf/rf_plot",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_rf)) dir.create(out_rf)
  

  set.seed(123) 
  #===== importance plot ========#
  #varImpPlot(rf, main = "Variable Importance")
  
  
  
  if (class(rf) == "list"){
    compute_importance <- function(model) {
      return(varImp(model, scale = FALSE)$importance)
    }
    
    # Compute variable importance for each model in the ensemble
    importance_list <- lapply(rf, compute_importance)
    # Calculate the average feature importance across all models in the ensemble
    average_importance <- Reduce('+', importance_list) / length(rf)
    
    # Create a dataframe for easier manipulation
    average_importance_df <- data.frame(Features = rownames(average_importance), Importance = rowSums(average_importance))
    
    # Order the dataframe based on importance
    importance_data <- average_importance_df[order(-average_importance_df$Importance),]
    
    ytitle <- "Average Importance (Gini)"

    importance <- ggplot(importance_data, aes(x = reorder(Features, Importance), y = Importance)) 
  }else{
    # Calculate variable importances
    importance_data <- data.frame(Feature = rownames(rf$importance),
                                  Importance = rf$importance[, "MeanDecreaseGini"])
    if (dim(importance_data)[1] == 0){
      importance_data <- data.frame(Feature = rownames(rf$finalModel$importance),
                                    Importance = rf$finalModel$importance[, "MeanDecreaseGini"])
    }
    
    # Sort the data by importance
    importance_data <- importance_data %>% arrange(desc(Importance))
    
    ytitle <- "Mean Decrease in Gini"
    
    importance <- ggplot(importance_data, aes(x = reorder(Features, Importance), y = Importance)) 
  }
  
  

  
  
  dim(importance_data)[1]
  
  # plot size ratio
  if (dim(importance_data)[1] < 6){
    num.subgroup <- 1
  }else{
    num.subgroup <- dim(importance_data)[1]*0.15
  }
  
  # Create the importance plot using ggplot2
  importance <- importance +
    geom_bar(stat = "identity", fill = "steelblue") + coord_flip() + 
    xlab("Features") + ylab(ytitle) + ggtitle("Variable Importance") + 
    theme(panel.grid = element_blank(),  panel.background = element_rect(fill = "white", colour = "Black",size = 0.7, linetype = "solid"), 
          aspect.ratio = num.subgroup/1) 
  
  
  pdf(sprintf("%s/rf.importancePlot.%s.%s%s.pdf", out_rf, 
              project, 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              format(Sys.Date(), "%y%m%d")), height = height, width = width)
  print(importance)
  dev.off()
}





