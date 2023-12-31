% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Go_rf_function_sets.R
\name{Go_rfcaret}
\alias{Go_rfcaret}
\title{Random Forest Model Training Using Caret Package}
\usage{
Go_rfcaret(trainSetDf, testSetDf, outcome, orders, formula)
}
\arguments{
\item{trainSetDf}{DataFrame containing the training dataset.}

\item{testSetDf}{DataFrame containing the test dataset.}

\item{outcome}{Name of the outcome variable in the dataset.}

\item{orders}{Vector of levels for the outcome variable.}

\item{formula}{The formula for the Random Forest model.}
}
\value{
A dataframe with the best ntree, nodesize, and Accuracy.
}
\description{
This function uses the caret package to perform Random Forest model training with hyperparameter tuning.
}
\details{
The function applies a cross-validation approach to tune hyperparameters like ntree and nodesize to optimize model performance.
}
\examples{
# Example usage:
Go_rfcaret(trainSetDf = trainData, testSetDf = testData, outcome = "Species", orders = c("setosa", "versicolor", "virginica"), formula = Species ~ .)

}
