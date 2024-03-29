% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.R
\name{model}
\alias{model}
\title{Test for effect of SNP on outcome variance using the LAD-BF model}
\usage{
model(data, x, y, covar1 = NULL, covar2 = NULL)
}
\arguments{
\item{data}{Dataframe of observations}

\item{x}{Name of SNP dosage}

\item{y}{Name of outcome}

\item{covar1}{Optional vector of covariates to include in the first-stage model}

\item{covar2}{Optional vector of covariates to include in the second-stage model}
}
\value{
Dataframe containing variance effect for SNP=1 (phi_x1) and SNP=2 (phi_x2) with SE and p and F-stat
}
\description{
Test for effect of SNP on outcome variance using the LAD-BF model
}
