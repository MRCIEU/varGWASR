#' Test for effect of SNP on outcome variance using the LAD-BF model
#'
#' @param data Dataframe of observations
#' @param x Name of SNP dosage
#' @param y Name of outcome
#' @param covar1 Optional vector of covariates to include in the first-stage model
#' @param covar2 Optional vector of covariates to include in the second-stage model
#' @return res Vector of results: var(Y|G==0), var(Y|G==1), var(Y|G==2) and test P-value
#' @export
model <- function(data, x, y, covar1=NULL, covar2=NULL){
    if (any(is.na(data))) stop("Dataframe contains NA values")
    if (!x %in% names(data)) stop(paste0(x, " was not in dataframe"))
    if (!y %in% names(data)) stop(paste0(y, " was not in dataframe"))
    if (!all(covar1 %in% names(data))) stop(paste0("Dataframe is missing some of: ", paste0(covar1, collapse=" ")))
    if (!all(covar2 %in% names(data))) stop(paste0("Dataframe is missing some of: ", paste0(covar2, collapse=" ")))
    if (any(data[[x]] > 2) | any(data[[x]] < 0)) stop("X contains values < 0 | > 2")
    if (!is.numeric(data[[x]])) stop("Genotype must be numeric")

    # prepare first-stage fit matrix
    if (!is.null(covar1)){
        X <- data %>% dplyr::select(!!x, !!covar1) %>% as.matrix
    } else {
        X <- data %>% dplyr::select(!!x) %>% as.matrix
    }
    # first-stage fit
    fit <- cqrReg::qrfit(X=X, y=data[[y]], tau=.5, method="mm")
    b <- rbind(fit$b, fit$beta)
    # predicted
    X <- cbind(rep(1, nrow(X)), X)
    fitted <- X %*% b
    # residual
    d <- data[[y]] - fitted
    # abs residual
    d <- abs(as.vector(d))
    # second-stage model
    data[[x]] <- as.factor(round(data[[x]]))
    if (!is.null(covar2)){
        X <- data %>% dplyr::select(!!x, !!covar2)
        fit2 <- lm(d ~ ., data=X)
        X <- data %>% dplyr::select(!!covar2)
        fit_null <- lm(d ~ ., data=X)
        p <- anova(fit_null, fit2) %>% broom::tidy(.) %>% dplyr::pull(p.value) %>% dplyr::nth(2)
    } else {
        X <- data %>% dplyr::select(!!x)
        fit2 <- lm(d ~ ., data=X)
        fit_null <- lm(d ~ 1, data=data)
        p <- anova(fit_null, fit2) %>% broom::tidy(.) %>% dplyr::pull(p.value) %>% dplyr::nth(2)
    }
    # extract coef
    b0 <- fit2 %>% broom::tidy(.) %>% dplyr::filter(term == "(Intercept)") %>% dplyr::pull("estimate")
    b1 <- fit2 %>% broom::tidy(.) %>% dplyr::filter(term == "x1") %>% dplyr::pull("estimate")
    b2 <- fit2 %>% broom::tidy(.) %>% dplyr::filter(term == "x2") %>% dplyr::pull("estimate")
    res <- c(
        b0^2/(2/pi), # var(Y|G==0)
        b0^2/(2/pi) + (2*b0*b1+b1^2)/(2/pi), # var(Y|G==1)
        b0^2/(2/pi) + (2*b0*b2+b2^2)/(2/pi), # var(Y|G==2)
        p # P-value
    )
    return(res)
}

#' Bootstrap function to obtain SE for LAD-BF effects
#'
#' Usage result <- boot::boot(data=data, statistic=varGWASR::model_bs, R=500, x="SNP", y="outcome", covar1=c("c1", "c2"), covar2=c("c3", "c4")) %>% broom::tidy
#'
#' @param data Dataframe of observations
#' @param indices Vector of nrow indices for bootstrapping SEs
#' @param x Name of SNP dosage
#' @param y Name of outcome
#' @param covar1 Optional vector of covariates to include in the first-stage model
#' @param covar2 Optional vector of covariates to include in the second-stage model
#' @return result Vector of betas
#' @export
model_bs <- function(data, indices, x, y, covar1=NULL, covar2=NULL){
    d <- data[indices,] # allows boot to select sample
    result <- varGWASR::model(d, x, y, covar1=covar1, covar2=covar2)[1:3]
    return(result)
}