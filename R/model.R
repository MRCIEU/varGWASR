#' Test for effect of SNP on outcome variance using the LAD-BF model
#'
#' @param data Dataframe of observations
#' @param x Name of SNP dosage
#' @param y Name of outcome
#' @param covar1 Optional vector of covariates to include in the first-stage model
#' @param covar2 Optional vector of covariates to include in the second-stage model
#' @return Dataframe containing variance effect for SNP=1 (phi_x1) and SNP=2 (phi_x2) with SE and p and F-stat
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
        f <- anova(fit_null, fit2) %>% broom::tidy(.) %>% dplyr::pull(statistic) %>% dplyr::nth(2)
    } else {
        X <- data %>% dplyr::select(!!x)
        fit2 <- lm(d ~ ., data=X)
        fit_null <- lm(d ~ 1, data=data)
        p <- anova(fit_null, fit2) %>% broom::tidy(.) %>% dplyr::pull(p.value) %>% dplyr::nth(2)
        f <- anova(fit_null, fit2) %>% broom::tidy(.) %>% dplyr::pull(statistic) %>% dplyr::nth(2)
    }

    # deltamethod
    v1 <- car::deltaMethod(fit2, "(2*b0*b1+b1^2)/(2/pi)", parameterNames=c("b0", "b1", "b2"))
    v2 <- car::deltaMethod(fit2, "(2*b0*b2+b2^2)/(2/pi)", parameterNames=c("b0", "b1", "b2"))

    res <- data.frame(
        phi_x1=v1$Estimate,
        se_x1=v1$SE,
        phi_x2=v2$Estimate,
        se_x2=v2$SE,
        phi_f=f,
        phi_p=p
    )

    return(res)
}