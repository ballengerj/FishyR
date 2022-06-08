#-------------------------------------------------------------------------------
#' Perform Model Selection Using Generalized Linear Models (GLMs)
#'
#' \code{GLM.Select} generates a vector that contains an \code{\link{AIC}}, an
#'   \code{\link[MuMIn]{AICc}}, a \code{\link{BIC}}, and a dispersion parameter
#'   estimate.  Intent is to use this function in a loop to extract model
#'   information criterion for a range of model formulations, facilitating the
#'   automated selection of *best* model.
#' @param formula character string representing a potential model formula.
#' @param data Data frame containing columns with names corresponding to the
#'   values in the formula.
#' @param family a description of the error distribution and link function to
#'   be used in the model (see \code{\link{glm}} and \code{\link{family}}.
#'   Outside of those values that can be used in \code{\link{glm}}, you can
#'   also specify it as *nb*.  This tells the model to use the
#'   \code{\link[MASS]{glm.nb}} function available in the package **MASS**.
#' @return vector that contains the \code{\link{AIC}},
#'   \code{\link[MuMIn]{AICc}}, \code{\link{BIC}} and dispersion parameter
#'   estimate for a given model structure.
#' @family Model Selection
#' @export
GLM.Select <- function(formula, data, family) {
    if(family == "nb") {
        tmp <- try(MASS::glm.nb(as.formula(formula), data = data))
        w <- ifelse(class(try(AIC(tmp))) == "numeric", AIC(tmp), NA)
        y <- ifelse(class(try(MuMIn::AICc(tmp))) == "numeric",
                    MuMIn::AICc(tmp),
                    NA)
        z <- ifelse(class(try(BIC(tmp))) == "numeric", BIC(tmp), NA)
        Disp <- ifelse(class(try(AIC(tmp))) == "numeric",
                       sum(resid(tmp, type="pearson") ^ 2)/ tmp$df.resid,
                       NA)
        c(w, y, z, Disp)
    } else {
        tmp <- try(glm(as.formula(formula), data = data,
                       family = family))
        w <- ifelse(class(try(AIC(tmp))) == "numeric", AIC(tmp), NA)
        y <- ifelse(class(try(MuMIn::AICc(tmp))) == "numeric",
                    MuMIn::AICc(tmp),
                    NA)
        z <- ifelse(class(try(BIC(tmp))) == "numeric", BIC(tmp), NA)
        Disp <- ifelse(class(try(AIC(tmp))) == "numeric",
                       sum(resid(tmp, type = "pearson") ^ 2)/ tmp$df.resid,
                       NA)
        c(w, y, z, Disp)
    }
}

#-------------------------------------------------------------------------------
#' Perform Model Selection Using Zero-Inflated Generalized Linear Models (GLMs)
#'
#' \code{ZI.Select} generates a vector that contains an \code{\link{AIC}}, an
#'   \code{\link[MuMIn]{AICc}}, a \code{\link{BIC}}, and a dispersion parameter
#'   estimate.  Intent is to use this function in a loop to extract model
#'   information criterion for a range of model formulations, facilitating the
#'   automated selection of *best* model.
#' @param formula character string representing a potential model formula
#' @param data data frame containing columns wiht names corresponding to the
#'   values in the formula
#' @param dist see help page for \code{\link[pscl]{zeroinfl}}
#' @param link see help page for \code{\link[pscl]{zeroinfl}}
#' @return vector that contains the \code{\link{AIC}},
#'   \code{\link[MuMIn]{AICc}}, \code{\link{BIC}} and dispersion parameter
#'   estimate for a given model structure
#' @family Model Selection
#' @export
ZI.Select <- function(formula, data, dist, link) {
    tmp <- try(pscl::zeroinfl(as.formula(formula), data = data,
                              dist = dist, link = link))
    w <- ifelse(class(try(AIC(tmp))) == "numeric", AIC(tmp), NA)
    y <- ifelse(class(try(MuMIn::AICc(tmp))) == "numeric", MuMIn::AICc(tmp), NA)
    z <- ifelse(class(try(AIC(tmp))) == "numeric", BIC(tmp), NA)
    Disp <- ifelse(class(try(AIC(tmp))) == "numeric",
                   sum(resid(tmp, type = "pearson") ^ 2)/ tmp$df.resid,
                   NA)
    c(w, y, z, Disp)
}

#-------------------------------------------------------------------------------
#' Perform Model Selection Using Generalized Linear Mixed Models (GLMMs)
#'
#' \code{GLMM.Select} generates a vector that contains an \code{\link{AIC}},
#' an \code{\link[MuMIn]{AICc}}, a \code{\link{BIC}}, and a dispersion
#' parameter estimate. Intent is to use this function in a loop to extract
#' model information criterion for a range of model formulations, facilitating
#' the automated selection of *best* model.
#' @param formula combined fixed and random effects formula, following lme4
#' syntax.
#' @param data data frame (tibbles are OK) containing model variables. Not
#' required, but strongly recommended; if \code{data} is not specified,
#' downstream methods such as prediction with new data
#' (\code{predict(fitted_model, newdata = ...)}) will fail. If it is necessary
#' to call \code{glmmTMB} with model variables taken from the environment
#' rather than from a data frame, specifying \code{data=NULL} will suppress
#' the warning message.
#' @param family a family function, a character string naming a family
#' function, or the result of a call to a family function (variance/link
#' function) information. See \code{\link{family}} for a generic discussion of
#' families or \code{\link[glmmTMB]{family_glmmTMB}} for details of
#' \code{glmmTMB}-specific families.
#' @param ziformula a \emph{one-sided} (i.e., no response variable) formula for
#' zero-inflation combining fixed and random effects: the default \code{~0}
#' specifies no zero-inflation. Specifying \code{~.} sets the zero-inflation
#' formula identical to the right-hand side of \code{formula} (i.e., the
#' conditional effects formula); terms can also be added or subtracted.
#' \strong{When using \code{~.} as the zero-inflation formula in models where
#' the conditional effects formula contains an offset term, the offset term
#' will automatically be dropped}. The zero-inflation model uses a logit link.
#' @param dispformula a \emph{one-sided} formula for dispersion containing only
#' fixed effects: the default \code{~1} specifies the standard dispersion given
#' any family. The argument is ignored for families that do not have a
#' dispersion parameter. For an explanation of the dispersion parameter for
#' each family, see \code{\link[glmmTMB]{sigma}}. The dispersion model uses a log link.
#' In Gaussian mixed models, \code{dispformula=~0} fixes the residual variance
#' to be 0 (actually a small non-zero value), forcing variance into the random
#' effects. The precise value can be controlled via
#' \code{control=glmmTMBControl(zero_dispval=...)}; the default value is
#' \code{sqrt(.Machine$double.eps)}.
#' @param weights weights, as in \code{glm}. Not automatically scaled to have
#' sum 1.
#' @param offset offset for conditional model (only).
#' @param contrasts an optional list, e.g., \code{list(fac1="contr.sum")}. See
#' the \code{contrasts.arg} of \code{\link{model.matrix.default}}.
#' @param na.action a function that specifies how to handle observations
#' containing \code{NA}s.  The default action (\code{na.omit}, inherited from
#' the 'factory fresh' value of \code{getOption("na.action")}) strips any
#' observations with any missing values in any variables. Using
#' \code{na.action = na.exclude} will similarly drop observations with missing
#' values while fitting the model, but will fill in \code{NA} values for the
#' predicted and residual values for cases that were excluded during the
#' fitting process because of missingness.
#' @param se whether to return standard errors.
#' @param verbose whether progress indication should be printed to the console.
#' @param doFit whether to fit the full model, or (if FALSE) return the
#' pre-processed data and parameter objects, without fitting the model.
#' @param control control parameters, see \code{\link[glmmTMB]{glmmTMBControl}}.
#' @param REML whether to use REML estimation rather than maximum likelihood.
#' @param start starting values, expressed as a list with possible components
#' \code{beta}, \code{betazi}, \code{betad} (fixed-effect parameters for
#' conditional, zero-inflation, dispersion models); \code{b}, \code{bzi}
#' (conditional modes for conditional and zero-inflation models); \code{theta},
#' \code{thetazi} (random-effect parameters, on the standard deviation/Cholesky
#' scale, for conditional and z-i models); \code{thetaf} (extra family
#' parameters, e.g., shape for Tweedie models).
#' @param map a list specifying which parameter values should be fixed to a
#' constant value rather than estimated. \code{map} should be a named list
#' containing factors corresponding to a subset of the internal parameter
#' names (see \code{start} parameter). Distinct factor values are fitted as
#' separate parameter values, \code{NA} values are held fixed: e.g.,
#' \code{map=list(beta=factor(c(1,2,3,NA)))} would fit the first three
#' fixed-effect parameters of the conditional model and fix the fourth
#' parameter to its starting value. In general, users will probably want to
#' use \code{start} to specify non-default starting values for fixed
#' parameters. See \code{\link[TMB]{MakeADFun}} for more details.
#' @param sparseX a named logical vector containing (possibly) elements named
#' "cond", "zi", "disp" to indicate whether fixed-effect model matrices for
#' particular model components should be generated as sparse matrices, e.g.
#' \code{c(cond=TRUE)}. Default is all \code{FALSE}
#' @importFrom stats gaussian binomial poisson nlminb as.formula terms model.weights
#' @importFrom lme4 subbars findbars mkReTrms nobars
#' @importFrom Matrix t
#' @importFrom TMB MakeADFun sdreport
#' @importFrom glmmTMB glmmTMB glmmTMBControl
#' @importFrom stats na.fail sigma
#' @details
#' See \code{\link[glmmTMB]{glmmTMB}}
#' @note
#' For more information about the \pkg{glmmTMB} package, see Brooks et al.
#' (2017) and the \code{vignette(package="glmmTMB")} collection. For the
#' underlying \pkg{TMB} package that performs the model estimation, see
#' Kristensen et al. (2016).
#' @references
#' Brooks, M. E., Kristensen, K., van Benthem, K. J., Magnusson, A., Berg,
#' C. W., Nielsen, A., Skaug, H. J., Mächler, M. and Bolker, B. M. (2017).
#' glmmTMB balances speed and flexibility among packages for zero-inflated
#' generalized linear mixed modeling. \emph{The R Journal}, \bold{9}(2),
#' 378--400.
#'
#' Kristensen, K., Nielsen, A., Berg, C. W., Skaug, H. and Bell, B. (2016).
#' TMB: Automatic differentiation and Laplace approximation. \emph{Journal of
#' Statistical Software}, \bold{70}, 1--21.
#'
#' Millar, R. B. (2011). \emph{Maximum Likelihood Estimation and Inference:
#' With Examples in R, SAS and ADMB.} Wiley, New York.
#' @export
#'
GLMM.Select <- function(formula, data, family = gaussian(), ziformula = ~0,
                        dispformula = ~1, weights = NULL, offset = NULL,
                        contrasts = NULL, na.action = na.fail, se = TRUE,
                        verbose = FALSE, doFit = TRUE,
                        control = glmmTMBControl(profile = TRUE, collect = TRUE),
                        REML = FALSE, start = NULL,
                        map = NULL, sparseX = NULL) {
    tmp <- try(glmmTMB(formula = as.formula(formula), data = data,
                       family = family, ziformula = ziformula,
                       dispformula = dispformula, weights = weights,
                       offset = offset, contrasts = contrasts,
                       na.action = na.action, se = se, verbose = verbose,
                       doFit = doFit, control = control, REML = REML,
                       start = start, map = map, sparseX = sparseX))
    w <- ifelse(class(try(AIC(tmp))) == "numeric", AIC(tmp),
                NA)
    y <- ifelse(class(try(MuMIn::AICc(tmp))) == "numeric",
                MuMIn::AICc(tmp), NA)
    z <- ifelse(class(try(BIC(tmp))) == "numeric", BIC(tmp),
                NA)
    Disp <- ifelse(class(try(AIC(tmp))) == "numeric",
                   sigma(tmp), NA)
    Fit <- tmp
    list(AIC = w, AICc = y, BIC = z, Dispersion = Disp, Fitted = Fit)
}

