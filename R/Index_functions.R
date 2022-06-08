#' Provide Names for a Series of Models
#'
#' \code{Model.Names} generates a character vector that can be used to represent
#'   names of individual model runs.  Intent is to quickly generate a set of
#'   unique model names that can be used to store specific model runs.
#' @param dat vector used to indicate how individual variables in the models are
#'   included in a specific model run.  Typically, this would represent a
#'   numeric vector the length of the maximum number of co-variates being
#'   considered during model development.  In this vector, a value of 0 would
#'   represent that variable not being included in that particular model run.
#'   Generally, a value of 1 represents that a covariate (whether continuous or
#'   discrete) is included in the model.  Any value >1 (e.g., 3) represents that
#'   the co-variate was included in the model as a continous variable and
#'   modeled using multiple polynomials to the order suggested by the number
#'   (e.g., 3 suggests we modeled the co-variate as a 3rd order polynomial).
#' @param names character vector representing the maximum number of co-variates
#'   being considered during model development.  For "two-stage' type models
#'   (e.g., zero-inflated models),their may be duplicate names as the same
#'   co-variate may be included in both stages of the model.  In the "two-stage"
#'   scenario all co-variates included in the 1st stage should be supplied
#'   first, followed by all co-variates included in the 2nd stage.
#' @param model character string denoting the type of model being specified.
#'   Currently, model can take 1 of 2 values: "ZINB" (default) or "GLM".  If
#'   another value is provided for model, function will stop, reporting the
#'   error "An inappropriate value was provided for 'model'".
#' @param structure value that represents the model structure (default =
#'   "ZINB").  The investigator is free to specify this in any way they deem
#'   appropriate.
#' @param count.only logical indicating whether you are only modeling the 'count
#'   sub-model' of a zero-inflated model.  If 'TRUE' (default), you are treating
#'   the 'zero-inflation sub-model' as an intercept only model.  If 'FALSE' you
#'   are specifying you are modeling both the 'count' and 'zero-inflation'
#'   sub-models of a zero-inflated model using co-variates.
#' @param count.n numeric value indicating the position of the last 'names'
#'   value that represents co-variates included in the 1st stage of a two-stage
#'   model.
#' @return Character string representing a model name which represents the proposed model structure.
#' @examples
#' # GLM Model
#' Model.Names(c(1,0,2,3), names = c("Yr", "D", "L", "T"), model = "GLM",
#' structure = "Pois")
#' # Zero-Inflated Model - Count Sub-model Only
#' Model.Names(c(1,0,2,3), names = c("Yr", "D", "L", "T"), count.n = 4)
#' # Zero-Inflated Model - Both Sub-Models
#' Model.Names(c(1,0,2,3,0,1,0,0), names = c("Yr","D","L","T","Yr","D","L","T"),
#' count.only = FALSE, count.n = 4)
#' @export
Model.Names <- function(dat, names, model = "ZINB", structure = "ZINB",
                        count.only = TRUE, count.n) {
    if (model == "GLM") {
        paste0(structure, paste0(names, dat, collapse = ""))
    } else {
        if (model == "ZINB") {
            if (count.only == TRUE) {
                paste0(structure, "_Count_", paste0(names[1 : count.n],
                                                    dat[1 : count.n],
                                                    collapse = ""),
                       "_ZI_Intercept_")
            } else {
                if (count.only == FALSE) {
                    paste0(structure,
                           "_Count_",
                           paste0(names[1 : count.n],
                                  dat[1 : count.n],
                                  collapse = ""),
                           "_ZI",
                           paste0(names[(count.n + 1) : length(dat)],
                                  dat[(count.n + 1) : length(dat)],
                                  collapse = ""))

                } else {
                    stop("'count.only' must take the value 'TRUE' or 'FALSE'")
                }
            }
        } else {
            stop("An inappropriate value was provided for 'model'")
        }
    }
}

# ------------------------------------------------------------------------------
#' Develop formulas for GAM models
#'
#' \code{GAM.Formulas} generates a character vector that can be used to
#'   represent formulas for use in individual GAM model runs.  The intent is to
#'   quickly generate a GAM model formula that incorporates information from
#'   several potential covariates, the type of covariate utilized, and the
#'   spline "basis" you specify.
#' @param Response Character string specifying the name of the response
#'   variable.  This should match a specific column name in your data frame.
#' @param names character vector representing the names of unique
#'   co-variates you want to include in your GAM model.  These should match
#'   specific column names in your data frame.
#' @param VarType character string denoting the type of co-variate (continuous
#'   ("C") or discrete ("F")) for each co-varaite listed in *predictors*.  This
#'   character string must be the same length as the character string specified
#'   by *predictors*.  If a value other than "C" or "F" is provided, the
#'   function will return an error.
#' @param basis a two letter character string indicating the (penalized)
#'   smoothing basis to use (e.g., "tp" for thin plate regression spline, "cr"
#'   for cubic regression spline).  see \code{\link[mgcv]{s}} and
#'   \code{\link[mgcv]{smooth.terms}} for more information.
#' @return Character string representing a potential GAM model formulat which
#'   represents the proposed model structure defined by the function inputs.
#' @family Model Formulas
#' @examples
#' data(iris)
#' form <- GAM.Formulas(Response = "Sepal.Length", names =
#'   c("Petal.Length", "Petal.Width", "Sepal.Width", "Species"), VarType =
#'   c(rep("C", 3), "F"), basis = "tp")
#' form
#' @export
GAM.Formulas <- function(Response, names, VarType, basis = "cr") {
    tmp <- names
    for(j in 1: length(tmp)) {
        if(VarType[j] == "F") {
           tmp[j] <- tmp[j]
        } else {
            if(VarType[j] == "C") {
                tmp[j] <- paste0("s(", tmp[j],", bs = '",
                                 basis, "')")
            } else {
                stop("Invalid input provided for 'VarType'")
            }
        }
    }
    paste0(Response, " ~ ", paste0(tmp, collapse = " + "))
}
#-------------------------------------------------------------------------------
#' Develop formulas for GLM models
#'
#' \code{GLM.Formulas} generates a character vector that can be used to
#'   represent formulas for use in individual GLM model runs.  The intent is to
#'   quickly generate a GLM model formula that incorporates information from
#'   several potential covariates and the type of covariate utilized.
#' @param Response Character string specifying the name of the response
#'   variable.  This should match a specific column name in your data frame.
#' @param Predictors numeric vector representing 'value' for each co-variate
#'   included in the model, where the definition of 'value' is different
#'   depending on the type of covariate in question.  For discrete covariates,
#'   a value of 0 indicates that co-variate should not be included in the
#'   model; all other values indicate the variable will be included.  For
#'   continuous covariates, a value of 0 indicates that co-variate should not
#'   be included in the model; positive values indicate the polynomial order
#'   for that co-variate.
#' @param offset character vector representing an "offset" variable to be used
#'   in the formula.
#' @param VarType character string denoting the type of co-variate (continuous
#'   ("C") or discrete ("F")) for each co-varaite listed in *predictors*.  This
#'   character string must be the same length as the character string specified
#'   by *predictors*.  If a value other than "C" or "F" is provided, the
#'   function will return an error.
#' @param names character vector representing the names of unique
#'   co-variates you want to include in your GLM model.  These should match
#'   specific column names in your data frame.
#' @return Character string representing a potential GLM model formula which
#'   represents the proposed model structure defined by the function inputs.
#' @family Model Formulas
#' @examples
#' data(iris)
#' form <- GLM.Formulas(Response = "Sepal.Length", Predictors = c(1, 2, 0, 3),
#' names = c("Petal.Length", "Petal.Width", "Sepal.Width", "Species"), VarType =
#' c(rep("C", 3), "F"))
#' form
#' m1 <- glm(as.formula(form), data = iris, family="gaussian")
#' summary(m1)
#' plot(m1, all.terms = TRUE, seWithMean = TRUE)
#' @export
GLM.Formulas <- function(Predictors, Response, offset = NULL, VarType, names) {
    tmp <- Predictors
    for(j in 1:length(tmp)) {
        if(VarType[j] == "F") {
            tmp[j] <- ifelse(tmp[j] == 0, NA, names[j])
        } else {
            if(VarType[j] == "C") {
                tmp[j] <- ifelse(tmp[j] == 0,
                                 NA,
                                 paste0("poly(", names[j],
                                        ", ",
                                        tmp[j],
                                        ", raw = TRUE)"))
            } else {
                stop("Invalid value for 'vartype' provided")
            }
        }
    }
    if(is.null(offset) == FALSE) {
        offset <- paste0("offset(",offset,") + ")
    } else {
        offset <- "NA + "
    }
    tmp <- paste0(Response," ~ ", offset, paste(tmp, collapse = " + "))
    tmp <- gsub("NA + ","", tmp, fixed = TRUE)
    gsub(" + NA","", tmp, fixed = TRUE)
}
#-------------------------------------------------------------------------------
#' Develop formulas for Zero-Inflation models
#'
#' \code{ZI.Formulas} generates a character vector that can be used to
#'   represent formulas for use in individual zero-inflation model runs.  The
#'   intent is to quickly generate a zero-inflation model formula that
#'   incorporates information from several potential covariates and the type of
#'   covariate utilized.  The same or different co-variates can be included in
#'   the different sub-models
#' @param Response Character string specifying the name of the response
#'   variable.  This should match a specific column name in your data frame.
#' @param Predictors numeric vector representing 'value' for each co-variate
#'   included in the model, where the definition of 'value' is different
#'   depending on the type of covariate in question.  For discrete covariates,
#'   a value of 0 indicates that co-variate should not be included in the
#'   model; all other values indicate the variable will be included.  For
#'   continuous covariates, a value of 0 indicates that co-variate should not
#'   be included in the model; positive values indicate the polynomial order
#'   for that co-variate.
#' @param offset character vector representing an "offset" variable to be used
#'   in the formula.
#' @param VarType character string denoting the type of co-variate (continuous
#'   ("C") or discrete ("F")) for each co-varaite listed in *predictors*.  This
#'   character string must be the same length as the character string specified
#'   by *predictors*.  If a value other than "C" or "F" is provided, the
#'   function will return an error.
#' @param names character vector representing the names of unique
#'   co-variates you want to include in your zero-inflation model.  These
#'   should match specific column names in your data frame.
#' @param count.n integer representing the number of covariates included in the
#'   *count* sub-model of the overall zero-inflation model.  If equal to (or
#'   greater) the length of *names* then all covariates are only included in
#'   the *count* sub-model.  If less than the length of *names*, then those
#'   co-variates at positions greater than *count.n* are included in the
#'   zero-inflation sub-model.
#' @return Character string representing a potential zero-inflation model
#'   formula which  represents the proposed model structure defined by the
#'   function inputs.
#' @family Model Formulas
#' @examples
#' data(iris)
#' form <- ZI.Formulas(Response = "Sepal.Length", Predictors = c(1, 2, 0, 3, 1,
#' 0, 0, 0), names = rep(c("Petal.Length", "Petal.Width", "Sepal.Width",
#' "Species"), 2), VarType = rep(c(rep("C", 3), "F"), 2), count.n = 4)
#' form
#' @export
ZI.Formulas <- function(Predictors, Response, offset = NULL, VarType, names,
                        count.n) {
    Count <- Predictors[1 : count.n]
    for(j in 1 : count.n) {
        if(VarType[j] == "F") {
            Count[j] <- ifelse(Count[j] == 0, NA, names[j])
        } else {
            if(VarType[j] == "C") {
                Count[j] <- ifelse(Count[j] == 0,
                                   NA,
                                   paste0("poly(",
                                          names[j],
                                          ", ", Count[j],
                                          ", raw = TRUE)"))
            } else {
                stop("Invalid value for 'VarType' provided")
            }
        }
    }
    if (count.n >= length(Predictors)) {
        ZI <- " 1"
    } else {
        ZI <- Predictors[(count.n + 1) : length(Predictors)]
        for(j in (count.n + 1) : length(Predictors)) {
            if(VarType[j] == "F") {
                ZI[j - count.n] <- ifelse(ZI[j - count.n] == 0,
                                          NA,
                                          names[j])
            } else {
                if(VarType[j] == "C") {
                    ZI[j - count.n] <- ifelse(ZI[j - count.n] ==0,
                                              NA,
                                              paste0("poly(",
                                                     names[j],
                                                     ", ",
                                                     ZI[j - count.n],
                                                     ", raw = TRUE)"))
                } else {
                    stop("Invalid value for 'VarType' provided")
                }
            }
        }
    }
    Y.Model <- paste0(Response," ~ ")
    C.Model <- paste0(Count, collapse = " + ")
    Z.Model <- paste0(ZI, collapse = " + ")
    O.Model <- paste0("offset(", offset, ") + ")
    if(is.null(offset) == TRUE) {
        tmp <- paste0(Y.Model, C.Model, " | ", Z.Model)
        tmp <- gsub("NA + ","", tmp, fixed = TRUE)
        gsub("+ NA", "", tmp, fixed = TRUE)
    } else {
        tmp <- paste0(Y.Model, O.Model, C.Model, " | ", O.Model, Z.Model)
        gsub("+ NA", "", tmp, fixed = TRUE)
    }
}
#-------------------------------------------------------------------------------
#' Develop formulas for a variety of different types of modles
#'
#' \code{Model.Formulas} generates a character vector that can be used to
#'   represent formulas for use in individual GAM, GLM, or zero-inflation
#'   model runs.  In essence it is a wrapper for the individual functions
#'   \code{\link{GAM.Formulas}}, \code{\link{GLM.Formulas}}, and
#'   \code{\link{ZI.Formulas}}.  The intent is to quickly generate model
#'   formulas that incorporates information from several potential covariates
#'   and the type of covariate utilized.
#' @param model character string specifiying the type of model formula to be
#'   constructed.  Currently only "GAM", "GLM", and "ZI" are acceptable values,
#'   with all others giving an error.
#' @param Response see help pages for \code{\link{GAM.Formulas}},
#'   \code{\link{GLM.Formulas}}, and \code{\link{ZI.Formulas}} for description.
#' @param Predictors see help pages for \\code{\link{GLM.Formulas}} and
#'   \code{\link{ZI.Formulas}} for description.
#' @param offset see help pages for \code{\link{GLM.Formulas}} and
#'   \code{\link{ZI.Formulas}} for description.
#' @param VarType see help pages for \code{\link{GAM.Formulas}},
#'   \code{\link{GLM.Formulas}}, and \code{\link{ZI.Formulas}} for description.
#' @param names see help pages for \code{\link{GAM.Formulas}},
#'   \code{\link{GLM.Formulas}}, and \code{\link{ZI.Formulas}} for description.
#' @param count.n see help page for \code{\link{ZI.Formulas}} for description.
#' @param basis see help pages for \code{\link{GAM.Formulas}} for description.
#' @return Character string representing a potential model formula.
#' @family Model Formulas
#' @examples
#' data(iris)
#' # GAM Model
#' form.GAM <- GAM.Formulas(Response = "Sepal.Length", names =
#'   c("Petal.Length", "Petal.Width", "Sepal.Width", "Species"), VarType =
#'   c(rep("C", 3), "F"), basis = "tp")
#' form.GAM
#' # GLM Model
#' form.GLM <- GLM.Formulas(Response = "Sepal.Length", Predictors =
#' c(1, 2, 0, 3), names = c("Petal.Length", "Petal.Width", "Sepal.Width",
#' "Species"), VarType = c(rep("C", 3), "F"))
#' form.GLM
#' # Zero-Inflation Model
#' form.ZI <- ZI.Formulas(Response = "Sepal.Length", Predictors =
#' c(1, 2, 0, 3, 1, 0, 0, 0), names = rep(c("Petal.Length", "Petal.Width",
#' "Sepal.Width", "Species"), 2), VarType = rep(c(rep("C", 3), "F"), 2),
#' count.n = 4)
#' form.ZI
#' @export
Model.Formulas <- function(Predictors, Response, model, offset = NULL, VarType,
                           names, count.n, basis) {
    models <- c("GAM", "GLM", "ZI")
    approp.basis <- c("tp", "ts", "ds", "cr", "cs", "cc", "sos", "ps", "cp",
                      "re", "mrf", "gp", "so")
    if(is.null(Response) == TRUE)
        stop("'Response' variable must be provided")
    if(is.null(names) == TRUE)
        stop("'names' for variable(s) must be provided")
    if(is.null(VarType) == TRUE)
        stop("'VarType(s)' must be specified for each predictor variable")
    if(model %in% models == FALSE)
        stop("An inappropriate model 'type' was provided to the function")
    if(model == "GAM") {
        if(basis %in% approp.basis == FALSE) {
            stop("An inappropriate basis for the smoother was provided: see help for 'smooth.terms' in the package 'mgcv'")
        } else return(GAM.Formulas(names, Response, VarType, basis))
    }
    if(model == "GLM") {
        if(length(names) != length(Predictors)) {
            stop("A list of variable 'names' the same length as the length of predictors must be provided")
        } else return(GLM.Formulas(Predictors, Response, offset, VarType,
                                   names))
    }
    if(model == "ZI") {
        if(length(names) != length(Predictors)) {
            stop("A list of variable 'names' the same length as the length of predictors must be provided")
        } else {
           if(is.null(count.n) == TRUE) {
               stop("A value for 'count.n' must be provided for a zero-inflation model formula")
           } else
               return(ZI.Formulas(Predictors, Response, offset, VarType, names,
                           count.n))
        }
    }
}

#-------------------------------------------------------------------------------
#' Summarize results of a standardized relative abundance index
#'
#' \code{Index.Summary}
#' @param Y numeric vector representing response variable.  In the case of a
#'   relative abundance index, this would represent your CPUE measure
#' @param X numeric vector representing the primary covariate of interest.  As
#'   such, this is the variable you are wanting to get a predicted "effect" for
#'   over a range of values this variable can take.
#' @param grid prediction grid.  This most commonly would be created via a call
#'   to \code{Index.Grid} if used in a relative abundance index development
#'   question
#' @param models character vector of named models for which you want to get an
#'   index summary for
#' @param model.name character vector defining model type.  All values of the
#'   character vector must take one of five defined values:
#'   \describe{
#'     \item{Poisson}{identifies use of a "Poisson GLM" for the development of
#'       the model}
#'     \item{NB}{identifies use of a "Negative Binomial GLM" for the development
#'       of the model}
#'     \item{ZIP}{identifies use of a "Zero-Inflated Poisson GLM" for the
#'       development of the model}
#'     \item{ZINB}{identifies use of a "Zero-Inflated Negative Binomial GLM" for
#'       the development of the model}
#'     \item{Nominal}{identifies the calculation of a nominal index}
#'   }
#' @return list the same length as \code{models} that contains summary
#'   statistics for each model explored
#' @family Model Evaluation
#' @export
Index.Summary <- function(Y, X, grid, models,
                          model.name = c("Poisson", "NB", "ZIP", "ZINB",
                                         "Nominal")) {
    model.list.names <- model.name
    model.list <- sapply(model.list.names, function(x) NULL)
    model.results.list <- sapply(model.list.names, function(x) NULL)
    for(i in 1 : length(models)) {
        if(model.name[i] == "Poisson") {
            Results <- predict(models[[i]], grid, type = "response",
                               se.fit = TRUE)
            model.list[[i]] <- Results
            Year = unique(grid[,1])
            Raw = aggregate(Results$fit ~ as.vector(grid[, 1]), FUN = mean)[, 2]
            Normalized = Raw / mean(Raw)
            Centered = Raw - mean(Raw)
            Scaled = scale(Raw)
            model.results.list[[i]] <- data.frame(Year = Year, Raw = Raw,
                                                  Normalized = Normalized,
                                                  Centered = Centered,
                                                  Scaled = Scaled)
        }
        if(model.name[i] == "NB") {
            Results <- predict(models[[i]], grid, type = "response",
                               se.fit = TRUE)
            model.list[[i]] <- Results
            Year = unique(grid[, 1])
            Raw = aggregate(Results$fit ~ grid[, 1], FUN = mean)[, 2]
            Normalized = Raw / mean(Raw)
            Centered = Raw - mean(Raw)
            Scaled = scale(Raw)
            model.results.list[[i]] <- data.frame(Year = Year, Raw = Raw,
                                                  Normalized = Normalized,
                                                  Centered = Centered,
                                                  Scaled = Scaled)
        }
        if(model.name[i] == "ZIP") {
            Results <- predict(models[[i]], grid, type = "response")
            model.list[[i]] <- Results
            Year = unique(grid[, 1])
            Raw = aggregate(Results ~ grid[, 1], FUN = mean)[, 2]
            Normalized = Raw / mean(Raw)
            Centered = Raw - mean(Raw)
            Scaled = scale(Raw)
            model.results.list[[i]] <- data.frame(Year = Year, Raw = Raw,
                                                  Normalized = Normalized,
                                                  Centered = Centered,
                                                  Scaled = Scaled)
        }
        if(model.name[i] == "ZINB") {
            Results <- predict(models[[i]], grid, type = "response")
            model.list[[i]] <- Results
            Year = unique(grid[, 1])
            Raw = aggregate(Results ~ grid[, 1], FUN = mean)[,2]
            Normalized = Raw / mean(Raw)
            Centered = Raw - mean(Raw)
            Scaled = scale(Raw)
            model.results.list[[i]] <- data.frame(Year = Year, Raw = Raw,
                                                  Normalized = Normalized,
                                                  Centered = Centered,
                                                  Scaled = Scaled)
        }
        if(model.name[i] == "Nominal") {
            model.list[[i]] <- NULL
            Year = unique(grid[,1])
            Raw = aggregate(Y ~ X, FUN = mean)[, 2]
            Normalized = Raw / mean(Raw)
            CV = (aggregate(Y ~ X, FUN = sd)[, 2] /
                      sqrt(aggregate(Y ~ X, FUN = length)[, 2])) / Raw
            Centered = Raw - mean(Raw)
            Scaled = scale(Raw)
            model.results.list[[i]] <- data.frame(Year = Year, Raw = Raw,
                                                  Normalized = Normalized,
                                                  CV = CV, Centered = Centered,
                                                  Scaled = Scaled)
        }
    }
    model.results <- model.list
    model.summary <- model.results.list
    Results <<- list(Models = model.results, Summary = model.summary)
}

#-------------------------------------------------------------------------------
#' Develop "grid" over which to predict CPUE from an index model based on
#'   covariate values
#'
#' \code{Index.Grid}
#' @param X numeric vector representing the primary covariate of interest.  As
#'   such, this is the variable you are wanting to get a predicted "effect" for
#'   over a range of values this variable can take.
#' @param X.type character vector defining structure of the \code{X} variable of
#'   interest.  See \code{type} for a description of the possible \code{X.type}s
#' @param X.name character string defining the name of the \code{X} variable of
#'   interest
#' @param covariates numeric list containing additional covariates included in
#'   the final models being considered
#' @param type character vector defining structure of each covariate.  It must
#'   be the same length as the \code{covariates} list.  Each covariate must be
#'   defined as one of four values:
#'   \describe{
#'     \item{factor}{discrete covariate where data fall in discrete bins}
#'     \item{mean}{any variable for which you want to use the mean value in the
#'       model prediction frame.  This would be commonly used to represent an
#'       offset variable or potentially continous covariates}
#'     \item{numeric}{continuous covariate, i.e., data fall along a continuous
#'       scale}
#'     \item{continuous}{continuous covariate, i.e., data fall along a
#'       continuous scale}
#'     \item{discrete}{numeric value representing a discrete value for which you
#'       want to use in a prediction frame for a given covariete.  For example,
#'       you may want to choose a specific depth or latitude}
#'     \item{median}{any variable for which you want to use the median value in
#'       the model prediction frame.  This would be commonly used to represent
#'       an offset variable or potentially continous covariates}
#'     \item{mode}{any variable for which you want to use the modal value in the
#'       model prediction frame.  This would be commonly used to represent an
#'       offset variable or potentially continous covariates}
#'   }
#' @param bin.length numeric value representing the number of discrete points of
#'   each continuous covariate to use in the prediction frame.  Recommendation
#'   is to use an odd number such that the median of the range of the contiuous
#'   variable is used.  Function used the value to choose n = \code{bin.length}
#'   equally spaced points from the minimum of each continuous covariate to the
#'   maximum of each continuous covariate
#' @param cont.length numeric value representing the number of discrete points
#'   of each continuous covariate to use in the prediction frame.
#'   Recommendation is to use an odd number such that the median of the range of
#'   the contiuous variable is used.  Function used the value to choose n =
#'   \code{cont.length} equally spaced points from the minimum of each
#'   continuous covariate to the maximum of each continuous covariate. It is
#'   analogous to \code{bin.length}, just allows for additional flexibility of
#'   having more precision for the primary variable of interest
#' @return list the same length as \code{models} that contains summary
#'   statistics for each model explored
#' @family Model Evaluation
#' @examples
#' X <- factor(sample(c(seq(1990, 2015, 1)), size = 10000, replace = TRUE))
#' df <- data.frame(depth = rnorm(n = 10000, mean = 0, sd = 1), temp = rnorm(
#' n = 10000, mean = 0, sd = 0.5), lat = rnorm(n = 10000, mean =0, sd = 3))
#' Index.Grid(X = X, X.type = "factor", X.name = "Year", covariates = df,
#' type = rep("mean", 3))
#' Index.Grid(X = X, X.type = "factor", X.name = "Year", covariates = df, type =
#' rep("numeric", 3), bin.length = 3)
#' @export
Index.Grid <- function(X, X.type, X.name, covariates, type, bin.length=17,
                       cont.length = 51) {
    tmp.list.names <- c(X.name, names(covariates))
    tmp.list <- sapply(tmp.list.names, function(x) NULL)

    if(X.type == "factor")
        tmp.list[[1]] <- levels(X)
    if(X.type == "continuous")
        tmp.list[[1]] <- seq(min(X, na.rm = TRUE), max(X, na.rm = TRUE),
                             length = cont.length)

    for(i in 1 : length(covariates)) {
        if(type[i] == "factor")
            tmp.list[[i + 1]] <- levels(covariates[[i]])
        if(type[i] == "mean")
            tmp.list[[i + 1]] <- mean(covariates[[i]])
        if(type[i] == "numeric")
            tmp.list[[i + 1]] <- seq(min(covariates[[i]]), max(covariates[[i]]),
                                     length = bin.length)
        if(type[i] == "continuous")
            tmp.list[[i + 1]] <- seq(min(covariates[[i]]), max(covariates[[i]]),
                                     length = cont.length)
        if(type[i] == "discrete")
            tmp.list[[i + 1]] <- covariates[[i]]
        if(type[i] == "median")
            tmp.list[[i + 1]] <- median(covariates[[i]])
        if(type[i] == "mode")
            tmp.list[[i + 1]] <- mode(covariates[[i]])
    }
    tmp <- expand.grid(tmp.list)
    return(tmp)
}

#-------------------------------------------------------------------------------
#' Calculate the dispersion estimate from a GLM model
#' \code{disp}
#' @param mod Name of a fitted model of various types from which Pearson
#'   residuals and a residual degrees of freedom can be extracted from the
#'   model object
#' @family Model Evaluation
#' @examples
#' X <- seq(1, 100)
#' Y <- 100 + 0.2 * X + rnorm(1, mean = 0, sd = 5)
#' mod <- glm(Y ~ X)
#' disp(mod) # Dispersion parameter estimated using functio
#' sigma(mod)^2 # Dispersion parameter estimated using the calculated sigma
#' # paramter of the model object
#' @export
disp = function(mod = NA){
    E = residuals(mod, type='pearson')
    d = sum(E^2)/mod$df.resid
    return(d)
}
