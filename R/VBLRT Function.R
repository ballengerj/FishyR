#' Von Bertalanffy Growth Function Likelihood Ratio Test
#'
#' \code{VBLRT} constructs a series of von Bertalanffy growth function models
#'   looking at likelihood ratio tests to determine if there are differences
#'   between defined groups.  Here we only allow for a single grouping variable
#'   with two levels (i.e., sex (Male and Female)).
#' @param len
#' @param age
#' @param group
#' @param error
#' @param select
#' @param Linf
#' @param k
#' @param t0
#' @return Individual von Bertalanffy growth function model fits, based on data
#'   and specifications defined, along with results of likelihood ratio tests
#'   and information criterion comparing the fits of the different models.
#' @examples
#' @export
VBLRT <- function(len = NULL, age = NULL, group = NULL, error = 1, select = 1,
                  Linf = NULL, k = NULL, t0 = NULL) {
    # Defining some error messages
    if (is.null(len)) stop("len is missing")
    if (is.null(age)) stop("age is missing")
    if (is.null(group)) stop("group is missing")
    if (length(age) != length(len)) stop("Vectors of different lengths")
    if (nlevels(as.factor(group)) > 2) stop("Only two groups allowed in data")
    if (select == 2 & (is.null(Linf) | is.null(k) | is.null(t0)))
        stop("User-specified values of Linf, K, and t0 are required")
    # Manipulationg "group" variable such that it has values of 0 and 1
    cat <- as.integer(as.factor(group)) - 1
    # Combining model data into a data frame, removing any non-complete cases
    x <- as.data.frame(cbind(len, age, cat))
    x <- na.exclude(x)
    # Setting initial wgt term to NA
    wgt <- NULL
    # Defining input parameter estimates
    if (select == 1) {
        m1xt <- NULL
        m2xt <- NULL
        mbxt <- NULL
        g1 <- with(x, aggregate(len ~ cat + age, FUN = mean))
        m1 <- g1[g1[, 1] == 0, ]
        m2 <- g1[g1[, 1] == 1, ]
        m1t <- m1[, c(2:3)]
        m1t[, 1] <- m1t[, 1] - 1
        m2t <- m2[, c(2:3)]
        m2t[, 1] <- m2t[, 1] - 1
        m1xt <- merge(m1, m1t, by.x = c("Group.2"), by.y = c("Group.2"))
        out1 <- lm(m1xt[, 4] ~ m1xt[, 3])
        k1 <- abs(log(coef(out1)[2]))
        L1 <- -coef(out1)[1] / (coef(out1)[2] - 1)
        dx1 <- as.data.frame(cbind(L1 - m1$x, m1[, 2]))
        dx1 <- dx1[dx1[, 1] > 0, ]
        t01 <- (coef(lm(log(dx1[, 1]) ~ dx1[, 2]))[1] - log(L1))/k1
        m2xt <- merge(m2, m2t, by.x = c("Group.2"), by.y = c("Group.2"))
        out2 <- lm(m2xt[, 4] ~ m2xt[, 3])
        k2 <- abs(log(coef(out2)[2]))
        L2 <- -coef(out2)[1] / (coef(out2)[2] - 1)
        dx2 <- as.data.frame(cbind(L2 - m2$x, m2[, 2]))
        dx2 <- dx2[dx2[, 1] > 0, ]
        t02 <- (coef(lm(log(dx2[, 1]) ~ dx2[, 2]))[1] - log(L2))/k2
        g2 <- with(x, aggregate(len ~ round(age, 0), FUN = mean))
        mbt <- g2
        mbt[, 1] <- mbt[, 1] - 1
        mbxt <- merge(g2, mbt, by.x = c("Group.1"), by.y = c("Group.1"))
        outboth <- lm(mbxt[, 3] ~ mbxt[, 2])
        kboth <- abs(log(coef(outboth)[2]))
        Lboth <- -coef(outboth)[1] / (coef(outboth)[2] - 1)
        dxb <- as.data.frame(cbind(Lboth - g2$x, g2[, 1]))
        dxb <- dxb[dxb[, 1] > 0, ]
        t0b <- (coef(lm(log(dxb[, 1]) ~ dxb[, 2]))[1] - log(Lboth))/kboth
        Ld <- L2 - L1
        kd <- k2 - k1
        td <- t02 - t01
    }
    if (select == 2) {
        L1 <- L2 <- Lboth <- Linf
        k1 <- k2 <- kboth <- k
        t01 <- t02 <- t0b <- t0
        Ld <- 0
        kd <- 0
        td <- 0
    }
    # Defining Error Structures
    # Constant variance for all lijs
    if (error == 1)
        x$wgt <- 1
    # Constant variance for all mean leangths at age
    if (error == 2) {
        x <- with(x, aggregate(len ~ cat + age, FUN = mean))
        names(x) <- c("cat", "age", "len")
        x$wgt <- 1
    }
    # Variance of lij varies with age #
    if (error == 3) {
        d4 <- merge(with(x, aggregate(len ~ cat + age, FUN = mean)),
                    with(x, aggregate(len ~ cat + age, FUN = var)),
                    by.y = c("Group.1","Group.2"),
                    by.x = c("Group.1", "Group.2"))
        d4 <- merge(d4, with(x , aggregate(len ~ cat + age, FUN = length)),
                    by.y = c("Group.1", "Group.2"),
                    by.x = c("Group.1","Group.2"))
        names(d4) <- c("cat", "age", "len", "s2", "n")
        x <- d4
        x$wgt <- x$n / x$s2
        if (any(is.na(x$wgt)))
            stop("At least one age has a single length observation.\nNeed at least two observations to calculate variance.")
    }
    # Models
    # All parameters are different
    Ho <- nls(len ~ (Linf + ls * cat) *
                  (1 - exp(-(k + ks * cat) * (age -(t0 + ts * cat)))),
              data = x,
              weights = wgt,
              start = list(Linf = L1, ls = Ld, k = k1, ks = kd, t0 = t01,
                           ts = td),
              control = list(maxiter = 100, minFactor = 1 / 2 ^ 30,
                             warnOnly = TRUE))
    resid0 <- residuals(Ho)
    # All parameters are the same
    H1 <- nls(len ~ Linf * (1 - exp(-k * (age - t0))), data = x,
              weights = wgt,
              start = list(Linf = Lboth, k = kboth, t0 = t0b),
              control = list(maxiter = 100, minFactor = 1 / 2 ^ 30,
                             warnOnly = TRUE))
    resid1 <- residuals(H1)
    # The Linf parameters are the same #
    H2 <- nls(len ~ Linf * (1 - exp(-(k + ks*cat) * (age - (t0 + ts*cat)))),
              data = x, weights = wgt,
              start = list(Linf = Lboth, k = k1, ks = kd, t0 = t01, ts = td),
              control = list(maxiter = 100, minFactor = 1 / 2 ^ 30,
                             warnOnly = TRUE))
    resid2 <- residuals(H2)
    # The k parameters are the same #
    H3 <- nls(len ~ (Linf + ls * cat) * (1 - exp(-k * (age - (t0 + ts*cat)))),
              data = x, weights = wgt,
              start = list(Linf = L1, ls = Ld, k= kboth, t0 = t01, ts = td),
              control = list(maxiter = 100, minFactor = 1 / 2 ^ 30,
                             warnOnly = TRUE))
    resid3 <- residuals(H3)
    # The t0 parameters are the same #
    H4 <- nls(len ~ (Linf + ls * cat) *
                  (1 - exp(-(k + ks * cat) * (age - t0))),
              data = x, weights = wgt,
              start = list(Linf = L1, ls= Ld, k = k1, ks = kd, t0 = t0b),
              control = list(maxiter = 100, minFactor = 1 / 2 ^ 30,
                             warnOnly = TRUE))
    resid4 <- residuals(H4)
    # The Linf and k parameters are the same #
    H5 <- nls(len ~ Linf * (1 - exp(-k * (age - (t0 + ts * cat)))),
              data = x, weights = wgt,
              start = list(Linf = Lboth, k = kboth, t0 = t01, ts = td),
              control = list(maxiter = 100, minFactor=1 / 2 ^ 30,
                             warnOnly = TRUE))
    resid5 <- residuals(H5)
    # The Linf and t0 parameters are the same #
    H6 <- nls(len ~ Linf * (1 - exp(-(k + ks * cat) * (age - t0))),
              data = x, weights = wgt,
              start = list(Linf = Lboth, k = k1, ks = kd, t0 = t0b),
              control = list(maxiter = 100, minFactor = 1 / 2 ^ 30,
                             warnOnly = TRUE))
    resid6 <- residuals(H6)
    # The k and t0 parameters are the same #
    H7 <- nls(len ~ (Linf + ls * cat) * (1 - exp(-k * (age - t0))),
              data = x, weights = wgt,
              start = list(Linf = L1, ls = Ld, k = kboth, t0 = t0b),
              control = list(maxiter = 100, minFactor = 1 / 2 ^ 30,
                             warnOnly = TRUE))
    resid7 <- residuals(H7)
    # SSR Residuals Vector
    RSS <- c(sum(resid0 ^ 2), sum(resid1 ^ 2), sum(resid2 ^ 2),
             sum(resid3 ^ 2), sum(resid4 ^ 2), sum(resid5 ^ 2),
             sum(resid6 ^ 2), sum(resid7 ^ 2))
    N <- length(resid0)
    # Calculating the chi-squared test statistic
    X <- round(c(-N * log(RSS[1] / RSS[2]), -N * log(RSS[1] / RSS[3]),
                 -N * log(RSS[1] / RSS[4]), -N * log(RSS[1] / RSS[5]),
                 -N * log(RSS[1] / RSS[6]), -N * log(RSS[1] / RSS[7]),
                 -N * log(RSS[1] / RSS[8])))
    # Determining the df for the model in question #
    df <- c(length(coef(Ho)) - length(coef(H1)),
            length(coef(Ho)) - length(coef(H2)),
            length(coef(Ho)) -length(coef(H3)),
            length(coef(Ho)) - length(coef(H4)),
            length(coef(Ho)) - length(coef(H5)),
            length(coef(Ho)) - length(coef(H6)),
            length(coef(Ho)) - length(coef(H7)))
    # Determing the p-value for the model in question #
    p <- round(1 - pchisq(X, df), 4)
    # Creating labels and data frames that summarizes the resulting data #
    labs <- c("Ho", "H1", "H2", "H3", "H4", "H5", "H6", "H7")
    hyp <- c("Linf1=Linf2, k1=k2, t01=t02", "Linf1=Linf2", "k1=k2", "t01=t02",
             "Linf1=Linf2, k1=k2","Linf1=Linf2, t01=t02", "k1=k2, t01=t02")
    hyp2 <- c("All Different", "Linf1=Linf2, k1=k2, t01=t02", "Linf1=Linf2",
              "k1=k2", "t01=t02", "Linf1=Linf2, k1=k2", "Linf1=Linf2, t01=t02",
              "k1=k2, t01=t02")
    labels <- c("Ho vs H1", "Ho vs H2", "Ho vs H3", "Ho vs H4", "Ho vs H5",
                "Ho vs H6", "Ho vs H7")
    aic <- c(AIC(Ho), AIC(H1), AIC(H2), AIC(H3), AIC(H4), AIC(H5), AIC(H6),
             AIC(H7))
    bic <- c(BIC(Ho), BIC(H1), BIC(H2), BIC(H3), BIC(H4), BIC(H5), BIC(H6),
             BIC(H7))
    compout <- data.frame(tests = labels, hypothesis = hyp, chisq = round(X,2),
                          df = df, p = p)
    AIC_BIC <- data.frame(tests = labs, hypothesis = hyp2, AIC = aic, BIC = bic)
    rss <- as.data.frame(cbind(labs, RSS))
    names(rss) <- c("model", "RSS")
    residuals_all <- as.data.frame(cbind(resid0, resid1, resid2, resid3,
                                         resid4, resid5, resid6, resid7))
    nlsout <- list(compout, summary(Ho), summary(H1), summary(H2), summary(H3),
                   summary(H4), summary(H5), summary(H6), summary(H7), rss,
                   residuals_all, AIC_BIC)
    names(nlsout) <- c("Results", c(paste("Model", labs)), "RSS","Residuals","AIC_BIC")
    return(nlsout)
}

