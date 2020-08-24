##' @title Plot method for bmgarch objects.
##' @param x bmgarch object.
##' @param type String (Default: "mean"). Whether to plot conditional means ("mean"), variance ("var"), or correlations ("cor"). 
##' @param askNewPage askNewPage Logical (Default: True). Whether to ask for new plotting page.
##' @param CrI CrI Numeric vector (Default: \code{c(.025, .975)}). Lower and upper bound of predictive credible interval.
##' @param ... Not used
##' @return List of ggplot objects (one per time series).
##' @author Stephen R. Martin
##' @import ggplot2
##' @importFrom graphics plot
##' @importFrom grDevices devAskNewPage
##' @export
plot.bmgarch <- function(x, type = "mean", askNewPage = TRUE, CrI = c(.025, .975), ...) {
    ## Bind period, L, U locally to plot.bmgarch to avoid CMD R check  note
    period <- L <- U <- NULL
    x.fitted <- fitted(x, CrI = CrI, digits = 4)
    x.observed <- x$RTS_full
    nt <- x$nt
    TS_names <- lapply(x.fitted$backcast, function(x) {dimnames(x)[[3]]})
    TS_length <- x$TS_length

    df <- as.data.frame(x.fitted)

    # rename % cols to L and U.
    LU <- paste0(CrI * 100, "%")
    colnames(df)[colnames(df) %in% LU] <- c("L", "U")

    plt <- list()

    # Subset by correct param type.
    if(!(type %in% c("mean","var","cor"))) {
        stop("'type' must be 'mean', 'var', or 'cor'.")
    }
    if(type == "cor" & x$param == "CCC") {
        stop("CCC does not model correlations over time.")
    }

    df <- df[df$param == type, ]

    for(i in TS_names[[type]]) {
        df.i <- df[df$TS == i,]

        plt[[i]] <- ggplot(data = df.i, aes(x = period, y = mean, ymin = L, ymax = U)) +
            geom_line() +
            geom_ribbon(alpha = .3)
        l <- labs(x = "Time Period",
                  y = switch(type,
                            mean = "Conditional Means",
                            var = "Conditional Variances",
                            cor = "Conditional Correlations",
                            NULL),
                  title = paste0(x$param, "-MGARCH"),
                  subtitle = i)
        plt[[i]] <- plt[[i]] + l
        print(plt[[i]])
        if(which(i %in% TS_names[[type]]) == 1) {
            devAskNewPage(ask = askNewPage)
        }
    }

    return(invisible(plt))
}
##' @title Plot method for forecast.bmgarch objects.
##' @param x forecast.bmgarch object. See \code{\link{forecast.bmgarch}}.
##' @param type String (Default: "mean"). Whether to plot conditional means ("mean"), variance ("var"), or correlations ("cor"). 
##' @param askNewPage Logical (Default: True). Whether to ask for new plotting page.
##' @param last_t Integer (Default: 100). Only show \code{last_t} observations in plot.
##' @param ... Not used
##' @return List of ggplot objects (one per time series).
##' @author Stephen R. Martin
##' @import ggplot2
##' @importFrom graphics plot
##' @importFrom grDevices devAskNewPage
##' @export
plot.forecast.bmgarch <- function(x, type = "mean", askNewPage = TRUE, last_t = 100, ...) {
    ## Bind period, L, U locally to plot.bmgarch to avoid CMD R check  note
    period <- L <- U <- NULL
    nt <- x$meta$nt
    TS_names <- lapply(x$forecast, function(x) {dimnames(x)[[3]]})
    x.observed <- as.data.frame(x$meta$RTS_full)
    x.observed$period <- seq_len(x$meta$TS_length)

    df <- as.data.frame(x, backcast = TRUE)

    # rename % cols to L and U.
    LU <- paste0(x$meta$CrI * 100, "%")
    colnames(df)[colnames(df) %in% LU] <- c("L", "U")

    # Subset by correct param type.
    if(!(type %in% c("mean","var","cor"))) {
        stop("'type' must be 'mean', 'var', or 'cor'.")
    }
    condCor <- any(sapply(x$meta_list, function(x) {x$param != "CCC"}))
    ## if(type == "cor" & x$meta$param == "CCC") {
    if(type == "cor" & !condCor) {
        stop("CCC does not model correlations over time.")
    }
    df <- df[df$param == type, ]
    df$type <- ifelse(df$type == "backcast", "Backcast", "Forecast")

    plt <- list()
    TS_length <- max(df$period)

    for(i in TS_names[[type]]) {
        df.i <- df[df$TS == i,]

        plt[[i]] <- ggplot(data = df.i, aes(x = period, y = mean, ymin = L, ymax = U, color = type, fill = type)) +
            geom_line() +
            geom_ribbon(alpha = .3)
        l <- labs(x = "Time Period",
                  y = switch(type,
                            mean = "Conditional Means",
                            var = "Conditional Variances",
                            cor = "Conditional Correlations",
                            NULL),
                  title = paste0(x$meta$param, "-MGARCH"),
                  subtitle = i,
                  color = "Type",
                  fill = "Type")
        plt[[i]] <- plt[[i]] + l
        plt[[i]] <- plt[[i]] + coord_cartesian(xlim = c(TS_length - last_t, TS_length))

        print(plt[[i]])
        if(which(i %in% TS_names[[type]]) == 1) {
            devAskNewPage(ask = askNewPage)
        }
    }

    return(invisible(plt))
}
