
#' Title
#'
#' @param x
#' @param ntop
#' @param choices
#' @param arrowColors
#' @param biplot
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
genepca <- function(x,ntop,choices=c(1,2),arrowColors = "steelblue", groupNames="group", biplot=FALSE,...) {
  # intgroup <- c("condition","tissue")

  rv <- rowVars(assay(x))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
  pca <- prcomp((assay(x)[select, ]))


  #ggbiplot(pca,choices = c(1,2))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  # if (!all(intgroup %in% names(colData(x)))) {
  # stop("the argument 'intgroup' should specify columns of colData(dds)")
  # }
  # intgroup.df <- as.data.frame(colData(x)[, intgroup, drop = FALSE])
  # group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
  if(biplot){
    ggbiplotFede(pca,arrowColors = arrowColors,choices=choices,groupNames=groupNames,...)
  } else {
    nobs.factor <- sqrt(nrow(pca$x) - 1)
    devs <- pca$sdev
    pcast <- pca
    pcast$x <- sweep(pca$x, 2, 1/(devs * nobs.factor), FUN = "*") * nobs.factor
    d <- data.frame(PC1 = pcast$x[, choices[1]], PC2 = pcast$x[, choices[2]],
                    # group = group,
                    # intgroup.df,
                    names = rownames((assay(x)[select, ])))

    # if (returnData) {
    #   attr(d, "percentVar") <- percentVar[1:2]
    #   return(d)
    # }

    ggplot(data = d, aes_string(x = "PC1", y = "PC2")) +
      geom_point(size = 3) +
      xlab(paste0("PC",choices[1],": ", round(percentVar[choices[1]] * 100), "% variance")) +
      ylab(paste0("PC",choices[2],": ", round(percentVar[choices[2]] * 100), "% variance")) +
      # geom_text(aes(label=names),hjust=0.25, vjust=-0.5, show.legend = F) +
      ggtitle("title") + theme_bw()
  }
}



#' Title
#'
#' @param pcobj
#' @param choices
#' @param scale
#' @param pc.biplot
#' @param obs.scale
#' @param var.scale
#' @param groups
#' @param ellipse
#' @param ellipse.prob
#' @param labels
#' @param labels.size
#' @param alpha
#' @param var.axes
#' @param circle
#' @param circle.prob
#' @param varname.size
#' @param varname.adjust
#' @param varname.abbrev
#' @param arrowColors
#' @param returnData
#' @param coordEqual
#' @param scaleArrow
#' @param useRownamesAsLabels
#' @param point_size
#' @param annotation
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
ggbiplotFede <- function (pcobj, choices = NULL, scale = 1, pc.biplot = TRUE,
                          obs.scale = 1 - scale, var.scale = scale, groups = NULL,
                          ellipse = FALSE, ellipse.prob = 0.68, labels = NULL, labels.size = 3,
                          alpha = 1, var.axes = TRUE, circle = FALSE, circle.prob = 0.69,
                          varname.size = 4, varname.adjust = 1.5, varname.abbrev = FALSE,
                          arrowColors = NULL, groupNames="group",
                          returnData=FALSE,coordEqual=FALSE, scaleArrow = 1,
                          useRownamesAsLabels=TRUE, point_size=2,annotation = NULL,
                          ...)
{
#   library("ggplot2")
#   library("plyr")
#   library("scales")
#   library("grid")

  stopifnot(length(choices) == 2)
  if (inherits(pcobj, "prcomp")) {
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$rotation
  }
  #   else if (inherits(pcobj, "princomp")) {
  #     nobs.factor <- sqrt(pcobj$n.obs)
  #     d <- pcobj$sdev
  #     u <- sweep(pcobj$scores, 2, 1/(d * nobs.factor), FUN = "*")
  #     v <- pcobj$loadings
  #   }
  #   else if (inherits(pcobj, "PCA")) {
  #     nobs.factor <- sqrt(nrow(pcobj$call$X))
  #     d <- unlist(sqrt(pcobj$eig)[1])
  #     u <- sweep(pcobj$ind$coord, 2, 1/(d * nobs.factor), FUN = "*")
  #     v <- sweep(pcobj$var$coord, 2, sqrt(pcobj$eig[1:ncol(pcobj$var$coord),
  #                                                   1]), FUN = "/")
  #   }
  #   else if (inherits(pcobj, "lda")) {
  #     nobs.factor <- sqrt(pcobj$N)
  #     d <- pcobj$svd
  #     u <- predict(pcobj)$x/nobs.factor
  #     v <- pcobj$scaling
  #     d.total <- sum(d^2)
  #   }
  #   else {
  #     stop("Expected a object of class prcomp, princomp, PCA, or lda")
  #   }
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale,
                              FUN = "*"))
  v <- sweep(v, 2, d^var.scale, FUN = "*")
  df.v <- as.data.frame(v[, choices])
  names(df.u) <- c("xvar", "yvar")
  names(df.v) <- names(df.u)
  if (pc.biplot) {
    df.u <- df.u * nobs.factor
  }


  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  v.scale <- rowSums(v^2)
  df.v <- r * df.v/sqrt(max(v.scale))
  if (obs.scale == 0) {
    u.axis.labs <- paste("standardized PC", choices, sep = "")
  } else {
    u.axis.labs <- paste("PC", choices, sep = "")
  }
  u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%% explained var.)",
                                            100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
  if (!is.null(labels)) {
    df.u$labels <- labels
  }
  if (!is.null(groups)) {
    df.u$groups <- groups
  }

  # additionally...
  df.u$ids <- rownames(df.u)
  if(!is.null(annotation)) {
    df.u$geneNames <- annotation$gene_name[match(df.u$ids,rownames(annotation))]
  } else {
    df.u$geneNames <- df.u$ids
  }
  if (varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  } else {
    df.v$varname <- rownames(v)
  }
  df.v$angle <- with(df.v, (180/pi) * atan(yvar/xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar))/2)

  if(returnData){
    return(df.u)
  }

  g <- ggplot(data = df.u, aes_string(x = "xvar", y = "yvar")) + xlab(u.axis.labs[1]) +
    ylab(u.axis.labs[2]) # + coord_equal() # REMOVED OTHERWISE BRUSH DOES NOT WORK PROPERLY
  if(coordEqual) g <- g + coord_equal()

  if (!is.null(df.u$labels)) {
    if (!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups),
                         size = labels.size)
    } else {
      g <- g + geom_text(aes(label = labels), size = labels.size)
    }
  } else {
    if (!is.null(df.u$groups)) {
      g <- g + geom_point(aes(color = groups), size= point_size,alpha = alpha)
    } else {
      g <- g + geom_point(size=point_size,alpha = alpha)
    }
  }

  if(useRownamesAsLabels) {

    g <- g + geom_text(aes_string(label = "geneNames"), size = labels.size,hjust=0.25, vjust=-0.75)
  }

  if (!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    ell <- ddply(df.u, "groups", function(x) {
      if (nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2,
                       mu, FUN = "+"), groups = x$groups[1])
    })
    names(ell)[1:2] <- c("xvar", "yvar")
    g <- g + geom_path(data = ell, aes(color = groups, group = groups))
  }
  # moved down to have the arrows drawn on top of the points and not vice versa
  if (var.axes) {
    if (circle) {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi,
                                                length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r *
                             sin(theta))
      g <- g + geom_path(data = circle, color = "steelblue",
                         size = 1/2, alpha = 1/3)
    }
    df.v$scaleArrow <- scaleArrow # quick fix for mapping scaling of the arrows
    arrowColors <-  as.factor(arrowColors)
    df.v$arrowColors <- arrowColors
    df.v$groupNames <- groupNames
    df.v$sca_x <- df.v$xvar * scaleArrow
    df.v$sca_y <- df.v$yvar * scaleArrow
    df.v$sta_x <- 0
    df.v$sta_y <- 0
    g <- g + geom_segment(data = df.v, aes_string(x = "sta_x", y = "sta_y", xend = "sca_x", yend ="sca_y", color = "arrowColors"),
                          arrow = arrow(length = unit(1/2, "picas"))) +
      scale_color_manual(values = levels(arrowColors),name="Group",labels=levels(groupNames))

  }

  if (var.axes) {
    g <- g + geom_text(data = df.v, aes_string(label = "varname",
                                        x = "sca_x", y = "sca_y",# angle = angle,
                                        hjust = "hjust"),
                       color = arrowColors, size = varname.size)
  }
  g <- g + theme_bw()
  return(g)
}

