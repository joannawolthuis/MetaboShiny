#' @title Generate before and after normalization plot
#' @description Function to generate ggplot or plotly plot for data normalization
#' @param mSet mSet object
#' @param cf Function to get plot colors from
#' @return GGPLOT or PLOTLY object(s)
#' @seealso 
#'  \code{\link[reshape2]{melt}}
#'  \code{\link[ggplot2]{ggplot}},\code{\link[ggplot2]{geom_density}},\code{\link[ggplot2]{aes}},\code{\link[ggplot2]{labs}},\code{\link[ggplot2]{geom_boxplot}},\code{\link[ggplot2]{geom_abline}},\code{\link[ggplot2]{scale_continuous}}
#' @rdname ggplotNormSummary
#' @export 
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot geom_density aes ylab xlab geom_boxplot geom_hline scale_y_continuous
ggplotNormSummary <- function(mSet,
                              cf){
  
  # load in original data (pre-normalization, post-filter)
  orig_data <- as.data.frame(mSet$dataSet$proc)
  # load in normalized data
  norm_data <- as.data.frame(mSet$dataSet$norm)
  
  # isolate which samples and mz values are available in both tables
  candidate.samps <- intersect(rownames(orig_data), rownames(norm_data))
  candidate.mzs <- intersect(colnames(orig_data), colnames(norm_data))
  
  # at random, pick 20 compounds and 20 samples to plot data from
  sampsize = if(nrow(norm_data) > 20) 20 else nrow(norm_data)
  which_cpds <- sample(candidate.mzs, sampsize, replace = FALSE, prob = NULL)
  which_samps <- sample(candidate.samps, sampsize, replace = FALSE, prob = NULL)
  
  # isolate these samples from original table and melt into long format (ggplot needs it!)
  orig_melt <- reshape2::melt(cbind(which_samps, orig_data[which_samps, which_cpds]))
  orig_melt[is.na(orig_melt)] <- 0 # replace NA values with 0
  
  # isolate these samples from normalized table and melt into long format
  norm_melt <- reshape2::melt(cbind(which_samps, norm_data[which_samps, which_cpds]))
  
  # create base plot with base theme and font size for original data
  plot <- ggplot2::ggplot(data=orig_melt)
  
  # first result plot: is a density plot of chosen 20 mz values with 20 samples
  RES1 <- plot + ggplot2::geom_density(ggplot2::aes(x=value,y=..scaled..), colour="blue", fill="blue", alpha=0.4) +
    ggplot2::ylab("density") + 
    ggplot2::xlab("intensity")
  
  # second result plot: shows the spread of the intensities before normalization
  RES2 <- plot + ggplot2::geom_boxplot(
    ggplot2::aes(y=value, x=variable),
    color=cf(sampsize),
    alpha=0.4) + ggplot2::geom_hline(ggplot2::aes(yintercept=median(value))) + 
    ggplot2::coord_flip() +
    ggplot2::xlab("m/z") + 
    ggplot2::ylab("intensity")
  
  # create base plot with base theme and font size for normalized data
  plot <- ggplot2::ggplot(data=norm_melt)
  
  # third result plot: a density plot of chosen 20 mz values post normalization
  RES3 <- plot + ggplot2::geom_density(ggplot2::aes(x=value,y=..scaled..), colour="pink", fill="pink", alpha=0.4) + 
    ggplot2::ylab("density") + 
    ggplot2::xlab("intensity")
  
  # fourth result plot: spread of intensities after normalization
  RES4 <- plot + ggplot2::geom_boxplot(
    ggplot2::aes(y=value, x=variable),
    color=cf(sampsize),
    alpha=0.4) + ggplot2::geom_hline(ggplot2::aes(yintercept=median(value))) + 
    ggplot2::coord_flip() +
    ggplot2::xlab("m/z") + 
    ggplot2::ylab("intensity")
  
  scaleFUN <- function(x) sprintf("%.1f", x)
  
  list(tl=RES1 + ggplot2::scale_y_continuous(labels=scaleFUN), 
       bl=RES2, 
       tr=RES3 + ggplot2::scale_y_continuous(labels=scaleFUN), 
       br=RES4)
  
}

#' @title Generate sample normalization before/after plots
#' @description Function to generate ggplot or plotly plot for sample intensity normalization
#' @param mSet mSet object
#' @param cf Function to get plot colors from
#' @return GGPLOT or PLOTLY object(s)
#' @seealso 
#'  \code{\link[reshape2]{melt}}
#'  \code{\link[ggplot2]{ggplot}},\code{\link[ggplot2]{geom_density}},\code{\link[ggplot2]{aes}},\code{\link[ggplot2]{labs}},\code{\link[ggplot2]{geom_boxplot}},\code{\link[ggplot2]{geom_abline}},\code{\link[ggplot2]{scale_continuous}}
#' @rdname ggplotSampleNormSummary
#' @export 
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot geom_density aes ylab xlab geom_boxplot geom_hline scale_y_continuous
ggplotSampleNormSummary <- function(mSet,
                                    cf){
  # 4 by 4 plot, based on random 20-30 picked
  orig_data <- as.data.frame(mSet$dataSet$proc)
  norm_data <- as.data.frame(mSet$dataSet$norm)
  
  candidate.samps <- intersect(rownames(orig_data), rownames(norm_data))
  candidate.mzs <- intersect(colnames(orig_data), colnames(norm_data))
  
  sampsize = if(nrow(norm_data) > 20) 20 else nrow(norm_data)
  
  which_samps <- sample(candidate.samps,
                        sampsize,
                        replace = FALSE,
                        prob = NULL)
  
  sumsOrig <- rowSums(orig_data[which_samps,])
  sumsNorm <- rowSums(norm_data[which_samps,])
  
  orig_data$Label <- rownames(orig_data)
  orig_melt <- reshape2::melt(orig_data[which_samps,],
                              id.vars = "Label")
  orig_melt_sums <- reshape2::melt(sumsOrig)
  orig_melt_sums$variable <- rownames(orig_melt_sums)
  
  norm_data$Label <- rownames(norm_data)
  norm_melt <- reshape2::melt(norm_data[which_samps,],
                              id.vars="Label")
  norm_melt_sums <- reshape2::melt(sumsNorm)
  norm_melt_sums$variable <- rownames(norm_melt_sums)
  
  RES1 <- ggplot2::ggplot(data=orig_melt_sums) +
    ggplot2::geom_density(ggplot2::aes(x=value,y=..scaled..), colour="blue", fill="blue", alpha=0.4) + 
    ggplot2::ylab("density") + 
    ggplot2::xlab("intensity")
  
  RES2 <- ggplot2::ggplot(data=orig_melt) +
    ggplot2::geom_boxplot(
      ggplot2::aes(y=value,x=Label),
      color=cf(sampsize),
      alpha=0.4) + ggplot2::geom_hline(ggplot2::aes(yintercept=median(value),text=Label)) + ggplot2::coord_flip() +
    ggplot2::xlab("m/z") + 
    ggplot2::ylab("intensity")
  
  RES3 <- ggplot2::ggplot(data=norm_melt_sums) +
    ggplot2::geom_density(ggplot2::aes(x=value, y=..scaled..), colour="pink", fill="pink", alpha=0.4) +
    ggplot2::ylab("density") + 
    ggplot2::xlab("intensity")
  
  RES4 <- ggplot2::ggplot(data=norm_melt) +
    ggplot2::geom_boxplot(
      ggplot2::aes(y=value,x=Label),
      color=cf(sampsize),
      alpha=0.4) + ggplot2::geom_hline(ggplot2::aes(yintercept=median(value),text=Label))+ggplot2::coord_flip() + 
    ggplot2::xlab("m/z") + 
    ggplot2::ylab("intensity")
  
  scaleFUN <- function(x) sprintf("%.2f", x)
  
  list(tl=RES1 + ggplot2::scale_y_continuous(labels=scaleFUN), 
       bl=RES2, 
       tr=RES3 + ggplot2::scale_y_continuous(labels=scaleFUN), 
       br=RES4)
  
}


#' @title Generate MEBA plot
#' @description Function to generate ggplot or plotly plot for MEBA analysis
#' @param mSet mSet object
#' @param cpd m/z value of interest
#' @param draw.average PARAM_DESCRIPTION, Default: T
#' @param cols Colors to use
#' @param cf Function to get plot colors from
#' @return GGPLOT or PLOTLY object(s)
#' @seealso 
#'  \code{\link[stringr]{str_match}}
#'  \code{\link[ggplot2]{ggplot}},\code{\link[ggplot2]{geom_path}},\code{\link[ggplot2]{aes}},\code{\link[ggplot2]{scale_x_discrete}},\code{\link[ggplot2]{scale_manual}},\code{\link[ggplot2]{stat_summary_bin}}
#'  \code{\link[Hmisc]{capitalize}}
#' @rdname ggplotMeba
#' @export 
#' @importFrom stringr str_match
#' @importFrom ggplot2 ggplot geom_line aes scale_x_discrete scale_color_manual stat_summary
#' @importFrom Hmisc capitalize
ggplotMebaSingle <- function(mSet, cpd, draw.average=T, cols,
                       cf){
  
  time.mode = mSet$settings$exp.type
  classes = unique(switch(time.mode, 
                          t1f=mSet$dataSet$facA, 
                          t=mSet$dataSet$sbj))
  spec.cols = cf(length(classes))
  cols <- if(is.null(cols)) spec.cols else{
    if(length(cols) < length(classes)){
      cols <- spec.cols
    }
    cols
  }
  
  profile <- getProfile(mSet, 
                        cpd, 
                        mode="multi")
  
  cpd = stringr::str_match(cpd, "(\\d+\\.\\d+)")[,2]
  
  profile$Individual <- mSet$dataSet$covars[match(profile$Sample,
                                                  table = mSet$dataSet$covars$sample),"individual"][[1]]
  profile$Color <- switch(time.mode, 
                          t1f=profile$GroupA, 
                          t=profile$Individual, 
                          group=profile$GroupA)
  
  p <- ggplot2::ggplot(data=profile) +
    ggplot2::geom_line(size=if(draw.average) 0.3 else 1, ggplot2::aes(x=GroupB, 
                                                                      y=Abundance, 
                                                                      group=Individual, 
                                                                      color=Color,
                                                                      text=Sample), alpha=0.4) +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    
    ggplot2::scale_color_manual(values=cols) +
    ggtitle(paste(cpd, "m/z")) + 
    xlab(Hmisc::capitalize(mSet$dataSet$facB.lbl)) + 
    ggplot2::labs(color = Hmisc::capitalize(switch(mSet$settings$exp.type, 
                                                   t1f=mSet$dataSet$facA.lbl,
                                                   t="Individual")))
  if(draw.average){
    p <- p + ggplot2::stat_summary(fun="mean", size=2, 
                                   geom="line", ggplot2::aes(x=GroupB, 
                                                             y=Abundance, 
                                                             color = Color, 
                                                             group = switch(time.mode, 
                                                                            t=c(1),
                                                                            t1f=Color)))
  }
  p
}

#' @title Generate black/white gradient
#' @description Function to generate black/white gradient 
#' @param n Number of colors to include in gradient
#' @return Vector of color values
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  bw.cols = blackwhite.colors(256)
#'  }
#' }
#' @rdname blackwhite.colors
#' @export 
blackwhite.colors <- function(n){
  gray.colors(n, start=0, end=1)
}

#' @title Generate intensity summary plot
#' @description Function to generate ggplot or plotly plot for the current m/z selected
#' @param mSet mSet object
#' @param cpd m/z value of interest
#' @param shape.fac Change shape based on this metadata column, Default: 'label'
#' @param cols Colors to use, Default: c("black", "pink")
#' @param cf Function to get plot colors from, Default: rainbow
#' @param mode Normal(nm), time series etc., Default: 'nm'
#' @param styles Which plot styles to apply (each adds a new layer), Default: c("box", "beeswarm")
#' @param add_stats Add statistics-based line in plot? And use what to do so?, Default: 'mean'
#' @param color.fac Change fill color based on this metadata column, Default: 'label'
#' @param text.fac Change hover text based on this metadata column, Default: 'label'
#' @return GGPLOT or PLOTLY object(s)
#' @seealso 
#'  \code{\link[stringr]{str_match}}
#'  \code{\link[ggplot2]{ggplot}},\code{\link[ggplot2]{geom_boxplot}},\code{\link[ggplot2]{geom_violin}},\code{\link[ggplot2]{geom_point}},\code{\link[ggplot2]{annotate}},\code{\link[ggplot2]{stat_summary_bin}},\code{\link[ggplot2]{scale_manual}},\code{\link[ggplot2]{labs}}
#'  \code{\link[Hmisc]{capitalize}}
#'  \code{\link[ggbeeswarm]{geom_beeswarm}}
#' @rdname ggplotSummary
#' @export 
#' @importFrom stringr str_match
#' @importFrom data.table data.table
#' @importFrom ggplot2 ggplot geom_boxplot geom_violin geom_point annotate stat_summary scale_color_manual scale_fill_manual xlab
#' @importFrom Hmisc capitalize
#' @importFrom ggbeeswarm geom_beeswarm
ggplotSummary <- function(mSet, cpd, shape.fac = "label", cols = c("black", "pink"),
                          cf = rainbow, 
                          mode = "nm", 
                          styles=c("box", "beeswarm"), add_stats = "mean",
                          color.fac = "label",
                          text.fac = "label",
                          fill.fac="label"){
  
  sourceTable = mSet$dataSet$norm
  
  if(length(styles) == 0){
    styles = c("beeswarm")
  }
  # - - -
  
  if(mSet$settings$exp.type %in% c("t","t1f", "2f")){
    mode = "multi"
  }
  
  profile <- getProfile(mSet, 
                        cpd, 
                        mode=if(mode == "nm") "stat" else "multi")
  
  cpd = stringr::str_match(cpd, "(\\d+\\.\\d+)")[,2]
  
  df_line <- data.table::data.table(x = c(1,2),
                                    y = rep(min(profile$Abundance - 0.1), 2))
  stars = ""
  
  # try({
  #   pval <- if(mode == "nm"){
  #     mSet$analSet$tt$sig.mat[which(rownames(mSet$analSet$tt$sig.mat) == cpd), "p.value"]
  #   }else{
  #     int.col <- grep("adj|Adj", colnames(mSet$analSet$aov2$sig.mat),value=T)
  #     int.col <- grep("int|Int", int.col, value=T)
  #     mSet$analSet$aov2$sig.mat[which(rownames(mSet$analSet$aov2$sig.mat) == cpd), int.col]
  #   }
  #   stars <- p2stars(pval)
  # })
  
  p <- ggplot2::ggplot()
  
  for(adj in c("color", "shape", "text", "fill")){
    adj.fac = switch(adj,
                     "shape" = shape.fac,
                     "color" = color.fac,
                     "text" = text.fac,
                     "fill" = fill.fac)
    profile[[Hmisc::capitalize(adj)]] <- if(adj.fac != "label") as.factor(mSet$dataSet$covars[[adj.fac]]) else {
      if(adj == "fill"){
        switch(mode, 
               multi = profile$GroupB, 
               nm = profile$Group)
      }else{
        switch(mode, 
               multi = profile$GroupA, 
               nm = profile$Group)  
      }}
  }
  
  nshape = length(unique(profile$Shape))
  
  if(nshape > 5){
    symbols = c(1:25)
    # fill > color
    print("Too many shapes - fill property only works with <6 shapes. Outlines removed.")
    profile$Color <- profile$Fill
  }else{
    symbols = c(21:25) 
  }
  
  profiles <- switch(mode,
                     multi = split(profile, f = profile$GroupA),
                     nm = list(x = data.table::data.table(profile)))
  
  if(length(cols) < length(levels(profile$Color))){
    cols <- cf(length(levels(profile$Color)))
  }
  
  i = 1
  
  suppressWarnings({
    
    for(prof in profiles){
      
      groupCols <- grep(x=colnames(prof), pattern="Group", value=T)
      
      for(groupCol in groupCols){
        prof[, (groupCol) := factor(get(groupCol), levels = {
          lvls = levels(get(groupCol))
          numconv = as.numeric(as.character(lvls))
          if(all(!is.na(numconv))){
            order = order(numconv)
          }else{
            order = order(as.character(lvls))
          }
          lvls[order]
        })] 
      }
      
      for(style in styles){
        switch(mode,
               nm = {
                 p <- switch(style,
                             box = p + ggplot2::geom_boxplot(data = prof, alpha=0.4, ggplot2::aes(x = Group,
                                                                                                  y = Abundance,
                                                                                                  shape = Shape,
                                                                                                  text = Text,
                                                                                                  color = Color,
                                                                                                  fill = Fill)),
                             violin = p + ggplot2::geom_violin(data = prof, alpha=0.4, position = "identity", ggplot2::aes(x = Group,
                                                                                                                           y = Abundance,
                                                                                                                           color = Color,
                                                                                                                           fill = Fill)),
                             beeswarm = p + ggbeeswarm::geom_beeswarm(data = prof, alpha=0.7, size = 2, 
                                                                      #position = position_dodge(width=.3), 
                                                                      ggplot2::aes(x = Group,
                                                                                   y = Abundance,
                                                                                   text = Text,
                                                                                   shape = Shape,
                                                                                   color = Color,
                                                                                   fill = Fill)),
                             scatter = p + ggplot2::geom_point(data = prof, alpha=0.7, size = 2, ggplot2::aes(x = Group,
                                                                                                              y = Abundance,
                                                                                                              text=Text,
                                                                                                              shape = Shape,
                                                                                                              color = Color,
                                                                                                              fill = Fill), 
                                                               position = position_jitterdodge())
                 )
               },
               multi = {
                 p <- switch(style,
                             box = p + ggplot2::geom_boxplot(data = prof, alpha=0.4, ggplot2::aes(x = GroupB,
                                                                                                  y = Abundance,
                                                                                                  text = Text,
                                                                                                  shape = Shape,
                                                                                                  color = Color,
                                                                                                  fill = Fill)),
                             violin = p + ggplot2::geom_violin(data = prof, alpha=0.4, position = "identity", ggplot2::aes(x = GroupB,
                                                                                                                           y = Abundance,
                                                                                                                           group = GroupB,
                                                                                                                           text = Text,
                                                                                                                           color = Color,
                                                                                                                           fill = Fill)),#GroupB)),
                             beeswarm = {
                               p + ggbeeswarm::geom_beeswarm(data = prof, alpha=0.7, size = 2, position = ggplot2::position_dodge(width=.3), ggplot2::aes(x = GroupB,
                                                                                                                                                          y = Abundance,
                                                                                                                                                          text=Text,
                                                                                                                                                          shape = Shape,
                                                                                                                                                          color = Color,
                                                                                                                                                          fill = Fill))},#GroupA))},
                             scatter = p + ggplot2::geom_point(data = prof, alpha=0.7, size = 2, position = ggplot2::position_jitterdodge(), ggplot2::aes(x = GroupB,
                                                                                                                                                          y = Abundance,
                                                                                                                                                          text=Text,
                                                                                                                                                          shape = Shape,
                                                                                                                                                          color = Color,
                                                                                                                                                          fill = Fill))#GroupA))
                 )
               })
      }
      
      p <- p + ggplot2::annotate("text",
                                 x = switch(mode, nm = 1.5,
                                            multi = max(as.numeric(as.factor(profile$GroupB)))/2 + .5),
                                 y = min(profile$Abundance - 0.3),
                                 label = stars, size = 8, col = "black") + ggplot2::ggtitle(paste(cpd, "m/z"))
      
      if(!("box" %in% styles)){
        p <- switch(add_stats,
                    median = {
                      p + ggplot2::stat_summary(data = prof,
                                                ggplot2::aes( x = if(mode == "nm") Group else GroupB,
                                                              y = Abundance,
                                                              color = if(mode == "nm") Group else GroupA),
                                                fun = median,
                                                fun.ymin = median,
                                                fun.ymax = median,
                                                geom = "crossbar",
                                                width = 0.5,
                                                color = switch(mode,
                                                               multi = cols[i],
                                                               nm = cols[1:length(unique(levels(profile$Group)))]))
                    },
                    mean = {
                      p + ggplot2::stat_summary(data = prof,
                                                ggplot2::aes(x = if(mode == "nm") Group else GroupB,
                                                             y = Abundance),
                                                fun = mean,
                                                fun.ymin = mean,
                                                fun.ymax = mean,
                                                geom = "crossbar",
                                                width = 0.5,
                                                color = switch(mode,
                                                               ts = cols[i],
                                                               nm = cols[1:length(unique(levels(profile$Group)))]))
                    },
                    none = {
                      p
                    }
        )
      }
      i <- i + 1
    }
    
    if(mode == "multi"){
      p <- p + ggplot2::scale_color_manual(values=cols)
      p <- p + ggplot2::scale_fill_manual(values=cols)
      if(mSet$settings$exp.type == "t"){
        p <- p + ggplot2::xlab("Time")
      }else{
        p <- p + ggplot2::xlab(Hmisc::capitalize(mSet$dataSet$facB.lbl))
      }
    }else{
      p <- p + ggplot2::scale_color_manual(values = cols) + 
               ggplot2::scale_fill_manual(values = cols)
      p <- p + ggplot2::xlab(Hmisc::capitalize(gsub(x=mSet$settings$cls.name, pattern = ":.*$", replacement="")))
    }
    p <- p + ggplot2::scale_shape_manual(values = symbols) + ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(shape = 21)),
                                                                             color = ggplot2::guide_legend(override.aes = list(shape = 21)))
    p  
  })
}

#' @export 
ggPlotASCA <- function(mSet, cf, n=20){
  
  profile = data.table::as.data.table(mSet$analSet$asca$sig.list$Model.ab, keep.rownames = T)
  
  if(nrow(profile)==0){
    shiny::showNotification("No significant hits")
    return(NULL)
  }
  
  profile$Peak <- c(1:nrow(profile))
  colnames(profile)[1:3] <- c("m/z", "Leverage", "SPE")
  scaleFUN <- function(x) sprintf("%.5f", x)
  scaleFUN2 <- function(x) sprintf("%.0f", x)
  p <- ggplot2::ggplot(data=profile) +
    ggplot2::geom_point(ggplot2::aes(y=SPE,
                                     x=Leverage,
                                     text=`m/z`,
                                     color=`SPE`, 
                                     key=`m/z`)) +
    ggplot2::scale_colour_gradientn(colours = cf(n)) + 
    ggplot2::coord_flip() +
    ggplot2::scale_x_log10(labels = scaleFUN) +
    ggplot2::scale_y_continuous(labels=scaleFUN2)
  p
}

#' @export 
ggPlotMeba <- function(mSet, cf, n=20, topn=NULL){
  profile = data.table::as.data.table(mSet$analSet$MB$stats, keep.rownames = T)
  if(nrow(profile)==0){
    shiny::showNotification("No significant hits")
    return(NULL)
  }
  
  if(!is.null(topn)){
    profile = profile[order(abs(V2), decreasing = T)]
    profile = profile[1:min(topn, nrow(profile)),]
  }
  
  profile$Peak <- c(1:nrow(profile))
  colnames(profile)[1:2] <- c("m/z", "Hotelling-T2")
  scaleFUN <- function(x) sprintf("%.1f", x)
  scaleFUN2 <- function(x) sprintf("%.0f", x)
  xaxis = seq(0,600, 50)
  p <- ggplot2::ggplot(data=profile) +
    ggplot2::geom_point(ggplot2::aes(y=Peak,
                                     x=`Hotelling-T2`,
                                     text=`m/z`,
                                     color=`Hotelling-T2`, 
                                     key=`m/z`)) +
    ggplot2::geom_segment(aes(y = Peak, 
                              yend = Peak,
                              color=`Hotelling-T2`,
                              x = 0,
                              xend = `Hotelling-T2`)) +
    ggplot2::scale_colour_gradientn(colours = cf(n)) + ggplot2::coord_flip() +
    ggplot2::scale_x_log10(labels = scaleFUN)+
    ggplot2::scale_y_continuous(labels=scaleFUN2)
  p
}


#' @title Generate ANOVA plot
#' @description Function to generate ggplot or plotly plot for ANOVA
#' @param mSet mSet object
#' @param cf Function to get plot colors from
#' @param n Amount of colors in gradient, Default: 20
#' @return GGPLOT or PLOTLY object(s)
#' @seealso 
#'  \code{\link[data.table]{as.data.table}}
#'  \code{\link[shiny]{showNotification}}
#'  \code{\link[ggplot2]{ggplot}},\code{\link[ggplot2]{geom_point}},\code{\link[ggplot2]{aes}},\code{\link[ggplot2]{scale_x_discrete}},\code{\link[ggplot2]{scale_colour_gradient}},\code{\link[ggplot2]{scale_continuous}}
#' @rdname ggPlotAOV
#' @export 
#' @importFrom data.table as.data.table
#' @importFrom shiny showNotification
#' @importFrom ggplot2 ggplot geom_point aes scale_x_discrete scale_colour_gradientn scale_y_continuous
ggPlotAOV <- function(mSet, cf, n=20, topn=NULL){
  
  which_aov = if(mSet$settings$exp.type %in% c("t","2f", "t1f")) "aov2" else "aov"
  
  profile <- if(which_aov == "aov"){
    data.table::as.data.table(mSet$analSet[[which_aov]]$p.log[mSet$analSet[[which_aov]]$inx.imp],
                              keep.rownames = T)
  }else{
    data.table::as.data.table(-log(mSet$analSet[[which_aov]]$sig.mat[,if(mSet$settings$exp.type == "t") "Adjusted P-val" else "Interaction(adj.p)"]),
                              keep.rownames = T)  
  }
  
  if(nrow(profile)==0){
    shiny::showNotification("No significant hits")
    return(NULL)
  }
  
  if(!is.null(topn)){
    profile = profile[order(abs(V2), decreasing = if(which_aov == "aov") T else F)]
    profile = profile[1:min(topn, nrow(profile)),]
  }
  
  #profile[,2] <- round(profile[,2], digits = 2)
  profile$Peak <- c(1:nrow(profile))
  colnames(profile)[1:2] <- c("m/z", "-log(p)")
  scaleFUN <- function(x) sprintf("%.2f", x)
  scaleFUN2 <- function(x) sprintf("%.0f", x)
  xaxis = seq(0,600, 50)
  p <- ggplot2::ggplot(data=profile) +
    ggplot2::geom_point(ggplot2::aes(y=Peak,
                                     x=`-log(p)`,
                                     text=`m/z`,
                                     color=`-log(p)`, 
                                     key=`m/z`)) +
    ggplot2::geom_segment(aes(y = Peak, 
                              yend = Peak,
                              color=`-log(p)`,
                              x = 0,
                              xend = `-log(p)`)) +
    ggplot2::scale_colour_gradientn(colours = cf(n)) + ggplot2::coord_flip() +
    ggplot2::scale_x_continuous(labels=scaleFUN) + 
    ggplot2::scale_y_continuous(labels=scaleFUN2) + 
    ggplot2::xlab(if(which_aov == "aov") "-log(p)" else if(mSet$settings$exp.type == "t") "-log(Adjusted P-val)" else "-log(Interaction(adj.p))")
  p
}

#' @title Generate T-TEST plot
#' @description Function to generate ggplot or plotly plot for T-TEST
#' @param mSet mSet object
#' @param cf Function to get plot colors from
#' @param n Number of colors in gradient, Default: 20
#' @return GGPLOT or PLOTLY object(s)
#' @seealso 
#'  \code{\link[data.table]{as.data.table}}
#'  \code{\link[shiny]{showNotification}}
#'  \code{\link[ggplot2]{ggplot}},\code{\link[ggplot2]{geom_point}},\code{\link[ggplot2]{aes}},\code{\link[ggplot2]{scale_colour_gradient}}
#' @rdname ggPlotTT
#' @export 
#' @importFrom data.table as.data.table
#' @importFrom shiny showNotification
#' @importFrom ggplot2 ggplot geom_point aes scale_colour_gradientn
ggPlotTT <- function(mSet, cf, n=20, topn=NULL){
  profile <- data.table::as.data.table(mSet$analSet$tt$p.log[mSet$analSet$tt$inx.imp],keep.rownames = T)
  
  if(nrow(profile)==0){
    shiny::showNotification("No significant hits")
    return(NULL)
  }
  
  if(!is.null(topn)){
    profile = profile[order(abs(V2), decreasing = T)]
    profile = profile[1:min(topn, nrow(profile)),]
  }
  
  profile[,2] <- round(profile[,2], digits = 2)
  profile$Peak <- c(1:nrow(profile))
  colnames(profile)[1:2] <- c("m/z", "-log(p)")
  profile[["-log(p)"]] <- as.numeric(sprintf("%.1f", profile[["-log(p)"]]))
  
  xaxis = seq(0,600, 50)
  
  # ---------------------------
  
  p = ggplot2::ggplot() +
    ggplot2::geom_point(data=profile,ggplot2::aes(y=Peak,
                                                  x=`-log(p)`,
                                                  text=`m/z`,
                                                  color=`-log(p)`, 
                                                  key=`m/z`),
                        size=2.5) +
    ggplot2::geom_segment(data=profile,aes(y = Peak,
                              yend = Peak,
                              color=`-log(p)`,
                              x = 0,
                              xend = `-log(p)`)) +
    ggplot2::scale_colour_gradientn(colours = cf(n)) +
    ggplot2::coord_flip()
  #ggplot2::scale_y_continuous()
  p
}

#' @title Generate pattern analysis plot
#' @description Function to generate ggplot or plotly plot for pattern analysis
#' @param mSet mSet object
#' @param cf Function to get plot colors from
#' @param n Number of colors in gradient, Default: 20
#' @return GGPLOT or PLOTLY object(s)
#' @seealso 
#'  \code{\link[data.table]{as.data.table}}
#'  \code{\link[shiny]{showNotification}}
#'  \code{\link[ggplot2]{ggplot}},\code{\link[ggplot2]{geom_bar}},\code{\link[ggplot2]{labs}},\code{\link[ggplot2]{coord_flip}},\code{\link[ggplot2]{scale_colour_gradient}},\code{\link[ggplot2]{scale_continuous}}
#' @rdname ggPlotPattern
#' @export 
#' @importFrom data.table as.data.table
#' @importFrom shiny showNotification
#' @importFrom ggplot2 ggplot geom_bar ggtitle coord_flip ylab xlab labs scale_colour_gradientn scale_fill_gradientn scale_y_continuous
ggPlotPattern <- function(mSet, cf, n=20){
  profile <- data.table::as.data.table(mSet$analSet$corr$cor.mat,keep.rownames = T)
  profile <- profile[1:n]
  
  if(nrow(profile)==0){
    shiny::showNotification("No significant hits")
    return(NULL)
  }
  
  colnames(profile)[1] <- c("m/z")
  #profile$Peak <- c(1:nrow(profile))
  scaleFUN <- function(x) sprintf("%.2f", x)
  
  profile$`m/z` <- reorder(x = profile$`m/z`, X = -profile$`p-value`)
  
  # ---------------------------
  p <- ggplot2::ggplot(data=profile) +
    ggplot2::geom_bar(mapping = ggplot2::aes(x = `m/z`, 
                                             y = correlation, 
                                             key = `m/z`, 
                                             text = `m/z`, 
                                             color = `p-value`, 
                                             fill = `p-value`), 
                      stat = "identity", alpha=0.5) + 
    ggplot2::ggtitle("Correlated m/z values") +
    ggplot2::coord_flip() +
    ggplot2::ylab("correlation") + 
    ggplot2::xlab("m/z") + 
    ggplot2::labs(fill="p-value", 
                  color="p-value") + 
    ggplot2::scale_colour_gradientn(colours = cf(n)) +
    ggplot2::scale_fill_gradientn(colours = cf(n)) +
    ggplot2::scale_y_continuous(labels=scaleFUN)
  
  p
}

#' @title Generate fold-change analysis plot
#' @description Function to generate ggplot or plotly plot for fold-change analysis
#' @param mSet mSet object
#' @param cf Function to get plot colors from
#' @param n Number of colors in gradient, Default: 20
#' @return GGPLOT or PLOTLY object(s)
#' @seealso 
#'  \code{\link[data.table]{as.data.table}}
#'  \code{\link[shiny]{showNotification}}
#'  \code{\link[ggplot2]{ggplot}},\code{\link[ggplot2]{geom_point}},\code{\link[ggplot2]{aes}},\code{\link[ggplot2]{geom_abline}},\code{\link[ggplot2]{scale_continuous}},\code{\link[ggplot2]{scale_colour_gradient}}
#' @rdname ggPlotFC
#' @export 
#' @importFrom data.table as.data.table
#' @importFrom shiny showNotification
#' @importFrom ggplot2 ggplot geom_point aes geom_vline scale_y_continuous scale_colour_gradientn
ggPlotFC <- function(mSet, cf, n=20, topn=NULL){
  profile <- data.table::as.data.table(mSet$analSet$fc$fc.log[mSet$analSet$fc$inx.imp],keep.rownames = T)
  
  if(nrow(profile)==0){
    shiny::showNotification("No significant hits")
    return(NULL)
  }
  
  if(!is.null(topn)){
    profile = profile[order(abs(V2), decreasing = T)]
    profile = profile[1:min(topn, nrow(profile)),]
  }
  
  colnames(profile) <- c("m/z", "log2fc")
  profile$Peak <- c(1:nrow(profile))
  scaleFUN <- function(x) sprintf("%.0f", x)
  # ---------------------------
  p <- ggplot2::ggplot(data=profile) +
    ggplot2::geom_point(ggplot2::aes(y=Peak, 
                                     x=log2fc, 
                                     color=log2fc, 
                                     key=`m/z`,
                                     text=`m/z`)) +
    ggplot2::geom_segment(aes(y = Peak, 
                              yend = Peak,
                              color=log2fc,
                              x = 0,
                              xend = log2fc)) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = 0)) +
    ggplot2::scale_y_continuous(labels=scaleFUN) +
    ggplot2::scale_colour_gradientn(colours = cf(n)) + 
    ggplot2::coord_flip()
  
  p
}

#' @importFrom ggplot2 ggplot geom_point aes scale_colour_gradientn
ggPlotCombi <- function(mSet,
                        cf,
                        n=20){
  
    dt = data.table::as.data.table(mSet$analSet$combi$sig.mat)
    anal1_trans = mSet$analSet$combi$trans$x
    anal2_trans = mSet$analSet$combi$trans$y
    
    anal1 = mSet$analSet$combi$source$x
    anal2 = mSet$analSet$combi$source$y
    
    anal1_col = colnames(dt)[2]
    anal2_col = colnames(dt)[3]
    
    colnames(dt) <- c("m/z", "x", "y")#anal1_col, anal2_col)
    dt$col <- abs(dt[,2]*dt[,3])
    dt$significant <- "YES"
    
    if(length(mSet$analSet$combi$all.vals$x) > 0 & length(mSet$analSet$combi$all.vals$y) > 0){
      x.all = mSet$analSet$combi$all.vals$x
      x.tbl = data.table::data.table("m/z" = names(x.all),
                                     x = x.all)
      y.all = mSet$analSet$combi$all.vals$y
      y.tbl = data.table::data.table("m/z" = names(y.all),
                                     x = y.all)
      dt.merged = merge(x.tbl, y.tbl, by.x="m/z", by.y="m/z")
      
      dt.all = data.table::data.table("m/z" = dt.merged$`m/z`,
                                         x = dt.merged$x.x,
                                         y = dt.merged$x.y,
                                         col = c(0),
                                         significant = "NO")
      dt.all[`m/z` %in% dt$`m/z`]$significant <- "YES"
      dt = dt.all
    }
    
    scaleFUN <- function(x) sprintf("%.2f", x)
    
    anal1_col = if(anal1_trans != "none") paste0(anal1_trans,"(", anal1_col,")") else anal1_col
    anal2_col = if(anal2_trans != "none") paste0(anal2_trans,"(", anal2_col,")") else anal2_col
    
    p <- ggplot2::ggplot() +
      ggplot2::geom_point(data=dt, ggplot2::aes(x=x,
                                                y=y,
                                                color=significant,#col,
                                                text=`m/z`,
                                                key=`m/z`),cex=3) +
                          #,alpha = sapply(dt$sig, function(x) ifelse(x=="YES",1,0.2))) +
      ggplot2::geom_segment(data=dt[significant=="YES"],aes(y = 0,
                                        yend = y,
                                        x = 0,
                                        xend = x,
                                        color = significant#col 
      ),alpha=0.2,linetype=6) +
      scale_color_manual(values=c("YES" = "red",
                                  "NO" = "darkgray")) +
      #ggplot2::scale_colour_gradientn(colours = cf(n),guide=FALSE) +
      ggplot2::scale_x_continuous(labels=scaleFUN) + 
      ggplot2::xlab(paste0(anal1, ": ", anal1_col)) + 
      ggplot2::ylab(paste0(anal2, ": ", anal2_col))
    p
}

#' @title Generate volcano plot
#' @description Function to generate ggplot or plotly plot for volcano plot
#' @param mSet mSet object
#' @param cf Function to get plot colors from
#' @param n Number of colors in gradient, Default: 20
#' @return GGPLOT or PLOTLY object(s)
#' @seealso 
#'  \code{\link[shiny]{showNotification}}
#'  \code{\link[ggplot2]{ggplot}},\code{\link[ggplot2]{geom_point}},\code{\link[ggplot2]{aes}},\code{\link[ggplot2]{scale_colour_gradient}}
#' @rdname ggPlotVolc
#' @export 
#' @importFrom shiny showNotification
#' @importFrom ggplot2 ggplot geom_point aes scale_colour_gradientn
ggPlotVolc <- function(mSet,
                       cf,
                       n=20){

  vcn<-mSet$analSet$volcano;
  
  if(nrow(vcn$sig.mat)==0){
    shiny::showNotification("No significant hits")
    return(NULL)
  }
  
  dt <- as.data.frame(vcn$sig.mat)[,c(2,4)]
  dt <- cbind(cpd = rownames(dt), dt)
  colnames(dt) <- c("m/z", "log2FC", "-log10P")
  dt$col <- with(dt, abs(log2FC*`-log10P`))
  
  scaleFUN <- function(x) sprintf("%.2f", x)
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data=dt, ggplot2::aes(x=log2FC,
                                              y=`-log10P`,
                                              color=col,
                                              text=`m/z`,
                                              key=`m/z`),cex=3) +
    ggplot2::geom_segment(data=dt,aes(y = 0,
                              yend = `-log10P`,
                              x = 0,
                              xend = `log2FC`,
                              color = col 
                              ),alpha=0.2,linetype=6) +
    ggplot2::scale_colour_gradientn(colours = cf(n),guide=FALSE) +
    ggplot2::scale_x_continuous(labels=scaleFUN)
  
  p
}

#' @title Generate PLS-DA classification plot
#' @description Function to generate ggplot or plotly plot for PLS-DA classification
#' @param mSet mSet object
#' @param pls.type PLSDA type, Default: 'plsda'
#' @param cf Function to get plot colors from
#' @param pcs PCs used to classify, Default: 3
#' @return GGPLOT or PLOTLY object(s)
#' @seealso 
#'  \code{\link[ggplot2]{geom_bar}},\code{\link[ggplot2]{ggtheme}},\code{\link[ggplot2]{scale_manual}}
#' @rdname ggPlotClass
#' @export 
#' @importFrom ggplot2 geom_bar theme_minimal scale_fill_manual
ggPlotClass <- function(mSet,
                        pls.type = "plsda",
                        cf,
                        pcs = 3){
  res <- mSet$analSet$plsda$fit.info
  colnames(res) <- 1:ncol(res)
  # best.num <- mSet$analSet$plsda$best.num
  # choice <- mSet$analSet$plsda$choice
  df <- reshape2::melt(res)
  df$Component <- paste0("PC",df$Component)
  colnames(df) <- c("Metric", "Component", "Value")
  
  scaleFUN <- function(x) sprintf("%.2f", x)
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x=Metric, y=Value, fill=Metric)) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::theme_minimal() +
    ggplot2::facet_grid(~Component) +
    ggplot2::scale_fill_manual(values=cf(pcs)) +
    ggplot2::scale_y_continuous(labels=scaleFUN)
  p
}

#' @title Generate PLS-DA permutation plot
#' @description Function to generate ggplot or plotly plot for PLS-DA permutation
#' @param mSet mSet object
#' @param pls.type PLS-DA type, Default: 'plsda'
#' @param cf Function to get plot colors from
#' @param pcs PCs used for permutation, Default: 3
#' @return GGPLOT or PLOTLY object(s)
#' @seealso 
#'  \code{\link[stringr]{str_match}}
#'  \code{\link[ggplot2]{geom_freqpoly}},\code{\link[ggplot2]{scale_manual}},\code{\link[ggplot2]{geom_segment}},\code{\link[ggplot2]{geom_label}}
#' @rdname ggPlotPerm
#' @export 
#' @importFrom stringr str_match
#' @importFrom ggplot2 geom_histogram scale_fill_manual geom_segment geom_text
ggPlotPerm <- function(mSet,
                       pls.type = "plsda",
                       cf,
                       pcs = 3){
  bw.vec <- mSet$analSet$plsda$permut
  len <- length(bw.vec)
  df <- reshape2::melt(bw.vec)
  colnames(df) = "acc"
  # round p value
  pval <- mSet$analSet$plsda$permut.p
  rounded <- round(as.numeric(stringr::str_match(pval, "0\\.\\d*")), digits = 3)
  pval <- gsub(pval, pattern = "(0\\.\\d*)", replacement=rounded)
  # - - -
  scaleFUN <- function(x) sprintf("%.2f", x)
  
  p <- ggplot2::ggplot(df) +
    ggplot2::geom_histogram(mapping=ggplot2::aes(x=acc, y=..count.., fill=factor(..count..)),
                            binwidth=0.01) +
    ggplot2::scale_fill_manual(values=cf(20)) +
    
    ggplot2::labs(x="Accuracy", y = "Permutations") +
    ggplot2::geom_segment(data=df,
                          color="black",
                          x=bw.vec[1],
                          xend=bw.vec[1],
                          y=0,
                          ggplot2::aes(yend=.1*nrow(df)),
                          size=1.5,
                          linetype=8) +
    ggplot2::geom_text(mapping = ggplot2::aes(x = bw.vec[1], y =  .11*nrow(df), label = pval), color = "black", size = 4)+
    ggplot2::scale_x_continuous(labels=scaleFUN)
  
  p
}

ggPlotMLMistakes <- function(predictions,
                             labels,
                             test_sampnames, 
                             cutoffs,
                             covars,
                             metadata_focus = c(), cf=rainbow, 
                             show_reps=F,smooth_line=T){
  # miss metadata plot
  metadata_with_sample = c("sample", metadata_focus)
  test_meta = covars[sample %in% test_sampnames, ..metadata_with_sample]
  uniqvars = unique(test_meta[[2]])
  
  lvls = levels(labels[[1]])
  if(length(lvls) > 2){
    data = data.frame(text = "Only available for data with 2 groups!")
    ggplot2::ggplot(data) + ggplot2::geom_text(ggplot2::aes(label = text), x = 0.5, y = 0.5, size = 10) +
      ggplot2::theme(text = ggplot2::element_text(family = lcl$aes$font$family)) + ggplot2::theme_bw()
  }else{
    pred <- ROCR::prediction(lapply(predictions, function(l) l[[1]]), labels)
    cutoffs = pred@cutoffs
    predictions = pred@predictions
    
    all_reps = pbapply::pblapply(1:length(predictions), function(rep){
      predict_test = predictions[[rep]]#[[1]]
      labels = labels[[rep]]
      cutoffs = cutoffs[[rep]]
      classes = levels(labels)
      wrong_hits = lapply(cutoffs, function(x){
        predicted_labels = sapply(predict_test, function(y) if(y < x) classes[1] else classes[2])
        correct = predicted_labels == labels
        incorrect = test_sampnames[which(!correct)]
        tbl = data.table(mistaken = incorrect,
                         meta_var = test_meta[sample %in% incorrect, ..metadata_with_sample][[2]],
                         cutoff = rep(x, length(incorrect)))
        for(unique_var in uniqvars){
          ntotal = sum(test_meta[[metadata_focus]] == unique_var)
          tbl[meta_var == unique_var, "var_total"] <- ntotal
          in_missing = nrow(tbl[meta_var == unique_var])
          miss_perc = in_missing/tbl[meta_var == unique_var,"var_total"][[1]][1]*100
          tbl[meta_var == unique_var, "wrong_perc_var"] <- miss_perc  
        }
        tbl
      })
      res = data.table::rbindlist(wrong_hits)
      res$rep = rep
      res[cutoff != Inf] 
    })
    res = data.table::rbindlist(all_reps)
    
    line_fun = if(smooth_line) geom_smooth else geom_line
    p = ggplot2::ggplot(data = res,aes(x = cutoff, 
                                       y = wrong_perc_var,
                                       text = meta_var,
                                       color = meta_var)) + 
      #geom_point() +
      line_fun(cex = 1,se = FALSE ) +
      ggplot2::scale_color_manual(name = metadata_focus, 
                                  values=cf(length(unique(res$meta_var)))) +
      ggplot2::xlab("Cutoff") + ggplot2::ylab("% of testing mistakes")
    if(show_reps) p + facet_grid("rep") else p
  }
}
#' @title Generate ROC/PrecRec plot
#' @description Function to generate ggplot or plotly ROC/PrecRec plot for machine learning
#' @param data Model data
#' @param attempts Number of models that were created, Default: 50
#' @param cf Function to get plot colors from
#' @param curve_type roc or precrec. Default: 'roc'
#' @param class_type Binary or multivariate? (b/m), Default: 'b'
#' @return GGPLOT or PLOTLY object(s)
#' @seealso
#'  \code{\link[shiny]{showNotification}}
#'  \code{\link[pbapply]{pbapply}}
#'  \code{\link[ggplot2]{geom_path}},\code{\link[ggplot2]{annotate}},\code{\link[ggplot2]{stat_summary_bin}},\code{\link[ggplot2]{scale_manual}},\code{\link[ggplot2]{coord_fixed}},\code{\link[ggplot2]{coord_cartesian}}
#' @rdname ggPlotCurves
#' @export 
#' @importFrom shiny showNotification
#' @importFrom pbapply pbsapply
#' @importFrom ggplot2 geom_path annotate stat_summary_bin scale_color_manual coord_fixed coord_cartesian
ggPlotCurves <- function(perf.long,
                         #labels,
                         #predictions,
                         attempts = 50,
                         cf,
                         curve_type="roc",#roc, precrec
                         class_type="b",
                         forceorigin=F){
  perf.long <- data.table::as.data.table(perf.long)
  perf.long <- perf.long[metric == curve_type, -"metric"]
  mean.auc = mean(perf.long$AUC_AVG)
  means.per.comp = perf.long[, lapply(.SD, mean), by = comparison]
  ncomp = length(unique(means.per.comp$comparison))
  if(ncomp > 2){
    cols = cf(ncomp + 1)
    class_type = "m"
    try({
      shiny::showNotification("Calculating AUCs per comparison...")
    })
    for(comp in unique(perf.long$comparison)){
      new.name = paste0(comp, " || avg. AUC=", round(means.per.comp[comparison == comp]$AUC_PAIR, digits=3), " ||")
      perf.long[comparison == comp]$comparison <- c(new.name)
    } 
  }else{
    cols = cf(max(as.numeric(perf.long$attempt)))
    class_type = "b"
  }
  
  #if(curve_type == "precrec") forceorigin <- T
  
  if(forceorigin){
      origin_rows = lapply(1:length(unique(perf.long$attempt)), function(i){
        data.table::data.table(attempt = i, 
                               y = if(curve_type == "precrec") 1 else 0,
                               x = if(curve_type == "precrec") 0 else 1, 
                               AUC_AVG = unique(perf.long[attempt == i]$AUC_AVG),
                               AUC_PAIR = unique(perf.long[attempt == i]$AUC_PAIR),
                               comparison = unique(perf.long[attempt == i]$comparison))
      })
      perf.long = rbind(data.table::rbindlist(origin_rows),perf.long)
  }
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_path(data = perf.long, 
                       alpha=.5,
                       cex=.7,
                       ggplot2::aes(x = x,
                                    y = y,
                                    color = if(class_type == "m") comparison else as.factor(attempt),
                                    text = if(class_type == "m") comparison else as.factor(attempt),
                                    key = if(class_type == "m") comparison else as.factor(attempt),
                                    group = attempt)) +
    ggplot2::labs(color = if(class_type == "m") "Class" else "Attempt",
                  text = if(class_type == "m") "Class" else "Attempt",
                  key = if(class_type == "m") "Class" else "Attempt") +
    ggplot2::stat_summary_bin(data = perf.long,
                              #alpha=.6,
                              ggplot2::aes(x,
                                           y,
                                           group = comparison,
                                           color = if(class_type == "m") comparison else NULL),
                              fun=mean, geom="line",
                              cex = 2.3) + #,color="black") +
    ggplot2::scale_color_manual(values = cols) +
    #scale_y_continuous(limits=c(0,1)) +
    ggplot2::coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE) +
    ggplot2::coord_cartesian(xlim = c(.04,.96), ylim = c(.04,.96)) +
    ggplot2::xlab(if(curve_type == "roc") "FPR" else "Recall") + 
    ggplot2::ylab(if(curve_type == "roc") "TPR" else "Precision") +
    ggplot2::annotate("text",
                      label = paste0("Average AUC: ",
                                     format(mean.auc,
                                            2,
                                            drop0trailing = TRUE,
                                            digits = 2)),
                      size = 8,
                      x = 0.77,
                      y = 0.03)
  
  p
}

#' @title Generate machine learning importance bar plot
#' @description Function to generate ggplot or plotly variable importance barplot for machine learning
#' @param data Model data
#' @param attempts Number of models in data, Default: 50
#' @param cf Function to get plot colors from
#' @param topn Top number of compounds to display in plot, Default: 50
#' @param ml_name ML name as defined by user
#' @param ml_type ML model type
#' @return GGPLOT or PLOTLY object(s)
#' @seealso 
#'  \code{\link[Rmisc]{group.CI}}
#'  \code{\link[ggplot2]{geom_bar}},\code{\link[ggplot2]{scale_colour_gradient}},\code{\link[ggplot2]{theme}},\code{\link[ggplot2]{margin}},\code{\link[ggplot2]{geom_label}}
#' @rdname ggPlotBar
#' @export 
#' @importFrom Rmisc group.CI
#' @importFrom ggplot2 geom_bar scale_fill_gradientn theme element_blank geom_text
ggPlotBar <- function(data,
                      attempts=50,
                      cf,
                      topn=50,
                      ml_name,
                      ml_type){
  
  if(ml_name != ""){
    lname = ml_name
  }else{
    lname <- "all"
  }
  
  if(ml_type == "glmnet"){
    colnames(data) = c("m/z", "importance.mean", "dummy")
    data.ordered <- data[order(data$importance, decreasing=T),1:2]
  }else{
    data.norep <- data#[,-3]
    colnames(data.norep)[1] <- "m/z"
    data.ci = Rmisc::group.CI(importance ~ `m/z`, data.norep)
    data.ordered <- data.ci[order(data.ci$importance.mean, decreasing = T),]      
  }
  
  data.subset <- data.ordered[1:topn,]    
  data.subset$`m/z` <- gsub("`|^X","",data.subset$`m/z`)
  data.subset$`m/z` <- gsub("\\.$","-",data.subset$`m/z`)
  
  p <- ggplot2::ggplot(data.subset, 
                       ggplot2::aes(x = reorder(`m/z`, -importance.mean),
                                    y = importance.mean,
                                    fill = importance.mean,
                                    colour = importance.mean,
                                    text = `m/z`,
                                    key = `m/z`)) +
    ggplot2::geom_bar(stat = "identity") +
    #ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))+
    ggplot2::scale_fill_gradientn(colors=cf(20)) +
    ggplot2::labs(x="Top hits (m/z)",y=if(ml_type == "glmnet") "Times included in final model" else "Relative importance (%)")
  
  if(topn <= 15){
    p <- p + ggplot2::geom_text(ggplot2::aes(x=`m/z`, 
                                             y=importance.mean + .02*max(importance.mean), 
                                             label=substr(as.character(`m/z`),1,3)
    ), size = 4) + ggplot2::expand_limits(y=max(data.subset$importance.mean) + max(data.subset$importance.mean)*0.1)
  }
  
  mzdata <- p$data

  list(mzdata = mzdata, plot = p)
}

#' @title Generate PCA/PLS-DA loadings plot
#' @description Function to generate ggplot or plotly loadings plot for PCA
#' @param mSet mSet object
#' @param cf Function to get plot colors from
#' @param pcx X axis PC to use
#' @param pcy Y axis PC to use
#' @param type pca or plsda?, Default: 'pca'
#' @return GGPLOT or PLOTLY object(s)
#' @seealso 
#'  \code{\link[data.table]{as.data.table}}
#'  \code{\link[ggplot2]{geom_point}},\code{\link[ggplot2]{scale_continuous}},\code{\link[ggplot2]{scale_colour_gradient}}
#'  \code{\link[gsubfn]{fn}}
#' @rdname plotPCAloadings.2d
#' @export 
#' @importFrom data.table as.data.table
#' @importFrom ggplot2 geom_point scale_x_continuous scale_y_continuous scale_colour_gradientn
#' @importFrom gsubfn fn
plotPCAloadings.2d <- function(mSet,
                               cf,
                               pcx,
                               pcy,
                               type = "pca"){
  
  pcx = as.numeric(pcx)
  pcy = as.numeric(pcy)
  
  switch(type,
         pca = {
           df <- mSet$analSet$pca$rotation
           x.var <- round(mSet$analSet$pca$variance[pcx] * 100.00, digits=1)
           y.var <- round(mSet$analSet$pca$variance[pcy] * 100.00, digits=1)
         }, plsda = {
           plsda.table <- data.table::as.data.table(round(mSet$analSet$plsr$Xvar
                                                          / mSet$analSet$plsr$Xtotvar
                                                          * 100.0,
                                                          digits = 2),
                                                    keep.rownames = T)
           colnames(plsda.table) <- c("PC", "var")
           plsda.table[, "PC"] <- paste0("PC", 1:nrow(plsda.table))
           
           x.var <- plsda.table[PC == paste0("PC",pcx)]$var
           y.var <- plsda.table[PC == paste0("PC",pcy)]$var
           
           # --- coordinates ---
           df <- mSet$analSet$plsr$loadings
           class(df) <- "matrix"
           colnames(df) <- paste0("PC", 1:ncol(df))
         })
  
  df = as.data.frame(df)

  df$extremity <- apply(df, 1, function(row) max(abs(c(row[[pcx]], 
                                                       row[[pcy]]))))
  
  scaleFUN <- function(x) sprintf("%.4s", x)
  
  prefix = switch(type, 
                  pca = "PC",
                  plsda = "Component ")
  
  p <- ggplot2::ggplot(df, ggplot2::aes(.data[[paste0("PC",pcx)]], .data[[paste0("PC",pcy)]])) +
    ggplot2::geom_point(ggplot2::aes(color = extremity,
                                     size = extremity,
                                     text = rownames(df),
                                     key = rownames(df)),
                        #pch=21, size = 2, stroke = 2,
                        #fill="white", 
                        alpha=0.7)+
    
    ggplot2::scale_size_area(max_size = 15) +
    ggplot2::scale_x_continuous(labels=scaleFUN,name=gsubfn::fn$paste("$prefix$pcx ($x.var%)")) +
    ggplot2::scale_y_continuous(labels=scaleFUN,name=gsubfn::fn$paste("$prefix$pcy ($y.var%)")) +
    ggplot2::scale_colour_gradientn(colors=cf(20))
  #ggplot2::scale_y_discrete(labels=scaleFUN) +
  #ggplot2::scale_x_discrete(labels=scaleFUN)
  p 
}

#' @title Generate PCA/PLS-DA loadings plot
#' @description Function to generate ggplot or plotly loadings plot for PCA
#' @param mSet mSet object
#' @param cf Function to get plot colors from
#' @param pcx X axis PC to use
#' @param pcy Y axis PC to use
#' @param pcz Z axis PC to use
#' @param font Font family to use in plotly plot
#' @param type pca or plsda?, Default: 'pca'
#' @return GGPLOT or PLOTLY object(s)
#' @seealso 
#'  \code{\link[data.table]{as.data.table}}
#'  \code{\link[gsubfn]{fn}}
#' @rdname plotPCAloadings.3d
#' @export 
#' @importFrom data.table as.data.table
#' @importFrom gsubfn fn
plotPCAloadings.3d <- function(mSet,
                               cf,
                               pcx,
                               pcy,
                               pcz,
                               font,
                               type = "pca"){
  
  pcx <- as.numeric(pcx)
  pcy <- as.numeric(pcy)
  pcz <- as.numeric(pcz)
  
  switch(type,
         pca = {
           df <- mSet$analSet$pca$rotation
           x.var <- round(mSet$analSet$pca$variance[pcx] * 100.00, digits=1)
           y.var <- round(mSet$analSet$pca$variance[pcy] * 100.00, digits=1)
           z.var <- round(mSet$analSet$pca$variance[pcz] * 100.00, digits=1)
         }, plsda = {
           plsda.table <- data.table::as.data.table(round(mSet$analSet$plsr$Xvar
                                                          / mSet$analSet$plsr$Xtotvar
                                                          * 100.0,
                                                          digits = 2),
                                                    keep.rownames = T)
           colnames(plsda.table) <- c("PC", "var")
           plsda.table[, "PC"] <- paste0("PC", 1:nrow(plsda.table))
           
           x.var <- plsda.table[PC == paste0("PC",pcx)]$var
           y.var <- plsda.table[PC == paste0("PC",pcy)]$var
           z.var <- plsda.table[PC == paste0("PC",pcz)]$var
           
           # --- coordinates ---
           df <- mSet$analSet$plsr$loadings
           class(df) <- "matrix"
           colnames(df) <- paste0("PC", 1:ncol(df))
         })
  
  df <- as.data.frame(df)
  
  df$extremity <- apply(df, 1, function(row) max(abs(c(row[[pcx]], 
                                                       row[[pcy]],
                                                       row[[pcz]]))))
  
  basic_scene = list(
    aspectmode="cube",
    aspectratio=list(x=1,y=1,z=1),
    hoverlabel = list(bgcolor = ~extremity),
    camera = list(
      eye = list(x=0, y=0, z= 2)
    ),
    xaxis = list(
      titlefont = list(size = font$ax.txt.size * 1.5),
      title = gsubfn::fn$paste("$pcx ($x.var%)")),
    yaxis = list(
      titlefont = list(size = font$ax.txt.size * 1.5),
      title = gsubfn::fn$paste("$pcy ($y.var%)")),
    zaxis = list(
      titlefont = list(size = font$ax.txt.size * 1.5),
      title = gsubfn::fn$paste("$pcz ($z.var%)"))) 
  
  cols = cf(8)
  bins = seq(0, max(df$extremity), length.out = 8)
  bins = bins/max(bins)
  colscale <- lapply(1:8, function(i) c(bins[i], cols[i]))
  
  p <- plot_ly(df, 
               x = df[,pcx], 
               y = df[,pcy], 
               z = df[,pcz],
               key = rownames(df),
               text = rownames(df),
               hoverinfo = "text", 
               marker = list(size = ~extremity * 800,
                             sizemode = "diameter",
                             sizes = c(5, 100),
                             symbol = "circle",
                             opacity = 1,
                             color = ~extremity,
                             colorscale = list(bins, cols),
                             line = list(width = 0.1,
                                         color = "white"),
                             showscale = FALSE)
  ) %>%
    add_markers() %>%
    layout(scene = basic_scene)
  p
}


#' @title Generate 3d scatter plot
#' @description Function to generate ggplot or plotly plot for PCA/PLS-DA/T-SNE
#' @param mSet mSet object
#' @param cols Colors to use
#' @param shape.fac Marker shape based on which metadata column?, Default: 'label'
#' @param pcx X component
#' @param pcy Y component
#' @param pcz Z component
#' @param type pca, plsda or tsne?, Default: 'pca'
#' @param font Font family to use in plot
#' @param col.fac Marker fill based on which metadata column?, Default: 'label'
#' @param mode normal or timeseries mode?, Default: 'normal'
#' @param cf Function to get plot colors from
#' @return GGPLOT or PLOTLY object(s)
#' @seealso 
#'  \code{\link[data.table]{as.data.table}}
#'  \code{\link[plotly]{plot_ly}}
#'  \code{\link[rgl]{ellipse3d}}
#'  \code{\link[gsubfn]{fn}}
#' @rdname plotPCA.3d
#' @export 
#' @importFrom data.table as.data.table
#' @importFrom plotly plot_ly
#' @importFrom rgl ellipse3d
#' @importFrom gsubfn fn
plotPCA.3d <- function(mSet,
                       cols,
                       shape.fac="label",
                       pcx, pcy, pcz,
                       type="pca",font,
                       col.fac = "label",
                       fill.fac = "label",
                       mode="normal",cf,
                       ellipse=T){
  
  pcx = as.numeric(pcx)
  pcy = as.numeric(pcy)
  pcz = as.numeric(pcz)
  
  switch(type,
         ica = {
           df = mSet$analSet$ica$S
           x.var = ""
           y.var = ""
           z.var = ""
         },
         umap = {
           df = mSet$analSet$umap$layout
           x.var = ""
           y.var = ""
           z.var = ""
         },
         tsne = {
           df = mSet$analSet$tsne$x
           x.var = ""
           y.var = ""
           z.var = ""
         },
         pca = {
           df <- mSet$analSet$pca$x
           x.var <- round(mSet$analSet$pca$variance[pcx] * 100.00, digits=1)
           y.var <- round(mSet$analSet$pca$variance[pcy] * 100.00, digits=1)
           z.var <- round(mSet$analSet$pca$variance[pcz] * 100.00, digits=1)
         }, plsda = {
           plsda.table <- data.table::as.data.table(round(mSet$analSet$plsr$Xvar
                                                          / mSet$analSet$plsr$Xtotvar
                                                          * 100.0,
                                                          digits = 2),
                                                    keep.rownames = T)
           colnames(plsda.table) <- c("PC", "var")
           plsda.table[, "PC"] <- paste0("PC", 1:nrow(plsda.table))
           
           x.var <- round(plsda.table[PC == pcx]$var, digits=1)
           y.var <- round(plsda.table[PC == pcy]$var, digits=1)
           z.var <- round(plsda.table[PC == pcz]$var, digits=1)
           
           # --- coordinates ---
           df <- mSet$analSet$plsr$scores
           class(df) <- "matrix"
           colnames(df) <- paste0("PC", 1:ncol(df))
         })
  
  df <- as.data.frame(df)
  rownames(df) <- rownames(mSet$dataSet$norm)
  
  if(mode != "normal"){
    fac.lvls <- length(levels(mSet$dataSet$facA))
    classes = mSet$dataSet$facA
    df_list <- split(df, mSet$dataSet$facB)
  }else{
    fac.lvls <- length(levels(mSet$dataSet$cls))
    classes = mSet$dataSet$cls
    df_list <- list(df)
  }
  
  cols <- if(is.null(cols)) cf(length(levels(classes))) else{
    if(length(cols) < length(levels(classes))){
      cols <- cf(levels(classes))
    }
    cols
  }
  
  cols <- if(length(unique(classes)) > length(cols)){
    cols <- cf(length(unique(classes)))
  }else{
    cols[c(1:length(unique(classes)))]
  }
  
  symbol.vec<-if(is.null(shape.fac)){
    rep('circle', times = length(classes))
  }else if(shape.fac == "label"){
    rep('circle', times = length(classes))
  }else{
    as.factor(mSet$dataSet$covars[, ..shape.fac][[1]])
  }
  
  col.vec <- if(is.null(col.fac)){
    classes
  }else if(shape.fac == "label"){
    classes
  }else{
    as.factor(mSet$dataSet$covars[, ..col.fac][[1]])
  }
  
  fill.vec <- if(is.null(fill.fac)){
    classes
  }else if(fill.fac == "label"){
    classes
  }else{
    as.factor(mSet$dataSet$covars[, ..fill.fac][[1]])
  }
  
  plots_facet <- lapply(1:length(df_list), function(i){
    
    df = df_list[[i]]
    
    orig_idx = match(rownames(df), rownames(mSet$dataSet$norm))
    
    plots <- plotly::plot_ly(scene = paste0("scene", if(i > 1) i else ""))
    
    if(ellipse){
      
      show.orbs <- c(1:length(levels(classes)))
      
      for(class in levels(classes)){
        samps <- rownames(mSet$dataSet$norm)[which(classes == class)]
        row = which(rownames(df) %in% samps)
        # ---------------------
        xc=df[row, pcx]
        yc=df[row, pcy]
        zc=df[row, pcz]
        
        # --- plot ellipse ---
        worked = F
        
        try({
          o <- rgl::ellipse3d(cov(cbind(xc,yc,zc)),
                              centre=c(mean(xc),
                                       mean(yc),
                                       mean(zc)),
                              level = 0.95)
          worked = T
        })
        
        if(worked){
          mesh <- c(list(x = o$vb[1, o$ib]/o$vb[4, o$ib],
                         y = o$vb[2, o$ib]/o$vb[4, o$ib],
                         z = o$vb[3, o$ib]/o$vb[4, o$ib]))
          plots = plots %>% add_mesh(
            x=mesh$x,
            y=mesh$y,
            z=mesh$z,
            type='mesh3d',
            alphahull = 0,
            opacity=0.1
          )
          
          adj_plot <- plotly_build(plots)
          rgbcols <- toRGB(cols[show.orbs])
          c = 1
          
          for(i in seq_along(adj_plot$x$data)){
            item = adj_plot$x$data[[i]]
            if(item$type == "mesh3d"){
              adj_plot$x$data[[i]]$color <- rgbcols[c]
              adj_plot$x$data[[i]]$visible <- TRUE
              #adj_plot$x$data[[i]]$hoverinfo <- "none"
              c = c + 1
            }
          }
        }else{
          adj_plot = plots
        }
        show.orbs <- c(show.orbs, worked)
      }       
    }else{
      adj_plot = plots
    }
    
    t <- list(family = font$family)
    
    df$shape <- symbol.vec[orig_idx]
    df$fill <- fill.vec[orig_idx]
    df$color <- col.vec[orig_idx]
    df$x <- df[,pcx]
    df$y <- df[,pcy]
    df$z <- df[,pcz]
    
    x = as.numeric(df$color)
    limits=range(x)
    pal=cf(200)
    outlineCols = pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
    
    # --- return ---
    pca_plot <- adj_plot %>% add_trace(
      data = df,
      hoverinfo = 'text',
      text = rownames(df),
      x = ~x,
      y = ~y,
      z = ~z,
      visible = rep(T, times=fac.lvls),
      type = "scatter3d",
      color = ~fill,
      colors = cols,
      opacity = 1,
      symbol = ~shape,
      symbols = c("circle", "square", "diamond", 
                  "cross", "x", "triangle-up", 
                  "pentagon", "hexagram", "star", 
                  "diamond", "hourglass", "bowtie",
                  "asterisk", "hash", "y","line"),
      marker = list(
        line = list(
          width = 1.5,
          color = outlineCols
        )
      )
    ) 
    # --- return ---
    pca_plot
  })
  
  title_prefix = switch(type,
                        tsne = "t-sne dimension ",
                        umap = "umap dimension ",
                        ica = "IC",
                        pca = "PC",
                        plsda = "Component ")
  
  basic_scene = list(
    aspectmode="cube",
    aspectratio=list(x=1,y=1,z=1),
    camera = list(
      eye = list(x=0, y=0, z= 2)
    ),
    xaxis = list(
      titlefont = list(size = font$ax.txt.size * 1.5),
      title = paste0(title_prefix, pcx, if(x.var != "") gsubfn::fn$paste("($x.var%)") else "")),
    yaxis = list(
      titlefont = list(size = font$ax.txt.size * 1.5),
      title =  paste0(title_prefix, pcy, if(y.var != "") gsubfn::fn$paste("($y.var%)") else "")),
    zaxis = list(
      titlefont = list(size = font$ax.txt.size * 1.5),
      title = paste0(title_prefix, pcz, if(z.var != "") gsubfn::fn$paste("($z.var%)") else ""))
  ) 
  if(mode == "normal"){
    plots_facet[[1]] %>% layout(font = t, 
                                scene = basic_scene) %>%
      config(toImageButtonOptions = list(format = "svg"))
  }else{
    maxrows = ceiling(length(plots_facet)/2)
    
    x_start = rep(c(0, 0.5), maxrows)
    x_end = rep(c(0.5, 1), maxrows)
    
    y_start = rev(c(0, 0, rep(sapply(1:(maxrows-1), function(i) c(1/maxrows)*i), each=2)))
    y_end = rev(rep(sapply(1:(maxrows), function(i) c(1/maxrows)*i), each=2))
    #y_end = c(0.5, 0.5,1, 1)
    
    domains = lapply(1:10, function(i){
      list(x = c(x_start[i], x_end[i]),
           y = c(y_start[i], y_end[i])) 
    })
    
    # TODO: make this a less ugly solution ; w;"
    subplot(plots_facet) %>% layout(font = t,
                                    scene = append(basic_scene,
                                                   list(domain=domains[[1]])),
                                    scene2 = append(basic_scene,
                                                    list(domain=domains[[2]])),
                                    scene3 = append(basic_scene,
                                                    list(domain=domains[[3]])),
                                    scene4 = append(basic_scene,
                                                    list(domain=domains[[4]])),
                                    scene5 = append(basic_scene,
                                                    list(domain=domains[[5]])),
                                    scene6 = append(basic_scene,
                                                    list(domain=domains[[6]])),
                                    scene7 = append(basic_scene,
                                                    list(domain=domains[[7]])),
                                    scene8 = append(basic_scene,
                                                    list(domain=domains[[8]])),
                                    scene9 = append(basic_scene,
                                                    list(domain=domains[[9]])),
                                    scene10 = append(basic_scene,
                                                     list(domain=domains[[10]]))) %>%
      config(toImageButtonOptions = list(format = "svg"))
  }
}



#' @title Generate 2d scatter plot
#' @description Function to generate ggplot or plotly plot for PCA/PLS-DA/T-SNE
#' @param mSet mSet object
#' @param cols Colors to use
#' @param shape.fac Marker shape based on which metadata column?, Default: 'label'
#' @param pcx X component
#' @param pcy Y component
#' @param type pca, plsda or tsne?, Default: 'pca'
#' @param col.fac Marker fill based on which metadata column?, Default: 'label'
#' @param mode normal or timeseries mode?, Default: 'normal'
#' @param cf Function to get plot colors from
#' @return GGPLOT or PLOTLY object(s)
#' @seealso 
#'  \code{\link[data.table]{as.data.table}}
#'  \code{\link[ggplot2]{ggplot}},\code{\link[ggplot2]{geom_point}},\code{\link[ggplot2]{stat_ellipse}},\code{\link[ggplot2]{scale_continuous}},\code{\link[ggplot2]{scale_manual}},\code{\link[ggplot2]{labs}}
#'  \code{\link[gsubfn]{fn}}
#'  \code{\link[Hmisc]{capitalize}}
#' @rdname plotPCA.2d
#' @export 
#' @importFrom data.table as.data.table
#' @importFrom ggplot2 ggplot geom_point stat_ellipse scale_x_continuous scale_y_continuous scale_fill_manual scale_color_manual ggtitle
#' @importFrom gsubfn fn
#' @importFrom Hmisc capitalize
plotPCA.2d <- function(mSet, shape.fac = "label", cols, col.fac = "label",  fill.fac = "label",
                       pcx, pcy, mode="normal", type="pca",
                       cf = rainbow, ellipse=T){
  
  classes <- if(mode == "ipca"){
    mSet$dataSet$facA
  }else{
    mSet$dataSet$cls
  }
  
  pcx = as.numeric(pcx)
  pcy = as.numeric(pcy)
  
  switch(type,
         ica = {
           df <- mSet$analSet$ica$S
           x.var <- ""
           y.var <- ""
           fac.lvls <- length(levels(mSet$dataSet$cls))
           
           xc=df[, pcx]
           yc=df[, pcy]
           
           dat_long <- data.table::data.table(variable = rownames(mSet$dataSet$norm),
                                              group = classes,
                                              x = xc,
                                              y = yc)
         },
         umap = {
           df <- mSet$analSet$umap$layout
           x.var <- ""
           y.var <- ""
           fac.lvls <- length(levels(mSet$dataSet$cls))
           
           xc=df[, pcx]
           yc=df[, pcy]
           
           dat_long <- data.table::data.table(variable = rownames(mSet$dataSet$norm),
                                              group = classes,
                                              x = xc,
                                              y = yc)
         },
         tsne = {
           df <- mSet$analSet$tsne$x
           x.var <- ""
           y.var <- ""
           fac.lvls <- length(levels(mSet$dataSet$cls))
           
           xc=df[, pcx]
           yc=df[, pcy]
           
           dat_long <- data.table::data.table(variable = rownames(mSet$dataSet$norm),
                                              group = classes,
                                              x = xc,
                                              y = yc)
         },
         pca = {
           df <- mSet$analSet$pca$x
           x.var <- round(mSet$analSet$pca$variance[pcx] * 100.00, digits=1)
           y.var <- round(mSet$analSet$pca$variance[pcy] * 100.00, digits=1)
           fac.lvls <- length(levels(mSet$dataSet$cls))
           
           xc=mSet$analSet$pca$x[, pcx]
           yc=mSet$analSet$pca$x[, pcy]
           
           dat_long <- data.table(variable = names(xc),
                                  group = classes,
                                  x = xc,
                                  y = yc)
         },
         plsda = {
           plsda.table <- data.table::as.data.table(round(mSet$analSet$plsr$Xvar
                                                          / mSet$analSet$plsr$Xtotvar
                                                          * 100.0,
                                                          digits = 2),
                                                    keep.rownames = T)
           colnames(plsda.table) <- c("PC", "var")
           plsda.table[, "PC"] <- paste0("PC", 1:nrow(plsda.table))
           
           x.var <- plsda.table[PC == paste0("PC", pcx)]$var
           y.var <- plsda.table[PC == paste0("PC", pcy)]$var
           
           # --- coordinates ---
           df <- mSet$analSet$plsr$scores
           colnames(df) <- paste0("PC", 1:ncol(df))
           rownames(df) <- rownames(mSet$dataSet$norm)
           
           xc=df[, pcx]
           yc=df[, pcy]
           
           dat_long <- data.table::data.table(variable = names(xc),
                                              group = classes,
                                              x = xc,
                                              y = yc)
         })
  
  if(mode == "ipca"){
    fac.lvls <- length(levels(mSet$dataSet$facA))
    dat_long$groupB <- mSet$dataSet$facB
  }
  
  dat_long$fill <- if(is.null(fill.fac)){
    dat_long$group
  } else if(fill.fac == "label"){
    dat_long$group
  }else{
    as.factor(mSet$dataSet$covars[,..fill.fac][[1]])
  }
  
  dat_long$color <- if(is.null(col.fac)){
    factor(1) # all same shape...
  } else if(col.fac == "label"){
    dat_long$group
  }else{
    as.factor(mSet$dataSet$covars[,..col.fac][[1]])
  }
  
  symbols = c(1:25)
  
  adjcols = cf(200)
  
  dat_long$shape <- if(is.null(shape.fac)){
    factor(1) # all same shape...
  } else if(shape.fac == "label"){
    dat_long$group
  }else{
    shapes = as.factor(mSet$dataSet$covars[,..shape.fac][[1]])
    lvls = unique(shapes)
    if(length(lvls) > 5){
      print(">5 shapes! This will sacrifice point outline as variable.")
      dat_long$color <- dat_long$fill
      adjcols = cols
    }else{
      symbols <- c(21:25) 
    }
    shapes
  }
  
  cols <- if(is.null(cols)) cf(length(levels(classes))) else{
    if(length(cols) < length(levels(classes))){
      cols <- cf(levels(classes))
    }
    cols
  }
  
  ggplot2::scale_shape_manual(values = symbols)
  
  title_prefix = switch(type,
                  tsne = "t-sne dimension ",
                  umap = "umap dimension ",
                  ica = "IC",
                  pca = "PC",
                  plsda = "Component ")
  
  p <- ggplot2::ggplot(dat_long, ggplot2::aes(x, y,group=group)) +
    ggplot2::geom_point(size=6, ggplot2::aes(
      shape=shape,
      text=variable,
      fill=fill,
      color=color), alpha=0.7,stroke = 1.5)+
    ggplot2::scale_x_continuous(name = paste0(title_prefix, pcx, if(x.var != "") gsubfn::fn$paste("($x.var%)") else ""))+
    ggplot2::scale_y_continuous(name = paste0(title_prefix, pcy, if(y.var != "") gsubfn::fn$paste("($y.var%)") else ""))+
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::scale_shape_manual(values = as.numeric(symbols)) +
    ggplot2::guides(fill = ggplot2::guide_legend(fill = ggplot2::guide_legend(override.aes = list(shape = 21)),
                                                 color = ggplot2::guide_legend(override.aes = list(shape = 21))))
  
  if(ellipse){
    p <- p + ggplot2::stat_ellipse(geom = "polygon", ggplot2::aes(fill=group), alpha = 0.3,level = .95, type = "norm",)
  }
  
  if(mode == "ipca"){
    p <- p + ggplot2::facet_wrap(~groupB,ncol = 2)
    p <- p + ggplot2::ggtitle(Hmisc::capitalize(mSet$dataSet$facB.lbl))
  }
  p
}

#' @title Generate Venn plot
#' @description Function to generate ggplot or plotly plot for Venn diagram
#' @param mSet mSet object
#' @param venn_yes Table with data-subsets to include in venn diagram
#' @param top Top x m/z per category chosen for intersection, Default: 100
#' @param cols Colors to use
#' @param cf Function to get plot colors from
#' @return GGPLOT or PLOTLY object(s)
#' @seealso 
#'  \code{\link[stringr]{str_split}}
#'  \code{\link[ggVennDiagram]{ggVennDiagram}}
#'  \code{\link[ggplot2]{lims}}
#' @rdname ggPlotVenn
#' @export 
#' @importFrom stringr str_split
#' @importFrom ggVennDiagram ggVennDiagram
#' @importFrom ggplot2 lims
ggPlotVenn <- function(mSet,
                       venn_yes,
                       top = 100,
                       cols,
                       cf){
  
  flattened <- getTopHits(mSet, unlist(venn_yes$now), top)
  
  parseFun = function(labels){
    sapply(labels, function(label){
      if(grepl(label, pattern=":")){
        split = trimws(stringr::str_split(label, ":")[[1]])
        stats_on = toupper(split[1])
        which_stats = split[length(split)]
        subset = grep("=", split, value=T)
        paste0(stats_on, "\n", if(length(subset)>0) paste0(subset, collapse="\n", "\n") else "", paste0(">> ", which_stats, " <<\n\n"))
      }else{
        label
      }
    })
  }
  
  p = ggVennDiagram::ggVennDiagram(flattened,
                                   label_alpha = 1, 
                                   cf = cf,
                                   parseFun = parseFun) + ggplot2::lims(x = switch(as.character(length(flattened)),
                                                                                   "2"= c(-7,11),
                                                                                   "3"= c(-5,9),
                                                                                   "4"= c(-.05,1.13)),
                                                                        y = switch(as.character(length(flattened)),
                                                                                   "2"= c(-5.5,7.5),
                                                                                   "3"= c(-4,10),
                                                                                   "4"= c(.15,1)))
  list(plot = p, info = flattened)
  
}

#' @title Generate PCA Scree plot
#' @description Function to generate ggplot or plotly scree plot forPCA
#' @param mSet mSet object
#' @param cf Function to get plot colors from
#' @param pcs Principal components displayed, Default: 20
#' @return GGPLOT or PLOTLY object(s)
#' @seealso 
#'  \code{\link[ggplot2]{ggplot}},\code{\link[ggplot2]{geom_path}},\code{\link[ggplot2]{scale_colour_gradient}}
#' @rdname ggPlotScree
#' @export 
#' @importFrom data.table data.table
#' @importFrom ggplot2 ggplot geom_line scale_colour_gradientn
ggPlotScree <- function(mSet, cf, pcs=20){
  df <- data.table::data.table(
    pc = 1:length(names(mSet$analSet$pca$variance)),
    var = round(mSet$analSet$pca$variance*100,digits = 1))
  p <- ggplot2::ggplot(data=df[1:20,]) + 
    ggplot2::geom_line(mapping = ggplot2::aes(x=pc, y=var), cex=1, color="black") +
    ggplot2::geom_point(mapping = ggplot2::aes(x=pc, y=var, color=var), cex=3) +
    ggplot2::scale_colour_gradientn(colours = cf(200)) + 
    ggplot2::xlab("Principal components") +
    ggplot2::ylab("% Variance")
  # - - - - -
  p
}

#' @title Generate wordcloud bar plot
#' @description Function to generate ggplot or plotly barplot for word cloud
#' @param wcdata Word cloud data
#' @param cf Function to get plot colors from
#' @param plotlyfy Convert plot to plotly object?, Default: T
#' @return GGPLOT or PLOTLY object(s)
#' @seealso 
#'  \code{\link[ggplot2]{geom_bar}},\code{\link[ggplot2]{coord_flip}},\code{\link[ggplot2]{scale_colour_gradient}},\code{\link[ggplot2]{theme}},\code{\link[ggplot2]{margin}}
#' @rdname ggPlotWordBar
#' @export 
#' @importFrom ggplot2 geom_bar coord_flip scale_fill_gradientn theme element_blank
ggPlotWordBar <- function(wcdata, cf, plotlyfy=T){
  g <- ggplot2::ggplot(wcdata, ggplot2::aes(y = freq, x = reorder(word, 
                                                                  freq, 
                                                                  sum)))
  g <- g + ggplot2::geom_bar(ggplot2::aes(fill = freq),
                             stat = "identity") +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_gradientn(colors=cf(256)) +
    ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                   axis.ticks.x=ggplot2::element_blank()) +
    ggplot2::labs(x="Word",y="Frequency")
  g
}

#' @title Generate power plot
#' @description Function to generate ggplot or plotly plot for power analysis
#' @param mSet mSet object
#' @param cf Function to get plot colors from
#' @param comparisons Which 1 vs 1 classes to use?
#' @param max_samples Max amount of samples simulated per group
#' @return GGPLOT or PLOTLY object(s)
#' @seealso 
#'  \code{\link[data.table]{rbindlist}},\code{\link[data.table]{data.table-package}}
#'  \code{\link[ggplot2]{geom_path}},\code{\link[ggplot2]{stat_summary_bin}},\code{\link[ggplot2]{coord_cartesian}}
#' @rdname ggPlotPower
#' @export 
#' @importFrom data.table rbindlist data.table
#' @importFrom ggplot2 geom_path stat_summary_bin coord_cartesian
ggPlotPower <- function(mSet,
                        cf,
                        comparisons,
                        max_samples){
  
  cols = cf(length(comparisons))
  
  data = data.table::rbindlist(lapply(comparisons, function(comp) data.table::data.table(samples = mSet$analSet$power[[comp]]$Jpred,
                                                                                         power = mSet$analSet$power[[comp]]$pwrD,
                                                                                         comparison = c(comp))))
  
  #data$comparison <- substr(gsub(data$comparison, pattern = " .*$", replacement = ""), 1, 10)
  
  if(ncol(data) == 1){
    stop("Something went wrong! Try other settings please :(")
  }else{
    p <- ggplot2::ggplot(data, ggplot2::aes(x=samples,y=power)) +
      ggplot2::geom_path(alpha=.5,
                         cex=.5,
                         ggplot2::aes(color = comparison, group = comparison)) +
      ggplot2::stat_summary_bin(#alpha=.6,
        ggplot2::aes(samples, 
                     power, 
                     group=comparison), 
        fun=mean, geom="line", 
        cex = 2.3,color="black") +
      ggplot2::stat_summary_bin(#alpha=.6,
        ggplot2::aes(samples, power, 
                     color=comparison
                     #,group=comparison
        ), 
        fun=mean, geom="line", 
        cex = 1.2) +
      ggplot2::stat_summary_bin(ggplot2::aes(samples,
                                             power), 
                                fun=mean, color="black", 
                                geom="line", cex = 2)# +
    # ggplot2::coord_cartesian(xlim = c(0,max_samples), 
    #                          ylim = c(.04,.96))
    p 
  }
}

#' @title Generate MUMMICHOG plot
#' @description Function to generate ggplot or plotly plot for MUMMICHOG
#' @param mum_mSet mSet object
#' @param anal.type Mummichog or GSEA? , Default: 'mummichog'
#' @param cf Function to get plot colors from
#' @return GGPLOT or PLOTLY object(s)
#' @seealso 
#'  \code{\link[ggplot2]{labs}},\code{\link[ggplot2]{scale_colour_gradient}}
#' @rdname ggPlotMummi
#' @export 
#' @importFrom ggplot2 ylab xlab scale_colour_gradientn
ggPlotMummi <- function(mSet, cf){
  
  anal.type = if("mummi.resmat" %in% names(mSet$analSet$enrich)) "mummichog" else "gsea"
  
  if (anal.type == "mummichog") {
    mummi.mat <- mSet$analSet$enrich$mummi.resmat
    y <- -log10(mummi.mat[, 5])
    x <- mummi.mat[, 3]/mummi.mat[, 4]
    pathnames <- rownames(mummi.mat)
  } else {
    gsea.mat <- mSet$analSet$enrich$mummi.gsea.resmat
    if(is.null(gsea.mat)) stop("No hits found.")
    y <- -log10(gsea.mat[, 3])
    x <- gsea.mat[, 2]/gsea.mat[, 1]
    pathnames <- rownames(gsea.mat)
  }
  inx <- order(y, decreasing = T)
  y <- y[inx]
  x <- x[inx]
  path.nms <- pathnames[inx]
  sqx <- sqrt(x)
  min.x <- min(sqx, na.rm = TRUE)
  max.x <- max(sqx, na.rm = TRUE)
  if (min.x == max.x) {
    max.x = 1.5 * max.x
    min.x = 0.5 * min.x
  }
  maxR <- (max.x - min.x)/40
  minR <- (max.x - min.x)/160
  radi.vec <- minR + (maxR - minR) * (sqx - min.x)/(max.x - 
                                                      min.x)
  bg.vec <- heat.colors(length(y))
  
  df <- data.frame(path.nms, x, y)
  
  scaleFUN <- function(x) sprintf("%.2f", x)
  
  p = ggplot2::ggplot(df) + ggplot2::geom_point(ggplot2::aes(y = y, 
                                                             x = x, 
                                                             size = `radi.vec`, 
                                                             color = `radi.vec`,
                                                             text = path.nms,
                                                             key = path.nms)) +
    # ggtitle("Enrichment Results") +
    ggplot2::ylab("-log10(p)") + 
    ggplot2::xlab("Significant/expected hits") +
    ggplot2::scale_colour_gradientn(colours = cf(20)) +
    ggplot2::scale_y_continuous(labels=scaleFUN)
  
  p
}
