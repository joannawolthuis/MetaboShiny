#' @export
#' Plot a summary of metaboanalystR normalization results. Takes the total of 20 m/z values and 20 samples before and after normalization and plots that distribution.
#' @param mSet input user mSet
#' @param plot.theme function for ggplot theme used.
#' @return list of four plots that fit in a 2x2 raster used in MetaboShiny.
ggplotNormSummary <- function(mSet,
                              plot.theme,
                              font,
                              cf){
  
  # load in original data (pre-normalization, post-filter)
  orig_data <- as.data.frame(mSet$dataSet$prenorm)
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
  plot <- ggplot2::ggplot(data=orig_melt) +
    plot.theme(base_size = 15) #+ facet_grid(. ~ variable)
  
  # first result plot: is a density plot of chosen 20 mz values with 20 samples
  RES1 <- plot + ggplot2::geom_density(ggplot2::aes(x=value), colour="blue", fill="blue", alpha=0.4)+
    ggplot2::theme(legend.position="none",
                   axis.text=ggplot2::element_text(size=font$ax.num.size),
                   axis.title=ggplot2::element_text(size=font$ax.txt.size),
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.line = ggplot2::element_line(colour = 'black', size = .5),
                   text = ggplot2::element_text(family = font$family))
  
  # second result plot: shows the spread of the intensities before normalization
  RES2 <- plot + ggplot2::geom_boxplot(
    ggplot2::aes(x=value,y=variable),
    color=cf(sampsize),
    alpha=0.4) + ggplot2::geom_vline(ggplot2::aes(xintercept=median(orig_melt$value)))+
    ggplot2::theme(legend.position="none",
                   axis.text=ggplot2::element_text(size=font$ax.num.size),
                   axis.title=ggplot2::element_text(size=font$ax.txt.size),
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.line = ggplot2::element_line(colour = 'black', size = .5),
                   text = ggplot2::element_text(family = font$family))
  
  # create base plot with base theme and font size for normalized data
  plot <- ggplot2::ggplot(data=norm_melt) +
    plot.theme(base_size = 15) #+ facet_grid(. ~ variable)
  
  # third result plot: a density plot of chosen 20 mz values post normalization
  RES3 <- plot + ggplot2::geom_density(ggplot2::aes(x=value), colour="pink", fill="pink", alpha=0.4)+
    ggplot2::theme(legend.position="none",
                   axis.text=ggplot2::element_text(size=font$ax.num.size),
                   axis.title=ggplot2::element_text(size=font$ax.txt.size),
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.line = ggplot2::element_line(colour = 'black', size = .5),
                   text = ggplot2::element_text(family = font$family))
  
  # fourth result plot: spread of intensities after normalization
  RES4 <- plot + ggplot2::geom_boxplot(
    ggplot2::aes(x=value,y=variable),
    color=cf(sampsize),
    alpha=0.4) + ggplot2::geom_vline(ggplot2::aes(xintercept=median(norm_melt$value))) +
    ggplot2::theme(legend.position="none",
                   axis.text=ggplot2::element_text(size=font$ax.num.size),
                   axis.title=ggplot2::element_text(size=font$ax.txt.size),
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.line = ggplot2::element_line(colour = 'black', size = .5),
                   text = ggplot2::element_text(family = font$family))
  
  # - - - - - - - - - - - - - - - - - - -
  
  list(tl=RES1, bl=RES2, tr=RES3, br=RES4)
}

#' @export
#' Plot a summary of metaboanalystR normalization results. Takes the total of 20 samples before and after normalization and plots that distribution.
#' @param mSet input user mSet
#' @param plot.theme function for ggplot theme used.
#' @return list of four plots that fit in a 2x2 raster used in MetaboShiny.
ggplotSampleNormSummary <- function(mSet,
                                    plot.theme,
                                    font,
                                    cf){
  # 4 by 4 plot, based on random 20-30 picked
  orig_data <- as.data.frame(mSet$dataSet$prenorm)
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
    plot.theme(base_size = 15) + ggplot2::geom_density(ggplot2::aes(x=value), colour="blue", fill="blue", alpha=0.4)+
    ggplot2::theme(legend.position="none",
                   axis.text=ggplot2::element_text(size=font$ax.num.size),
                   axis.title=ggplot2::element_text(size=font$ax.txt.size),
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.line = ggplot2::element_line(colour = 'black', size = .5),
                   text = ggplot2::element_text(family = font$family))
  
  RES2 <- ggplot2::ggplot(data=orig_melt) +
    plot.theme(base_size = 15) + ggplot2::geom_boxplot(
      ggplot2::aes(x=value,y=Label),
      color=cf(sampsize),
      alpha=0.4) + ggplot2::geom_vline(ggplot2::aes(xintercept=median(orig_melt$value),text=Label))+
    ggplot2::theme(legend.position="none",
                   axis.text=ggplot2::element_text(size=font$ax.num.size),
                   axis.title=ggplot2::element_text(size=font$ax.txt.size),
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.line = ggplot2::element_line(colour = 'black', size = .5),
                   text = ggplot2::element_text(family = font$family))
  
  RES3 <- ggplot2::ggplot(data=norm_melt_sums) +
    plot.theme(base_size = 15) + ggplot2::geom_density(ggplot2::aes(x=value), colour="pink", fill="pink", alpha=0.4)+
    ggplot2::theme(legend.position="none",
                   axis.text=ggplot2::element_text(size=font$ax.num.size),
                   axis.title=ggplot2::element_text(size=font$ax.txt.size),
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.line = ggplot2::element_line(colour = 'black', size = .5),
                   text = ggplot2::element_text(family = font$family))
  
  RES4 <- ggplot2::ggplot(data=norm_melt) +
    plot.theme(base_size = 15) + ggplot2::geom_boxplot(
      ggplot2::aes(x=value,y=Label),
      color=cf(sampsize),
      alpha=0.4) + ggplot2::geom_vline(ggplot2::aes(xintercept=median(norm_melt$value),text=Label))+
    ggplot2::theme(legend.position="none",
                   axis.text=ggplot2::element_text(size=font$ax.num.size),
                   axis.title=ggplot2::element_text(size=font$ax.txt.size),
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.line = ggplot2::element_line(colour = 'black', size = .5),
                   text = ggplot2::element_text(family = font$family))
  
  list(tl=RES1, bl=RES2, tr=RES3, br=RES4)
}


#' @export
ggplotMeba <- function(mSet, cpd, draw.average=T, cols,
                       cf,
                       plot.theme,
                       plotlyfy=TRUE,font){
  
  time.mode = mSet$dataSet$exp.type
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
  
  profile <- MetaboShiny::getProfile(mSet, 
                                     cpd, 
                                     mode="multi")
  profile$Individual <- mSet$dataSet$covars[match(mSet$dataSet$covars$sample, 
                                                  table = profile$Sample),"individual"][[1]]
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
    plot.theme(base_size = 15) +
    ggplot2::scale_color_manual(values=cols) +
    ggplot2::theme(legend.position=switch(time.mode, t="none", t1f="right"),
                   axis.text=ggplot2::element_text(size=font$ax.num.size),
                   axis.title=ggplot2::element_text(size=font$ax.txt.size),
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.line = ggplot2::element_line(colour = 'black', size = .5),
                   text = ggplot2::element_text(family = font$family)) +
    ggtitle(paste(cpd, "m/z")) + 
    xlab(Hmisc::capitalize(mSet$dataSet$facB.lbl)) + 
    labs(color = Hmisc::capitalize(switch(mSet$dataSet$exp.type, 
                                          t1f=mSet$dataSet$facA.lbl,
                                          t="Individual")))
  if(draw.average){
    p <- p + ggplot2::stat_summary(fun.y="mean", size=2, 
                                   geom="line", ggplot2::aes(x=GroupB, 
                                                             y=Abundance, 
                                                             color = Color, 
                                                             group = switch(time.mode, 
                                                                            t=c(1),
                                                                            t1f=Color)))
  }
  
  if(plotlyfy){
    plotly::ggplotly(p, tooltip="Individual", originalData=T)
  }else{
    p
  }
}

blackwhite.colors <- function(n){
  gray.colors(n, start=0, end=1)
}

#' @export
ggplotSummary <- function(mSet, cpd, shape.fac = "label", cols = c("black", "pink"),
                          cf = rainbow, plot.theme,
                          mode = "nm", plotlyfy = TRUE,
                          styles=c("box", "beeswarm"), add_stats = "mean",
                          col.fac = "label", txt.fac = "label",font){
  
  sourceTable = mSet$dataSet$norm
  
  if(length(styles) == 0){
    styles = c("beeswarm")
  }
  # - - -
  
  if(mSet$dataSet$exp.type %in% c("t1f", "2f")){
    mode = "multi"
  }
  
  profile <- MetaboShiny::getProfile(mSet, cpd, mode=if(mode == "nm") "stat" else "multi")

  df_line <- data.table::data.table(x = c(1,2),
                                    y = rep(min(profile$Abundance - 0.1), 2))
  stars = ""
  
  try({
    pval <- if(mode == "nm"){
      mSet$analSet$tt$sig.mat[which(rownames(mSet$analSet$tt$sig.mat) == cpd), "p.value"]
    }else{
      int.col <- grep("adj|Adj", colnames(mSet$analSet$aov2$sig.mat),value=T)
      int.col <- grep("int|Int", int.col, value=T)
      mSet$analSet$aov2$sig.mat[which(rownames(mSet$analSet$aov2$sig.mat) == cpd), int.col]
    }
    stars <- p2stars(pval)
  })
  
  profile$Shape <- if(shape.fac == "label"){
    as.factor(mSet$dataSet$cls)
  }else{
    as.factor(mSet$dataSet$covars[,..shape.fac][[1]])
  }
  
  profile$Color <- if(col.fac == "label"){
    as.factor(mSet$dataSet$cls)
  }else{
    as.factor(mSet$dataSet$covars[,..col.fac][[1]])
  }
  
  profile$Text <- if(txt.fac == "label"){
    as.factor(mSet$dataSet$cls)
  }else{
    as.factor(mSet$dataSet$covars[,..txt.fac][[1]])
  }
  
  p <- ggplot2::ggplot() + plot.theme(base_size = 15)
  
  profiles <- switch(mode,
                     multi = split(profile, f = profile$GroupA),
                     nm = list(x = data.table::data.table(profile)))
  
  if(length(cols) < length(levels(profile$Color))){
    cols <- cf(length(levels(profile$Color)))
  }
  
  i = 1
  
  for(prof in profiles){
    
    for(style in styles){
      
      switch(mode,
             nm = {
               p <- switch(style,
                           box = p + ggplot2::geom_boxplot(data = prof, alpha=0.4, aes(x = Group,
                                                                                       y = Abundance,
                                                                                       text = Text,
                                                                                       color = Group,
                                                                                       fill = Group)),
                           violin = p + ggplot2::geom_violin(data = prof, alpha=0.4, position = "identity", aes(x = Group,
                                                                                                                y = Abundance,
                                                                                                                color = Group,
                                                                                                                fill = Group)),
                           beeswarm = p + ggbeeswarm::geom_beeswarm(data = prof, alpha=0.7, size = 2, position = position_dodge(width=.3), aes(x = Group,
                                                                                                                                               y = Abundance,
                                                                                                                                               text= Text,
                                                                                                                                               color = Group,
                                                                                                                                               fill = Color)),
                           scatter = p + ggplot2::geom_point(data = prof, alpha=0.7, size = 2, aes(x = Group,
                                                                                                   y = Abundance,
                                                                                                   text=Text,
                                                                                                   color = Group,
                                                                                                   fill = Color), position = position_jitterdodge())
               )
             },
             multi = {
               p <- switch(style,
                           box = p + ggplot2::geom_boxplot(data = prof, alpha=0.4, aes(x = GroupB,
                                                                                       y = Abundance,
                                                                                       text = Text,
                                                                                       color = GroupA,
                                                                                       fill = GroupA)),
                           violin = p + ggplot2::geom_violin(data = prof, alpha=0.4, position = "identity", aes(x = GroupB,
                                                                                                                y = Abundance,
                                                                                                                group = GroupB,
                                                                                                                text = Text,
                                                                                                                color = GroupA,
                                                                                                                fill = GroupA)),
                           beeswarm = p + ggbeeswarm::geom_beeswarm(data = prof, alpha=0.7, size = 2, position = position_dodge(width=.3), aes(x = GroupB,
                                                                                                                                               y = Abundance,
                                                                                                                                               text=Text,
                                                                                                                                               color = GroupA,
                                                                                                                                               fill = GroupA)),
                           scatter = p + ggplot2::geom_point(data = prof, alpha=0.7, size = 2, position = position_jitterdodge(), aes(x = GroupB,
                                                                                                                                      y = Abundance,
                                                                                                                                      text=Text,
                                                                                                                                      color = GroupA,
                                                                                                                                      fill = GroupA))
               )
             })
    }
    
    p <- p +ggplot2::theme(legend.position="none",
                           plot.title = ggplot2::element_text(hjust = 0.5),
                           axis.text=ggplot2::element_text(size=font$ax.num.size),
                           axis.title=ggplot2::element_text(size=font$ax.txt.size),
                           axis.line = ggplot2::element_line(colour = 'black', size = .5),
                           text = ggplot2::element_text(family = font$family))+
      ggplot2::annotate("text",
                        x = switch(mode, nm = 1.5,
                                   multi = max(as.numeric(as.factor(profile$GroupB)))/2 + .5),
                        y = min(profile$Abundance - 0.3),
                        label = stars, size = 8, col = "black") + ggtitle(paste(cpd, "m/z"))
    
    if(!("box" %in% styles)){
      p <- switch(add_stats,
                  median = {
                    p + ggplot2::stat_summary(data = prof,
                                              aes( x = if(mode == "nm") Group else GroupB,
                                                   y = Abundance,
                                                   color = if(mode == "nm") Group else GroupA),
                                              fun.y = median,
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
                                              aes(x = if(mode == "nm") Group else GroupB,
                                                  y = Abundance),
                                              fun.y = mean,
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
    p <- p + ggplot2::scale_color_manual(values=cols[1:length(unique(profile$GroupA))])
    p <- p + ggplot2::scale_fill_manual(values = cols[1:length(unique(profile$GroupA))])
    if(mSet$dataSet$exp.type == "t"){
      p <- p + ggplot2::xlab("Time")
    }else{
      p <- p + ggplot2::xlab(Hmisc::capitalize(mSet$dataSet$facB.lbl))
    }
  }else{
    p <- p + ggplot2::scale_color_manual(values=cols[1:length(unique(profile$Color))])
    if(all(as.character(profile$Color) == as.character(profile$Group))){
      p <- p + ggplot2::scale_fill_manual(values = cols[1:length(unique(profile$Group))])
    }else{
      ncols = length(levels(profile$Color))
      scale = c(cols[1:length(unique(levels(profile$Group)))],cf(ncols))
      names(scale) <- c(levels(profile$Group),levels(profile$Color))
      p <- p + ggplot2::scale_fill_manual(values = scale)
    }  
    p <- p + ggplot2::xlab(Hmisc::capitalize(gsub(x=mSet$dataSet$cls.name, pattern = ":.*$", replacement="")))
  }
  # ---------------
  if(plotlyfy){
    plotly::ggplotly(p, tooltip = "Text")#, originalData=T)
  }else{
    p
  }
}

ggPlotAOV <- function(mSet, cf, n=20,
                      plot.theme, plotlyfy=TRUE,font){
  
  which_aov = if(mSet$dataSet$exp.type %in% c("t", "2f", "t1f")) "aov2" else "aov"
  
  profile <- data.table::as.data.table(mSet$analSet[[which_aov]]$p.log[mSet$analSet[[which_aov]]$inx.imp],
                                       keep.rownames = T)
  
  if(nrow(profile)==0){
    shiny::showNotification("No significant hits")
    return(NULL)
  }
  
  colnames(profile) <- c("cpd", "-logp")
  profile[,2] <- round(profile[,2], digits = 2)
  profile$Peak <- c(1:nrow(profile))
  colnames(profile)[1:2] <- c("m/z", "-log(p)")
  scaleFUN <- function(x) sprintf("%.0f", x)
  
  xaxis = seq(0,600, 50)
  # ---------------------------
  p <- ggplot2::ggplot(data=profile) +
    ggplot2::geom_point(ggplot2::aes(x=`m/z`, y=`-log(p)`,
                                     text=`m/z`, color=`-log(p)`, 
                                     key=`m/z`)) +
    plot.theme(base_size = 15) +
    ggplot2::scale_x_discrete(breaks = xaxis, labels=as.character(xaxis)) + 
    ggplot2::theme(legend.position="none",
                   axis.text=ggplot2::element_text(size=font$ax.num.size),
                   axis.title=ggplot2::element_text(size=font$ax.txt.size),
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.line = ggplot2::element_line(colour = 'black', size = .5),
                   text = ggplot2::element_text(family = font$family))+
    ggplot2::scale_colour_gradientn(colours = cf(n)) +
    #ggplot2::scale_y_log10()+
    ggplot2::scale_y_continuous(labels=scaleFUN)
  if(plotlyfy){
    plotly::ggplotly(p, tooltip="m/z")
  }else{
    p
  }
}
ggPlotTT <- function(mSet, cf, n=20,
                     plot.theme, plotlyfy=TRUE,font){
  profile <- data.table::as.data.table(mSet$analSet$tt$p.log[mSet$analSet$tt$inx.imp],keep.rownames = T)
  
  if(nrow(profile)==0){
    shiny::showNotification("No significant hits")
    return(NULL)
  }
  
  colnames(profile) <- c("cpd", "-logp")
  profile[,2] <- round(profile[,2], digits = 2)
  profile$Peak <- c(1:nrow(profile))
  colnames(profile)[1:2] <- c("m/z", "-log(p)")
  scaleFUN <- function(x) sprintf("%.0f", x)
  
  xaxis = seq(0,600, 50)
  # ---------------------------
  p <- ggplot2::ggplot(data=profile) +
    ggplot2::geom_point(ggplot2::aes(x=`m/z`, y=`-log(p)`,
                                     text=`m/z`, color=`-log(p)`, 
                                     key=`m/z`)) +
    plot.theme(base_size = 15) +
    ggplot2::scale_x_discrete(breaks = xaxis, labels=as.character(xaxis)) + 
    ggplot2::theme(legend.position="none",
                   axis.text=ggplot2::element_text(size=font$ax.num.size),
                   axis.title=ggplot2::element_text(size=font$ax.txt.size),
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.line = ggplot2::element_line(colour = 'black', size = .5),
                   text = ggplot2::element_text(family = font$family))+
    ggplot2::scale_colour_gradientn(colours = cf(n)) +
    #ggplot2::scale_y_log10()+
    ggplot2::scale_y_continuous(labels=scaleFUN)
  if(plotlyfy){
    plotly::ggplotly(p, tooltip="m/z")
  }else{
    p
  }
}

ggPlotPattern <- function(mSet, cf, n=20,
                          plot.theme,
                          plotlyfy=TRUE,font){
  profile <- data.table::as.data.table(mSet$analSet$corr$cor.mat,keep.rownames = T)
  profile <- profile[1:n]
  
  if(nrow(profile)==0){
    shiny::showNotification("No significant hits")
    return(NULL)
  }
  
  colnames(profile)[1] <- c("cpd")
  #profile$Peak <- c(1:nrow(profile))
  scaleFUN <- function(x) sprintf("%.2f", x)
  
  # ---------------------------
  p <- ggplot2::ggplot(data=profile) +
    ggplot2::geom_bar(mapping = aes(x = reorder(cpd, -`p-value`), 
                                    y = correlation, 
                                    color = `p-value`, 
                                    text = reorder(cpd, -`p-value`),
                                    fill = `p-value`), 
                      stat = "identity", alpha=0.5) + 
    ggplot2::ggtitle(paste("Associated with pattern", 
                           mSet$analSet$corr$pattern)) +
    ggplot2::coord_flip() +
    ggplot2::ylab("correlation") + 
    ggplot2::xlab("m/z") + 
    ggplot2::labs(fill="p-value", 
                  color="p-value")+
    plot.theme(base_size = 15) +
    ggplot2::theme(legend.position="none",
                   axis.text=ggplot2::element_text(size=font$ax.num.size),
                   axis.title=ggplot2::element_text(size=font$ax.txt.size),
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5),
                   axis.line = ggplot2::element_line(colour = 'black', size = .5),
                   text = ggplot2::element_text(family = font$family)) +
    ggplot2::scale_colour_gradientn(colours = cf(n)) +
    ggplot2::scale_fill_gradientn(colours = cf(n)) +
    ggplot2::scale_y_continuous(labels=scaleFUN)
  
  if(plotlyfy){
    plotly::ggplotly(p, tooltip="text")
  }else{
    p
  }
}

ggPlotFC <- function(mSet, cf, n=20,
                     plot.theme,
                     plotlyfy=TRUE,font){
  profile <- data.table::as.data.table(mSet$analSet$fc$fc.log[mSet$analSet$fc$inx.imp],keep.rownames = T)
  
  if(nrow(profile)==0){
    shiny::showNotification("No significant hits")
    return(NULL)
  }
  
  colnames(profile) <- c("cpd", "log2fc")
  profile$Peak <- c(1:nrow(profile))
  # ---------------------------
  p <- ggplot2::ggplot(data=profile) +
    ggplot2::geom_point(ggplot2::aes(x=Peak, y=log2fc, text=log2fc, color=log2fc, key=cpd)) +
    ggplot2::geom_abline(ggplot2::aes(intercept = 0, slope = 0)) +
    plot.theme(base_size = 15) +
    ggplot2::theme(legend.position="none",
                   axis.text=ggplot2::element_text(size=font$ax.num.size),
                   axis.title=ggplot2::element_text(size=font$ax.txt.size),
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.line = ggplot2::element_line(colour = 'black', size = .5),
                   text = ggplot2::element_text(family = font$family))+
    ggplot2::scale_colour_gradientn(colours = cf(n))
  
  if(plotlyfy){
    plotly::ggplotly(p, tooltip="log2fc")
  }else{
    p
  }
}

ggPlotVolc <- function(mSet,
                       cf,
                       n=20,
                       plot.theme,
                       plotlyfy=TRUE,
                       font ){
  
  vcn<-mSet$analSet$volcano;
  
  if(nrow(vcn$sig.mat)==0){
    shiny::showNotification("No significant hits")
    return(NULL)
  }
  
  dt <- data.table::as.data.table(vcn$sig.mat[,c(2,4)],keep.rownames = T)
  colnames(dt) <- c("cpd", "log2FC", "-log10P")
  p <- ggplot2::ggplot() +
    #ggplot2::geom_point(data=dt[!imp.inx], ggplot2::aes(x=log2FC, y=minlog10P)) +
    ggplot2::geom_point(data=dt, ggplot2::aes(x=log2FC,
                                              y=`-log10P`,
                                              text=cpd,
                                              color=abs(log2FC*`-log10P`),
                                              key=cpd)) +
    plot.theme(base_size = 15) +
    ggplot2::theme(legend.position="none",
                   axis.text=ggplot2::element_text(size=font$ax.num.size),
                   axis.title=ggplot2::element_text(size=font$ax.txt.size),
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.line = ggplot2::element_line(colour = 'black', size = .5),
                   text = ggplot2::element_text(family = font$family))+
    ggplot2::scale_colour_gradientn(colours = cf(n),guide=FALSE)
  
  if(plotlyfy){
    plotly::ggplotly(p, tooltop="cpd")
  }else{
    p
  }
}

ggPlotClass <- function(mSet,
                        pls.type = "plsda",
                        cf,
                        pcs = 3,
                        plot.theme,
                        plotlyfy=TRUE,font){
  res <- mSet$analSet$plsda$fit.info
  colnames(res) <- 1:ncol(res)
  # best.num <- mSet$analSet$plsda$best.num
  # choice <- mSet$analSet$plsda$choice
  df <- melt(res)
  df$Component <- paste0("PC",df$Component)
  colnames(df) <- c("Metric", "Component", "Value")
  p <- ggplot(df, aes(x=Metric, y=Value, fill=Metric)) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::theme_minimal() +
    plot.theme(base_size = 15) +
    ggplot2::theme(legend.position="none",
                   axis.text=ggplot2::element_text(size=font$ax.num.size),
                   axis.title=ggplot2::element_text(size=font$ax.txt.size),
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.line = ggplot2::element_line(colour = 'black', size = .5),
                   text = ggplot2::element_text(family = font$family))+
    facet_grid(~Component) +
    ggplot2::scale_fill_manual(values=cf(pcs))
  if(plotlyfy){
    plotly::ggplotly(p)
  }else{
    p
  }
}

ggPlotPerm <- function(mSet,
                       pls.type = "plsda",
                       cf,
                       pcs = 3,
                       plot.theme,
                       plotlyfy=TRUE,font){
  bw.vec <- mSet$analSet$plsda$permut
  len <- length(bw.vec)
  df <- melt(bw.vec)
  colnames(df) = "acc"
  # round p value
  pval <- mSet$analSet$plsda$permut.p
  rounded <- round(as.numeric(stringr::str_match(pval, "0\\.\\d*")), digits = 3)
  pval <- gsub(pval, pattern = "(0\\.\\d*)", replacement=rounded)
  # - - -
  p <- ggplot(df) +
    ggplot2::geom_histogram(mapping=aes(x=acc, y=..count.., fill=factor(..count..)),
                            binwidth=0.01) +
    ggplot2::scale_fill_manual(values=cf(20)) +
    plot.theme(base_size = 15) +
    ggplot2::theme(legend.position="none",
                   axis.text=ggplot2::element_text(size=font$ax.num.size),
                   axis.title=ggplot2::element_text(size=font$ax.txt.size),
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.line = ggplot2::element_line(colour = 'black', size = .5),
                   text = ggplot2::element_text(family = font$family))+
    labs(x="Accuracy", y = "Permutations") +
    ggplot2::geom_segment(data=df,
                          color="black",
                          x=bw.vec[1],
                          xend=bw.vec[1],
                          y=0,
                          aes(yend=.1*nrow(df)),
                          size=1.5,
                          linetype=8) +
    ggplot2::geom_text(mapping = aes(x = bw.vec[1], y =  .11*nrow(df), label = pval), color = "black", size = 4)
  
  if(plotlyfy){
    plotly::ggplotly(p)
  }else{
    p
  }
}

ggPlotROC <- function(data,
                      attempts = 50,
                      cf,
                      plot.theme,
                      plotlyfy=TRUE,font,
                      class_type="b"){
  
  mean.auc <- data$m_auc
  perf.long <- data$perf
  
  means.per.comp=perf.long[, lapply(.SD, mean), by = .(comparison)]
  
  if(length(unique(means.per.comp$comparison))>2){
    class_type = "m"
    shiny::showNotification("Calculating AUCs per comparison...")
    perf.long$comparison <- pbapply::pbsapply(perf.long$comparison,
                                              function(comp){
                                                paste0(comp, " || avg. AUC=", round(means.per.comp[comparison == comp]$AUC, digits=3), " ||")
                                              })  
  }else{
    class_type = "b"
  }
  
  cols = cf(attempts)
  
  p <- ggplot(perf.long, aes(FPR,TPR)) +
    ggplot2::geom_path(alpha=.5,
                       cex=.5,
                       aes(color = comparison, group = attempt)) +
    ggplot2::annotate("text",
                      label = paste0("Average AUC: ",
                                     format(mean.auc,
                                            2,
                                            drop0trailing = TRUE,
                                            digits = 2)),
                      size = 8,
                      x = 0.77,
                      y = 0.03) +
    plot.theme(base_size = 10) +
    ggplot2::stat_summary_bin(#alpha=.6,
      aes(FPR, TPR, 
          group=comparison), 
      fun.y=mean, geom="line", 
      cex = 2.3,color="black")+
    ggplot2::stat_summary_bin(#alpha=.6,
      aes(FPR, TPR, 
          color=comparison, 
          group=comparison), 
      fun.y=mean, geom="line", 
      cex = 1.2) +
    ggplot2::stat_summary_bin(aes(FPR, TPR), 
                              fun.y=mean, color="black", 
                              geom="line", cex = 2) +
    
    #ggplot2::scale_color_gradientn(colors = cols) +
    ggplot2::theme(legend.position= if(class_type == "b") "none" else "right",
                   axis.text=ggplot2::element_text(size=font$ax.num.size),
                   axis.title=ggplot2::element_text(size=font$ax.txt.size),
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.line = ggplot2::element_line(colour = 'black', size = .5),
                   text = ggplot2::element_text(family = font$family))+
    #scale_y_continuous(limits=c(0,1)) +
    ggplot2::coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE) +
    ggplot2::coord_cartesian(xlim = c(.04,.96), ylim = c(.04,.96))
  
  if(plotlyfy){
    plotly::ggplotly(p)
  }else{
    p
  }
}

ggPlotBar <- function(data,
                      attempts=50,
                      cf,
                      topn=50,
                      plot.theme,
                      ml_name,
                      ml_type,
                      plotlyfy=TRUE,
                      font){
  
  if(ml_name != ""){
    lname = ml_name
  }else{
    lname <- "all"
  }
  
  if(ml_type == "glmnet"){
    colnames(data) = c("mz", "importance.mean", "dummy")
    data.ordered <- data[order(data$importance, decreasing=T),1:2]
  }else{
    data.norep <- data[,-3]
    data.ci = Rmisc::group.CI(importance ~ mz, data.norep)
    
    data.ordered <- data.ci[order(data.ci$importance.mean, decreasing = T),]
  }
  
  data.subset <- data.ordered[1:topn,]    
  
  p <- ggplot(data.subset, aes(x = reorder(mz,-importance.mean),
                               y = importance.mean,
                               label = mz)) +
    ggplot2::geom_bar(stat = "identity",
                      aes(fill = importance.mean)) +
    ggplot2::scale_fill_gradientn(colors=cf(20)) +
    plot.theme() +
    ggplot2::theme(legend.position="none",
                   axis.text=ggplot2::element_text(size=font$ax.num.size),
                   axis.title=ggplot2::element_text(size=font$ax.txt.size),
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.line = ggplot2::element_line(colour = 'black', size = .5),
                   axis.text.x=ggplot2::element_blank(),
                   axis.ticks.x=ggplot2::element_blank(),
                   text = ggplot2::element_text(family = font$family))+
    labs(x="Top hits (m/z)",y=if(ml_type == "glmnet") "Times included in final model" else "Relative importance (%)")
  
  if(topn <= 15){
    p <- p + ggplot2::geom_text(aes(x=mz, y=importance.mean, label=sapply(mz, function(x){
      if(is.na(as.numeric(as.character(x)))){
        if(grepl(x, pattern = "_")){
          as.character(gsub(x, pattern = "_", replacement = "\n"))
        }else{
          as.character(x)
        }
      }else{
        round(as.numeric(as.character(x)),digits=1)
      }
    })), size = 4, vjust = -.5, lineheight = .6)
    plotlyfy=F
  }
  
  mzdata <- p$data
  mzdata$mz <- gsub(mzdata$mz, pattern = "`|'", replacement="")
  
  if(plotlyfy){
    list(mzdata = mzdata, plot = plotly::ggplotly(p, tooltip="label"))
  }else{
    list(mzdata = mzdata, plot = p)
  }
}

plotPCA.3d <- function(mSet,
                       cols ,
                       shape.fac="label",
                       pcx, pcy, pcz,
                       mode="pca",font,
                       cf){
  
  switch(mode,
         pca = {
           df <- mSet$analSet$pca$x
           x.var <- round(mSet$analSet$pca$variance[pcx] * 100.00, digits=1)
           y.var <- round(mSet$analSet$pca$variance[pcy] * 100.00, digits=1)
           z.var <- round(mSet$analSet$pca$variance[pcz] * 100.00, digits=1)
         }, ipca = {
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
           
           x.var <- plsda.table[PC == pcx]$var
           y.var <- plsda.table[PC == pcy]$var
           z.var <- plsda.table[PC == pcz]$var
           
           # --- coordinates ---
           df <- mSet$analSet$plsr$scores
           class(df) <- "matrix"
           colnames(df) <- paste0("PC", 1:ncol(df))
         })
  
  if(mode == "ipca"){
    fac.lvls <- length(levels(mSet$dataSet$exp.fac))
    classes = mSet$dataSet$exp.fac
  }else{
    fac.lvls <- length(levels(mSet$dataSet$cls))
    classes = mSet$dataSet$cls
  }
  
  df <- as.data.frame(df)
  rownames(df) <- rownames(mSet$dataSet$norm)
  df_list <- list(df)
  
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
  
  # --- add ellipses ---
  
  symbols = c('circle',
              'diamond',
              'square',
              'x',
              'o')
  
  symbol.vec<-if(is.null(shape.fac)){
    rep('circle', times = length(classes))
  }else if(shape.fac == "label"){
    rep('circle', times = length(classes))
  }else{
    as.factor(mSet$dataSet$covars[,..shape.fac][[1]])
  }
  
  plots_facet <- lapply(1:length(df_list), function(i){
    
    df = df_list[[i]]
    
    orig_idx = which(rownames(df) %in% rownames(mSet$dataSet$norm))
    
    plots <- plotly::plot_ly(showlegend = F)
    
    show.orbs <- c()
    
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
        plots = plots %>% add_trace(
          x=mesh$x,
          y=mesh$y,
          z=mesh$z,
          type='mesh3d',
          alphahull=0,
          opacity=0.1,
          hoverinfo="none"
        )
      }
      
      show.orbs <- c(show.orbs, worked)
    }
    
    adj_plot <- plotly_build(plots)
    rgbcols <- toRGB(cols[show.orbs])
    c = 1
    
    for(i in seq_along(adj_plot$x$data)){
      item = adj_plot$x$data[[i]]
      if(item$type == "mesh3d"){
        adj_plot$x$data[[i]]$color <- rgbcols[c]
        adj_plot$x$data[[i]]$visible <- TRUE
        c = c + 1
      }
    }
    t <- list(
      family = font$family
    )
    
    # --- return ---
    pca_plot <- adj_plot %>% add_trace(
      hoverinfo = 'text',
      text = rownames(df),
      x = df[,pcx],
      y = df[,pcy],
      z = df[,pcz],
      visible = rep(T, times=fac.lvls),
      type = "scatter3d",
      opacity=1,
      color = classes[orig_idx],
      colors=cols,
      symbol = symbol.vec[orig_idx],
      symbols = c('circle',
                  'diamond',
                  'square',
                  'x',
                  'o')
    ) %>%  layout(font = t,
                  scene = list(
                    aspectmode="cube",
                    xaxis = list(
                      titlefont = list(size = font$ax.txt.size * 1.5),
                      title = gsubfn::fn$paste("$pcx ($x.var %)")),
                    yaxis = list(
                      titlefont = list(size = font$ax.txt.size * 1.5),
                      title = gsubfn::fn$paste("$pcy ($y.var %)")),
                    zaxis = list(
                      titlefont = list(size = font$ax.txt.size * 1.5),
                      title = gsubfn::fn$paste("$pcz ($z.var %)"))))
    # --- return ---
    pca_plot
  })
  subplot(plots_facet,shareX = F, shareY = F)
}

plotPCA.2d <- function(mSet, shape.fac = "label", cols,
                       pcx, pcy, mode="pca", plot.theme,
                       plotlyfy="T", font, cf){
  
  classes <- switch(mode,
                    ipca = mSet$dataSet$facA,
                    pca = mSet$dataSet$cls,
                    plsda = mSet$dataSet$cls)
  
  symbols = c("16",#'circle',
              "18",#'diamond',
              "15",#'square',
              "4",#'x',
              "1"#'o'
  )
  
  
  switch(mode,
         pca={
           df <- mSet$analSet$pca$x
           #x =1;y=2;z=3
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
         ipca = {
           df <- mSet$analSet$pca$x
           #x =1;y=2;z=3
           x.var <- round(mSet$analSet$pca$variance[pcx] * 100.00, digits=1)
           y.var <- round(mSet$analSet$pca$variance[pcy] * 100.00, digits=1)
           fac.lvls <- length(levels(mSet$dataSet$facA))
           
           xc=mSet$analSet$pca$x[, pcx]
           yc=mSet$analSet$pca$x[, pcy]
           
           dat_long <- data.table(variable = names(xc),
                                  group = classes,
                                  groupB = mSet$dataSet$facB,
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
           
           x.var <- plsda.table[PC == pcx]$var
           y.var <- plsda.table[PC == pcy]$var
           
           # --- coordinates ---
           df <- mSet$analSet$plsr$scores
           colnames(df) <- paste0("PC", 1:ncol(df))
           rownames(df) <- rownames(mSet$dataSet$norm)
           
           xc=df[, pcx]
           yc=df[, pcy]
           
           dat_long <- data.table(variable = names(xc),
                                  group = classes,
                                  x = xc,
                                  y = yc)
         })
  # - - - - - - - - -
  
  dat_long$shape <- if(is.null(shape.fac)){
    factor(1) # all same shape...
  } else if(shape.fac == "label"){
    dat_long$group
  }else{
    as.factor(mSet$dataSet$covars[,..shape.fac][[1]])
  }
  
  cols <- if(is.null(cols)) cf(length(levels(classes))) else{
    if(length(cols) < length(levels(classes))){
      cols <- cf(levels(classes))
    }
    cols
  }
  
  p <- ggplot(dat_long, aes(x, y,group=group)) +
    ggplot2::geom_point(size=5, aes(shape=shape,
                                    text=variable,
                                    fill=group,
                                    color=group), alpha=0.7)+
    ggplot2::stat_ellipse(geom = "polygon", aes(fill=group), alpha = 0.3,level = .95) +
    plot.theme() +
    ggplot2::theme(legend.position="none",
                   axis.text=ggplot2::element_text(size=font$ax.num.size),
                   axis.title=ggplot2::element_text(size=font$ax.txt.size),
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.line = ggplot2::element_line(colour = 'black', size = .5),
                   text = ggplot2::element_text(family = font$family))+
    ggplot2::scale_x_continuous(name=gsubfn::fn$paste("$pcx ($x.var%)")) +
    ggplot2::scale_y_continuous(name=gsubfn::fn$paste("$pcy ($y.var%)")) +
    ggplot2::scale_fill_manual(values=cols) +
    ggplot2::scale_color_manual(values = cols)
  
  if(mode == "ipca"){
    p <- p + facet_wrap(~groupB,)
    p <- p + ggplot2::ggtitle(Hmisc::capitalize(mSet$dataSet$facB.lbl))
  }
  
  if(plotlyfy){
    plotly::ggplotly(p)
  }else{
    p
  }
  
}

ggPlotVenn <- function(mSet,
                       venn_yes,
                       top = 100,
                       cols,
                       cf,
                       plotlyfy=TRUE,font){
  
  experiments <- stringr::str_match(unlist(venn_yes$now), 
                                    pattern = "\\(.*\\)")[,1]
  
  experiments <- unique(gsub(experiments, pattern = "\\(\\s*(.+)\\s*\\)", replacement="\\1"))
  
  table_list <- lapply(experiments, function(experiment){
    
    analysis = mSet$storage[[experiment]]$analysis
    
    rgx_exp <- gsub(experiment, pattern = "\\(", replacement = "\\\\(")
    rgx_exp <- gsub(rgx_exp, pattern = "\\)", replacement = "\\\\)")
    rgx_exp <- gsub(rgx_exp, pattern = "\\-", replacement = "\\\\-")
    rgx_exp <- gsub(rgx_exp, pattern = "\\+", replacement = "\\\\+")
    
    categories = grep(unlist(venn_yes$now),
                      pattern = paste0("\\(",rgx_exp, "\\)"), value = T)
    
    categories = gsub(categories, pattern = " \\(\\s*(.+)\\s*\\)", replacement = "")
    
    # go through the to include analyses
    
    tables <- lapply(categories, function(name){
      
      base_name <- search_name <- gsub(name, pattern = " -.*$| ", replacement="")
      
      if(base_name %in% gbl$constants$ml.models){
        search_name <- "ml"
      }
      
      # fetch involved mz values
      tbls <- switch(search_name,
                     ml = {
                       which.ml <- gsub(name, pattern = "^.*- | ", replacement="")
                       mzvals = analysis$ml[[base_name]][[which.ml]]$bar[order(analysis$ml[[base_name]][[which.ml]]$bar$importance,
                                                                               decreasing = T),]$mz
                       mzvals <- type.convert(gsub(mzvals, pattern = "'|`", replacement=""))
                       res <- list(mzvals)
                       names(res) <- paste0(which.ml, " (", base_name, ")")
                       # - - -
                       res
                     },
                     aov = {
                       res = list(as.numeric(rownames(analysis$aov$sig.mat[order(analysis$aov$sig.mat[,2],
                                                                                 decreasing = F),])))
                       names(res) = base_name
                       res
                     },
                     aov2 = {
                       res = list(as.numeric(rownames(analysis$aov2$sig.mat[order(analysis$aov2$sig.mat[,"Interaction(adj.p)"],
                                                                                  decreasing = F),])))
                       names(res) = base_name
                       res
                     },
                     asca = {
                       res = list(as.numeric(rownames(analysis$asca$sig.list$Model.ab[order(analysis$asca$sig.list$Model.ab[,1],
                                                                                            decreasing = T),])))
                       names(res) = base_name
                       res
                     },
                     MB = {
                       res = list(as.numeric(rownames(analysis$MB$stats))[order(analysis$MB$stats[,1],
                                                                                decreasing = T)])
                       names(res) = base_name
                       res
                     },
                     tt = {
                       res = list(as.numeric(rownames(analysis$tt$sig.mat[order(analysis$tt$sig.mat[,2],
                                                                                decreasing = F),])))
                       names(res) = base_name
                       res
                     },
                     fc = {
                       res = list(as.numeric(rownames(analysis$fc$sig.mat[order(abs(analysis$fc$sig.mat[,2]),
                                                                                decreasing = F),])))
                       names(res) = base_name
                       res
                     },
                     volcano = {
                       res = list(as.numeric(rownames(analysis$volcano$sig.mat)))
                       names(res) = base_name
                       res
                     },
                     plsda = {
                       
                       which.plsda <- gsub(name, pattern = "^.*- | ", replacement="")
                       
                       compounds_pc <- data.table::as.data.table(analysis$plsda$vip.mat,keep.rownames = T)
                       colnames(compounds_pc) <- c("rn", paste0("PC", 1:(ncol(compounds_pc)-1)))
                       ordered_pc <- setorderv(compounds_pc, which.plsda, -1)
                       
                       res <- list(ordered_pc$rn)
                       names(res) <- paste0(which.plsda, " (PLS-DA)")
                       # - - -
                       res
                     },
                     pca = {
                       which.pca <- gsub(name, pattern = "^.*- | ", replacement="")
                       
                       compounds_pc <- data.table::as.data.table(analysis$pca$rotation,keep.rownames = T)
                       ordered_pc <- setorderv(compounds_pc, which.pca, -1)
                       res <- list(ordered_pc$rn)
                       names(res) <- paste0(which.pca, " (PCA)")
                       # - - -
                       res
                     },
                     volc = {
                       res <- list(rownames(analysis$volcano$sig.mat))
                       names(res) = base_name
                       res
                     },
                     {MetaboShiny::metshiAlert("Not currently supported...")
                       return(NULL)})
      
      if(is.null(tbls)) return(NULL)
      
      # user specified top hits only
      tbls_top <- lapply(tbls, function(tbl){
        if(length(tbl) < top){
          tbl
        }else{
          tbl[1:top]
        }
      })
      names(tbls_top) <- paste0(experiment, ": ", names(tbls_top))
      tbls_top
    })
    
    # unnest the nested lists
    flattened <- flattenlist(tables)
    
    # remove NAs
    flattened <- lapply(flattened, function(x) x[!is.na(x)])
    
    #rename and remove regex-y names
    names(flattened) <- gsub(x = names(flattened), pattern = "(.*\\.)(.*$)", replacement = "\\2")
    
    # return
    flattened
  })
  
  flattened <- flattenlist(table_list)
  names(flattened) <- gsub(x = names(flattened), pattern = "(.*\\.)(.*$)", replacement = "\\2")
  flattened <- lapply(flattened, function(x) x[!is.na(x)])
  
  # how many circles need to be plotted? (# of included analysis)
  circles = length(flattened)
  
  # generate the initial plot - the POLYGONS
  venn.plot <- VennDiagram::venn.diagram(x = flattened,
                                         filename = NULL)
  
  # split the plots into its individual elements
  items <- strsplit(as.character(venn.plot), split = ",")[[1]]
  
  # get which are circles
  circ_values <- data.frame(
    id = 1:length(grep(items, pattern="polygon"))
    #,value = c(3, 3.1, 3.1, 3.2, 3.15, 3.5)
  )
  
  # get which are text
  txt_values <- data.frame(
    id = grep(items, pattern="text"),
    value = unlist(lapply(grep(items, pattern="text"), function(i) venn.plot[[i]]$label))
  )
  
  # TODO: figure out what i did here again...
  txt_values$value <- gsub(x = txt_values$value, pattern = "(.*\\.)(.*$)", replacement = "\\2")
  #categories <- c(categories, input$rf_choice, input$ls_choice, input$plsda_choice)
  
  # get x and y values for circles
  x_c = unlist(lapply(grep(items, pattern="polygon"), function(i) venn.plot[[i]]$x))
  y_c = unlist(lapply(grep(items, pattern="polygon"), function(i) venn.plot[[i]]$y))
  
  # get x and y values for text
  x_t = unlist(lapply(grep(items, pattern="text"), function(i) venn.plot[[i]]$x))
  y_t = unlist(lapply(grep(items, pattern="text"), function(i)venn.plot[[i]]$y))
  
  # table with positions and ids for circles
  positions_c <- data.frame(
    id = rep(circ_values$id, each = length(x_c)/length(circ_values$id)),
    x = x_c,
    y = y_c
  )
  
  # table with positions and ids for text
  positions_t <- data.frame(
    id = rep(txt_values$id, each = length(x_t)/length(txt_values$id)),
    x = x_t,
    y = y_t
  )
  
  # merge them together for use in ggplot
  datapoly <- merge(circ_values, positions_c, by=c("id"))
  datatxt <- merge(txt_values, positions_t, by=c("id"))
  
  # make sure only the wanted analyses are in there
  numbers <- datatxt[!(datatxt$value %in% names(flattened)),]
  headers <- datatxt[(datatxt$value %in% names(flattened)),]
  
  # move numbers slightly if there are only 2 circles
  if(circles == 2){
    occur <- table(numbers$y)
    newy <- names(occur[occur == max(occur)])
    # - - -
    numbers$y <- as.numeric(c(newy))
  }
  
  # generate plot with ggplot
  p <- ggplot(datapoly,
              aes(x = x,
                  y = y)) + ggplot2::geom_polygon(colour="black", alpha=0.5, aes(fill=id, group=id)) +
    ggplot2::geom_text(mapping = aes(x=x-.02, y=y, label=value), data = numbers, size = 4, hjust = 0) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position="none",
                   text=ggplot2::element_text(#size=font$ax.num.size,
                     family = font$family),
                   panel.grid = ggplot2::element_blank()) +
    ggplot2::scale_fill_gradientn(colours =
                                    cf(circles)) +
    ggplot2::coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE)
  
  # - - text with white outline - -
  p <- p + shadowtext::geom_shadowtext(mapping = aes(x=x, y=y, label=gsub(value,
                                                                          pattern = ":",
                                                                          replacement = "\n")),
                                       color = "black",
                                       bg.color = "white",
                                       data = headers,
                                       size = 5, hjust = 0.5,
                                       fontface = 2) +
    ggplot2::scale_x_continuous(expand = c(.1, .1)) +
    ggplot2::scale_y_continuous(expand = c(.1, .1))
  if(plotlyfy){
    list(plot = plotly::ggplotly(p), info = flattened)
  }else{
    list(plot = p, info = flattened)
  }
}

ggPlotScree <- function(mSet, plot.theme, cf, font, pcs=20){
  df <- data.table::data.table(
    pc = 1:length(names(mSet$analSet$pca$variance)),
    var = round(mSet$analSet$pca$variance*100,digits = 1))
  p <- ggplot2::ggplot(data=df[1:20,]) + ggplot2::geom_line(mapping = aes(x=pc, y=var, colour=var), cex=3) +
    plot.theme(base_size = 10) +
    ggplot2::scale_colour_gradientn(colours = cf(20)) +
    ggplot2::theme(legend.position="none",
                   axis.text=ggplot2::element_text(size=font$ax.num.size),
                   axis.title=ggplot2::element_text(size=font$ax.txt.size),
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.line = ggplot2::element_line(colour = 'black', size = .5),
                   text = ggplot2::element_text(family = font$family))
  # - - - - -
  p
}

ggPlotWordBar <- function(wcdata, plot.theme, cf, font, plotlyfy=T){
  g <- ggplot(wcdata, aes(y = freq, x = reorder(word, 
                                                freq, 
                                                sum)))
  g <- g + ggplot2::geom_bar(aes(fill = freq),
                             stat = "identity") +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_gradientn(colors=cf(256)) +
    plot.theme() +
    ggplot2::theme(legend.position="none",
                   axis.text=ggplot2::element_text(size=font$ax.num.size),
                   axis.title=ggplot2::element_text(size=font$ax.txt.size),
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.line = ggplot2::element_line(colour = 'black', size = .5),
                   axis.text.x=ggplot2::element_blank(),
                   axis.ticks.x=ggplot2::element_blank(),
                   text = ggplot2::element_text(family = font$family)) +
    labs(x="Word",y="Frequency")
  if(plotlyfy){
    plotly::ggplotly(g,tooltip = "freq")
  }else{
    g
  }
}

ggPlotPower <- function(mSet,
                        cf,
                        plot.theme,
                        plotlyfy = TRUE,
                        font,
                        comparisons,
                        max_samples){
  
  cols = cf(length(comparisons))
  
  data = data.table::rbindlist(lapply(comparisons, function(comp) data.table::data.table(samples = mSet$analSet$power[[comp]]$Jpred,
                                                                                         power = mSet$analSet$power[[comp]]$pwrD,
                                                                                         comparison = c(comp))))
  
  print(head(data))
  
  #data$comparison <- substr(gsub(data$comparison, pattern = " .*$", replacement = ""), 1, 10)
    
  if(ncol(data) == 1){
    stop("Something went wrong! Try other settings please :(")
  }else{
    p <- ggplot(data, aes(samples,power)) +
      ggplot2::geom_path(alpha=.5,
                         cex=.5,
                         aes(color = comparison, group = comparison)) +
      plot.theme(base_size = 10) +
      ggplot2::stat_summary_bin(#alpha=.6,
        aes(samples, 
            power, 
            group=comparison), 
        fun.y=mean, geom="line", 
        cex = 2.3,color="black") +
      ggplot2::stat_summary_bin(#alpha=.6,
        aes(samples, power, 
            color=comparison, 
            group=comparison), 
        fun.y=mean, geom="line", 
        cex = 1.2) +
      ggplot2::stat_summary_bin(aes(samples,
                                    power), 
                                fun.y=mean, color="black", 
                                geom="line", cex = 2) +
      ggplot2::theme(legend.position= "none",
                     axis.text=ggplot2::element_text(size=font$ax.num.size),
                     axis.title=ggplot2::element_text(size=font$ax.txt.size),
                     plot.title = ggplot2::element_text(hjust = 0.5),
                     axis.line = ggplot2::element_line(colour = 'black', size = .5),
                     text = ggplot2::element_text(family = font$family, size = 15)) +
      ggplot2::coord_cartesian(xlim = c(0,max_samples), ylim = c(.04,.96))
    
    if(plotlyfy){
      plotly::ggplotly(p, tooltip="comparison")
    }else{
      p
    } 
  }
}