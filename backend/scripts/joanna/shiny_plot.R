#' @export
#' Plot a summary of metaboanalystR normalization results. Takes the total of 20 m/z values and 20 samples before and after normalization and plots that distribution.
#' @param mSet input user mSet
#' @param colmap character vector of colours used in plotting
#' @param plot.theme function for ggplot theme used.
#' @return list of four plots that fit in a 2x2 raster used in MetaboShiny.
ggplotNormSummary <- function(mSet, colmap = global$vectors$mycols, plot.theme = global$functions$plot.themes[[getOptions("user_options.txt")$gtheme]]){
  # load in original data (pre-normalization, post-filter)
  orig_data <- as.data.frame(mSet$dataSet$procr)
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
    theme(legend.position="none",
          axis.text=element_text(size=global$constants$font.aes$ax.num.size),
          axis.title=element_text(size=global$constants$font.aes$ax.txt.size),
          legend.title.align = 0.5,
          axis.line = element_line(colour = 'black', size = .5),
          text = element_text(family = global$constants$font.aes$font))
  
  # second result plot: shows the spread of the intensities before normalization
  RES2 <- plot + ggplot2::geom_boxplot(
    ggplot2::aes(x=value,y=variable),
    color=rainbow(sampsize), 
    alpha=0.4) + ggplot2::geom_vline(ggplot2::aes(xintercept=median(orig_melt$value)))+
    theme(legend.position="none",
          axis.text=element_text(size=global$constants$font.aes$ax.num.size),
          axis.title=element_text(size=global$constants$font.aes$ax.txt.size),
          legend.title.align = 0.5,
          axis.line = element_line(colour = 'black', size = .5),
          text = element_text(family = global$constants$font.aes$font))
  
  # create base plot with base theme and font size for normalized data
  plot <- ggplot2::ggplot(data=norm_melt) +
    plot.theme(base_size = 15) #+ facet_grid(. ~ variable)
  
  # third result plot: a density plot of chosen 20 mz values post normalization
  RES3 <- plot + ggplot2::geom_density(ggplot2::aes(x=value), colour="pink", fill="pink", alpha=0.4)+
    theme(legend.position="none",
          axis.text=element_text(size=global$constants$font.aes$ax.num.size),
          axis.title=element_text(size=global$constants$font.aes$ax.txt.size),
          legend.title.align = 0.5,
          axis.line = element_line(colour = 'black', size = .5),
          text = element_text(family = global$constants$font.aes$font))
  
  # fourth result plot: spread of intensities after normalization
  RES4 <- plot + ggplot2::geom_boxplot(
    ggplot2::aes(x=value,y=variable),
    color=rainbow(sampsize), 
    alpha=0.4) + ggplot2::geom_vline(ggplot2::aes(xintercept=median(norm_melt$value))) +
    theme(legend.position="none",
          axis.text=element_text(size=global$constants$font.aes$ax.num.size),
          axis.title=element_text(size=global$constants$font.aes$ax.txt.size),
          legend.title.align = 0.5,
          axis.line = element_line(colour = 'black', size = .5),
          text = element_text(family = global$constants$font.aes$font))
  
  # - - - - - - - - - - - - - - - - - - -
  
  list(tl=RES1, bl=RES2, tr=RES3, br=RES4)
} 

#' @export
#' Plot a summary of metaboanalystR normalization results. Takes the total of 20 samples before and after normalization and plots that distribution.
#' @param mSet input user mSet
#' @param plot.theme function for ggplot theme used.
#' @return list of four plots that fit in a 2x2 raster used in MetaboShiny.
ggplotSampleNormSummary <- function(mSet, plot.theme = global$functions$plot.themes[[getOptions("user_options.txt")$gtheme]]){
  # 4 by 4 plot, based on random 20-30 picked 
  orig_data <- as.data.frame(mSet$dataSet$procr)
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
    theme(legend.position="none",
          axis.text=element_text(size=global$constants$font.aes$ax.num.size),
          axis.title=element_text(size=global$constants$font.aes$ax.txt.size),
          legend.title.align = 0.5,
          axis.line = element_line(colour = 'black', size = .5),
          text = element_text(family = global$constants$font.aes$font))
  
  RES2 <- ggplot2::ggplot(data=orig_melt) +
    plot.theme(base_size = 15) + ggplot2::geom_boxplot(
      ggplot2::aes(x=value,y=Label),
      color=rainbow(sampsize), 
      alpha=0.4) + ggplot2::geom_vline(ggplot2::aes(xintercept=median(orig_melt$value),text=Label))+
    theme(legend.position="none",
          axis.text=element_text(size=global$constants$font.aes$ax.num.size),
          axis.title=element_text(size=global$constants$font.aes$ax.txt.size),
          legend.title.align = 0.5,
          axis.line = element_line(colour = 'black', size = .5),
          text = element_text(family = global$constants$font.aes$font))
  
  RES3 <- ggplot2::ggplot(data=norm_melt_sums) +
    plot.theme(base_size = 15) + ggplot2::geom_density(ggplot2::aes(x=value), colour="pink", fill="pink", alpha=0.4)+
    theme(legend.position="none",
          axis.text=element_text(size=global$constants$font.aes$ax.num.size),
          axis.title=element_text(size=global$constants$font.aes$ax.txt.size),
          legend.title.align = 0.5,
          axis.line = element_line(colour = 'black', size = .5),
          text = element_text(family = global$constants$font.aes$font))
  
  RES4 <- ggplot2::ggplot(data=norm_melt) +
    plot.theme(base_size = 15) + ggplot2::geom_boxplot(
      ggplot2::aes(x=value,y=Label),
      color=rainbow(sampsize), 
      alpha=0.4) + ggplot2::geom_vline(ggplot2::aes(xintercept=median(norm_melt$value),text=Label))+
    theme(legend.position="none",
          axis.text=element_text(size=global$constants$font.aes$ax.num.size),
          axis.title=element_text(size=global$constants$font.aes$ax.txt.size),
          legend.title.align = 0.5,
          axis.line = element_line(colour = 'black', size = .5),
          text = element_text(family = global$constants$font.aes$font))
  
  list(tl=RES1, bl=RES2, tr=RES3, br=RES4)
} 


#' @export
ggplotMeba <- function(cpd, draw.average=T, cols=global$vectors$mycols, cf = global$functions$color.functions[[getOptions("user_options.txt")$gspec]], plot.theme = global$functions$plot.themes[[getOptions("user_options.txt")$gtheme]], plotlyfy=TRUE){
  cols <- if(is.null(cols)) cf(length(levels(mSet$dataSet$cls))) else(cols)
  profile <- getProfile(cpd, mode="time")
  p <- if(draw.average){
    ggplot2::ggplot(data=profile) +
      ggplot2::geom_line(size=0.3, ggplot2::aes(x=Time, y=Abundance, group=Sample, color=Group, text=Sample), alpha=0.4) +
      stat_summary(fun.y="mean", size=1.5, geom="line", ggplot2::aes(x=Time, y=Abundance, color=Group, group=Group)) +
      ggplot2::scale_x_discrete(expand = c(0, 0)) +
      plot.theme(base_size = 15) +
      ggplot2::scale_color_manual(values=cols) +
      theme(legend.position="none",
            axis.text=element_text(size=global$constants$font.aes$ax.num.size),
            axis.title=element_text(size=global$constants$font.aes$ax.txt.size),
            legend.title.align = 0.5,
            axis.line = element_line(colour = 'black', size = .5),
            text = element_text(family = global$constants$font.aes$font)) +
      ggtitle(cpd)
  } else{
    ggplot2::ggplot(data=profile) +
      ggplot2::geom_line(size=0.7, ggplot2::aes(x=Time, y=Abundance, group=Sample, color=Group, text=Sample)) +
      ggplot2::scale_x_discrete(expand = c(0, 0)) +
      plot.theme(base_size = 15) +
      ggplot2::scale_color_manual(values=cols) +
      theme(legend.position="none",
            axis.text=element_text(size=global$constants$font.aes$ax.num.size),
            axis.title=element_text(size=global$constants$font.aes$ax.txt.size),
            legend.title.align = 0.5,
            axis.line = element_line(colour = 'black', size = .5),
            text = element_text(family = global$constants$font.aes$font)) + 
      ggtitle(cpd)
  }
  if(plotlyfy){
    ggplotly(p, tooltip="Sample", originalData=T)
  }else{
    p
  }
}

blackwhite.colors <- function(n){
  gray.colors(n, start=0, end=1)
}

#' @export
ggplotSummary <- function(cpd = curr_cpd, shape.fac = "label", cols = c("black", "pink"), sourceTable = mSet$dataSet$norm, 
                          cf = rainbow, plot.theme = global$functions$plot.themes[[getOptions("user_options.txt")$gtheme]], 
                          mode = "nm", plotlyfy = TRUE, 
                          styles=c("box", "beeswarm"), add_stats = "mean",
                          col.fac = "label", txt.fac = "label"){
  
  if(length(styles) == 0){
    styles = c("beeswarm")
  }
  # - - -
  
  profile <- getProfile(cpd, mode=if(mode == "nm") "stat" else "time", sourceTable = sourceTable)
  df_line <- data.table(x = c(1,2),
                        y = rep(min(profile$Abundance - 0.1),2))
  stars = ""
  
  try({
    pval <- if(mode == "nm"){
      mSet$analSet$tt$sig.mat[which(rownames(mSet$analSet$tt$sig.mat) == curr_cpd), "p.value"]
    }else{
      mSet$analSet$aov2$sig.mat[which(rownames(mSet$analSet$aov2$sig.mat) == curr_cpd), 'Interaction(adj.p)']
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
                     ts = split(profile, f = profile$Group),
                     nm = list(x = data.table(profile)))

  i = 1
  
  for(prof in profiles){
    
    for(style in styles){
      
      switch(mode,
             nm = {
               p <- p + switch(style,
                               box = ggplot2::geom_boxplot(data = prof, alpha=0.4, aes(x = Group,
                                                                                       y = Abundance,
                                                                                       text = Text,
                                                                                       color = Group, 
                                                                                       fill = Group)),
                               violin = ggplot2::geom_violin(data = prof, alpha=0.4, position = "identity", aes(x = Group,
                                                                                                                y = Abundance,
                                                                                                                color = Group,
                                                                                                                fill = Group)),
                               beeswarm = ggbeeswarm::geom_beeswarm(data = prof, alpha=0.7, size = 2, position = position_dodge(width=.3), aes(x = Group,
                                                                                                                                               y = Abundance,
                                                                                                                                               text= Text,
                                                                                                                                               color = Group, 
                                                                                                                                               fill = Color)),
                               scatter = ggplot2::geom_point(data = prof, alpha=0.7, size = 2, aes(x = Group, 
                                                                                                   y = Abundance, 
                                                                                                   text=Text, 
                                                                                                   color = Group, 
                                                                                                   fill = Color), position = position_jitterdodge())
               )
             },
             ts = {
               p <- p + switch(style,
                               box = ggplot2::geom_boxplot(data = prof, alpha=0.4, aes(x = Time,
                                                                                       y = Abundance,
                                                                                       text = Text,
                                                                                       color = Group, 
                                                                                       fill = Group)),
                               violin = ggplot2::geom_violin(data = prof, alpha=0.4, position = "identity", aes(x = Time,
                                                                                                                y = Abundance,
                                                                                                                group = Time,
                                                                                                                color = Group, 
                                                                                                                fill = Group)),
                               beeswarm = ggbeeswarm::geom_beeswarm(data = prof, alpha=0.7, size = 2, position = position_dodge(width=.3), aes(x = Time,
                                                                                                                                               y = Abundance,
                                                                                                                                               text=Text,
                                                                                                                                               color = Group, 
                                                                                                                                               fill = Color)),
                               scatter = ggplot2::geom_point(data = prof, alpha=0.7, size = 2, position = position_jitterdodge(), aes(x = Time,
                                                                                                                                      y = Abundance,
                                                                                                                                      text=Text,
                                                                                                                                      color = Group, 
                                                                                                                                      fill = Color))
               )
             })

      }
    
    p <- p + theme(legend.position="none",
                   axis.text=element_text(size=global$constants$font.aes$ax.num.size),
                   axis.title=element_text(size=global$constants$font.aes$ax.txt.size),
                   legend.title.align = 0.5,
                   axis.line = element_line(colour = 'black', size = .5),
                   text = element_text(family = global$constants$font.aes$font))+ 
      ggplot2::annotate("text", 
                        x = switch(mode, nm = 1.5,
                                   ts = max(as.numeric(profile$Time))/2 + .5), 
                        y = min(profile$Abundance - 0.3), 
                        label = stars, size = 8, col = "black") + ggtitle(cpd)
    
    if(!("box" %in% styles)){
      p <- switch(add_stats,
                  median = {
                    p + stat_summary(data = prof,
                                     aes( x = if(mode == "nm") Group else Time,
                                          y = Abundance,
                                          color = Group),
                                     fun.y = median, 
                                     fun.ymin = median, 
                                     fun.ymax = median,
                                     geom = "crossbar", 
                                     width = 0.5, 
                                     color = switch(mode, 
                                                    ts = cols[i],
                                                    nm = cols[1:length(unique(levels(profile$Group)))]))
                  },
                  mean = {
                    p + stat_summary(data = prof,
                                     aes( x = if(mode == "nm") Group else Time,
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
  
  p <- p + ggplot2::scale_color_manual(values=cols[1:length(unique(levels(profile$Group)))])
  
  # print(cpd)
  # print(col.fac)
  # print(shape.fac)
  # print(txt.fac)
  # print(styles)

  if(all(as.character(profile$Color) == as.character(profile$Group))){
    p <- p + ggplot2::scale_fill_manual(values = cols[1:length(unique(levels(profile$Group)))])
  }else{
    ncols = length(levels(profile$Color))
    scale = c(cols[1:length(unique(levels(profile$Group)))],cf(ncols))
    names(scale) <- c(levels(profile$Group),levels(profile$Color))
    p <- p + ggplot2::scale_fill_manual(values = scale)
  }
  # ---------------
  if(plotlyfy){
    ggplotly(p, tooltip = "Text")#, originalData=T)
  }else{
    p
  }
}

ggPlotTT <- function(cf=global$functions$color.functions[[getOptions("user_options.txt")$gspec]], n, plot.theme = global$functions$plot.themes[[getOptions("user_options.txt")$gtheme]], plotlyfy=TRUE){
  profile <- as.data.table(mSet$analSet$tt$p.log[mSet$analSet$tt$inx.imp],keep.rownames = T)
  colnames(profile) <- c("cpd", "p")
  profile$Peak <- c(1:nrow(profile)) 
  # ---------------------------
  p <- ggplot2::ggplot(data=profile) +
    ggplot2::geom_point(ggplot2::aes(x=Peak, y=p,text=cpd, color=p, key=cpd)) +
    plot.theme(base_size = 15) +
    theme(legend.position="none",
          axis.text=element_text(size=global$constants$font.aes$ax.num.size),
          axis.title=element_text(size=global$constants$font.aes$ax.txt.size),
          legend.title.align = 0.5,
          axis.line = element_line(colour = 'black', size = .5),
          text = element_text(family = global$constants$font.aes$font))+
    ggplot2::scale_colour_gradientn(colours = cf(n)) +
    ggplot2::scale_y_log10()
  if(plotlyfy){
    ggplotly(p, tooltip="cpd")
  }else{
    p
  }
}

ggPlotFC <- function(cf=global$functions$color.functions[[getOptions("user_options.txt")$gspec]], n, plot.theme = global$functions$plot.themes[[getOptions("user_options.txt")$gtheme]], plotlyfy=TRUE){
  profile <- as.data.table(mSet$analSet$fc$fc.log[mSet$analSet$fc$inx.imp],keep.rownames = T)
  profile
  colnames(profile) <- c("cpd", "log2fc")
  profile$Peak <- c(1:nrow(profile)) 
  # ---------------------------
  p <- ggplot2::ggplot(data=profile) +
    ggplot2::geom_point(ggplot2::aes(x=Peak, y=log2fc, text=log2fc, color=log2fc, key=cpd)) +
    ggplot2::geom_abline(ggplot2::aes(intercept = 0, slope = 0)) +
    plot.theme(base_size = 15) +
    theme(legend.position="none",
          axis.text=element_text(size=global$constants$font.aes$ax.num.size),
          axis.title=element_text(size=global$constants$font.aes$ax.txt.size),
          legend.title.align = 0.5,
          axis.line = element_line(colour = 'black', size = .5),
          text = element_text(family = global$constants$font.aes$font))+
    ggplot2::scale_colour_gradientn(colours = cf(n))
  
  if(plotlyfy){
    ggplotly(p, tooltip="log2fc")
  }else{
    p
  }
}

ggPlotVolc <- function(cf=global$functions$color.functions[[getOptions("user_options.txt")$gspec]], n=256, plot.theme = global$functions$plot.themes[[getOptions("user_options.txt")$gtheme]], plotlyfy=TRUE){
  vcn<-mSet$analSet$volcano;
  dt <- as.data.table(vcn$sig.mat[,c(2,4)],keep.rownames = T)
  colnames(dt) <- c("cpd", "log2FC", "-log10P")
  p <- ggplot2::ggplot() +
    #ggplot2::geom_point(data=dt[!imp.inx], ggplot2::aes(x=log2FC, y=minlog10P)) +
    ggplot2::geom_point(data=dt, ggplot2::aes(x=log2FC, 
                                              y=`-log10P`,
                                              text=cpd,
                                              color=abs(log2FC*`-log10P`), 
                                              key=cpd)) +
    plot.theme(base_size = 15) +
    theme(legend.position="none",
          axis.text=element_text(size=global$constants$font.aes$ax.num.size),
          axis.title=element_text(size=global$constants$font.aes$ax.txt.size),
          legend.title.align = 0.5,
          axis.line = element_line(colour = 'black', size = .5),
          text = element_text(family = global$constants$font.aes$font))+
    ggplot2::scale_colour_gradientn(colours = cf(n),guide=FALSE)
  
  if(plotlyfy){
    ggplotly(p, tooltop="cpd")
  }else{
    p
  }
}

ggPlotPCApairs <- function(cols = c("black", "pink"), 
                           pc.num, 
                           type = "pca", 
                           cf = rainbow, 
                           plot.theme = global$functions$plot.themes[[getOptions("user_options.txt")$gtheme]],
                           plotlyfy=TRUE){
  
  cols <- if(is.null(cols)) cf(length(levels(mSet$dataSet$cls))) else(cols)
  
  ggplot <- function(...) ggplot2::ggplot(...) + scale_color_brewer(palette="Set1")
  unlockBinding("ggplot",parent.env(asNamespace("GGally")))
  assign("ggplot",ggplot,parent.env(asNamespace("GGally")))
  
  switch(type,
         pca = {
           pclabels <- paste0("PC", 1:pc.num, " (", round(100 * mSet$analSet$pca$variance[1:pc.num], 
                                                          1), "%)")
           pca_matr <- as.data.table(mSet$analSet$pca$x[, 1:pc.num])
         },
         plsca = {
           pclabels <- paste0("PC", 1:ncol(mSet$analSet$plsr$scores), " (", round(mSet$analSet$plsr$Xvar 
                                                                                  / mSet$analSet$plsr$Xtotvar 
                                                                                  * 100.0,
                                                                                  digits = 2), "%)")
           pca_matr <- as.data.table(mSet$analSet$plsr$scores[, 1:pc.num])
         },
         oplsda = {
           # pclabels <- paste0("PC", 1:ncol(mSet$analSet$oplsda), " (", round(mSet$analSet$plsr$Xvar 
           #                                                                        / mSet$analSet$plsr$Xtotvar 
           #                                                                        * 100.0,
           #                                                                        digits = 2), "%)")
           # pca_matr <- as.data.table(mSet$analSet$plsr$scores[, 1:pc.num])
         },
         splsda = {
           # pclabels <- paste0("PC", 1:ncol(mSet$analSet$splsr$variates$X), " (", round(mSet$analSet$splsr$explained_variance
           #                                                                   * 100.0,
           #                                                                   digits = 2), "%)")
           # pca_matr <- as.data.table(mSet$analSet$plsr$scores[, 1:pc.num])
           # 
           # mSet$analSet$splsr$variates
           # 
         }
  )
  colnames(pca_matr)  <- pclabels
  pca_matr$class <- mSet$dataSet$cls
  p <- GGally::ggpairs(pca_matr, 
                       columns = c(1:(ncol(pca_matr)-1)),
                       aes(colour = class, 
                           alpha = 0.4),
                       upper = list(continuous = wrap("density", alpha = 0.5), combo = "box"),
                       lower = list(continuous = wrap("points", alpha = 0.3), combo = wrap("dot", alpha = 0.4)),
                       diag = list(continuous = wrap("densityDiag"))
  )
  
  for(i in 1:p$nrow) {
    for(j in 1:p$ncol){
      p[i,j] <- p[i,j] + 
        scale_fill_manual(values=c(cols)) +
        scale_color_manual(values=c(cols))  
    }
  }
  
  if(plotlyfy){
    ggplotly(p)
  }else{
    p
  }
}

ggPlotClass <- function(pls.type = "plsda", 
                        cf = global$functions$color.functions[[getOptions("user_options.txt")$gspec]], 
                        pcs = 3, 
                        plot.theme = global$functions$plot.themes[[getOptions("user_options.txt")$gtheme]],
                        plotlyfy=TRUE){
  res <- mSet$analSet$plsda$fit.info
  colnames(res) <- 1:ncol(res)
  # best.num <- mSet$analSet$plsda$best.num
  # choice <- mSet$analSet$plsda$choice
  df <- melt(res)
  df$Component <- paste0("PC",df$Component)
  colnames(df) <- c("Metric", "Component", "Value")
  p <- ggplot(df, aes(x=Metric, y=Value, fill=Metric)) +
    geom_bar(stat="identity") + 
    theme_minimal() + 
    plot.theme(base_size = 15) +
    theme(legend.position="none",
          axis.text=element_text(size=global$constants$font.aes$ax.num.size),
          axis.title=element_text(size=global$constants$font.aes$ax.txt.size),
          legend.title.align = 0.5,
          axis.line = element_line(colour = 'black', size = .5),
          text = element_text(family = global$constants$font.aes$font))+
    facet_grid(~Component) + 
    scale_fill_manual(values=cf(pcs))
  if(plotlyfy){
    ggplotly(p)
  }else{
    p
  }
}

ggPlotPerm <- function(pls.type = "plsda", 
                       cf = global$functions$color.functions[[getOptions("user_options.txt")$gspec]], 
                       pcs = 3, 
                       plot.theme = global$functions$plot.themes[[getOptions("user_options.txt")$gtheme]],
                       plotlyfy=TRUE){
  bw.vec <- mSet$analSet$plsda$permut
  len <- length(bw.vec)
  df <- melt(bw.vec)
  colnames(df) = "acc"
  # round p value
  pval <- mSet$analSet$plsda$permut.p
  rounded <- round(as.numeric(str_match(pval, "0\\.\\d*")), digits = 3)
  pval <- gsub(pval, pattern = "(0\\.\\d*)", replacement=rounded)
  # - - -
  p <- ggplot(df) +
    geom_histogram(mapping=aes(x=acc, y=..count.., fill=factor(..count..)),
                   binwidth=0.01) +
    scale_fill_manual(values=cf(20)) +
    plot.theme(base_size = 15) + 
    theme(legend.position="none",
          axis.text=element_text(size=global$constants$font.aes$ax.num.size),
          axis.title=element_text(size=global$constants$font.aes$ax.txt.size),
          legend.title.align = 0.5,
          axis.line = element_line(colour = 'black', size = .5),
          text = element_text(family = global$constants$font.aes$font))+
    labs(x="Accuracy", y = "Permutations") +
    geom_segment(data=df,
                 color="black",
                 x=bw.vec[1],
                 xend=bw.vec[1],
                 y=0,
                 aes(yend=.1*nrow(df)),
                 size=1.5,
                 linetype=8) +
    geom_text(mapping = aes(x = bw.vec[1], y =  .11*nrow(df), label = pval), color = "black", size = 4)
  
  if(plotlyfy){
    ggplotly(p)
  }else{
    p
  }
}

plot.many <- function(res.obj = models, 
                      which_alpha = 1, 
                      plot.theme = global$functions$plot.themes[[getOptions("user_options.txt")$gtheme]],
                      plotlyfy=TRUE){
  
  predictions <- if(length(res.obj) > 1) do.call("cbind", lapply(res.obj, function(x) x$prediction)) else data.frame(res.obj[[1]]$prediction)
  
  colnames(predictions) <- if(length(res.obj) > 1) sapply(res.obj, function(x) x$alpha) else res.obj[[1]]$alpha
  testY = res.obj[[1]]$labels
  
  if(length(unique(testY)) > 2){
    
    return("not supported yet")
    # https://stats.stackexchange.com/questions/112383/roc-for-more-than-2-outcome-categories
    #
    # for (type.id in length(unique(testY))) {
    #   type = as.factor(iris.train$Species == lvls[type.id])
    #   
    #   nbmodel = NaiveBayes(type ~ ., data=iris.train[, -5])
    #   nbprediction = predict(nbmodel, iris.test[,-5], type='raw')
    #   
    #   score = nbprediction$posterior[, 'TRUE']
    #   actual.class = iris.test$Species == lvls[type.id]
    #   
    #   pred = prediction(score, actual.class)
    #   nbperf = performance(pred, "tpr", "fpr")
    #   
    #   roc.x = unlist(nbperf@x.values)
    #   roc.y = unlist(nbperf@y.values)
    #   lines(roc.y ~ roc.x, col=type.id+1, lwd=2)
    #   
    #   nbauc = performance(pred, "auc")
    #   nbauc = unlist(slot(nbauc, "y.values"))
    #   aucs[type.id] = nbauc
    # }
  }else{
    data <- data.frame(D = as.numeric(as.factor(testY))-1,
                       D.str = testY)
    data <- cbind(data, predictions)
    if(length(res.obj) > 1){
      roc_coord <- plotROC::melt_roc(data, "D", m = 3:ncol(data))
    }else{
      roc_coord <- data.frame(D = rep(data[, "D"], length(3)),
                              M = data[, 3], 
                              name = rep(names(data)[3], each = nrow(data)), 
                              stringsAsFactors = FALSE)
    }
  }
  
  names(roc_coord)[which(names(roc_coord) == "name")] <- "alpha"
  
  roc_coord <- roc_coord[roc_coord$alpha %in% which_alpha,]
  # plot
  p <- ggplot2::ggplot(roc_coord, 
                       ggplot2::aes(d = D, m = M, color = alpha)) + 
    plotROC::geom_roc(labelsize=0,show.legend = TRUE) + 
    plotROC::style_roc() + 
    theme(axis.text=element_text(size=19),
          axis.title=element_text(size=19,face="bold"),
          legend.title=element_text(size=19),
          legend.text=element_text(size=19)) +
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE)
  
  if(plotlyfy){
    ggplotly(p)
  }else{
    p
  }
}


ggPlotROC <- function(data, 
                      attempts = 50, 
                      cf = global$functions$color.functions[[getOptions("user_options.txt")$gspec]], 
                      plot.theme = global$functions$plot.themes[[getOptions("user_options.txt")$gtheme]],
                      plotlyfy=TRUE){
  
  require(ROCR)
  require(ggplot2)
  require(data.table)
  
  
  pred <- ROCR::prediction(data$predictions, data$labels)
  perf <- ROCR::performance(pred, "tpr", "fpr")
  perf_auc <- ROCR::performance(pred, "auc")
  
  # - - - old method - - -
  perf.long <- rbindlist(lapply(1:length(perf@x.values), function(i){
    xvals <- perf@x.values[[i]]
    yvals <- perf@y.values[[i]]
    aucs <- signif(perf_auc@y.values[[i]][[1]], digits = 2)
    
    res <- data.table::data.table(attempt = c(i),
                                  FPR = xvals,
                                  TPR = yvals,
                                  AUC = aucs)
    res
  }))
  
  aucs <- signif(unlist(perf_auc@y.values), digits = 2)
  
  cols = cf(attempts)
  
  p <- ggplot(perf.long, aes(FPR,TPR)) +
    stat_summary_bin(aes(FPR, TPR), fun.y=mean, geom="line", colour="black", cex = 2) +
    geom_path(alpha=.5,
              cex=.5,
              aes(color = attempt, group = attempt)) +
    geom_text(label = paste0("Average AUC: ", mean(aucs)), size = 8, mapping = aes(x = 0.82, y = 0.03)) + 
    plot.theme(base_size = 10) +
    ggplot2::scale_color_gradientn(colors = cols) +
    theme(legend.position="none",
          axis.text=element_text(size=global$constants$font.aes$ax.num.size),
          axis.title=element_text(size=global$constants$font.aes$ax.txt.size),
          legend.title.align = 0.5,
          axis.line = element_line(colour = 'black', size = .5),
          text = element_text(family = global$constants$font.aes$font))+
    #scale_y_continuous(limits=c(0,1)) +
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE) + 
    coord_cartesian(xlim = c(.04,.96), ylim = c(.04,.96))
  
  if(plotlyfy){
    ggplotly(p)
  }else{
    p
  }
}

ggPlotBar <- function(data, 
                      attempts=50, 
                      cf = global$functions$color.functions[[getOptions("user_options.txt")$gspec]], 
                      topn=50, 
                      plot.theme = global$functions$plot.themes[[getOptions("user_options.txt")$gtheme]],
                      ml_name, 
                      ml_type,
                      plotlyfy=TRUE){
  
  if(ml_name != ""){
    lname = ml_name
  }else{
    lname <- paste0(input$ml_train_regex,"|", input$ml_test_regex)
  }
  
  if(ml_type == "ls"){
    p <- ggplot(data[1:topn,], aes(mz,count)) + geom_bar(stat = "identity", aes(fill = count)) +
      scale_fill_gradientn(colors=cf(20)) +
      plot.theme() + 
      theme(legend.position="none",
            axis.text=element_text(size=global$constants$font.aes$ax.num.size),
            axis.title=element_text(size=global$constants$font.aes$ax.txt.size),
            legend.title.align = 0.5,
            axis.line = element_line(colour = 'black', size = .5),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            text = element_text(family = global$constants$font.aes$font))+
      labs(x="Top hits",y="Times chosen")
    
  }else{
    p <- ggplot(data[1:topn,], aes(mz,mda)) + geom_bar(stat = "identity", aes(fill = mda)) +
      scale_fill_gradientn(colors=cf(20)) +
      plot.theme() + 
      theme(legend.position="none",
            axis.text=element_text(size=global$constants$font.aes$ax.num.size),
            axis.title=element_text(size=global$constants$font.aes$ax.txt.size),
            legend.title.align = 0.5,
            axis.line = element_line(colour = 'black', size = .5),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            text = element_text(family = global$constants$font.aes$font))+
      labs(x="Top hits",y="Mean Decrease Accuracy")
  }
  if(plotlyfy){
    ggplotly(p)
  }else{
    p
  }
}

plotPCA.3d <- function(mSet, 
                       cols = global$vectors$mycols, 
                       shape.fac="label", 
                       pcx, pcy, pcz, 
                       mode="pca"){
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
           plsda.table <- as.data.table(round(mSet$analSet$plsr$Xvar 
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
    #df_split_idx <- split(1:nrow(df), f = sapply(strsplit(rownames(df), split = "_T"), function(x) x[[2]]))
    #df_list <- lapply(df_split_idx, function(idx_list) as.data.frame(df[idx_list,]))
  }else{
    fac.lvls <- length(levels(mSet$dataSet$cls))
    classes = mSet$dataSet$cls
  }
  
  df <- as.data.frame(df)
  rownames(df) <- rownames(mSet$dataSet$norm)
  df_list <- list(df)
  
  cols <- cols[c(1:length(unique(classes)))]
  
  # --- add ellipses ---
  
  symbols = c('circle',
              'diamond', 
              'square', 
              'x',
              'o')
  
  symbol.vec<-if(is.null(shape.fac)){
    rep('circle', times = length(classes))
  }
  else if(shape.fac == "label"){
    rep('circle', times = length(classes))
  }else{
    as.factor(mSet$dataSet$covars[,..shape.fac][[1]])
  }
  
  plots_facet <- lapply(1:length(df_list), function(i){
    
    df = df_list[[i]]
    
    orig_idx = which(rownames(df) %in% rownames(mSet$dataSet$norm))
    
    plots <- plotly::plot_ly(showlegend = F)
    
    for(class in levels(classes)){
      
      samps <- rownames(mSet$dataSet$norm)[which(classes == class)]
      row = which(rownames(df) %in% samps)
      # ---------------------
      xc=df[row, pcx]
      yc=df[row, pcy]
      zc=df[row, pcz]
      
      # --- plot ellipse ---
      o <- rgl::ellipse3d(cov(cbind(xc,yc,zc)), 
                          centre=c(mean(xc), 
                                   mean(yc), 
                                   mean(zc)), 
                          level = 0.95)
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
    adj_plot <<- plotly_build(plots)
    rgbcols <- toRGB(cols)
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
      family = global$constants$font.aes$font
    )
    # --- return ---
    pca_plot <<- adj_plot %>% add_trace(
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
    ) %>%  layout(font=t,
                  scene = list(
                    aspectmode="cube",
                    xaxis = list(
                      titlefont = list(size = global$constants$font.aes$ax.txt.size * 1.5),
                      title = gsubfn::fn$paste("$pcx ($x.var %)")),
                    yaxis = list(
                      titlefont = list(size = global$constants$font.aes$ax.txt.size * 1.5),
                      title = gsubfn::fn$paste("$pcy ($y.var %)")),
                    zaxis = list(
                      titlefont = list(size = global$constants$font.aes$ax.txt.size * 1.5),
                      title = gsubfn::fn$paste("$pcz ($z.var %)")))) 
    # --- return ---
    pca_plot
  })
  subplot(plots_facet,shareX = F, shareY = F)
}

plotPCA.2d <- function(mSet, shape.fac = "label", cols = global$vectors$mycols, pcx, pcy, mode="pca", plot.theme = global$functions$plot.themes[[getOptions("user_options.txt")$gtheme]], plotlyfy="T"){
  
  classes <- switch(mode, 
                    ipca = mSet$dataSet$exp.fac,
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
           fac.lvls <- length(levels(mSet$dataSet$exp.fac))
           
           xc=mSet$analSet$pca$x[, pcx]
           yc=mSet$analSet$pca$x[, pcy]
           
           dat_long <- data.table(variable = names(xc),
                                  group = classes,
                                  time = mSet$dataSet$time.fac,
                                  x = xc, 
                                  y = yc)
           
         },
         plsda = {
           plsda.table <- as.data.table(round(mSet$analSet$plsr$Xvar 
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
  
  p <- ggplot(dat_long, aes(x, y,group=group)) +
    geom_point(size=5, aes(shape=shape,
                           text=variable,
                           fill=group, 
                           color=group), alpha=0.7)+
    stat_ellipse(geom = "polygon", aes(fill=group), alpha = 0.3) +
    plot.theme() +
    theme(legend.position="none",
          axis.text=element_text(size=global$constants$font.aes$ax.num.size),
          axis.title=element_text(size=global$constants$font.aes$ax.txt.size),
          legend.title.align = 0.5,
          axis.line = element_line(colour = 'black', size = .5),
          text = element_text(family = global$constants$font.aes$font))+
    scale_x_continuous(name=gsubfn::fn$paste("$pcx ($x.var%)")) +
    scale_y_continuous(name=gsubfn::fn$paste("$pcy ($y.var%)")) +
    scale_fill_manual(values=cols) +
    scale_color_manual(values = cols)
  
  if(mode == "ipca") p <- p + facet_wrap(~time)
  
  if(plotlyfy){
    ggplotly(p)
  }else{
    p
  }
  
}

ggPlotVenn <- function(mSet, 
                       venn_yes, 
                       top=100, 
                       cols=global$vectors$mycols, 
                       cf = global$functions$color.functions[[getOptions("user_options.txt")$gspec]], 
                       plotlyfy=TRUE){
  
  #venn_yes <<- isolate({venn_yes})
  
  experiments <- str_match(unlist(venn_yes$now), pattern = "\\(.*\\)")[,1]
  
  experiments <- unique(gsub(experiments, pattern = "\\(\\s*(.+)\\s*\\)", replacement="\\1"))
  
  table_list <- lapply(experiments, function(experiment){
    
    analysis = mSet$storage[[experiment]]$analysis
    #data = mSet$storage[[experiment]]$dataset
    
    rgx_exp <- gsub(experiment, pattern = "\\(", replacement = "\\\\(")
    rgx_exp <- gsub(rgx_exp, pattern = "\\)", replacement = "\\\\)")
    
    categories = grep(unlist(venn_yes$now), 
                      pattern = rgx_exp, value = T)
    
    categories = gsub(categories, pattern = "\\(\\s*(.+)\\s*\\)", replacement = "")
    
    # go through the to include analyses
    tables <- lapply(categories, function(name){
      
      base_name <- gsub(name, pattern = " -.*$| ", replacement="")
      
      # fetch involved mz values 
      tbls <- switch(base_name,
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
                     ls = {
                       which.ls <- gsub(name, pattern = "^.*- | ", replacement="")
                       res <- list(as.character(analysis$ml$ls[[which.ls]]$bar[order(analysis$ml$ls[[which.ls]]$bar$count, 
                                                                                     decreasing = T),]$mz))
                       names(res) <- paste0(which.ls, " (LS)")
                       # - - - 
                       res
                       
                     },
                     rf = {
                       which.rf <- gsub(name, pattern = "^.*- | ", replacement="")
                       res <- list(as.character(analysis$ml$rf[[which.rf]]$bar[order(analysis$ml$rf[[which.rf]]$bar$mda, 
                                                                                     decreasing = T),]$mz))
                       names(res) <- paste0(which.rf, " (RF)")
                       # - - -
                       res
                     },
                     plsda = {
                       
                       which.plsda <- gsub(name, pattern = "^.*- | ", replacement="")
                       
                       compounds_pc <- as.data.table(analysis$plsda$vip.mat,keep.rownames = T)
                       colnames(compounds_pc) <- c("rn", paste0("PC", 1:(ncol(compounds_pc)-1)))
                       ordered_pc <- setorderv(compounds_pc, which.plsda, -1)
                       
                       res <- list(ordered_pc$rn)
                       names(res) <- paste0(which.plsda, " (PLS-DA)")
                       # - - -
                       res
                     },
                     volc = {
                       res <- list(rownames(analysis$volcano$sig.mat))
                       names(res) = base_name
                       res
                     })
      
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
  
  # TODO: this is sloppy... IDK
  global$vectors$venn_lists <<- flattened
  
  # how many circles need to be plotted? (# of included analysis)
  circles = length(flattened)
  
  # generate the initial plot - the POLYGONS
  venn.plot <- VennDiagram::venn.diagram(x = flattened,
                                         filename = NULL)
  
  # split the plots into its individual elements
  items <- strsplit(as.character(venn.plot), split = ",")[[1]]
  
  # get which are circles
  circ_values <<- data.frame(
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
                  y = y)) + geom_polygon(colour="black", alpha=0.5, aes(fill=id, group=id)) +
    geom_text(mapping = aes(x=x-.02, y=y, label=value), data = numbers, size = 6, hjust = 0) +
    theme_void() +
    theme(legend.position="none",
          text=element_text(size=),
          panel.grid = element_blank()) +
    scale_fill_gradientn(colours = 
                           cf(circles)) +
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE)
  
  # - - text with white outline - - 
  p <- p + shadowtext::geom_shadowtext(mapping = aes(x=x, y=y, label=gsub(value, 
                                                                          pattern = ":", 
                                                                          replacement = "\n")), 
                                       color = "black",
                                       bg.color = "white",
                                       data = headers, 
                                       size = 7, hjust = 0.5, 
                                       fontface = 2) +
    scale_x_continuous(expand = c(.1, .1)) + 
    scale_y_continuous(expand = c(.1, .1))
  if(plotlyfy){
    ggplotly(p)
  }else{
    p
  }
}