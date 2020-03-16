#' @export
#' Plot a summary of metaboanalystR normalization results. Takes the total of 20 m/z values and 20 samples before and after normalization and plots that distribution.
#' @param mSet input user mSet
#' @return list of four plots that fit in a 2x2 raster used in MetaboShiny.
ggplotNormSummary <- function(mSet,
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
  plot <- ggplot2::ggplot(data=orig_melt)
  
  # first result plot: is a density plot of chosen 20 mz values with 20 samples
  RES1 <- plot + ggplot2::geom_density(ggplot2::aes(x=value), colour="blue", fill="blue", alpha=0.4)
  
  # second result plot: shows the spread of the intensities before normalization
  RES2 <- plot + ggplot2::geom_boxplot(
    ggplot2::aes(x=value,y=variable),
    color=cf(sampsize),
    alpha=0.4) + ggplot2::geom_vline(ggplot2::aes(xintercept=median(orig_melt$value)))
  
  # create base plot with base theme and font size for normalized data
  plot <- ggplot2::ggplot(data=norm_melt)
  
  # third result plot: a density plot of chosen 20 mz values post normalization
  RES3 <- plot + ggplot2::geom_density(ggplot2::aes(x=value), colour="pink", fill="pink", alpha=0.4)
  # fourth result plot: spread of intensities after normalization
  RES4 <- plot + ggplot2::geom_boxplot(
    ggplot2::aes(x=value,y=variable),
    color=cf(sampsize),
    alpha=0.4) + ggplot2::geom_vline(ggplot2::aes(xintercept=median(norm_melt$value)))

    list(tl=RES1, bl=RES2, 
         tr=RES3, br=RES4)
  
}

#' @export
#' Plot a summary of metaboanalystR normalization results. Takes the total of 20 samples before and after normalization and plots that distribution.
#' @param mSet input user mSet
#' @return list of four plots that fit in a 2x2 raster used in MetaboShiny.
ggplotSampleNormSummary <- function(mSet,
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
     ggplot2::geom_density(ggplot2::aes(x=value), colour="blue", fill="blue", alpha=0.4)
  
  RES2 <- ggplot2::ggplot(data=orig_melt) +
     ggplot2::geom_boxplot(
      ggplot2::aes(x=value,y=Label),
      color=cf(sampsize),
      alpha=0.4) + ggplot2::geom_vline(ggplot2::aes(xintercept=median(orig_melt$value),text=Label))
  
  RES3 <- ggplot2::ggplot(data=norm_melt_sums) +
     ggplot2::geom_density(ggplot2::aes(x=value), colour="pink", fill="pink", alpha=0.4)
  
  RES4 <- ggplot2::ggplot(data=norm_melt) +
     ggplot2::geom_boxplot(
      ggplot2::aes(x=value,y=Label),
      color=cf(sampsize),
      alpha=0.4) + ggplot2::geom_vline(ggplot2::aes(xintercept=median(norm_melt$value),text=Label))
  
    list(tl=RES1, bl=RES2, 
         tr=RES3, br=RES4)
  
}


#' @export
ggplotMeba <- function(mSet, cpd, draw.average=T, cols,
                       cf){
  
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
  p
}

blackwhite.colors <- function(n){
  gray.colors(n, start=0, end=1)
}

#' @export
ggplotSummary <- function(mSet, cpd, shape.fac = "label", cols = c("black", "pink"),
                          cf = rainbow, 
                          mode = "nm", 
                          styles=c("box", "beeswarm"), add_stats = "mean",
                          color.fac = "label",
                          text.fac = "label"){
  
  sourceTable = mSet$dataSet$norm
  
  if(length(styles) == 0){
    styles = c("beeswarm")
  }
  # - - -
  
  if(mSet$dataSet$exp.type %in% c("t1f", "2f")){
    mode = "multi"
  }
  
  profile <- MetaboShiny::getProfile(mSet, 
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
  
  for(adj in c("color", "shape", "text")){
    adj.fac = switch(adj,
                     "shape" = shape.fac,
                     "color" = color.fac,
                     "text" = text.fac)
    profile[[Hmisc::capitalize(adj)]] <- if(adj.fac != "label") as.factor(mSet$dataSet$covars[[adj.fac]]) else switch(mode, multi = profile$GroupA, nm = profile$Group)
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
                             box = p + ggplot2::geom_boxplot(data = prof, alpha=0.4, aes(x = Group,
                                                                                         y = Abundance,
                                                                                         shape = Shape,
                                                                                         text = Text,
                                                                                         color = Group,
                                                                                         fill = Color)),
                             violin = p + ggplot2::geom_violin(data = prof, alpha=0.4, position = "identity", aes(x = Group,
                                                                                                                  y = Abundance,
                                                                                                                  color = Group,
                                                                                                                  fill = Color)),
                             beeswarm = p + ggbeeswarm::geom_beeswarm(data = prof, alpha=0.7, size = 2, 
                                                                      #position = position_dodge(width=.3), 
                                                                      aes(x = Group,
                                                                          y = Abundance,
                                                                          text = Text,
                                                                          shape = Shape,
                                                                          color = Group,
                                                                          fill = Color)),
                             scatter = p + ggplot2::geom_point(data = prof, alpha=0.7, size = 2, aes(x = Group,
                                                                                                     y = Abundance,
                                                                                                     text=Text,
                                                                                                     shape = Shape,
                                                                                                     color = Group,
                                                                                                     fill = Color), position = position_jitterdodge())
                 )
               },
               multi = {
                 p <- switch(style,
                             box = p + ggplot2::geom_boxplot(data = prof, alpha=0.4, aes(x = GroupB,
                                                                                         y = Abundance,
                                                                                         text = Text,
                                                                                         shape = Shape,
                                                                                         color = Color,
                                                                                         fill = GroupA)),
                             violin = p + ggplot2::geom_violin(data = prof, alpha=0.4, position = "identity", aes(x = GroupB,
                                                                                                                  y = Abundance,
                                                                                                                  group = GroupB,
                                                                                                                  text = Text,
                                                                                                                  color = Color,
                                                                                                                  fill = GroupA)),
                             beeswarm = {
                               p + ggbeeswarm::geom_beeswarm(data = prof, alpha=0.7, size = 2, position = position_dodge(width=.3), aes(x = GroupB,
                                                                                                                                        y = Abundance,
                                                                                                                                        text=Text,
                                                                                                                                        shape = Shape,
                                                                                                                                        color = Color,
                                                                                                                                        fill = GroupA))},
                             scatter = p + ggplot2::geom_point(data = prof, alpha=0.7, size = 2, position = position_jitterdodge(), aes(x = GroupB,
                                                                                                                                        y = Abundance,
                                                                                                                                        text=Text,
                                                                                                                                        shape = Shape,
                                                                                                                                        color = Color,
                                                                                                                                        fill = GroupA))
                 )
               })
      }
      
      p <- p + ggplot2::annotate("text",
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
    p  
  })
}

ggPlotAOV <- function(mSet, cf, n=20){
  
  which_aov = if(mSet$dataSet$exp.type %in% c("t", "2f", "t1f")) "aov2" else "aov"
  
  profile <- data.table::as.data.table(mSet$analSet[[which_aov]]$p.log[mSet$analSet[[which_aov]]$inx.imp],
                                       keep.rownames = T)
  
  if(nrow(profile)==0){
    shiny::showNotification("No significant hits")
    return(NULL)
  }
  
  profile[,2] <- round(profile[,2], digits = 2)
  profile$Peak <- c(1:nrow(profile))
  colnames(profile)[1:2] <- c("m/z", "-log(p)")
  scaleFUN <- function(x) sprintf("%.1f", x)
  
  xaxis = seq(0,600, 50)
  # ---------------------------
  p <- ggplot2::ggplot(data=profile) +
    ggplot2::geom_point(ggplot2::aes(x=`m/z`, y=`-log(p)`,
                                     text=`m/z`, 
                                     color=`-log(p)`, 
                                     key=`m/z`)) +
    
    ggplot2::scale_x_discrete(breaks = xaxis, labels=as.character(xaxis)) + 
    ggplot2::scale_colour_gradientn(colours = cf(n)) +
    #ggplot2::scale_y_log10()+
    ggplot2::scale_y_continuous(labels=scaleFUN)
  p
}

ggPlotTT <- function(mSet, cf, n=20){
  profile <- data.table::as.data.table(mSet$analSet$tt$p.log[mSet$analSet$tt$inx.imp],keep.rownames = T)
  
  if(nrow(profile)==0){
    shiny::showNotification("No significant hits")
    return(NULL)
  }
  
  profile[,2] <- round(profile[,2], digits = 2)
  profile$Peak <- c(1:nrow(profile))
  colnames(profile)[1:2] <- c("m/z", "-log(p)")
  scaleFUN <- function(x) sprintf("%.1f", x)
  
  xaxis = seq(0,600, 50)
  # ---------------------------
  p <- ggplot2::ggplot(data=profile) +
    ggplot2::geom_point(ggplot2::aes(x=`m/z`, y=`-log(p)`,
                                     text=`m/z`,
                                     color=`-log(p)`, 
                                     key=`m/z`)) +
    ggplot2::scale_x_discrete(breaks = xaxis, labels=as.character(xaxis)) +
    ggplot2::scale_colour_gradientn(colours = cf(n)) +
    #ggplot2::scale_y_log10()+
    ggplot2::scale_y_continuous(labels=scaleFUN)
  p
}

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
    ggplot2::geom_bar(mapping = aes(x = `m/z`, 
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
                  color="p-value")+
    
    ggplot2::scale_colour_gradientn(colours = cf(n)) +
    ggplot2::scale_fill_gradientn(colours = cf(n)) +
    ggplot2::scale_y_continuous(labels=scaleFUN)
  
  p
}

ggPlotFC <- function(mSet, cf, n=20){
  profile <- data.table::as.data.table(mSet$analSet$fc$fc.log[mSet$analSet$fc$inx.imp],keep.rownames = T)
  
  if(nrow(profile)==0){
    shiny::showNotification("No significant hits")
    return(NULL)
  }
  
  colnames(profile) <- c("m/z", "log2fc")
  profile$Peak <- c(1:nrow(profile))
  scaleFUN <- function(x) sprintf("%.0f", x)
  # ---------------------------
  p <- ggplot2::ggplot(data=profile) +
    ggplot2::geom_point(ggplot2::aes(x=Peak, 
                                     y=log2fc, 
                                     color=log2fc, 
                                     key=`m/z`,
                                     text=`m/z`)) +
    ggplot2::geom_abline(ggplot2::aes(intercept = 0, slope = 0)) +
    ggplot2::scale_y_continuous(labels=scaleFUN) +
    ggplot2::scale_colour_gradientn(colours = cf(n))
  
  p
}

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
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data=dt, ggplot2::aes(x=log2FC,
                                              y=`-log10P`,
                                              color=col,
                                              text=`m/z`,
                                              key=`m/z`)) +
    
    ggplot2::scale_colour_gradientn(colours = cf(n),guide=FALSE) 
  
  p
}

ggPlotClass <- function(mSet,
                        pls.type = "plsda",
                        cf,
                        pcs = 3){
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
    
    facet_grid(~Component) +
    ggplot2::scale_fill_manual(values=cf(pcs))
  p
}

ggPlotPerm <- function(mSet,
                       pls.type = "plsda",
                       cf,
                       pcs = 3){
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
  
  p
}

ggPlotROC <- function(data,
                      attempts = 50,
                      cf,
                      class_type="b"){
  
  mean.auc <- data$m_auc
  perf.long <- data$perf
  
  means.per.comp=perf.long[, lapply(.SD, mean), by = comparison]
  
  ncomp = length(unique(means.per.comp$comparison))
  if(ncomp > 2){
    cols = cf(ncomp)+1
    class_type = "m"
    try({
      shiny::showNotification("Calculating AUCs per comparison...")
    })
    perf.long$comparison <- pbapply::pbsapply(perf.long$comparison,
                                              function(comp){
                                                paste0(comp, " || avg. AUC=", round(means.per.comp[comparison == comp]$AUC_PAIR, digits=3), " ||")
                                              })  
  }else{
    cols = cf(attempts)
    class_type = "b"
  }
  
  p <- ggplot(perf.long, aes(FPR,TPR,
                             key = attempt,
                             text = attempt)) +
    ggplot2::geom_path(alpha=.5,
                       cex=.7,
                       aes(color = if(class_type == "m") comparison else as.factor(attempt),
                           text = if(class_type == "m") comparison else as.factor(attempt),
                           key = if(class_type == "m") comparison else as.factor(attempt),
                           group = attempt)) +
    ggplot2::annotate("text",
                      label = paste0("Average AUC: ",
                                     format(mean.auc,
                                            2,
                                            drop0trailing = TRUE,
                                            digits = 2)),
                      size = 8,
                      x = 0.77,
                      y = 0.03) +
    labs(color = if(class_type == "m") "Comparison" else "Attempt",
         text = if(class_type == "m") "Comparison" else "Attempt",
         key = if(class_type == "m") "Comparison" else "Attempt") +
    
    ggplot2::stat_summary_bin(#alpha=.6,
      aes(FPR, TPR, 
          group = comparison), 
      fun.y=mean, geom="line", 
      cex = 2.3,color="black") +
    ggplot2::scale_color_manual(values = cols) +
    #scale_y_continuous(limits=c(0,1)) +
    ggplot2::coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE) +
    ggplot2::coord_cartesian(xlim = c(.04,.96), ylim = c(.04,.96))
  
  p
}

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
    data.norep <- data[,-3]
    colnames(data.norep)[1] <- "m/z"
    data.ci = Rmisc::group.CI(importance ~ `m/z`, data.norep)
    data.ordered <- data.ci[order(data.ci$importance.mean, decreasing = T),]
  }
  
  data.subset <- data.ordered[1:topn,]    
  data.subset$`m/z` <- reorder(x = data.subset$`m/z`, X = -data.subset$importance.mean)
  
  p <- ggplot(data.subset, aes(x = `m/z`,
                               y = importance.mean,
                               text = `m/z`,
                               key = `m/z`)) +
    ggplot2::geom_bar(stat = "identity",
                      aes(fill = importance.mean)) +
    ggplot2::scale_fill_gradientn(colors=cf(20)) +
    
    ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                   axis.ticks.x=ggplot2::element_blank())+
    labs(x="Top hits (m/z)",y=if(ml_type == "glmnet") "Times included in final model" else "Relative importance (%)")
  
  if(topn <= 15){
    p <- p + ggplot2::geom_text(aes(x=`m/z`, y=importance.mean, label=sapply(`m/z`, function(x){
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
  }
  
  mzdata <- p$data
  mzdata$`m/z` <- gsub(mzdata$`m/z`, pattern = "`|'", replacement="")

  list(mzdata = mzdata, plot = p)
  
}

plotPCAloadings.2d <- function(mSet,
                               cf,
                               pcx,
                               pcy,
                               type = "pca"){
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
           
           x.var <- plsda.table[PC == pcx]$var
           y.var <- plsda.table[PC == pcy]$var
           
           # --- coordinates ---
           df <- mSet$analSet$plsr$loadings
           class(df) <- "matrix"
           colnames(df) <- paste0("PC", 1:ncol(df))
         })
  
  df = as.data.frame(df)
  df$extremity <- apply(df, 1, function(row) max(abs(c(row[[pcx]], 
                                                       row[[pcy]]))))
  
  #df <- df[order(df$extremity, decreasing=T),]#[1:2000,]
  
  scaleFUN <- function(x) sprintf("%.4s", x)
  
  p <- ggplot(df, aes(df[[pcx]], df[[pcy]])) +
    ggplot2::geom_point(aes(color = extremity,
                            size = extremity,
                            text = rownames(df),
                            key = rownames(df)),
                        #pch=21, size = 2, stroke = 2,
                        #fill="white", 
                        alpha=0.7)+
    
    scale_size_area(max_size = 15) +
    ggplot2::scale_x_continuous(labels=scaleFUN,name=gsubfn::fn$paste(if(type != "tsne") "$pcx ($x.var%)" else "t-sne dimension 1")) +
    ggplot2::scale_y_continuous(labels=scaleFUN,name=gsubfn::fn$paste(if(type != "tsne") "$pcy ($y.var%)" else "t-sne dimension 2")) +
    ggplot2::scale_colour_gradientn(colors=cf(20))
  #scale_y_discrete(labels=scaleFUN) +
  #scale_x_discrete(labels=scaleFUN)
  p 
}

plotPCAloadings.3d <- function(mSet,
                               cf,
                               pcx,
                               pcy,
                               pcz,
                               font,
                               type = "pca"){
  
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
           
           x.var <- plsda.table[PC == pcx]$var
           y.var <- plsda.table[PC == pcy]$var
           z.var <- plsda.table[PC == pcz]$var
           
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
                             showscale = FALSE,
                             showlegend = FALSE)
  ) %>%
    add_markers() %>%
    layout(scene = basic_scene)
  p
}


plotPCA.3d <- function(mSet,
                       cols,
                       shape.fac="label",
                       pcx, pcy, pcz,
                       type="pca",font,
                       col.fac = "label",
                       mode="normal",cf){
  
  print(type)
  
  switch(type,
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
    as.factor(mSet$dataSet$covars[, ..shape.fac][[1]])
  }
  
  col.vec <- if(is.null(col.fac)){
    classes
  }else if(shape.fac == "label"){
    classes
  }else{
    as.factor(mSet$dataSet$covars[, ..col.fac][[1]])
  }
  
  plots_facet <- lapply(1:length(df_list), function(i){
    
    df = df_list[[i]]
    
    orig_idx = match(rownames(df), rownames(mSet$dataSet$norm))
    
    plots <- plotly::plot_ly(showlegend = F, 
                             scene = paste0("scene", if(i > 1) i else ""))
    
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
          #type='mesh3d',
          alphahull = 0,
          opacity=0.1,
          hoverinfo="none"
        )
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
      }else{
        adj_plot = plots
      }
      show.orbs <- c(show.orbs, worked)
    }
    
    t <- list(family = font$family)
    
    # --- return ---
    pca_plot <- adj_plot %>% add_trace(
      hoverinfo = 'text',
      text = rownames(df),
      x = df[, pcx],
      y = df[, pcy],
      z = df[, pcz],
      visible = rep(T, times=fac.lvls),
      type = "scatter3d",
      color = classes,
      colors = cols,
      opacity = 1,
      marker = list(
        line = list(#color = col.vec[orig_idx],
          #colors = cols,
          width = 2),
        symbol = symbol.vec[orig_idx],
        symbols = c('circle',
                    'diamond',
                    'square',
                    'x',
                    'o'))
    ) 
    # --- return ---
    pca_plot
  })
  
  basic_scene = list(
    aspectmode="cube",
    aspectratio=list(x=1,y=1,z=1),
    camera = list(
      eye = list(x=0, y=0, z= 2)
    ),
    xaxis = list(
      titlefont = list(size = font$ax.txt.size * 1.5),
      title = if(type != "tsne") gsubfn::fn$paste("$pcx ($x.var%)") else "t-sne dimension 1"),
    yaxis = list(
      titlefont = list(size = font$ax.txt.size * 1.5),
      title = if(type != "tsne") gsubfn::fn$paste("$pcy ($y.var%)") else "t-sne dimension 2"),
    zaxis = list(
      titlefont = list(size = font$ax.txt.size * 1.5),
      title = if(type != "tsne") gsubfn::fn$paste("$pcz ($z.var%)") else "t-sne dimension 3"))
  
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



plotPCA.2d <- function(mSet, shape.fac = "label", cols, col.fac = "label",
                       pcx, pcy, mode="normal", type="pca",
                       cf = rainbow){
  
  classes <- if(mode == "ipca"){
    mSet$dataSet$facA
  }else{
    mSet$dataSet$cls
  }
  
  symbols = c("16",#'circle',
              "18",#'diamond',
              "15",#'square',
              "4",#'x',
              "1"#'o'
  )
  
  switch(type,
         tsne = {
           df <- mSet$analSet$tsne$x
           x.var <- ""
           y.var <- ""
           fac.lvls <- length(levels(mSet$dataSet$cls))
           
           xc=mSet$analSet$tsne$x[, pcx]
           yc=mSet$analSet$tsne$x[, pcy]
           
           dat_long <- data.table(variable = rownames(mSet$dataSet$norm),
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
  
  if(mode == "ipca"){
    fac.lvls <- length(levels(mSet$dataSet$facA))
    dat_long$groupB <- mSet$dataSet$facB
  }
  
  # - - - - - - - - -
  dat_long$color <- if(is.null(col.fac)){
    factor(1) # all same shape...
  } else if(col.fac == "label"){
    dat_long$group
  }else{
    as.factor(mSet$dataSet$covars[,..col.fac][[1]])
  }
  
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
  
  p <- ggplot2::ggplot(dat_long, aes(x, y,group=group)) +
    ggplot2::geom_point(size=5, aes(shape=shape,
                                    text=variable,
                                    fill=group,
                                    color=color), alpha=0.7)+
    ggplot2::stat_ellipse(geom = "polygon", aes(fill=group), alpha = 0.3,level = .95) +
    ggplot2::scale_x_continuous(name=gsubfn::fn$paste(if(type != "tsne") "$pcx ($x.var%)" else "t-sne dimension 1")) +
    ggplot2::scale_y_continuous(name=gsubfn::fn$paste(if(type != "tsne") "$pcy ($y.var%)" else "t-sne dimension 2")) +
    ggplot2::scale_fill_manual(values=cols) +
    ggplot2::scale_color_manual(values = cols)
  
  if(mode == "ipca"){
    p <- p + facet_wrap(~groupB,ncol = 2)
    p <- p + ggplot2::ggtitle(Hmisc::capitalize(mSet$dataSet$facB.lbl))
  }
  p
}

ggPlotVenn <- function(mSet,
                       venn_yes,
                       top = 100,
                       cols,
                       cf){
  
  flattened <- getTopHits (mSet, unlist(venn_yes$now), top)
  
  parseFun = function(labels){
    sapply(labels, function(label){
      if(grepl(label, pattern=":")){
        split = trimws(stringr::str_split(label, ":")[[1]])
        stats_on = toupper(split[1])
        which_stats = split[length(split)]
        subset = grep("=", split, value=T)
        print(subset)
        paste0(stats_on, "\n", if(length(subset)>0) paste0(subset, collapse="\n", "\n") else "", paste0(">> ", which_stats, " <<\n\n"))
      }else{
        label
      }
    })
  }
  
  #labels = c("species: aov", "species: PC1 (PLS-DA)", "species_type:tissue=distal intestine: tt")
  #cat(parseFun(labels))
  
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

ggPlotScree <- function(mSet, cf, pcs=20){
  df <- data.table::data.table(
    pc = 1:length(names(mSet$analSet$pca$variance)),
    var = round(mSet$analSet$pca$variance*100,digits = 1))
  p <- ggplot2::ggplot(data=df[1:20,]) + ggplot2::geom_line(mapping = aes(x=pc, y=var, colour=var), cex=3) +
    
    ggplot2::scale_colour_gradientn(colours = cf(20))
  # - - - - -
  p
}

ggPlotWordBar <- function(wcdata, cf, plotlyfy=T){
  g <- ggplot(wcdata, aes(y = freq, x = reorder(word, 
                                                freq, 
                                                sum)))
  g <- g + ggplot2::geom_bar(aes(fill = freq),
                             stat = "identity") +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_gradientn(colors=cf(256)) +
    ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                   axis.ticks.x=ggplot2::element_blank()) +
    labs(x="Word",y="Frequency")
  g
}

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
    p <- ggplot(data, aes(samples,power)) +
      ggplot2::geom_path(alpha=.5,
                         cex=.5,
                         aes(color = comparison, group = comparison)) +
      ggplot2::stat_summary_bin(#alpha=.6,
        aes(samples, 
            power, 
            group=comparison), 
        fun.y=mean, geom="line", 
        cex = 2.3,color="black") +
      ggplot2::stat_summary_bin(#alpha=.6,
        aes(samples, power, 
            color=comparison
            #,group=comparison
            ), 
        fun.y=mean, geom="line", 
        cex = 1.2) +
      ggplot2::stat_summary_bin(aes(samples,
                                    power), 
                                fun.y=mean, color="black", 
                                geom="line", cex = 2) +
      ggplot2::coord_cartesian(xlim = c(0,max_samples), ylim = c(.04,.96))
    
    p 
  }
}

ggPlotMummi <- function(mum_mSet, anal.type = "mummichog", cf){
  if (anal.type == "mummichog") {
    mummi.mat <- mum_mSet$mummi.resmat
    y <- -log10(mummi.mat[, 5])
    x <- mummi.mat[, 3]/mummi.mat[, 4]
    pathnames <- rownames(mummi.mat)
  }
  else {
    gsea.mat <- mum_mSet$mummi.gsea.resmat
    print(colnames(gsea.mat))
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
  dat = as.data.frame(cbind(pathway = rownames(mummi.mat),
                            y = as.numeric(y), 
                            x = as.numeric(x), 
                            as.numeric(radi.vec)))
  
  scaleFUN <- function(x) sprintf("%.4s", x)
  
  p = ggplot(dat) + geom_point(aes(y = y, 
                                   x = x, 
                                   size = `radi.vec`, 
                                   color = `radi.vec`,
                                   text = pathway)) +
     ggtitle("Enrichment Results") +
    ggplot2::ylab(switch(anal.type, 
                         mummichog = "-log10(FET)", 
                         gsea = "-log10(p.value)")) + 
    ggplot2::xlab("significant hits / all pathway hits") +
    ggplot2::scale_colour_gradientn(colours = cf(20)) +
    scale_y_discrete(labels=scaleFUN) +
    scale_x_discrete(labels=scaleFUN)
  
  p
}
