#' @export
ggplotNormSummary <- function(mSet){
  # 4 by 4 plot, based on random 20-30 picked 
  orig_data <- mSet$dataSet$procr
  norm_data <- mSet$dataSet$norm

  sampsize = if(nrow(norm_data) > 20) 20 else nrow(norm_data)
  which_cpds <- sample(colnames(norm_data), sampsize, replace = FALSE, prob = NULL)
  which_samps <- sample(rownames(norm_data), sampsize, replace = FALSE, prob = NULL)
  
  orig_melt <- reshape2::melt(orig_data[which_samps, which_cpds])
  orig_melt[is.na(orig_melt)] <- 0
  
  norm_melt <- reshape2::melt(norm_data[which_samps, which_cpds])

  plot <- ggplot2::ggplot(data=orig_melt) +
    plot.theme(base_size = 15) #+ facet_grid(. ~ variable)
  
  RES1 <- plot + ggplot2::geom_density(ggplot2::aes(x=value), colour="blue", fill="blue", alpha=0.4)

  RES2 <- plot + ggplot2::geom_boxplot(alpha=0.4,
                              ggplot2::aes(x=value,y=variable),
                              color=rainbow(sampsize), 
                              alpha=0.4) + ggplot2::geom_vline(ggplot2::aes(xintercept=median(orig_melt$value)))
  
  plot <- ggplot2::ggplot(data=norm_melt) +
    plot.theme(base_size = 15) #+ facet_grid(. ~ variable)
  
  RES3 <- plot + ggplot2::geom_density(ggplot2::aes(x=value), colour="pink", fill="pink", alpha=0.4)
  
  RES4 <- plot + ggplot2::geom_boxplot(alpha=0.4,
                      ggplot2::aes(x=value,y=variable),
                      color=rainbow(sampsize), 
                      alpha=0.4) + ggplot2::geom_vline(ggplot2::aes(xintercept=median(norm_melt$value)))
  
  list(tl=RES1, bl=RES2, tr=RES3, br=RES4)
} 

ggplotSampleNormSummary <- function(mSet){
  # 4 by 4 plot, based on random 20-30 picked 
  orig_data <- mSet$dataSet$orig
  norm_data <- mSet$dataSet$norm
  
  sampsize = if(nrow(norm_data) > 20) 20 else nrow(norm_data)
  
  which_samps <- sample(rownames(orig_data), 
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
    plot.theme(base_size = 15) + ggplot2::geom_density(ggplot2::aes(x=value), colour="blue", fill="blue", alpha=0.4)
  
  RES2 <- ggplot2::ggplot(data=orig_melt) +
    plot.theme(base_size = 15) + ggplot2::geom_boxplot(alpha=0.4,
                      ggplot2::aes(x=value,y=Label),
                      color=rainbow(sampsize), 
                      alpha=0.4) + ggplot2::geom_vline(ggplot2::aes(xintercept=median(orig_melt$value),text=Label))

  RES3 <- ggplot2::ggplot(data=norm_melt_sums) +
    plot.theme(base_size = 15) + ggplot2::geom_density(ggplot2::aes(x=value), colour="pink", fill="pink", alpha=0.4)
  
  RES4 <- ggplot2::ggplot(data=norm_melt) +
    plot.theme(base_size = 15) + ggplot2::geom_boxplot(alpha=0.4,
                      ggplot2::aes(x=value,y=Label),
                      color=rainbow(sampsize), 
                      alpha=0.4) + ggplot2::geom_vline(ggplot2::aes(xintercept=median(norm_melt$value),text=Label))
   
  list(tl=RES1, bl=RES2, tr=RES3, br=RES4)
} 


#' @export
ggplotMeba <- function(cpd, draw.average=T, cols=NULL, cf){
  cols <- if(is.null(cols)) cf(length(levels(mSet$dataSet$cls))) else(cols)
  profile <- getProfile(cpd, mode="time")
  plot <- if(draw.average){
    ggplot2::ggplot(data=profile) +
      ggplot2::geom_line(size=0.3, ggplot2::aes(x=Time, y=Abundance, group=Sample, color=Group, text=Sample), alpha=0.4) +
      stat_summary(fun.y="mean", size=1.5, geom="line", ggplot2::aes(x=Time, y=Abundance, color=Group, group=Group)) +
      ggplot2::scale_x_discrete(expand = c(0, 0)) +
      plot.theme(base_size = 15) +
      ggplot2::scale_color_manual(values=cols)
  } else{
    ggplot2::ggplot(data=profile) +
      ggplot2::geom_line(size=0.7, ggplot2::aes(x=Time, y=Abundance, group=Sample, color=Group, text=Sample)) +
      ggplot2::scale_x_discrete(expand = c(0, 0)) +
      plot.theme(base_size = 15) +
      ggplot2::scale_color_manual(values=cols)
  }
  # ---------------
  plotly::ggplotly(plot, tooltip="Sample")
}

blackwhite.colors <- function(n){
  gray.colors(n, start=0, end=1)
}

#' @export
ggplotSummary <- function(cpd = curr_cpd, cols=c("black", "pink"), sourceTable = mSet$dataSet$norm, cf=rainbow){
  cols <- if(is.null(cols)) cf(length(levels(mSet$dataSet$cls))) else(cols)
  if(substring(mSet$dataSet$format,4,5)!="ts"){
    # --- ggplot ---
    profile <- getProfile(cpd, mode="stat", sourceTable = sourceTable)
    df_line <- data.table(x = c(1,2),
                          y = rep(min(profile$Abundance - 0.1),2))
    stars = ""
    try({
      pval <- mSet$analSet$tt$sig.mat[which(rownames(mSet$analSet$tt$sig.mat) == curr_cpd), "p.value"]
      if(length(pval) == 0){
        stars <- ""
      }else{
        if(pval > 0.05) stars <- "n.s."
        else if(pval < 0.05 & pval > 0.01) stars <- "*"
        else if(pval < 0.01 & pval > 0.001) stars <- "***"
        else stars <- "****"
      } 
    })
    # -----------
    # ggplot
    plot <- ggplot2::ggplot(data=profile, ggplot2::aes(x=Group,y=Abundance, fill=Group, color=Group)) +
      ggplot2::geom_boxplot(alpha=0.4) +
      ggplot2::geom_point(ggplot2::aes(text=Sample),alpha=0.4, size = 2, shape = 1, position = position_dodge(width=0.1)) +
      plot.theme(base_size = 15) +
      ggplot2::scale_fill_manual(values=cols) +
      ggplot2::scale_color_manual(values=cols)
    plot <- plot + ggplot2::annotate("text", x = 1.5, y = min(profile$Abundance - 0.3), label = stars, size = 8, col = "black")
    # ---------------
    plotly::ggplotly(plot, tooltip="Sample")
    #
  }else if(mSet$dataSet$design.type =="time"){ # time trend within phenotype
    cpd = curr_cpd
    profile <- getProfile(cpd, mode="time")
    # -----------
    print(profile)
    # ggplot
    # p <- plot_ly(profile, 
    #              x = ~Time, 
    #              y = ~Abundance, 
    #              color = ~Group, 
    #              type = "box", 
    #              colors=cols, 
    #              boxpoints = 'all'
    #              ) %>%
    #   layout(boxmode = "group")
    # p
    plot <- ggplot2::ggplot(data=profile, ggplot2::aes(x=Time, y=Abundance, group=interaction(Time, Group),fill=Group, color=Group)) +
      ggplot2::geom_boxplot(alpha=0.4) +
      ggplot2::geom_point(ggplot2::aes(text=Sample),alpha=0.4, size = 2, shape = 1, position = position_dodge(width=0.1)) +
      plot.theme(base_size = 15) +
      ggplot2::scale_fill_manual(values=cols) +
      ggplot2::scale_color_manual(values=cols)+
      geom_text(x = 1.5, y = 1, label = "***")
    # ---------------
    plotly::ggplotly(plot, tooltip="Sample") %>% layout(boxmode = "group")
  }
  #
}

ggPlotTT <- function(cf, n){
  profile <- as.data.table(mSet$analSet$tt$p.log[mSet$analSet$tt$inx.imp],keep.rownames = T)
  colnames(profile) <- c("cpd", "p")
  profile$Peak <- c(1:nrow(profile)) 
  # ---------------------------
  plot <- ggplot2::ggplot(data=profile) +
    ggplot2::geom_point(ggplot2::aes(x=Peak, y=p,text=cpd, color=p, key=cpd)) +
    plot.theme(base_size = 15) +
    theme(axis.text=element_text(size=12)) +
    ggplot2::scale_colour_gradientn(colours = cf(n)) +
    ggplot2::scale_y_log10()
  plotly::ggplotly(plot, tooltip="cpd")
}

ggPlotFC <- function(cf, n){
  profile <- as.data.table(mSet$analSet$fc$fc.log[mSet$analSet$fc$inx.imp],keep.rownames = T)
  profile
  colnames(profile) <- c("cpd", "log2fc")
  profile$Peak <- c(1:nrow(profile)) 
  # ---------------------------
  plot <- ggplot2::ggplot(data=profile) +
    ggplot2::geom_point(ggplot2::aes(x=Peak, y=log2fc, text=log2fc, color=log2fc, key=cpd)) +
    ggplot2::geom_abline(ggplot2::aes(intercept = 0, slope = 0)) +
    plot.theme(base_size = 15) +
    ggplot2::scale_colour_gradientn(colours = cf(n))
  plotly::ggplotly(plot, tooltip="log2fc")
}

ggPlotVolc <- function(cf=rainbow, n=256){
    vcn<-mSet$analSet$volcano;
    dt <- as.data.table(vcn$sig.mat[,c(2,4)],keep.rownames = T)
    colnames(dt) <- c("cpd", "log2FC", "-log10P")
    plot <- ggplot2::ggplot() +
      #ggplot2::geom_point(data=dt[!imp.inx], ggplot2::aes(x=log2FC, y=minlog10P)) +
      ggplot2::geom_point(data=dt, ggplot2::aes(x=log2FC, 
                                       y=`-log10P`,
                                       text=cpd,
                                       color=abs(log2FC*`-log10P`), 
                                       key=cpd)) +
      plot.theme(base_size = 15) +
      ggplot2::scale_colour_gradientn(colours = cf(n),guide=FALSE)
    plotly::ggplotly(plot, tooltip="cpd")
}

ggPlotPCApairs <- function(cols = c("black", "pink"), pc.num, type = "pca", cf = rainbow){
  
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
  
    p
    
    plotly::ggplotly(p)
  # }
  # else {
  #   pairs(mSet$analSet$pca$x[, 2:pc.num], labels = pclabels)
  # }
}

ggPlotClass <- function(pls.type = "plsda", cf = "rainbow"){
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
    facet_grid(~Component) + 
    scale_fill_manual(values=cf(3))
  return(p)
  }

ggPlotPerm<- function(pls.type = "plsda", cf = "rainbow"){
  bw.vec <- mSet$analSet$plsda$permut
  len <- length(bw.vec)
  df <- melt(bw.vec)
  colnames(df) = "acc"
  p <- ggplot(df) +
    geom_histogram(mapping=aes(x=acc, y=..count.., fill=factor(..count..)),
                   binwidth=0.01) +
    scale_fill_manual(values=cf(20)) +
    theme_minimal() + 
    plot.theme(base_size = 15) + 
    theme(legend.position="none",legend.title.align = 0.5)+
    labs(title="PLSDA predictive value",x="Accuracy", y = "Permutations") +
    geom_segment(data=df,
                 color="black",
                 x=bw.vec[1],
                 xend=bw.vec[1],
                 y=0,aes(yend=.1*nrow(df)),
                 size=1.5,
                 linetype=8)+
    annotate("text",x=bw.vec[1],y=.11*nrow(df),label = mSet$analSet$plsda$permut.p, color="black",size=4)
 return(p)
}
