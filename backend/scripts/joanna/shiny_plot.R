

#' @export
ggplotNormSummary <- function(mSet){
  # 4 by 4 plot, based on random 20-30 picked 
  orig_data <- mSet$dataSet$orig
  norm_data <- mSet$dataSet$norm
  
  which_cpds <- sample(colnames(orig_data[4:ncol(orig_data)]), 20, replace = FALSE, prob = NULL)
  which_samps <- sample(rownames(orig_data), 20, replace = FALSE, prob = NULL)
  
  orig_melt <- melt(orig_data[, which_cpds])
  norm_melt <- melt(norm_data[, which_cpds])

  plot <- ggplot(data=orig_melt) +
    theme_minimal(base_size = 10) #+ facet_grid(. ~ variable)
  
  RES1 <- plot + geom_density(aes(x=value), colour="blue", fill="blue", alpha=0.4)

  RES2 <- plot + geom_boxplot(alpha=0.4,
                              aes(x=value,y=variable),
                              color=rainbow(20), 
                              alpha=0.4) + geom_vline(aes(xintercept=median(orig_melt$value)))
  
  plot <- ggplot(data=norm_melt) +
    theme_minimal(base_size = 10) #+ facet_grid(. ~ variable)
  
  RES3 <- plot + geom_density(aes(x=value), colour="pink", fill="pink", alpha=0.4)
  
  RES4 <- plot + geom_boxplot(alpha=0.4,
                      aes(x=value,y=variable),
                      color=rainbow(20), 
                      alpha=0.4) + geom_vline(aes(xintercept=median(norm_melt$value)))
  
  list(tl=RES1, bl=RES2, tr=RES3, br=RES4)
  
} 

ggplotSampleNormSummary <- function(mSet){
  # 4 by 4 plot, based on random 20-30 picked 
  orig_data <- mSet$dataSet$orig
  norm_data <- mSet$dataSet$norm
  
  which_samps <- sample(rownames(orig_data), 
                        20, 
                        replace = FALSE, 
                        prob = NULL)
  
  sumsOrig <- rowSums(orig_data[which_samps,])
  sumsNorm <- rowSums(norm_data[which_samps,])

  orig_data$Label <- rownames(orig_data)
  orig_melt <- melt(orig_data[which_samps,],
                    id.vars = "Label")
  orig_melt_sums <- melt(sumsOrig)
  orig_melt_sums$variable <- rownames(orig_melt_sums)
  
  norm_data$Label <- rownames(norm_data)
  norm_melt <- melt(norm_data[which_samps,],
                    id.vars="Label")
  norm_melt_sums <- melt(sumsNorm)
  norm_melt_sums$variable <- rownames(norm_melt_sums)
  
  RES1 <- ggplot(data=orig_melt_sums) +
    theme_minimal(base_size = 10) + geom_density(aes(x=value), colour="blue", fill="blue", alpha=0.4)
  
  RES2 <- ggplot(data=orig_melt) +
    theme_minimal(base_size = 10) + geom_boxplot(alpha=0.4,
                      aes(x=value,y=Label),
                      color=rainbow(20), 
                      alpha=0.4) + geom_vline(aes(xintercept=median(orig_melt$value),text=Label))

  RES3 <- ggplot(data=norm_melt_sums) +
    theme_minimal(base_size = 10) + geom_density(aes(x=value), colour="pink", fill="pink", alpha=0.4)
  
  RES4 <- ggplot(data=norm_melt) +
    theme_minimal(base_size = 10) + geom_boxplot(alpha=0.4,
                      aes(x=value,y=Label),
                      color=rainbow(20), 
                      alpha=0.4) + geom_vline(aes(xintercept=median(norm_melt$value),text=Label))
   
  list(tl=RES1, bl=RES2, tr=RES3, br=RES4)
} 


#' @export
ggplotMeba <- function(cpd, draw.average=T, cols=c("Red", "Green")){
  cols <- if(is.null(cols))c("Blue","Pink") else(cols)
  profile <- getProfile(cpd, mode="time")
  plot <- if(draw.average){
    ggplot(data=profile) +
      geom_line(size=0.3, aes(x=Time, y=Abundance, group=Sample, color=Group, text=Sample), alpha=0.4) +
      stat_summary(fun.y="mean", size=1.5, geom="line", aes(x=Time, y=Abundance, color=Group, group=Group)) +
      scale_x_discrete(expand = c(0, 0)) +
      theme_minimal(base_size = 10) +
      scale_color_manual(values=cols)
  } else{
    ggplot(data=profile) +
      geom_line(size=0.7, aes(x=Time, y=Abundance, group=Sample, color=Group, text=Sample)) +
      scale_x_discrete(expand = c(0, 0)) +
      theme_minimal(base_size = 10) +
      scale_color_manual(values=cols)
  }
  # ---------------
  ggplotly(plot, tooltip="Sample")
}

blackwhite.colors <- function(n){
  gray.colors(n, start=0, end=1)
}

#' @export
ggplotSummary <- function(cpd = curr_cpd, cols=c("Blue", "Pink")){
  cols <- if(is.null(cols))c("Blue","Pink") else(cols)
  if(substring(mSet$dataSet$format,4,5)!="ts"){
    # --- ggplot ---
    profile <- getProfile(cpd, mode="stat")
    # -----------
    # ggplot
    plot <- ggplot(data=profile, aes(x=Group,y=Abundance, fill=Group, color=Group)) +
      geom_boxplot(alpha=0.4) +
      geom_point(aes(text=Sample),alpha=0.4, size = 2, shape = 1, position = position_dodge(width=0.1)) +
      theme_minimal(base_size = 10) +
      scale_fill_manual(values=cols) +
      scale_color_manual(values=cols)
    # ---------------
    ggplotly(plot, tooltip="Sample")
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
    plot <- ggplot(data=profile, aes(x=Time, y=Abundance, group=interaction(Time, Group),fill=Group, color=Group)) +
      geom_boxplot(alpha=0.4) +
      geom_point(aes(text=Sample),alpha=0.4, size = 2, shape = 1, position = position_dodge(width=0.1)) +
      theme_minimal(base_size = 10) +
      scale_fill_manual(values=cols) +
      scale_color_manual(values=cols) 
    # ---------------
    ggplotly(plot, tooltip="Sample") %>% layout(boxmode = "group")
  }
  #
}

ggPlotTT <- function(cf, n){
  profile <- as.data.table(mSet$analSet$tt$p.log[mSet$analSet$tt$inx.imp],keep.rownames = T)
  colnames(profile) <- c("cpd", "p")
  profile$Peak <- c(1:nrow(profile)) 
  # ---------------------------
  plot <- ggplot(data=profile) +
    geom_point(aes(x=Peak, y=p,text=cpd, color=p, key=cpd)) +
    theme_minimal(base_size = 10) +
    scale_colour_gradientn(colours = cf(n)) +
    scale_y_log10()
  ggplotly(plot, tooltip="cpd")
}

ggPlotFC <- function(cf, n){
  profile <- as.data.table(mSet$analSet$fc$fc.log[mSet$analSet$fc$inx.imp],keep.rownames = T)
  profile
  colnames(profile) <- c("cpd", "log2fc")
  profile$Peak <- c(1:nrow(profile)) 
  # ---------------------------
  plot <- ggplot(data=profile) +
    geom_point(aes(x=Peak, y=log2fc, text=log2fc, color=log2fc, key=cpd)) +
    geom_abline(aes(intercept = 0, slope = 0)) +
    theme_minimal(base_size = 10) +
    scale_colour_gradientn(colours = cf(n))
  ggplotly(plot, tooltip="log2fc")
}

ggPlotVolc <- function(cf, n){
    vcn<-mSet$analSet$volcano;
    dt <- as.data.table(cbind(vcn$fc.log, vcn$p.log), keep.rownames=T)
    imp.inx<-(vcn$inx.up | vcn$inx.down) & vcn$inx.p;
    colnames(dt) <- c("cpd", "log2FC", "minlog10P")
    plot <- ggplot() +
      #geom_point(data=dt[!imp.inx], aes(x=log2FC, y=minlog10P)) +
      geom_point(data=dt[imp.inx], aes(x=log2FC, 
                                       y=minlog10P,
                                       text=cpd,
                                       color=abs(log2FC*minlog10P), 
                                       key=cpd)) +
      theme_minimal(base_size = 10) +
      scale_colour_gradientn(colours = cf(n),guide=FALSE)
    ggplotly(plot, tooltip="cpd")
    #abline (v = vcn$max.xthresh, lty=3);
    #abline (v = vcn$min.xthresh, lty=3);
    #abline (h = vcn$thresh.y, lty=3);
}


