#' @export
ggplotMeba <- function(mz, draw.average, cols=c("Red", "Green")){
  cols <- if(is.null(cols)) c("Red", "Green") else(cols)
  print(cols)
  profile <- getProfile(mz, mode="time")
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

#' @export
ggplotSummary <- function(mz = curr_mz, cols=c("Red", "Green")){
  cols <- if(is.null(cols)) c("Red", "Green") else(cols)
  print(cols)
  if(substring(dataSet$format,4,5)!="ts"){
    # --- ggplot ---
    profile <- getProfile(mz, mode="stat")
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
  }else if(dataSet$design.type =="time"){ # time trend within phenotype
    profile <- getProfile(mz, mode="time")
    # -----------
    print(profile)
    # ggplot
    plot <- ggplot(data=profile, aes(x=Time, y=Abundance, fill=Group, color=Group)) +
      geom_boxplot(alpha=0.4, aes(group=Group)) +
      geom_point(aes(text=Sample),alpha=0.4, size = 2, shape = 1, position = position_dodge(width=0.1)) +
      theme_minimal(base_size = 10) +
      scale_fill_manual(values=cols) +
      scale_color_manual(values=cols) + facet_wrap(~Group)
    # ---------------
    ggplotly(plot, tooltip="Sample")
  }
  #
}

ggPlotTT <- function(cf, n){
  profile <- as.data.table(analSet$tt$p.log[analSet$tt$inx.imp],keep.rownames = T)
  colnames(profile) <- c("mz", "p")
  profile$Peak <- c(1:nrow(profile)) 
  # ---------------------------
  plot <- ggplot(data=profile) +
    geom_point(aes(x=Peak, y=p,text=mz, color=p, key=mz)) +
    theme_minimal(base_size = 10) +
    scale_colour_gradientn(colours = cf(n)) +
    scale_y_log10()
  ggplotly(plot, tooltip="mz")
}

ggPlotFC <- function(cf, n){
  profile <- as.data.table(analSet$fc$fc.log[analSet$fc$inx.imp],keep.rownames = T)
  profile
  colnames(profile) <- c("mz", "log2fc")
  profile$Peak <- c(1:nrow(profile)) 
  # ---------------------------
  plot <- ggplot(data=profile) +
    geom_point(aes(x=Peak, y=log2fc, text=log2fc, color=log2fc, key=mz)) +
    geom_abline(aes(intercept = 0, slope = 0)) +
    theme_minimal(base_size = 10) +
    scale_colour_gradientn(colours = cf(n))
  ggplotly(plot, tooltip="log2fc")
}

ggPlotVolc <- function(cf, n){
    vcn<-analSet$volcano;
    dt <- as.data.table(cbind(vcn$fc.log, vcn$p.log), keep.rownames=T)
    imp.inx<-(vcn$inx.up | vcn$inx.down) & vcn$inx.p;
    colnames(dt) <- c("mz", "log2FC", "minlog10P")
    plot <- ggplot() +
      #geom_point(data=dt[!imp.inx], aes(x=log2FC, y=minlog10P)) +
      geom_point(data=dt[imp.inx], aes(x=log2FC, y=minlog10P,text=mz,color=abs(log2FC*minlog10P), key=mz)) +
      theme_minimal(base_size = 10) +
      scale_colour_gradientn(colours = cf(n),guide=FALSE)
    ggplotly(plot, tooltip="mz")
    #abline (v = vcn$max.xthresh, lty=3);
    #abline (v = vcn$min.xthresh, lty=3);
    #abline (h = vcn$thresh.y, lty=3);
}


