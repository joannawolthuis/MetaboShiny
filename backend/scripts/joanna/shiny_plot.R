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
    if(vcn$paired){
      xlim<-c(-nrow(dataSet$norm)/2, nrow(dataSet$norm)/2)*1.2;
      
      # merge fc.all two rows into one, bigger one win
      fc.all <- apply(vcn$fc.all, 2, function(x){ if(x[1] > x[2]){return(x[1])}else{return(-x[2])}})
      hit.inx <- vcn$inx.p & (vcn$inx.up | vcn$inx.down);
      print(cbind(fc.all, vcm$p.log))
      plot(fc.all, vcn$p.log, xlim=xlim, pch=20, cex=ifelse(hit.inx, 1.2, 0.8),
           col = ifelse(hit.inx, "red", "black"),
           xlab="Count of Significant Pairs", ylab="-log10(p)");
      
    }else{
      imp.inx<-(vcn$inx.up | vcn$inx.down) & vcn$inx.p;
      print(head(cbind(vcn$fc.log, vcn$p.log)))
      
      plot(vcn$fc.log, vcn$p.log, pch=20, cex=ifelse(imp.inx, 1.2, 0.7),
           col = ifelse(imp.inx, "red", "black"),
           xlab="log2 (FC)", ylab="-log10(p)");
    }
    
    abline (v = vcn$max.xthresh, lty=3);
    abline (v = vcn$min.xthresh, lty=3);
    abline (h = vcn$thresh.y, lty=3);
    axis(4); # added by Beomsoo
}


# library(manhattanly)
# testy <- as.data.table(analSet$volcano$sig.mat, keep.rownames = T)
# testy2 <- as.data.table(cbind(analSet$volcano$fc.log, analSet$volcano$p.log), keep.rownames=T)
# testy2
# 
# colnames(testy) <- c("m/z", "a", "EFFECTSIZE", "P", "b")
# colnames(testy2) <- c("m/z", "EFFECTSIZE", "log10p")
# 
# volcanorObject <- volcanor(testy2, snp = "m/z")
# volcanoly(testy, varName = "m/z",  genomewideline = 1, effect_size_line = FALSE) %>% layout(title = "", xaxis=list(title="Log10 fold change"))
# 
# ggPlotVolc(rainbow, 10)
