#' @export
ggplotMeba <- function(mz, draw.average, cols=c("Red", "Green")){
  cols <- if(is.null(cols)) c("Red", "Green") else(cols)
  print(cols)
  profile <- getProfile(mz)
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
ggplotSummary <- function(mz, cols=c("Red", "Green")){
  cols <- if(is.null(cols)) c("Red", "Green") else(cols)
  print(cols)
  if(substring(dataSet$format,4,5)!="ts"){
    
    par(mar=c(4,4,2,2), mfrow = c(1,2), oma=c(0,0,2,0));
    
    mns <- by(as.numeric(dataSet$proc[, mz]), dataSet$proc.cls, mean, na.rm=T);
    sds <- by(as.numeric(dataSet$proc[, mz]), dataSet$proc.cls, sd, na.rm=T);
    
    ups <- mns + sds;
    dns <- mns - sds;
    
    # all concentration need start from 0
    y <- c(0, dns, mns, ups);
    
    rg <- range(y) + 0.05 * diff(range(y)) * c(-1, 1)
    pt <- pretty(y)
    
    axp=c(min(pt), max(pt[pt <= max(rg)]),length(pt[pt <= max(rg)]) - 1);
    
    # ymk <- pretty(c(0,ymax));
    x <- barplot(mns, col= color.vec(), las=2, yaxp=axp, ylim=range(pt));
    arrows(x, dns, x, ups, code=3, angle=90, length=.1);
    axis(1, at=x, col="white", col.tick="black", labels=F);
    box();
    mtext("Original Conc.", line=1);
    
    boxplot(dataSet$norm[, mz]~dataSet$cls,las=2, col=color.vec());
    mtext("Normalized Conc.", line=1);
    title(main=cmpdNm, out=T);
    #
  }else if(dataSet$design.type =="time0"){
    #
    plotProfile(mz);
    #
  }else{
    if(dataSet$design.type =="time"){ # time trend within phenotype
      profile <- getProfile(mz)
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
}