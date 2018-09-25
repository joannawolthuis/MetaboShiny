# GET DATA

library(data.table)
library(gsubfn)
library(plotly)
library(ggplot2)

# SETTINGS:
# norm: sum, fillmiss = rf, missperc = 1%, rest default
# chicken_mset <- mSet
# swine_mset <- mSet

t <- list(
  family = "Trebuchet MS",
  size = 40,
  color = toRGB("black"))

#save(chicken_mset, swine_mset, file="~/data_twopager.RData")

load("~/data_twopager.RData")

# PLS-DA BICOLOURED PLOT
mSet <- chicken_mset

msets <- list(chicken_mset, swine_mset)
# get table
plsda.tables <- lapply(msets, function(mSet){
  plsda.table <- as.data.table(round(mSet$analSet$plsr$Xvar 
                                     / mSet$analSet$plsr$Xtotvar 
                                     * 100.0,
                                     digits = 2),
                               keep.rownames = T)
  colnames(plsda.table) <- c("Principal Component", "% variance")
  plsda.table[, "Principal Component"] <- paste0("PC", 1:nrow(plsda.table))
  # - - - - -
  plsda.table
})

coord.tables <- lapply(list(chicken_mset, swine_mset), function(mSet){
  coords <- as.data.frame.array(mSet$analSet$plsr$scores)
  colnames(coords) <- paste0("PC", 1:ncol(coords))
  # - - - - - 
  as.data.frame(coords)
})

joined.plsda.tables <- rbindlist(plsda.tables)
concat.plsda.tables <- aggregate(`% variance` ~ `Principal Component`, joined.plsda.tables, paste, collapse = " / ")

joined.coord.tables <- rbindlist(coord.tables)

joined.cls <- unlist(lapply(msets, function(mSet){
  cls <- mSet$dataSet$cls
}))


# --- coordinates ---

# --- vip table ---
# colnames(mSet$analSet$plsda$vip.mat) <- paste0("PC", 1:ncol(mSet$analSet$plsda$vip.mat))
# compounds_pc <- as.data.table(mSet$analSet$plsda$vip.mat,keep.rownames = T)
# ordered_pc <- setorderv(compounds_pc, input$plsda_vip_cmp, -1)
# plsda_tab <<- cbind(ordered_pc[1:50, c("rn")], 
#                     Rank=c(1:50))

mSet <- swine_mset

plsda.table <- as.data.table(round(mSet$analSet$plsr$Xvar 
                                   / mSet$analSet$plsr$Xtotvar 
                                   * 100.0,
                                   digits = 2),
                             keep.rownames = T)
colnames(plsda.table) <- c("Principal Component", "% variance")
plsda.table[, "Principal Component"] <- paste0("PC", 1:nrow(plsda.table))

coords <- as.data.frame.array(mSet$analSet$plsr$scores)
colnames(coords) <- paste0("PC", 1:ncol(coords))

# PLOT


x <- "PC1"
y <- "PC2"
z <- "PC3"

x.var <- plsda.table[`Principal Component` == x, 
                     `% variance`]
y.var <- plsda.table[`Principal Component` == y, 
                     `% variance`]
z.var <- plsda.table[`Principal Component` == z, 
                     `% variance`]
fac.lvls <- unique(mSet$dataSet$cls)

chosen.colors <- c('#DF2935',
                   '#3772FF')

#chosen.colors <- if(fac.lvls == length(color.vec())) color.vec() else rainbow(length(fac.lvls))
# --- add ellipses ---
classes <- mSet$dataSet$cls
plots <- plotly::plot_ly(showlegend=F) 
#%>% layout(paper_bgcolor = 'rgb(235, 235, 235)',
#  plot_bgcolor = 'rgb(235, 235, 235)')

for(class in levels(classes)){
  row = which(classes == class)
  # ---------------------
  xc=coords[row, x]
  yc=coords[row, y]
  zc=coords[row, z]
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
rgbcols <- toRGB(chosen.colors)
c = 1
for(i in seq_along(adj_plot$x$data)){
  item = adj_plot$x$data[[i]]
  if(item$type == "mesh3d"){
    adj_plot$x$data[[i]]$color <- rgbcols[c]
    adj_plot$x$data[[i]]$visible <- TRUE
    adj_plot$x$data[[i]]$hoverinfo <- "none"
    c = c + 1
  }
}
# ---------------

countrysymbols <- c("Brazil"= "circle", 
                    "Brasil" = "circle",
                    "Spain"="square", 
                    "Netherlands"="diamond",
                    "Italy" = "cross")

plotA <- adj_plot %>%
  add_trace(
    x = coords[,x], 
    y = coords[,y], 
    z = coords[,z], 
    type = "scatter3d",
    marker = list(
      color = sapply(as.numeric(mSet$dataSet$cls), function(i) chosen.colors[i]),
      size = 9,
      symbol = sapply(as.character(as.factor(mSet$dataSet$covars$Country)),function(i) countrysymbols[[i]]),
      line = list(
        color = "black",
        width = 2
      )
    ),
    opacity=1,
    hoverinfo = 'text',
    text = rownames(mSet$dataSet$norm)
  ) %>%  layout(scene = list(
    aspectmode="cube",
    xaxis = list(
      titlefont = t,
      title = fn$paste("$x ($x.var %)"),
      gridcolor = 'rgb(255, 255, 255)'),
    yaxis = list(
      titlefont = t,
      title = fn$paste("$y ($y.var %)"),
      gridcolor = 'rgb(255, 255, 255)'),
    zaxis = list(
      titlefont = t,
      title = fn$paste("$z ($z.var %)"),
      gridcolor = 'rgb(255, 255, 255)'))
  )

plotA

# LEGEND

ax <- list(
  title = "",
  zeroline = FALSE,
  showline = FALSE,
  showticklabels = FALSE,
  showgrid = FALSE
)

data <- as.data.table(countrysymbols[-1],keep.rownames = T)

plot_ly(x  = c(-1,1,3,1, 1,1,1,0, 2),
        y = c(6,6,6,5, 4, 3, 2, 1, 1),
        text = c(".",".",".",data$V1, "Unhealthy", "Healthy"),#~V1,
        textposition= 'top',
        marker = list(
          color = c("white","white","white","white","white","white","white",chosen.colors),
          size = c(20,20,20,20,20,20,20,40,40),
          symbol = c("cross","cross","cross", data$V2, "circle", "circle"),
          line = list(
            color = c("black","black","black","black","black","black","black","white","white"),
            width = 2
          )
        )) %>%
  add_text(textfont = t, textposition = "top") %>%
  layout(xaxis = ax, yaxis = ax) 

# HISTOGRAM


df_both <- rbind(df1, df2)


df1 = {
  mSet <- swine_mset
  bw.vec1 <- mSet$analSet$plsda$permut
  len <- length(bw.vec)
  df <- melt(bw.vec)
  colnames(df) = "acc"
  df$animal = "swine"
  df
}
df2 = {
  mSet <- chicken_mset
  bw.vec2 <- mSet$analSet$plsda$permut
  len <- length(bw.vec)
  df <- melt(bw.vec)
  colnames(df) = "acc"
  df$animal = "chicken"
  df
}
df_both <- rbind(df1, df2)

chick_pval <- gsub(chicken_mset$analSet$plsda$permut.p, pattern = "0\\.00333333333333333", replacement = "0.003")
swine_pval <- gsub(swine_mset$analSet$plsda$permut.p, pattern = "0\\.00333333333333333", replacement = "0.003")

my_binwidth = (max(df_both$acc)-min(df_both$acc))/30;

ggplot(df_both,aes(x=acc)) + 
  #plot.theme(base_size = 15) +
  theme(text=element_text(size=16,  family="Trebuchet MS", hjust = 0.5)) +
  labs(title="PLSDA predictive value",x="Accuracy", y = "Permutations") +
  #geom_histogram(data=subset(df_both,animal == 'chicken'), fill = "yellow", color="black", alpha = 0.7,linetype=8, cex=.7) +
  geom_density(data=subset(df_both,animal == 'chicken'),aes(y = 0.01*..count..), fill = "yellow", color="black", alpha = 0.7,linetype=8, cex=.7) +
  geom_density(data=subset(df_both,animal == 'swine'), aes(y = 0.01*..count..), fill = "pink", color="black", alpha = 0.7,linetype=1, cex=.7) + 
  geom_segment(data=df_both,
               #color="yellow",
               x=bw.vec1[1],
               xend=bw.vec1[1],
               y=0,
               aes(yend=.05*nrow(df_both)),
               size=.7, 
               #cex=2
               linetype=8
  ) +
  geom_segment(data=df_both,
               #color="pink",
               x=bw.vec2[1],
               xend=bw.vec2[1],
               y=0,
               aes(yend=.05*nrow(df_both)),
               size=.7,
               #cex=2,
               linetype=1)+
  geom_label(aes(x=bw.vec1[1]-0.015,y=.05*nrow(df_both),label = chick_pval), color="black", label.size = 1, size=5) +
  geom_label(aes(x=bw.vec2[1]-0.015,y=.05*nrow(df_both),label = swine_pval), color="black", label.size = 1, size=5)

# ACCURACY CHARTS

require(ggplot2)

res <- chicken_mset$analSet$plsda$fit.info
colnames(res) <- 1:ncol(res)
# best.num <- mSet$analSet$plsda$best.num
# choice <- mSet$analSet$plsda$choice
df1 <- reshape2::melt(res)
df1$animal = "chicken"

res <- swine_mset$analSet$plsda$fit.info
colnames(res) <- 1:ncol(res)
# best.num <- mSet$analSet$plsda$best.num
# choice <- mSet$analSet$plsda$choice
df2 <- reshape2::melt(res)
df2$animal = "swine"

df <- rbind(df1, df2)

df$Component <- paste0("PC",df$Component)
colnames(df) <- c("Metric", "Component", "Value", "animal")
p <- ggplot(df, aes(x=Metric, y=Value, fill=animal)) +
  geom_bar(stat="identity",position='dodge', colour="black") + 
  theme_minimal() + 
  theme_gray() +
  facet_grid(~Component) +
  theme(text=element_text(size=15,  family="Trebuchet MS", hjust = 0.5)) +
  scale_fill_manual(values=c("yellow", "pink"))
p

# ROC
xvals <- chicken_mset$analSet$ml$ls$lasso_all$roc
xvals <- swine_mset$analSet$ml$ls$lasso_pig$roc

require(ROCR)
require(ggplot2)
require(data.table)

pred <- ROCR::prediction(xvals$predictions, xvals$labels)
perf <- ROCR::performance(pred, "tpr", "fpr")

# - - - old method - - -
# plot(perf,col="grey82",lty=3, cex.lab=1.3)
# plot(perf,lwd=3,avg="vertical",spread.estimate="boxplot",add=TRUE)

perf.long <<- rbindlist(lapply(1:length(perf@x.values), function(i){
  xvals <- perf@x.values[[i]]
  yvals <- perf@y.values[[i]]
  
  # - - - - - - - - -
  
  res <- data.table::data.table(attempt = c(i),
                                FPR = xvals,
                                TPR = yvals)
  res
}))
perf.long$animal = "pig"

perf2 = perf.long


perf.joined <- rbind(perf1, perf2)
cols = rainbow(50)

#nbins <- length(perf.long$TPR)

my_ci <- function(x){
  res = data.frame(
    y=mean(x), 
    ymin=mean(x) - 2 * sd(x), 
    ymax=mean(x) + 2 * sd(x))
  res
}

perf.joined[animal=="pig", animal:="swine"]
perf.joined[animal=="broilers", animal:="poultry"]

ggplot(perf.joined, aes(FPR, TPR, color=animal)) +
  #geom_smooth(fullrange=TRUE, method="loess", alpha= 0.8,color=perf.joined$col)+
  #geom_point()+
  stat_summary(fun.data="my_ci", geom="smooth", cex=2) +
  #stat_summary(fun.data="my_ci", geom="smooth", cex=2, color="black") +
  geom_abline(intercept=seq(0, 1, 25),
              slope=1,
              colour="gray",
              alpha = 0.5) +
  plot.theme(base_size = 10) +
  #ggplot2::scale_color_gradientn(colors = cols) +
  ggplot2::theme_gray() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=19,face="bold"),
        legend.title=element_text(size=15, face="bold"),
        legend.text=element_text(size=12),
        text=element_text(size=16,  family="Trebuchet MS", hjust = 0.5)
        #legend.position="none"
  ) +
  scale_color_manual(name = "Animal", values = c("yellow", "pink")) +
  coord_cartesian(xlim = c(.04,.96), ylim = c(.04,.96))

# VENN PREP

repeats <- chicken_mset$analSet$ml$ls$lasso_all$bar
repeats <- swine_mset$analSet$ml$ls$lasso_pig$bar

ml_type = repeats[[1]]$type

print(ml_type)

data <- switch(ml_type,
               rf = {
                 res <- aggregate(. ~ rn, rbindlist(lapply(repeats, function(x) as.data.table(x$feats, keep.rownames=T))), mean)
                 data <- res[order(res$MDA, decreasing = TRUE),]
                 colnames(data) <- c("mz", "mda")
                 # - - -
                 data
               },
               ls = {
                 feat_count <- lapply(repeats, function(x){
                   beta <- x$model$beta
                   feats <- which(beta[,1] > 0)
                   names(feats)
                 })
                 feat_count_tab <- table(unlist(feat_count))
                 feat_count_dt <- data.table::data.table(feat_count_tab)
                 colnames(feat_count_dt) <- c("mz", "count")
                 data <- feat_count_dt[order(feat_count_dt$count, decreasing = T)]
                 # - - -
                 data
               })
# data$mz <- round(as.numeric(data$mz), digits = 2)
data$mz <- factor(data$mz, levels=data$mz)

ml_bar_tab <<- data

venn2 <- data

p <- ggplot(data[1:topn,], aes(mz,count))
p + geom_bar(stat = "identity", aes(fill = count)) +
  geom_hline(aes(yintercept=attempts)) + 
  scale_fill_gradientn(colors=rainbow(20)) +
  theme(legend.position="none",
        axis.text=element_text(size=8),
        axis.title=element_text(size=13,face="bold"),
        #axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(x="Top hits",y="Times chosen")

# - - - - - VENN DIAGRAM - - - - - -

tophits = 100
chk.ls <- venn1[1:tophits,]
chk.tt <- chicken_mset$analSet$tt$sig.mat[order(chicken_mset$analSet$tt$sig.mat[,'p.value'], decreasing = F),][1:tophits,]
sw.ls <- venn2[1:tophits,]
sw.tt <- swine_mset$analSet$tt$sig.mat[order(swine_mset$analSet$tt$sig.mat[,'p.value'], decreasing = F),]

tables <- list(venn1[order(venn1$count, decreasing = T),]$mz, # lasso chicken
               venn2[order(venn2$count, decreasing = T),]$mz) # lasso swine

tables <- list(venn1[order(venn1$count, decreasing = T),]$mz, # lasso chicken
               rownames(chicken_mset$analSet$tt$sig.mat)[order(chicken_mset$analSet$tt$sig.mat[,"p.value"], decreasing = T)]) # t-test chicken


tables <- list("ml poultry" = chk.ls$mz, "t-test poultry" = rownames(chk.tt),
               "ml swine" = sw.ls$mz, "t-test swine" = rownames(sw.tt))

#names(tables) <- c("chicken", "pig")
# - - unlist - -

flattened <- flattenlist(tables)
names(flattened) <- gsub(x = names(flattened), pattern = "(.*\\.)(.*$)", replacement = "\\2")

circles = length(flattened)

# - - - - - - --

venn.plot <- VennDiagram::venn.diagram(x = flattened,
                                       filename = NULL)

# - - - - - - - - -

items <- strsplit(as.character(venn.plot), split = ",")[[1]]

circ_values <<- data.frame(
  id = 1:length(grep(items, pattern="polygon"))
  #,value = c(3, 3.1, 3.1, 3.2, 3.15, 3.5)
)

txt_values <- data.frame(
  id = grep(items, pattern="text"),
  value = unlist(lapply(grep(items, pattern="text"), function(i) venn.plot[[i]]$label))
)

txt_values$value <- gsub(x = txt_values$value, pattern = "(.*\\.)(.*$)", replacement = "\\2")
#categories <- c(categories, input$rf_choice, input$ls_choice, input$plsda_choice)

x_c = unlist(lapply(grep(items, pattern="polygon"), function(i) venn.plot[[i]]$x))
y_c = unlist(lapply(grep(items, pattern="polygon"), function(i) venn.plot[[i]]$y))

x_t = unlist(lapply(grep(items, pattern="text"), function(i) venn.plot[[i]]$x))
y_t = unlist(lapply(grep(items, pattern="text"), function(i)venn.plot[[i]]$y))

positions_c <- data.frame(
  id = rep(circ_values$id, each = length(x_c)/length(circ_values$id)),
  x = x_c,
  y = y_c
)

positions_t <- data.frame(
  id = rep(txt_values$id, each = length(x_t)/length(txt_values$id)),
  x = x_t,
  y = y_t
)

datapoly <- merge(circ_values, positions_c, by=c("id"))
datatxt <- merge(txt_values, positions_t, by=c("id"))

numbers <- datatxt[!(datatxt$value %in% names(flattened)),]
headers <- datatxt[(datatxt$value %in% names(flattened)),]

if(circles == 2){
  occur <- table(numbers$y)
  newy <- names(occur[occur == max(occur)])
  # - - -
  numbers$y <- as.numeric(c(newy))
}

datapoly <- as.data.table(datapoly)
datapoly[id==1, col:="yellow"]
datapoly[id==2, col:="#FC9E4F"]
datapoly[id==3, col:="pink"]
datapoly[id==4, col:="#DA4167"]

p <- ggplot(datapoly, 
            aes(x = x, 
                y = y)) + 
  geom_polygon(colour="black", alpha=0.3, aes(fill=id, group=id),fill=datapoly$col) +
  geom_text(mapping = aes(x=x, y=y, label=value), data = numbers, size = 7,family="Trebuchet MS") +
  geom_text(mapping = aes(x=x, y=y, label=value), data = headers, size = 9,family="Trebuchet MS") +
  theme_void() +
  theme(legend.position="none") + 
  scale_fill_gradientn(colours = rainbow(circles)) +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE)
p 

# Hypergeometric testing for venn diagrams?

require(gmp)

enrich_pvalue <- function(N, A, B, k)
{
  m <- A + k
  n <- B + k
  i <- k:min(m,n)
  # - - - 
  as.numeric( sum(gmp::chooseZ(m,i) * gmp::chooseZ(N-m,n-i)) / gmp::chooseZ(N,n) )
}

N = ncol(swine_mset$dataSet$norm) # pool of possible options
A = 93 # circle A
B = 6 # circle B
k = 3# overlap

{
  pval <- enrich_pvalue(N, A, B, k)
  
  stars <- if(pval < 0.0001){
   " (****)" 
  }else if(pval < 0.001){
    " (***)"
  }else if(pval < 0.01){
    " (**)"
  }else if(pval < 0.05){
    " (*)"
  }else{
    ""
  }
  
  print(paste0(pval, stars))
}

# COUNT ALL CPDS IN DB

cpd.count <- pbapply::pblapply(db_list, function(db){
  db_path <- paste0(options$db_dir, "/", db, ".full.db")
  print(db_path)
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), db_path) # change this to proper var later
  res <- DBI::dbGetQuery(conn, "SELECT count(*) FROM extended")
  # - - -
  res
})

sum(unlist(cpd.count))
