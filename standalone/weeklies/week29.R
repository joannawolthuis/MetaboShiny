# !!! local 

# - - - SEPERATE TESTS - - -

csv_loc <- "/Users/jwolthuis/Analysis/SP/BrazilAndSpain_W.csv"

csv <- fread(csv_loc, header = T)

csv[1:10,1:10]

require(data.table)

split.by.country <- split(csv, by = "Country")
split.by.farm <- split(csv, by = c("Country","Farm"))

current.exp <- split.by.country

results <- pbapply::pblapply(1:length(current.exp), function(i){
  name <- names(current.exp)[i]
  table <- current.exp[[i]]
  file.name <- paste0(name, "_", basename(csv_loc))
  new.loc <- file.path(dirname(csv_loc), file.name)
  print(new.loc)
  # - get QCs -
  batches <- unique(table$Batch)
  print(batches)
  QCs <- csv[Sample %like% "^QC" & Batch %in% batches]
  result <- rbind(table, QCs)
  
  print(result[,1:20])
  # - - -
  fwrite(x = result, file = new.loc)
})

# - - - store metaboshiny results - - -

week29 <- list()

#week29$Brazil <- mSet$analSet
#week29$Spain <- mSet$analSet
week29$Both <- mSet$analSet

save(week29, file="week29_spainbrazilboth.RData")
#week29$farm1 <- mSet$analSet

# - - CROSS COUNTRY - -

tophits = 200

flattened <- list(sp_tt = rownames(week29$Spain$tt$sig.mat), br_tt = rownames(week29$Brazil$tt$sig.mat), both_tt = rownames(week29$Both$tt$sig.mat))
flattened <- list(sp_fc = rownames(week29$Spain$fc$sig.mat), br_fc = rownames(week29$Brazil$fc$sig.mat), both_fc = rownames(week29$Both$fc$sig.mat))
flattened <- list(sp_vc = rownames(week29$Spain$volcano$sig.mat), br_vc = rownames(week29$Brazil$volcano$sig.mat), both_vc = rownames(week29$Both$volcano$sig.mat))
flattened <- list(sp_ls = week29$Spain$ml$ls[[1]]$mz, br_ls = week29$Brazil$ml$ls[[1]]$mz, both=week29$Both$ml$ls[[1]]$mz)
flattened <- list(sp_rf = week29$Spain$ml$rf[[1]]$mz[1:tophits], br_ls = week29$Brazil$ml$rf[[1]]$mz[1:tophits], both=week29$Both$ml$rf[[1]]$mz[1:tophits])

flattened <- list(sp_vc = rownames(week29$Spain$volcano$sig.mat), 
                  br_vc = rownames(week29$Brazil$volcano$sig.mat),
                  sp_ls = week29$Spain$ml$ls[[1]]$mz, 
                  br_ls = week29$Brazil$ml$ls[[1]]$mz )

circles <- length(flattened)

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
  newy <- rep(min(as.numeric(names(occur))), times = 3)
  # - - -
  numbers$y <- as.numeric(c(newy))
}

p <- ggplot(datapoly, 
            aes(x = x, 
                y = y)) + 
  geom_polygon(colour="black", alpha=0.5, aes(fill=id, group=id)) +
  geom_text(mapping = aes(x=x, y=y, label=value), data = numbers, size = 5) +
  geom_text(mapping = aes(x=x, y=y, label=value), data = headers, fontface="bold", size = 7) +
  theme_void() +
  theme(legend.position="none") + 
  scale_fill_gradientn(colours = rainbow(circles)) +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE)

ggplotly(p, tooltip = "label") %>% layout(plot_bgcolor='white')

# - - whodunnit? - -

cpds <- Reduce(intersect, x = flattened)

# use the combined one as patdb...
patdb <- "/Users/jwolthuis/Analysis/SP/BrazilAndSpain_W.db"

get_matches(cpd = cpds[1], file.path(options$db_dir, "hmdb.full.db"),searchid = "mz")
get_matches(cpd = cpds[2], file.path(options$db_dir, "hmdb.full.db"),searchid = "mz")


# - - - checking subdirectory search - - -

mountDir <- "/Users/jwolthuis/MountPoint/Data/Metabolomics/DSM/"

all_dirs <- list.dirs(path = mountDir, recursive = T)
res_dirs <- grep(all_dirs, pattern = "results$",value = T)
