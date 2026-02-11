init_global_state <- function(env){
  evalq(
    {
      # === THE BELOW LIST CONTAINS ALL GLOBAL VARIABLES THAT METABOSHINY CALLS UPON LATER ===
      gbl <- list(
        constants = MetaboShiny:::build_gbl_constants(),
        functions = MetaboShiny:::build_gbl_functions(),
        paths = MetaboShiny:::build_gbl_paths(),
        vectors = MetaboShiny:::build_gbl_vectors(adducts)
      )
      
      try({
        gbl$vectors$kegg_pathways <- KEGGREST::keggList("pathway")
      }, silent = T)
      gbl$vectors$db_categories$all <- gbl$vectors$db_list
      gbl$vectors$example_sizes <- file.size(list.files(gbl$paths$volumes[["Examples"]],full.names = T))
      gbl$vectors$example_md5s <- tools::md5sum(list.files(gbl$paths$volumes[["Examples"]],full.names = T))
      gbl$vectors$wordcloud$filters <- list(
        stopwords = unique(c(tidytext::stop_words$word, 
                             qdapDictionaries::Top200Words,
                             tm::stopwords("english"))),
        metabolomics = c("metabolism", "metabolic",
                         "metabolomic", "metabolomics",
                         "biochemical", "mass", "spectrometry", 
                         "nmr", "direct", "infusion","exposome","papers",
                         "compounds","compound"),
        default = c(gbl$vectors$db_list, "exposome",
                    "Synonyms", "synonyms"))
      
      data(isotopes, package = "enviPat")
      
      radioTooltip <- function(id, choice, title, placement = "bottom", trigger = "hover", options = NULL){
        options = shinyBS:::buildTooltipOrPopoverOptionsList(title, placement, trigger, options)
        options = paste0("{'", paste(names(options), options, sep = "': '", collapse = "', '"), "'}")
        bsTag <- shiny::tags$script(shiny::HTML(paste0("
    $(document).ready(function() {
      setTimeout(function() {
        $('input', $('#", id, "')).each(function(){
          if(this.getAttribute('value') == '", choice, "') {
            opts = $.extend(", options, ", {html: true});
            $(this.parentElement).tooltip('destroy');
            $(this.parentElement).tooltip(opts);
          }
        })
      }, 500)
    });
  ")))
        htmltools::attachDependencies(bsTag, shinyBS:::shinyBSDep)
      }
      
      # interleave for sorting later ...
      add_idx <- order(c(seq_along(gbl$vectors$pos_adducts$Name), seq_along(gbl$vectors$neg_adducts$Name)))
      sort_order <<- unlist(c(gbl$vectors$pos_adducts$Name, gbl$vectors$neg_adducts$Name))[add_idx]
      
      session_cl <- parallel::makeCluster(1)
      debug_mSet <- list()
      debug_lcl <- list()
      debug_input <- list()
      
      msg.vec <- c()
      
      #  cliff delta
      USE_CLIFF_DELTA_R <- TRUE
      tryCatch(
        expr = {
          Rcpp::cppFunction('double cliffd(NumericVector x, NumericVector y){
         int len_x = x.size();
         int len_y = y.size();

         int sign_sum = 0;
         for(int i = 0; i < len_x; i++){
            for(int j = 0; j < len_y; j++){
               if(x[i] < y[j]){
                  sign_sum--;
               } else if(x[i] > y[j]){
                  sign_sum++;
               }
            }
         }

         return sign_sum / ((double)len_x*len_y);
      }
      ')
        },
        error = function(e){
          warning('C++ version of cliff delta could not be compiled. Using R version')
          USE_CLIFF_DELTA_R <<- TRUE
        }
      )
    },
    envir = env
  )
}
