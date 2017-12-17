#' @export
get.mw <- function(mol.formula, charge){
  # --------------------------------------
  #library(CIAAWconsensus)
  #library(CHNOSZ)
  # --------------------------------------
  atom.masses <- CIAAWconsensus::ciaaw.mass
  sodium.row <- data.frame(
    isotope = '22Na',
    element = 'sodium',
    mass = 22.98976928,
    uncertainty = 2.0e-9
  ) # sodium value from ciaaw website
  # ------
  electron.mass <- 0.00054858
  proton.mass <- 1.00727646677
  # ------
  atom.masses <- rbind(atom.masses, sodium.row)
  atom.masses$symbol <- gsub('\\d', '', atom.masses$isotope)
  # --------------------------------------
  if(mol.formula == '2Na-H'){
    mol.formula <- 'Na2-H'
  }
  # --------------------------------------
  atom.df <- makeup(mol.formula)
  atom.makeup <- names(atom.df)
  masses <- lapply(1:length(atom.makeup), FUN=function(atom.index, atom.df, atom.makeup){
    atom.id <- atom.makeup[[atom.index]]
    atom.count <- atom.df[atom.id]
    atom.weight <- min(atom.masses[atom.masses$symbol == atom.id, 'mass'])
    total.weight <- atom.count * atom.weight + -(charge)*electron.mass
  }, atom.df=atom.df, atom.makeup=atom.makeup)
  Reduce("+", masses)[[1]]
}
