## calculate relative abundancies from theoretical mass and composition

#library(Rdisop)

#theor.table <- read.table(file="C:/Users/mraves/Metabolomics/TheoreticalMZ_Negative_composition.txt", sep="\t", header=TRUE)
options(digits=16)

#library(OrgMassSpecR)

library(lattice)
# The following list was copied from Rdisop elements.R and corrected for C, H, O, Cl, S according to NIST
Ba <- list(name= 'Ba', mass=130, isotope=list(mass=c(-0.093718, 0, -0.094958, 0, -0.095514, -0.094335, -0.095447, -0.094188, -0.094768),                                      abundance=c(0.00106, 0, 0.00101, 0, 0.02417, 0.06592, 0.07854, 0.1123, 0.717)))
Br <- list(name= 'Br', mass=79,  isotope=list(mass=c(-0.0816639, 0, -0.083711),                                                                                               abundance=c(0.5069, 0, 0.4931)))
C <-  list(name= 'C',  mass=12,  isotope=list(mass=c(0, 0.003354838, 0.003241989),                                                                                            abundance=c(0.9893, 0.0107, 0)))
Ca <- list(name= 'Ca', mass=40,  isotope=list(mass=c(-0.0374094, 0, -0.0413824, -0.0412338, -0.0445194, 0, -0.046311, 0, -0.047467),                                          abundance=c(0.96941, 0, 0.00647, 0.00135, 0.02086, 0, 4e-05, 0, 0.00187)))
Cl <- list(name= 'Cl', mass=35,  isotope=list(mass=c(-0.03114732, 0, -0.03409741),                                                                                            abundance=c(0.7576, 0, 0.2424)))
Cr <- list(name= 'Cr', mass=50,  isotope=list(mass=c(-0.0539536, 0, -0.0594902, -0.0593487, -0.0611175),                                                                      abundance=c(0.04345, 0, 0.83789, 0.09501, 0.02365)))
Cu <- list(name= 'Cu', mass=63,  isotope=list(mass=c(-0.0704011, 0, -0.0722071),                                                                                              abundance=c(0.6917, 0, 0.3083)))
F <-  list(name= 'F',  mass=19,  isotope=list(mass=c(-0.00159678),                                                                                                            abundance=c(1)))
Fe <- list(name= 'Fe', mass=54,  isotope=list(mass=c(-0.0603873, 0, -0.0650607, -0.0646042, -0.0667227),                                                                      abundance=c(0.058, 0, 0.9172, 0.022, 0.0028)))
H <-  list(name= 'H',  mass=1,   isotope=list(mass=c(0.00782503207, 0.014101778, 0.01604928),                                                                                 abundance=c(0.999885, 0.000115, 0)))
Hg <- list(name= 'Hg', mass=196, isotope=list(mass=c(-0.034193, 0, -0.033257, -0.031746, -0.0317, -0.029723, -0.029383, 0, -0.026533),                                        abundance=c(0.0015, 0, 0.0997, 0.1687, 0.231, 0.1318, 0.2986, 0, 0.0687)))
I <-  list(name= 'I',  mass=127, isotope=list(mass=c(-0.095527),                                                                                                              abundance=c(1)))
K <-  list(name= 'K',  mass=39,  isotope=list(mass=c(-0.0362926, -0.0360008, -0.0381746),                                                                                     abundance=c(0.932581, 0.000117, 0.067302)))
Li <- list(name= 'Li', mass=6,   isotope=list(mass=c(0.0151214, 0.016003),                                                                                                    abundance=c(0.075, 0.925)))
Mg <- list(name= 'Mg', mass=24,  isotope=list(mass=c(-0.0149577, -0.0141626, -0.0174063),                                                                                     abundance=c(0.7899, 0.1, 0.1101)))
Mn <- list(name= 'Mn', mass=55,  isotope=list(mass=c(-0.0619529),                                                                                                             abundance=c(1)))
N <-  list(name= 'N',  mass=14,  isotope=list(mass=c(0.003074002, 0.00010897),                                                                                                abundance=c(0.99634, 0.00366)))
Na <- list(name= 'Na', mass=23,  isotope=list(mass=c(-0.0102323),                                                                                                             abundance=c(1)))
Ni <- list(name= 'Ni', mass=58,  isotope=list(mass=c(-0.0646538, 0, -0.0692116, -0.0689421, -0.0716539, 0, -0.0720321),                                                       abundance=c(0.68077, 0, 0.26223, 0.0114, 0.03634, 0, 0.00926)))
O <-  list(name= 'O',  mass=16,  isotope=list(mass=c(-0.00508538044, -0.0008683, -0.0008397),                                                                                 abundance=c(0.99757, 0.000381, 0.00205)))
P <-  list(name= 'P',  mass=31,  isotope=list(mass=c(-0.026238),                                                                                                              abundance=c(1)))
S <-  list(name= 'S',  mass=32,  isotope=list(mass=c(-0.027929, -0.02854124, -0.0321331, 0, -0.03291924),                                                                     abundance=c(0.9499, 0.0075, 0.0425, 0, 1e-04)))
Se <- list(name= 'Se', mass=74,  isotope=list(mass=c(-0.0775254, 0, -0.080788, -0.0800875, -0.0826924, 0, -0.0834804, 0, -0.0833022),                                         abundance=c(0.0089, 0, 0.0936, 0.0763, 0.2378, 0, 0.4961, 0, 0.0873)))
Si <- list(name= 'Si', mass=28,  isotope=list(mass=c(-0.0230729, -0.0235051, -0.0262293),                                                                                     abundance=c(0.9223, 0.0467, 0.031)))
Sn <- list(name= 'Sn', mass=112, isotope=list(mass=c(-0.095174, 0, -0.097216, -0.096652, -0.098253, -0.097044, -0.098391, -0.09669, -0.0978009, 0, -0.0965596, 0, -0.0947257),abundance=c(0.0097, 0, 0.0065, 0.0034, 0.1453, 0.0768, 0.2423, 0.0859, 0.3259, 0, 0.0463, 0, 0.0579)))
Zn <- list(name= 'Zn', mass=64,  isotope=list(mass=c(-0.0708552, 0, -0.0739653, -0.0728709, -0.0751541, 0, -0.074675),                                                        abundance=c(0.486, 0, 0.279, 0.041, 0.188, 0, 0.006)))

NH4 <- list(name= "NH4", mass=18, isotope=list(mass=c(0.03437, 0.03141, -0.95935)), abundance=c(0.995, 0.004, 0.001)) # SISweb: 18.03437 100 19.03141 0.4 19.04065 0.1
Ac <- list(name= "Ac", mass=60, isotope=list(mass=c(0.02114, 0.02450, 0.02538)), abundance=c(0.975, 0.021, 0.004))   # SISweb: 60.02114 100 61.02450 2.2 62.02538 0.4
NaCl <- list(name= "NaCl", mass=58, isotope=list(mass=c(-0.04137, 0, -0.04433)), abundance=c(0.755, 0, 0.245))   # SISweb: 57.95862 100 59.95567 32.4
NaCl2 <- list(name= "NaCl2", mass=116, isotope=list(mass=c(-0.08274, 0, -0.08866)), abundance=c(0.755, 0, 0.245))   # SISweb: 57.95862 100 59.95567 32.4
NaCl3 <- list(name= "NaCl3", mass=174, isotope=list(mass=c(-0.12411, 0, -0.13299)), abundance=c(0.755, 0, 0.245))   # SISweb: 57.95862 100 59.95567 32.4
NaCl4 <- list(name= "NaCl4", mass=232, isotope=list(mass=c(-0.16548, 0, -0.17732)), abundance=c(0.755, 0, 0.245))   # SISweb: 57.95862 100 59.95567 32.4
NaCl5 <- list(name= "NaCl5", mass=290, isotope=list(mass=c(-0.20685, 0, -0.22165)), abundance=c(0.755, 0, 0.245))   # SISweb: 57.95862 100 59.95567 32.4
For <- list(name= "For", mass=45, isotope=list(mass=c(-0.00233, 0.00103)), abundance=c(0.989, 0.011))   # SISweb: 46.00549 100 47.00885 1.1 (47.0097 0.1) 48.00973 0.4
Na2 <- list(name= "2Na-H", mass=46, isotope=list(mass=c(-1.0282896)), abundance=c(1))  # SISweb for Na2: 45.97954 100 # minus 1 H !
Met <- list(name= "CH3OH", mass=32, isotope=list(mass=c(1.034045,1.037405)), abundance=c(0.989,0.011)) # SISweb: 32.02622 100 33.02958 1.1 33.0325 0.1 34.03046 0.2
CH3OH <- list(name= "CH3OH", mass=32, isotope=list(mass=c(1.034045,1.037405)), abundance=c(0.989,0.011)) # SISweb: 32.02622 100 33.02958 1.1 33.0325 0.1 34.03046 0.2
Na3 <- list(name= "3Na-2H", mass=69, isotope=list(mass=c(-2.0463469)), abundance=c(1))  # SISweb for Na2: 45.97954 100 # minus 1 H !
KCl <- list(name= "KCl", mass=74, isotope=list(mass=c(-0.06744,0.92961,-0.06744,0.92772)), abundance=c(0.7047,0.2283,0.0507,0.0162))   # SISweb: 73.93256 100 75.92961 32.4 75.93067 7.2 77.92772 2.3
H2PO4 <- list(name= "H2PO4", mass=97, isotope=list(mass=c(-0.03091)), abundance=c(1))
HSO4 <- list(name= "HSO4", mass=97, isotope=list(mass=c(-0.04042,0,-0.04462)), abundance=c(0.96,0,0.04))
Met2 <- list(name= "Met2", mass=64, isotope=list(mass=c(1.060265,1.013405)), abundance=c(0.978,0.022))
Met3 <- list(name= "Met3", mass=96, isotope=list(mass=c(1.086485,1.089845)), abundance=c(0.969,0.031))
Met4 <- list(name= "Met4", mass=128, isotope=list(mass=c(1.112705,1.116065)), abundance=c(0.959,0.041))
Met5 <- list(name= "Met5", mass=160, isotope=list(mass=c(1.20935,1.142285)), abundance=c(0.949,0.051))
NaminH <- list(name= "Na-H", mass=21, isotope=list(mass=c(-0.02571416)), abundance=c(1))
KminH <- list(name= "K-H", mass=37, isotope=list(mass=c(-0.05194,0.94617)), abundance=c(0.9328,0.0672))
H2O <- list(name= "H2O", mass=-19, isotope=list(mass=c(-0.01894358)), abundance=c(1))
NaK <- list(name= "NaK-H", mass=61, isotope=list(mass=c(-0.054345,0.943765)), abundance=c(0.9328,0.0672))
min2H <- list(name= "min2H", mass=-2, isotope=list(mass=c(-0.0151014)), abundance=c(1))
plus2H <- list(name= "plus2H", mass=2, isotope=list(mass=c(0.0151014)), abundance=c(1))
plus2Na <- list(name= "plus2Na", mass=46,  isotope=list(mass=c(-0.02046), abundance=c(1)))
plusNaH <- list(name= "plusNaH", mass=24,  isotope=list(mass=c(-0.00295588), abundance=c(1)))
plusKH <- list(name= "plusKH", mass=40, isotope=list(mass=c(-0.029008,0.969101)), abundance=c(0.9328,0.0672))
plusHNH <- list(name= "plusHNH", mass=19, isotope=list(mass=c(0.04164642)), abundance=c(1))
min3H <- list(name= "min3H", mass=-3, isotope=list(mass=c(-0.02182926)), abundance=c(1))
plus3H <- list(name= "plus3H", mass=3, isotope=list(mass=c(0.02182926)), abundance=c(1))
plus3Na <- list(name= "plus3Na", mass=68, isotope=list(mass=c(0.96931)), abundance=c(1))
plus2NaH <- list(name= "plus2NaH", mass=47, isotope=list(mass=c(0.2985712)), abundance=c(1))
plusNa2H <- list(name= "plusNa2H", mass=25, isotope=list(mass=c(0.00432284)), abundance=c(1))

Cl37 <- list(name= "Cl37", mass=37,  isotope=list(mass=c(-0.03409741)),  abundance=c(1)) 

allelements <- list(Ba, Br, C, Ca, Cl, Cr, Cu, F, Fe, H, Hg, I, K, Li, Mg, Mn, N, Na, Ni, O, P, S, Se, Si, Sn, Zn)
allAdducts <- list(Ba, Br, Ca, Cl, Cl37, Cr, Cu, Fe, Hg, I, K, Li, Mg, Mn, Na, Ni, Se, Si, Sn, Zn, NH4, Ac, NaCl, For,Na2,CH3OH,NaCl2,NaCl3,NaCl4,NaCl5,Na3,KCl,H2PO4,HSO4,Met2,Met3,Met4,Met5,NaminH,KminH,H2O,NaK,min2H,plus2H,plus2Na,plusNaH,plusKH,min3H,plus3H,plusHNH,plus3Na,plus2NaH,plusNa2H)

#atomsinuse <- c("P",     "O",         "N",           "C",     "H",          "S",        "Cl")
#atomicWts <- c(30.97376163, 15.99491463, 14.0030740052, 12.0000, 1.0078250321, 31.9720707, 34.968852721)
#electron <- 0.00054858
#Mol.comp <- c(0,4,0,0,1,1,0)
#Mol.exact <- sum(Mol.comp * atomicWts) + electron
#Mol.corr <- Mol.exact - 0.0022 + 0.000007*Mol.exact  # mass as found in peak group list
#getMass(Mol) + electron - 0.0022 + 0.000007*getMass(Mol)

#atomsinuse <- c("P",     "O",         "N",           "C",     "H",          "S",        "Cl",          "D",           "34S",          "18O")
#atomicWts <- c(30.97376163, 15.99491463, 14.0030740052, 12.0000, 1.0078250321, 31.9720707, 34.968852721, 2.0141017778, 33.96786690, 17.9991610)
#electron <- 0.00054858
#Mol.comp <- c(0,4,0,0,1,1,0,0,0,0)  # main peak HSO4
#Mol.comp <- c(0,4,0,0,0,1,0,1,0,0)  # deuterated
#Mol.comp <- c(0,4,0,0,1,0,0,0,1,0)  # 34S
#Mol.comp <- c(0,3,0,0,1,1,0,0,0,1)  # 18O
#Mol.exact <- sum(Mol.comp * atomicWts) + electron
#Mol.corr <- Mol.exact - 0.0022 + 0.000007*Mol.exact  # mass as found in peak group list

atomsinuse <- c("P",     "O",         "N",           "C",     "H",          "S",        "Cl",          "D",           "13C",       "34S",          "18O",     "37Cl")
atomicWts <- c(30.97376163, 15.99491463, 14.0030740052, 12.0000, 1.0078250321, 31.9720707, 34.968852721, 2.0141017778, 13.0033548378, 33.96786690, 17.9991610, 36.96590259)
electron <- 0.00054858
############# P O N C H  S Cl D 13C 34S 18O 37Cl
#Mol.comp <- c(0,6,0,6,12,0,1, 0, 0, 0,  0,  0)  # main peak Galactose HCl negative ion
#Mol.comp <- c(0,6,0,5,12,0,1, 0, 1, 0,  0,  0)  # 13C Galactose HCl negative ion
#Mol.comp <- c(0,6,0,6,11,0,1, 1, 0, 0,  0,  0)  # D Galactose HCl negative ion
#Mol.comp <- c(0,5,0,6,12,0,1, 0, 0, 0,  1,  0)  # 18O Galactose HCl negative ion
#Mol.comp <- c(0,6,0,6,12,0,0, 0, 0, 0,  0,  1)  # 18O Galactose HCl negative ion
#Mol.exact <- sum(Mol.comp * atomicWts) + electron
#Mol.corr <- Mol.exact - 0.0022 + 0.000007*Mol.exact  # mass as found in peak group list

Hmass <- H$mass + H$isotope$mass[1]
Dmass <- H$mass + 1 + H$isotope$mass[2]
Tmass <- H$mass + 2 + H$isotope$mass[3]
C13mass <- C$mass + 1 + C$isotope$mass[2]
N15mass <- N$mass + 1 + N$isotope$mass[2]
