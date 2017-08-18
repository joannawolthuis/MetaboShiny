multiform.joanna <- function (formula1, fact)
{
  formulas <- c()
  for (i in 1:length(formula1)) {
    formula1[i] <- gsub("^D(?![A-z])", "[2]H", formula1[i], perl=TRUE)
    ende1 <- nchar(formula1[i])
    element1 <- c()
    number1 <- c()
    j <- c(1)
    while (j <= ende1) {
      if (substr(formula1[i], j, j) == c("[")) {
        b <- j
        while (any(substr(formula1[i], j, j) == c("]")) !=
               TRUE) {
          j <- c(j + 1)
        }
        k <- j
        while (any(substr(formula1[i], j, j) == c("0",
                                                  "1", "2", "3", "4", "5", "6", "7", "8", "9")) !=
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        element1 <- c(element1, substr(formula1[i],
                                       b, m))
      }
      if (any(substr(formula1[i], j, j) == c("0", "1",
                                             "2", "3", "4", "5", "6", "7", "8", "9")) !=
          TRUE) {
        k <- j
        while (any(substr(formula1[i], j, j) == c("0",
                                                  "1", "2", "3", "4", "5", "6", "7", "8", "9")) !=
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        j <- c(j - 1)
        element1 <- c(element1, substr(formula1[i],
                                       k, m))
      }
      if (any(substr(formula1[i], j, j) == c("0", "1",
                                             "2", "3", "4", "5", "6", "7", "8", "9")) ==
          TRUE) {
        k <- j
        while (any(substr(formula1[i], j, j) == c("0",
                                                  "1", "2", "3", "4", "5", "6", "7", "8", "9")) ==
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        j <- c(j - 1)
        number1 <- c(number1, as.numeric(substr(formula1[i],
                                                k, m)))
      }
      j <- j + 1
    }
    number1 <- fact * number1
    formula_fin <- ""
    for (i in 1:length(element1)) {
      formula_fin <- paste(formula_fin, element1[i], number1[i],
                           sep = "")
    }
    formulas <- c(formulas, formula_fin)
  }
  return(formulas)
}

mergeform.joanna <- function (formula1, formula2)
{
  formula2 <- gsub("^D(?![A-z])", "[2]H", formula2, perl=TRUE)
  ende2 <- nchar(formula2)
  element2 <- c()
  number2 <- c()
  j <- c(1)
  while (j <= ende2) {
    if (substr(formula2, j, j) == c("[")) {
      b <- j
      while (any(substr(formula2, j, j) == c("]")) !=
             TRUE) {
        j <- c(j + 1)
      }
      k <- j
      while (any(substr(formula2, j, j) == c("0", "1",
                                             "2", "3", "4", "5", "6", "7", "8", "9")) !=
             TRUE) {
        j <- c(j + 1)
      }
      m <- c(j - 1)
      element2 <- c(element2, substr(formula2, b, m))
    }
    if (any(substr(formula2, j, j) == c("0", "1", "2", "3",
                                        "4", "5", "6", "7", "8", "9")) != TRUE) {
      k <- j
      while (any(substr(formula2, j, j) == c("0", "1",
                                             "2", "3", "4", "5", "6", "7", "8", "9")) !=
             TRUE) {
        j <- c(j + 1)
      }
      m <- c(j - 1)
      j <- c(j - 1)
      element2 <- c(element2, substr(formula2, k, m))
    }
    if (any(substr(formula2, j, j) == c("0", "1", "2", "3",
                                        "4", "5", "6", "7", "8", "9")) == TRUE) {
      k <- j
      while (any(substr(formula2, j, j) == c("0", "1",
                                             "2", "3", "4", "5", "6", "7", "8", "9")) ==
             TRUE) {
        j <- c(j + 1)
      }
      m <- c(j - 1)
      j <- c(j - 1)
      number2 <- c(number2, as.numeric(substr(formula2,
                                              k, m)))
    }
    j <- j + 1
  }
  formulas <- c()
  for (i in 1:length(formula1)) {
    formula1[i] <- gsub("^D(?![A-z])", "[2]H", formula1[i], perl=TRUE)
    ende1 <- nchar(formula1[i])
    element1 <- c()
    number1 <- c()
    j <- c(1)
    while (j <= ende1) {
      if (substr(formula1[i], j, j) == c("[")) {
        b <- j
        while (any(substr(formula1[i], j, j) == c("]")) !=
               TRUE) {
          j <- c(j + 1)
        }
        k <- j
        while (any(substr(formula1[i], j, j) == c("0",
                                                  "1", "2", "3", "4", "5", "6", "7", "8", "9")) !=
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        element1 <- c(element1, substr(formula1[i],
                                       b, m))
      }
      if (any(substr(formula1[i], j, j) == c("0", "1",
                                             "2", "3", "4", "5", "6", "7", "8", "9")) !=
          TRUE) {
        k <- j
        while (any(substr(formula1[i], j, j) == c("0",
                                                  "1", "2", "3", "4", "5", "6", "7", "8", "9")) !=
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        j <- c(j - 1)
        element1 <- c(element1, substr(formula1[i],
                                       k, m))
      }
      if (any(substr(formula1[i], j, j) == c("0", "1",
                                             "2", "3", "4", "5", "6", "7", "8", "9")) ==
          TRUE) {
        k <- j
        while (any(substr(formula1[i], j, j) == c("0",
                                                  "1", "2", "3", "4", "5", "6", "7", "8", "9")) ==
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        j <- c(j - 1)
        number1 <- c(number1, as.numeric(substr(formula1[i],
                                                k, m)))
      }
      j <- j + 1
    }
    both <- unique(c(element1, element2))
    counts <- c()
    for (i in 1:length(both)) {
      if (any(element1 == both[i])) {
        it1 <- c(number1[element1 == both[i]])
      }
      else {
        it1 <- c(0)
      }
      if (any(element2 == both[i])) {
        it2 <- c(number2[element2 == both[i]])
      }
      else {
        it2 <- c(0)
      }
      counts <- c(counts, it1 + it2)
    }
    formula_all <- ""
    for (i in 1:length(both)) {
      formula_all <- paste(formula_all, both[i], counts[i],
                           sep = "")
    }
    formulas <- c(formulas, formula_all)
  }
  return(formulas)
}

subform.joanna <- function (formula1, formula2)
{
  formula2 <- gsub("^D(?![A-z])", "[2]H", formula2, perl=TRUE)
  ende2 <- nchar(formula2)
  element2 <- c()
  number2 <- c()
  j <- c(1)
  while (j <= ende2) {
    if (substr(formula2, j, j) == c("[")) {
      b <- j
      while (any(substr(formula2, j, j) == c("]")) !=
             TRUE) {
        j <- c(j + 1)
      }
      k <- j
      while (any(substr(formula2, j, j) == c("0", "1",
                                             "2", "3", "4", "5", "6", "7", "8", "9")) !=
             TRUE) {
        j <- c(j + 1)
      }
      m <- c(j - 1)
      element2 <- c(element2, substr(formula2, b, m))
    }
    if (any(substr(formula2, j, j) == c("0", "1", "2", "3",
                                        "4", "5", "6", "7", "8", "9")) != TRUE) {
      k <- j
      while (any(substr(formula2, j, j) == c("0", "1",
                                             "2", "3", "4", "5", "6", "7", "8", "9")) !=
             TRUE) {
        j <- c(j + 1)
      }
      m <- c(j - 1)
      j <- c(j - 1)
      element2 <- c(element2, substr(formula2, k, m))
    }
    if (any(substr(formula2, j, j) == c("0", "1", "2", "3",
                                        "4", "5", "6", "7", "8", "9")) == TRUE) {
      k <- j
      while (any(substr(formula2, j, j) == c("0", "1",
                                             "2", "3", "4", "5", "6", "7", "8", "9")) ==
             TRUE) {
        j <- c(j + 1)
      }
      m <- c(j - 1)
      j <- c(j - 1)
      number2 <- c(number2, as.numeric(substr(formula2,
                                              k, m)))
    }
    j <- j + 1
  }
  formulas <- c()
  for (i in 1:length(formula1)) {
    formula1[i] <- gsub("^D(?![A-z])", "[2]H", formula1[i], perl=TRUE)
    ende1 <- nchar(formula1[i])
    element1 <- c()
    number1 <- c()
    j <- c(1)
    while (j <= ende1) {
      if (substr(formula1[i], j, j) == c("[")) {
        b <- j
        while (any(substr(formula1[i], j, j) == c("]")) !=
               TRUE) {
          j <- c(j + 1)
        }
        k <- j
        while (any(substr(formula1[i], j, j) == c("0",
                                                  "1", "2", "3", "4", "5", "6", "7", "8", "9")) !=
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        element1 <- c(element1, substr(formula1[i],
                                       b, m))
      }
      if (any(substr(formula1[i], j, j) == c("0", "1",
                                             "2", "3", "4", "5", "6", "7", "8", "9")) !=
          TRUE) {
        k <- j
        while (any(substr(formula1[i], j, j) == c("0",
                                                  "1", "2", "3", "4", "5", "6", "7", "8", "9")) !=
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        j <- c(j - 1)
        element1 <- c(element1, substr(formula1[i],
                                       k, m))
      }
      if (any(substr(formula1[i], j, j) == c("0", "1",
                                             "2", "3", "4", "5", "6", "7", "8", "9")) ==
          TRUE) {
        k <- j
        while (any(substr(formula1[i], j, j) == c("0",
                                                  "1", "2", "3", "4", "5", "6", "7", "8", "9")) ==
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        j <- c(j - 1)
        number1 <- c(number1, as.numeric(substr(formula1[i],
                                                k, m)))
      }
      j <- j + 1
    }
    formula_all <- TRUE
    for (i in 1:length(element2)) {
      if (any(element2[i] == element1) == FALSE) {
        formula_all <- paste(element2[i], " from formula 2 not part of formula1",
                             sep = "")
      }
      else {
        if (number2[i] > number1[element2[i] == element1]) {
          formula_all <- paste("Atom number of ", element2[i],
                               " from formula 2 not fully subset of formula1 atom number",
                               sep = "")
        }
      }
    }
    if (formula_all == TRUE) {
      formula_all <- ""
      counts <- c()
      for (i in 1:length(element1)) {
        if (any(element2 == element1[i])) {
          counts <- c(counts, number1[i] - (number2[element2 ==
                                                      element1[i]]))
        }
        else {
          counts <- c(counts, number1[i])
        }
      }
      element1 <- element1[counts > 0]
      counts <- counts[counts > 0]
      for (i in 1:length(counts)) {
        formula_all <- paste(formula_all, element1[i],
                             counts[i], sep = "")
      }
    }
    formulas <- c(formulas, formula_all)
  }
  return(formulas)
}

check.chemform.joanna <- function (isotopes, chemforms, get_sorted = FALSE, get_list = FALSE)
{
  multif <- function(formula1, fact, numbers) {
    formulas <- c()
    formula1 <- gsub("^D(?![A-z])", "[2]H", formula1, perl=TRUE)
    ende1 <- nchar(formula1)
    element1 <- c()
    number1 <- c()
    j <- c(1)
    while (j <= ende1) {
      if (substr(formula1, j, j) == c("[")) {
        b <- j
        while (any(substr(formula1, j, j) == c("]")) !=
               TRUE) {
          j <- c(j + 1)
        }
        k <- j
        while (any(substr(formula1, j, j) == numbers) !=
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        element1 <- c(element1, substr(formula1, b,
                                       m))
      }
      if (any(substr(formula1, j, j) == numbers) != TRUE) {
        k <- j
        while (any(substr(formula1, j, j) == numbers) !=
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        j <- c(j - 1)
        element1 <- c(element1, substr(formula1, k,
                                       m))
      }
      if (any(substr(formula1, j, j) == numbers) == TRUE) {
        k <- j
        while (any(substr(formula1, j, j) == numbers) ==
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        j <- c(j - 1)
        number1 <- c(number1, as.numeric(substr(formula1,
                                                k, m)))
      }
      j <- j + 1
    }
    number1 <- fact * number1
    formula_fin <- ""
    for (p in 1:length(element1)) {
      formula_fin <- paste(formula_fin, element1[p], number1[p],
                           sep = "")
    }
    formulas <- c(formulas, formula_fin)
    return(formulas)
  }
  capitals <- c("[", LETTERS)
  numbers <- as.character(0:10)
  allpossible <- c(capitals, numbers, "(", ")", "]", letters)
  masses <- c()
  warn <- c()
  elem <- unique(as.character(isotopes[, 1]))
  isotopes2 <- matrix(nrow = length(elem), ncol = 2)
  isotopes2[, 1] <- elem
  for (i in 1:length(elem)) {
    intermed <- isotopes[isotopes[, 1] == elem[i], ]
    if (is.vector(intermed) == TRUE) {
      isotopes2[, 2][isotopes2[, 1] == elem[i]] <- intermed[3]
    }
    else {
      isotopes2[, 2][isotopes2[, 1] == elem[i]] <- intermed[,
                                                            3][as.numeric(intermed[, 4]) == max(as.numeric(intermed[,
                                                                                                                    4]))]
    }
  }
  info <- isotopes2
  if (get_list)
    listed <- vector("list", length(chemforms))
  for (i in 1:length(chemforms)) {
    mass <- c(0)
    warnit <- FALSE
    if (chemforms[i] == "") {
      warnit <- TRUE
    }
    if (warnit == FALSE) {
      formel <- as.character(chemforms[i])
      formel <- strsplit(formel, " ")[[1]]
      m <- strsplit(formel, as.character())
      m <- m[[1]]
    }
    if (warnit == FALSE) {
      for (j in 1:length(m)) {
        if (any(allpossible == m[j]) == FALSE) {
          warnit <- TRUE
        }
      }
    }
    if (warnit == FALSE) {
      if (any(m == "[") || any(m == "]") || any(m == "(") ||
          any(m == ")")) {
        getit1 <- c(0)
        getit2 <- c(0)
        a <- c(1)
        while ((a) <= length(m)) {
          if (m[a] == "[") {
            getit1 <- getit1 + 1
          }
          if (m[a] == "]") {
            getit1 <- getit1 - 1
          }
          if (m[a] == "(") {
            getit2 <- getit2 + 1
          }
          if (m[a] == ")") {
            getit2 <- getit2 - 1
          }
          if (getit1 > 0 & (any(numbers == m[a]) ==
                            FALSE & m[a] != "[" & m[a] != "]")) {
            warnit <- TRUE
          }
          if (getit1 < 0 || getit2 < 0) {
            warnit <- TRUE
          }
          a <- a + 1
        }
        if (getit1 != 0 || getit2 != 0) {
          warnit <- TRUE
        }
      }
    }
    if (warnit == FALSE) {
      if (any(any(capitals == m[1]) == FALSE & m[1] !=
              "(" & m[1] != ")")) {
        warnit <- TRUE
      }
      if (length(m) == 1) {
        m <- c(m, "1")
        formel <- paste(formel, "1", sep = "")
      }
    }
    if (warnit == FALSE) {
      for (k in 2:length(m)) {
        if (m[k - 1] == "(" & m[k] == ")") {
          warnit <- TRUE
        }
      }
    }
    if (warnit == FALSE) {
      if (any(m == "(")) {
        m2 <- c()
        for (j in 1:(length(m) - 1)) {
          if (m[j] == ")" & any(m[j + 1] == numbers) ==
              FALSE) {
            m2 <- c(m2, m[j], "1")
          }
          else {
            m2 <- c(m2, m[j])
          }
        }
        m2 <- c(m2, m[length(m)])
        if (m[length(m)] == ")") {
          m2 <- c(m2, "1")
        }
        m <- m2
      }
      m2 <- m[1]
      for (j in 2:length(m)) {
        if ((any(m[j] == capitals) || m[j] == ")" ||
             m[j] == "(") & all(m[j - 1] != numbers) &
            m[j - 1] != "(" & m[j - 1] != "]") {
          m2 <- c(m2, "1", m[j])
        }
        else {
          m2 <- c(m2, m[j])
        }
      }
      if (all(m[length(m)] != numbers)) {
        m2 <- c(m2, "1")
      }
      m <- m2
      formel <- ""
      for (k in 1:length(m)) {
        formel <- paste(formel, m[k], sep = "")
      }
    }
    if (warnit == FALSE) {
      while (any(m == "(")) {
        a <- c(1)
        getit1 <- 1
        getit2 <- 1
        while (getit1 != 0 & getit2 != 0 & a <= length(m)) {
          if (m[a] == "(") {
            getit1 <- 2
            from <- a
          }
          if (m[a] == ")") {
            getit2 <- 2
            to <- a
          }
          if (getit1 == 2 & getit2 == 2) {
            b <- a + 1
            count <- ""
            while (any(m[b] == numbers & b <= length(m))) {
              count <- paste(count, m[b], sep = "")
              b <- b + 1
            }
            count <- as.numeric(count)
            m2 <- ""
            for (k in (from + 1):(to - 1)) {
              m2 <- paste(m2, m[k], sep = "")
            }
            m2 <- multif(m2, count, numbers)
            m2 <- strsplit(m2, as.character())[[1]]
            m3 <- c()
            doneit <- FALSE
            for (z in 1:length(m)) {
              if (z < from || z >= b) {
                m3 <- c(m3, m[z])
              }
              else {
                if (doneit == FALSE & (z >= from ||
                                       z < b)) {
                  m3 <- c(m3, m2)
                  doneit <- TRUE
                }
              }
            }
            m <- m3
            getit1 <- 0
            getit2 <- 0
          }
          a <- a + 1
        }
      }
      formel <- ""
      for (k in 1:length(m)) {
        formel <- paste(formel, m[k], sep = "")
      }
    }
    if (warnit == FALSE) {
      element1 <- c()
      number1 <- c()
      j <- c(1)
      while (j <= nchar(formel)) {
        if (substr(formel, j, j) == c("[")) {
          b <- j
          while (any(substr(formel, j, j) == c("]")) !=
                 TRUE & j <= nchar(formel)) {
            j <- c(j + 1)
          }
          k <- j
          while (any(substr(formel, j, j) == numbers) !=
                 TRUE) {
            j <- c(j + 1)
          }
          z <- c(j - 1)
          element1 <- c(element1, substr(formel, b,
                                         z))
        }
        if (any(substr(formel, j, j) == numbers) !=
            TRUE) {
          k <- j
          while (any(substr(formel, j, j) == numbers) !=
                 TRUE & j <= nchar(formel)) {
            j <- c(j + 1)
          }
          z <- c(j - 1)
          j <- c(j - 1)
          element1 <- c(element1, substr(formel, k,
                                         z))
        }
        if (any(substr(formel, j, j) == numbers) ==
            TRUE) {
          k <- j
          while (any(substr(formel, j, j) == numbers) ==
                 TRUE & j <= nchar(formel)) {
            j <- c(j + 1)
          }
          z <- c(j - 1)
          j <- c(j - 1)
          number1 <- c(number1, as.numeric(substr(formel,
                                                  k, z)))
        }
        j <- j + 1
      }
    }
    if (warnit == FALSE) {
      for (j in 1:length(element1)) {
        if (any(element1[j] == as.character(isotopes[,
                                                     1])) == FALSE) {
          warnit <- TRUE
        }
      }
      if (length(element1) != length(number1)) {
        warnit <- TRUE
      }
    }
    if (warnit == FALSE) {
      element2 <- c()
      number2 <- c()
      doneit <- rep(FALSE, length(element1))
      for (j in 1:length(element1)) {
        if (doneit[j] == FALSE) {
          doneit[element1 == element1[j]] <- TRUE
          element2 <- c(element2, element1[element1 ==
                                             element1[j]][1])
          number2 <- c(number2, as.character(sum(as.numeric(number1[element1 ==
                                                                      element1[j]]))))
        }
      }
      element1 <- element2
      rm(element2)
      number1 <- number2
      rm(number2)
      if (get_sorted) {
        this <- order(match(element1, info))
        number1 <- number1[this]
        element1 <- element1[this]
      }
      formel <- ""
      for (k in 1:length(element1)) {
        formel <- paste(formel, element1[k], number1[k],
                        sep = "")
        mass <- mass + (as.numeric(info[info[, 1] ==
                                          element1[k], 2][1]) * as.numeric(number1[k]))
      }
    }
    if (warnit == FALSE) {
      if (!get_list) {
        warn <- c(warn, FALSE)
        masses <- c(masses, mass)
        chemforms[i] <- formel
      }
      else {
        number1 <- as.numeric(number1)
        names(number1) <- element1
        listed[[i]] <- number1
        names(listed)[i] <- chemforms[i]
      }
    }
    else {
      if (!get_list) {
        warn <- c(warn, TRUE)
        masses <- c(masses, -9999)
      }
      else {
        listed[[i]] <- numeric()
        names(listed)[i] <- "invalid formula"
      }
    }
  }
  if (!get_list) {
    checked <- data.frame(warn, chemforms, masses)
    names(checked) <- c("warning", "new_formula", "monoisotopic_mass")
    checked[, 2] <- as.character(checked[, 2])
    return(checked)
  }
  else {
    return(listed)
  }
}

check.ded.joanna <- function(formulas, deduct)
{
  warn <- c()
  formel2 <- as.character(deduct)
  formel2 <- gsub("^D(?![A-z])", "[2]H", formel2, perl=TRUE)
  ende2 <- nchar(formel2)
  element2 <- c()
  number2 <- c()
  j <- c(1)
  while (j <= ende2) {
    if (substr(formel2, j, j) == c("[")) {
      b <- j
      while (any(substr(formel2, j, j) == c("]")) != TRUE) {
        j <- c(j + 1)
      }
      k <- j
      while (any(substr(formel2, j, j) == c("0", "1",
                                            "2", "3", "4", "5", "6", "7", "8", "9")) !=
             TRUE) {
        j <- c(j + 1)
      }
      m <- c(j - 1)
      element2 <- c(element2, substr(formel2, b, m))
    }
    if (any(substr(formel2, j, j) == c("0", "1", "2", "3",
                                       "4", "5", "6", "7", "8", "9")) != TRUE) {
      k <- j
      while (any(substr(formel2, j, j) == c("0", "1",
                                            "2", "3", "4", "5", "6", "7", "8", "9")) !=
             TRUE) {
        j <- c(j + 1)
      }
      m <- c(j - 1)
      j <- c(j - 1)
      element2 <- c(element2, substr(formel2, k, m))
    }
    if (any(substr(formel2, j, j) == c("0", "1", "2", "3",
                                       "4", "5", "6", "7", "8", "9")) == TRUE) {
      k <- j
      while (any(substr(formel2, j, j) == c("0", "1",
                                            "2", "3", "4", "5", "6", "7", "8", "9")) ==
             TRUE) {
        j <- c(j + 1)
      }
      m <- c(j - 1)
      j <- c(j - 1)
      number2 <- c(number2, as.numeric(substr(formel2,
                                              k, m)))
    }
    j <- j + 1
  }
  element3 <- c()
  number3 <- c()
  getit <- as.character(levels(as.factor(element2)))
  for (j in 1:length(getit)) {
    element3 <- c(element3, getit[j])
    number3 <- c(number3, sum(number2[element2 == getit[j]]))
  }
  element2 <- element3
  number2 <- number3
  for (i in 1:length(formulas)) {
    warn <- c(warn, FALSE)
    element1 <- c()
    number1 <- c()
    formel1 <- as.character(formulas[i])
    formel1 <- gsub("^D(?![A-z])", "[2]H", formel1, perl=TRUE)
    ende1 <- nchar(formel1)
    j <- c(1)
    while (j <= ende1) {
      if (substr(formel1, j, j) == c("[")) {
        b <- j
        while (any(substr(formel1, j, j) == c("]")) !=
               TRUE) {
          j <- c(j + 1)
        }
        k <- j
        while (any(substr(formel1, j, j) == c("0", "1",
                                              "2", "3", "4", "5", "6", "7", "8", "9")) !=
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        element1 <- c(element1, substr(formel1, b, m))
      }
      if (any(substr(formel1, j, j) == c("0", "1", "2",
                                         "3", "4", "5", "6", "7", "8", "9")) != TRUE) {
        k <- j
        while (any(substr(formel1, j, j) == c("0", "1",
                                              "2", "3", "4", "5", "6", "7", "8", "9")) !=
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        j <- c(j - 1)
        element1 <- c(element1, substr(formel1, k, m))
      }
      if (any(substr(formel1, j, j) == c("0", "1", "2",
                                         "3", "4", "5", "6", "7", "8", "9")) == TRUE) {
        k <- j
        while (any(substr(formel1, j, j) == c("0", "1",
                                              "2", "3", "4", "5", "6", "7", "8", "9")) ==
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        j <- c(j - 1)
        number1 <- c(number1, as.numeric(substr(formel1,
                                                k, m)))
      }
      j <- j + 1
    }
    element3 <- c()
    number3 <- c()
    getit <- as.character(levels(as.factor(element1)))
    for (j in 1:length(getit)) {
      element3 <- c(element3, getit[j])
      number3 <- c(number3, sum(number1[element1 == getit[j]]))
    }
    element1 <- element3
    number1 <- number3
    for (j in 1:length(element2)) {
      if (any(element2[j] == element1) == FALSE) {
        warn[i] <- c(TRUE)
      }
      else {
        if (as.numeric(number2[element2 == element2[j]]) >
            as.numeric(number1[element1 == element2[j]])) {
          warn[i] <- c(TRUE)
        }
      }
    }
  }
  return(warn)
}
