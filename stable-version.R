#      ___
#     /\  \
#    /::\  \
#   /:/\:\__\
#  /:/ /:/  /
# /:/_/:/__/___
# \:\/:::::/  /
#  \::/~~/~~~~
#   \:\~~\
#    \:\__\
#     \/__/
#################################################
# global POSIX character patterns
pattern_MRN = "[[:digit:]]{8}-[[:digit:]]{1}"
# pattern_BiopsyNum = "[[:upper:]]{2}-[[:digit:]]{2}-[[:digit:]]{5}"
pattern_BiopsyNum = "SD-[[:digit:]]{2}-[[:digit:]]{5}"
pattern_Loc = "[[:upper:]]{3}"
pattern_Date = "[[:digit:]]{2}/[[:digit:]]{2}/[[:digit:]]{2}"
pattern_Diagnosis = ""
pattern_Specimen = ""
pattern_Histology = ""
pattern_Ulceration = ""
pattern_Depth = ""
pattern_Regression = ""
pattern_SunDam = ""
pattern_PeriMarg = ""
pattern_DeepMarg = ""
pattern_TNM = "pT[[:digit:]]{1}x"

# preliminary functions
EncodeMRN <- function(MRN) {
# MRN_1 = (EMRN_1 - b)*p^-1 %% 10^8
# MRN_2 = (EMRN_2 - b)*p^-1 %% 10^1
# for p = 17, p^-1 mod 10^8 is 5882353
# for p = 17, p^-1 mod 10^1 is 3
  # Args:
  #   MRN: an MRN
  # Returns:
  #   encoded MRN
  a_1 = as.integer(substr(MRN, 1, 8))  # the first 8 numbers before "-"
  a_2 = as.integer(substr(MRN, 10,10))  # the last number after "-"
  p = 17
  b = 5
  n_1 = 10^8
  n_2 = 10^1
  MRN_1 = as.integer((p*a_1+b) %% n_1)
  MRN_2 = as.integer((p*a_2+b) %% n_2)
  return(sprintf("%s-%s", toString(MRN_1), toString(MRN_2)))
}
BiopsyNum <- function(raw) {
  # Args:
  #   raw: a list of lists of words
  # Returns:
  #   the biopsy number corresponding to a MRN
  biopsynum = ""
  if (length(grep(pattern_BiopsyNum, raw)) > 0) {
    i = grep(pattern_BiopsyNum, raw)[1]  # line index
    j = grep(pattern_BiopsyNum, raw[[i]])[1]  # word index within line
    biopsynum = raw[[i]][j]
  }
  return(biopsynum)
}
Age <- function(raw) {
  # Args:
  #   raw: a list of lists of words
  # Returns:
  #   the Age (= DOB - Date Collected) at time of biopsy
  date_of_birth = ""
  date_collected = ""
  # DOB always occurrs in the same row as the MRN
  date_of_birth_index = grep(pattern_Date, raw[[1]])[1]  # index of MM/DD/YY in the first line of raw
  date_of_birth = raw[[1]][date_of_birth_index]
  # Date Collected : MM/DD/YY may not occurr in the second row
  rows = grep(pattern_Date, raw)
  for (i in 1:length(rows)) {
    cols = grep(pattern_Date, raw[[rows[i]]])
    for (j in 1:length(cols)) {
      if (raw[[rows[i]]][cols[j]-1] == ":" && raw[[rows[i]]][cols[j]-2] == "Collected") {
        date_collected = raw[[rows[i]]][cols[j]]
      }
    }
  }
  if (date_of_birth == "" || date_collected == "") {
    return("")
  } else {
    # compute age = date_collected - date_of_birth
    # might need to change if current year is not 2015
    yy_1 = as.numeric(substr(date_of_birth, 7, 8))
    yy_2 = as.numeric(substr(date_collected, 7, 8))
    if (yy_1 > 15) {
      yy_1 = 1900+yy_1
    } else {
      yy_1 = 2000+yy_1
    }
    if (yy_2 > 15) {
      yy_2 = 1900+yy_2
    } else {
      yy_2 = 2000+yy_2
    }
    return(yy_2 - yy_1)
  }
}
Sex <- function(raw) {
  # Args:
  #   raw: a list of lists of words
  # Returns:
  #   the sex of the patient
  sex_ind = grep("[[:upper:]]{1}", raw[[1]])  # indices of uppercase letters in the first row of raw
  sex = ""
  # ensure that sex has a MM/DD/YY before it
  i = 1
  while ((sex == "") & (i <= length(sex_ind))) {
    if (grepl(pattern_Date, raw[[1]][sex_ind[i]-1])) {
      sex = raw[[1]][sex_ind[i]]
    }
    i = i+1
  }
  return(sex)
}
Sex_Ind <- function(raw) {
  # Args:
  #   raw: a list of lists of words
  # Returns:
  #   the index of the sex of the patient
  sex_ind = grep("[[:upper:]]{1}", raw[[1]])  # indices of uppercase letters in the first row of raw
  ind = 0
  # ensure that sex has a MM/DD/YY before it
  i = 1
  while ((ind == 0) & (i <= length(sex_ind))) {
    if (grepl(pattern_Date, raw[[1]][sex_ind[i]-1])) {
      ind = sex_ind[i]
    }
    i = i+1
  }
  return(ind)
}
Loc <- function(raw) {
  # Args:
  #   raw: a list of lists of words
  # Returns:
  #   loc immediately after sex
  loc_ind = Sex_Ind(raw)+1
  loc = ""
  if (grepl(pattern_Loc, raw[[1]][loc_ind])) {
    loc = raw[[1]][loc_ind]
  }
  return(loc)
}
Physician <- function(raw) {
  # Args:
  #   raw: a list of lists of words
  # Returns:
  #   the requesting physician in form last_name, first_initial
  if (Loc(raw) == "") {
    phys_ind = Sex_Ind(raw)+1
  } else {
    phys_ind = Sex_Ind(raw)+2
  }
  physician = raw[[1]][phys_ind]
  physician = paste(physician, substr(raw[[1]][phys_ind+1], 1, 1), sep = " ")
  return(physician)
}
#NextWord <- function(i, j, raw) {
  # Args:
  #   i: row of raw
  #   j: col of raw
  #   raw: text of the form raw[[i]][j]
  # Returns:
  #   list (row, col) of the "next word" accounting for line breaks
#  row = i
#  col = j
#  if (i == length(raw) && j == length(raw[[i]])) {
#    return("ERROR")
#  }
#  if (j == length(raw[[i]])) {
#    row = i+1
#    col = 1
#  } else {
#    col = j+1
#  }
#  pair = c(row, col)
#  return(pair)
#}
Diagnosis <- function(raw) {
  # Args:
  #   raw: a list of lists of words
  # Returns:
  #   for now either "Melanoma", "Melanoma-in-situ", or "Residual-malignant-melanoma".
  diagnosis = ""
  if (length(grep("313993.00", raw)) > 0) {
    start = grep("313993.00", raw)[1]  # start of diagnosis section
    end = length(raw)  # all the way to the end
    # reassign raw to only include the diagnosis section
    raw = raw[start:end]
    # check for "obviously melanoma"
    if (grepl("Melanoma", raw[[3]][1])) {  # 3rd line of diagnosis section starts with "Melanoma"
      diagnosis = "Melanoma"
    }
    # check around each occurrence of "melanoma"
    row = grep("melanoma", raw, ignore.case = TRUE)  # rows where "melanoma" occurrs in diagnosis
    if (length(row) > 0) {
      diagnosis = "Melanoma"
      for (i in 1:length(row)) {
        col = grep("melanoma", raw[[row[i]]], ignore.case = TRUE)
        for (j in 1:length(col)) {
          # here raw[[row[i]]][col[j]] = "melanoma" or "Melanoma" etc.
          # so here is the proper place to check for indicators...
          # CAUTION: line breaks :P
          # ALSO CAUTION: naive thus far...
          if (grepl("melanoma-in-situ", raw[[row[i]]][col[j]], ignore.case = TRUE)) {
            diagnosis = "Melanoma-in-situ"
          }
          if (grepl("situ", raw[[row[i]]][col[j]+1], ignore.case = TRUE)) {
            diagnosis = "Melanoma-in-situ"
          }
          if (grepl("situ", raw[[row[i]]][col[j]+2], ignore.case = TRUE)) {
            diagnosis = "Melanoma-in-situ"
          }
        }
      }
    }
  }
  return(diagnosis)
}
Specimen <- function(raw) {
  # Args:
  #   raw: a list of lists of words
  # Returns:
  #   all "Skin, ..." specimens
  specimen = ""
  return(specimen)
}
Histology <- function(raw) {
  # Args:
  #   raw: a list of lists of words
  # Returns:
  #   the first Histology occurring in raw
  Histology = ""
  return(Histology)
}
Ulceration <- function(raw) {
  # Args:
  #   raw: a list of lists of words
  # Returns:
  #   the first Ulceration occurring in raw
  Ulceration = ""
  return(Ulceration)
}
Depth <- function(raw) {
  # Args:
  #   raw: a list of lists of words
  # Returns:
  #   the first Depth occurring in raw
  Depth = ""
  return(Depth)
}
Regression <- function(raw) {
  # Args:
  #   raw: a list of lists of words
  # Returns:
  #   the first Regression occurring in raw
  Regression = ""
  return(Regression)
}
SunDam <- function(raw) {
  # Args:
  #   raw: a list of lists of words
  # Returns:
  #   the first SunDam occurring in raw
  SunDam = ""
  return(SunDam)
}
PeriMarg <- function(raw) {
  # Args:
  #   raw: a list of lists of words
  # Returns:
  #   the first PeriMarg occurring in raw
  PeriMarg = ""
  return(PeriMarg)
}
DeepMarg <- function(raw) {
  # Args:
  #   raw: a list of lists of words
  # Returns:
  #   the first DeepMarg occurring in raw
  DeepMarg = ""
  return(DeepMarg)
}
TNM <- function(raw) {
  # Args:
  #   raw: a list of lists of words
  # Returns:
  #   the first TNM occurring in raw
  TNM = ""
  return(TNM)
}

# main function
parsefile <- function(txtfilename) {
  # Args:
  #   txtfilename: the name (as a string) of a (terrible) txt file
  # Returns:
  #   data.frame
  rows = scan(txtfilename, character(), sep = "\n", strip.white = TRUE)
  # raw[[i]][j] is the jth "word" in the ith nonempty row of sample.txt
  raw = list()
  for (i in 1:length(rows)) {
    Sys.setlocale('LC_ALL','C')  # suppress warning "input string 1 is invalid in this locale"
    raw[i] = strsplit(rows[i], split = " +")  # split each row by varying "whitespace"
  }
  # create empty data.frame
  df = data.frame(MRN = c(), BiopsyNum = c(), Age = c(), Sex = c(), Loc = c(), Physician = c(), Diagnosis = c(), Specimen = c(), Histology = c(), Ulceration = c(), Depth = c(), Regression = c(), SunDam = c(), PeriMarg = c(), DeepMarg = c(), TNM = c())
  # partition txt file by MRN
  index = grep(pattern_MRN, raw)  # (line) indices of med numbers
  # we only want MRNs with SD-XX-XXXXX in the same row
  new_index = list()
  for (i in 1:length(index)) {
    if (length(grep(pattern_BiopsyNum, raw[index[i]])) > 0) {
      if (length(new_index) == 0) {
        new_index = c(index[i])
      } else {
        new_index = c(new_index, index[i])
      }
    }
  }
  index = new_index
  num_MRN = length(index)  # number of MRNs in text file
  for (i in 1:num_MRN) {
    print(i)
    print(raw[[index[i]]][1])
    new_MRN = EncodeMRN(raw[[index[i]]][1])
    # create raw_temp a subset of raw to search
    raw_start = index[i]
    if (i == num_MRN) {
      raw_end = length(raw)  # then this is the last MRN so go to the end of the text file
    } else {
      raw_end = index[i+1] - 1  # then raw_end should be the line before the next MRN
    }
    raw_temp = raw[raw_start:raw_end]
    # call functions to search raw_temp for specific information
    new_BiopsyNum = BiopsyNum(raw_temp)
    print(new_BiopsyNum)
    new_Age = Age(raw_temp)
    new_Sex = Sex(raw_temp)
    new_Loc = Loc(raw_temp)
    new_Physician = Physician(raw_temp)
    new_Diagnosis = Diagnosis(raw_temp)
    new_Specimen = Specimen(raw_temp)
    if (new_Diagnosis == "Melanoma") {
      new_Histology = Histology(raw_temp)
      new_Ulceration = Ulceration(raw_temp)
      new_Depth = Depth(raw_temp)
      new_Regression = Regression(raw_temp)
      new_SunDam = SunDam(raw_temp)
      new_PeriMarg = PeriMarg(raw_temp)
      new_DeepMarg = DeepMarg(raw_temp)
      new_TNM = TNM(raw_temp)
    } else {
      new_Histology = ""
      new_Ulceration = ""
      new_Depth = ""
      new_Regression = ""
      new_SunDam = ""
      new_PeriMarg = ""
      new_DeepMarg = ""
      new_TNM = ""
    }
    # create new row to put in data.frame
    new_df = data.frame(MRN = c(new_MRN), BiopsyNum = c(new_BiopsyNum), Age = c(new_Age), Sex = c(new_Sex), Loc = c(new_Loc), Physician = c(new_Physician), Diagnosis = c(new_Diagnosis), Specimen = c(new_Specimen), Histology = c(new_Histology), Ulceration = c(new_Ulceration), Depth = c(new_Depth), Regression = c(new_Regression), SunDam = c(new_SunDam), PeriMarg = c(new_PeriMarg), DeepMarg = c(new_DeepMarg), TNM = c(new_TNM))
    # append the new row to the global data.frame
    df = rbind(df, new_df)
    # update raw_start
  }
  return(df)
}
