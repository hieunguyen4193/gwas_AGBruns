

QC.ChiTest100_plate1 <- unlist(lapply(snp.full.info$ChiTest100_plate1, function(x){
  if (x >= thres_hwe){
    return("PASS")
  } else {
    return("FAIL")
  }
}))

snp.full.info[["QC.ChiTest100_plate1"]] <- QC.ChiTest100_plate1

QC.Call.Freq_plate1 <- unlist(lapply(snp.full.info$Call.Freq_plate1, function(x){
  if (x >= thres_call.freq){
    return("PASS")
  } else {
    return("FAIL")
  }
}))

snp.full.info[["QC.Call.Freq_plate1"]] <- QC.Call.Freq_plate1

QC.Minor.Freq_plate1 <- unlist(lapply(snp.full.info$Minor.Freq_plate1, function(x){
  if (x >= thres_maf){
    return("PASS")
  } else {
    return("FAIL")
  }
}))

snp.full.info[["QC.Minor.Freq_plate1"]] <- QC.Minor.Freq_plate1

QC.ChiTest100_plate2 <- unlist(lapply(snp.full.info$ChiTest100_plate2, function(x){
  if (x >= thres_hwe){
    return("PASS")
  } else {
    return("FAIL")
  }
}))

snp.full.info[["QC.ChiTest100_plate2"]] <- QC.ChiTest100_plate2

QC.Call.Freq_plate2 <- unlist(lapply(snp.full.info$Call.Freq_plate2, function(x){
  if (x >= thres_call.freq){
    return("PASS")
  } else {
    return("FAIL")
  }
}))

snp.full.info[["QC.Call.Freq_plate2"]] <- QC.Call.Freq_plate2

QC.Minor.Freq_plate2 <- unlist(lapply(snp.full.info$Minor.Freq_plate2, function(x){
  if (x >= thres_maf){
    return("PASS")
  } else {
    return("FAIL")
  }
}))

snp.full.info[["QC.Minor.Freq_plate2"]] <- QC.Minor.Freq_plate2