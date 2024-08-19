.onAttach =
  function(libname, pkgname){
    msg = paste0("Loading required package: ", pkgname,"\n\n",
                 "If you are using DEprot in your work please cite:\n\n",
                 "HDAC1 and PRC2 mediate combinatorial control in SPI1/PU.1-dependent gene repression in murine erythroleukaemia.\n",
                 "Gregoricchio S. et al., Nucleic Acid Research (2022)\n",
                 "doi: 10.1093/nar/gkac613")
    packageStartupMessage(msg)
  }
