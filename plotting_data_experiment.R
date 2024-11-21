###############################################################################
# Loading data experiment ----
###############################################################################
plotting_data_experiment <- function(
    data,
    parameters,
    model=NULL,
    model.sto=NULL
    ) {
  
  # jpeg(("FREP1.jpeg"), width = 7, height = 7, units = 'in', res = 250)
  tiff(("FREP1.tiff"), width = 7, height = 7, units = 'in', res = 200)
  
  # aux for the length of the x-axis
  aux.x      = length(model$LD)-1
  aux.x.data = length(data[,c("data_screened_cage1")])-1
  
  # plotting data stratified by sex
  par(oma = c(0, 0, 0, 0), # outer margins
      mai = c(0.32, 0.4, 0.3, 0.2), # inner margins for each panel
      mgp = c(1.5, 0.5, 0))
  layout(matrix(c(1,2,3,4,5,6,7,8,9), nrow = 3, ncol = 3, byrow = F))
  
  ##############################################################################
  # LDonly ----
  plot(
    0,100*data[1,c("data_phenotype_LDonly_cage1")]/data$data_screened_cage1[1],
    ylim = 100*c(0,1),
    xlim = c(0,aux.x),
    ylab = "Homozygous RFP-gRNA-Q (%)",
    xlab = "Generations",
    main = "",
    pch  = 19,
    col = "magenta"
  )
  
  mtext("A", side = 3, adj = -0.2, line = 1.0)
  
  if (!is.null(model.sto)) {
    for (j in 1:length(model.sto)) {
      lines(c(0:aux.x),100*model.sto[[j]]$LDonly,type = "l",col = "lightblue")  
    }
  }  
  
  lines(c(0:aux.x.data),100*data[,c("data_phenotype_LDonly_cage1")]/data$data_screened_cage1,type = "l",col = "magenta")
  lines(c(0:aux.x.data),100*data[,c("data_phenotype_LDonly_cage2")]/data$data_screened_cage2,type = "l",col = "magenta")
  lines(c(0:aux.x.data),100*data[,c("data_phenotype_LDonly_cage3")]/data$data_screened_cage3,type = "l",col = "magenta")
  points(
    0,100*data[1,c("data_phenotype_LDonly_cage1")]/data$data_screened_cage1[1],
    pch  = 19,
    col = "magenta"
  )
  
  if (!is.null(model)) {
    lines(c(0:aux.x),100*model$LDonly,type = "l",col = "blue")
  }
  
  legend(
    "topright",
    c("experimental data","model sim. (deterministic)","model sim. (stochastic)"),
    lty = c(1,1),
    col = c("magenta","blue","lightblue"),
    bty = "n"
  )
  
  ##############################################################################
  # GFPLonly ----
  plot(
    0,100*data[1,c("data_phenotype_GFPLonly_cage1")]/data$data_screened_cage1[1],
    ylim = 100*c(0,1),
    xlim = c(0,aux.x),
    ylab = "Homozygous GFP (%)",
    xlab = "",
    main = "",
    pch  = 19,
    col = "magenta"
  )
  
  mtext("B", side = 3, adj = -0.2, line = 1.0)
  
  if (!is.null(model.sto)) {
    for (j in 1:length(model.sto)) {
      lines(c(0:aux.x),100*model.sto[[j]]$GFPLonly,type = "l",col = "lightblue")  
    }
  }  
  
  lines(c(0:aux.x.data),100*data[,c("data_phenotype_GFPLonly_cage1")]/data$data_screened_cage1,type = "l",col = "magenta")
  lines(c(0:aux.x.data),100*data[,c("data_phenotype_GFPLonly_cage2")]/data$data_screened_cage2,type = "l",col = "magenta")
  lines(c(0:aux.x.data),100*data[,c("data_phenotype_GFPLonly_cage3")]/data$data_screened_cage3,type = "l",col = "magenta")
  points(
    0,100*data[1,c("data_phenotype_GFPLonly_cage3")]/data$data_screened_cage1[1],
    pch  = 19,
    col = "magenta"
  )
  
  if (!is.null(model)) {
    lines(c(0:aux.x),100*model$GFPLonly,type = "l",col = "blue")
  }
  
  ##############################################################################
  # LDGFPL ----
  plot(
    0,100*data[1,c("data_phenotype_LDGFPL_cage1")]/data$data_screened_cage1[1],
    ylim = 100*c(0,1),
    xlim = c(0,aux.x),
    ylab = "RFP-gRNA-Q/GFP (%)",
    xlab = "",
    main = "",
    pch  = 19,
    col = "magenta"
  )
  
  mtext("C", side = 3, adj = -0.2, line = 1.0)
  
  if (!is.null(model.sto)) {
    for (j in 1:length(model.sto)) {
      lines(c(0:aux.x),100*model.sto[[j]]$LDGFPL,type = "l",col = "lightblue")  
    }
  }  
  
  lines(c(0:aux.x.data),100*data[,c("data_phenotype_LDGFPL_cage1")]/data$data_screened_cage1,type = "l",col = "magenta")
  lines(c(0:aux.x.data),100*data[,c("data_phenotype_LDGFPL_cage2")]/data$data_screened_cage2,type = "l",col = "magenta")
  lines(c(0:aux.x.data),100*data[,c("data_phenotype_LDGFPL_cage3")]/data$data_screened_cage3,type = "l",col = "magenta")
  points(
    0,100*data[1,c("data_phenotype_LDGFPL_cage1")]/data$data_screened_cage1[1],
    pch  = 19,
    col = "magenta"
  )
  
  if (!is.null(model)) {
    lines(c(0:aux.x),100*model$LDGFPL,type = "l",col = "blue")
  }

  ##############################################################################
  # L224 ----
  plot(
    0,100*data[1,c("data_allelefreq_L224_cage1")]/data$data_genotyped_cage1[1],
    ylim = 100*c(0,1),
    xlim = c(0,aux.x),
    ylab = "L224 allele freq (%)",
    xlab = "",
    main = "",
    pch  = 19,
    col = "magenta"
  )
  
  mtext("D", side = 3, adj = -0.2, line = 1.0)
  
  if (!is.null(model.sto)) {
    for (j in 1:length(model.sto)) {
      lines(c(0:aux.x),100*model.sto[[j]]$L224,type = "l",col = "lightblue")  
    }
  }  
  
  lines(c(0:aux.x.data),100*data[,c("data_allelefreq_L224_cage1")]/data$data_genotyped_cage1,type = "l",col = "magenta")
  lines(c(0:aux.x.data),100*data[,c("data_allelefreq_L224_cage2")]/data$data_genotyped_cage2,type = "l",col = "magenta")
  lines(c(0:aux.x.data),100*data[,c("data_allelefreq_L224_cage3")]/data$data_genotyped_cage3,type = "l",col = "magenta")
  points(
    0,100*data[1,c("data_allelefreq_L224_cage1")]/data$data_genotyped_cage1[1],
    pch  = 19,
    col = "magenta"
  )
  
  if (!is.null(model)) {
    lines(c(0:aux.x),100*model$L224,type = "l",col = "blue")
  }
  
  ##############################################################################
  # Q224 ----
  plot(
    0,100*data[1,c("data_allelefreq_Q224_cage1")]/data$data_genotyped_cage1[1],
    ylim = 100*c(0,1),
    xlim = c(0,aux.x),
    ylab = "Q224 allele freq (%)",
    xlab = "",
    main = "",
    pch  = 19,
    col = "magenta"
  )
  
  mtext("E", side = 3, adj = -0.2, line = 1.0)
  
  if (!is.null(model.sto)) {
    for (j in 1:length(model.sto)) {
      lines(c(0:aux.x),100*model.sto[[j]]$Q224,type = "l",col = "lightblue")  
    }
  }  
  
  lines(c(0:aux.x.data),100*data[,c("data_allelefreq_Q224_cage1")]/data$data_genotyped_cage1,type = "l",col = "magenta")
  lines(c(0:aux.x.data),100*data[,c("data_allelefreq_Q224_cage2")]/data$data_genotyped_cage2,type = "l",col = "magenta")
  lines(c(0:aux.x.data),100*data[,c("data_allelefreq_Q224_cage3")]/data$data_genotyped_cage3,type = "l",col = "magenta")
  points(
    0,100*data[1,c("data_allelefreq_Q224_cage1")]/data$data_genotyped_cage1[1],
    pch  = 19,
    col = "magenta"
  )
  
  if (!is.null(model)) {
    lines(c(0:aux.x),100*model$Q224,type = "l",col = "blue")
  }
  
  ##############################################################################
  # NHEJ ----
  plot(
    0,100*data[1,c("data_allelefreq_NHEJ_cage1")]/data$data_genotyped_cage1[1],
    ylim = 100*c(0,1),
    xlim = c(0,aux.x),
    ylab = "NHEJ allele freq (%)",
    xlab = "",
    main = "",
    pch  = 19,
    col = "magenta"
  )
  
  mtext("F", side = 3, adj = -0.2, line = 1.0)
  
  if (!is.null(model.sto)) {
    for (j in 1:length(model.sto)) {
      lines(c(0:aux.x),100*model.sto[[j]]$NHEJ,type = "l",col = "lightblue")  
    }
  }  
  
  lines(c(0:aux.x.data),100*data[,c("data_allelefreq_NHEJ_cage1")]/data$data_genotyped_cage1,type = "l",col = "magenta")
  lines(c(0:aux.x.data),100*data[,c("data_allelefreq_NHEJ_cage2")]/data$data_genotyped_cage2,type = "l",col = "magenta")
  lines(c(0:aux.x.data),100*data[,c("data_allelefreq_NHEJ_cage3")]/data$data_genotyped_cage3,type = "l",col = "magenta")
  points(
    0,100*data[1,c("data_allelefreq_NHEJ_cage1")]/data$data_genotyped_cage1[1],
    pch  = 19,
    col = "magenta"
  )
  
  if (!is.null(model)) {
    lines(c(0:aux.x),100*model$NHEJ,type = "l",col = "blue")
  }
  
  # ##############################################################################
  # LDonly + LDGFPL ----
  plot(
    0,100*(data[1,c("data_phenotype_LDonly_cage1")]+data[1,c("data_phenotype_LDGFPL_cage1")])/data$data_screened_cage1[1],
    ylim = 100*c(0,1),
    xlim = c(0,aux.x),
    ylab = "Mosquitoes with linked drive (%)",
    xlab = "",
    main = "",
    pch  = 19,
    col = "magenta"
  )
  
  mtext("G", side = 3, adj = -0.2, line = 1.0)
  
  if (!is.null(model.sto)) {
    for (j in 1:length(model.sto)) {
      lines(c(0:aux.x),100*(model.sto[[j]]$LDonly+model.sto[[j]]$LDGFPL),type = "l",col = "lightblue")  
    }
  }  
  
  lines(c(0:aux.x.data),100*(data[,c("data_phenotype_LDonly_cage1")]+data[,c("data_phenotype_LDGFPL_cage1")])/data$data_screened_cage1,type = "l",col = "magenta")
  lines(c(0:aux.x.data),100*(data[,c("data_phenotype_LDonly_cage2")]+data[,c("data_phenotype_LDGFPL_cage2")])/data$data_screened_cage2,type = "l",col = "magenta")
  lines(c(0:aux.x.data),100*(data[,c("data_phenotype_LDonly_cage3")]+data[,c("data_phenotype_LDGFPL_cage3")])/data$data_screened_cage3,type = "l",col = "magenta")
  points(
    0,100*(data[1,c("data_phenotype_LDonly_cage1")]+data[1,c("data_phenotype_LDGFPL_cage1")])/data$data_screened_cage1[1],
    pch  = 19,
    col = "magenta"
  )
  
  if (!is.null(model)) {
    lines(c(0:aux.x),100*(model$LDonly+model$LDGFPL),type = "l",col = "blue")
  }
  
  # ##############################################################################
  # LDonly + LDGFPL ----
  plot(
    0,100*(data[1,c("data_phenotype_GFPLonly_cage1")]+data[1,c("data_phenotype_LDGFPL_cage1")])/data$data_screened_cage1[1],
    ylim = 100*c(0,1),
    xlim = c(0,aux.x),
    ylab = "Mosquitoes with GFP (%)",
    xlab = "",
    main = "",
    pch  = 19,
    col = "magenta"
  )
  
  mtext("H", side = 3, adj = -0.2, line = 1.0)
  
  if (!is.null(model.sto)) {
    for (j in 1:length(model.sto)) {
      lines(c(0:aux.x),100*(model.sto[[j]]$GFPLonly+model.sto[[j]]$LDGFPL),type = "l",col = "lightblue")  
    }
  }  
  
  lines(c(0:aux.x.data),100*(data[,c("data_phenotype_GFPLonly_cage1")]+data[,c("data_phenotype_LDGFPL_cage1")])/data$data_screened_cage1,type = "l",col = "magenta")
  lines(c(0:aux.x.data),100*(data[,c("data_phenotype_GFPLonly_cage2")]+data[,c("data_phenotype_LDGFPL_cage2")])/data$data_screened_cage2,type = "l",col = "magenta")
  lines(c(0:aux.x.data),100*(data[,c("data_phenotype_GFPLonly_cage3")]+data[,c("data_phenotype_LDGFPL_cage3")])/data$data_screened_cage3,type = "l",col = "magenta")
  points(
    0,100*(data[1,c("data_phenotype_GFPLonly_cage1")]+data[1,c("data_phenotype_LDGFPL_cage1")])/data$data_screened_cage1[1],
    pch  = 19,
    col = "magenta"
  )
  
  if (!is.null(model)) {
    lines(c(0:aux.x),100*(model$GFPLonly+model$LDGFPL),type = "l",col = "blue")
  }
  
  
  # ##############################################################################
  # V9 ----
  plot(
    0,100*data[1,c("data_phenotype_V9_cage1")]/data$data_screened_cage1[1],
    ylim = 100*c(0,1),
    xlim = c(0,aux.x),
    ylab = "vasa-Cas9 (%)",
    xlab = "",
    main = "",
    pch  = 19,
    col = "magenta"
  )
  
  mtext("I", side = 3, adj = -0.2, line = 1.0)
  
  if (!is.null(model.sto)) {
    for (j in 1:length(model.sto)) {
      lines(c(0:aux.x),100*model.sto[[j]]$V9,type = "l",col = "lightblue")  
    }
  }  
  
  lines(c(0:aux.x.data),100*data[,c("data_phenotype_V9_cage1")]/data$data_screened_cage1,type = "l",col = "magenta")
  lines(c(0:aux.x.data),100*data[,c("data_phenotype_V9_cage2")]/data$data_screened_cage2,type = "l",col = "magenta")
  lines(c(0:aux.x.data),100*data[,c("data_phenotype_V9_cage3")]/data$data_screened_cage3,type = "l",col = "magenta")
  points(
    0,100*data[1,c("data_phenotype_V9_cage1")]/data$data_screened_cage1[1],
    pch  = 19,
    col = "magenta"
  )
  
  if (!is.null(model)) {
    lines(c(0:aux.x),100*model$V9,type = "l",col = "blue")
  }
  
  
  dev.off()
  return()
}
