########################################################################################################################

rm(list = ls())                                                 
library(rvif)      # for the detection of multicolinearity and Klein and Goldberger data
library(faraway)   # for Cheddar data
library(xtable)    # to export tables in LaTeX format
library(jtools)    # for 'export_summs'

########################################################################################################################

# 1.- The data

  data(cheddar)
  head(cheddar, n=5)

# 2.- OLS estimation
  
  taste = cheddar$taste
  acetic = cheddar$Acetic
  h2s = cheddar$H2S
  lactic = cheddar$Lactic
  model = lm(taste~acetic+h2s+lactic)
  summary(model)
      
# 3.- Detection of multicolinearity

  X = model.matrix(model)  
  RdetR(X)
  CNs(X)
  ki(X)
  vifs = c(NA, VIF(X))
  cvs = c(NA, CVs(X))
  kis1 = ki(X)[[1]]
  kis2 = c(NA, ki(X)[[2]])
  kis3 = c(NA, ki(X)[[3]])
  table1 = data.frame(vifs, cvs, kis1, kis2, kis3)
  colnames(table1) = c("VIF", "CV", "k_i", "% essential", "% non-essential")
  rownames(table1) = c("Intercept", "Acetic", "H2S", "Lactic")
  table1
  
  # % of essential and non-essential multicollinearity for centered data
  
  Xcentered <- cbind(X[, 1], scale(X[, -1], center = TRUE, scale = FALSE))
  ki(Xcentered)

# 4.- Raise regression: stewart index calculation
  
  # 'acetic' is the variable to be raised since is the variable with the highest VIF and k_i 
  
  model_aux = lm(acetic~h2s+lactic)  # expression (2)
  e = model_aux$residuals
  
  top = 0.7
  step = 0.02
  landas = seq(0, top, step)
  Xraised = X
  k1 = numeric()
  k2 = numeric()
  k3 = numeric()
  k4 = numeric()
  i = 0
  for (landa in landas){
    Xraised[,2] = X[,2] + landa*e                                   # model design matrix (6)
    kis_raised = ki(Xraised)
    i = i + 1
    k1[i] = kis_raised[[1]][1]
    k2[i] = kis_raised[[1]][2]
    k3[i] = kis_raised[[1]][3]
    k4[i] = kis_raised[[1]][4]
  }
  ki_raised = data.frame(landas, k1, k2, k3, k4)
  ki_raised 
  table2 <- xtable(ki_raised, digits = 2)
  print(table2, type = "latex", file = "example2 ki_raised.tex", include.rownames = FALSE)
  
  # joint representation
  
  pdf("example2 graph-joint.pdf") 
    matplot(landas, ki_raised[,-1], type = "l", lwd = 2, pch = 19, lty = 1, col = c("blue", "cyan", "darkgreen", "brown"), xlab = expression(lambda), ylab = "Stewart index")
    points(c(0, 0, 0, 0), kis1, pch = 19, col = c("blue", "cyan", "darkgreen", "brown"), cex = 1.2)
    legend("topright", legend = expression(k[1], k[2], k[3], k[4]), col = c("blue", "cyan", "darkgreen", "brown"), lty = 1, lwd = 2)
    abline(h=1, col="black")
    abline(h=64.4293, col="black", lty=2)
    abline(h=134.0091, col="black", lty=3)
  dev.off()
      
  # separate representation
  
  pdf("example2 graph-separate.pdf") 
    par(mfrow = c(2, 2))
    for (i in 2:5) {
      plot(landas, ki_raised[, i],
         type = "l", lwd = 2,
         main = rownames(table1)[i-1],
         xlab = expression(lambda),
         ylab = "Stewart index",
         col="lightblue")
      points(0, kis1[i-1], pch = 19, col = "blue", cex = 1.2)
      abline(h=1, col="black")
      abline(h=64.4293, col="black", lty=2)
      abline(h=134.0091, col="black", lty=3)
    }
    par(mfrow = c(1, 1))
  dev.off()
  
  # asymptotic behavior
  
  top = 150
  step = 1
  landas = seq(0, top, step)
  Xraised = X
  k1 = numeric()
  k2 = numeric()
  k3 = numeric()
  k4 = numeric()
  i = 0
  for (landa in landas){
    Xraised[,2] = X[,2] + landa*e                                   # model design matrix (6)
    kis_raised = ki(Xraised)
    i = i + 1
    k1[i] = kis_raised[[1]][1]
    k2[i] = kis_raised[[1]][2]
    k3[i] = kis_raised[[1]][3]
    k4[i] = kis_raised[[1]][4]
  }
  ki_raised = data.frame(landas, k1, k2, k3, k4)
  
  # horizontal asymptotes
  
  Xasymptote = X[,-2]
  ki_a = ki(Xasymptote)
  ki_a[[1]]
  
  # joint representation
  
  pdf("example2 asymptotic behavior.pdf") 
      matplot(landas, ki_raised[,-1], type = "l", lwd = 2, pch = 19, lty = 1, col = c("blue", "cyan", "darkgreen", "brown"), xlab = expression(lambda), ylab = "Stewart index")
      points(c(0, 0, 0, 0), kis1, pch = 19, col = c("blue", "cyan", "darkgreen", "brown"), cex = 1.2)
      legend("topright", legend = expression(k[1], k[2], k[3], k[4]), col = c("blue", "cyan", "darkgreen", "brown"), lty = 1, lwd = 2)
      abline(h=1, col="black")
      abline(h=ki_a[[1]][1], col="black", lty=2)
      abline(h=ki_a[[1]][2], col="black", lty=2)
      abline(h=ki_a[[1]][3], col="black", lty=2)
  dev.off()
  
# 5.- Raise regression: estimation for the selected raising factors
  
  # for landa = 0.16 all k_{i}^{2} < 134.0091
  
    landa = 0.16
    Xraised = X
    Xraised[,2] = X[,2] + landa*e
    acetic = Xraised[,2]
    model_raised1 = lm(taste~acetic+h2s+lactic)
    summary(model_raised1)
    
      RdetR(Xraised)
      CNs(Xraised)
      ki(Xraised)
      vifs = c(NA, VIF(Xraised))
      cvs = c(NA, CVs(Xraised))
      kis = ki(Xraised)[[1]]
      table3 = data.frame(vifs, cvs, kis)
      colnames(table3) = c("VIF", "CV", "k_i")
      rownames(table3) = c("Intercept", "Wage income", "Non farm income", "Farm income")
      table3
    
  # for landa = 0.68 all k_{i}^{2} < 64.4293
      
    landa = 0.68
    Xraised = X
    Xraised[,2] = X[,2] + landa*e
    acetic = Xraised[,2]
    model_raised2 = lm(taste~acetic+h2s+lactic)
    summary(model_raised2)
      
      RdetR(Xraised)
      CNs(Xraised)
      ki(Xraised)
      vifs = c(NA, VIF(Xraised))
      cvs = c(NA, CVs(Xraised))
      kis = ki(Xraised)[[1]]
      table4 = data.frame(vifs, cvs, kis)
      colnames(table4) = c("VIF", "CV", "k_i")
      rownames(table4) = c("Intercept", "Wage income", "Non farm income", "Farm income")
      table4

  # for landa = 150 
      
      landa = 150
      Xraised = X
      Xraised[,2] = X[,2] + landa*e
      acetic = Xraised[,2]
      model_raised3 = lm(taste~acetic+h2s+lactic)
      summary(model_raised3)
      
      RdetR(Xraised)
      CNs(Xraised)
      ki(Xraised)
      vifs = c(NA, VIF(Xraised))
      cvs = c(NA, CVs(Xraised))
      kis = ki(Xraised)[[1]]
      table5 = data.frame(vifs, cvs, kis)
      colnames(table5) = c("VIF", "CV", "k_i")
      rownames(table5) = c("Intercept", "Wage income", "Non farm income", "Farm income")
      table5      
      
########################################################################################################################  
      
# centered model
      
      acetic = cheddar$Acetic - mean(cheddar$Acetic)
      centered = lm(taste~acetic+h2s+lactic)
      summary(centered)
      
      Xc = model.matrix(centered)  
      RdetR(Xc)
      CNs(Xc)
      ki(Xc)
      vifs = c(NA, VIF(Xc))
      cvs = c(NA, CVs(Xc))
      kis1 = ki(Xc)[[1]]
      kis2 = c(NA, ki(Xc)[[2]])
      kis3 = c(NA, ki(Xc)[[3]])
      table6 = data.frame(vifs, cvs, kis1, kis2, kis3)
      colnames(table6) = c("VIF", "CV", "k_i", "% essential", "% non-essential")
      rownames(table6) = c("Intercept", "Acetic", "H2S", "Lactic")
      table6
            
########################################################################################################################        

# summary
  
  export_summs(model, model_raised1, model_raised2, model_raised3, centered, model.names = c("OLS", "Raised (1)", "Raised (2)", "Raised (3)", "Centrado"), number_format = "%.3f") 
  table1
  table3
  table4
  table5
  table6

########################################################################################################################