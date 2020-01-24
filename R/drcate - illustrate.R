### Illustration for function "drcate"

require("locpol")
require("KernSmooth")
require("ggplot2")
require("haven")

source("drcate.R")

df = read_dta("http://www.stata-press.com/data/r13/cattaneo2.dta")

drcate(Y = df$bweight, X = df$mage, V = cbind(df$medu, df$fage), D = df$alcohol, alpha = 0.05, ps = "logit")