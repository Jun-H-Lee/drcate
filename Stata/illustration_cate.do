*** Illustration code for lowcate.ado
clear

webuse set http://www.stata-press.com/data/r13
webuse cattaneo2

drcate bweight alcohol mage medu fage, ps("logit") bwidth(.781) graph("on") ci("on") ate("on")

ereturn list
/* varlist [1] : dependent variable
   varlist [2] : treatment variable
   varlist [3] : covariate of interest
   varlist [4::end] : other covariates
   ps : "logit" or "probit"
   bwidth : user-specified bandwidth # defalut = Plugin bandwidth of X */
