setwd("../../../Dropbox/Smoking MR/")
data = read.sas7bdat("data/k06113_waldenbergerf4_tra290713.sas7bdat")

metabolites = colnames(data)[27:163]
confounders = c("utalter","ucsex","utbmi", "my.alkkon")
##
# variable explaination:
#   utgewi: weight
#   utsmkreg: smoking status (1:5, different from utcigreg in that it set different classes of fomer smoker)
#   utlipi: lipid lowering drug intake.
#   utphys: active or inactive in physical cativities
#   utglukfast_a: fasting plasma glucose (mg/dl)
#   utglukrand_a: random bblood glucose (mg/dl)
#   uh_ins: insulin in serum (ulU/ml)
#   uc101: Pregancy (1 (pregant),2 (no),3 (maybe not)

#Exclusions
# •  Individuals on lipid-lowering medication (utlipid), on hypertensive treatment (utantihy) or pregnant (uc101)
data = subset(data, subset = utlipi==2&utantihy==2&uc101!=1)
# •	Individuals with extreme metabolite levels (± 4 standard deviations) 
data[,metabolites] = sapply(data[,metabolites], 
             function(x){
               x[x>mean(x,na.rm=T)+4*sd(x,na.rm=T)|x<mean(x,na.rm=T)-4*sd(x,na.rm=T)]=NA
               return(x)
                } 
             )
# •  Exclusion of non-European ethnicity
# •  Exclude the samples without genotype measurements.


# 5) Phenotypes
# a) Smoking status: never/former/current smoker (uctigreg)
data$my.cigreg = data$utcigreg
data$my.cigreg[which(data$utcigreg<3)] = 2
data$my.cigreg[which(data$utcigreg==3)] = 1
data$my.cigreg[which(data$utcigreg==4)] = 0
# b) Serum metabolites from high-throughput NMR spectrometer 
metabolties
# c) Potential confounding factors: age (utalter), sex (ucsex), alcohol usage, (BMI)(utbmi)
data$my.alkkon = rep(0, nrow(data))
data$my.alkkon[which(data$utalkkon >=40 & data$ucsex==1 )] = 1
data$my.alkkon[which(data$utalkkon >=20 & data$ucsex==2 )] = 1


# 6) Genotypes
# a) Extract from GWAS data, rs1051730/rs16969968 in the nicotine acetylcholine receptor gene cluster (CHRNA5-CHRNA3-CHRNB4). SNPs are in perfect linkage disequilibrium, so either could be used.
# b) Investigate frequency of minor alleles. 
table(data$rs1051730_probAA, useNA = "always")
table(data$rs1051730_probAB, useNA = "always")
table(data$rs1051730_probBB, useNA = "always")
table(data$rs16969968_probAA, useNA = "always")
table(data$rs16969968_probAB, useNA = "always")
table(data$rs16969968_probBB, useNA = "always")

data$rs1051730 = NA
data$rs1051730[which(data$rs1051730_probAA==1)] = 0
data$rs1051730[which(data$rs1051730_probAB==1)] = 1
data$rs1051730[which(data$rs1051730_probBB==1)] = 2
data$rs16969968 = NA
data$rs16969968[which(data$rs16969968_probAA==1)] = 0
data$rs16969968[which(data$rs16969968_probAB==1)] = 1
data$rs16969968[which(data$rs16969968_probBB==1)] = 2

# c) Test assumptions of Hardy-Weinberg Equilibrium e.g. using chi-squared test. 


# *** In all analyses, use 2 models: a) with adjustment for sex and age; b) with adjustment for sex, age and confounding factors listed above***
# 7) Analysis 
# a) Characteristics of study population (where available): gender, age, BMI, prevalence of overweight and obesity, systolic blood pressure, total cholesterol, HDL/LDL cholesterol, triglycerides, glucose, insulin, physical activity measures, dietary measures, alcohol intake, smoking prevalence. 
data$obes = 0
data$obes[which(data$utbmi>=25&data$utbmi<30)] = 1
data$obes[which(data$utbmi>=25&data$utbmi>=30)] = 2
feature = colnames(data)[1:25]
rst = characteristics(data[,feature], factor=as.factor(data$my.cigreg), 2, na.rm = T)


# b) Investigate observational associations between smoking status and overall metabolic traits, using linear regression with indicator variables.  
rst = NULL
for(i in metabolites){
  data$m = data[,i]
  model = lm(m ~ as.factor(my.cigreg)
          +utalter+ucsex
             #+ utbmi + my.alkkon
             #+ as.factor(obes) + utsysmm + ul_cholal+ul_tria+utglukfast_a+uh_ins+as.factor(utphact)
     , data)
  rst = rbind(rst, c(summary(model)$coef[2,],summary(model)$coef[3,]))
}
rownames(rst) = metabolites
write.csv(rst, file = "associations between smoking status and overall metabolic traits_model1.csv")


# c) Investigate associations between genotype and smoking status assuming an additive model and using logistic/multinomial regression.
model = glm( as.factor(my.cigreg) ~ rs16969968
           +utalter+ucsex
          #+ utbmi + my.alkkon
           #+ as.factor(obes) + utsysmm + ul_cholal+ul_tria+utglukfast_a+uh_ins+as.factor(utphact)
            , family=binomial
           , subset = which(data$my.cigreg!=2)
             , data)
rst = summary(model)$coef
write.csv(rst, file = "associations between genotype and smoking status_FSvsNS_model1.csv")

# d) Investigate associations between genotype and possible confounding factors stated above.
rst=NULL
for( i in c(confounders, "utalkkon")){
  data$conf = data[,i]
  model = glm( conf~ rs16969968, data=data)
  rst = rbind(rst, summary(model)$coef[2,])
}
rownames(rst) = c(confounders, "utalkkon")
write.csv(rst, "associations between rs16969968 and confounders.csv")

# e) Test for interaction between genotype and smoking status in linear regression of metabolite levels, using a likelihood ratio test.
rst = NULL
for(i in metabolites){
  data$m = data[,i]
  model.full = lm(m ~ as.factor(my.cigreg)+as.factor(my.cigreg):rs16969968
             +utalter+ucsex
             + utbmi + my.alkkon
             #+ as.factor(obes) + utsysmm + ul_cholal+ul_tria+utglukfast_a+uh_ins+as.factor(utphact)
             , data)
  model.reduced = update(model.full, .~.-+as.factor(my.cigreg):rs16969968, subset = which(!is.na(data$rs16969968)))
  tst = lrtest(model.full, model.reduced)
  rst = rbind(rst, tst[2,])
}
rownames(rst) = metabolites
write.csv(rst, file = "Test for interaction between genotype and smoking status_model1.csv")


# f) i) Stratify individuals based on smoking status: never smoker, former smoker, current smoker
# ii) Within each stratum, investigate the associations of rs1051730/rs16969968 and all metabolite traits, and estimate continuous effects using linear regression (assuming an additive model)
rst = NULL
for(i in metabolites){
  data$m = data[,i]
  model = lm(m ~ rs16969968
                  +utalter+ucsex
                  #+ utbmi + my.alkkon
                  #+ as.factor(obes) + utsysmm + ul_cholal+ul_tria+utglukfast_a+uh_ins+as.factor(utphact)
                  , subset = which(data$my.cigreg==0)
                  , data)
  rst = rbind(rst, summary(model)$coef[2,])
}
rownames(rst) = metabolites
write.csv(rst, file = "associations of rs16969968 and metabolite traits_NS_model1.csv")


# g) Use instrumental variable (IV) analysis to obtain estimates of the causal associations between tobacco exposure (cotinine levels) and metabolite levels using two-sample methods i.e. taking effect estimate of rs1051730/rs16969968  cotinine association in independent samples (from JNCI paper ) and effect estimates of rs1051730/rs16969968  metabolites associations. 
# rst=NULL
# require(sem)
# for(i in metabolites){
#   data$m = data[,i]
#   model = tsls(m ~ my.cigreg , ~rs16969968, data, subset = which(!is.na(data$rs16969968)))
#   coef = cbind(model$coef, confint.default(model))
#   pvalue = 2*pnorm(-abs(coef[,1]*1.96/(coef[,1]-coef[,2])))
#   coef = cbind(coef, pvalue)
#   rst = rbind(rst, coef[2,])
# }
# rownames(rst) = metabolites
# write.csv(rst, file = "IV analysis of rs16969968_model.csv")
# 
# rst=NULL
# require(sem)
# for(i in metabolites){
#   data$m = data[,i]
#   model = tsls(m ~ rs16969968 
#                +utalter+ucsex 
#                #+ utbmi + my.alkkon
#                , data
#                , subset = which(!is.na(data$rs16969968))
#                )
#   rst = rbind(rst, summary(model)$coef[2,])
# }
# rownames(rst) = metabolites

## Taylor series expansion
estimate.taylor = function(theta1, theta2, var1, var2, cov12=0){
  return (theta1/theta2+theta1*var2/theta2^3 -cov12/theta2^2)
}

ci.Fieller = function(theta1, theta2, var1, var2, cov12=0){
  a = (theta1/theta2)/(1-1.96^2*var2/theta2^2)
  b = 1.96*sqrt(var1/theta1^2+var2/theta2^2-1.96*(var1/theta1^2)*(var2/theta2^2))
  rst = data.frame(a*(1+b), a*(1-b))
  return(rst)
}

variation.taylor = function(theta1, theta2, var1, var2, cov12=0){
  return (theta1^2*var2/theta2^4 + var1/theta2^2 - 2*theta1*cov12/theta2^3)
}

snpmetabo = read.csv("IV analysis/associations of rs1051730 and metabolite traits_S_model2.csv",row.names = 1)
theta2 = 138.72; var2 = ((138.72-97.91)/1.96)^2*2932
rst = apply(snpmetabo, 1, 
      function(x){
        beta = estimate.taylor(theta1=x[1], theta2, var1 = x[2]^2*177, var2 = var2)
        beta.var = variation.taylor(theta1=x[1], theta2, var1 = x[2]^2*177, var2 = var2)
        beta.up = beta + 1.96*sqrt(beta.var)
        beta.lo = beta - 1.96*sqrt(beta.var)
        return(c(beta, beta.lo, beta.up))
      }
    )
rst = t(rst)
colnames(rst) = c("beta", "beta.lo", "beta.up")

# h) Produce graphs of change in metabolite levels per-allele increase of rs1051730/ rs16969968 against the corresponding change in metabolites per-unit increase in tobacco exposure (cotinine levels) respectively, where the “unit” corresponds to the per-allele effect of rs1051730/rs16969968 on cotinine levels. 


plotchange <- function(x,y, erx, ery){
  xmax = max(x+erx);xmin = min(x-erx)
  ymax = max(y+ery);ymin = min(y-ery)
  
  plot(x,y, xlim=c(xmin, xmax) , ylim = c(ymin, ymax), 
       ylab = "metabolite change per allel", xlab = "metabolite change per 138.72 unite cotinine") 
  arrows(x,y,x+erx,y, angle=90, length=.1) 
  arrows(x,y,x-erx,y, angle=90, length=.1) 
  arrows(x,y,x,y-ery, angle=90, length=.1) 
  arrows(x,y,x,y+ery, angle=90, length=.1) 
}

metabo.percontine = read.csv("IV analysis of rs1051730_model2.csv", row.names = 1)
metabo.perallele = read.csv("IV analysis/associations of rs1051730 and metabolite traits_S_model2.csv", row.names = 1)

pdf("metabolite change per-allele against per unit continine_model2.pdf", width = 6, height = 9)
plotchange( metabo.percontine[,1],metabo.perallele[,1], erx=metabo.percontine[,3]-metabo.percontine[,1],ery = metabo.perallele[,2]*1.96)
abline(h = 0, col = "red", lty = 2)
abline(v = 0, col = "red", lty = 2)
dev.off()