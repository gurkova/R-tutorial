rm(list=ls())

library("memisc")
library("dplyr")
library("psych")
library("lmtest")
library("sjPlot")
library("sgof")
library("stargazer") 
library("sandwich")
library("ggplot2")
library("foreign")
library("car")
library("haven")
library("plm")
library("hexbin")
library("devtools")
library("rlms")
library("jtools")
library("aod")
library("labelled")

#Function cse is useful to apply robust standard errors in case of heteroscedasticity
#Applied with OLS
cse = function(reg) {
  rob = sqrt(diag(vcovHC(reg, type = "HC1")))
  return(rob)
}


#clse function for correct standard errors
#Applied in case of panel data (with FE method)
#clustered SEs, clustered on "group"... could also cluster on "time" 
clse = function(reg) { 
  # index(reg, "id") returns the id or entity variable vector 
  G = length(unique(index(reg,"id")))
  N = length(index(reg,"id"))
  dfa = (G/(G - 1))   # note Bluhm multiplies this by finite-sample df adjustment
  rob = sqrt(diag(dfa*vcovHC(reg, method="arellano", type = "HC1", 
                             cluster = "group")))
  return(rob)
}

#ivse function is applied for robust standard errors with TSLS method
ivse = function(reg) {
  rob = robust.se(reg)[,2]
  return(rob)
}

memory.limit()
memory.limit(size=100000)


D <- rlms_read("USER_RLMS.sav")
saveRDS(D, "Data_rlms.rds")

#Chossing required variables
H <- select(D, id_w, idind, year, psu, status, marst, occup08, h5, h6, age, i2, i3, i8, j8, j10, j21.1.3, j21.1.4, 
            j21.1.5, j21.1.6, j21.1.7, j21.1.8, j21.1.9, j23, j24, j25, j26, j72.18a, j72.171, j72.172, j72.173, 
            j216ac, j216bc, j217a, j217b, l5, m20.7, m20.8)

saveRDS(H, "Data_dirty.rds")

#Rename variables
colnames(H)[colnames(H)=="h6"] <- "birth_year"
colnames(H)[colnames(H)=="i3"] <- "birth_place"
colnames(H)[colnames(H)=="i2"] <- "birth_country"
colnames(H)[colnames(H)=="i8"] <- "rus_year"
colnames(H)[colnames(H)=="j8"] <- "hours_month"
colnames(H)[colnames(H)=="j10"] <- "wage_month"
colnames(H)[colnames(H)=="m20.7"] <- "invalid"
colnames(H)[colnames(H)=="m20.8"] <- "inv_group"
colnames(H)[colnames(H)=="occup08"] <- "worktype"
colnames(H)[colnames(H)=="idind"] <- "id"
colnames(H)[colnames(H)=="marst"] <- "married"
colnames(H)[colnames(H)=="year"] <- "year_w"
colnames(H)[colnames(H)=="j72.18a"] <- "educ"
colnames(H)[colnames(H)=="j72.171"] <- "kids"
colnames(H)[colnames(H)=="j72.172"] <- "n_kids"
colnames(H)[colnames(H)=="j72.173"] <- "kids18"
colnames(H)[colnames(H)=="j216ac"] <- "father_work"
colnames(H)[colnames(H)=="j216bc"] <- "mother_work"
colnames(H)[colnames(H)=="j217a"] <- "father_educ"
colnames(H)[colnames(H)=="j217b"] <- "mother_educ"
colnames(H)[colnames(H)=="l5"] <- "health_problems"
colnames(H)[colnames(H)=="h5"] <- "female"
colnames(H)[colnames(H)=="psu"] <- "region"

#Recode variables
H$female[D$h5 == 2] <- 1
H$female[D$h5 == 1] <- 0

H$benefits <- NA
H$benefits[H$j21.1.3 == 1 | H$j21.1.4 == 1 | H$j21.1.5 == 1 | H$j21.1.6 == 1 | H$j21.1.7 | H$j21.1.8 == 1 | H$j21.1.9 ==1] <- 1
H$benefits[H$j21.1.3 == 2 | H$j21.1.4 == 2 | H$j21.1.5 == 2 | H$j21.1.6 == 2 | H$j21.1.7 | H$j21.1.8 == 2 | H$j21.1.9 ==2] <- 0

H$invalid[H$invalid == 1] <- 1
H$invalid[H$invalid == 2 | H$invalid == 5] <- 0

H$inv_group <- as.factor(H$inv_group)

H$worktype <- as.factor(H$worktype)

H$married[D$marst == 2 | D$marst == 3 | D$marst == 7] <- 1
H$married[D$marst == 1 | D$marst == 4 | D$marst == 5 | D$marst == 6] <- 0

H$gov <- NA
H$gov[H$j23 == 1 ] <- 1
H$gov[H$j23 == 2] <- 0

H$kids[H$kids == 2] <- 0

H$school <- NA
H$school[H$educ == 1 | H$educ == 2] <- 1
H$school[H$educ > 2 & H$educ < 18] <- 0

H$ptu <- NA
H$ptu[H$educ > 2 & H$educ < 7] <- 1
H$ptu[H$educ >= 7 & H$educ < 18] <- 0
H$ptu[H$educ == 15] <- 1
H$ptu[H$educ < 3] <- 0

H$university <- NA
H$university[H$educ == 7 | H$educ == 10 | H$educ == 11 | H$educ == 12 | H$educ == 16] <- 1
H$university[H$educ < 7 | H$educ == 8 | H$educ == 9 | H$educ == 13 | H$educ == 14 | H$educ == 15 
             | H$educ == 17] <- 0

H$ucheny <- NA
H$ucheny[H$educ == 8 | H$educ == 9 | H$educ == 8 | H$educ == 13 | H$educ == 14 | H$educ == 17] <- 1
H$ucheny[H$educ < 8 | H$educ == 10 | H$educ == 11 | H$educ == 12 | H$educ == 15 | H$educ == 16] <- 0

H$educ <- as.factor(H$educ)

H$health_problems[H$health_problems == 2] <- 0

H$father_work[H$father_work == 0110] <- 0
H$father_work[H$father_work >= 1110 & H$father_work < 2111] <- 1
H$father_work[H$father_work >= 2111 & H$father_work < 3111] <- 2
H$father_work[H$father_work >= 3111 & H$father_work < 4111] <- 3
H$father_work[H$father_work >= 4111 & H$father_work < 5111] <- 4
H$father_work[H$father_work >= 5111 & H$father_work < 6111] <- 5
H$father_work[H$father_work >= 6111 & H$father_work < 7111] <- 6
H$father_work[H$father_work >= 7111 & H$father_work < 8111] <- 7
H$father_work[H$father_work >= 8111 & H$father_work < 9111] <- 8
H$father_work[H$father_work >= 9111] <- 9

H$father_work <- as.factor(H$father_work)

H$mother_work[H$mother_work == 0110] <- 0
H$mother_work[H$mother_work >= 1110 & H$mother_work < 2111] <- 1
H$mother_work[H$mother_work >= 2111 & H$mother_work < 3111] <- 2
H$mother_work[H$mother_work >= 3111 & H$mother_work < 4111] <- 3
H$mother_work[H$mother_work >= 4111 & H$mother_work < 5111] <- 4
H$mother_work[H$mother_work >= 5111 & H$mother_work < 6111] <- 5
H$mother_work[H$mother_work >= 6111 & H$mother_work < 7111] <- 6
H$mother_work[H$mother_work >= 7111 & H$mother_work < 8111] <- 7
H$mother_work[H$mother_work >= 8111 & H$mother_work < 9111] <- 8
H$mother_work[H$mother_work >= 9111] <- 9

H$mother_work <- as.factor(H$mother_work)

H$f_noeduc <- NA
H$f_noeduc[H$father_educ == 1] <- 1
H$f_noeduc[H$father_educ > 1 & H$father_educ < 13] <- 0
summary(H$f_noeduc)

H$f_school <- NA
H$f_school[H$father_educ == 2 | H$father_educ == 3 | H$father_educ == 7 | H$father_educ == 12] <- 1
H$f_school[H$father_educ == 1 | H$father_educ == 4 | H$father_educ == 5 | H$father_educ == 6 
           | H$father_educ == 8 | H$father_educ == 9 | H$father_educ == 10 | H$father_educ == 11] <- 0

H$f_ptu <- NA
H$f_ptu[H$father_educ == 4 | H$father_educ == 5 | H$father_educ == 6 | H$father_educ == 8] <- 1
H$f_ptu[H$father_educ == 1 | H$father_educ == 2 | H$father_educ == 3 | H$father_educ == 7 
        | H$father_educ == 9 | H$father_educ == 10 | H$father_educ == 11 | H$father_educ == 12] <- 0

H$f_university <- NA
H$f_university[H$father_educ == 9 | H$father_educ == 10 | H$father_educ == 11] <- 1
H$f_university[H$father_educ < 9 | H$father_educ == 12] <- 0


H$m_noeduc <- NA
H$m_noeduc[H$mother_educ == 1] <- 1
H$m_noeduc[H$mother_educ > 1 & H$mother_educ < 13] <- 0


H$m_school <- NA
H$m_school[H$mother_educ == 2 | H$mother_educ == 3 | H$mother_educ == 7 | H$mother_educ == 12] <- 1
H$m_school[H$mother_educ == 1 | H$mother_educ == 4 | H$mother_educ == 5 | H$mother_educ == 6 
           | H$mother_educ == 8 | H$mother_educ == 9 | H$mother_educ == 10 | H$mother_educ == 11] <- 0

H$m_ptu <- NA
H$m_ptu[H$mother_educ == 4 | H$mother_educ == 5 | H$mother_educ == 6 | H$mother_educ == 8] <- 1
H$m_ptu[H$mother_educ == 1 | H$mother_educ == 2 | H$mother_educ == 3 | H$mother_educ == 7 
        | H$mother_educ == 9 | H$mother_educ == 10 | H$mother_educ == 11 | H$mother_educ == 12] <- 0

H$m_university <- NA
H$m_university[H$mother_educ == 9 | H$mother_educ == 10 | H$mother_educ == 11] <- 1
H$m_university[H$mother_educ < 9 | H$mother_educ == 12] <- 0


H$age2 <- H$age*H$age

H$status <- as.factor(H$status)

H$birth_country <- as.factor(H$birth_country)
H$birth_place <- as.factor(H$birth_place)

#Recode wages to make log of hourly wages
H$realwage <- H$wage_month/H$hours_month


#Estimate real wages with data on inflation from IMF
H$inflation <- 1

H$inflation[H$id_w==5] <- 17.97
H$inflation[H$id_w==6] <- 15.66
H$inflation[H$id_w==7] <- 12.23
H$inflation[H$id_w==8] <- 6.60
H$inflation[H$id_w==9] <- 5.47
H$inflation[H$id_w==10] <- 4.50
H$inflation[H$id_w==11] <- 3.89
H$inflation[H$id_w==12] <- 3.42
H$inflation[H$id_w==13] <- 3.08
H$inflation[H$id_w==14] <- 2.74
H$inflation[H$id_w==15] <- 2.50
H$inflation[H$id_w==16] <- 2.29
H$inflation[H$id_w==17] <- 2.01
H$inflation[H$id_w==18] <- 1.8
H$inflation[H$id_w==19] <- 1.68
H$inflation[H$id_w==20] <- 1.55
H$inflation[H$id_w==21] <- 1.48
H$inflation[H$id_w==22] <- 1.38
H$inflation[H$id_w==23] <- 1.28
H$inflation[H$id_w==24] <- 1.11
H$inflation[H$id_w==25] <- 1.04

H$realwage < H$realwage*H$inflation
H$log.real.wage <- log(H$realwage)

#Subset dataframe 
K <- subset(H, age %in% c(14:72))
K <- subset(K, wage_month %in% c(8:350000))


#Selest new dataframe (with NAs)
U <- select(K, id, id_w, year_w, female, region, status, married, worktype, age, age2, birth_country, birth_place, 
            log.real.wage, educ, school, ptu, university, ucheny, kids, n_kids, kids18, father_work, mother_work, 
            father_educ, mother_educ, f_noeduc, f_school, f_ptu, f_university, m_noeduc, m_school, m_ptu, m_university, health_problems, 
            invalid, inv_group, benefits, gov)

#Remove labels (discriptiona of variables)
U <- remove_labels(U)
summary(U)

saveRDS(U, "Data_inter.rds")


#Subset required dataset
N <- subset(U, worktype %in% c(0:9))
N <- subset(N, married %in% c(0:1))
N <- subset(N, educ %in% c(1:17))
N <- subset(N, invalid %in% c(0:1))

#Subset another datasets to be able to compare
K <- subset(N, father_work %in% c(0:9))
K <- subset(K, mother_work %in% c(0:9))
K <- subset(K, father_educ %in% c(1:12))
K <- subset(K, mother_educ %in% c(1:12))

L <- subset(K, health_problems %in% c(1:10000))

d <- density(L$wage_month) # returns the density data 
plot(d)

#Model influence of different factors on real wages of individual
model1 <- plm(log.real.wage ~ female + married + worktype + age + age2 + educ + kids + invalid, 
              data = S, index = c("id","id_w"), model="pooling")
summary(model1)


model2 <- plm(log.real.wage ~ female + married + worktype + age + age2 + ptu + university + ucheny 
              + kids + kids18 + n_kids + invalid + f_school + f_ptu + f_university + m_school + m_ptu 
              + m_university, data = K, index = c("id","id_w"), model="pooling")
summary(model2)

model3 <- plm(log.real.wage ~ female + married + worktype + age + age2 + ptu + university + ucheny 
              + kids + kids18 + n_kids + invalid + f_ptu + f_university + m_ptu 
              + m_university, data = K, index = c("id","id_w"), model="pooling")
summary(model3)

model4 <- plm(log.real.wage ~ female + married + worktype + age + age2 + ptu + university + ucheny 
              + kids + kids18 + n_kids + invalid + f_ptu + f_university + m_ptu 
              + m_university + father_work + mother_work, data = K, index = c("id","id_w"), model="pooling")
summary(model4)

model5 <- plm(log.real.wage ~ female + married + worktype + age + age2 + ptu + university + ucheny 
              + kids + kids18 + n_kids + invalid + f_ptu + f_university + m_ptu 
              + m_university + status, data = K, index = c("id","id_w"), model="pooling")
summary(model5)


model6 <- plm(log.real.wage ~ female + married + worktype + age + age2 + ptu + university + ucheny 
              + kids + kids18 + n_kids + inv_group + f_ptu + f_university + m_ptu 
              + m_university + status, data = K, index = c("id","id_w"), model="pooling")

summary(model6)

#Fixed effects model
model7 <- plm(log.real.wage ~ worktype + age + age2 + ptu + university + ucheny 
              + kids + kids18 + n_kids + invalid + f_ptu + f_university + m_ptu 
              + m_university, data = K, , effect="individual", model = "within", index = c("id","id_w"))
summary(model7)

model8 <- plm(log.real.wage ~ age, data = K, effect="individual", model = "within", index = c("id","id_w"))
summary(model7)
