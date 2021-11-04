#################################################################################
#
# Example code for Trial construct 2, effect size = 1.85, sample size = 240
#
#################################################################################

library(nlme)
library(lme4)
library(lmerTest)
library(tidyr)
library(lsmeans)
library(dplyr)

load("mcidata.Rdata")

slp1 <- c(1.85/18^2)

############ accrual pattern
ptno <- unique(mcidata$ptno)
ptcomp <- mciwide1$ptno[!is.na(mciwide1$m18)]
pt03 <- mciwide1$ptno[is.na(mciwide1$m03)] 
pt03 <- pt03[!pt03 %in% ptcomp]
pt06 <- mciwide1$ptno[is.na(mciwide1$m06)]
pt06 <- pt06[!pt06 %in% c(pt03, ptcomp)]
pt12 <- mciwide1$ptno[is.na(mciwide1$m12)]
pt12 <- pt12[!pt12 %in% c(pt03, pt06, ptcomp)]
pt18 <- mciwide1$ptno[is.na(mciwide1$m18)]
pt18 <- pt18[!pt18 %in% c(pt03, pt06, pt12, ptcomp)]

  
tpva <- lmepva <- glspva <- tpvva <- lmepvva <- glspvva <- tpvvva <- lmepvvva <- glspvvva <- 
  tpvvvva <- lmepvvvva <- glspvvvva <- NULL
for (i in 1:334) {
  set.seed(i)
  
  ## Stratified sampling from a real randomized trial
  simid <- c(sample(ptcomp, 210, replace = T), sample(pt03, 22, replace = T), sample(pt06, 13, replace = T),
             sample(pt12, 14, replace = T), sample(pt18, 21, replace = T))
  simid <- sample(simid)
  idx <- sample(c(1:280), 140, replace = F)
  txid <- simid[idx]
  ctid <- simid[-idx]
  
  simdat1 <- mciwide1[match(ctid, mciwide1$ptno),]
  simdat2 <- mciwide1[match(txid, mciwide1$ptno),]
  simdat1$id <- 1:140
  simdat2$id <- 141:280
  simdat1$arm2 <- "control"
  simdat2$arm2 <- "tx"
  
  simdat1 <- simdat1 %>% mutate(m03chg = m03 - sc,
                                m06chg = m06 - sc,
                                m12chg = m12 - sc,
                                m18chg = m18 - sc)
  
  simdat2 <- simdat2 %>% mutate(m03chg = m03 - sc - slp1*3^2,
                                m06chg = m06 - sc - slp1*6^2,
                                m12chg = m12 - sc - slp1*12^2,
                                m18chg = m18 - sc - slp1*18^2)
  
  simdat <- rbind(simdat1, simdat2)
  
  
  newid <- simdat$id[!is.na(simdat$m18chg)]
  newid2 <- sample(newid, 34)
  
  idm12 <- simdat$id[!is.na(simdat$m12chg)]
  idm12 <- idm12[!idm12 %in% newid2]
  idm12 <- sample(idm12, 36)
  
  idm06 <- simdat$id[!is.na(simdat$m06chg)]
  idm06 <- idm06[!idm06 %in% c(newid2, idm12)]
  idm06 <- sample(idm06, 70)
  
  idm03 <- simdat$id[!is.na(simdat$m03chg)]
  idm03 <- idm03[!idm03 %in% c(newid2, idm12, idm06)]
  idm03 <- sample(idm03, 84) 
  
  idrest <- simdat$id[!simdat$id %in% c(newid2, idm12, idm06, idm03)] 
  
  newid3 <- simdat$id[!simdat$id %in% newid2]
  newid4 <- sample(newid3, 180)
  newid5 <- c(newid2, newid4)
  
  ######################################################################
  #
  # Scenario 0 - follow original protocol
  #
  ######################################################################
  
  simdatS0 <- simdat
  
  mcilong <- simdatS0 %>% pivot_longer(m03chg:m18chg, names_to = "viscode", values_to = "adas") %>% mutate(
    viscode = factor(viscode, levels = c("m03chg","m06chg","m12chg","m18chg")),
    arm2 = factor(arm2, levels = c("control", "tx")),
    vis = case_when(
      viscode == 'm03chg' ~ 3,
      viscode == 'm06chg' ~ 6,
      viscode == 'm12chg' ~ 12,
      viscode == 'm18chg' ~ 18
    )
  ) %>% arrange(id, vis)
  
  mcinomiss <- mcilong[!is.na(mcilong$adas),]
  
  ### Slope model
  mod1 <- lme(adas ~ sc + sc*vis +
                arm2 + vis + vis * arm2, random = ~ 1 | id, data = mcinomiss)
  
  lsm1 <-  lsmeans(mod1, c("vis", "arm2"), at = list(vis=18) )  
  con1 <-  lsm1 %>% lsmeans::contrast("trt.vs.ctrl1", adjust= "none")  %>% summary  %>% as.data.frame
  lmepva <- c(lmepva, con1[1,6])
  
  ### Categorical model
  mod2 <- gls(adas ~ sc + sc*viscode + arm2 * viscode,
              weights = varIdent(form ~ 1 | viscode),
              data = mcinomiss, 
              correlation = corSymm(form = ~ 1 | id))
  lsm2 <-  lsmeans(mod2, c("viscode", "arm2"), by = "viscode" )  
  con2 <-  lsm2 %>% lsmeans::contrast("trt.vs.ctrl1", adjust= "none")  %>% summary  %>% as.data.frame
  glspva <- c(glspva, con2[4,7])
  
  ### T test
  ttest <- t.test(simdat$m18chg[simdat$arm2 == "control"], simdat$m18chg[simdat$arm2 == "tx"])
  tpva <- c(tpva, ttest$p.value)
  
  
  ######################################################################
  #
  # Scenario 1 - Study end when COVID hit
  #
  ######################################################################
  
  simdatS1 <- simdat
  
  newdata <- simdatS1 %>% filter(id %in% newid5) %>% mutate(
    m18chg = case_when(
      !id %in% newid2 ~ -9999,
      TRUE ~ m18chg
    ),
    m12chg = case_when(
      !id %in% c(idm12,newid2) ~ -9999,
      TRUE ~ m12chg
    ),
    m06chg = case_when(
      !id %in% c(idm06, idm12, newid2) ~ -9999,
      TRUE ~ m06chg
    ),
    m03chg = case_when(
      !id %in% c(idm03, idm06, idm12, newid2) ~ -9999,
      TRUE ~ m03chg
    )
  )
  newdata[newdata == -9999] <- NA
  
  mcilong <- newdata %>% pivot_longer(m03chg:m18chg, names_to = "viscode", values_to = "adas") %>% mutate(
    viscode = factor(viscode, levels = c("m03chg","m06chg","m12chg","m18chg")),
    arm2 = factor(arm2, levels = c("control", "tx")),
    vis = case_when(
      viscode == 'm03chg' ~ 3,
      viscode == 'm06chg' ~ 6,
      viscode == 'm12chg' ~ 12,
      viscode == 'm18chg' ~ 18
    )
  ) %>% arrange(id, vis)
  
  mcinomiss <- mcilong[!is.na(mcilong$adas),]
  
  ### Slope model
  mod1 <- lme(adas ~ sc + sc * arm2 +
                arm2 + vis + vis * arm2, random = ~ 1 | id, data = mcinomiss)
  
  lsm1 <-  lsmeans(mod1, c("vis", "arm2"), at = list(vis=18) )  
  con1 <-  lsm1 %>% lsmeans::contrast("trt.vs.ctrl1", adjust= "none")  %>% summary  %>% as.data.frame
  lmepvva <- c(lmepvva, con1[1,6])
  
  ### Categorical model
  mod2 <- gls(adas ~ sc + sc * viscode + arm2 + viscode + arm2 * viscode,
              varIdent(form ~ 1 | viscode),
              data = mcinomiss,
              correlation = corSymm(form = ~ 1 | id))
  
  lsm2 <-  lsmeans(mod2, c("viscode", "arm2"), by = "viscode", mode = "df.error" )
  
  tryCatch({
    lsm2 <-  lsmeans(mod2, c("viscode", "arm2"), by = "viscode" )
  },error=function(e){"error"})
  
  con2 <-  lsm2 %>% lsmeans::contrast("trt.vs.ctrl1", adjust= "none")  %>% summary  %>% as.data.frame
  glspvva <- c(glspvva, con2[4,7])
  
  ### T test
  ttest <- t.test(newdata$m18chg[newdata$arm2 == "control"], newdata$m18chg[newdata$arm2 == "tx"])
  tpvva <- c(tpvva, ttest$p.value)
  
  
  
  ######################################################################
  #
  # Scenario 2 - Study resume after COVID (6 month period)
  #
  ######################################################################
  
  ## When Covid hit, treatment group will miss 6 months tx effects
  simdat2 <- simdat2 %>% mutate(m03chg = case_when(
    !id %in% c(idm03, idm06, idm12, newid2) ~ m03chg + slp1*3^2,
    TRUE ~ m03chg
  ), m18chg = case_when(
    !id %in% newid2 ~ m18chg + slp1*6^2,
    TRUE ~ m18chg
  ), m12chg = case_when(
    !id %in% c(idm12, newid2) ~ m12chg + slp1*6^2,
    TRUE ~ m12chg
  ), m06chg = case_when(
    !id %in% c(idm03, idm06, idm12, newid2) ~ m06chg + slp1*6^2,
    id %in% idm03 ~ m06chg + slp1*3^2,
    TRUE ~ m06chg
  )
  )
  
  simdatS2 <- rbind(simdat1, simdat2)
  
  newdata <- simdatS2 %>% filter(id %in% newid5) %>% mutate(
    m18chg = case_when(
      id %in% idm12 ~ -9999,
      TRUE ~ m18chg
    ),
    m12chg = case_when(
      id %in% idm06 ~ -9999,
      TRUE ~ m12chg
    ),
    m06chg = case_when(
      id %in% c(idm03, idrest) ~ -9999,
      TRUE ~ m06chg
    ),
    m03chg = case_when(
      id %in% idrest ~ -9999,
      TRUE ~ m03chg
    )
  )
  newdata[newdata == -9999] <- NA
  
  mcilong <- newdata %>% pivot_longer(m03chg:m18chg, names_to = "viscode", values_to = "adas") %>% mutate(
    viscode = factor(viscode, levels = c("m03chg","m06chg","m12chg","m18chg")),
    arm2 = factor(arm2, levels = c("control", "tx")),
    vis = case_when(
      viscode == 'm03chg' ~ 3,
      viscode == 'm06chg' ~ 6,
      viscode == 'm12chg' ~ 12,
      viscode == 'm18chg' ~ 18
    )
  ) %>% arrange(id, vis)
  
  mcinomiss <- mcilong[!is.na(mcilong$adas),]
  
  ### Slope model
  mod1 <- lme(adas ~ sc + sc * arm2 +
                arm2 + vis + vis * arm2, random = ~ 1 | id, data = mcinomiss)
  lsm1 <-  lsmeans(mod1, c("vis", "arm2"), at = list(vis=18) )  
  con1 <-  lsm1 %>% lsmeans::contrast("trt.vs.ctrl1", adjust= "none")  %>% summary  %>% as.data.frame
  lmepvvva <- c(lmepvvva, con1[1,6])
  
  ### Categorical model
  mod2 <- gls(adas ~ sc + sc * viscode + arm2 + viscode + arm2 * viscode,
              varIdent(form ~ 1 | viscode),
              data = mcinomiss,
              correlation = corSymm(form = ~ 1 | id))
  lsm2 <-  lsmeans(mod2, c("viscode", "arm2"), by = "viscode", mode = "df.error" )
  
  tryCatch({
    lsm2 <-  lsmeans(mod2, c("viscode", "arm2"), by = "viscode" )
  },error=function(e){"error"})
  
  con2 <-  lsm2 %>% lsmeans::contrast("trt.vs.ctrl1", adjust= "none")  %>% summary  %>% as.data.frame
  #sandpv2a <- c(sandpv2a, coef_test(mod2, vcov = "CR0", coefs = 12)[,5])
  glspvvva <- c(glspvvva, con2[4,7])
  
  ### T test
  ttest <- t.test(newdata$m18chg[newdata$arm2 == "control"], newdata$m18chg[newdata$arm2 == "tx"])
  tpvvva <- c(tpvvva, ttest$p.value)
  
  
  ######################################################################
  #
  # Scenario 3 - Remote assessments during COVID (6 month period)
  #
  ######################################################################
  
  simdatS3 <- rbind(simdat1, simdat2)
  
  ## By remotely assessments, ADAS should shift by 0.5 points
  newdata2 <- simdatS3 %>% filter(id %in% newid5) %>% mutate(
    m18chg = case_when(
      id %in% idm12 ~ m18chg+.5,
      TRUE ~ m18chg
    ),
    m12chg = case_when(
      id %in% idm06 ~ m12chg+.5,
      TRUE ~ m12chg
    ),
    m06chg = case_when(
      id %in% c(idm03, idrest) ~ m06chg+.5,
      TRUE ~ m06chg
    ),
    m03chg = case_when(
      id %in% idrest ~ m03chg+.5,
      TRUE ~ m03chg
    )
  )
  
  
  mcilong <- newdata %>% pivot_longer(m03chg:m18chg, names_to = "viscode", values_to = "adas") %>% mutate(
    viscode = factor(viscode, levels = c("m03chg","m06chg","m12chg","m18chg")),
    arm2 = factor(arm2, levels = c("control", "tx")),
    vis = case_when(
      viscode == 'm03chg' ~ 3,
      viscode == 'm06chg' ~ 6,
      viscode == 'm12chg' ~ 12,
      viscode == 'm18chg' ~ 18
    )
  ) %>% arrange(id, vis)
  
  mcinomiss <- mcilong[!is.na(mcilong$adas),]
  
  ### Slope model
  mod1 <- lme(adas ~ sc + sc * arm2 +
                arm2 + vis + vis * arm2, random = ~ 1 | id, data = mcinomiss)
  lsm1 <-  lsmeans(mod1, c("vis", "arm2"), at = list(vis=18) )  
  con1 <-  lsm1 %>% lsmeans::contrast("trt.vs.ctrl1", adjust= "none")  %>% summary  %>% as.data.frame
  lmepvvvva <- c(lmepvvvva, con1[1,6])
  
  ### Categorical model
  mod2 <- gls(adas ~ sc + sc * viscode + arm2 + viscode + arm2 * viscode,
              varIdent(form ~ 1 | viscode),
              data = mcinomiss,
              correlation = corSymm(form = ~ 1 | id))
  lsm2 <-  lsmeans(mod2, c("viscode", "arm2"), by = "viscode", mode = "df.error" )
  
  tryCatch({
    lsm2 <-  lsmeans(mod2, c("viscode", "arm2"), by = "viscode" )
  },error=function(e){"error"})
  
  con2 <-  lsm2 %>% lsmeans::contrast("trt.vs.ctrl1", adjust= "none")  %>% summary  %>% as.data.frame
  #sandpv2a <- c(sandpv2a, coef_test(mod2, vcov = "CR0", coefs = 12)[,5])
  glspvvvva <- c(glspvvvva, con2[4,7])
  
  ### T test
  ttest <- t.test(newdata$m18chg[newdata$arm2 == "control"], newdata$m18chg[newdata$arm2 == "tx"])
  tpvvvva <- c(tpvvvva, ttest$p.value)
  
  
  print(i)
  
}

adas1d1.85 <- data.frame(lmepvva, glspvva, tpvva)
