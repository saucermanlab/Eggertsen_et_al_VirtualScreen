library(tidyverse)
library(emmeans)

T_1 <- c(2657.19175627240, 3366.69921875000, 2132.78521617852, 3332.97484276730)
N_1 <- c(3004.41935483871, 2199.72285251216, 3453.30935251799, 2422.92448979592)
TM20_1 <- c(2260.71309192201, 2339.00819672131, 2813.93406593407, 3390.01716738197)
TM80_1 <- c(2191.17534246575, 1948.83681214421, 2058.32318501171, 2109.80837004405)
TM320_1 <- c(1863.81992337165, 1625.52090032154, 1359.59550561798, 1215.93734335840)
NM20_1 <- c(1898.57563587684, 2312.04672897196, 3459.46494464945, 2780.76063829787)
NM80_1 <- c(2464.00844594595, 2534.65816326531, 2994.96009975062, 2568.39232409382)
NM320_1 <- c(2311.32213438735, 2323.44482173175, 2695.08064516129, 2402.77701149425)
M_1 <- c(2296.53638814016, 1806.12622950820, 1990.10359408034, 2324.74598070740)
T_2 <- c(4153.15873015873, 2751.77209302326, 3429.43809523810, 4199.79381443299)
N_2 <- c(3877.84536082474, 3177.06611570248, 2748.67639902676, 3124.56022408964)
TM20_2 <- c(3205.05511811024, 3305.13452914798, 2929.63924050633, 3642.15167095116)
TM80_2 <- c(3007.02016129032, 2396.78801843318, 2774.29337539432, NA) #486.841269841270
TM320_2 <- c(1769.20895522388, 2089.20300751880, 1681.67757009346, 1645.79467680608)
NM20_2 <- c(3051.77973568282, 3665.76744186047, 3069.60496613995, 3520.99112426036)
NM80_2 <- c(3555.90551181102, 2599.49180327869, 3152.53145336226, NA) #337.050847457627
NM320_2 <- c(3498.53030303030, 2695.72164948454, 3643.93540051680, NA) #4691.18934911243
M_2 <- c(2419.74761904762, 2418.605678233439, 3305.6052631578950, NA) #214.985915492958

Traw <- data.frame(isolation=c(rep("1",20),rep("2",20)),
                    conditions=rep(c(rep("TGFb",4),
                                     rep("TGFb + 20nM Midostaurin",4),
                                     rep("TGFb + 80nM Midostaurin",4),
                                     rep("TGFb + 320nM Midostaurin",4),
                                     rep("320nM Midostaurin",4)),2),
                    cellarea=c(T_1,TM20_1,TM80_1,TM320_1,M_1,T_2,TM20_2,TM80_2,TM320_2,M_2))

Nraw <- data.frame(isolation=c(rep("1",20),rep("2",20)),
                   conditions=rep(c(rep("NE",4),
                                    rep("NE + 20nM Midostaurin",4),
                                    rep("NE + 80nM Midostaurin",4),
                                    rep("NE + 320nM Midostaurin",4),
                                    rep("320nM Midostaurin",4)),2),
                   cellarea=c(N_1,NM20_1,NM80_1,NM320_1,M_1,N_2,NM20_2,NM80_2,NM320_2,M_2))

# releveling conditions to reference T_1
Traw <- Traw %>% 
  mutate(conditions=factor(conditions,levels=c("TGFb",
                                               "TGFb + 20nM Midostaurin",
                                               "TGFb + 80nM Midostaurin",
                                               "TGFb + 320nM Midostaurin",
                                               "320nM Midostaurin")))

Nraw <- Nraw %>% 
  mutate(conditions=factor(conditions,levels=c("NE",
                                               "NE + 20nM Midostaurin",
                                               "NE + 80nM Midostaurin",
                                               "NE + 320nM Midostaurin",
                                               "320nM Midostaurin")))

modelT <- lm(cellarea~isolation+conditions,data=Traw)
summary(modelT)  #these are RAW t-tests with RAW p-values

modelN <- lm(cellarea~isolation+conditions,data=Nraw)
summary(modelN)  #these are RAW t-tests with RAW p-values

# post-hoc
Tmeans <- emmeans(modelT,specs= trt.vs.ctrl ~ conditions) # dunnet's, compare treatments to control
Tmeans$contrasts

Nmeans <- emmeans(modelN,specs= trt.vs.ctrl ~ conditions) # dunnet's, compare treatments to control
Nmeans$contrasts


### Cell Counts ###

cT_1 <- c(680,398,786,345)
cN_1 <- c(441,702,545,562)
cTM20_1 <- c(783,570,645,366)
cTM80_1 <- c(421,575,500,504)
cTM320_1 <- c(284,356,587,434)
cNM20_1 <- c(809,731,473,476)
cNM80_1 <- c(667,457,473,536)
cNM320_1 <- c(525,605,401,461)
cM_1 <- c(391,620,493,343)
cT_2 <- c(126,215,525,388)
cN_2 <- c(97,242,411,357)
cTM20_2 <- c(127,223,158,389)
cTM80_2 <- c(248,217,317,NA)
cTM320_2 <- c(201,266,214,263)
cNM20_2 <- c(227,258,443,338)
cNM80_2 <- c(127,183,461,NA)
cNM320_2 <- c(198,291,387,NA)
cM_2 <- c(210,317,190,NA)

Traw <- data.frame(isolation=c(rep("1",20),rep("2",20)),
                   conditions=rep(c(rep("TGFb",4),
                                    rep("TGFb + 20nM Midostaurin",4),
                                    rep("TGFb + 80nM Midostaurin",4),
                                    rep("TGFb + 320nM Midostaurin",4),
                                    rep("320nM Midostaurin",4)),2),
                   cellarea=c(cT_1,cTM20_1,cTM80_1,cTM320_1,cM_1,cT_2,cTM20_2,cTM80_2,cTM320_2,cM_2))

Nraw <- data.frame(isolation=c(rep("1",20),rep("2",20)),
                   conditions=rep(c(rep("NE",4),
                                    rep("NE + 20nM Midostaurin",4),
                                    rep("NE + 80nM Midostaurin",4),
                                    rep("NE + 320nM Midostaurin",4),
                                    rep("320nM Midostaurin",4)),2),
                   cellarea=c(cN_1,cNM20_1,cNM80_1,cNM320_1,cM_1,cN_2,cNM20_2,cNM80_2,cNM320_2,cM_2))

# releveling conditions to reference T_1
Traw <- Traw %>% 
  mutate(conditions=factor(conditions,levels=c("TGFb",
                                               "TGFb + 20nM Midostaurin",
                                               "TGFb + 80nM Midostaurin",
                                               "TGFb + 320nM Midostaurin",
                                               "320nM Midostaurin")))

Nraw <- Nraw %>% 
  mutate(conditions=factor(conditions,levels=c("NE",
                                               "NE + 20nM Midostaurin",
                                               "NE + 80nM Midostaurin",
                                               "NE + 320nM Midostaurin",
                                               "320nM Midostaurin")))

modelT <- lm(cellarea~isolation+conditions,data=Traw)
summary(modelT)  #these are RAW t-tests with RAW p-values

modelN <- lm(cellarea~isolation+conditions,data=Nraw)
summary(modelN)  #these are RAW t-tests with RAW p-values

# post-hoc
Tmeans <- emmeans(modelT,specs= trt.vs.ctrl ~ conditions) # dunnet's, compare treatments to control
Tmeans$contrasts

Nmeans <- emmeans(modelN,specs= trt.vs.ctrl ~ conditions) # dunnet's, compare treatments to control
Nmeans$contrasts
