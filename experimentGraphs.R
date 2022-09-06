library(ggplot2)
library(ggsignif)
library(stringr)

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
TM80_2 <- c(3007.02016129032, 2396.78801843318, 2774.29337539432) #486.841269841270
TM320_2 <- c(1769.20895522388, 2089.20300751880, 1681.67757009346, 1645.79467680608)
NM20_2 <- c(3051.77973568282, 3665.76744186047, 3069.60496613995, 3520.99112426036)
NM80_2 <- c(3555.90551181102, 2599.49180327869, 3152.53145336226) #337.050847457627
NM320_2 <- c(3498.53030303030, 2695.72164948454, 3643.93540051680) #4691.18934911243
M_2 <- c(2419.74761904762, 2418.605678233439, 3305.6052631578950) #214.985915492958

cellcount_T1 <- c(552.250000000000,591,500,415.250000000000,461.750000000000)
cellstd_T1 <- c(214.246548007975,174.017240525185,62.9338276816742,129.854726521602,122.646850754514)
cellcount_T2 <- c(313.500000000000,224.250000000000,260.666666666667,236,239)
cellstd_T2 <- c(44.5208284588985,29.2228328252641,17.0630638558305,8.33916462642792,22.7620541935721)
cellcount_N1 <- c(562.500000000000,622.250000000000,533.250000000000,498,461.750000000000)
cellstd_N1 <- c(107.283114545891,173.557627317269,95.4650896052234,87.4757109145162,122.646850754514)
cellcount_N2 <- c(276.750000000000,316.500000000000,257,292,239)
cellstd_N2 <- c(34.7553203361250,24.1078306780183,59.6247525035627,31.5013227235514,22.7620541935721)

se <- function(x) sqrt(var(x)/length(x))

theme_Bryan <- function () {
  theme_bw(base_size=14*96/72) %+replace%
    theme(
      #text = element_text(family = "sans"),
      plot.title = element_text(size = 7*96/72, hjust = 0.5), #hjust = 0.5,
      strip.text = element_text(size = 7*96/72, hjust = 0.5), #hjust = 0.5,
      strip.text.x = element_text(size = 7*96/72, hjust = 0.5), #hjust = 0.5,
      strip.text.y = element_text(size = 7*96/72, hjust = 0.5), #hjust = 0.5, #10 == 2.57mm
      strip.background = element_blank(),
      axis.title =element_text(size=6*96/72), #, face = "bold"
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.text = element_text(size = 5*96/72),
      legend.title=element_text(size=6*96/72),
      legend.text=element_text(size=5*96/72),
      axis.ticks.length =  unit(2.75, "pt"))
}

dfraw <- data.frame(isolation=c(rep("1",20),rep("2",20)),
                    conditions=rep(c("TGFb","TGFb + 20nM Midostaurin","TGFb + 80nM Midostaurin","TGFb + 320nM Midostaurin","320nM Midostaurin"),10))
pT <- c('0.4782','0.0009','<.0001','0.0010')
pN <- c('0.9961','0.9160','0.8321','0.0366')

# TGFb data
df <- data.frame(isolation=c(rep("1",5),rep("2",5)),
                 conditions=c("TGFb","TGFb + 20nM Midostaurin","TGFb + 80nM Midostaurin","TGFb + 320nM Midostaurin","320nM Midostaurin",
                         "TGFb","TGFb + 20nM Midostaurin","TGFb + 80nM Midostaurin","TGFb + 320nM Midostaurin","320nM Midostaurin"),
                 means=c(mean(T_1), mean(TM20_1), mean(TM80_1), mean(TM320_1), mean(M_1), 
                         mean(T_2), mean(TM20_2), mean(TM80_2), mean(TM320_2), mean(M_2)),
                 sm=c(se(T_1), se(TM20_1), se(TM80_1), se(TM320_1), se(M_1), 
                      se(T_2), se(TM20_2), se(TM80_2), se(TM320_2), se(M_2)))

ggplot(df, aes(x=conditions, y=means, fill=isolation)) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin=means-sm, ymax=means+sm), width=0.2, position=position_dodge(0.9)) +
  geom_text(x="TGFb + 80nM Midostaurin",y=3200,label="**", size = 8) +
  geom_text(x="TGFb + 320nM Midostaurin",y=2300,label="***", size = 8) +
  geom_text(x="320nM Midostaurin",y=3200,label="**", size = 8) +
  scale_x_discrete(limits = c("TGFb","TGFb + 20nM Midostaurin","TGFb + 80nM Midostaurin","TGFb + 320nM Midostaurin","320nM Midostaurin"),
                   labels = function(x) str_wrap(x,width = 10)) +
  scale_fill_grey() + 
  theme_Bryan()

# NE data
df <- data.frame(isolation=c(rep("1",5),rep("2",5)),
                 conditions=c("NE","NE + 20nM Midostaurin","NE + 80nM Midostaurin","NE + 320nM Midostaurin","320nM Midostaurin",
                              "NE","NE + 20nM Midostaurin","NE + 80nM Midostaurin","NE + 320nM Midostaurin","320nM Midostaurin"),
                 means=c(mean(N_1), mean(NM20_1), mean(NM80_1), mean(NM320_1), mean(M_1), 
                         mean(N_2), mean(NM20_2), mean(NM80_2), mean(NM320_2), mean(M_2)),
                 sm=c(se(N_1), se(NM20_1), se(NM80_1), se(NM320_1), se(M_1), 
                      se(N_2), se(NM20_2), se(NM80_2), se(NM320_2), se(M_2)))

ggplot(df, aes(x=conditions, y=means, fill=isolation)) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin=means-sm, ymax=means+sm), width=0.2, position=position_dodge(0.9)) +
  geom_text(x="320nM Midostaurin",y=3200,label="*", size = 8) +
  scale_x_discrete(limits = c("NE","NE + 20nM Midostaurin","NE + 80nM Midostaurin","NE + 320nM Midostaurin","320nM Midostaurin"),
                   labels = function(x) str_wrap(x,width = 10)) +
  scale_fill_grey() + 
  theme_Bryan()

# TGFb data cell count
df <- data.frame(isolation=c(rep("1",5),rep("2",5)),
                 conditions=c("TGFb","TGFb + 20nM Midostaurin","TGFb + 80nM Midostaurin","TGFb + 320nM Midostaurin","320nM Midostaurin",
                              "TGFb","TGFb + 20nM Midostaurin","TGFb + 80nM Midostaurin","TGFb + 320nM Midostaurin","320nM Midostaurin"),
                 means=c(cellcount_T1,cellcount_T2),
                 sm=c(cellstd_T1,cellstd_T2))

ggplot(df, aes(x=conditions, y=means, fill=isolation)) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin=means-sm, ymax=means+sm), width=0.2, position=position_dodge(0.9)) +
  geom_text(x="TGFb + 80nM Midostaurin",y=3200,label="**", size = 8) +
  geom_text(x="TGFb + 320nM Midostaurin",y=2300,label="***", size = 8) +
  geom_text(x="320nM Midostaurin",y=3200,label="**", size = 8) +
  scale_x_discrete(limits = c("TGFb","TGFb + 20nM Midostaurin","TGFb + 80nM Midostaurin","TGFb + 320nM Midostaurin","320nM Midostaurin"),
                   labels = function(x) str_wrap(x,width = 10)) +
  scale_fill_grey() + 
  theme_Bryan()

# NE data cell count
df <- data.frame(isolation=c(rep("1",5),rep("2",5)),
                 conditions=c("NE","NE + 20nM Midostaurin","NE + 80nM Midostaurin","NE + 320nM Midostaurin","320nM Midostaurin",
                              "NE","NE + 20nM Midostaurin","NE + 80nM Midostaurin","NE + 320nM Midostaurin","320nM Midostaurin"),
                 means=c(cellcount_N1,cellcount_N2),
                 sm=c(cellstd_N1,cellstd_N2))

ggplot(df, aes(x=conditions, y=means, fill=isolation)) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin=means-sm, ymax=means+sm), width=0.2, position=position_dodge(0.9)) +
  geom_text(x="320nM Midostaurin",y=3200,label="*", size = 8) +
  scale_x_discrete(limits = c("NE","NE + 20nM Midostaurin","NE + 80nM Midostaurin","NE + 320nM Midostaurin","320nM Midostaurin"),
                   labels = function(x) str_wrap(x,width = 10)) +
  scale_fill_grey() + 
  theme_Bryan()

## Predictions

NE_p <- c(0.943165541889592,0.941051500986458,0.939126263117417,0.935983825071283,0.159300000000000)
TGFb_p <- c(0.753618077319266,0.645558168173154,0.525232134717231,0.302167803139185,0.159300000000000)

# TGFb prediction
df <- data.frame(conditions=c("TGFb","TGFb + low Midostaurin","TGFb + medium Midostaurin","TGFb + high Midostaurin","Midostaurin"),
                 predicted=TGFb_p)

ggplot(df, aes(x=conditions,y=predicted)) +
  geom_bar(width=0.8, position=position_dodge(), stat="identity", colour='black', fill='#3182bd') +
  scale_x_discrete(limits = c("TGFb","TGFb + low Midostaurin","TGFb + medium Midostaurin","TGFb + high Midostaurin","Midostaurin"),
                   labels = function(x) str_wrap(x,width = 10)) +
  scale_fill_grey() + 
  theme_Bryan()

# NE prediction
df <- data.frame(conditions=c("NE","NE + low Midostaurin","NE + medium Midostaurin","NE + high Midostaurin","Midostaurin"),
                 predicted=NE_p)

ggplot(df, aes(x=conditions,y=predicted)) +
  geom_bar(width=0.8, position=position_dodge(), stat="identity", colour='black', fill='#3182bd') +
  scale_x_discrete(limits = c("NE","NE + low Midostaurin","NE + medium Midostaurin","NE + high Midostaurin","Midostaurin"),
                   labels = function(x) str_wrap(x,width = 10)) +
  scale_fill_grey() + 
  theme_Bryan()

#
#
#
#
##cleaning it up##
#
#
#
#


