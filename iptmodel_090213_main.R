# Model of IPT in Brazil
# David Dowdy
# Revision date: Sept. 3, 2013

############################################################################
# DEFINE USER PARAMETERS
# Entries & TB transmission
extrapars <-c(TP = 100000,     # total population size
              LTBI = 0.08,     # prevalence of latent TB in new entries
              MDR = 0.005,     # prevalence of MDR-TB (among all TB) in new entries
              B = 0.00006544,   # transmissibility of TB
              RI = (0.43+0.149*0.57)/(0.63+0.149*0.37), # relative infectiousness of HIV-pos TB
              RIR = 0.647,       # relative transmissibility of DR-TB
              LP = 0.56)       # protection against reinfection if latently infected

# Initial population -> make sure this agrees w/ extrapars
total_pop <- 100000    # size of initial population
latent_prev <- 0.08    # prevalence of latent TB on entering the population
mdr_prev <- 0.005      # prevalence of MDR-TB among active TB entering the population

# TB progression
prop_rapid <- 0.0866   # proportion of infections leading to rapid progression
prop_rapidh <- 0.998    # same, but for HIV, CD4>350
lat_prog <- 0.0005      # progression of latent TB, HIV-neg
lat_progh <- 0.100       # progression of latent TB, HIV, CD4>350
art_prot <- 0.998         # relative reduction in TB incidence from ART (Williams PNAS)
ipt_prot <- (1-0.67)    # proportion reduction in TB progression, IPT

# Mortality
m_notb <- 0.026        # mortality and exit from the population, no TB
m_h1 <- 0.0306            # mortality from HIV, CD4>350 (50% protection)
m_h2 <- 0.0612            # mortality from HIV, CD4<350
m_h1a <- 0.0306          # mortality for HIV on ART, CD4>350 (50% protection)
m_h2a <- 0.0306          # mortality for HIV on ART, CD4<350
m_tb <- 0.0376         # mortality from smear-positive TB, no HIV
m_tbh1 <- 0.0624         # mortality from smear-positive TB, CD4>350 (halfway between)
m_tbh2 <- 0.0872         # mortality from smear-positive TB, CD4<350

# Drug Resistance
pdr_ipt <- 0.00442        # proportion of IPT regimens that generate DR-TB
pdr_tx  <- 0.00442        # proportion of treatments that generate DR-TB

# HIV
hiv_inc <- 0.0003800    # annual incidence of HIV
hiv_prog12 <- 0.2      # progression from CD4 >350 to <350
art_inith1 <- 0        # ART initiation, CD4 >350, set to 0
art_inith2 <- 1.035    # ART initiation, CD4 <350, fit to ART prev
art_fail   <- 0.2      # ART failure, assume 20% of those starting

# TB Treatment, Including IPT
ipt_rate <- 0          # rate of placing HIV-negative pts on IPT
ipt_rateh <- 1/3.93      # IPT rate, HIV, baseline 1/3.93
tx_rate <- 0.87             # treatment rate smear-positive TB, no HIV
rl <- 0.01320             # rate of relapse, treated TB
rlr <- rl             # rate of relapse, treated MDR-TB



#####################################################################
# Set initial matrices and population size (will change this after coming to equilib)
# Reminder of the Matrix Structure:
# matrix = 9 x 5
# rows: 9 TB states (1=S, 2=L, 3=Lr, 4=I, 5=Ir, 6=A, 7=Ar, 8=T, 9=Tr)
# cols: 5 HIV states (1=no HIV, 2=1, 3=1a, 4=2, 5=2a)

# Start all the various matrices off as empty:
statemat <- array(total_pop, dim=c(9,5)) # matrix of initial population sizes
mortmat <- array(0, dim=c(9,5))       # matrix of mortality rates
entrymat <- array(0, dim=c(9,5))      # matrix of entries into the population
mastermat <- array(0, dim=c(9,5,9,5)) # master transmission matrix, 1st 2 dimensions: loss, next 2:gain

# Set the initial population states:
# statemat[1,] <- statemat[1,]*(0.85-0.003-0.001-0.00003-0.0015)
# statemat[2,] <- statemat[2,]*0.15
# statemat[3,] <- statemat[3,]*0.003
# statemat[4,] <- statemat[4,]*0
# statemat[5,] <- statemat[5,]*0
# statemat[6,] <- statemat[6,]*0.001
# statemat[7,] <- statemat[7,]*0.00003
# statemat[8,] <- statemat[8,]*0.001
# statemat[9,] <- statemat[9,]*0.0005
# statemat[,1] <- statemat[,1]*0.994           # 0.6% adult HIV prev, per UNAIDS
# statemat[,2] <- statemat[,2]*0.006*0.2*0.8   # 20% not on ART, 80% by CD4 count
# statemat[,3] <- statemat[,3]*0.006*0.8*0.8
# statemat[,4] <- statemat[,4]*0.006*0.2*0.2
# statemat[,5] <- statemat[,5]*0.006*0.8*0.2

# now, prepare them as a list for entry into ode:
# statename<-paste("S",1:45,sep="")
# statevect <-c(statemat)
# names(statevect)<-statename

################################################
# Mortality
# matrix = 9 x 5
# rows: 9 TB states (1=S, 2=L, 3=Lr, 4=I, 5=Ir, 6=A, 7=Ar, 8=T, 9=Tr)
# cols: 5 HIV states (1=no HIV, 2=1, 3=1a, 4=2, 5=2a)
mortmat[,] <- -m_notb
mortmat[,2] <- mortmat[,2] - m_h1
mortmat[,3] <- mortmat[,3] - m_h1a
mortmat[,4] <- mortmat[,4] - m_h2
mortmat[,5] <- mortmat[,5] - m_h2a
mortmat[6:7,1] <- mortmat[6:7,1] - m_tb
mortmat[6:7,2:3] <- mortmat[6:7,2:3] - m_tbh1
mortmat[6:7,4:5] <- mortmat[6:7,4:5] - m_tbh2

mortname<-paste("m",1:45,sep="")
mortvect <-c(mortmat)
names(mortvect)<-mortname

i <- rep(1:5,each=9)
j <- rep(1:9,5)
k <- array(c(j,i),dim=c(45,2))
m <- array(c(k,k),dim=c(45,4))
mastermat[m]<-mastermat[m]+mortvect


################################################
# Entry
# matrix = 9 x 5
# rows: 9 TB states (1=S, 2=L, 3=Lr, 4=I, 5=Ir, 6=A, 7=Ar, 8=T, 9=Tr)
# cols: 5 HIV states (1=no HIV, 2=1, 3=1a, 4=2, 5=2a)

# will need to multiply these equations by total mortality, not population size:
entrymat[1,1] <- (1-latent_prev)
entrymat[2,1] <- latent_prev*(1-mdr_prev)
entrymat[3,1] <- latent_prev*mdr_prev

entryname<-paste("E",1:45,sep="")
entryvect <-c(entrymat)
names(entryvect)<-entryname

# mastermat[1,1,1,1]<-mastermat[1,1,1,1]+entrymat[1,1]
# mastermat[2,1,2,1]<-mastermat[2,1,2,1]+entrymat[2,1]
# mastermat[3,1,3,1]<-mastermat[3,1,3,1]+entrymat[3,1]



################################################
# TB Progression (Latent -> Active)
# matrix = 9 x 5
# rows: 9 TB states (1=S, 2=L, 3=Lr, 4=I, 5=Ir, 6=A, 7=Ar, 8=T, 9=Tr)
# cols: 5 HIV states (1=no HIV, 2=1, 3=1a, 4=2, 5=2a)

mastermat[c(2,4),1,6,1]<-lat_prog               # from L/I to A
mastermat[c(2,4),2,6,2]<-lat_progh
mastermat[c(2,4),3,6,3]<-lat_progh*art_prot
mastermat[c(2,4),4,6,4]<-lat_progh
mastermat[c(2,4),5,6,5]<-lat_progh*art_prot

mastermat[c(3,5),1,7,1]<-lat_prog               # from Lr/Ir to Ar
mastermat[c(3,5),2,7,2]<-lat_progh
mastermat[c(3,5),3,7,3]<-lat_progh*art_prot
mastermat[c(3,5),4,7,4]<-lat_progh
mastermat[c(3,5),5,7,5]<-lat_progh*art_prot

# add in IPT protection
i <- array(c(rep(4,5),(1:5),rep(6,5),(1:5)),dim=c(5,4))
mastermat[i]<-mastermat[i]*ipt_prot


################################################
# Treatment and Drug Resistance
# matrix = 9 x 5
# rows: 9 TB states (1=S, 2=L, 3=Lr, 4=I, 5=Ir, 6=A, 7=Ar, 8=T, 9=Tr)
# cols: 5 HIV states (1=no HIV, 2=1, 3=1a, 4=2, 5=2a)

mastermat[6,1,8,1]<-tx_rate*(1-pdr_tx)                # from A to T
mastermat[6,2,8,2]<-tx_rate*(1-pdr_tx)
mastermat[6,3,8,3]<-tx_rate*(1-pdr_tx)
mastermat[6,4,8,4]<-tx_rate*(1-pdr_tx)
mastermat[6,5,8,5]<-tx_rate*(1-pdr_tx)

mastermat[6,1,7,1]<-tx_rate*pdr_tx                    # from A to Ar
mastermat[6,2,7,2]<-tx_rate*pdr_tx
mastermat[6,3,7,3]<-tx_rate*pdr_tx
mastermat[6,4,7,4]<-tx_rate*pdr_tx
mastermat[6,5,7,5]<-tx_rate*pdr_tx

mastermat[7,1,9,1]<-tx_rate                          # from Ar to Tr
mastermat[7,2,9,2]<-tx_rate
mastermat[7,3,9,3]<-tx_rate
mastermat[7,4,9,4]<-tx_rate
mastermat[7,5,9,5]<-tx_rate

# Relapse:
mastermat[8,1,6,1]<-rl                      # from T to A
mastermat[8,2,6,2]<-rl
mastermat[8,3,6,3]<-rl*art_prot
mastermat[8,4,6,4]<-rl
mastermat[8,5,6,5]<-rl*art_prot

mastermat[9,1,7,1]<-rlr                      # from Tr to Ar
mastermat[9,2,7,2]<-rlr
mastermat[9,3,7,3]<-rlr*art_prot
mastermat[9,4,7,4]<-rlr
mastermat[9,5,7,5]<-rlr*art_prot

# IPT:
mastermat[2,1,4,1]<-ipt_rate*(1-pdr_ipt)              # from L to I
mastermat[2,2,4,2]<-ipt_rateh*(1-pdr_ipt)
mastermat[2,3,4,3]<-ipt_rateh*(1-pdr_ipt)
mastermat[2,4,4,4]<-ipt_rateh*(1-pdr_ipt)
mastermat[2,5,4,5]<-ipt_rateh*(1-pdr_ipt)

mastermat[2,1,5,1]<-ipt_rate*pdr_ipt                  # from L to Ir
mastermat[2,2,5,2]<-ipt_rateh*pdr_ipt
mastermat[2,3,5,3]<-ipt_rateh*pdr_ipt
mastermat[2,4,5,4]<-ipt_rateh*pdr_ipt
mastermat[2,5,5,5]<-ipt_rateh*pdr_ipt

mastermat[3,1,5,1]<-ipt_rate                          # from Lr to Ir
mastermat[3,2,5,2]<-ipt_rateh
mastermat[3,3,5,3]<-ipt_rateh
mastermat[3,4,5,4]<-ipt_rateh
mastermat[3,5,5,5]<-ipt_rateh


################################################
# HIV Epi & Treatment
# matrix = 9 x 5
# rows: 9 TB states (1=S, 2=L, 3=Lr, 4=I, 5=Ir, 6=A, 7=Ar, 8=T, 9=Tr)
# cols: 5 HIV states (1=no HIV, 2=1, 3=1a, 4=2, 5=2a)

# HIV infection & progression:
i1 <- array(c((1:9),rep(1,9),(1:9),rep(2,9)),dim=c(9,4)) 
mastermat[i1]<-hiv_inc                        # from no HIV to HIV1
i2 <- array(c((1:9),rep(2,9),(1:9),rep(4,9)),dim=c(9,4))
mastermat[i2]<-hiv_prog12                     # from HIV1 to HIV2

# Starting ART:
i4 <- array(c((1:9),rep(2,9),(1:9),rep(3,9)),dim=c(9,4)) 
mastermat[i4]<-art_inith1                                     # from HIV1 to HIV1a
i5 <- array(c((1:9),rep(4,9),(1:9),rep(5,9)),dim=c(9,4)) 
mastermat[i5]<-art_inith2                                     # from HIV2 to HIV2a

# ART Failure:
i7 <- array(c(rep((1:9),2),rep(3,9),rep(5,9),rep((1:9),2),rep(2,9),rep(4,9)),dim=c(18,4)) 
mastermat[i7]<-art_fail


###################################################################
# Function for input into the differential equation solver:
# first, put mastermat into a useable format
mastername<-paste("M",1:2025,sep="")  # 45x45 parameters
mastervect <-c(mastermat)
names(mastervect)<-mastername
mastervect<-c(mastervect,extrapars)

TBdx<-function(t, state, parameters) {
        x<-as.list(c(state,parameters))
        with(x,{
             statev <- unlist(x[1:45])
             masterv <- unlist(x[46:(46+45*45-1)])

             # force of infection
             force <- B*(S6)+B*RI*(S15+S24+S33+S42)
             forcedr <- B*RIR*(S7)+B*RI*RIR*(S16+S25+S34+S43)

             # now, create an "infection matrix"
             infectmat <- array(0, dim=c(9,5,9,5)) 
             infectmat[1,1,2,1]<-(1-prop_rapid)*force           # from S to L
             infectmat[1,2,2,2]<-(1-prop_rapidh)*force
             infectmat[1,3,2,3]<-(1-prop_rapidh*art_prot)*force
             infectmat[1,4,2,4]<-(1-prop_rapidh)*force
             infectmat[1,5,2,5]<-(1-prop_rapidh*art_prot)*force

             infectmat[1,1,3,1]<-(1-prop_rapid)*forcedr           # from S to Lr
             infectmat[1,2,3,2]<-(1-prop_rapidh)*forcedr
             infectmat[1,3,3,3]<-(1-prop_rapidh*art_prot)*forcedr
             infectmat[1,4,3,4]<-(1-prop_rapidh)*forcedr
             infectmat[1,5,3,5]<-(1-prop_rapidh*art_prot)*forcedr

             infectmat[c(1:5,8,9),1,6,1]<-prop_rapid*force       # from S/L/I/T to A
             infectmat[c(1:5,8,9),2,6,2]<-prop_rapidh*force
             infectmat[c(1:5,8,9),3,6,3]<-prop_rapidh*art_prot*force
             infectmat[c(1:5,8,9),4,6,4]<-prop_rapidh*force
             infectmat[c(1:5,8,9),5,6,5]<-prop_rapidh*art_prot*force

             infectmat[c(1:5,8,9),1,7,1]<-prop_rapid*forcedr       # from S/L/I/T to Ar
             infectmat[c(1:5,8,9),2,7,2]<-prop_rapidh*forcedr
             infectmat[c(1:5,8,9),3,7,3]<-prop_rapidh*art_prot*forcedr
             infectmat[c(1:5,8,9),4,7,4]<-prop_rapidh*forcedr
             infectmat[c(1:5,8,9),5,7,5]<-prop_rapidh*art_prot*forcedr

             # add in protection from LTBI
             i <- array(c(rep(2:5,each=2),rep(1,8),rep(6:7,4),rep(1,8)),dim=c(8,4))
             infectmat[i]<-infectmat[i]*LP

             # convert infectmat into a vector
             infectv<-c(infectmat)

             # add infectv to masterv
             masterv3<-masterv+infectv

             # sum deaths, so entries can be added
             i3<-array(1:45,dim=c(45,2))
             i4<-1:45

             # create the change matrix
             mastermatr <-array(masterv3,dim=c(45,45))
             deaths<-sum(mastermatr[i3]*statev[i4])

             # rate of change
             delta <- rep(0,45)
             delta[1] <- sum(mastermatr[,1]*statev) - sum(mastermatr[1,]*statev[1]) + mastermatr[1,1]*statev[1] -
                   (1-LTBI)*deaths
             delta[2] <- sum(mastermatr[,2]*statev) - sum(mastermatr[2,]*statev[2]) + mastermatr[2,2]*statev[2] -
                   (LTBI*(1-MDR))*deaths
             delta[3] <- sum(mastermatr[,3]*statev) - sum(mastermatr[3,]*statev[3]) + mastermatr[3,3]*statev[3] -
                   (LTBI*MDR*deaths)
             delta[4] <- sum(mastermatr[,4]*statev) - sum(mastermatr[4,]*statev[4]) + mastermatr[4,4]*statev[4]
             delta[5] <- sum(mastermatr[,5]*statev) - sum(mastermatr[5,]*statev[5]) + mastermatr[5,5]*statev[5]
             delta[6] <- sum(mastermatr[,6]*statev) - sum(mastermatr[6,]*statev[6]) + mastermatr[6,6]*statev[6]
             delta[7] <- sum(mastermatr[,7]*statev) - sum(mastermatr[7,]*statev[7]) + mastermatr[7,7]*statev[7]
             delta[8] <- sum(mastermatr[,8]*statev) - sum(mastermatr[8,]*statev[8]) + mastermatr[8,8]*statev[8]
             delta[9] <- sum(mastermatr[,9]*statev) - sum(mastermatr[9,]*statev[9]) + mastermatr[9,9]*statev[9]
             delta[10] <- sum(mastermatr[,10]*statev) - sum(mastermatr[10,]*statev[10]) + mastermatr[10,10]*statev[10]
             delta[11] <- sum(mastermatr[,11]*statev) - sum(mastermatr[11,]*statev[11]) + mastermatr[11,11]*statev[11]
             delta[12] <- sum(mastermatr[,12]*statev) - sum(mastermatr[12,]*statev[12]) + mastermatr[12,12]*statev[12]
             delta[13] <- sum(mastermatr[,13]*statev) - sum(mastermatr[13,]*statev[13]) + mastermatr[13,13]*statev[13]
             delta[14] <- sum(mastermatr[,14]*statev) - sum(mastermatr[14,]*statev[14]) + mastermatr[14,14]*statev[14]
             delta[15] <- sum(mastermatr[,15]*statev) - sum(mastermatr[15,]*statev[15]) + mastermatr[15,15]*statev[15]
             delta[16] <- sum(mastermatr[,16]*statev) - sum(mastermatr[16,]*statev[16]) + mastermatr[16,16]*statev[16]
             delta[17] <- sum(mastermatr[,17]*statev) - sum(mastermatr[17,]*statev[17]) + mastermatr[17,17]*statev[17]
             delta[18] <- sum(mastermatr[,18]*statev) - sum(mastermatr[18,]*statev[18]) + mastermatr[18,18]*statev[18]
             delta[19] <- sum(mastermatr[,19]*statev) - sum(mastermatr[19,]*statev[19]) + mastermatr[19,19]*statev[19]
             delta[20] <- sum(mastermatr[,20]*statev) - sum(mastermatr[20,]*statev[20]) + mastermatr[20,20]*statev[20]
             delta[21] <- sum(mastermatr[,21]*statev) - sum(mastermatr[21,]*statev[21]) + mastermatr[21,21]*statev[21]
             delta[22] <- sum(mastermatr[,22]*statev) - sum(mastermatr[22,]*statev[22]) + mastermatr[22,22]*statev[22]
             delta[23] <- sum(mastermatr[,23]*statev) - sum(mastermatr[23,]*statev[23]) + mastermatr[23,23]*statev[23]
             delta[24] <- sum(mastermatr[,24]*statev) - sum(mastermatr[24,]*statev[24]) + mastermatr[24,24]*statev[24]
             delta[25] <- sum(mastermatr[,25]*statev) - sum(mastermatr[25,]*statev[25]) + mastermatr[25,25]*statev[25]
             delta[26] <- sum(mastermatr[,26]*statev) - sum(mastermatr[26,]*statev[26]) + mastermatr[26,26]*statev[26]
             delta[27] <- sum(mastermatr[,27]*statev) - sum(mastermatr[27,]*statev[27]) + mastermatr[27,27]*statev[27]
             delta[28] <- sum(mastermatr[,28]*statev) - sum(mastermatr[28,]*statev[28]) + mastermatr[28,28]*statev[28]
             delta[29] <- sum(mastermatr[,29]*statev) - sum(mastermatr[29,]*statev[29]) + mastermatr[29,29]*statev[29]
             delta[30] <- sum(mastermatr[,30]*statev) - sum(mastermatr[30,]*statev[30]) + mastermatr[30,30]*statev[30]
             delta[31] <- sum(mastermatr[,31]*statev) - sum(mastermatr[31,]*statev[31]) + mastermatr[31,31]*statev[31]
             delta[32] <- sum(mastermatr[,32]*statev) - sum(mastermatr[32,]*statev[32]) + mastermatr[32,32]*statev[32]
             delta[33] <- sum(mastermatr[,33]*statev) - sum(mastermatr[33,]*statev[33]) + mastermatr[33,33]*statev[33]
             delta[34] <- sum(mastermatr[,34]*statev) - sum(mastermatr[34,]*statev[34]) + mastermatr[34,34]*statev[34]
             delta[35] <- sum(mastermatr[,35]*statev) - sum(mastermatr[35,]*statev[35]) + mastermatr[35,35]*statev[35]
             delta[36] <- sum(mastermatr[,36]*statev) - sum(mastermatr[36,]*statev[36]) + mastermatr[36,36]*statev[36]
             delta[37] <- sum(mastermatr[,37]*statev) - sum(mastermatr[37,]*statev[37]) + mastermatr[37,37]*statev[37]
             delta[38] <- sum(mastermatr[,38]*statev) - sum(mastermatr[38,]*statev[38]) + mastermatr[38,38]*statev[38]
             delta[39] <- sum(mastermatr[,39]*statev) - sum(mastermatr[39,]*statev[39]) + mastermatr[39,39]*statev[39]
             delta[40] <- sum(mastermatr[,40]*statev) - sum(mastermatr[40,]*statev[40]) + mastermatr[40,40]*statev[40]
             delta[41] <- sum(mastermatr[,41]*statev) - sum(mastermatr[41,]*statev[41]) + mastermatr[41,41]*statev[41]
             delta[42] <- sum(mastermatr[,42]*statev) - sum(mastermatr[42,]*statev[42]) + mastermatr[42,42]*statev[42]
             delta[43] <- sum(mastermatr[,43]*statev) - sum(mastermatr[43,]*statev[43]) + mastermatr[43,43]*statev[43]
             delta[44] <- sum(mastermatr[,44]*statev) - sum(mastermatr[44,]*statev[44]) + mastermatr[44,44]*statev[44]
             delta[45] <- sum(mastermatr[,45]*statev) - sum(mastermatr[45,]*statev[45]) + mastermatr[45,45]*statev[45]

             # return the rate of change
             list(c(delta))
             })
}


outgen<-function(statevect, mastervect) {
        x<-as.list(c(statevect,mastervect))
        with(x,{
             statev <- unlist(x[1:45])
		 rl <- unlist(x[46])
		 rlr <- unlist(x[47])
		 tx_rate <- unlist(x[48])
             masterv <- unlist(x[49:(49+45*45-1)])
             statename<-paste("S",1:45,sep="")
             names(statev)<-statename

             # force of infection
             force <- B*(S6)+B*RI*(S15+S24+S33+S42)
             forcedr <- B*RIR*(S7)+B*RI*RIR*(S16+S25+S34+S43)

             # now, create an "infection matrix"
             # now, create an "infection matrix"
             infectmat <- array(0, dim=c(9,5,9,5)) 
             infectmat[1,1,2,1]<-(1-prop_rapid)*force           # from S to L
             infectmat[1,2,2,2]<-(1-prop_rapidh)*force
             infectmat[1,3,2,3]<-(1-prop_rapidh*art_prot)*force
             infectmat[1,4,2,4]<-(1-prop_rapidh)*force
             infectmat[1,5,2,5]<-(1-prop_rapidh*art_prot)*force

             infectmat[1,1,3,1]<-(1-prop_rapid)*forcedr           # from S to Lr
             infectmat[1,2,3,2]<-(1-prop_rapidh)*forcedr
             infectmat[1,3,3,3]<-(1-prop_rapidh*art_prot)*forcedr
             infectmat[1,4,3,4]<-(1-prop_rapidh)*forcedr
             infectmat[1,5,3,5]<-(1-prop_rapidh*art_prot)*forcedr

             infectmat[c(1:5,8,9),1,6,1]<-prop_rapid*force       # from S/L/I/T to A
             infectmat[c(1:5,8,9),2,6,2]<-prop_rapidh*force
             infectmat[c(1:5,8,9),3,6,3]<-prop_rapidh*art_prot*force
             infectmat[c(1:5,8,9),4,6,4]<-prop_rapidh*force
             infectmat[c(1:5,8,9),5,6,5]<-prop_rapidh*art_prot*force

             infectmat[c(1:5,8,9),1,7,1]<-prop_rapid*forcedr       # from S/L/I/T to Ar
             infectmat[c(1:5,8,9),2,7,2]<-prop_rapidh*forcedr
             infectmat[c(1:5,8,9),3,7,3]<-prop_rapidh*art_prot*forcedr
             infectmat[c(1:5,8,9),4,7,4]<-prop_rapidh*forcedr
             infectmat[c(1:5,8,9),5,7,5]<-prop_rapidh*art_prot*forcedr

             # add in protection from LTBI
             i <- array(c(rep(2:5,each=2),rep(1,8),rep(6:7,4),rep(1,8)),dim=c(8,4))
             infectmat[i]<-infectmat[i]*LP

             # convert infectmat into a vector
             infectv<-c(infectmat)

             # add infectv to masterv
             masterv3<-masterv+infectv

             # sum deaths, so entries can be added
             i3<-array(1:45,dim=c(45,2))
             i4<-1:45

             # create the change matrix
             mastermatr <-array(masterv3,dim=c(45,45))
             deaths<-sum(mastermatr[i3]*statev[i4])

             # remove losses to the population for death calc (not in main code)
             mastermatr[i3]<-mastermatr[i3]+0.022

		 i3<-c(seq(6,46,by=9),seq(7,46,by=9))
		 incidence <- sum(mastermatr[,i3]*statev)

		 i3b<-c(seq(15,46,by=9),seq(16,46,by=9))
		 hivtbinc <- sum(mastermatr[,i3b]*statev)

		 i3c<-array(c(i3,i3),dim=c(10,2))
		 TBmort <--sum(mastermatr[i3c]*statev[i3])

		 i4<-array(c(6:7),dim=c(2,2))
		 i5<-c(6:7)
 		 nonHIVTBmort<--sum(mastermatr[i4]*statev[i5])
		 HIVTBmort <-TBmort-nonHIVTBmort

		 i6<-c(37:45)
		 i7<-c(28:45)
		 ARTprop<-sum(statev[i6])/sum(statev[i7])*100

		 i8<-c(seq(6,46,by=9),seq(7,46,by=9))
		 TBprev<-sum(statev[i8])

		 i9<-c(seq(7,46,by=9))
		 MDRinc <- sum(mastermatr[,i9]*statev)

		 i10<-c(seq(10,45))
		 HIVprev<-sum(statev[i10])
		 i10b<-array(c(i10,i10),dim=c(36,2))
		 HIVmort <--sum(mastermatr[i10b]*statev[i10])

		 i11<-c(seq(8,46,by=9),seq(9,46,by=9))
		 i11b<-array(c(rep(seq(8,46,by=9),2),rep(seq(9,46,by=9),2),rep(i3,2)),dim=c(20,2))
		 relapse<- sum(mastermatr[i11b]*statev[i11])/incidence*100*(tx_rate-0.2*(incidence-hivtbinc)/incidence)/tx_rate
		 i12<-c(seq(8,46,by=9),seq(9,46,by=9))
		 i12b<-array(c(i12,i9,i9),dim=c(10,2))
		 i13<-c(seq(6,46,by=9))
		 i12c<-array(c(i13,i9),dim=c(5,2))
		 mdrrelapse <-sum(mastermatr[i12b]*statev[i12])/MDRinc*100*(tx_rate-0.2*(incidence-hivtbinc)/incidence)/tx_rate + 
			sum(mastermatr[i12c]*statev[i13])/MDRinc*100
		 i50<-c("TB inc", "HIVTB inc", "TB mort", "HIVTB mort",
		        "ART prop", "TB prev", "MDRTB inc", "HIV prev", "HIV mort",
		        "prop retreat", "prop retreat MDR",
		       incidence, hivtbinc, TBmort, HIVTBmort,
		       ARTprop, TBprev, MDRinc, HIVprev, HIVmort, relapse, mdrrelapse,
		       122, 15.74, 6.5, 1.54, 80, 132, 1.87, 600, 24.6, 27.4, 51.8)
		 i51<-array(i50,dim=c(11,3))
		 i51

		 })
}



# Times at which the model will be evaluated (0.01-year steps, for 5 years):
times <- seq(0, 5, by = 0.1)

# Bring to equilib: statevect as last output
 infile<-(read.csv("C:\\Users\\Family\\Documents\\R\\ipt_output_apr222012.csv"))
 statevect <-c(infile[51,])
 statevect <-c(as.numeric(statevect[3:47]))
 statename<-paste("S",1:45,sep="")
 names(statevect)<-statename
library(deSolve)



out <- ode(statevect, times, TBdx, mastervect)
out2 <- c(out[50,2:46],rl,rlr,tx_rate)
outvect <- outgen(out2,mastervect)
outdata <-as.data.frame(out)

T <- out[,1]
S <- out[,seq(2,46,by=9)]
L <- out[,c(seq(3,46,by=9),seq(4,46,by=9),seq(5,46,by=9),seq(6,46,by=9))]
A <- out[,c(seq(7,46,by=9),seq(8,46,by=9))]
R <- out[,c(seq(9,46,by=9),seq(10,46,by=9))]
H <- out[,c(seq(11,46))]
M <- out[,c(seq(8,46,by=9))]

total <-array(c(rowSums(out[,2:46])))

Ss <-array(c(rowSums(S)))
Ls <-array(c(rowSums(L)))
As <-array(c(rowSums(A)))
Rs <-array(c(rowSums(R)))
Hs <-array(c(rowSums(H)))
Ms <-array(c(rowSums(M)))

output <- as.data.frame(array(c(T,Ss,Ls,As,Rs,Hs,Ms),dim=c(length(T),7)))
write.csv(outdata,file="C:\\Users\\Family\\Documents\\R\\ipt_output_apr22eq2012.csv")

#############################################################
# Plot results

par(mfcol=c(3,2),
mfrow=c(3,2))
plot(T,Ss)
plot(T,Ls)
plot(T,As)
plot(T,Rs)
plot(T,Hs)
plot(T,Ms)


##############################################################
# outmat <-array(0,dim=c(50,10))
# sensanal <-array(0,dim=c(200,2))

# SENSITIVITY ANALYSES:
# Determining the range over which each parameter needs to be varied
# First, sweep over 200 parameter values:

# for (z in 1:40) {
# 	extrapara <-c(TP = 100000,     # total population size
#       	        LTBI = 0.08,     # prevalence of latent TB in new entries
#             	  MDR = 0.005,     # prevalence of MDR-TB (among all TB) in new entries
# 	              B = 0.0000955,   # transmissibility of TB
#       	        RI = 0.149,      # relative infectiousness of smear-negatives
#             	  RIR = 0.749,       # relative transmissibility of DR-TB
# 	              LP = 0.56 + (z)*0.01)       # protection against reinfection if latently infected
# 
# # Insert parameter here, then make sure to change sensanal output:
# 
# 
# 	# Bring to equilib: statevect as last output
# 	infile<-(read.csv("C:\\Users\\Family\\Documents\\R\\ipt_output_apr52012.csv"))
# 	statevect <-c(infile[51,])
# 	statevect <-c(as.numeric(statevect[3:79]))
# 	statename<-paste("S",1:77,sep="")
# 	names(statevect)<-statename
# 	library(deSolve)
# 
# 
#       mastername<-paste("M",1:5929,sep="")  # 77x77 parameters
#       mastervect <-c(mastermat)
#       names(mastervect)<-mastername
#       mastervect<-c(mastervect,extrapara)
# 
#       outsa <- ode(statevect, times, TBdx, mastervect)
# 	outsa2 <- c(outsa[50,2:78],rl,rlr,tx_rate)
#       outvect <- outgen(outsa2,mastervect)
#       sensanal[z,] <- c(outvect[1],extrapara[7])
# }
# 
# # give the results of the parameter sweep:
# sensanal[1:40,]

# # write.csv(outdata,file="C:\\Users\\Family\\Documents\\R\\ipt_output_apr192012.csv")

outvect