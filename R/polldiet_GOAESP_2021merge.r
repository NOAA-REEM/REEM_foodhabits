setwd("C:/src/REEM_foodhabits")

ppfile       <- "data/GOA_raw.csv.gz"                    # raw Pred/prey diet file (RACEBASE query result)
pplookfile   <- "lookups/Alaska_PreyLookup_MASTER.csv"  # lookup file for grouping prey
preylook_col <- "GOA_MED"                      # Column name in pplookfile to use to group prey
min_sample   <- 5                               # Minimum sample size 

# List of predator-specific values.  
preds <- list(
  "W.pollock"  = list(nodc="8791030701", A_L=0.00553096, B_L=3.044172,   LCLASS=c(0,10,25,40,55,999) ),
  "P.cod"      = list(nodc="8791030401", A_L=0.00411781, B_L=3.25325765, LCLASS=c(0,10,30,60,85,999) ),
  "Arrowtooth" = list(nodc=c("8857040100", "8857040102"), # ATF includes Atheresthes sp. unid
                                         A_L=0.00443866, B_L=3.19894001, LCLASS=c(0,10,30,50,999) ),
  "P.halibut"  = list(nodc="8857041901", A_L=0.01093251, B_L=3.24,       LCLASS=c(0,10,50,70,999) )
)

predators <- c("W.pollock")
# Read in lookup and raw data and ensure 10-digit nodc numbers are read as text keys
preylooktable <- read.csv(pplookfile)
  preylook <- preylooktable[,preylook_col]
  names(preylook)<-sprintf("%010.0f",preylooktable$NODC_CODE)
  fullprey<-unique(preylooktable[,preylook_col])
rawdat   <- read.csv(ppfile)
  rawdat$PREY_NODC <-sprintf("%010.0f",rawdat$PREY_NODC)
  rawdat$PRED_NODC <-sprintf("%010.0f",rawdat$PRED_NODC)  
  rawdat$preyguild <- preylook[rawdat$PREY_NODC]


  # CENT_WEST
  stratblock <- list("WEST_CENTRAL" = c(10,11,12,13,20,21,22,30,31,32,33,35,
                        110,111,112,120,121,122,130,131,132,133,134,
                        210,220,221,230,231,232,310,320,330,410,420,430,510,520,530)
  )
  
  yearlist <- c(1990, 1993, 1996, 1999, 2001, 2003, 2005, 2007, 2009, 2011, 2013, 2017, 2019, 2021)  
  
  
outdat <- e_logvals <- e_sdvals <- e_vals <- o_vals <- cp_vals <- cper_dat <- NULL   

# Start Calcs here

for (PRED in predators){ 
ThisPred <- preds[[PRED]] 

preddat <- rawdat[rawdat$PRED_NODC %in% ThisPred$nodc,] 

# Create a crosstab table SURV and list of prey allprey using the preylook_col as the crosstab
gtab <- tapply(preddat$TWT,list(as.character(preddat$PREDJOIN),preddat$preyguild),sum)
gtab[is.na(gtab)]<-0
predtab <- unique(preddat[,c("PREDJOIN","PRED_NODC","PRED_NAME","PRED_LEN","HAULJOIN","REGION","STRATUM","YEAR","MONTH","DAY","GEAR_TEMP",
                             "GEAR_DEPTH","BOTTOM_DEPTH","SURFACE_TEMP","STATIONID","CRUISE_TYPE","RLAT","RLONG")])
preycross <- gtab[as.character(predtab$PREDJOIN),]
# Bin by pred LCLASS (in cm)
lbin  <- cut(predtab$PRED_LEN,ThisPred$LCLASS,right=F)
GDAT    <-cbind(predtab,lbin,preycross)
allprey <- colnames(preycross)
alllen  <- unique(as.character(lbin))

SURV <- GDAT[GDAT$CRUISE_TYPE=="Race_Groundfish" & GDAT$MONTH%in%6:8, ]

# add any additional filters for the data here:     
for (STRAT in names(stratblock)){
  cat(PRED,STRAT,"\n"); flush.console()
  stratcodes <- unlist(stratblock[STRAT])
  
for (LL in alllen){
  lencodes <- LL 
  
  for (YY in yearlist){
    yearcodes <- YY
    
    SELPRED <- SURV[ SURV$STRATUM %in% stratcodes &
                       SURV$lbin %in% lencodes    &
                       SURV$YEAR %in% yearcodes   ,]      
    
    if (length(SELPRED[,1])>=min_sample){
      allfood <- SELPRED[,allprey]
      goodprey <- allfood #allfood[,colSums(allfood)>0] to remove colums with 0 for dirichlet  
      sptot    <- colSums(goodprey)
      spadd <- 0.0 #spadd    <- DETECT*sptot/sum(sptot)
      SCI <- t(t(goodprey)+spadd)/ (ThisPred$A_L * (SELPRED$PRED_LEN ^ ThisPred$B_L))
      sp_prey  <- colnames(goodprey)
      Pfull <- colSums(goodprey>0)  
      Pnum  <- colSums(goodprey>=0)
      elog  <- log((Pfull+0.5)/(Pnum-Pfull+0.5))
      e_sd  <- sqrt(1.0/(Pfull+0.5) + 1.0/(Pnum-Pfull+0.5))
      
      
      # Make matrix of individual stomachs
      sampmat <- matrix(as.numeric(unlist(SCI)),length(SCI[,1]),length(sp_prey))
      colnames(sampmat)<-sp_prey
      Nind <- length(sampmat[,1])    
      IND  <- 1:Nind #IND <- sample.int(Nind, Nind, replace = T)  
      tot_wt    <- sum(rowSums(sampmat[IND,]))
      cperw     <- sum(rowSums(sampmat[IND,]))/Nind
      pathsum <- colSums(sampmat[IND,]) #+ spadd
      dietprop <- pathsum/sum(pathsum)
      cperfull  <- pathsum/Pfull
      # Saving the data     
      Nsamp <- Nind
      Nfull <- sum(rowSums(sampmat[IND,])>0)
      SCIperN <- cperw 
      #outdat    <- rbind(outdat,data.frame(PRED,STRAT,YY,LL,Nsamp,Nfull,SCIperN,t(dietprop)))
      #cper_dat  <- rbind(cper_dat,data.frame(PRED,STRAT,YY,LL,Nsamp,Nfull,SCIperN,t(cperfull)))
      #e_logvals <- rbind(e_logvals,data.frame(PRED,STRAT,YY,LL,Nsamp,Nfull,SCIperN,t(elog) ))
      #e_sdvals  <- rbind(e_sdvals, data.frame(PRED,STRAT,YY,LL,Nsamp,Nfull,SCIperN,t(e_sd) ))
      o_vals    <- rbind(o_vals,data.frame(PRED,STRAT,YY,LL,Nsamp,Nfull,SCIperN,sp_prey,dietprop,cperfull,elog,e_sd))
      #cp_vals   <- rbind(cp_vals,data.frame(PRED,STRAT,YY,LL,Nsamp,Nfull,SCIperN,sp_prey,cperfull))
      #e_vals    <- rbind(e_vals,data.frame(PRED,STRAT,YY,LL,Nsamp,Nfull,SCIperN,sp_prey,elog,e_sd))
      #latmat[latind,sp_prey] <- dietprop[sp_prey]
    } # end of if sample size conditional 
    
    
    #cat(stratcodes," ",lencodes,"\n") 
  } # end of yearlist
} # end of lenblock
} #end of stratblock      
}  # end predator loop 

o_vals$cperfull[is.nan(o_vals$cperfull)]<-0        

# Add transformed logit values for frequency of occurrence 
o_vals$meanper <- exp(o_vals$elog)/(1+exp(o_vals$elog))       
o_vals$lo95per <- exp(o_vals$elog-1.96*o_vals$e_sd)/(1+exp(o_vals$elog-1.96*o_vals$e_sd))
o_vals$hi95per <- exp(o_vals$elog+1.96*o_vals$e_sd)/(1+exp(o_vals$elog+1.96*o_vals$e_sd))

write.csv(o_vals,"results/out_diets_allGOA.csv",row.names=F)





