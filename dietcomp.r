
ppfile       <- "BS_raw.csv"                    # raw Pred/prey diet file (RACEBASE query result)
pplookfile   <- "Alaska_PreyLookup_MASTER.csv"  # lookup file for grouping prey
preylook_col <- "EBS_CRAB"                      # Column name in pplookfile to use to group prey
min_sample   <- 5                               # Minimum sample size 

# List of predator-specific values.  
preds <- list(
  "W.pollock"  = list(nodc="8791030701", A_L=0.00553096, B_L=3.044172,   LCLASS=c(0,10,25,40,55,999) ),
  "P.cod"      = list(nodc="8791030401", A_L=0.00411781, B_L=3.25325765, LCLASS=c(0,10,30,60,85,999) ),
  "Arrowtooth" = list(nodc=c("8857040100", "8857040102"), # ATF includes Atheresthes sp. unid
                                         A_L=0.00443866, B_L=3.19894001, LCLASS=c(0,10,30,50,999) ),
  "P.halibut"  = list(nodc="8857041901", A_L=0.01093251, B_L=3.24,       LCLASS=c(0,10,50,70,999) )
)

# Predators to calculate (all must be on above list)
predators <- c("P.cod","P.halibut")
  
# Years to output
yearlist <- 1985:2019

# Group strata as desired
stratblock <- list(
    "SE_inner"  = c(10),
    "SE_middle" = c(31),
    "SE_outer"  = c(50),
    "Pribs"     = c(32,42),
    "NW_inner"  = c(20),
    "NW_middle" = c(41),
    "NW_outer"  = c(61),
    "StMatt"    = c(43,62),
    "NW_corner" = c(82,90),
    "NBS"       = c(70,71,81)
  )

# Read in lookup and raw data and ensure 10-digit nodc numbers are read as text keys
preylooktable <- read.csv(pplookfile)
  preylook <- preylooktable[,preylook_col]
  names(preylook)<-sprintf("%010.0f",preylooktable$NODC_CODE)
  fullprey<-unique(preylooktable[,preylook_col])
rawdat   <- read.csv(ppfile)
  rawdat$PREY_NODC <-sprintf("%010.0f",rawdat$PREY_NODC)
  rawdat$PRED_NODC <-sprintf("%010.0f",rawdat$PRED_NODC)  
  rawdat$preyguild <- preylook[rawdat$PREY_NODC]

    
# Start main loop here
outdat <- e_logvals <- e_sdvals <- e_vals <- o_vals <- cp_vals <- cper_dat <- NULL     

for (PRED in predators){    
  # Pick the predator 
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
  
  # add any additional filters for the data here:     
  SURV <- GDAT[GDAT$CRUISE_TYPE=="Race_Groundfish" & GDAT$MONTH%in%6:8,]
  
  
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
        
write.csv(o_vals,"out_diets.csv",row.names=F)


raw_original <- read.csv("EBS_Pcod_length_cpue.csv")
PRED <- "P.cod"
ThisPred <- preds[[PRED]] 
preddat <- rawdat[rawdat$PRED_NODC %in% ThisPred$nodc,] 
#Bioenergetics calculations
 source("bioen_ftc_fun.R")
 load("bioen_parUSE.Rdata")
bioen <- parUSE[parUSE$X=="Pcod_kir1",] # for cod


# Start table creation here
raw_cpue <- raw_original
haul_tab <- tapply(raw_cpue$frequency,list(raw_cpue$haul_join),sum)
raw_cpue$tot_freq <- haul_tab[as.character(raw_cpue$haul_join)]

# Calculate daily max consumption (g/predator/day) based on pred length
wt_g <- ThisPred$A_L * ((raw_cpue$fish_length/10)^ThisPred$B_L)
raw_cpue$Cmax <- bioen$CA * (wt_g^(1+bioen$CB))

# need to replace gear_temp NAs with stratum averages
Tfix <- raw_cpue[!is.na(raw_cpue$gear_temp),c("haul_num","stratum","gear_temp")]
tdat <- aggregate(Tfix$gear_temp,list(Tfix$haul_num,Tfix$stratum),mean)
tave <- tdat$x; names(tave) <- paste(tdat$Group.1,tdat$Group.2)
raw_cpue$tclean <- ifelse(is.na(raw_cpue$gear_temp),tave[paste(raw_cpue$haul_num,raw_cpue$stratum)],raw_cpue$gear_temp)

tmp<-list(TempC=raw_cpue$tclean)
raw_cpue$fTC<-fTC_fun(bioen,tmp)

raw_cpue$demand_g_km2_day <- raw_cpue$num_cpue * (raw_cpue$frequency/raw_cpue$tot_freq) * raw_cpue$Cmax * raw_cpue$fTC 


#Diet table lookups
stratlook <- sub("[0-9]","",names(unlist(stratblock)))
names(stratlook) <- unlist(stratblock)
raw_cpue$bigstrat <- stratlook[as.character(raw_cpue$stratum)]

raw_cpue$lbin <- cut(raw_cpue$fish_length/10,ThisPred$LCLASS,right=F)

dietkey <- paste(PRED,raw_cpue$bigstrat, raw_cpue$haul_num, raw_cpue$lbin)

prey <- o_vals[o_vals$sp_prey=="Opilio",]
row.names(prey) <- paste(prey$PRED,prey$STRAT, prey$YY,prey$LL)
raw_cpue$Opilio <- prey[dietkey,"dietprop"]

prey <- o_vals[o_vals$sp_prey=="Bairdi",]
row.names(prey) <- paste(prey$PRED,prey$STRAT, prey$YY,prey$LL)
raw_cpue$Bairdi <- prey[dietkey,"dietprop"]

prey <- o_vals[o_vals$sp_prey=="Unid.Chion",]
row.names(prey) <- paste(prey$PRED,prey$STRAT, prey$YY,prey$LL)
raw_cpue$Unid.Chion <- prey[dietkey,"dietprop"]

prey <- o_vals[o_vals$sp_prey=="W.pollock",]
row.names(prey) <- paste(prey$PRED,prey$STRAT, prey$YY,prey$LL)
raw_cpue$Pollock <- prey[dietkey,"dietprop"]

raw_cpue$Other <- 1-raw_cpue$Opilio-raw_cpue$Bairdi-raw_cpue$Unid.Chion-raw_cpue$Pollock

raw_cpue$Missing <- ifelse(is.na(raw_cpue$Other),1,0)

raw_cpue$Opilio[is.na(raw_cpue$Other)] <- 0
raw_cpue$Bairdi[is.na(raw_cpue$Other)] <- 0
raw_cpue$Unid.Chion[is.na(raw_cpue$Other)] <- 0
raw_cpue$Pollock[is.na(raw_cpue$Other)] <- 0
raw_cpue$Other[is.na(raw_cpue$Other)] <- 0

raw_cpue$Opilio_cons <- raw_cpue$Opilio * raw_cpue$demand_g_km2_day
raw_cpue$Bairdi_cons <- raw_cpue$Bairdi * raw_cpue$demand_g_km2_day
raw_cpue$Unid.Chion_cons <- raw_cpue$Unid.Chion * raw_cpue$demand_g_km2_day
raw_cpue$Pollock_cons <- raw_cpue$Pollock * raw_cpue$demand_g_km2_day
raw_cpue$Other_cons <- raw_cpue$Other * raw_cpue$demand_g_km2_day
raw_cpue$Missing_cons <- raw_cpue$Missing * raw_cpue$demand_g_km2_day

haul_cons <- aggregate(cbind(demand_g_km2_day, Opilio_cons, Bairdi_cons, Unid.Chion_cons, Pollock_cons, Other_cons, Missing_cons) ~ haul_join+haul_num+stratum, 
                       data = raw_cpue, sum)

strat_cons <- aggregate(cbind(demand_g_km2_day, Opilio_cons, Bairdi_cons, Unid.Chion_cons, Pollock_cons, Other_cons, Missing_cons) ~ haul_num+stratum, 
                        data = haul_cons, mean)

write.csv(raw_cpue,"raw_cpue2.csv",row.names=F)

