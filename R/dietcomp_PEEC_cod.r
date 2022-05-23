output_file  <- "results/out_diets_2021.csv"
ppfile       <- "data/BS_raw.csv.gz"                    # raw Pred/prey diet file (RACEBASE query result)
pplookfile   <- "lookups/Alaska_PreyLookup_MASTER.csv"  # lookup file for grouping prey
stratlookfile<- "lookups/EBS_strata.csv" 
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
predators <- c("P.cod") #,"P.halibut")
  
# Years to output
yearlist <- 1985:2021

# Get Stratum data and strata combining information for making larger domains (bigstrats)
  strat_table <- read.csv(stratlookfile, row.names="Stratum")
  bigstrat <- aggregate(Area_km2~Domain,strat_table,sum)
  bigstrat_area<-bigstrat$Area_km2; names(bigstrat_area)<-bigstrat$Domain
# Results should be strat_table, bigstrat, and bigstrat_area
  
# Get Haul data from all stations (including those with no stomachs)
# formerly in source("stationcounts.r") 
  bts_datdir   <- "data/bts_data" # Data files (bottom-trawl public format, from RACE)
  rawfiles <- c("ebs1982_1984.csv",        "ebs1985_1989.csv",
              "ebs1990_1994_fixed.csv",  "ebs1995_1999.csv",
              "ebs2000_2004_fixed.csv",  "ebs2005_2008.csv",
              "ebs2009_2012.csv",        "ebs2013_2016.csv",
              "ebs2017.csv", "ebs2018.csv", "ebs2019.csv","ebs2021.csv",
              "nbs1982_2019.csv","nbs2021.csv");
  EBS_RAW  <- NULL		
  for (f in rawfiles){cat(f,"\n"); flush.console(); 
    EBS_RAW <- rbind(EBS_RAW,read.csv(paste(bts_datdir,f,sep="/"),stringsAsFactors=F))
  }
  hauls <- unique(EBS_RAW[,c("LATITUDE","LONGITUDE","STATION","STRATUM","YEAR","DATETIME",
                             "BOT_DEPTH","BOT_TEMP","SURF_TEMP","VESSEL","CRUISE","HAUL")])

  hauls$bigstrat <- strat_table[as.character(hauls$STRATUM),"Domain"]
  haulcounts <- aggregate(STATION ~ bigstrat+YEAR, hauls, length)
  row.names(haulcounts)<-paste(haulcounts$bigstrat,haulcounts$YEAR,sep='_')
# Results should be hauls table and haulcounts table 

  
# Read in lookup and raw data and ensure 10-digit nodc numbers are read as text keys
  preylooktable <- read.csv(pplookfile)
  preylook <- preylooktable[,preylook_col]
  names(preylook)<-sprintf("%010.0f",preylooktable$NODC_CODE)
  fullprey<-unique(preylooktable[,preylook_col])
  rawdat   <- read.csv(ppfile)
  rawdat$PREY_NODC <-sprintf("%010.0f",rawdat$PREY_NODC)
  rawdat$PRED_NODC <-sprintf("%010.0f",rawdat$PRED_NODC)  
  rawdat$preyguild <- preylook[rawdat$PREY_NODC]
# Result should be rawdat, containing all stomach data read in with NODC name lookups applied
    
  
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
  
  
  for (STRAT in bigstrat$Domain){  # Was stratblock
    cat(PRED,STRAT,"\n"); flush.console()
    stratcodes <- row.names(strat_table)[which(strat_table$Domain==STRAT)]  #unlist(stratblock[STRAT])
    
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
          # THe below calculate the empirical logit transform for presence/absence
          # (frequency of occurrence) data.
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
        
write.csv(o_vals,output_file,row.names=F)

##




raw_original <- read.csv("data/EBS_Pcod_length_cpue.csv")
raw_original <- rbind(raw_original,read.csv("data/NBS_Pcod_length_cpue.csv"))

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
raw_cpue$wt_g <- ThisPred$A_L * ((raw_cpue$fish_length/10)^ThisPred$B_L)
raw_cpue$Cmax_g_day_fish <- bioen$CA * (raw_cpue$wt_g^(1+bioen$CB))

# need to replace gear_temp NAs with stratum averages
Tfix <- raw_cpue[!is.na(raw_cpue$gear_temp),c("haul_num","stratum","gear_temp")]
tdat <- aggregate(Tfix$gear_temp,list(Tfix$haul_num,Tfix$stratum),mean)
tave <- tdat$x; names(tave) <- paste(tdat$Group.1,tdat$Group.2)
raw_cpue$tclean <- ifelse(is.na(raw_cpue$gear_temp),tave[paste(raw_cpue$haul_num,raw_cpue$stratum)],raw_cpue$gear_temp)

tmp<-list(TempC=raw_cpue$tclean)
raw_cpue$fTC<-fTC_fun(bioen,tmp)

raw_cpue$numcpue_at_len   <- raw_cpue$num_cpue * (raw_cpue$frequency/raw_cpue$tot_freq)
raw_cpue$NtimesW          <- raw_cpue$numcpue_at_len * raw_cpue$wt_g
raw_cpue$Cmax_Ncpue       <- raw_cpue$numcpue_at_len * raw_cpue$Cmax_g_day_fish
raw_cpue$demand_g_km2_day <- raw_cpue$Cmax_Ncpue * raw_cpue$fTC 


#Diet table lookups
stratlook <- sub("[0-9]","",names(unlist(stratblock)))
names(stratlook) <- unlist(stratblock)
raw_cpue$bigstrat <- stratlook[as.character(raw_cpue$stratum)]

raw_cpue$lbin <- cut(raw_cpue$fish_length/10,ThisPred$LCLASS,right=F)

dietkey <- paste(PRED,raw_cpue$bigstrat, raw_cpue$haul_num, raw_cpue$lbin)

prey <- o_vals[o_vals$sp_prey=="Opilio",]
row.names(prey) <- paste(prey$PRED,prey$STRAT, prey$YY,prey$LL)
raw_cpue$SCIperN <- prey[dietkey,"SCIperN"]
raw_cpue$Opilio_cperfull <- prey[dietkey,"cperfull"]
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

prey <- o_vals[o_vals$sp_prey=="Opilio",]
write.csv(prey,"results/Opilio_dc_all.csv",row.names=F)
write.csv(raw_cpue,"results/haul_cpue_cons_all.csv",row.names=F)


strat_len_cons <- aggregate(cbind(numcpue_at_len, NtimesW, Cmax_Ncpue, demand_g_km2_day) ~ bigstrat+haul_num+lbin, 
                           data = raw_cpue, sum)
strat_len_cons$scount <- haulcounts[paste(strat_len_cons$bigstrat,strat_len_cons$haul_num,sep='_'),"STATION"]
strat_len_cons$sarea <- bigstrat_area[strat_len_cons$bigstrat]

strat_len_cons[,c("numcpue_millions_extrap","NtimesW_mt_extrap","Cmax_mt_day_extrap","demand_mt_day_extrap")] <- 
  strat_len_cons[,c("numcpue_at_len","NtimesW","Cmax_Ncpue","demand_g_km2_day")]/strat_len_cons$scount * strat_len_cons$sarea/1e6 # from g to mt
names(strat_len_cons)[2]<-"year"

dietkey <- paste(PRED,strat_len_cons$bigstrat, strat_len_cons$year, strat_len_cons$lbin)

  prey <- o_vals[o_vals$sp_prey=="Opilio",]
  row.names(prey) <- paste(prey$PRED,prey$STRAT, prey$YY,prey$LL)
  strat_len_cons$SCIperN <- prey[dietkey,"SCIperN"]
  strat_len_cons$Opilio_cperfull <- prey[dietkey,"cperfull"]
  strat_len_cons$Opilio <- prey[dietkey,"dietprop"]

  prey <- o_vals[o_vals$sp_prey=="Bairdi",]
  row.names(prey) <- paste(prey$PRED,prey$STRAT, prey$YY,prey$LL)
  strat_len_cons$Bairdi <- prey[dietkey,"dietprop"]

  prey <- o_vals[o_vals$sp_prey=="Unid.Chion",]
  row.names(prey) <- paste(prey$PRED,prey$STRAT, prey$YY,prey$LL)
  strat_len_cons$Unid.Chion <- prey[dietkey,"dietprop"]

  prey <- o_vals[o_vals$sp_prey=="W.pollock",]
  row.names(prey) <- paste(prey$PRED,prey$STRAT, prey$YY,prey$LL)
  strat_len_cons$Pollock <- prey[dietkey,"dietprop"]

  strat_len_cons$Other <- 1-strat_len_cons$Opilio-strat_len_cons$Bairdi-strat_len_cons$Unid.Chion-strat_len_cons$Pollock

  strat_len_cons$Missing <- ifelse(is.na(strat_len_cons$Other),1,0)

  strat_len_cons$Opilio[is.na(strat_len_cons$Other)] <- 0
  strat_len_cons$Bairdi[is.na(strat_len_cons$Other)] <- 0
  strat_len_cons$Unid.Chion[is.na(strat_len_cons$Other)] <- 0
  strat_len_cons$Pollock[is.na(strat_len_cons$Other)] <- 0
  strat_len_cons$Other[is.na(strat_len_cons$Other)] <- 0


  strat_len_cons$Opilio_cons <- strat_len_cons$Opilio * strat_len_cons$demand_mt_day_extrap
  strat_len_cons$Bairdi_cons <- strat_len_cons$Bairdi * strat_len_cons$demand_mt_day_extrap
  strat_len_cons$Unid.Chion_cons <- strat_len_cons$Unid.Chion * strat_len_cons$demand_mt_day_extrap
  strat_len_cons$Pollock_cons <- strat_len_cons$Pollock * strat_len_cons$demand_mt_day_extrap
  strat_len_cons$Other_cons <- strat_len_cons$Other * strat_len_cons$demand_mt_day_extrap
  strat_len_cons$Missing_cons <- strat_len_cons$Missing * strat_len_cons$demand_mt_day_extrap
  
#strat_len_cons$num_avg <- strat_len_cons$numcpue_at_len/strat_len_cons$scount
#strat_len_cons$num_extrap <- strat_len_cons$num_avg * strat_len_cons$sarea

write.csv(strat_len_cons,"results/strat_len_cons_all.csv",row.names=F)


#################################################################

haul_cons <- aggregate(cbind(demand_g_km2_day, Opilio_cons, Bairdi_cons, Unid.Chion_cons, Pollock_cons, Other_cons, Missing_cons) ~ haul_join+haul_num+stratum, 
                       data = raw_cpue, sum)



strat_cons <- aggregate(cbind(demand_g_km2_day, Opilio_cons, Bairdi_cons, Unid.Chion_cons, Pollock_cons, Other_cons, Missing_cons) ~ haul_num+stratum, 
                        data = haul_cons, mean)

strat_tons <- strat_cons[,3:9] * stratarea_km[as.character(strat_cons$stratum)]/1e6
strat_out <- cbind(strat_cons[,1:2],strat_tons)
names(strat_out) <- c("Year","Stratum","Demand_mt_day","Opilio_mt","Bairdi_mt","Unid.Chion_mt","Pollock_mt","Other_mt","Missing_mt")



write.csv(strat_out,"Pcod_consumption.csv",row.names=F)


## LENGTH STRATA STUFF







