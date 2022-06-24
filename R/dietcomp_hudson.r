
# output file for final diet table 
  output_file  <- "results/out_diets_AI.csv"

# raw data file (in gz compressed format)
  ppfile       <- "data/BS_raw.csv.gz"                    # raw Pred/prey diet file (RACEBASE query result)

# Lookup for prey codes
  pplookfile   <- "lookups/Alaska_PreyLookup_MASTER.csv"  # lookup file for grouping prey

# Column to use
  preylook_col <- "EBS_CRAB"                      # Column name in pplookfile to use to group prey

# Lookup for stratum groupings
  stratlookfile<- "lookups/EBS_strata.csv"

# Minimum sample size - strata/years without at least this many samples are not reported 
  min_sample   <- 5                               

# List of predator-specific values.
# nodc is the list of nodc codes (or sole nodc code)
# A and B are length weight regression params (cm to grams)
# Lclass is length class divisions to use (in cm)
#
preds <- list(
  "W.pollock"     = list(nodc="8791030701", A_L=0.00553096, B_L=3.044172,   LCLASS=c(0,10,25,40,55,999) ),
  "P.cod"         = list(nodc="8791030401", A_L=0.00411781, B_L=3.25325765, LCLASS=c(0,10,30,60,85,999) ),
  "Arrowtooth"    = list(nodc="8857040102", A_L=0.00443866, B_L=3.19894001, LCLASS=c(0,10,30,50,999)    ),
  "P.halibut"     = list(nodc="8857041901", A_L=0.01093251, B_L=3.24,       LCLASS=c(0,10,50,70,999)    ),
  "POP"           = list(nodc="8826010102", A_L=0.0078,     B_L=3.15,       LCLASS=c(0,20,999)          ),
  "Atka.mackerel" = list(nodc="8827010501", A_L=0.0088,     B_L=3.273,      LCLASS=c(0,20,999)          )
)

# Predators to calculate (all must be on above list)
predators <- c("P.cod") #,"P.halibut")
  
# Years to output
yearlist <- 1985:2021  # can also use list yearlist <- c(1985,1987,1989) etc


# Leaving Biomass queries (from RACE Biomass data) out for now
# (this is written for EBS)
# # Get Stratum data and strata combining information for making larger domains (bigstrats)
#   strat_table <- read.csv(stratlookfile, row.names="Stratum")
#   bigstrat <- aggregate(Area_km2~Domain,strat_table,sum)
#   bigstrat_area<-bigstrat$Area_km2; names(bigstrat_area)<-bigstrat$Domain
# # Results should be strat_table, bigstrat, and bigstrat_area
#   
# # Get Haul data from all stations (including those with no stomachs)
# # formerly in source("stationcounts.r") 
#   bts_datdir   <- "data/bts_data" # Data files (bottom-trawl public format, from RACE)
#   rawfiles <- c("ebs1982_1984.csv",        "ebs1985_1989.csv",
#               "ebs1990_1994_fixed.csv",  "ebs1995_1999.csv",
#               "ebs2000_2004_fixed.csv",  "ebs2005_2008.csv",
#               "ebs2009_2012.csv",        "ebs2013_2016.csv",
#               "ebs2017.csv", "ebs2018.csv", "ebs2019.csv","ebs2021.csv",
#               "nbs1982_2019.csv","nbs2021.csv");
#   EBS_RAW  <- NULL		
#   for (f in rawfiles){cat(f,"\n"); flush.console(); 
#     EBS_RAW <- rbind(EBS_RAW,read.csv(paste(bts_datdir,f,sep="/"),stringsAsFactors=F))
#   }
#   hauls <- unique(EBS_RAW[,c("LATITUDE","LONGITUDE","STATION","STRATUM","YEAR","DATETIME",
#                              "BOT_DEPTH","BOT_TEMP","SURF_TEMP","VESSEL","CRUISE","HAUL")])
# 
#   hauls$bigstrat <- strat_table[as.character(hauls$STRATUM),"Domain"]
#   haulcounts <- aggregate(STATION ~ bigstrat+YEAR, hauls, length)
#   row.names(haulcounts)<-paste(haulcounts$bigstrat,haulcounts$YEAR,sep='_')
# # Results should be hauls table and haulcounts table 

  
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


