# R SCRIPT FOR GENERATING GUILD DATA FROM GROUNDFISH SURVEY BIOMASS
# Kerim Aydin Sept 2021

setwd("C:/src/REEM_foodhabits")

# Lookup table files
  species_lookup <- "lookups/EBS_RACE_fit_lookup_agg3_pre2021_cleanup.csv" # maps RACE codes to ECOPATH groups
  q_lookup       <- "lookups/GroupQ_pre2021.csv"        # survey Q for each PATH group, maps PATH group to guilds
  strata_lookup  <- "lookups/EBS_strata.csv"                    # Strata categories and surface areas

# RACE data directory and files
  datdir         <- "data/bts_data"  
  rawfiles <- c("ebs1982_1984.csv",        "ebs1985_1989.csv",
                "ebs1990_1994_fixed.csv",  "ebs1995_1999.csv",
                "ebs2000_2004_fixed.csv",  "ebs2005_2008.csv",
	   		        "ebs2009_2012.csv",        "ebs2013_2016.csv",
	   		        "ebs2017.csv", "ebs2018.csv", "ebs2019.csv", "nbs1982_2019.csv");  #

# Script below this line can probably can stay the same year to year.
################################################################################ 

# Read in lookup tables
  # Using RACE codes as primary key, lookup for "Ecopath" species groupings
    race_look  <- read.csv(species_lookup, row.names="RACE")
  # Using Ecopath group as primary key, lookup for catchability, data quality, and guild
    q_look     <- read.csv(q_lookup, row.names="GROUP")
  # Using RACE stratum# as primary key, lookup for "domain" definitions and surface area (km^2)
    strat_look <- read.csv(strata_lookup, row.names="Stratum")
  
# Load by-station data from list of rawfiles
  EBS_RAW  <- NULL
  for (f in rawfiles){cat(f,"\n"); flush.console(); 
    EBS_RAW <- rbind(EBS_RAW,read.csv(paste(datdir,f,sep="/")))
  }

# Apply lookup tables to raw data
  EBS_RAW$ebs_name   <- race_look[as.character(EBS_RAW$SID), "ECOPATH_Name"]
  # Detect if RACE has added any codes to the raw data that need to be put in the RACE lookup
    if (length(EBS_RAW$ebs_name[is.na(EBS_RAW$ebs_name)]) != 0){
      cat("Missing Codes Detected\n");
      unique(data.frame(EBS_RAW$SID,EBS_RAW$SCIENTIFIC,EBS_RAW$COMMON)[is.na(EBS_RAW$ebs_name),])
    }
  EBS_RAW$qq         <- q_look[EBS_RAW$ebs_name,"QQ"]
  EBS_RAW$guild      <- q_look[EBS_RAW$ebs_name,"GUILD"]
  EBS_RAW$ebs_domain <- strat_look[as.character(EBS_RAW$STRATUM),"Domain"]

# List of groups with no guild assignment
  #unique(EBS_RAW$ebs_name[is.na(EBS_RAW$guild)])

# Final selection of data to use (can apply other criteria here)
  SURV       <- EBS_RAW[!is.na(EBS_RAW$guild),]  #EBS_RAW[EBS_RAW$STRATUM<70, ] 

# Adjusted Biomass Calculation including Ecopath Q and
# multiplying by 0.1 to convert CPUE from kg/hectare to t/km^2
  SURV$qq_cpue <- 0.1 * SURV$WTCPUE * SURV$qq  
  
# List of guilds that have biomass entries in race database
  #GUILDS     <- unique(SURV$guild) #levels(as.factor(SURV$guild))[summary(SURV$guild)>0]  

# Sum survey area by domains specified in lookup
  dom_area <- tapply(strat_look$Area_km2, strat_look$Domain, sum)

# Haul summary tables - summarize by domain and year
  hauls   <- unique(SURV[,c("LATITUDE", "LONGITUDE", "STATION", "STRATUM", "YEAR", "DATETIME", 
                        "BOT_DEPTH", 
                        "BOT_TEMP", "SURF_TEMP", "VESSEL", "CRUISE", "HAUL", "ebs_domain")])

  haultots <- aggregate(STATION ~ YEAR+ebs_domain, data=hauls,length)
  area_km2 <- dom_area[haultots$ebs_domain]
  haultots <- cbind(haultots,area_km2)
  row.names(haultots) <- paste(haultots$YEAR,haultots$ebs_domain)
  
# Sum biomass by domain and add domain multipliers (#stations, area)
  qq_cpue_sum <- aggregate(qq_cpue ~ YEAR + ebs_domain + ebs_name, data=SURV,sum)
  qq_cpue_sum$guild      <- q_look[qq_cpue_sum$ebs_name,"GUILD"]  
  qq_cpue_sum$stat_count <- haultots[paste(qq_cpue_sum$YEAR,qq_cpue_sum$ebs_domain),"STATION"]
  qq_cpue_sum$area_km2   <- haultots[paste(qq_cpue_sum$YEAR,qq_cpue_sum$ebs_domain),"area_km2"]

# Final biomass in tons (can be summed across domains) and density averages 
  qq_cpue_sum$ave_bio_tkm2 <- qq_cpue_sum$qq_cpue/qq_cpue_sum$stat_count
  qq_cpue_sum$tot_bio_tons <- qq_cpue_sum$ave_bio_tkm2 * qq_cpue_sum$area_km2
  
  write.csv(qq_cpue_sum,"cpue_sum.csv",row.names=F)

  
  
  
  
# find the number of unique stations in each year (yeartot).  Can't just 
# take straight average because 0s are missing from the dataset
  #statlist  <- unique(paste(SURV$YEAR, SURV$VESSEL, SURV$CRUISE, SURV$HAUL))    
  yeartot   <- tapply(statlist,list(SURV$ebs_domain,SURV$YEAR),length)
  nyears    <- length (yeartot)

# Add Temperature and Cold Pool # NOTE THIS IS NOT USED IN GUILDS_CHAPTER
#  statkey  <-  paste(SURV$YEAR, SURV$VESSEL, SURV$CRUISE, SURV$HAUL)
#  haultemp  <- aggregate(data.frame(SURV$BOT_TEMP,SURV$YEAR,SURV$LATITUDE,SURV$LONGITUDE),list(statkey,SURV$DATETIME),mean)
#  # Nested regular expression for dates that seems to work with all date formats in data
#    hauldate <- as.Date(gsub(" .*","",(as.factor(gsub("-","/",as.character(haultemp$Group.2))))),format="%m/%d/%Y")  
#  t1 <- cut(haultemp$SURV.BOT_TEMP,c(-10000,-5,0,1,2,999),right=F)
# Tcats <- tapply(haultemp$SURV.BOT_TEMP,list(haultemp$SURV.YEAR,t1),length)
  


  G_TOTS <- matrix(0,nyears,length(GUILDS))
  rownames(G_TOTS)<-names(yeartot)
  colnames(G_TOTS)<-GUILDS   
# Multiplying by 0.1 converts CPUE from kg/hectare to t/km^2, sum by year and 
# RACE code.  Divide by # of stations in the year (yearmat) then multiply by 
# set area of the BS shelf (495,000km^2) for total tons.
  
  GOUT_LIST <- list()
  
  for (g in GUILDS){
    # g<-GUILDS[1]
    cat(g,"\n"); flush.console()
    GG <- SURV[SURV$guild==g,]
    # List of Ecopath species in the Guild
      fullG <- as.character(q_look$GROUP[q_look$GUILD==g])
    # Q value for that species (0 means don't use survey data)
      fullQ <-                 q_look$QQ[q_look$GUILD==g]
    # Ecopath biomasses to substitute for Q=0 species
     fullB <-                 q_look$BB[q_look$GUILD==g]
        
    # Crosstab sum of all groups from RACE data, 0.1 for kg/hec to t/km2
      CPUE_SUM <- tapply(0.1 * GG$qq * GG$WTCPUE,list(GG$YEAR,GG$ebs_name),sum)
        
    # Select which values are reliable from the survey, based on Q
    # (Q=0 means use the ECOPATH number   
      SurveyList <- (fullG %in% colnames(CPUE_SUM))&(fullQ>0)&!is.na(fullG)
      NotSurvey  <- !SurveyList&!is.na(fullG)
       
    # Make matrix of constant ecopath values
      UNSURVEYED <- t(matrix(rep(fullB[NotSurvey],nyears),length(fullG[NotSurvey]),nyears))
      rownames(UNSURVEYED)<-names(yeartot)
      colnames(UNSURVEYED)<-fullG[NotSurvey]
          
    # Get final values from survey data, put surveyed and unsurveyed together  
    # in a single matrix
      if (length(fullG[SurveyList])>0){ 
              CPUE_SF1  <- CPUE_SUM[,fullG[SurveyList]]
              CPUE_SF   <- CPUE_SF1[match(names(yeartot),rownames(CPUE_SF1)),]
              rownames(CPUE_SF)<-names(yeartot)
              CPUE_SF[is.na(CPUE_SF)]<-0 
              yearmat  <- matrix(yeartot,dim(CPUE_SF)[1],dim(CPUE_SF)[2])       
              CPUE_AVG <- CPUE_SF/yearmat
              GUILD_OUT <- cbind(CPUE_AVG,UNSURVEYED)
              #GUILD_OUT<-CPUE_AVG[match(names(yeartot),rownames(CPUE_AVG)),]
              #rownames(GUILD_OUT)<-names(yeartot)
              #GUILD_OUT[is.na(GUILD_OUT)]<-0
       }
       else {
              GUILD_OUT <- UNSURVEYED
       }
        
    # List by guild and species
      GOUT_LIST[[g]] <- GUILD_OUT
       
    # Convert to 1000s of tons for table of sums 
      G_TOTS[,g] <- 495000 * rowSums(GUILD_OUT)/1000
}
  
write.csv(G_TOTS,paste(outdir,"Guild_totals_2019.csv",sep="/"))

for (g in GUILDS){
    write.csv(GOUT_LIST[[g]],paste(outdir,
                              paste("Guild_2019_",sub(" ","_",g),".csv",sep=''),sep="/"))
}

