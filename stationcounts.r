# Data files (bottom-trawl public format, from RACE)
  datdir   <- "data/bts_data"
  rawfiles <- c("ebs1982_1984.csv",        "ebs1985_1989.csv",
                "ebs1990_1994_fixed.csv",  "ebs1995_1999.csv",
                "ebs2000_2004_fixed.csv",  "ebs2005_2008.csv",
	   		    "ebs2009_2012.csv",        "ebs2013_2016.csv",
	   		    "ebs2017.csv", "ebs2018.csv", "ebs2019.csv");
    
# Load by-station data
  EBS_RAW  <- NULL		
  for (f in rawfiles){cat(f,"\n"); flush.console(); 
    EBS_RAW <- rbind(EBS_RAW,read.csv(paste(datdir,f,sep="/"),stringsAsFactors=F))
  }

hauls <- unique(EBS_RAW[,c("LATITUDE","LONGITUDE","STATION","STRATUM","YEAR","DATETIME",
                           "BOT_DEPTH","BOT_TEMP","SURF_TEMP","VESSEL","CRUISE","HAUL")])




