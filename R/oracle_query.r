#install.packages("DBI")
#install.packages("odbc")

library(odbc)

# Using the kaydin@afsc oracle account (password redacted)
# Note: ODBC drivers were installed by IT as part of sqlplus installation 
# done May 2021 - some other driver versions (on others' computers) didn't
# work with the library.  
# Important: You must be logged into VPN for this to work.

con <- dbConnect(odbc(), "AFSC", uid="kaydin", pwd="[REDACTED - insert PW]")

setwd("C:/e/src/queries")

regions <- c("BS","GOA","AI")
RAWDAT  <- list()

for (RR in regions){
  qtest <- paste(
    "SELECT",
    "NPRED.nodc as PRED_NODC, NPRED.name as PRED_NAME,",
	"HH.REGION, HH.VESSEL, HH.CRUISE, HH.HAUL, HH.STRATUM, HH.YEAR, HH.MONTH, HH.DAY,",
	"HH.GEAR_TEMP, HH.GEAR_DEPTH, HH.BOTTOM_DEPTH, HH.SURFACE_TEMP, HH.START_HOUR,",
	"HH.HAULJOIN, HH.STATIONID, HH.CRUISE_TYPE,HH.RLAT, HH.RLONG,",
	"PP.PREDJOIN, PP.PRED_LEN,",
	"NPREY.nodc as PREY_NODC, NPREY.name as PREY_NAME,", #"NPREY.ECOPATH_Prey,", 
	"sum(PP.PREY_TWT) as TWT, sum(PP.PREY_CNT) as CNT", 
	"FROM",
	"foodlab.predprey PP, foodlab.haul HH, foodlab.nodc NPRED, foodlab.nodc NPREY",
	"WHERE",
	"(PP.VESSEL = HH.VESSEL AND PP.CRUISE = HH.CRUISE AND PP.HAUL = HH.HAUL)", 
	"AND (PP.PRED_NODC = NPRED.nodc AND PP.PREY_NODC = NPREY.nodc)", 
	"AND (HH.REGION=", paste("'", RR ,"'), ",sep=''),
	"AND (HH.CRUISE_TYPE='Race_Groundfish')",
	"AND (HH.YEAR <= 2015)",
	"AND (HH.MONTH >= 5 AND HH.MONTH <= 8)",
	"GROUP BY",
	"PP.PREDJOIN, NPRED.nodc, NPRED.name, HH.REGION,",
	"HH.VESSEL, HH.CRUISE, HH.HAUL, HH.STRATUM, HH.YEAR, HH.MONTH, HH.DAY,",
	"HH.GEAR_TEMP, HH.GEAR_DEPTH, HH.BOTTOM_DEPTH, HH.SURFACE_TEMP,",
	"HH.START_HOUR, HH.HAULJOIN, HH.STATIONID, HH.CRUISE_TYPE, HH.RLAT,",
	"HH.RLONG, PP.PRED_LEN,",
	"NPREY.nodc, NPREY.name", #"NPREY.ECOPATH_Prey",
    "",
  sep=' ')

  RAWDAT[[RR]] <- dbGetQuery(con, qtest)
  write.csv(RAWDAT[[RR]], paste(RR,"_raw.csv",sep=''), row.names=F)
}


nodc <- dbGetQuery(con, "select * from foodlab.nodc")

dbDisconnect(con)

