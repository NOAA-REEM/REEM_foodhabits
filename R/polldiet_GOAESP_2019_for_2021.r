
setwd("C:/src/REEM_foodhabits")

require(VGAM)
   
#rawdat   <- rbind(read.csv("GOA_allPrey_out.csv"),read.csv("BS_allPrey_out.csv"))
rawdat   <- read.csv("data/GOA_raw.csv.gz")
preylook <- read.csv("lookups/GOA_med_diet_lookup_2019.csv")
preycol  <- "GOA_MED"
#preds <- list("W.Pollock" = "W. Pollock",
#              "P.Cod" = "P. Cod",
#              "Arrowtooth" = c("Arrow or Kam","Arrowtooth","Kamchat fl"),
#              "P.Halibut"  = "P. Halibut"
#               )

preds   <- list("W.Pollock" = "Theragra chalcogramma (walleye pollock)")
LCLASS  <- c(0,10,25,40,60,999) 
lenlist <- levels(cut(seq(0,999),LCLASS,right=F))
#lenlist <- c("[0,10)", "[10,25)", "[25,40)", "[40,60)", "[60,999)")
#preds <- list("P.Halibut"  = "P. Halibut") 
##LCLASS <- c(0,20,40,60,80,100,999)
#LCLASS <- c(0,10,20,30,40,50,60,70,80,90,100,999)
          
# Halibut mm to cm regression for a (6.291e-6)*10^3.24 =  0.01093251
A_L <- c("W.Pollock"=0.00553096, "P.Cod"=0.00411781, "Arrowtooth"=0.00443866,"P.Halibut"=0.01093251)
B_L <- c("W.Pollock"=3.044172,   "P.Cod"=3.25325765, "Arrowtooth"=3.19894001,"P.Halibut"=3.24)

#lenlist <-  c(0,10,20,30,40,50,60,70,80,90,100)
# Build a crosstab-query (per-predator) using names in GOA_MED column
  guild <- as.character(preylook[,preycol])
  names(guild) <- as.character(preylook$PREYNAME)
  glist <- guild[as.character(rawdat$ECOPATH_PREY)]     
  gdat <- aggregate(rawdat$TWT,list(as.character(rawdat$PREDJOIN),glist),sum)
  gtab <- tapply(rawdat$TWT,list(as.character(rawdat$PREDJOIN),glist),sum)
  gtab[is.na(gtab)]<-0
  predtab <- unique(rawdat[,1:23])  
  preycross <- gtab[as.character(predtab$PREDJOIN),]
  allprey   <- colnames(preycross)
  lbin  <- cut(predtab$PRED_LEN,LCLASS,right=F)
  #lbin <- ifelse(predtab$PRED_LEN>=100,100,10*floor(predtab$PRED_LEN/10))
  GDAT<-cbind(predtab,lbin,preycross)


yrs <- c(1990, 1993, 1996, 1999, 2001, 2003, 2005, 2007, 2009, 2011, 2013, 2017, 2019)
#DOMS <- list("cent_west" = c("Centralshelf","Centralgully","Centralslope",
#                             "Westshelf","Westgully","Westslope"))
DOMS <- list("cent_west" = c(10,11,12,13,20,21,22,30,31,32,33,35,
                     110,111,112,120,121,122,130,131,132,133,134,
            210,220,221,230,231,232,310,320,330,410,420,430,510,520,530))

#DOMS <- list("central" = c("Centralshelf","Centralgully","Centralslope"),
#             "west"    = c("Westshelf","Westgully","Westslope"),
#             "east"    = c("Eastshelf","Eastgully","Eastslope"))

#yrs <- 1985:2014
#DOMS <- list("Inner"=c("InnerSE","InnerNW"),
#             "Middle"=c("MiddleSE","MiddleNW"),
#             "Outer"=c("OuterSE","OuterNW")
#            )

DETECT <- 0.0001 # Grams of empty prey for proportioning
SAMPLES <- 10000; 

#FinalEst  <- NULL
Nests     <- length(yrs)*length(preds)*length(lenlist)*length(DOMS)
alphalist <- matrix(NA,Nests,length(allprey))
colnames(alphalist)<-allprey
Nstomachs <- Nresamp <- Wdetect <- twtMean  <- twtSD <- cperwMean <- cperwSD  <- rep(NA,Nests) 
ThisSp <- ThisDom <- ThisLen <- ThisYr <- rep(NA,Nests) 

YIND <- 0
for (SPN in 1:(length(preds))){
  sp <- names(preds)[SPN] 
  spcodes <- preds[[sp]]
        TotWt0 <- as.numeric(rowSums(GDAT[,allprey]))
        SURV   <- GDAT[GDAT$CRUISE_TYPE=="Race_Groundfish" & 
                  GDAT$MONTH%in%6:8 & 
                  GDAT$PRED_NAME%in%spcodes &
                  TotWt0>0.0, ] 
  for (yr in yrs){
    for (dn in 1:(length(DOMS))){
      dom <- DOMS[[dn]]
      domname <- names(DOMS)[dn]
    #for (dom in DOMS){
      for (sbin in lenlist){
        YIND <- YIND + 1
        cat(sp,yr,domname,sbin,YIND,"of",Nests,"\n"); flush.console() 
        SELPRED <- SURV[SURV$lbin%in%sbin &
                        SURV$STRATUM%in%dom & 
                        SURV$YEAR%in%yr, ]        
        Nind <- length(SELPRED[,1])
        if (Nind>=10){
            cat(Nind,"\n"); flush.console()
            allfood <- SELPRED[,allprey]
            goodprey <- allfood[,colSums(allfood)>0]
            #Breaks if only one prey type is >0 (vector conversion); TODO fix
            sptot    <- colSums(goodprey)
            spadd    <- DETECT*sptot/sum(sptot)
            SCI <- t(t(goodprey)+spadd)/ (A_L[sp] * (SELPRED$PRED_LEN ^ B_L[sp]))
            sp_prey  <- colnames(goodprey)
            
            # Make matrix of individual stomachs
              sampmat <- matrix(as.numeric(unlist(SCI)),length(SCI[,1]),length(sp_prey))
              colnames(sampmat)<-sp_prey
              #Nind <- length(sampmat[,1])
      
           # Resample sums 
             resampmat <- matrix(NA,SAMPLES,length(sp_prey)) 
             colnames(resampmat) <- sp_prey
             tot_wt <- cperw <- rep(NA,SAMPLES)
             for (i in 1:SAMPLES){
               IND <- sample.int(Nind, Nind, replace = T)
               tot_wt[i]    <- sum(rowSums(sampmat[IND,]))
               cperw[i]     <- sum(rowSums(sampmat[IND,]))/Nind
               pathsum <- colSums(sampmat[IND,]) #+ spadd
               resampmat[i,]<-pathsum
             }
      
           # Fit resampled matrix to dirichlet
             sampprop <- resampmat/rowSums(resampmat)  
             tfit <- vglm(sampprop~1,dirichlet)
             alphas <- Coef(tfit)

# Final accumulation for loop
  alphalist[YIND,sp_prey] <- alphas
  Nstomachs[YIND]  <- Nind
  Nresamp[YIND]    <- SAMPLES
  Wdetect[YIND]    <- DETECT
  twtMean[YIND]    <- mean(tot_wt) 
  cperwMean[YIND]  <- mean(cperw)
  twtSD[YIND]      <- sd(tot_wt) 
  cperwSD[YIND]    <- sd(cperw)
 }
 else{
  #alphalist[YIND,sp_prey] <- alphas
  Nstomachs[YIND]  <- length(SELPRED[,1])
  Nstomachs[YIND]  <- Nind
  Nresamp[YIND]    <- NA
  Wdetect[YIND]    <- NA
  twtMean[YIND]    <- NA 
  cperwMean[YIND]  <- NA
  twtSD[YIND]      <- NA 
  cperwSD[YIND]    <- NA    
 }
 ThisSp[YIND] <- sp
 ThisYr[YIND] <- yr
 ThisDom[YIND] <- domname
 ThisLen[YIND] <- sbin
}}}}

alphalist[is.na(alphalist)]<-0
alpha0 <- rowSums(alphalist)
OutEst <- data.frame(ThisSp,ThisYr,ThisDom,ThisLen,Nstomachs,Nresamp,Wdetect,twtMean,twtSD,cperwMean,cperwSD,alpha0,alphalist)

CleanEst <- OutEst[!is.na(OutEst$Nresamp),]

dcol <- "Euphausiids"

dat <- CleanEst
   aa  <- dat[,dcol]
   az  <- dat$alpha0
   mn <- aa/az
   up <- qbeta(0.975,aa,az-aa)
   dn <- qbeta(0.025,aa,az-aa)
  
res <- data.frame(dat$ThisSp,dat$ThisYr,dat$ThisDom,dat$ThisLen,dat$Nstomachs,mn,up,dn)
colnames(res) <- c("Predator","Year","Domains","PredLength","Samples","Euphausiid_prop","upper95CI","lower95CI")  

#write.csv(OutEst,fout,row.names=F)
} # end LCLASS loop

} # end SPECIES loop



GOAfinal <- FinalEst
EBSfinal <- FinalEst
allfinal <- rbind(GOAfinal,EBSfinal)


cleanfin <- allfinal[!is.na(allfinal$Nresamp),]

pdoms <- c("central","west","Inner","Middle","Outer")
pnames <- c("Central GOA", "West GOA", "Inner BS", "Middle BS", "Outer BS")
psym  <- c(15,16,18,19,4)
pcol  <- c("blue","cyan","red","purple","black")
dset <- c( "Mesozooplankton", "Motile.epifauna","Pelagic.foragers","Benthic.foragers")

#"Infauna", ,"Other" "Structural.epifauna",
par(mfrow=c(2,2))
for (dcol in dset){
plot(seq(0,100,10),seq(0,1,0.1),ylim=c(0,1),type='n',ylab="Proportion in diet",xlab="Halibut fork length")
title(dcol)
if (dcol=="Mesozooplankton"){legend(75,1,pnames,pch=psym,col=pcol,bty="n")}
#dcol<-"Motile.epifauna"
k<-0
for (d in 1:(length(pdoms))){
   dd <- pdoms[d]
   dat <- cleanfin[cleanfin$domain==dd,]
   xx  <- dat$Lbin
   aa  <- dat[,dcol]
   az  <- dat$alpha0
   mn <- aa/az
   up <- qbeta(0.975,aa,az-aa)
   dn <- qbeta(0.025,aa,az-aa)
   segments(xx+k,y0=up,y1=dn,col=pcol[d])
   points(xx+k,mn,ylim=c(0,1),pch=psym[d],col=pcol[d])
   #text(xx+k,mn,k,ylim=c(0,1))
   k<-k+1 
}
}

dat <- FinalEst$Motile.epifauna
err <- FinalEst$alpha0

pdom <- "central"
aa   <- dat[FinalEst$domain==pdom]
az   <- err[FinalEst$domain==pdom]
mn <- aa/az
up <- qbeta(0.975,aa,az-aa)
dn <- qbeta(0.025,aa,az-aa)
segments(seq(0,100,10),y0=up,y1=dn,col="red")
points(seq(0,100,10),mn,ylim=c(0,1),pch=15,col="red")
pdom <- "west"
aa   <- dat[FinalEst$domain==pdom]
az   <- err[FinalEst$domain==pdom]
mn <- aa/az
up <- qbeta(0.975,aa,az-aa)
dn <- qbeta(0.025,aa,az-aa)
segments(seq(0,100,10)+1,y0=up,y1=dn,col="blue")
points(seq(0,100,10)+1,mn,ylim=c(0,1),pch=15,col="blue")



write.csv(FinalEst,"IPHC_finalest.csv",row.names=F)

pmean <-  alphalist/alpha0 
loCI  <- qbeta(0.025,alphalist,alpha0 - alphalist)
hiCI  <- qbeta(0.975,alphalist,alpha0 - alphalist)

xx <-as.numeric(rownames(pmean))
py <-"Pandalidae" 
plot(xx,pmean[,py],ylim=c(0,max(hiCI[,py])))
segments(xx,y0=lowCI[,py],y1=hiCI[,py])


  names(alphas)<- sp_prey
  alpha0 <- sum(alphas)

  xx <- seq(0,1,0.001)
  plot(xx,dbeta(xx,alphas["Euphausiids"],alpha0))
  
  
#TotWt <- as.numeric(rowSums(SURV[,allprey]))






#PREDS <- list("AK Plaice", 
#           c("AK Skate", "Unid Bathyraja","Unid Rajidae"), 
#           c("Aleutian Skate","Bering Skate", "Big Skate","Black Skate","Commander Skate","WhtBlotch Skate","Mud Skate"), 
#           c("Arrow or Kam","Arrowtooth","Kamchat fl"),  
#           "Capelin",   
#           "Eelpout", 
#           "Eulachon", 
#           "FH Sole", 
#           c("Gen Rock Sole","N Rock Sole","S Rock Sole"),  
#           c("Giant Grenadier","Macrouridae","Pacific Grenadier"), 
#           "Gr. Turbot", 
#           "Greenlings", 
#           "Herring", "Lg Sculpin", 
#           "Managed Forage", "Misc. Flatfish", 
#           "Oth pel. smelt", "P. Cod", "P. Halibut",  
#           "POP", "Rex sole",  "Sablefish",  "Sandlance", 
#           "Sculpin","W. Pollock",  "YF Sole")