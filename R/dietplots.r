library("ggplot2")
setwd("c:/src/REEM_foodhabits")

dat <- read.csv("results/out_diets_5cm.csv")

sdat <- dat[!(dat$LL %in% c("[5,10)","[100,999)")),c("sp_prey","LL","dietprop","STRAT")]
sdat$LL <- as.factor(sdat$LL)
sdat$sp_prey <- as.factor(sdat$sp_prey)

for(DOM in bigstrat$Domain){

sdat_filter <- sdat[sdat$STRAT==DOM,]  

agg <- aggregate(dietprop~LL+sp_prey,data=sdat_filter,mean)


sp_order <- c("Amphipods","Euphausiid","Shrimp","Polychaete","Epifauna",
              "Forage.fish","Opilio","Bairdi","Unid.Chion","King.crab","Octopus",
              "Flatfish","W.pollock")

not_shown <- c("Copepods","Oth.nonfish","Rockfish","Discards","Infauna")

png("core_all_lines.png",width=1500,height=1000)
#png(paste("core_",DOM,".png",sep=''),width=1500,height=1000)
yscale <- 5.5; datscale <-5
plot(NULL,xlim=c(-1,nlevels(agg$LL)+1), ylim=c(0,yscale*length(sp_order)+datscale),axes=F,xlab='',ylab='')

for (ss in 1:length(sp_order)){
  d <- agg[agg$sp_prey==sp_order[ss],]
  polygon(c(as.numeric(d$LL),rev(as.numeric(d$LL))),
          c(yscale*ss+d$dietprop*datscale,rev(yscale*ss-d$dietprop*datscale)), 
          col="forestgreen",border="forestgreen")
  text(0.9,yscale*ss,sp_order[ss],pos=2,cex=1.8)
  text(1:nlevels(agg$LL),0,substr(levels(agg$LL),2,3),cex=1.8)
  lines(c(nlevels(agg$LL)+0.5,nlevels(agg$LL)+0.5),c(7*yscale+datscale,7*yscale-datscale))
  lines(c(nlevels(agg$LL)+0.4,nlevels(agg$LL)+0.6),c(7*yscale+datscale,7*yscale+datscale))
  lines(c(nlevels(agg$LL)+0.4,nlevels(agg$LL)+0.6),c(7*yscale-datscale,7*yscale-datscale))
    text(nlevels(agg$LL)+0.5,7*yscale,"100%",pos=4,cex=1.8)
   mtext("Pacific cod fork length (cm)",side=1,cex=1.8)
   mtext("%diet by weight",side=2,cex=2)
   mtext("Pacific cod - summer BTS",side=3,line=0,cex=2,adj=0.1)
   #mtext(DOM,side=3,line=-2,cex=2.0,adj=0.18)   
}
abline(v=5)
abline(v=11)
abline(v=16,lty=2,col="grey40")
dev.off()
}


#####################################################################################
dat <- read.csv("results/out_diets_cod.csv")

#sdat <- dat[!(dat$LL %in% c("[5,10)","[100,999)")),c("sp_prey","LL","dietprop","STRAT")]

# "[30,60)"  "[60,85)"  "[10,30)"  "[85,999)" "[0,10)"

LEN <- "[60,85)"; lname <- "piscivores"
sdat <- dat[dat$LL %in% LEN , c("sp_prey","YY","LL","dietprop","STRAT")]
sdat$YY <- as.factor(sdat$YY)
sdat$sp_prey <- as.factor(sdat$sp_prey)

for(DOM in bigstrat$Domain){
  
  sdat_filter <- sdat[sdat$STRAT==DOM,]  
  
  if (length(sdat_filter[,1])){
  agg <- aggregate(dietprop~YY+sp_prey,data=sdat_filter,mean)
  
  
  sp_order <- c("Amphipods","Euphausiid","Shrimp","Polychaete","Epifauna",
                "Forage.fish","Opilio","Bairdi","Unid.Chion","King.crab","Octopus",
                "Flatfish","W.pollock")
  
  not_shown <- c("Copepods","Oth.nonfish","Rockfish","Discards","Infauna")
  
  #png("core_all_lines.png",width=1500,height=1000)
  png(paste("core",lname,DOM,".png",sep=''),width=1500,height=1000)
  yscale <- 5.5; datscale <-5
  plot(NULL,xlim=c(-2,length(1985:2019)+2), ylim=c(0,yscale*length(sp_order)+datscale),axes=F,xlab='',ylab='')
  
  for (ss in 1:length(sp_order)){
    d <- agg[agg$sp_prey==sp_order[ss],]
    
    for (yy in 1985:2019){ 
      t <- which(yy==d$YY)
      if (length(t)){
        x1 <- yy-1984
        polygon(c( x1-0.5, x1-0.5, x1+0.5, x1+0.5),
             c(yscale*ss-d$dietprop[t]*datscale, yscale*ss+d$dietprop[t]*datscale,
               yscale*ss+d$dietprop[t]*datscale, yscale*ss-d$dietprop[t]*datscale), 
             col="darkred",border="darkred")
      }
    }
    text(0,yscale*ss,sp_order[ss],pos=2,cex=1.8)
    text(1:length(1985:2019),0,1985:2019,cex=1.8,srt=90)
    lines(c(length(1985:2019)+1.0,length(1985:2019)+1.0),c(7*yscale+datscale,7*yscale-datscale))
    lines(c(length(1985:2019)+1.1,length(1985:2019)+0.9),c(7*yscale+datscale,7*yscale+datscale))
    lines(c(length(1985:2019)+1.1,length(1985:2019)+0.9),c(7*yscale-datscale,7*yscale-datscale))
    text(length(1985:2019)+1.0,7*yscale,"100%",pos=4,cex=1.8)
    #mtext("Pacific cod fork length (cm)",side=1,cex=1.8)
    mtext("%diet by weight",side=2,cex=2,line=1.1)
    #mtext("Pacific cod - summer BTS",side=3,line=0,cex=2,adj=0.1)
    mtext(paste(DOM,lname),side=3,line=-1.0,cex=2.0,adj=0.1)   
  }
  #abline(v=5)
  #abline(v=12)
  #abline(v=16,lty=2,col="grey40")
  dev.off()
  }
}
