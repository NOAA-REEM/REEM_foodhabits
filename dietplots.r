library("ggplot2")
setwd("c:/src/REEM_foodhabits")

dat <- read.csv("results/out_diets_5cm.csv")

sdat <- dat[!(dat$LL %in% c("[5,10)","[100,999)")),c("sp_prey","LL","dietprop")]
sdat$LL <- as.factor(sdat$LL)
sdat$sp_prey <- as.factor(sdat$sp_prey)

agg <- aggregate(dietprop~LL+sp_prey,data=sdat,mean)


sp_order <- c("Amphipods","Euphausiid","Shrimp","Polychaete","Epifauna",
              "Forage.fish","Opilio","Bairdi","Unid.Chion","King.crab","Octopus",
              "Flatfish","W.pollock")

not_shown <- c("Copepods","Oth.nonfish","Rockfish","Discards""Infauna",)

png("out.png",width=1500,height=1000)
yscale <- 5.5; datscale <-5
plot(NULL,xlim=c(-1,nlevels(agg$LL)+1), ylim=c(0,yscale*length(sp_order)+datscale),axes=F,xlab='',ylab='')

for (ss in 1:length(sp_order)){
  d <- agg[agg$sp_prey==sp_order[ss],]
  polygon(c(as.numeric(d$LL),rev(as.numeric(d$LL))),
          c(yscale*ss+d$dietprop*datscale,rev(yscale*ss-d$dietprop*datscale)), 
          col="blue")
  text(0.9,yscale*ss,sp_order[ss],pos=2,cex=1.8)
  text(1:nlevels(agg$LL),0,substr(levels(agg$LL),2,3),cex=1.8)
  lines(c(nlevels(agg$LL)+0.5,nlevels(agg$LL)+0.5),c(7*yscale+datscale,7*yscale-datscale))
  lines(c(nlevels(agg$LL)+0.4,nlevels(agg$LL)+0.6),c(7*yscale+datscale,7*yscale+datscale))
  lines(c(nlevels(agg$LL)+0.4,nlevels(agg$LL)+0.6),c(7*yscale-datscale,7*yscale-datscale))
    text(nlevels(agg$LL)+0.5,7*yscale,"100%",pos=4,cex=1.8)
   mtext("Pacific cod fork length (cm)",side=1,cex=1.8)
   mtext("%diet by weight",side=2,cex=2)
   mtext("Pacific cod - summer BTS",side=3,line=0,cex=2,adj=0.1)
   mtext("all EBS strata 1985-2019",side=3,line=-2,cex=2,adj=0.1)   
}

dev.off()
