library(ggplot2)
library(ggridges)

rawdat <- read.csv("data/all_PL_raw.csv")

# "6187010300" # Unid.chion
# "6187010301" # Opilio
# "6187010302" # Bairdi

prey <- c("6187010301")
pdat <- rawdat[as.character(rawdat$PRED_NODC)=="8791030401" & 
               rawdat$CRUISE_TYPE=="Race_Groundfish" &
               as.character(rawdat$PREY_NODC)%in%prey &
               rawdat$YEAR>=1985 &
               rawdat$MONTH>=6 & rawdat$MONTH<=8
               ,]


ggplot(pdat, aes(x = PREY_SIZE_MM, y = YEAR, group=YEAR)) +
  geom_density_ridges() +
  theme_ridges() + 
  xlim(0,100)
  theme(legend.position = "none")
