library(INLA)
library(readxl)
library(tidyverse)
library(rgdal)
library(raster)
library(ggpubr)
library(stringr)

source("Datawrangling.R")

######################
# modelling
######################

DDat.k <- DDat

DDat   <- DDat.k 

coords = cbind(DDat$LONGITUDE ,DDat$LATITUDE)
mesh <- inla.mesh.2d(coords, max.edge = c(2, 2)) 

k <- DDat$YEAR %>% unique() %>% length()

DDat$YEAR2 <- ((DDat$YEAR)-1996) 
k=DDat$YEAR2  %>% unique()%>% length()
DDat$iddd = 1:length(DDat$YEAR2)

spde <- inla.spde2.pcmatern(mesh = mesh, 
                            prior.range = c(5, 0.01), # P(range < 5) = 0.01
                            prior.sigma = c(1, 0.01))   # P(sigma > 1) = 0.01

iseti <- inla.spde.make.index('i', n.spde = spde$n.spde,
                              n.group = k)
isetj <- inla.spde.make.index('j', n.spde = spde$n.spde,
                              n.group = k)
isetk <- inla.spde.make.index('k', n.spde = spde$n.spde,
                              n.group = k)

A <- inla.spde.make.A(mesh = mesh,
                      loc = cbind(DDat$LONGITUDE, DDat$LATITUDE), group = DDat$YEAR2) 

###############
# Create the stack object for lead
stk.zero <- inla.stack(
  data = list(Y = cbind(DDat$binaryFat, NA)),
  A = list(A, 1),
  effects = list(iseti,
                 data.frame(Intercept.zero = 1,
                            Ez=log(DDat$Population),
                            SubEventz=DDat$SubEventID,
                            Summerz=DDat$Summer,
                            AASpringz=DDat$AASpring,
                            Winterz=DDat$Winter,
                            OLAz=DDat$OLA,
                            TPLFz=DDat$TPLF,
                            ONLFz=DDat$ONLF,
                            EUFFz=DDat$EUFF,
                            SPLAz=DDat$SPLA,
                            EDFz=DDat$EDF,
                            TPDMz=DDat$TPDM,
                            AfarRebelsz=DDat$AfarRebels,
                            Ginbot7z=DDat$Ginbot7,
                            EritSoldiersz=DDat$EritSoldiers,
                            GumuzMovz=DDat$GumuzMov,
                            YEARz=DDat$YEAR,
                            YEAR2z=DDat$YEAR2,
                            iddz=DDat$iddd)),
  tag = "Zero")

# Create the stack object for zinc
stk.count <- inla.stack(
  data = list(Y= cbind(NA, DDat$countFat)),
  A = list(A, A, 1),
  effects = list(
    isetj, isetk,
    data.frame(Intercept.count = 1, 
               E=log(DDat$Population),
               SubEvent=DDat$SubEventID,
               Summer=DDat$Summer,
               AASpring=DDat$AASpring,
               Winter=DDat$Winter,
               OLA=DDat$OLA,
               TPLF=DDat$TPLF,
               ONLF=DDat$ONLF,
               EUFF=DDat$EUFF,
               SPLA=DDat$SPLA,
               EDF=DDat$EDF,
               TPDM=DDat$TPDM,
               AfarRebels=DDat$AfarRebels,
               Ginbot7=DDat$Ginbot7,
               EritSoldiers=DDat$EritSoldiers,
               GumuzMov=DDat$GumuzMov,
               YEAR=DDat$YEAR,
               YEAR2=DDat$YEAR2,
               idd=DDat$iddd)),
  tag = "count")

join.stack <- inla.stack(
  stk.zero, stk.count)


h.spec <- list(rho = list(prior = 'pc.cor1', param = c(0, 0.9)))

# Model With Joint modeling as in Asmarian et. al.2019 

formulae <- Y ~ -1+Intercept.count+Intercept.zero+OLA+TPLF+ONLF+EUFF+SPLA+EDF+TPDM+
  Ginbot7+EritSoldiers+Summer+AASpring+Winter+E+
  
  OLAz+TPLFz+ONLFz+EUFFz+SPLAz+EDFz+TPDMz+
  Ginbot7z+EritSoldiersz+Summerz+AASpringz+
  Winterz+ f(SubEvent, model="iid") +f(SubEventz, model="iid")


formulae <- update(formulae,  ~. +
                     
  f(i, copy = "j")+
  f(j, model = spde, group = j.group, 
    control.group = list(model = 'ar1', hyper = h.spec))
           )


# Model fitting
resincnocount_3nbinomial <- inla(formulae,  
                                 family = c("binomial","nbinomial"),
                                 data = inla.stack.data(join.stack), 
                                 control.predictor = list(compute = TRUE,
                                                          A = inla.stack.A(join.stack)), 
                                 control.fixed = list(expand.factor.strategy = 'inla'),
                                 control.compute = list(config = TRUE,dic = TRUE,waic=TRUE,cpo=TRUE)
                                 ,verbose = TRUE)


###############
# Sequential Approach
###############

sdat0 <- inla.stack(
  data = list(y = DDat$binaryFat), 
  A = list(A,1), 
  effects = list(c(iseti,list(Intercept=1)),
                 list(E=log(DDat$Population),
                      SubEvent=DDat$SubEventID,
                      Summer=DDat$Summer,
                      AASpring=DDat$AASpring,
                      Winter=DDat$Winter,
                      OLA=DDat$OLA,
                      TPLF=DDat$TPLF,
                      ONLF=DDat$ONLF,
                      EUFF=DDat$EUFF,
                      SPLA=DDat$SPLA,
                      EDF=DDat$EDF,
                      TPDM=DDat$TPDM,
                      AfarRebels=DDat$AfarRebels,
                      Ginbot7=DDat$Ginbot7,
                      EritSoldiers=DDat$EritSoldiers,
                      GumuzMov=DDat$GumuzMov,
                      YEAR=DDat$YEAR,
                      YEAR2=DDat$YEAR2,
                      idd=DDat$iddd)),
  tag = 'stdata0')

h.spec <- list(rho = list(prior = 'pc.cor1', param = c(0, 0.9)))

# Model formula
formulae <- y ~ 0+OLA+TPLF+ONLF+EUFF+SPLA+EDF+TPDM+
  Ginbot7+EritSoldiers+Summer+AASpring+
  Winter+
  f(i, model = spde, group = i.group, 
 control.group = list(model = 'ar1', hyper = h.spec))

# Model fitting
resincnozero <- inla(formulae,  
                     family = "binomial",
                     data = inla.stack.data(sdat0), 
                     control.predictor = list(compute = TRUE,
                                              A = inla.stack.A(sdat0)), 
                     control.fixed = list(expand.factor.strategy = 'inla'),control.compute = list(config = TRUE,dic = TRUE,waic=TRUE,cpo=TRUE)
                     ,verbose = TRUE)

#####################
#Determine optimal c
#####################
n = nrow(DDat)
fittedProb =  1-resincnozero$summary.fitted.values$mean[1:n]

DIC=WAIC =NULL
KK = seq(0,1,length = 40)
for (k in KK) {
  
  aux = ifelse(fittedProb < k,TRUE,FALSE)
  Y = DDat$countFat
  aux.1 = which(is.na(Y))
  aux.sub = aux[aux.1]
  Y[aux.1][aux.sub] = 0
  
  sdat1 <- inla.stack(
    data = list(y =Y ), 
    A = list(A,1), 
    effects = list(c(iseti,list(Intercept=1)),
                   list(E=log(DDat$Population),
                        SubEvent=DDat$SubEventID,
                        Summer=DDat$Summer,
                        AASpring=DDat$AASpring,
                        Winter=DDat$Winter,
                        OLA=DDat$OLA,
                        TPLF=DDat$TPLF,
                        ONLF=DDat$ONLF,
                        EUFF=DDat$EUFF,
                        SPLA=DDat$SPLA,
                        EDF=DDat$EDF,
                        TPDM=DDat$TPDM,
                        AfarRebels=DDat$AfarRebels,
                        Ginbot7=DDat$Ginbot7,
                        EritSoldiers=DDat$EritSoldiers,
                        GumuzMov=DDat$GumuzMov,
                        #GLF=DDat$GLF,
                        YEAR=DDat$YEAR,
                        YEAR2=DDat$YEAR2,
                        idd=DDat$iddd)),
    tag = 'stdata1')
  
  resincnozero <- inla(formulae,  
                       family = "nbinomial",
                       data = inla.stack.data(sdat0), 
                       control.predictor = list(compute = TRUE,
                                                A = inla.stack.A(sdat0)), 
                       control.fixed = list(expand.factor.strategy = 'inla'),control.compute = list(config = TRUE,dic = TRUE,waic=TRUE,cpo=TRUE)
                       ,verbose = TRUE)
  
  DIC = c(DIC,resincnozero$dic$dic)
  WAIC = c(WAIC,resincnozero$waic$waic)
}
par(mfrow=c(1,2))
plot(KK,DIC, xlab= "c",ylab="DIC",type="b",ylim=c(10370,10370.4))
plot(KK,WAIC, xlab= "c",ylab="WAIC",type="b",ylim= c(10349.6,10350.1))


#####################
#Determine optimal c through optim
#####################
n = nrow(DDat)
fittedProb =  1-resincnozero$summary.fitted.values$mean[1:n]

GD = function(k) {
  
  aux = ifelse(fittedProb < k,TRUE,FALSE)
  Y = DDat$countFat
  aux.1 = which(is.na(Y))
  aux.sub = aux[aux.1]
  Y[aux.1][aux.sub] = 0
  
  sdat1 <- inla.stack(
    data = list(y =Y ), 
    A = list(A,1), 
    effects = list(c(iseti,list(Intercept=1)),
                   list(E=log(DDat$Population),
                        SubEvent=DDat$SubEventID,
                        Summer=DDat$Summer,
                        AASpring=DDat$AASpring,
                        Winter=DDat$Winter,
                        OLA=DDat$OLA,
                        TPLF=DDat$TPLF,
                        ONLF=DDat$ONLF,
                        EUFF=DDat$EUFF,
                        SPLA=DDat$SPLA,
                        EDF=DDat$EDF,
                        TPDM=DDat$TPDM,
                        AfarRebels=DDat$AfarRebels,
                        Ginbot7=DDat$Ginbot7,
                        EritSoldiers=DDat$EritSoldiers,
                        GumuzMov=DDat$GumuzMov,
                        YEAR=DDat$YEAR,
                        YEAR2=DDat$YEAR2,
                        idd=DDat$iddd)),
    tag = 'stdata1')
  
  resincno <- inla(formulae,  
                       family = "nbinomial",
                       data = inla.stack.data(sdat1), 
                       control.predictor = list(compute = TRUE,
                                                A = inla.stack.A(sdat1)), 
                       control.fixed = list(expand.factor.strategy = 'inla'),control.compute = list(config = TRUE,dic = TRUE,waic=TRUE,cpo=TRUE)
                       ,verbose = TRUE)

 m2 <- inla.posterior.sample(100, resincno)
  predprob3=NULL
  no=20   # number of deaths negative binomial parameter
  xi=1/resincno$summary.hyperpar[1,1]
  lppd = NULL
  for (i in 1:100) {
    Apred = exp(m2[[i]]$latent[grep("APredictor",rownames(m2[[i]]$latent )),])
    rex=dnbinom(rep(1,length(Apred)),mu=Apred,size = xi)
    lppd=cbind( lppd, rex)
  }
  lppd1 = rowMeans(lppd) %>% log
  lppd1 = sum( lppd1[(sdat1$data$data!=0 & !is.na(sdat1$data$data))])
  lapp2 = log(lppd)
  lapp2 = apply(lapp2,1,function(x) var(x))
  lapp2 = sum(lapp2)
  Waic = -2*(lppd1 - lapp2)
  return(Waic)
}
res = nlminb(start = 0.01, objective = GD, lower = c(0), 
             upper = c(1))

resincnocountAllYear$iset = iset
resincnocountAllYear$meshh = mesh

resincno <- resincnocountAllYear
summary(resincno)


############## Plot #######################
library(viridis)
path <- read_sf(dsn="Ethiopia.kml")
path = st_coordinates(path)
path = path[,c(1,2)]
############# MAP #############

# A map that leads to country (Check repository to download "shp_df.Rdata")
shp <- st_read(dsn = "eth_admbnda_adm1_csa_bofedb_2021.shp", stringsAsFactors = F)
subshp <- shp
shp_df <- st_geometry(subshp, region = "ADM1_EN")

# Projection of Spatial Field
boundary <- path

stepsize <- 0.08 

k = 1
nxy <- round(
  c(diff(range( boundary[, 1])), 
    diff(range( boundary[, 2]))) / stepsize)
projgrid <- inla.mesh.projector(
  mesh, xlim = range(boundary[, 1]), 
  ylim = range(boundary[, 2]), dims = nxy)

xmean <- list()



for (j in 1:k){
  xmean[[j]] <- inla.mesh.project(
    projgrid, resincno$summary.random$j$mean[isetj$j.group == j])
}

df <-  expand.grid(x = projgrid$x, y = projgrid$y)

datresult= NULL
datresult2 = NULL
X = 1996
for (i in 1:k) {
  
  df$mean_s <- as.vector(xmean[[i]])
  
  ind <- point.in.polygon(
    df$x, df$y,
    boundary[, 1], boundary[, 2]
  )
  
  dff <- df[which(ind == 1), ]
  dff$id=i+X
  datresult = rbind(datresult,dff)
}

datresult$idd  <- datresult$id
datresult2 = rbind(datresult2,datresult)
datresult3 = datresult2 
datresult3$mean_s[datresult3$mean_s< -1]=-1

subshp =  st_coordinates(shp_df) %>%as.data.frame()
#####
p <- ggplot(datresult3, aes(x=x, y=y)) +
  geom_tile(aes(fill = mean_s)) +

  scale_fill_viridis_c(option = "H", direction = 1,name="", guide=guide_colorbar(
    barheight = unit(60, units = "mm"),
    barwidth = unit(1, units = "mm"),
  ))+
  labs(title ="",#paste0(i+X),
       y = "",x="") +
  theme_light()+scale_size(guide = guide_legend(direction = "vertical"))+
  geom_point(data=subshp,aes(x=X,y=Y),color="black",size=0.1)
p+facet_wrap(~id,scales = "free")+
  theme(strip.background =element_rect(fill="grey"))+
  theme(strip.text = element_text(colour = 'black'))+
  theme(legend.position = "right")
#####


ggsave("AllModel2countOnlySpat.pdf", width = 5.23, height = 4.2, units = "in")
### 
