library(INLA)
library(readxl)
library(tidyverse)
library(rgdal)
library(raster)
library(ggpubr)
library(stringr)

# Data manipulation
setwd("C:\\Users\\OsafuEgbon\\OneDrive\\Documents\\R\\Rwork\\CrimeEthiopia")
setwd("/Users/egbonoa/Documents/Documents01_09_2023_From_USPLaptop/R/Rwork/CrimeEthiopia")
# Import and filter West and Central Africa
Dat=read_excel("Ethiopia_1997-2022_Jun24.xlsx")
Dat$Vgroup=0
rexgroup <- c("OLA","Oromo Liberation Army","Oromo Liberation Front","ONLF/A","EUFF","OLF","OLF/Shane","SPLA","ONLF","EDF","TPDM","EPPF","OFDM","Ogaden Liberation Front","TPLF","Tigray Peoples Liberation Front","Ogaden rebels", "ONLA", "Afar rebels","Ginbot 7","Eritrean soldiers","Gumuz People's Democratic Movement","Gambela Liberation Front (GLF)","Shane fighter")
rexgroupid <- 1:length(rexgroup)
for (i in 1:nrow(Dat)) {
  aux = str_which(Dat$NOTES[i],rexgroup)
  if(length(aux)==0){
    Dat$Vgroup[i]=0
  }else{
    Dat$Vgroup[i]=sample(aux,1)
  }
    
}
Dat$Vgroup
Dat$OLA=Dat$TPLF=Dat$ONLF=0
Dat$EUFF=Dat$SPLA=Dat$EDF=Dat$TPDM=Dat$AfarRebels=Dat$Ginbot7=Dat$EritSoldiers=0
Dat$GumuzMov=Dat$GLF=0

Dat$OLA[Dat$Vgroup==1|Dat$Vgroup==2|Dat$Vgroup==3 |Dat$Vgroup==5|Dat$Vgroup==6|
          Dat$Vgroup==7|Dat$Vgroup==24] = 1
Dat$TPLF[Dat$Vgroup==15 | Dat$Vgroup==16] = 1
Dat$ONLF[Dat$Vgroup==4 |Dat$Vgroup==9|Dat$Vgroup==14|
    Dat$Vgroup==17|Dat$Vgroup==18]=1
Dat$EUFF[Dat$Vgroup==5]=1
Dat$SPLA[Dat$Vgroup==8]=1
Dat$EDF[Dat$Vgroup==10]=1
Dat$TPDM[Dat$Vgroup==11]=1
#Dat$AfarRebels[Dat$Vgroup==19]=1
Dat$Ginbot7[Dat$Vgroup==20]=1
Dat$EritSoldiers[Dat$Vgroup==21]=1
#Dat$GumuzMov[Dat$Vgroup==22]=1
#Dat$GLF[Dat$Vgroup==23]=1


DDat <- Dat  #%>% filter(YEAR>=2021 & YEAR <= 2022)


# Create sub-event ID to be modeled non-linearly
DDat <- DDat %>% filter(EVENT_TYPE=="Explosions/Remote violence"|
                        EVENT_TYPE=="Violence against civilians"|
                        EVENT_TYPE=="Battles"|
                        EVENT_TYPE=="Riots")
SubEvent <- DDat$SUB_EVENT_TYPE %>%unique()
SubEvent <- data.frame(Event=SubEvent,id=1:length(SubEvent))
DDat$SubEventID <- 0
for (i in 1:nrow(SubEvent)) {
  DDat$SubEventID[DDat$SUB_EVENT_TYPE==SubEvent$Event[i]]=SubEvent$id[i]
}
# population 
Pop <- read.csv("API_SP.POP.TOTL_DS2_en_csv_v2_4019998.csv",skip = 4)
Cname.who <- Pop$Country.Name %>% unique
CnameDat <- DDat$COUNTRY %>% unique()

Pop$id=0
for (i in 1 : length(CnameDat)) {
  Pop$id[Pop$Country.Name==CnameDat[i]] <-1
}
Pop <- Pop %>% filter(id==1)
Pop <- Pop[,-(ncol(Pop))]
############ 
POP <- list()
#
for(k in 1: length(CnameDat)){
  PopNG <- Pop %>% filter(Country.Name==CnameDat[k])%>% unlist(use.names=FALSE)
  PopNG <- PopNG[42:length(PopNG)] 
  PopNG[which(is.na(PopNG))] <- PopNG[min(which(is.na(PopNG))-1)]
  POP[[k]] <- as.numeric(PopNG)
}

DDat$Population=0

for (i in 1:length(CnameDat)) { # ID
  
  for (j in 1997:2021) { # YEAR
    
    DDat$Population[(DDat$COUNTRY==CnameDat[i])&(DDat$YEAR==j)]=POP[[i]][[j-1996]]
    
  }
  
}
DDat$Population[(DDat$COUNTRY==CnameDat[i])&(DDat$YEAR==2022)] <- DDat$Population[(DDat$COUNTRY==CnameDat[i])&(DDat$YEAR==2021)][1]+2884858 # yearly change 

# temperate season


DDat$EVENT_DATE <- DDat$EVENT_DATE %>% as.Date()
DDat$Month <- months(DDat$EVENT_DATE)
DDat$Season <- 0
DDat$Season[DDat$Month=="janeiro" | DDat$Month=="fevereiro" |DDat$Month=="dezembro"] <- "Winter"
DDat$Season[DDat$Month=="mar?o"   | DDat$Month=="abril"     |DDat$Month=="maio"] <- "AASpring"
DDat$Season[DDat$Month=="junho"   | DDat$Month=="julho"     |DDat$Month=="agosto"] <- "Summer"
DDat$Season[DDat$Month=="setembro"   | DDat$Month=="outubro"     |DDat$Month=="novembro"] <- "Autumn"
DDat$AASpring= 0 
DDat$AASpring[DDat$Season=="AASpring"] <- 1
DDat$Summer= 0 
DDat$Summer[DDat$Season=="Summer"] <- 1
DDat$Winter= 0 
DDat$Winter[DDat$Season=="Winter"] <- 1

DDat$binaryFat = ifelse(DDat$FATALITIES>0,1,0)
DDat$countFat = ifelse(DDat$FATALITIES==0,NA,DDat$FATALITIES)

######################
# modelling
######################

DDat.k <- DDat

DDat   <- DDat.k #%>% filter(YEAR > 2016) 

coords = cbind(DDat$LONGITUDE ,DDat$LATITUDE)
mesh <- inla.mesh.2d(coords, max.edge = c(2, 2)) 
#mesh <- inla.mesh.2d(coords,offset = c(2, 2), 
#max.edge = c(2, 2), cutoff =0.1)
k <- DDat$YEAR %>% unique() %>% length()

DDat$YEAR2 <- ((DDat$YEAR)-1996) #%>% cut( breaks=c(-Inf, 10, 13, 18,Inf)) %>% as.factor() %>% as.numeric()
k=DDat$YEAR2  %>% unique()%>% length()
DDat$iddd = 1:length(DDat$YEAR2)

spde <- inla.spde2.pcmatern(mesh = mesh, 
                            prior.range = c(5, 0.01), # P(range < 0.05) = 0.01
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
                            #GLF=DDat$GLF,
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
               #GLF=DDat$GLF,
               YEAR=DDat$YEAR,
               YEAR2=DDat$YEAR2,
               idd=DDat$iddd)),
  tag = "count")

join.stack <- inla.stack(
    stk.zero, stk.count)
  
  
#stk.zinc.pr, stk.lead.pr,
  #stk.shared, stk.zinc.spec)

h.spec <- list(rho = list(prior = 'pc.cor1', param = c(0, 0.9)))
# Model formula
formulae1 <- Y ~ -1+Intercept.count+Intercept.zero+OLA+TPLF+ONLF+EUFF+SPLA+EDF+TPDM+
  Ginbot7+EritSoldiers+Summer+AASpring+
  Winter+E+OLAz+TPLFz+ONLFz+EUFFz+SPLAz+EDFz+TPDMz+
  Ginbot7z+EritSoldiersz+Summerz+AASpringz+
  Winterz+
  f(SubEvent, model="iid") +f(SubEventz, model="iid")+  f(i, copy = "j", fixed = F)+  f(j, model = spde)

formulae2 <- Y ~ -1+Intercept.count+Intercept.zero+OLA+TPLF+ONLF+EUFF+SPLA+EDF+TPDM+
  Ginbot7+EritSoldiers+Summer+AASpring+
  Winter+E+OLAz+TPLFz+ONLFz+EUFFz+SPLAz+EDFz+TPDMz+
  Ginbot7z+EritSoldiersz+Summerz+AASpringz+
  Winterz+
  f(SubEvent, model="iid") +f(SubEventz, model="iid")+  f(i, copy = "j")+  f(j, model = spde) + f(k, model = spde)

formulae3 <- Y ~ -1+Intercept.count+Intercept.zero+OLA+TPLF+ONLF+EUFF+SPLA+EDF+TPDM+
  Ginbot7+EritSoldiers+Summer+AASpring+
  Winter+E+OLAz+TPLFz+ONLFz+EUFFz+SPLAz+EDFz+TPDMz+
  Ginbot7z+EritSoldiersz+Summerz+AASpringz+
  Winterz+
  f(SubEvent, model="iid") +f(SubEventz, model="iid")+  f(i, copy = "j")+  f(j, model = spde) + f(k, model = spde)+
  f(YEAR2,model = "ar1")+ f(YEAR2z,model = "ar1")

formulae <- Y ~ -1 +
   f(i, model = spde, group = i.group, 
    control.group = list(model = 'ar1', hyper = h.spec))+
  f(j, model = spde, group = j.group, 
    control.group = list(model = 'ar1', hyper = h.spec))+
  f(k, copy = "i",fixed = T)

formulae <- Y ~ -1 +
  f(i, copy = "j", fixed = FALSE)+
  f(j, model = spde, group = j.group,
    control.group = list(model = 'ar1', hyper = h.spec)
    )

formulae <- Y ~ -1 +
  f(i, model = spde, group = i.group, 
    control.group = list(model = 'ar1', hyper = h.spec))+
  f(j, model = spde, group = j.group, 
    control.group = list(model = 'ar1', hyper = h.spec))+
  f(k, copy = "i",fixed = T)+f(YEAR2, model = "ar1")+f(YEAR2z, model = "ar1")

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
                      #GLF=DDat$GLF,
                      YEAR=DDat$YEAR,
                      YEAR2=DDat$YEAR2,
                      idd=DDat$iddd)),
  tag = 'stdata0')

h.spec <- list(rho = list(prior = 'pc.cor1', param = c(0, 0.9)))

# Model formula
formulae <- y ~ 0+OLA+TPLF+ONLF+EUFF+SPLA+EDF+TPDM+
  Ginbot7+EritSoldiers+Summer+AASpring+
  Winter+E+
  f(i, model = spde)#, group = i.group, 
                    #control.group = list(model = 'ar1', hyper = h.spec))

# Model fitting
resincnozero <- inla(formulae,  
                     #   E = E,
                     family = "binomial",
                     #family = "nbinomial",
                     #family = "gpoisson",
                     #family = "poisson",
                     data = inla.stack.data(sdat0), 
                     control.predictor = list(compute = TRUE,
                                              A = inla.stack.A(sdat0)), 
                     control.fixed = list(expand.factor.strategy = 'inla'),control.compute = list(config = TRUE,dic = TRUE,waic=TRUE,cpo=TRUE)
                     ,verbose = TRUE)

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
                       #   E = E,
                       #family = "binomial",
                       family = "nbinomial",
                       #family = "gpoisson",
                       #family = "poisson",
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

resincnocountAllYear$iset = iset
resincnocountAllYear$meshh = mesh

resincno <- resincnocountAllYear
#resincno <- resincnobin2010below
summary(resincno)


############## Plot #######################
library(viridis)
path <- read_sf(dsn="Ethiopia.kml")
path = st_coordinates(path)
path = path[,c(1,2)]
############# MAP #############
# A map that leads to country
shp <- st_read(dsn = "/Users/egbonoa/Documents/Documents01_09_2023_From_USPLaptop/R/Rwork/CrimeEthiopia/eth_adm_csa_bofedb_2021_shp/eth_admbnda_adm1_csa_bofedb_2021.shp", stringsAsFactors = F)
subshp <- shp
#shp_df <- broom::tidy(subshp, region = "ADM1_EN")
shp_df <- st_geometry(subshp, region = "ADM1_EN")
# Projection of Spatial Field
boundary <- path

stepsize <- 0.08#4 * 1 / 111

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

#plt <- list()
datresult= NULL
datresult2 = NULL
X = 2010
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
datresult$idcol <- ifelse(datresult$mean_s<0.0,1,1)

datresult2 = rbind(datresult2,datresult)
datresult3 = datresult2 
datresult3$mean_s[datresult3$mean_s< -1]=-1

subshp =  st_coordinates(shp_df) %>%as.data.frame()
#####
p <- ggplot(datresult3, aes(x=x, y=y)) +
  geom_tile(aes(fill = mean_s), alpha=datresult2$idcol) +
  #scale_fill_viridis_c(option = "A", direction = 1,name="") +
  #scale_fill_gradient(low="white", high="darkred")+
  scale_fill_viridis_c(option = "H", direction = 1,name="", guide=guide_colorbar(
    #label.position = "left",
    #direction = "vertical",
    barheight = unit(60, units = "mm"),
    barwidth = unit(1, units = "mm"),
    #title.position = "top",
    #title.vjust=0.5,
    #label.vjust = 0.5
  ))+
  labs(title ="",#paste0(i+X),
       y = "",x="") +
  theme_light()+scale_size(guide = guide_legend(direction = "vertical"))+
  geom_point(data=subshp,aes(x=X,y=Y),color="black",size=0.1)
 # geom_polygon(data = shp_df, aes(x = long, y = lat, group = group), colour = "black", fill = NA)
 p+facet_wrap(~id,scales = "free")+
   theme(strip.background =element_rect(fill="grey"))+
    theme(strip.text = element_text(colour = 'black'))+
    theme(legend.position = "right")
   #####
 

ggsave("AllModel2countOnlySpat.pdf", width = 5.23, height = 4.2, units = "in")
### 
################################################
### Compute uncertainty for the predicted values
###############################################
 m <- inla.posterior.sample(100, resincnobin)
 predprob=NULL
 for (i in 1:100) {
   Apred = m[[i]]$latent[1:4544,]
   predprob=cbind( predprob, exp(Apred)/(1+exp(Apred)))
 }
 predprob = rowMeans(predprob)
 DDat$predprob = predprob
 datprob = DDat %>% select(ADMIN1, predprob)
 datprob2 =datprob %>%group_by(ADMIN1)%>% summarise(Mean=mean(predprob))
 
 datprob2$ADMIN1[datprob2$ADMIN1=="Benshangul/Gumuz"] <- "Benishangul Gumz"
 datprob2$ADMIN1[datprob2$ADMIN1=="South West"] <- "South West Ethiopia"

 cloc3=NULL
 for(i in 1:length(subshp@data$ADM1_EN)){
   cloc3 <- c(cloc3,which(datprob2$ADMIN1==subshp@data$ADM1_EN[i]))
 }
 
 est.dataF <- data.frame(id=subshp@data$ADM1_EN,
                         mean=datprob2$Mean[cloc3])
 
 
 m2 <- inla.posterior.sample(100, resincnocount)
 predprob3=NULL
 no=20   # number of deaths negative binomial parameter
 xi=1.733 # estimated number
 for (i in 1:100) {
   Apred = m2[[i]]$latent[1:3677,]
   rex = exp(Apred)/(1+exp(Apred))
   mu= (rex* xi)/(1-rex)
   rex=rnegbin(rep(1,length(mu)),mu,xi)
   rex= rex>no
   predprob3=cbind( predprob3, rex)
 }
 predprob3 = rowMeans(predprob3)
 DDat$predprob3 = predprob3
 predprob3 = DDat %>% dplyr::select(ADMIN1, predprob3)
 datprob4 =predprob3 %>%group_by(ADMIN1)%>% summarise(Mean=mean(predprob3))
 
 datprob4$ADMIN1[datprob4$ADMIN1=="Benshangul/Gumuz"] <- "Benishangul Gumz"
 datprob4$ADMIN1[datprob4$ADMIN1=="South West"] <- "South West Ethiopia"
 
 cloc3=NULL
 for(i in 1:length(subshp@data$ADM1_EN)){
   cloc3 <- c(cloc3,which(datprob4$ADMIN1==subshp@data$ADM1_EN[i]))
 }
 
 est.dataF2 <- data.frame(id=subshp@data$ADM1_EN,
                         mean=datprob4$Mean[cloc3])
 #est.dataF$component = 1
 #est.dataF2$component = 2
 #Est.data= rbind( est.dataF, est.dataF2)
# Est.data$component = factor(Est.data$component,levels = c(1,2),labels = c("Binary","Count"))
 # A map that leads to country
 shp_df1 = shp_df2 = shp_df
 shp_df1$mean=0
 for (i in 1:nrow(est.dataF)) {
   shp_df1$mean[shp_df1$id==est.dataF$id[i]]<-est.dataF$mean[i]
 }
 shp_df2$mean=0
 for (i in 1:nrow(est.dataF2)) {
   shp_df2$mean[shp_df2$id==est.dataF2$id[i]]<-est.dataF2$mean[i]
 }
 shp_df1$component=1
 shp_df2$component=2
 shp_df = rbind( shp_df1, shp_df2)
 shp_df$component=  factor(shp_df$component,levels = c(1,2),labels = c("Binary","Count"))
p= ggplot() + geom_polygon(data = shp_df, aes(x = long, y = lat, group = group, fill = mean), colour = "black") + 
   theme_void()+
   scale_fill_viridis_c(option = "H", direction = 1,name="", guide=guide_colorbar(
     #label.position = "left",
     #direction = "vertical",
     barheight = unit(60, units = "mm"),
     barwidth = unit(1, units = "mm"),
     #title.position = "top",
     #title.vjust=0.5,
     #label.vjust = 0.5
   ))+
   labs(title = "",x="",y="")+theme_bw()
 p+facet_wrap(~component)
 ggsave("uncertainty.pdf", width = 8.23, height = 3.9, units = "in")
 ###
library(xtable)

tab=data.frame (Event=SubEvent$Event,
                mean=round(resincno$summary.random$SubEvent$mean,3),
                lw=round(resincno$summary.random$SubEvent$`0.025quant`,3),
                up=round(resincno$summary.random$SubEvent$`0.975quant`,3))


tab <- tab %>% arrange(desc(mean))
xtable(tab, type="latex",digits = 3)

###  plot the event
tab$id=0
tab_bin = tab
tab=rbind(tab_bin,tab_count)
tab$sig =unlist(lapply(1:nrow(tab),function(i) findInterval(0,c(tab$lw[i],tab$up[i]))))
tab$sig2=tab$sig
tab$sig2[tab$sig==1]=0.3
tab$sig2[tab$sig!=1]=1
tab$id = factor(tab$id, levels = c(0,1),labels = c("Binary", "Count"))
tab = tab[order(tab$mean),]
tab$mean = round(tab$mean,digits = 2)
ggplot(tab,aes(x=Event, y=mean, label=mean)) + 
  geom_point( aes(color=id), size=4, alpha=tab$sig2)  +
  #geom_segment(aes(y = 0, 
  #                 x = Event, 
  #                 yend = mean, 
  #                 xend = Event), 
  #             color = "black") +
  geom_text(color="black", size=1) +geom_hline(yintercept = 0,alpha=0.2)+
  labs(title=" ") + 
  coord_flip()+
  theme_bw() + theme(legend.title=element_blank())+labs(y="Effect",x="Event type")

ggsave("eventtype.pdf", width = 7.23, height = 5.2, units = "in")

###
# Fixed Effect
tab=resincno$summary.fixed[c(1,3,5)] %>% round(digits = 3)
tab <- tab %>% arrange(desc(mean))
xtable(tab, type="latex",digits = 3)
### plot actors ###
# count
tab_count= tab[-c(7,6,8,12),]
a = rownames(tab_count)
a[3] = "Eritrea Army"
rownames(tab_count) = a
a= rownames(tab_count)
tab_count$name=a
#binary
tab_bin= tab[-c(7,10,8,11),]
a = rownames(tab_bin)
a[2] = "Eritrea Army"
rownames(tab_bin) = a
a =rownames(tab_bin)
tab_bin$name=a
###
tab_bin$id=1
tab_count$id=2
tab=rbind(tab_bin,tab_count)
tab$id= factor(tab$id,levels = c(1,2),labels = c("Binary", "Count"))
tab$x= c(1:9,1:9)
tab %>% ggplot()+
  geom_line(aes(x=x,y=mean,color=id))+
  geom_line(aes(x=x,y=`0.025quant`,color=id),linetype="dashed",alpha=0.3)+
  geom_point(aes(x=x,y=`0.025quant`,color=id),linetype="dashed",alpha=0.1)+
  geom_line(aes(x=x,y=`0.975quant`,color=id),linetype="dashed",alpha=0.3)+
  geom_point(aes(x=x,y=`0.975quant`,color=id),linetype="dashed",alpha=0.1)+
  geom_point(aes(x=x,y=mean,color=id))+theme_bw() +labs(x="Rank position",y="Effect")+
  annotate("text", x = tab$x, y=tab$mean, label = tab$name,size=3.5)+
  theme(legend.title=element_blank())
ggsave("actors.pdf", width = 7.23, height = 4.2, units = "in")
#non-linear
tab=data.frame(id=resincno$summary.random$YEAR$ID,mean=resincno$summary.random$YEAR$mean,
           lw=resincno$summary.random$YEAR$`0.025quant`,up=resincno$summary.random$YEAR$`0.975quant`)

tab %>% ggplot()+
  geom_line(aes(x=id,y=mean))+geom_point(aes(x=id,y=mean))+
geom_line(aes(x=id,y=lw),color="red",linetype="dashed")+
  geom_line(aes(x=id,y=up),color="red",linetype="dashed")+
  theme_bw()+labs(x="Year",y="Effect on fatality")

ggsave("nonlinearbin.pdf", width = 7.23, height = 4.2, units = "in")


### Hyperparameters

xtable(resincnobin2010above$summary.hyperpar[-1,c(1,3,5)], type="latex",digits = 3)
######### Descriptive ########

### Plot descriptive

DDAT <- DDat%>% group_by(ADMIN1)

FAT <- DDAT %>% summarise(Fat= sum(FATALITIES), Event=n())
aux <- DDat %>% filter(YEAR==2022)
aux <- aux %>% group_by(ADMIN1)
rex <- aux %>% summarise(Pop=mean(Population))
FAT <- FAT %>% mutate(Pop=rex$Pop)
FAT <- FAT %>% mutate(rate= (Fat/Pop)*100000, eventRate=(Event/Pop)*100000,motarate=Fat/Event)
FAT
FAT$ADMIN1[FAT$ADMIN1=="Benshangul/Gumuz"] <- "Benishangul Gumz"
FAT$ADMIN1[FAT$ADMIN1=="South West"] <- "South West Ethiopia"



cloc3=NULL
for(i in 1:length(subshp@data$ADM1_EN)){
  cloc3 <- c(cloc3,which(FAT$ADMIN1==subshp@data$ADM1_EN[i]))
}

est.dataF <- data.frame(id=subshp@data$ADM1_EN,
                        mean=FAT$motarate[cloc3])

# A map that leads to country
shp_df$mean=0
for (i in 1:nrow(est.dataF)) {
  shp_df$mean[shp_df$id==est.dataF$id[i]]<-est.dataF$mean[i]
}


plt[[1]] <- ggplot() + geom_polygon(data = shp_df, aes(x = long, y = lat, group = group, fill = mean), colour = "black") + 
  theme_void()+scale_fill_gradient(name="",low="white", high="red")+
  labs(title = "",x="",y="")+theme_bw()

## Descriptive per year

DDAT <- DDat%>% group_by(YEAR)

FAT <- DDAT %>% summarise(Fat= sum(FATALITIES), Event=n(),Pop=mean(Population))
FAT <- FAT %>% mutate(rate= (Fat/Pop)*100000, eventRate=(Event/Pop)*100000,motarate=Fat/Event)

FAT%>% ggplot()+
  geom_line(aes(x=YEAR,y=motarate),size=1.0)+
  geom_point(aes(x=YEAR,y=motarate),size=1.0,color="red")+
  geom_smooth(aes(x=YEAR,y=motarate),size=1.0, color="orange",alpha=0.2)+
  labs(x="Year", y="Fatality per event")+
  theme_bw()

### Barplot

DDat %>%filter(FATALITIES<50)%>% ggplot(aes(x=(FATALITIES)))+
  geom_bar( width=0.7, fill="steelblue",size=1.5)+
  theme_bw()+labs(y="Count",x="Fatality")


## Word plot
#https://monkeylearn.com/word-cloud
a= paste(DDat$NOTES[sample(1:nrow(Dat),5000)] , collapse = '')
write.csv(a,"note.txt")

a= paste(DDat$ACTOR1[sample(1:nrow(Dat),5000)] , collapse = '')
write.csv(a,"Actor2.txt")

# Plot Hist by yearly group.

DDat$YERRGRP=0
DDat$YERRGRP[DDat$YEAR<=2005]=1
DDat$YERRGRP[DDat$YEAR>2005 & DDat$YEAR<=2016]=2
DDat$YERRGRP[DDat$YEAR>2016]=3
DDat$YERRGRP = factor(DDat$YERRGRP,levels = c(1,2,3),labels = c("1997-2005","2006-2016","2017-2022"))
p= DDat%>%filter(FATALITIES<200) %>%
  ggplot( aes(x=FATALITIES)) +
  geom_histogram(color="white",fill="black", alpha=0.5, position="identity")+labs(x="Fatality count per event",y= "Count",title = "")+theme_bw()
p+facet_wrap(~YERRGRP,scales = "free")


ggsave("barplot.pdf", width = 6.23, height = 2.9, units = "in")


# Plot Hist by yearly group overlay poison fit.

DDat$YERRGRP=0
DDat$YERRGRP[DDat$YEAR<=2005]=1
DDat$YERRGRP[DDat$YEAR>2005 & DDat$YEAR<=2016]=2
DDat$YERRGRP[DDat$YEAR>2016]=3
DDat$YERRGRP = factor(DDat$YERRGRP,levels = c(1,2,3),labels = c("1997-2005","2006-2016","2017-2022"))
p= DDat%>%filter(FATALITIES<200) %>%
  ggplot( aes(x=FATALITIES)) +
  geom_histogram(color="white",fill="black", alpha=0.5, position="identity")+labs(x="Fatality count per event",y= "Count",title = "")+theme_bw()
p+facet_wrap(~YERRGRP,scales = "free")

an ="2017-2022"
histDist(DDat$FATALITIES[DDat$YERRGRP==an&DDat$FATALITIES<200],family = "PO",main = an)
ggsave("barplot.pdf", width = 6.23, height = 2.9, units = "in")
