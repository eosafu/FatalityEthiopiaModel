#library(INLA)
library(readxl)
library(tidyverse)
#library(rgdal)
#library(raster)
#library(ggpubr)
#library(stringr)


# Import and filter West and Central Africa
Dat=read_excel("Ethiopia_1997-2022_Jun24.xlsx")
Dat$Vgroup=0
rexgroup <- c("OLA","Oromo Liberation Army","Oromo Liberation Front",
              "ONLF/A","EUFF","OLF","OLF/Shane","SPLA","ONLF","EDF",
              "TPDM","EPPF","OFDM","Ogaden Liberation Front","TPLF",
              "Tigray Peoples Liberation Front","Ogaden rebels", "ONLA",
              "Afar rebels","Ginbot 7","Eritrean soldiers",
              "Gumuz People's Democratic Movement",
              "Gambela Liberation Front (GLF)",
              "Shane fighter")

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
Dat$Ginbot7[Dat$Vgroup==20]=1
Dat$EritSoldiers[Dat$Vgroup==21]=1


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
