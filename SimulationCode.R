##############################################################################
#### This code shows an example of estimating a zero-inflated 
#### count model with R-INLA using joint and sequential estimation approaches
##############################################################################
library(INLA)
library(fields)
library(ggplot2)
library(viridisLite)
library(sf)

rm(list=ls())
set.seed(2010)
inla.seed = sample.int(n=1E7, size=1)
options(width=70, digits=3)

##############################################
####(1) Simulate Binary and count spatial fields
##############################################

  ##########################
  ######## Binary component
  ##########################
  # Marginal standard deviation
  sigma.u_0 = 0.2
  # - range
  range_0 = 1

  # Import Ethiopia boundary coordinates
  path <- read_sf(dsn="Ethiopia.kml")
  path = st_coordinates(path)
  path = path[,c(1,2)]
  
  # > c(min(DDat$LONGITUDE),min(DDat$LATITUDE),max(DDat$LONGITUDE),max(DDat$LONGITUDE),
  # +                             min(DDat$LONGITUDE), max(DDat$LONGITUDE), max(DDat$LONGITUDE), min(DDat$LATITUDE))
  # [1] 33.086  3.450 45.529 45.529 33.086 45.529 45.529  3.450
  #
  sample.locations = matrix(c(33.086,  3.450, 45.529, 45.529, 33.086, 45.529, 45.529,  3.450), nrow = 4, byrow = T)
  mesh.sim = inla.mesh.2d(loc = sample.locations, max.edge=c(2, 2))
  
  plot(mesh.sim)
  axis(1); axis(2)
  # Sample from SPDE model
  spde = inla.spde2.pcmatern(mesh.sim, prior.range = c(.5, .5), prior.sigma = c(.5, .5))
  
  Qu = inla.spde.precision(spde, theta=c(log(range_0), log(sigma.u_0)))
  u0 = inla.qsample(n=1, Q=Qu, seed = inla.seed)
  u0 = u0[ ,1]
  
  ##########################
  ######## Count component
  ##########################
  
  sigma.u_1 = 1.0
  # - the marginal standard deviation of the spatial field
  range_1 = 4
  
  Qu = inla.spde.precision(spde, theta=c(log(range_1), log(sigma.u_1)))
  u1 = inla.qsample(n=1, Q=Qu, seed = inla.seed)
  u1 = u1[ ,1]
  ###################################
  #### Plot simulates spatial fields
  ##################################
  boundary <- path
  
  stepsize <- 0.08 
  
  nxy <- round(
    c(diff(range( boundary[, 1])), 
      diff(range( boundary[, 2]))) / stepsize)
  projgrid <- inla.mesh.projector(
    mesh.sim, xlim = range(boundary[, 1]), 
    ylim = range(boundary[, 2]), dims = nxy)
  
## Project Spatial effect for binary component
  
   xmean <- inla.mesh.project(
             projgrid, u0)

  
  df <-  expand.grid(x = projgrid$x, y = projgrid$y)
  
    
    df$mean_s <- as.vector(xmean)
    
    ind <- point.in.polygon(
      df$x, df$y,
      boundary[, 1], boundary[, 2]
    )
    
    dff <- df[which(ind == 1), ]
##############################################
#### Import Ethiopia shp file for plotting
##############################################
subshp <- st_read(dsn = "/Users/egbonoa/Documents/Documents01_09_2023_From_USPLaptop/R/Rwork/CrimeEthiopia/eth_adm_csa_bofedb_2021_shp/eth_admbnda_adm1_csa_bofedb_2021.shp", stringsAsFactors = F)

shp_df <- st_geometry(subshp, region = "ADM1_EN")
    
subshp =  st_coordinates(shp_df) %>%as.data.frame()

  #####
  p0 <- ggplot(dff, aes(x=x, y=y)) +
    geom_tile(aes(fill = mean_s)) +
     
    scale_fill_viridis_c(option = "H", direction = 1,name="", guide=guide_colorbar(
      
      barheight = unit(60, units = "mm"),
      barwidth = unit(1, units = "mm"),
       
    ))+
    labs(title ="Count component",#paste0(i+X),
         y = "",x="") +
    theme_light()+scale_size(guide = guide_legend(direction = "vertical"))+
    geom_point(data=subshp,aes(x=X,y=Y),color="black",size=0.1)

## Project Spatial effect for count component

xmean <- inla.mesh.project(
  projgrid, u1)


df <-  expand.grid(x = projgrid$x, y = projgrid$y)


df$mean_s <- as.vector(xmean)

ind <- point.in.polygon(
  df$x, df$y,
  boundary[, 1], boundary[, 2]
)

dff <- df[which(ind == 1), ]


#####
p1 <- ggplot(dff, aes(x=x, y=y)) +
  geom_tile(aes(fill = mean_s)) +
  
  scale_fill_viridis_c(option = "H", direction = 1,name="", guide=guide_colorbar(
    
    barheight = unit(60, units = "mm"),
    barwidth = unit(1, units = "mm"),
    
  ))+
  labs(title ="Count component",
       y = "",x="") +
  theme_light()+scale_size(guide = guide_legend(direction = "vertical"))+
  geom_point(data=subshp,aes(x=X,y=Y),color="black",size=0.1)

library(patchwork)
  p0|p1

#############################
#### Generate response Data##
#############################
#  > min(DDat$LONGITUDE);max(DDat$LONGITUDE)
#  [1] 33.086
#  [1] 45.529  
#  > min(DDat$LATITUDE);max(DDat$LATITUDE)
#  [1] 3.45
#  [1] 14.695
  
# generate locations
long  = seq(33.086,45.529,length = 100)
lat   = seq(3.45,14.695,length = 100)
loc.data = expand.grid(long,lat) %>%as.matrix()
n = nrow(loc.data)

# - Derive projection matrix A
  
A = inla.spde.make.A(mesh=mesh.sim, loc=loc.data)
# Derive spatial effects from spatial field
ux0 = drop(A %*% u0)
ux1 = drop(A %*% u1)
  

# -Generate covaariates
set.seed(43567)
  x = rnorm(n,0,1)
  w = rbinom(n,1,0.5)
  # - Fix coefficients for binary and count components
  beta0 = c(1, -1.5, 0.5) 
  beta1 = c(2, 0.5,-1.0) 
  
# - calculate linear predictors
#############################
#############################
#  Simulation scheme 1 (where the same spatial effect generates both binary and count components)
  
  u = ux1
  lin.pred0  = beta0[1] + beta0[2]*x  +beta0[3]*w  +  0.4*u 
  lin.pred1 =  beta1[1] + beta1[2]*x+  beta1[3]*w  + u 
#############################
#############################
  #  Simulation scheme 2 (where the different spatial effects generate the binary and count components)(Uncomment to run simulation scheme 2)
  

 # lin.pred0  = beta0[1] + beta0[2]*x  +beta0[3]*w  +   ux0
 # lin.pred1 = beta0[1] + beta0[2]*x+  beta0[3]*w  +   ux1
  #############################
  #############################


  pix = exp(lin.pred0)/(1+exp(lin.pred0))
  lamdax = exp(lin.pred1)
  Y = NULL
  for(i in 1:n ){
    z0 =  rbinom(1,1,prob = pix[i])
    y=ifelse(z0==1,rpois(1,lambda =lamdax[i]),0 )
    Y = c(Y,y)
  }
  
  
  ##################
  #Estimation begins
  ##################
  #############################
  # (2) Estimation using a joint model (Asmarian, 2019)
  ############################# 
  
  mesh = inla.mesh.2d(loc = loc.data, max.edge=c(2, 2))
  A = inla.spde.make.A(mesh=mesh, loc=loc.data)
  
  plot(mesh, asp=1)
  points(loc.data, col='red', lwd=.1)
  
  spde = inla.spde2.pcmatern(mesh, 
                             prior.range = c(1, 0.1), 
                             prior.sigma = c(1, 0.1))
  
  
  iseti <- inla.spde.make.index('i', n.spde = spde$n.spde)
  isetj <- inla.spde.make.index('j', n.spde = spde$n.spde)
  isetk <- inla.spde.make.index('k', n.spde = spde$n.spde)
  
  z0 = ifelse(Y >0,1,0)
  z1 = ifelse(Y >0,Y,NA)
  
stack0 = inla.stack(tag='est0',

                      data=list(Y=cbind(z0,NA)),
                      effects=list(
                        iseti, 
                        data.frame(intercept0=1, x0=x,w0 = w)),
                      A=list(A, 1)
)
stack1 = inla.stack(tag='est1',
                      data=list(Y=cbind(NA,z1)),
                      effects=list(
                        isetj, isetk,
                        data.frame(intercept1=1, x1=x,w1=w)),
                      A=list(A,A, 1)
  )
  join.stack <- inla.stack(
    stack0, stack1)
  

  formulae <- Y ~ -1 +intercept0+intercept1+x0+x1+w0+w1+
    f(i, copy = "j",fixed=F)+
    f(j, model = spde)
  
  sim1 <- inla(formulae,  
               family = c("binomial","poisson"),
               data = inla.stack.data(join.stack), 
               control.predictor = list(compute = TRUE,
                                        A = inla.stack.A(join.stack)), 
               control.fixed = list(expand.factor.strategy = 'inla'),
               control.compute = list(config = TRUE,
                                      dic = TRUE,
                                      waic=TRUE,
                                      cpo=TRUE)
               ,verbose = TRUE) 
  
  summary(sim1)


  ######################################
  #(3) Estimation using sequential framework
  ######################################
  

  spde0 = inla.spde2.pcmatern(mesh, 
                             prior.range = c(1, .01), 
                             prior.sigma = c(1, .01))
  
  spde1 = inla.spde2.pcmatern(mesh, 
                              prior.range = c(1, .07), 
                              prior.sigma = c(1, .07))
  
  iseti <- inla.spde.make.index('i', n.spde = spde$n.spde )
  isetj <- inla.spde.make.index('j', n.spde = spde$n.spde )
  
  z0 = ifelse(Y >0,1,0)
  z1 = ifelse(Y >0,Y,NA)
  
  stack0 = inla.stack(tag='est0',
                      # - Name of the stack
                      data=list(Y=z0),
                      effects=list(
                        # - The spde Components
                        iseti, 
                        # - Data frame for the linear and non-linear effects as the case may be
                        data.frame(intercept0=1, x0=x,w0=w)),
                      # - The projection matrix A or C in the article
                      A=list(A, 1)
                      # - 1 here represents intercept.
  )
  

  formulae <- Y ~ -1 +intercept0+x0+w0+
    f(i, model = spde0)
  
  sim20 <- inla(formulae,  
               family = c("binomial"),
               data = inla.stack.data(stack0), 
               control.predictor = list(compute = TRUE,
                                        A = inla.stack.A(stack0)), 
               control.fixed = list(expand.factor.strategy = 'inla'),
               control.compute = list(config = TRUE,dic = TRUE,waic=TRUE,cpo=TRUE)
               ,verbose = TRUE) 
  
  
 summary(sim20)

 #### Optimize to obtain c
  
  n = length(Y)
  fittedProb =  1-sim20$summary.fitted.values$mean[1:n]
  
  GD = function(k) {
    z =  z1
    aux = ifelse(fittedProb < k,TRUE,FALSE)
    aux.1 = which(is.na(z1))
    aux.sub = aux[aux.1]
    z[aux.1][aux.sub] = 0
    
    stack1 = inla.stack(tag='est1',
                        # - Name of the stack
                        data=list(Y=z),
                        effects=list(
                          # - The spde Components
                          isetj,
                          # - Dataframe for the linear and non-linear effects as the case may be
                          data.frame(intercept1=1, x1= x,w1=w)),
                        # - The projection matrix A or C in the article
                        A=list(A, 1)
                        # - 1 here represents intercept.
    )
    formulae <- Y ~ -1 +intercept1+x1+w1+
      f(j, model = spde)
    
    sim21 <- inla(formulae,  
                  family = c("poisson"),
                  data = inla.stack.data(stack1), 
                  control.predictor = list(compute = TRUE,
                                           A = inla.stack.A(stack1)), 
                  control.fixed = list(expand.factor.strategy = 'inla'),
                  control.compute = list(config = TRUE,dic = TRUE,waic=TRUE,cpo=TRUE)
                  ,verbose = TRUE) 
    
    return(sim21$waic$waic)
  }
  res = nlminb(start = 0.01, objective = GD, lower = c(0), 
               upper = c(1))
  
  
  aux = ifelse(fittedProb < res$par,TRUE,FALSE)
  aux.1 = which(is.na(z1))
  aux.sub = aux[aux.1]
  z1[aux.1][aux.sub] = 0
  
  stack1 = inla.stack(tag='est1',
                      # - Name of the stack
                      data=list(Y=z1),
                      effects=list(
                        # - The spde Components
                        isetj,
                        # - Dataframe for the linear and non-linear effects as the case may be
                        data.frame(intercept1=1, x1= x,w1=w)),
                      # - The projection matrix A or C in the article
                      A=list(A, 1)
                      # - 1 here represents intercept.
  )
  
  formulae <- Y ~ -1 +intercept1+x1+w1+
                      f(j, model = spde)
  
  sim21 <- inla(formulae,  
               family = c("poisson"),
               data = inla.stack.data(stack1), 
               control.predictor = list(compute = TRUE,
                                        A = inla.stack.A(stack1)), 
               control.fixed = list(expand.factor.strategy = 'inla'),
               control.compute = list(config = TRUE,dic = TRUE,waic=TRUE,cpo=TRUE)
               ,verbose = TRUE) 
  
  
 summary(sim21)
 
 sim1$summary.fixed
 sim20$summary.fixed
 sim21$summary.fixed
 
################
# PLOT
################
boundary <- path

stepsize <- 0.08 

k = 1
nxy <- round(
  c(diff(range( boundary[, 1])), 
    diff(range( boundary[, 2]))) / stepsize)
projgrid <- inla.mesh.projector(
  mesh, xlim = range(boundary[, 1]), 
  ylim = range(boundary[, 2]), dims = nxy)



uhat1 = sim1$summary.random$i$mean
uhat2 = sim20$summary.random$i$mean
uhat3 = sim21$summary.random$j$mean

  xmean1  <- inla.mesh.project(
    projgrid,  (uhat1))
  xmean2  <- inla.mesh.project(
    projgrid,  (uhat2))
  xmean3  <- inla.mesh.project(
    projgrid,  (uhat3))


df <-  expand.grid(x = projgrid$x, y = projgrid$y)


  
  df$mean_s1 <- as.vector(xmean1)
  df$mean_s2 <- as.vector(xmean2)
  df$mean_s3 <- as.vector(xmean3)
  
  ind <- point.in.polygon(
    df$x, df$y,
    boundary[, 1], boundary[, 2]
  )
  
  dff <- df[which(ind == 1), ]

subshp =  st_coordinates(shp_df) %>%as.data.frame()
#####
p_fitJoint <- ggplot(dff, aes(x=x, y=y)) +
  geom_tile(aes(fill = (mean_s1) )) +
  scale_fill_viridis_c(option = "H", direction = 1,name="", guide=guide_colorbar(
    barheight = unit(60, units = "mm"),
    barwidth = unit(1, units = "mm"),
  ))+
  labs(title ="  ",
       y = "",x="") +
  theme_light()+scale_size(guide = guide_legend(direction = "vertical"))+
  geom_point(data=subshp,aes(x=X,y=Y),color="black",size=0.1)
p_fitSeq1 <- ggplot(dff, aes(x=x, y=y)) +
  geom_tile(aes(fill = (mean_s2) )) +
  scale_fill_viridis_c(option = "H", direction = 1,name="", guide=guide_colorbar(
    barheight = unit(60, units = "mm"),
    barwidth = unit(1, units = "mm"),
  ))+
  labs(title =" ",
       y = "",x="") +
  theme_light()+scale_size(guide = guide_legend(direction = "vertical"))+
  geom_point(data=subshp,aes(x=X,y=Y),color="black",size=0.1)
p_fitSeq2 <- ggplot(dff, aes(x=x, y=y)) +
  geom_tile(aes(fill = (mean_s3) )) +
  scale_fill_viridis_c(option = "H", direction = 1,name="", guide=guide_colorbar(
    barheight = unit(60, units = "mm"),
    barwidth = unit(1, units = "mm"),
  ))+
  labs(title =" ",
       y = "",x="") +
  theme_light()+scale_size(guide = guide_legend(direction = "vertical"))+
  geom_point(data=subshp,aes(x=X,y=Y),color="black",size=0.1)


library(patchwork)
(p1)/(p_fitJoint|p_fitSeq2)
