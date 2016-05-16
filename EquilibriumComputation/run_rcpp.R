
library("Rcpp")
library("RcppArmadillo")
library("MASS")
library("rbenchmark")
library("plot3D")
library("pracma")
library("lattice")
library("gridExtra")
library("nleqslv")

setwd("/home/david/Dropbox/Documents/Doctorado/Research/Didris/cppcode/tthetabarra4/UniversityFP/UniversityFP")

rm(list=ls())
set.seed(90)

Sys.setenv("PKG_CXXFLAGS"="-std=c++11 -fopenmp")
Sys.setenv("OMP_NUM_THREADS"=8)
Sys.setenv("PKG_LIBS"="-lgsl -lgslcblas -lnlopt -fopenmp")

sourceCpp("households.cpp")


#=============================================================#
#===                     PARAMETERS                       ====#
#=============================================================#

# 1. Grids
bbmin = 1;
bbmax = 10;
ttmin = 0;
ttmax = 1;
SIZEGRIDb = 200;
SIZEGRIDt = 200

ABARMIN = 0;
ABARMAX = 4;
SIZEABARGRID = 10;

# 2. Parameters of universities
Endowment = 0;
Fixedcost = 100;
aalpha1 = 0.25;
aalpha2 = 0.25;
aalpha3 = 1-aalpha1-aalpha2;


# 3. Parameters of individuals
bbeta = 0.971;
ssigma = 2;
ggamma = 0;

# 4. Parameters of economy
wl = 2;
r = 0.03;
Abar = 1;

# 5. Tax policies
tthetaRthreshold = 0;
bRthreshold = 20;
tax = 0;
Rsubsidy = 0;

ttgrid = seq(ttmin, ttmax, len = SIZEGRIDt);
bbgrid = seq(bbmin, bbmax, len = SIZEGRIDb);

abargrid = seq(ABARMIN, ABARMAX, len = SIZEABARGRID);


w = 4;
zh = 2.5;
zl = 2;
Pl = 0.5;
Ph = 1.1;
tthetaadmissionl = 0;
tthetaadmissionh = 0;

grid = meshgrid(bbgrid, ttgrid)

bpoints = matrix(grid$X, SIZEGRIDt*SIZEGRIDb,1)
tpoints = matrix(grid$Y, SIZEGRIDt*SIZEGRIDb,1)

decision(bpoints[99], tpoints[99], ssigma, bbeta, w, w*(1+zl), w*(1+zh), Pl, Ph, tthetaadmissionl, tthetaadmissionh, ggamma, r, Abar, tthetaRthreshold, bRthreshold, Rsubsidy, tax)

UNCSTUDY(bpoints[99], tpoints[99], ssigma, bbeta, w*(1+zh), Ph, ggamma, r, Abar, tthetaadmissionh, tthetaRthreshold, bRthreshold, Rsubsidy, tax);
CSTUDY(bpoints[99], tpoints[99], ssigma, bbeta, w*(1+zh), Ph, ggamma, r, Abar, tthetaadmissionh, tthetaRthreshold, bRthreshold, Rsubsidy, tax);

decisiones = matrix(0,SIZEGRIDt*SIZEGRIDb, 1)
n = SIZEGRIDt*SIZEGRIDb

for(i in 1:n){
  x = decision(bpoints[i], tpoints[i], ssigma, bbeta, w, w*(1+zl), w*(1+zh), Pl, Ph, tthetaadmissionl, tthetaadmissionh, ggamma, r, Abar, tthetaRthreshold, bRthreshold, Rsubsidy, tax)
  if (x[1] == 1){
    if(bpoints[i] >= Ph + (bbeta*(1+r))^(-1/ssigma)*w*tpoints[i]*(1+zh)){
      decisiones[i,1] = 6
    } else if(bpoints[i] < Ph + (bbeta*(1+r))^(-1/ssigma)*w*tpoints[i]*(1+zh) && bpoints[i] > Pl + (bbeta*(1+r))^(-1/ssigma)*w*tpoints[i]*(1+zl)){
      decisiones[i,1] = 5
    } else{
      decisiones[i,1] = 4
    }
  } else if(x[2] == 1){
    if(bpoints[i] >= Pl + (bbeta*(1+r))^(-1/ssigma)*w*tpoints[i]*(1+zl)){
      decisiones[i,1] = 3
    } else {
      decisiones[i,1] = 2
    }
  } else{
    if(bpoints[i] >= Pl + (bbeta*(1+r))^(-1/ssigma)*w*tpoints[i]*(1+zl)){
      decisiones[i,1] = 1
    }
  }
  if(tpoints[i] >= ((1+r)/w)*((Ph-Pl)/(zh-zl))){
    decisiones[i,1] = decisiones[i,1] + 0.2
  }
  if(tpoints[i] >= ((1+r)/w)*(Pl/(zl-1-r))){
    decisiones[i,1] = decisiones[i,1] + 0.2
  }
}

dat = data.frame(tt = tpoints, bb = bpoints, dd = decisiones)



levelplot(dd ~ tt*bb, data = dat,
          xlab = "Ability", ylab = "Wealth",
          main = "Model",
          col.regions = gray(100:10/100),
          colorkey = FALSE
)



implic = function(x){
  
  theta = x[1]; 
  wealth = x[2];
  
  Phi = (1/(1-ssigma))*((1/(1+(bbeta*(1+r))^(-1/ssigma)*(1+r)))^(1-ssigma))*((bbeta*(1+r))^((ssigma-1)/ssigma)+1)
  
  x[1] = (1/(1-ssigma))*(wealth-Pl+Abar)^(1-ssigma)+(bbeta/(1-ssigma))*(w*theta*(1+zl)-(1+r)*Abar)^(1-ssigma)-Phi*(w*theta*(2+r)+wealth*(1+r))^(1-ssigma)
  
  A1 = (bbeta*(1+zl))^(1/ssigma)
  A2 = (Phi*(1-ssigma)*(2+r))^(1/ssigma)
  
  x[2] = -(w*theta*(2+r)*A1+(1+r)*wealth*A1+(1+r)*Abar*A2-w*theta*(1+zl)*A2)
  
 return(x)
}













ttt = seq(0.05, 0.99, len = 50);
y = seq(0, 0, len = 50);
for(i in 1:50){
  
  func = function(x){
    
    theta = ttt[i]; 
    wealth = x;
    
    Phi = (1/(1-ssigma))*((1/(1+(bbeta*(1+r))^(-1/ssigma)*(1+r)))^(1-ssigma))*((bbeta*(1+r))^((ssigma-1)/ssigma)+1)
    
    x = (1/(1-ssigma))*(wealth-Pl+Abar)^(1-ssigma)+(bbeta/(1-ssigma))*(w*theta*(1+zl)-(1+r)*Abar)^(1-ssigma)-Phi*(w*theta*(2+r)+wealth*(1+r))^(1-ssigma)
    
    return(x)
  }
  
  x = nleqslv(2, func, control=list(btol=.01))
  y[i] = x$x
}

plot(ttt,y)















av = aaverages(zh, zl, bbgrid, ttgrid, ssigma, bbeta, w, Pl, Ph, ggamma, Abar, r, tthetaadmissionl, tthetaadmissionh, tthetaRthreshold, bRthreshold, Rsubsidy, tax)


decision(bbgrid[50], ttgrid[700], ssigma, bbeta, w, w*(1+zl), w*(1+zh), Pl, Ph, tthetaadmissionl, tthetaadmissionh, ggamma, r, Abar, tthetaRthreshold, bRthreshold, Rsubsidy, tax)


for (ib in 1 : 50){
  for (it in 1:700){
    x = decision(bbgrid[ib], ttgrid[it], ssigma, bbeta, w, w*(1+zl), w*(1+zh), Pl, Ph, tthetaadmissionl, tthetaadmissionh, ggamma, r, Abar, tthetaRthreshold, bRthreshold, Rsubsidy, tax)
    print(x[1])
  }
}



x = errZ(aalpha1, aalpha2, 0, 0, zh, zl, bbgrid, ttgrid, ssigma, bbeta, w, Pl, Ph, ggamma, Abar, r, tthetaadmissionl, tthetaadmissionh, tthetaRthreshold, bRthreshold, Rsubsidy, tax)



for (i in 1:50){
  for (j in 1:50){
    x = decision(bbgrid[i], ttgrid[j], ssigma, bbeta, w, w*(1+zl), w*(1+zh), Pl, Ph, tthetaadmissionl, tthetaadmissionh, ggamma, r, Abar, tthetaRthreshold, bRthreshold, Rsubsidy, tax)
    
  }
}









#--------------------------------------------------------#
#        Sample state space - Education decision         #
#--------------------------------------------------------#

# Empirical state space

A = matrix(0, 6, 10)

A[1,] = c(0.0146290381859574, 0.801522408792775, 0.618433417289759, 0.535430028884792, 0.490744371059606, 0.438139354854088, 0.429871578769211, 0.408807757785548, 0.398331582576988, 0.423754992140827)
A[2,] = c(0.0164736216768754, 0.53339545404938, 0.640760676198583, 0.717642449325497, 0.854870174374406, 0.944700158839534, 1.16409425396352, 1.31477196795604, 1.61889360722717, 2.12116881863381)
A[3,] = c(0.00654760308134761, 0.329532531364584, 0.506607859314152, 0.701124789775198, 0.96533763450481, 1.25919151159392, 1.79263508546744, 2.39751277479549, 3.33753375356538, 4.92900965714056)
A[4,] = c(0.000511643733936552, 0.186032893509745, 0.328716907678749, 0.588009382215278, 0.901313684702479, 1.4369214215677, 2.24101828021421, 3.65753724815139, 6.2425060921694, 12.056840506004)
A[5,] = c(0.000465497020671419, 0.141213179702542, 0.316334910237, 0.542294315355845, 0.847713558460114, 1.42606376452557, 2.67558836797225, 3.89520328858822, 7.48565255914089, 16.1196654250318)
A[6,] = c(0.000310972811137805, 0.0922789374309673, 0.198750803125325, 0.406185408893991, 0.777079250796986, 1.25494326548444, 2.23454952688146, 4.17973758403486, 7.87323274434798, 18.1305155222317)

K = matrix(A,60,1)
ES = matrix(seq(1,6,len=6),60,1)
DE = matrix(t(matrix(matrix(seq(1,10,len=10),1,60),10,6)),60,1)

dat = data.frame(x = ES, y = DE, z = K)

elevation.loess = loess(z ~ x*y, data = dat,
                        degree = 2, span = 0.25)

elevation.fit = expand.grid(list(x = seq(1, 6, 1), y = seq(1, 10, 1)))

z = predict(elevation.loess, newdata = elevation.fit)

elevation.fit$Height = as.numeric(z)

plot1 = levelplot(K ~ DE*ES, data = dat,
          xlab = list("Ability", cex = 1), ylab = list("Wealth", cex = 1),
          main = "Empirical",
          col.regions = gray(100:100/100),
          contour = TRUE,
          cuts = 20,
          at = c(0.9,1.5,3,5,10),
          colorkey = FALSE
)


# Model state space
Ph = 2.4
z = 2.535817
tthetaadmission = 3.286504e-06
x = aaverages(z, bbgrid, ttgrid, ssigma, bbeta, wl, Ph, ggamma, Abar, r, Endowment, Fixedcost, aalpha1, aalpha2, tthetaadmission,tthetaRthreshold,bRthreshold, Rsubsidy,tax)

grid = meshgrid(bbgrid, ttgrid)

bpoints = matrix(grid$X, SIZEGRID^2,1)
tpoints = matrix(grid$Y, SIZEGRID^2,1)

decisiones = matrix(0,SIZEGRID^2, 1)
for(i in 1:SIZEGRID^2){
  decisiones[i] = decision(bpoints[i], tpoints[i], ssigma, bbeta, wl, wl*(1+z), Ph, ggamma, r, Abar, tthetaadmission, tthetaRthreshold, bRthreshold, Rsubsidy, tax)
}

dat = data.frame(tt = tpoints, bb = bpoints, dd = decisiones)

plot2 = levelplot(dd ~ tt*bb, data = dat,
               xlab = "Ability", ylab = "Wealth",
               main = "Model",
               col.regions = gray(90:35/100),
               colorkey = FALSE
)

grid.arrange(plot2)
# pdf('/home/david/Dropbox/Documents/Doctorado/Research/Didris/Codigo DZ/model_state_space.pdf', width = 10, height = 5)
# grid.arrange(plot1,plot2, ncol=2, heights=unit(1, "npc"))
# dev.off();

#--------------------------------------------------------#
#   Encuentro policies de unviersidad de equilibrio      #
#--------------------------------------------------------#

subsidios = seq(0, 1, len = 10)
PP = matrix(0,1,length(subsidios))
TT = matrix(0,1,length(subsidios))
ZZ = matrix(0,1,length(subsidios))
TTp = matrix(0,1,length(subsidios))
BB = matrix(0,1,length(subsidios))
NN = matrix(0,1,length(subsidios))


for(ij in 1:length(abargrid)){
  
  # Rsubsidy = subsidios[ij];
  Abar = abargrid[ij]
  
  # Function to estimate surplus as function of tax - subsidy fixed
  surplus = function(x){
    surplus = FiscalSurpluss(Abar, bbgrid, ttgrid, ssigma, bbeta, wl, ggamma, r, Endowment, Fixedcost, aalpha1, aalpha2, tthetaRthreshold, bRthreshold, Rsubsidy, x)
    return(surplus)
  }
  
  # Tax that makes budget balanced in equilibrium
  x = uniroot(surplus, c(0, 0.1), tol = 0.000001)
  tax = x$root
  
  # Recovering equilibrium at this tax
  x = PriceQualitySolver(Abar, bbgrid, ttgrid, ssigma, bbeta, wl, ggamma, r, Endowment, Fixedcost, aalpha1, aalpha2, tthetaRthreshold, bRthreshold, Rsubsidy, tax)
  
  Ph = x[1]; tthetaadmission = x[2];
  
  z = FPsolversquarederror(bbgrid,ttgrid,ssigma,bbeta,wl, Ph,ggamma,Abar,r, Endowment,Fixedcost, aalpha1,aalpha2,tthetaadmission,tthetaRthreshold,bRthreshold,Rsubsidy,tax)
  
  x = aaverages(z, bbgrid, ttgrid, ssigma, bbeta, wl, Ph, ggamma, Abar, r, Endowment, Fixedcost, aalpha1, aalpha2, tthetaadmission,tthetaRthreshold,bRthreshold, Rsubsidy,tax)
  
  grid = meshgrid(bbgrid, ttgrid)
  
  bpoints = matrix(grid$X, SIZEGRID^2,1)
  tpoints = matrix(grid$Y, SIZEGRID^2,1)
  
  decisiones = matrix(0,SIZEGRID^2, 1)
  for(i in 1:SIZEGRID^2){
    decisiones[i] = decision(bpoints[i], tpoints[i], ssigma, bbeta, wl, wl*(1+z), Ph, ggamma, r, Abar, tthetaadmission, tthetaRthreshold, bRthreshold, Rsubsidy, tax)
  }
  
  dat = data.frame(tt = tpoints, bb = bpoints, dd = decisiones)
  
  levelplot(dd ~ tt*bb, data = dat,
            xlab = "\theta", ylab = "b",
            main = "State space",
            col.regions = gray(90:35/100)
  )
  
  
  PP[ij] = Ph
  TT[ij] = tthetaadmission
  ZZ[ij] = z
  TTp[ij] = x[1]
  BB[ij] = x[2]
  NN[ij] = x[3]
}



















