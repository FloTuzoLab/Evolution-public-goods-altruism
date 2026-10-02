##Note: only the first case is explained in detail as further scenarios
##reuse the same functions and algorithms.

library(nleqslv)##importing non linear equation systems solver package
library(viridis)##importing colourblind friendly package
ncol=256
palet<-viridis(ncol)
pal<-list(palet)
setwd("~/")#Setting the directory where to find relevant scripts
source("AP11.Plotting_multiple_images.R")##Importing scripts to make figures with multiple plots

##0.Dichothomy algorithm to find the minimum of a function (without change of concavity)
Dichotomy_solver_0<-function(f,tol){#Finding the zero of a function
  p_c<-1-1e-2
  p_init<-1e-2
  if(f(p_c)>0 && f(p_init)>0){
    return(0)
  }else if(f(p_c)<0 && f(p_init)<0){
    return(1)
  }else{
    p_inter<-(p_c+p_init)/2
    while((p_inter-p_init)>tol || f(p_inter)>tol){
      if(f(p_inter)>0){
        p_c<-p_inter
      }
      else{
        p_init<-p_inter
      }
      p_inter<-(p_c+p_init)/2
    }
  }
  return(p_inter)
}


##1.Rescaled average fitness of genotypes in a local area with frequency x ()
#1.1.Fitness function for cheaters
##Below: group-level relatedness R_GL (but results would also be valid replaceing R_GL by R=1/G+(1-1/G)R_GL); b,c,: benefit and cost of altruism; I: level of investment
w_C<-function(b,c,I,R_GL,x){#x denotes the frequency of altruists in a local area
  num<-(1+b*I*(1-R_GL)*x)
  den<-(1+(b-c)*I*x-b*c*I^2*x*(R_GL+(1-R_GL)*x))#Average fitness in a local area
  return(num/den)##function plugged for cheaters in equation (S.21)
}#Equation 

#1.2.Fitness function for altruists
w_A<-function(b,c,I,R_GL,x){#x denotes the frequency of altruists in a local area
  num<-(1-c*I)*(1+b*I*(R_GL+(1-R_GL)*x))
  den<-(1+(b-c)*I*x-b*c*I^2*x*(R_GL+(1-R_GL)*x))#Average fitness in a local area
  return(num/den)##function plugged for altruists in equation (S.21)
}

##2.Summary of outcomes
#2.0.Additive approximation (I->0, weak selection)

#Set of parameters
R_LT_set<-c(1/4,1/2,3/4)#Equivalent to R_LT, the relatedness due to processes above the competitive constraint
c<-1#Cost c
N_reso=250#Resolution for two plot parameters below
R_GL_set<-seq(0.01,0.99,length=N_reso)#Affinity bias in groups, equals R_GL as we consider large groups
b_set<-seq(1,10,length=N_reso)#Benefit b (which makes b/c a variable that only depends on c)

pA_out<-list()#Frequency of altruists at equilibrium using numerical solution
pA_analytical1<-list()#Frequency of altruists at equilibrium using analytical solution
diff_pA<-list()#Difference between above solutions


##2.0.A.Algorithm to generate the results for along the range of the three parameters (R_GL, R_LT and b)
for(f in 1:length(R_LT_set)){##Different values for R_LT parameter
  R_LT<-R_LT_set[f]##Setting R_LT
  theta<-1/R_LT-1##Conversion to pop. parameter theta
  pA_out[[f]]<-matrix(nrow=N_reso,ncol=N_reso)##Matrix of results for a value of R_LT
  pA_analytical1[[f]]<-matrix(nrow=N_reso,ncol=N_reso)##Matrix of results for a value of R_LT
  diff_pA[[f]]<-matrix(nrow=N_reso,ncol=N_reso)##Matrix of results for a value of R_LT
  for(i in 1:N_reso){##Different values for b, running along the b_set range
    b<-b_set[i]##Setting b as the second parameter current value
    print(c(f,i))##Printing the current state of the algorithm
    for(j in 1:N_reso){##Different values for affinity, running along the R_GL range
      R_GL<-R_GL_set[j]##Setting current value of R_GL
      Rtot<-R_GL+(1-R_GL)*R_LT##Determining total relatedness value
      I<-0.01#value set to a low level of investment (vanishing trade-off)
      DeltaW_C<-function(p){##Function determining the selection difference in favour of cheaters (against the mean population fitness)
        integrand_DWC<-function(x){w_C(b,c,I,R_GL,x)*x^(theta*p-1)*(1-x)^(theta*(1-p))/(1-p)}#Value of cheater's fitness in a local area weighted by how often cheaters see each type of local areas (composition) 
        return(integrate(integrand_DWC,0,1)$value/(beta(theta*p,theta*(1-p)))-1)##Return the difference between cheater mean fitness (averaged over all local areas) and total mean fitness (that equals 1) -> if positive, cheaters have the edge
      }
      pA_analytical1[[f]][i,j]=((R_GL*b-c)*(1+R_LT)-b*c*I*(R_GL+R_LT))/(I*b*c*(1-Rtot))#Analytical approximation for the frequency of altruists at evolutionary steady-state 
      if(pA_analytical1[[f]][i,j]>1){##Setting an upper bound of 1 for the frequency of altruists in the population (where altruists go to fixation)
        pA_analytical1[[f]][i,j]=1
      }
      if(R_GL*b-c<0){#Adding value for Hamilton's rule to speed up algorithm
        pA_out[[f]][i,j]<--0.1#When Hamilton's rule is not met: value out of the range for the frequency of altruists to set it aside as a grey area
        pA_analytical1[[f]][i,j]<--0.1#When Hamilton's rule is not met: value out of the range for the frequency of altruists to set it aside as a grey area
      }else{
        pA_out[[f]][i,j]<-Dichotomy_solver_0(DeltaW_C,1e-3)#Solving for the evolutionary steady-state (determining the frequency of altruists at equilibrium when cheaters and altruists have same fitness)
      }
      if(pA_out[[f]][i,j]>0 || pA_analytical1[[f]][i,j]>0){##Calculating the difference between solutions where altruism invades (when rare)
        diff_pA[[f]][i,j]<-abs(pA_out[[f]][i,j]-pA_analytical1[[f]][i,j])/abs(pA_out[[f]][i,j])
      }
    }
  }
}

##2.0.B.Figures with the results
R_var<-c(0,1,0.2)##Values of assortment bias to write on the axis
b_var<-c(1,10,3)##Values of b (b/c) to write on the axis
addtxt<-list(l=-0.05,h=1.1,txt=c("S6A.","S6B.","S6C."),srt = 0,font=2,col="black")##Letters to identify plots

##Figure S6: Plotting the frequency of altruists at steady-state (gradient of colours) using numerical solutions
multiplePlot("","","","",c(""),c("","","")
             ,ncol=800,b_var,R_var,pA_out,
             abs="b/c",ord=expression(paste("Assortment bias ",R[GL])),scale=c(0,1),palette=pal,cextext=1.5,TEXT_to_Add=addtxt,sub=c("R_LT=0","R_LT=1/4","R_LT=1/2","R_LT=3/4"),
             image=TRUE,pcex=1,subcex=2,labcex=1.75,axcex=0.8,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")

##Plotting the frequency of altruists at steady-state (gradient of colours) using analytical solutions
multiplePlot("","","","",c(""),c("","","")
             ,ncol=800,b_var,R_var,pA_analytical1,
             abs="b/c",ord=expression(paste("Assortment bias ",R[GL])),scale=c(0,1),palette=pal,cextext=1.5,TEXT_to_Add=addtxt,sub=c("R_LT=0","R_LT=1/4","R_LT=1/2","R_LT=3/4"),
             image=TRUE,pcex=1,subcex=2,labcex=1.75,axcex=0.8,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")

##Contrasting the results of the two solutions in the relevant part of the parameter space (dark blue implies good accuracy of analytical approximation)
multiplePlot("","","","",c(""),c("","","")
             ,ncol=800,b_var,R_var,diff_pA,
             abs="b/c",ord=expression(paste("Assortment bias ",R[GL])),scale=c(0,1),palette=pal,cextext=1.5,sub=c("R_LT=0","R_LT=1/4","R_LT=1/2","R_LT=3/4"),
             image=TRUE,pcex=1,subcex=2,labcex=1.75,axcex=0.8,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")

#2.1.Case of global competition (R_LT->0)

#Set of parameters
R_LT_set<-c(0.01)##Setting R_LT to a low value (0.01<<1)
I_set<-c(1/5,2/5,3/5)##Varying the level of investment

for(f in 1:length(I_set)){
  R_LT<-R_LT_set[1]##Only change compared to previous script: setting R_LT to a low level to mimic global competition
  theta<-1/R_LT-1#Corresponding value of Theta used in the Beta distribution
  I<-I_set[f]#Setting the current value of the investment level
  pA_out[[f]]<-matrix(nrow=N_reso,ncol=N_reso)#Frequency of altruists at evo equilibrium (numerical solution)
  pA_analytical1[[f]]<-matrix(nrow=N_reso,ncol=N_reso)#Frequency of altruists at evo equilibrium (analytical approximation)
  diff_pA[[f]]<-matrix(nrow=N_reso,ncol=N_reso)#Difference between the two solutions
  for(i in 1:N_reso){
    b<-b_set[i]
    print(c(f,i))
    for(j in 1:N_reso){
      R_GL<-R_GL_set[j]
      Rtot<-R_GL+(1-R_GL)*R_LT
      DeltaW_C<-function(p){
        integrand_DWC<-function(x){w_C(b,c,I,R_GL,x)*x^(theta*p-1)*(1-x)^(theta*(1-p))/(1-p)}
        return(integrate(integrand_DWC,0,1)$value/(beta(theta*p,theta*(1-p)))-1)
      }
      pA_analytical1[[f]][i,j]=((R_GL*b-c)*(1+R_LT)-b*c*I*(R_GL+R_LT))/(I*b*c*(1-Rtot))
      if(pA_analytical1[[f]][i,j]>1){
        pA_analytical1[[f]][i,j]=1
      }
      if(R_GL*b-c<1e-2){
        pA_out[[f]][i,j]<--0.1
        pA_analytical1[[f]][i,j]<--0.1
      }else{
        pA_out[[f]][i,j]<-Dichotomy_solver_0(DeltaW_C,1e-3)
      }
      if(pA_out[[f]][i,j]>0 || pA_analytical1[[f]][i,j]>0){
        diff_pA[[f]][i,j]<-abs(pA_out[[f]][i,j]-pA_analytical1[[f]][i,j])/abs(pA_out[[f]][i,j])
      }
    }
  }
}

addtxt<-list(l=-0.05,h=1.1,txt=c("S7A.","S7B.","S7C."),srt = 0,font=2,col="black")##Letters to identify plots
##Figure S7--case of local competition: Plotting the frequency of altruists at steady-state (gradient of colours) using numerical solutions
multiplePlot("","","","",c(""),c("","","")
             ,ncol=800,b_var,R_var,pA_out,
             abs="b/c",ord=expression(paste("Assortment bias ",R[GL])),scale=c(0,1),palette=pal,cextext=1.5,TEXT_to_Add=addtxt,sub=c("I=1/5","I=2/5","I=3/5"),
             image=TRUE,pcex=1,subcex=2,labcex=1.75,axcex=0.8,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")

##Plotting the frequency of altruists at steady-state (gradient of colours) using analytical solutions
multiplePlot("","","","",c(""),c("","","")
             ,ncol=800,b_var,R_var,pA_analytical1,
             abs="b/c",ord=expression(paste("Assortment bias ",R[GL])),scale=c(0,1),palette=pal,cextext=1.5,TEXT_to_Add=addtxt,sub=c("I=1/5","I=2/5","I=3/5"),
             image=TRUE,pcex=1,subcex=2,labcex=1.75,axcex=0.8,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")

##Contrasting the results of the two solutions in the relevant part of the parameter space (dark blue implies good accuracy of analytical approximation)
multiplePlot("","","","",c(""),c("","","")
             ,ncol=800,b_var,R_var,diff_pA,
             abs="b/c",ord=expression(paste("Assortment bias ",R[GL])),scale=c(0,1),palette=pal,cextext=1.5,sub=c("I=1/5","I=2/5","I=3/5"),
             image=TRUE,pcex=1,subcex=2,labcex=1.75,axcex=0.8,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")


#2.2.A.Trade-off scenario: influence of local competition (I & R_LT variable, starting with I=1/5)

#Set of parameters
R_LT_set<-c(1/4,1/2,3/4)##Varying R_LT=R_LT
I_set<-c(1/5,2/5,3/5)##Varying the investment level of altruism

for(f in 1:length(R_LT_set)){
  R_LT<-R_LT_set[f]##Varying the parameters R_LT
  theta<-1/R_LT-1
  pA_out[[f]]<-matrix(nrow=N_reso,ncol=N_reso)
  pA_analytical1[[f]]<-matrix(nrow=N_reso,ncol=N_reso)
  diff_pA[[f]]<-matrix(nrow=N_reso,ncol=N_reso)
  for(i in 1:N_reso){
    b<-b_set[i]
    print(c(f,i))
    for(j in 1:N_reso){
      R_GL<-R_GL_set[j]
      Rtot<-R_GL+(1-R_GL)*R_LT
      I<-I_set[1]#fixed value for I to I=1/5
      DeltaW_C<-function(p){
        integrand_DWC<-function(x){w_C(b,c,I,R_GL,x)*x^(theta*p-1)*(1-x)^(theta*(1-p))/(1-p)}
        return(integrate(integrand_DWC,0,1)$value/(beta(theta*p,theta*(1-p)))-1)
      }
      pA_analytical1[[f]][i,j]=((R_GL*b-c)*(1+R_LT)-b*c*I*(R_GL+R_LT))/(I*b*c*(1-Rtot))
      if(pA_analytical1[[f]][i,j]>1){
        pA_analytical1[[f]][i,j]=1
      }
      if(R_GL*b-c<1e-2){
        pA_out[[f]][i,j]<--0.1
        pA_analytical1[[f]][i,j]<--0.1
      }else{
        pA_out[[f]][i,j]<-Dichotomy_solver_0(DeltaW_C,1e-3)
      }
      if(pA_out[[f]][i,j]>0 || pA_analytical1[[f]][i,j]>0){
        diff_pA[[f]][i,j]<-abs(pA_out[[f]][i,j]-pA_analytical1[[f]][i,j])/abs(pA_out[[f]][i,j])
      }
    }
  }
}

addtxt<-list(l=-0.05,h=1.1,txt=c("S8A.","S8B.","S8C."),srt = 0,font=2,col="black")##Letters to identify plots
##Figure S8/ABC--case of local competition: Plotting the frequency of altruists at steady-state (gradient of colours) using numerical solutions
multiplePlot("","","","",c(""),c("","","")
             ,ncol=800,b_var,R_var,pA_out,
             abs="b/c",ord=expression(paste("Assortment bias ",R[GL])),scale=c(0,1),palette=pal,cextext=1.5,TEXT_to_Add=addtxt,sub=c(expression(paste(R[LT]," = 1/4")),expression(paste(R[LT]," = 1/2")),expression(paste(R[LT]," = 3/4"))),
             image=TRUE,pcex=1,subcex=2,labcex=1.75,axcex=0.8,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")

##Plotting the frequency of altruists at steady-state (gradient of colours) using analytical solutions
multiplePlot("","","","",c(""),c("","","")
             ,ncol=800,b_var,R_var,pA_analytical1,
             abs="b/c",ord=expression(paste("Assortment bias ",R[GL])),scale=c(0,1),palette=pal,cextext=1.5,TEXT_to_Add=addtxt,sub=c(expression(paste(R[LT]," = 1/4")),expression(paste(R[LT]," = 1/2")),expression(paste(R[LT]," = 3/4"))),
             image=TRUE,pcex=1,subcex=2,labcex=1.75,axcex=0.8,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")

##Contrasting the results of the two solutions in the relevant part of the parameter space (dark blue implies good accuracy of analytical approximation)
multiplePlot("","","","",c(""),c("","","")
             ,ncol=800,b_var,R_var,diff_pA,
             abs="b/c",ord=expression(paste("Assortment bias ",R[GL])),scale=c(0,1),palette=pal,cextext=1.5,sub=c(expression(paste(R[LT]," = 1/4")),expression(paste(R[LT]," = 1/2")),expression(paste(R[LT]," = 3/4"))),
             image=TRUE,pcex=1,subcex=2,labcex=1.75,axcex=0.8,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")


#2.2.B.Trade-off scenario: influence of local competition (I & R_LT variable, starting with I=2/5)

for(f in 1:length(R_LT_set)){
  R_LT<-R_LT_set[f]##Varying the parameters R_LT
  theta<-1/R_LT-1
  pA_out[[f]]<-matrix(nrow=N_reso,ncol=N_reso)
  pA_analytical1[[f]]<-matrix(nrow=N_reso,ncol=N_reso)
  diff_pA[[f]]<-matrix(nrow=N_reso,ncol=N_reso)
  for(i in 1:N_reso){
    b<-b_set[i]
    print(c(f,i))
    for(j in 1:N_reso){
      R_GL<-R_GL_set[j]
      Rtot<-R_GL+(1-R_GL)*R_LT
      I<-I_set[2]#fixed value to I=2/5
      DeltaW_C<-function(p){
        integrand_DWC<-function(x){w_C(b,c,I,R_GL,x)*x^(theta*p-1)*(1-x)^(theta*(1-p))/(1-p)}
        return(integrate(integrand_DWC,0,1)$value/(beta(theta*p,theta*(1-p)))-1)
      }
      pA_analytical1[[f]][i,j]=((R_GL*b-c)*(1+R_LT)-b*c*I*(R_GL+R_LT))/(I*b*c*(1-Rtot))
      if(pA_analytical1[[f]][i,j]>1){
        pA_analytical1[[f]][i,j]=1
      }
      if(R_GL*b-c<1e-2){
        pA_out[[f]][i,j]<--0.1
        pA_analytical1[[f]][i,j]<--0.1
      }else{
        pA_out[[f]][i,j]<-Dichotomy_solver_0(DeltaW_C,1e-3)
      }
      if(pA_out[[f]][i,j]>0 || pA_analytical1[[f]][i,j]>0){
        diff_pA[[f]][i,j]<-abs(pA_out[[f]][i,j]-pA_analytical1[[f]][i,j])/abs(pA_out[[f]][i,j])
      }
    }
  }
}

addtxt<-list(l=-0.05,h=1.1,txt=c("S8D.","S8E.","S8F."),srt = 0,font=2,col="black")##Letters to identify plots
##Figure S8/DEF--case of local competition: Plotting the frequency of altruists at steady-state (gradient of colours) using numerical solutions
multiplePlot("","","","",c(""),c("","","")
             ,ncol=800,b_var,R_var,pA_out,
             abs="b/c",ord=expression(paste("Assortment bias ",R[GL])),scale=c(0,1),palette=pal,cextext=1.5,TEXT_to_Add=addtxt,sub=c(expression(paste(R[LT]," = 1/4")),expression(paste(R[LT]," = 1/2")),expression(paste(R[LT]," = 3/4"))),
             image=TRUE,pcex=1,subcex=2,labcex=1.75,axcex=0.8,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")

##Plotting the frequency of altruists at steady-state (gradient of colours) using analytical solutions
multiplePlot("","","","",c(""),c("","","")
             ,ncol=800,b_var,R_var,pA_analytical1,
             abs="b/c",ord=expression(paste("Assortment bias ",R[GL])),scale=c(0,1),palette=pal,cextext=1.5,TEXT_to_Add=addtxt,sub=c(expression(paste(R[LT]," = 1/4")),expression(paste(R[LT]," = 1/2")),expression(paste(R[LT]," = 3/4"))),
             image=TRUE,pcex=1,subcex=2,labcex=1.75,axcex=0.8,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")

##Contrasting the results of the two solutions in the relevant part of the parameter space (dark blue implies good accuracy of analytical approximation)
multiplePlot("","","","",c(""),c("","","")
             ,ncol=800,b_var,R_var,diff_pA,
             abs="b/c",ord=expression(paste("Assortment bias ",R[GL])),scale=c(0,1),palette=pal,cextext=1.5,sub=c(expression(paste(R[LT]," = 1/4")),expression(paste(R[LT]," = 1/2")),expression(paste(R[LT]," = 3/4"))),
             image=TRUE,pcex=1,subcex=2,labcex=1.75,axcex=0.8,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")


#2.2.C.Trade-off scenario: influence of local competition (I & R_LT variable, starting with I=3/5)

for(f in 1:length(R_LT_set)){
  R_LT<-R_LT_set[f]
  theta<-1/R_LT-1
  pA_out[[f]]<-matrix(nrow=N_reso,ncol=N_reso)
  pA_analytical1[[f]]<-matrix(nrow=N_reso,ncol=N_reso)
  diff_pA[[f]]<-matrix(nrow=N_reso,ncol=N_reso)
  for(i in 1:N_reso){
    b<-b_set[i]
    print(c(f,i))
    for(j in 1:N_reso){
      R_GL<-R_GL_set[j]
      Rtot<-R_GL+(1-R_GL)*R_LT
      I<-I_set[3]#fixed value I=3/5
      DeltaW_C<-function(p){
        integrand_DWC<-function(x){w_C(b,c,I,R_GL,x)*x^(theta*p-1)*(1-x)^(theta*(1-p))/(1-p)}
        return(integrate(integrand_DWC,0,1)$value/(beta(theta*p,theta*(1-p)))-1)
      }
      pA_analytical1[[f]][i,j]=((R_GL*b-c)*(1+R_LT)-b*c*I*(R_GL+R_LT))/(I*b*c*(1-Rtot))
      if(pA_analytical1[[f]][i,j]>1){
        pA_analytical1[[f]][i,j]=1
      }
      if(R_GL*b-c<1e-2){
        pA_out[[f]][i,j]<--0.1
        pA_analytical1[[f]][i,j]<--0.1
      }else{
        pA_out[[f]][i,j]<-Dichotomy_solver_0(DeltaW_C,1e-3)
      }
      if(pA_out[[f]][i,j]>0 || pA_analytical1[[f]][i,j]>0){
        diff_pA[[f]][i,j]<-abs(pA_out[[f]][i,j]-pA_analytical1[[f]][i,j])/abs(pA_out[[f]][i,j])
      }
    }
  }
}

addtxt<-list(l=-0.05,h=1.1,txt=c("S8G.","S8H.","S8I."),srt = 0,font=2,col="black")##Letters to identify plots
##Figure S8/GHI--case of local competition: Plotting the frequency of altruists at steady-state (gradient of colours) using numerical solutions
multiplePlot("","","","",c(""),c("","","")
             ,ncol=800,b_var,R_var,pA_out,
             abs="b/c",ord=expression(paste("Assortment bias ",R[GL])),scale=c(0,1),palette=pal,cextext=1.5,TEXT_to_Add=addtxt,sub=c(expression(paste(R[LT]," = 1/4")),expression(paste(R[LT]," = 1/2")),expression(paste(R[LT]," = 3/4"))),
             image=TRUE,pcex=1,subcex=2,labcex=1.75,axcex=0.8,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")

##Plotting the frequency of altruists at steady-state (gradient of colours) using analytical solutions
multiplePlot("","","","",c(""),c("","","")
             ,ncol=800,b_var,R_var,pA_analytical1,
             abs="b/c",ord=expression(paste("Assortment bias ",R[GL])),scale=c(0,1),palette=pal,cextext=1.5,TEXT_to_Add=addtxt,sub=c(expression(paste(R[LT]," = 1/4")),expression(paste(R[LT]," = 1/2")),expression(paste(R[LT]," = 3/4"))),
             image=TRUE,pcex=1,subcex=2,labcex=1.75,axcex=0.8,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")

##Contrasting the results of the two solutions in the relevant part of the parameter space (dark blue implies good accuracy of analytical approximation)
multiplePlot("","","","",c(""),c("","","")##Grey area denotes difference higher than 10%
             ,ncol=800,b_var,R_var,diff_pA,
             abs="b/c",ord=expression(paste("Assortment bias ",R[GL])),scale=c(0,0.1),palette=pal,cextext=1.5,sub=c(expression(paste(R[LT]," = 1/4")),expression(paste(R[LT]," = 1/2")),expression(paste(R[LT]," = 3/4"))),
             image=TRUE,pcex=1,subcex=2,labcex=1.75,axcex=0.8,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")


#3.3.A.Boundary for monomorphism approximation: numerical and analytical results at the analytical threshold

for(f in 1:length(R_LT_set)){
  R_LT<-R_LT_set[f]
  theta<-1/R_LT-1
  pA_out[[f]]<-matrix(nrow=N_reso,ncol=N_reso)
  pA_analytical1[[f]]<-matrix(nrow=N_reso,ncol=N_reso)
  diff_pA[[f]]<-matrix(nrow=N_reso,ncol=N_reso)
  for(i in 1:N_reso){
    b<-b_set[i]
    print(c(f,i))
    for(j in 1:N_reso){
      R_GL<-R_GL_set[j]
      Rtot<-R_GL+(1-R_GL)*R_LT
      I<-(R_GL*b-c)*(1+R_LT)/(b*c)/((1+R_GL*R_LT))##Investment level for monomorphic altruism: analytical solution
      DeltaW_C<-function(p){
        integrand_DWC<-function(x){w_C(b,c,I,R_GL,x)*x^(theta*p-1)*(1-x)^(theta*(1-p))/(1-p)}
        return(integrate(integrand_DWC,0,1)$value/(beta(theta*p,theta*(1-p)))-1)
      }
      pA_analytical1[[f]][i,j]=((R_GL*b-c)*(1+R_LT)-b*c*I*(R_GL+R_LT))/(I*b*c*(1-Rtot))
      if(pA_analytical1[[f]][i,j]>1){
        pA_analytical1[[f]][i,j]=1
      }
      if(R_GL*b-c<1e-2){
        pA_out[[f]][i,j]<--0.1
        pA_analytical1[[f]][i,j]<--0.1
      }else{
        pA_out[[f]][i,j]<-Dichotomy_solver_0(DeltaW_C,1e-3)
      }
      if(pA_out[[f]][i,j]>0 || pA_analytical1[[f]][i,j]>0){
        diff_pA[[f]][i,j]<-abs(pA_out[[f]][i,j]-pA_analytical1[[f]][i,j])/abs(pA_out[[f]][i,j])
      }
    }
  }
}

##Plotting the frequency of altruists at steady-state (gradient of colours) using numerical solutions
multiplePlot("","","","",c(""),c("","","")
             ,ncol=800,b_var,R_var,pA_out,
             abs="b/c",ord=expression(paste("Assortment bias ",R[GL])),scale=c(0,1),palette=pal,cextext=1.5,sub=c(expression(paste(R[LT]," = 1/4")),expression(paste(R[LT]," = 1/2")),expression(paste(R[LT]," = 3/4"))),
             image=TRUE,pcex=1,subcex=2,labcex=1.75,axcex=0.8,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")

##Plotting the frequency of altruists at steady-state (gradient of colours) using analytical solutions
multiplePlot("","","","",c(""),c("","","")
             ,ncol=800,b_var,R_var,pA_analytical1,
             abs="b/c",ord=expression(paste("Assortment bias ",R[GL])),scale=c(0,1),palette=pal,cextext=1.5,sub=c(expression(paste(R[LT]," = 1/4")),expression(paste(R[LT]," = 1/2")),expression(paste(R[LT]," = 3/4"))),
             image=TRUE,pcex=1,subcex=2,labcex=1.75,axcex=0.8,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")


#3.3.B.Boundary for polymorphism approximation: numerical and analytical results at the analytical threshold

for(f in 1:length(R_LT_set)){
  R_LT<-R_LT_set[f]
  theta<-1/R_LT-1
  pA_out[[f]]<-matrix(nrow=N_reso,ncol=N_reso)
  pA_analytical1[[f]]<-matrix(nrow=N_reso,ncol=N_reso)
  diff_pA[[f]]<-matrix(nrow=N_reso,ncol=N_reso)
  for(i in 1:N_reso){
    b<-b_set[i]
    print(c(f,i))
    for(j in 1:N_reso){
      R_GL<-R_GL_set[j]
      Rtot<-R_GL+(1-R_GL)*R_LT
      I<-(R_GL*b-c)*(1+R_LT)/(b*c*(R_GL+R_LT))##Investment level for polymorphic altruism: analytical solution
      DeltaW_C<-function(p){
        integrand_DWC<-function(x){w_C(b,c,I,R_GL,x)*x^(theta*p-1)*(1-x)^(theta*(1-p))/(1-p)}
        return(integrate(integrand_DWC,0,1)$value/(beta(theta*p,theta*(1-p)))-1)
      }
      pA_analytical1[[f]][i,j]=((R_GL*b-c)*(1+R_LT)-b*c*I*(R_GL+R_LT))/(I*b*c*(1-Rtot))
      if(pA_analytical1[[f]][i,j]>1){
        pA_analytical1[[f]][i,j]=1
      }
      if(R_GL*b-c<1e-2){
        pA_out[[f]][i,j]<--0.1
        pA_analytical1[[f]][i,j]<--0.1
      }else{
        pA_out[[f]][i,j]<-Dichotomy_solver_0(DeltaW_C,1e-3)
        if(pA_out[[f]][i,j]>-0.02 && pA_out[[f]][i,j]<0.01){##Setting proportion of altruists at equilibrium at 0 when the algorithm lead to slightly lower than 0 values
          pA_out[[f]][i,j]=0
        }
      }
      if(pA_out[[f]][i,j]>0 || pA_analytical1[[f]][i,j]>-0.01){
        diff_pA[[f]][i,j]<-abs(pA_out[[f]][i,j]-pA_analytical1[[f]][i,j])/abs(pA_out[[f]][i,j])
      }else{
        diff_pA[[f]][i,j]<-0
      }
    }
  }
}

##Case of local competition: Plotting the frequency of altruists at steady-state (gradient of colours) using numerical solutions
multiplePlot("","","","",c(""),c("","","")
             ,ncol=800,b_var,R_var,pA_out,lev=c(0.1,0.3,0.5),##difference more important with greater local competition and b/c values, as expected from the exact numerical series of the additive case
             abs="b/c",ord=expression(paste("Assortment bias ",R[GL])),scale=c(-0.01,1),palette=pal,cextext=1.5,sub=c(expression(paste(R[LT]," = 1/4")),expression(paste(R[LT]," = 1/2")),expression(paste(R[LT]," = 3/4"))),
             image=TRUE,pcex=1,subcex=2,labcex=1.75,axcex=1.25,globcex=0.65,legcex=1,contourlab=TRUE,meth="edge",contcex=1,colorkey="COMMON")

##Plotting the frequency of altruists at steady-state (gradient of colours) using analytical solutions
multiplePlot("","","","",c(""),c("","","")##logically returns no altruism, even in a polymorphism
             ,ncol=800,b_var,R_var,pA_analytical1,
             abs="b/c",ord=expression(paste("Assortment bias ",R[GL])),scale=c(0,1),palette=pal,cextext=1.5,sub=c(expression(paste(R[LT]," = 1/4")),expression(paste(R[LT]," = 1/2")),expression(paste(R[LT]," = 3/4"))),
             image=TRUE,pcex=1,subcex=2,labcex=1.75,axcex=0.8,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")

##Calculating the average frequency of altruists (while expectation is 0)
mean(pA_out[[1]][pA_out[[1]]>-0.001])
mean(pA_out[[2]][pA_out[[2]]>-0.001])
mean(pA_out[[3]][pA_out[[3]]>-0.001])

#3.3.C.Boundary for monomorphism approximation
##Here, the aim is to determine the difference between the analytical and the numerical solutions for both thresholds

#Decreasing resolution to speeed up resolution
N_reso=100#Resolution for two plot parameters below
R_GL_set<-seq(0.01,0.99,length=N_reso)#Affinity bias in groups, equals R_GL as we consider large groups
b_set<-seq(1,10,length=N_reso)#Benefit b (which makes b/c a variable that only depends on c)

pA_out<-list()#Frequency of altruists at equilibrium using numerical solution
pA_analytical1<-list()#Frequency of altruists at equilibrium using analytical solution
diff_I_fix<-list()#Difference between above solutions
diff_I_poly<-list()#Difference between above solutions

for(f in 1:length(R_LT_set)){
  R_LT<-R_LT_set[f]
  theta<-1/R_LT-1
  pA_out[[f]]<-matrix(nrow=N_reso,ncol=N_reso)
  pA_analytical1[[f]]<-matrix(nrow=N_reso,ncol=N_reso)
  diff_I_fix[[f]]<-matrix(nrow=N_reso,ncol=N_reso)
  for(i in 1:N_reso){
    b<-b_set[i]
    print(c(f,i))
    for(j in 1:N_reso){
      pA_out[[f]][i,j]=1
      R_GL<-R_GL_set[j]
      Rtot<-R_GL+(1-R_GL)*R_LT
      I_thr<-(R_GL*b-c)*(1+R_LT)/(b*c)/((1+R_GL*R_LT))##Investment level for monomorphic altruism
      compteur<-0
      while(pA_out[[f]][i,j]>0.99){#Finding the threshold for the numerical solution with precision 1e-3 on the investment level
        I<-I_thr*(1+compteur/1000)
        DeltaW_C<-function(p){
          integrand_DWC<-function(x){w_C(b,c,I,R_GL,x)*x^(theta*p-1)*(1-x)^(theta*(1-p))/(1-p)}
          return(integrate(integrand_DWC,0,1)$value/(beta(theta*p,theta*(1-p)))-1)
        }
        pA_analytical1[[f]][i,j]=((R_GL*b-c)*(1+R_LT)-b*c*I*(R_GL+R_LT))/(I*b*c*(1-Rtot))
        if(pA_analytical1[[f]][i,j]>1){
          pA_analytical1[[f]][i,j]=1
        }
        if(R_GL*b-c<1e-2){
          pA_out[[f]][i,j]<--0.1
          pA_analytical1[[f]][i,j]<--0.1
        }else{
          pA_out[[f]][i,j]<-Dichotomy_solver_0(DeltaW_C,1e-3)
        }
        compteur<-compteur+1
      }
      diff_I_fix[[f]][i,j]<-(compteur/1000)
    }
  }
}

ncol=256
palet<-magma(ncol)
pal<-list(palet)
addtxt<-list(l=-0.05,h=1.1,txt=c("S9A.","S9B.","S9C."),srt = 0,font=2,col="black")##Letters to identify plots

#Figure S9/ABC:Plotting the difference between the numerical and the analytical investment threshold for altruism fixation/monomorphism
multiplePlot("","","","",c(""),c("","","")
             ,ncol=800,b_var,R_var,diff_I_fix,
             abs="b/c",ord=expression(paste("Assortment bias ",R[GL])),lev=c(0.02,0.04,0.06,0.08),scale=c(0,0.16),palette=pal,cextext=1.5,TEXT_to_Add=addtxt,sub=c(expression(paste(R[LT]," = 1/4")),expression(paste(R[LT]," = 1/2")),expression(paste(R[LT]," = 3/4"))),
             image=TRUE,pcex=1,subcex=2,labcex=2.75,axcex=1,globcex=0.55,legcex=1,contourlab=TRUE,meth="edge",contcex=1.3,colorkey="COMMON")


#3.3.D.Boundary for polymorphism approximation

for(f in 1:length(R_LT_set)){
  R_LT<-R_LT_set[f]
  theta<-1/R_LT-1
  pA_out[[f]]<-matrix(nrow=N_reso,ncol=N_reso)
  pA_analytical1[[f]]<-matrix(nrow=N_reso,ncol=N_reso)
  diff_I_poly[[f]]<-matrix(nrow=N_reso,ncol=N_reso)
  for(i in 1:N_reso){
    b<-b_set[i]
    print(c(f,i))
    for(j in 1:N_reso){
      pA_out[[f]][i,j]=1
      R_GL<-R_GL_set[j]
      Rtot<-R_GL+(1-R_GL)*R_LT
      I_thr<-(R_GL*b-c)*(1+R_LT)/(b*c*(R_GL+R_LT))##Investment level for monomorphic altruism
      compteur<-0
      while(pA_out[[f]][i,j]>0.01){
        I<-I_thr*(1+compteur/1000)
        DeltaW_C<-function(p){
          integrand_DWC<-function(x){w_C(b,c,I,R_GL,x)*x^(theta*p-1)*(1-x)^(theta*(1-p))/(1-p)}
          return(integrate(integrand_DWC,0,1)$value/(beta(theta*p,theta*(1-p)))-1)
        }
        pA_analytical1[[f]][i,j]=((R_GL*b-c)*(1+R_LT)-b*c*I*(R_GL+R_LT))/(I*b*c*(1-Rtot))
        if(pA_analytical1[[f]][i,j]>1){
          pA_analytical1[[f]][i,j]=1
        }
        if(R_GL*b-c<1e-2){
          pA_out[[f]][i,j]<--0.1
          pA_analytical1[[f]][i,j]<--0.1
        }else{
          pA_out[[f]][i,j]<-Dichotomy_solver_0(DeltaW_C,1e-3)
        }
        compteur<-compteur+1
      }
      diff_I_poly[[f]][i,j]<-(compteur/1000)
    }
  }
}

addtxt<-list(l=-0.05,h=1.1,txt=c("S9D.","S9E.","S9F."),srt = 0,font=2,col="black")##Letters to identify plots
##FigS9/DEF:Plotting the difference between the numerical and the analytical investment threshold for cheater-altruism polymorphism
multiplePlot("","","","",c(""),c("","","")
             ,ncol=800,b_var,R_var,diff_I_poly,##difference more important with greater local competition and b/c values, as expected from the exact numerical series of the additive case
             abs="b/c",ord=expression(paste("Assortment bias ",R[GL])),lev=c(0.02,0.04,0.06,0.08,0.1),scale=c(0,0.16),palette=pal,cextext=1.5,TEXT_to_Add=addtxt,sub=c(expression(paste(R[LT]," = 1/4")),expression(paste(R[LT]," = 1/2")),expression(paste(R[LT]," = 3/4"))),
             image=TRUE,pcex=1,subcex=2,labcex=2.75,axcex=1.75,globcex=0.55,legcex=1,contourlab=TRUE,meth="edge",contcex=1.3,colorkey="COMMON")

