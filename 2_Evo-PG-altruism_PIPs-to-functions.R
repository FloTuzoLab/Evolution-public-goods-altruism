setwd("~/")#Setting the directory where to find relevant scripts
source("AP11.Plotting_multiple_images.R")##Importing scripts to make figures with multiple plots

#0.Definition of fitness functions

DeltaW<-function(I1,I2,b,c,R){#Invasion fitness as a function of investment levels - Multiplicative linear trade-off
  Wres<-(1-c*I1)*(1+b*I1)#Resident fitness
  Wmut<-(1-c*I2)*(1+b*R*I2+b*(1-R)*I1)#Mutant fitness
  return(Wmut-Wres)#'Invasion fitness' measure used to determine the sign of the invasion exponent as S(mut,res)=Wmut/Wres, with invasion happening when log S>0 <-> S>1
}

DeltaW_supra_cost_F2<-function(I1,I2,b,c,R){#Invasion fitness as a function of investment levels - Case 1 of high order accelerating costs (Fig. S11 A.-D.)
  Wres<-(1-c*(I1)^4.2-2*c*(I1)^7.5)*(1+b*I1)
  Wmut<-(1-c*(I2)^4.2-2*c*(I2)^7.5)*(1+b*R*I2+b*(1-R)*I1)
  return(Wmut-Wres)
}

DeltaW_supra_cost_F3<-function(I1,I2,b,c,R){#Invasion fitness as a function of investment levels - Case 2 of high order accelerating costs: truncated exponential costs (Fig. S11 E.-H.)
  Wres<-(1-c*I1-c/2*(I1)^2-c/6*(I1)^3)*(1+b*I1)
  Wmut<-(1-c*I2-c/2*(I2)^2-c/6*(I2)^3)*(1+b*R*I2+b*(1-R)*I1)
  return(Wmut-Wres)
}

DeltaW_sublinear_cost<-function(I1,I2,b,c,R){#Invasion fitness as a function of investment levels - Case of concave cost function (Fig. S12)
  Wres<-(1-c/2*(I1)^0.5-c/2*(I1)^1.5)*(1+b*I1)
  Wmut<-(1-c/2*(I2)^0.5-c/2*(I2)^1.5)*(1+b*R*I2+b*(1-R)*I1)
  return(Wmut-Wres)
}

DeltaW_supralinear_benefit5<-function(I1,I2,b,c,R){#Invasion fitness as a function of investment levels - Case of accelerating benefit functions (Fig.S13 A.-H.)
  Wres<-(1-c*I1)*(1+10*b*I1^5+b*I1)
  Wmut<-(1-c*I2)*(1+b*R*I2+10*b*R*I2^5+b*(1-R)*I1+10*b*(1-R)*I1^5)
  return(Wmut-Wres)
}


##1.A.Example of a PIP for the level of investment

##Set of parameters and variables
b<-3#benefit
c<-1#cost
R=3/4#Relatedness (scenario without local competition, R=R_GL)
N_reso=500#Resolution of plots
I_1<-seq(0,1,length=N_reso)#Resident investment level
I_2<-seq(0,1,length=N_reso)#Mutant investment level

##Invasion outcome matrix
invasion_matrix<-matrix(nrow=N_reso , ncol=N_reso)

##Algorithm determining outcome for each point on the grid
for (i in 1:N_reso){#Looping over resident strategies
  for (j in 1:N_reso){#Looping over mutant strategies
    IR<-I_1[i]#Setting the current resident investment level
    IM<-I_2[j]#Setting the current mutant investment level
    if(i==j){
      invasion_matrix[i,j]=0#Invasion fitness equals 0 along the diagonal for which mutant and investment levels are identical
    }
    else{
      invasion_matrix[i,j]=DeltaW(IR,IM,b,c,R)#Invasion fitness 
    }
  }
}

##Plot settings
x_var<-c(0,1,0.2)##Normalised values of investment to write on the axes
ncol=800
jet.colors <- colorRampPalette(c(rep("white",400),rep("black",400)))#Scale used to mask the intensity of selection (and focus on the sign)
palet<-jet.colors(ncol)
pal<-list(palet)

##Pairwise invasibility plot
multiplePlot("","","","",c(""),c("")
             ,ncol=ncol,x_var,x_var,list(invasion_matrix),
             abs="Resident log10 [Enzyme R1R2]",ord="Mutant log10 [Enzyme R1M1]",
             scale=c(-5,5),palette=pal,cextext=2,image=TRUE,pcex=1.25,
             subcex=2.5,labcex=1.75,axcex=1,globcex=0.625,legcex=0.75,
             contourlab=TRUE,meth="edge",contcex=0.75,colorkey="TRUE")

##1.B.Example of a PIP with mutual invasion and fixation areas

##Invasion outcome matrix
invasion_matrix_detail<-matrix(nrow=N_reso , ncol=N_reso)

##Algorithm determining outcome for each point on the grid
for (i in 1:N_reso){
  for (j in 1:N_reso){
    IR<-I_1[i]
    IM<-I_2[j]
    if(i==j){##Diagonal where mutant=resident
      invasion_matrix_detail[i,j]=0
    }
    else if(DeltaW(IR,IM,b,c,R)>0 && DeltaW(IM,IR,b,c,R)<0){##Mutual invasion IM/IR occurs
      invasion_matrix_detail[i,j]=2
    }
    else if(DeltaW(IR,IM,b,c,R)>0 && DeltaW(IM,IR,b,c,R)>0){##Full invasion/Fixation occurs
      invasion_matrix_detail[i,j]=1
    }
    else{##Loss of the mutant
      invasion_matrix_detail[i,j]=-1
    }
  }
}

##Plot settings
ncol=1200
jet.colors <- colorRampPalette(c(rep("white",400),##White=mutant loss
                                 rep("coral",400),##coral=mutual invasion
                                 rep("blue",400)))##blue=full invasion
palet<-jet.colors(ncol)
pal<-list(palet)

##Pairwise invasibility plot
multiplePlot("","","","",c(""),c("")
             ,ncol=ncol,x_var,x_var,list(invasion_matrix_detail),
             abs="Resident investment I1",ord="Mutant investment I2",
             scale=c(-1,3),palette=pal,cextext=2,image=TRUE,pcex=1.25,
             subcex=2.5,labcex=1.75,axcex=1,globcex=0.625,legcex=0.75,
             contourlab=TRUE,meth="edge",contcex=0.75,colorkey="FALSE")

##2.Sensitivity study of fixation and invasion areas relative to R for multiplicative-linear trade-off

##Parameters
b_set<-c(8)#benefit
R_set<-c(1/4,1/2,3/4,1)#relatedness values

##Invasion outcome matrix
invasion_matrix_detail<-matrix(nrow=N_reso , ncol=N_reso)
list_inv_mat<-list()

##Algorithm determining outcome for each point on the grid
for (bi in 1:length(b_set)){## has no influence once investment is scaled against the ESS investment level, focus on a single value
  b<-b_set[bi]
  for(Ri in 1:length(R_set)){
    R<-R_set[Ri]
    I_1<-seq(0,(R*b-c)/(R*b*c),length=N_reso)
    I_2<-seq(0,(R*b-c)/(R*b*c),length=N_reso)
    for (i in 1:N_reso){
      for (j in 1:N_reso){
        IR<-I_1[i]
        IM<-I_2[j]
        if(IR>((R*b-c)/(R*b*c)) || IM>((R*b-c)/(R*b*c))){
          invasion_matrix_detail[i,j]=-2
        }
        else if(i==j){
          invasion_matrix_detail[i,j]=0
        }
        else if(DeltaW(IR,IM,b,c,R)>0 && DeltaW(IM,IR,b,c,R)<0){
          invasion_matrix_detail[i,j]=2
        }
        else if(DeltaW(IR,IM,b,c,R)>0 && DeltaW(IM,IR,b,c,R)>0){
          invasion_matrix_detail[i,j]=1
        }
        else{
          invasion_matrix_detail[i,j]=-1
        }
      }
    }
    list_inv_mat[[Ri+(bi-1)*length(R_set)]]<-invasion_matrix_detail
  }
}

##Plot settings
addtxt<-list(l=-0.05,h=1.1,txt=c("S5A.","S5B.","S5C.","S5D."),srt = 0,font=2,col="black")##Letters to identify plots
ncol=1200
jet.colors <- colorRampPalette(c(rep("white",400),##White=mutant loss
                                 rep("coral",400),##coral=mutual invasion
                                 rep("blue",400)))##blue=full invasion
palet<-jet.colors(ncol)
pal<-list(palet)

##Set of PIPs of Figure S5
multiplePlot("","","","",c(""),c("","","","")
             ,ncol=ncol,x_var,x_var,list_inv_mat,TEXT_to_Add=addtxt,
             abs="Resident",ord="Mutant",scale=c(-1,3),palette=pal,cextext=2,
             image=TRUE,pcex=1.25,subcex=2.5,labcex=1.75,axcex=1,globcex=0.625,
             legcex=0.75,contourlab=TRUE,meth="edge",contcex=0.75,colorkey="COMMON")

##3.Influence of accelerating cost functions

##Parameters
b_set<-c(5)#benefit
R_set<-c(1/4,2/4,3/4,1)#range of relatedness R=R_LT values considered

##Invasion outcomes matrix
invasion_matrix_detail<-matrix(nrow=N_reso , ncol=N_reso)
list_inv_mat<-list()

##3.A.Case 1 of accelerating costs considered in Figure S11 A.-D.
##Algorithm determining outcome for each point on the grid
for (bi in 1:length(b_set)){
  b<-b_set[bi]#Current benefit
  for(Ri in 1:length(R_set)){
    R<-R_set[Ri]#Current relatedness
    I_1<-seq(0,b/6,length=N_reso)#Resident range of investment level
    I_2<-seq(0,b/6,length=N_reso)#Mutant range of investment level
    for (i in 1:N_reso){
      for (j in 1:N_reso){
        IR<-I_1[i]
        IM<-I_2[j]
        if(i==j){
          invasion_matrix_detail[i,j]=0
        }
        else if(DeltaW_supra_cost_F2(IR,IM,b,c,R)>0 && DeltaW_supra_cost_F2(IM,IR,b,c,R)<0){
          invasion_matrix_detail[i,j]=2
        }
        else if(DeltaW_supra_cost_F2(IR,IM,b,c,R)>0 && DeltaW_supra_cost_F2(IM,IR,b,c,R)>0){
          invasion_matrix_detail[i,j]=1
        }
        else{
          invasion_matrix_detail[i,j]=-1
        }
      }
    }
    list_inv_mat[[Ri+(bi-1)*length(R_set)]]<-invasion_matrix_detail
  }
}

dev.off()

##Plot settings
addtxt<-list(l=-0.05,h=1.1,txt=c("S11A.","S11B.","S11C.","S11D."),srt = 0,font=2,col="black")##Letters to identify plots
ncol=1200
jet.colors <- colorRampPalette(c(rep("white",400),##White=mutant loss
                                 rep("coral",400),##coral=mutual invasion
                                 rep("blue",400)))##blue=full invasion
palet<-jet.colors(ncol)
pal<-list(palet)

##Set of PIPs of Figure S11 A.-D.
multiplePlot("","","","",c(""),c("","","","")
             ,ncol=ncol,x_var,x_var,list_inv_mat,TEXT_to_Add=addtxt,
             abs="Resident",ord="Mutant",scale=c(-1,3),palette=pal,cextext=2,
             image=TRUE,pcex=1.25,subcex=2.5,labcex=1.75,axcex=1,globcex=0.625,
             legcex=0.75,contourlab=TRUE,meth="edge",contcex=0.75,colorkey="COMMON")

##3.B.Case 2 of accelerating costs considered in Figure S11 E.-H.

##Invasion outcomes matrix
invasion_matrix_detail<-matrix(nrow=N_reso , ncol=N_reso)
list_inv_mat<-list()

##Algorithm determining outcome for each point on the grid
for (bi in 1:length(b_set)){
  b<-b_set[bi]#Current benefit
  for(Ri in 1:length(R_set)){
    R<-R_set[Ri]#Current relatedness
    I_1<-seq(0,b/9,length=N_reso)#Resident range of investment level
    I_2<-seq(0,b/9,length=N_reso)#Mutant range of investment level
    for (i in 1:N_reso){
      for (j in 1:N_reso){
        IR<-I_1[i]
        IM<-I_2[j]
        if(i==j){
          invasion_matrix_detail[i,j]=0
        }
        else if(DeltaW_supra_cost_F3(IR,IM,b,c,R)>0 && DeltaW_supra_cost_F3(IM,IR,b,c,R)<0){
          invasion_matrix_detail[i,j]=2
        }
        else if(DeltaW_supra_cost_F3(IR,IM,b,c,R)>0 && DeltaW_supra_cost_F3(IM,IR,b,c,R)>0){
          invasion_matrix_detail[i,j]=1
        }
        else{
          invasion_matrix_detail[i,j]=-1
        }
      }
    }
    list_inv_mat[[Ri+(bi-1)*length(R_set)]]<-invasion_matrix_detail
  }
}

dev.off()

##Plot settings
addtxt<-list(l=-0.05,h=1.1,txt=c("S11E.","S11F.","S11G.","S11H."),srt = 0,font=2,col="black")##Letters to identify plots
ncol=1200
jet.colors <- colorRampPalette(c(rep("white",400),##White=mutant loss
                                 rep("coral",400),##coral=mutual invasion
                                 rep("blue",400)))##blue=full invasion
palet<-jet.colors(ncol)
pal<-list(palet)

##Set of PIPs of Figure S11 E.-H.
multiplePlot("","","","",c(""),c("","","",""),ncol=ncol,
             x_var,x_var,list_inv_mat,TEXT_to_Add=addtxt,
             abs="Resident",ord="Mutant",scale=c(-1,3),palette=pal,cextext=2,
             image=TRUE,pcex=1.25,subcex=2.5,labcex=1.75,axcex=1,globcex=0.625,
             legcex=0.75,contourlab=TRUE,meth="edge",contcex=0.75,colorkey="COMMON")

##3.C.Case of ombined costs including a sub-linear component of Figure S12

##Parameters
b_set<-c(9)
R_set<-c(1/4,2/4,3/4,1)#range of relatedness R=R_LT values considered

##Invasion outcome matrix
invasion_matrix_detail<-matrix(nrow=N_reso , ncol=N_reso)
list_inv_mat<-list()

##Algorithm determining outcome for each point on the grid
for (bi in 1:length(b_set)){
  b<-b_set[bi]#Current benefit
  for(Ri in 1:length(R_set)){
    R<-R_set[Ri]#Current relatedness
    I_1<-seq(0,b/12,length=N_reso)#Resident range of investment level
    I_2<-seq(0,b/12,length=N_reso)#Mutant range of investment level
    for (i in 1:N_reso){
      for (j in 1:N_reso){
        IR<-I_1[i]
        IM<-I_2[j]
        if(i==j){
          invasion_matrix_detail[i,j]=0
        }
        else if(DeltaW_sublinear_cost(IR,IM,b,c,R)>0 && DeltaW_sublinear_cost(IM,IR,b,c,R)<0){
          invasion_matrix_detail[i,j]=2
        }
        else if(DeltaW_sublinear_cost(IR,IM,b,c,R)>0 && DeltaW_sublinear_cost(IM,IR,b,c,R)>0){
          invasion_matrix_detail[i,j]=1
        }
        else{
          invasion_matrix_detail[i,j]=-1
        }
      }
    }
    list_inv_mat[[Ri+(bi-1)*length(R_set)]]<-invasion_matrix_detail
  }
}

dev.off()
##Plot settings
addtxt<-list(l=-0.05,h=1.1,txt=c("S12A.","S12B.","S12C.","S12D."),srt = 0,font=2,col="black")##Letters to identify plots
ncol=1200
jet.colors <- colorRampPalette(c(rep("white",400),##White=mutant loss
                                 rep("coral",400),##coral=mutual invasion
                                 rep("blue",400)))##blue=full invasion
palet<-jet.colors(ncol)
pal<-list(palet)

##Set of PIPs of Figure S12
multiplePlot("","","","",c(""),c("","","","")
             ,ncol=ncol,x_var,x_var,list_inv_mat,TEXT_to_Add=addtxt,
             abs="Resident",ord="Mutant",scale=c(-1,3),palette=pal,cextext=2,
             image=TRUE,pcex=1.25,subcex=2.5,labcex=1.75,axcex=1,globcex=0.625,
             legcex=0.75,contourlab=TRUE,meth="edge",contcex=0.75,colorkey="COMMON")

##3.D.Case of supra-linear benefits of Figure S13

##Parameters
b_set<-c(5,10)#Range of benefits considered
R_set<-c(0.25,0.5,0.75,1)#Range of relatedness R=R_LT

##Invasion outcome matrix
invasion_matrix_detail<-matrix(nrow=N_reso , ncol=N_reso)
list_inv_mat<-list()

##Algorithm determining outcome for each point on the grid
for (bi in 1:length(b_set)){
  b<-b_set[bi]#Current benefit
  for(Ri in 1:length(R_set)){
    R<-R_set[Ri]#Current benefit
    I_1<-seq(0,b/4,length=N_reso)#Resident investment level
    I_2<-seq(0,b/4,length=N_reso)#Mutant investment level
    for (i in 1:N_reso){
      for (j in 1:N_reso){
        IR<-I_1[i]
        IM<-I_2[j]
        if(i==j){
          invasion_matrix_detail[i,j]=0
        }
        else if(DeltaW_supralinear_benefit5(IR,IM,b,c,R)>0 && DeltaW_supralinear_benefit5(IM,IR,b,c,R)<0){
          invasion_matrix_detail[i,j]=2
        }
        else if(DeltaW_supralinear_benefit5(IR,IM,b,c,R)>0 && DeltaW_supralinear_benefit5(IM,IR,b,c,R)>0){
          invasion_matrix_detail[i,j]=1
        }
        else{
          invasion_matrix_detail[i,j]=-1
        }
      }
    }
    list_inv_mat[[Ri+(bi-1)*length(R_set)]]<-invasion_matrix_detail
  }
}

dev.off()

##Plot settings
addtxt<-list(l=-0.05,h=1.1,txt=c("S13A.","S13B.","S13C.","S13D.",
                                 "S13E.","S13F.","S13G.","S13H."),srt = 0,font=2,col="black")##Letters to identify plots
ncol=1200
jet.colors <- colorRampPalette(c(rep("white",400),##White=mutant loss
                                 rep("coral",400),##coral=mutual invasion
                                 rep("blue",400)))##blue=full invasion
palet<-jet.colors(ncol)
pal<-list(palet)

##Set of PIPs of Figure S13 (b=5:A.-D.; b=10:E.-H.)
multiplePlot("","","","",c("",""),c("","","","")
             ,ncol=ncol,x_var,x_var,list_inv_mat,TEXT_to_Add=addtxt,
             abs="Resident",ord="Mutant",scale=c(-1,3),palette=pal,cextext=2,
             image=TRUE,pcex=1.25,subcex=2.5,labcex=1.75,axcex=1,globcex=0.625,
             legcex=0.75,contourlab=TRUE,meth="edge",contcex=0.75,colorkey="COMMON")

                              ###END###

