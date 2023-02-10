rm(list=ls())
library(ggplot2)
library(INLA)
library(readxl)
library(RColorBrewer)
library(scales)
library(sf)
library(tmap)


#################################################################################
## Figure 1: Standardized mortality ratio (SMR) for the different cancer sites ##
##           computed by year for the whole of Spanish provinces               ##
#################################################################################
load("CancerData_SpainPROV.Rdata")

T <- length(unique(data[[1]]$Year))
T.from <- min(data[[1]]$Year)
T.to <- max(data[[1]]$Year)

Table <- lapply(data, function(xx){
  aux <- data.frame(O=aggregate(xx$O, by=list(Year=xx$Year), sum)$x,
                    E=aggregate(xx$E, by=list(Year=xx$Year), sum)$x)
  aux$SMR <- round(aux$O/aux$E,3)
  aux
})

Table.SMR <- data.frame(Year=rep(seq(1,T),length(data)),
                        Cancer=rep(as.factor(names(data)),each=T),
                        SMR=as.numeric(unlist(lapply(Table, function(x) x$SMR))))

Figure1 <- ggplot(data=Table.SMR, aes(x=Year, y=SMR, group=Cancer)) +
  geom_line(aes(color=Cancer), lwd=1.1) +
  geom_hline(yintercept=1, linetype="dashed") + 
  xlab("Year") + ylab("SMR") + 
  scale_x_continuous(breaks=seq(1,T,2), labels=seq(T.from,T.to,2)) +
  theme(legend.position=c(0.9,0.85),
        legend.text=element_text(size=14),
        legend.title=element_text(size=16),
        axis.text=element_text(size=14),
        axis.title=element_text(size=14, face="bold"))
print(Figure1)

ggsave("Figure1.eps", Figure1, width=10, height=6)


##########################################################################################################
## Table 4: Summary statistics (mean, standard deviation, minimum, quantiles and maximum)               ##
##          of standardized mortality ratio (SMR) for the different cancer sites in some selected years ##
##########################################################################################################

## Lung cancer ##
O.T1 <- data$Lung[data$Lung$Year==2006,"SMR"]
O.T2 <- data$Lung[data$Lung$Year==2009,"SMR"]
O.T3 <- data$Lung[data$Lung$Year==2013,"SMR"]
O.T4 <- data$Lung[data$Lung$Year==2017,"SMR"]
O.T5 <- data$Lung[data$Lung$Year==2020,"SMR"]

aux <- do.call(rbind,lapply(list("2006"=O.T1,"2009"=O.T2,"2013"=O.T3,"2017"=O.T4,"2020"=O.T5), function(x){
  data.frame(mean=mean(x), sd=sd(x), min=min(x), q1=quantile(x,0.25), q2=quantile(x,0.5), q3=quantile(x,0.75), max=max(x))
}))
round(aux,2)   


## Colorectal cancer ##
O.T1 <- data$Colorectal[data$Colorectal$Year==2006,"SMR"]
O.T2 <- data$Colorectal[data$Colorectal$Year==2009,"SMR"]
O.T3 <- data$Colorectal[data$Colorectal$Year==2013,"SMR"]
O.T4 <- data$Colorectal[data$Colorectal$Year==2017,"SMR"]
O.T5 <- data$Colorectal[data$Colorectal$Year==2020,"SMR"]

aux <- do.call(rbind,lapply(list("2006"=O.T1,"2009"=O.T2,"2013"=O.T3,"2017"=O.T4,"2020"=O.T5), function(x){
  data.frame(mean=mean(x), sd=sd(x), min=min(x), q1=quantile(x,0.25), q2=quantile(x,0.5), q3=quantile(x,0.75), max=max(x))
}))
round(aux,2)    


## Stomach cancer ##
O.T1 <- data$Stomach[data$Stomach$Year==2006,"SMR"]
O.T2 <- data$Stomach[data$Stomach$Year==2009,"SMR"]
O.T3 <- data$Stomach[data$Stomach$Year==2013,"SMR"]
O.T4 <- data$Stomach[data$Stomach$Year==2017,"SMR"]
O.T5 <- data$Stomach[data$Stomach$Year==2020,"SMR"]

aux <- do.call(rbind,lapply(list("2006"=O.T1,"2009"=O.T2,"2013"=O.T3,"2017"=O.T4,"2020"=O.T5), function(x){
  data.frame(mean=mean(x), sd=sd(x), min=min(x), q1=quantile(x,0.25), q2=quantile(x,0.5), q3=quantile(x,0.75), max=max(x))
}))
round(aux,2)  


## LOCP cancer ##
O.T1 <- data$LOCP[data$LOCP$Year==2006,"SMR"]
O.T2 <- data$LOCP[data$LOCP$Year==2009,"SMR"]
O.T3 <- data$LOCP[data$LOCP$Year==2013,"SMR"]
O.T4 <- data$LOCP[data$LOCP$Year==2017,"SMR"]
O.T5 <- data$LOCP[data$LOCP$Year==2020,"SMR"]

aux <- do.call(rbind,lapply(list("2006"=O.T1,"2009"=O.T2,"2013"=O.T3,"2017"=O.T4,"2020"=O.T5), function(x){
  data.frame(mean=mean(x), sd=sd(x), min=min(x), q1=quantile(x,0.25), q2=quantile(x,0.5), q3=quantile(x,0.75), max=max(x))
}))
round(aux,2)  


###########################################################################
## Table 5: Model selection criteria and computational time (in minutes) ##
###########################################################################
compute.DIC <- function(x){
  data.frame(mean.deviance=x$dic$mean.deviance,
             p.eff=x$dic$p.eff,
             DIC=x$dic$dic,
             WAIC=x$waic$waic,
             LS=-sum(log(x$cpo$cpo)),
             Time=round(x$cpu.used[4])/60)
}

## NOTE: The models must be previously fitted by using the `Fit_models.R` script ##
load("INLAmodels_MCAR_Leroux_M2.Rdata")

Table <- do.call(rbind,lapply(list(Type.I=LCAR.t1,Type.II=LCAR.t2,Type.III=LCAR.t3,Type.IV=LCAR.t4), compute.DIC))
round(Table,1)


###################################################################################################
## Figure 2: Maps of posterior median estimates of $\exp{\theta_{ij}}$ (top) and                 ##
##           posterior exceedence probabilities for $Pr(\exp{\theta_{ij}}>1 | {\bf O})$ (bottom) ##
###################################################################################################
Model <- LCAR.t2

S <- nrow(carto)
J <- length(data)

## Posterior mean estimates (top) ##
spatial.mean <- lapply(Model$marginals.random$idx, function(x) inla.emarginal(function(y) exp(y),x))
spatial.mean <- matrix(unlist(spatial.mean),S,J,byrow=F)
colnames(spatial.mean) <- c("Lung","Colorectal","Stomach","LOCP")

carto.aux <- cbind(carto,spatial.mean)
paleta <- brewer.pal(8,"YlOrRd")
values <- c(-Inf,0.77,0.83,0.91,1,1.1,1.20,1.30,Inf)

Map1 <- tm_shape(carto.aux) +
  tm_polygons(col=c("Lung","Colorectal","Stomach","LOCP"), id="name", palette=paleta,
              title="", legend.show=T, legend.reverse=T,
              style="fixed", breaks=values, interval.closure="left") +
  tm_layout(main.title="", main.title.position=0.35, panel.label.size=1.5,
            panel.labels=paste(c("Lung","Colorectal","Stomach","LOCP"),"cancer"),
            legend.outside=T, legend.outside.position="right", legend.frame=F,
            legend.outside.size=0.2, outer.margins=c(0.02,0.01,0.02,0.01)) + 
  tm_facets(nrow=2, ncol=2)

tmap_mode("plot")
tmap_save(Map1, file="Figure2.eps")


## Posterior exceedence probabilities ##
spatial.prob <- lapply(Model$marginals.random$idx, function(x) 1-inla.pmarginal(log(1),x))
spatial.prob <- matrix(unlist(spatial.prob),S,J,byrow=F)
colnames(spatial.prob) <- c("Lung","Colorectal","Stomach","LOCP")

carto.aux <- cbind(carto,spatial.prob)
paleta <- brewer.pal(6,"Blues")[-1]
values <- c(0,0.1,0.2,0.8,0.9,1)

Map2 <- tm_shape(carto.aux) + 
  tm_polygons(col=c("Lung","Colorectal","Stomach","LOCP"), id="name", palette=paleta,
              title="", legend.show=T, legend.reverse=T,
              style="fixed", breaks=values, interval.closure="left",
              labels=c("[0-0.1)","[0.1-0.2)","[0.2-0.8)","[0.8-0.9)","[0.9-1]")) + 
  tm_layout(main.title="", main.title.position=0.35, panel.label.size=1.5,
            panel.labels=paste(c("Lung","Colorectal","Stomach","LOCP"),"cancer"),
            legend.outside=T, legend.outside.position="right", legend.frame=F,
            legend.outside.size=0.2, outer.margins=c(0.02,0.01,0.02,0.01)) + 
  tm_facets(nrow=2, ncol=2)

tmap_mode("plot")
tmap_save(Map2, file="Figure2b.eps")


#############################################################################################
## Table 6: Estimated between-disease spatial correlations (posterior medians and 95% CI). ##
##          Significant correlations are highlighted in bold.                              ##
#############################################################################################
round(Model$summary.cor[,c("0.5quant","0.025quant","0.975quant")],3)


###########################################################################
## Figure 3: Posterior median estimates of year-specific mortality risks ##
##           $\exp(\gamma_{tj})$ and 95\% credible intervals             ##
###########################################################################
marginals <- list(Model$marginals.random$idy.1,Model$marginals.random$idy.2,
                  Model$marginals.random$idy.3,Model$marginals.random$idy.4)

x <- 1:T
title <- paste(c("Lung","Colorectal","Stomach","LOCP"),"cancer")

graphics.off()
postscript("Figure3.eps", width=10, height=7, horizontal=FALSE)
par(mfrow=c(2,2), pty="m")
k <- 1

for(i in marginals){
  aux <- lapply(i, function(x) inla.tmarginal(function(y) exp(y),x))
  temporal.mean <- unlist(lapply(aux, function(x) inla.emarginal(function(y) y,x)))
  q1 <- unlist(lapply(aux, function(x) inla.qmarginal(0.025,x)))
  q2 <- unlist(lapply(aux, function(x) inla.qmarginal(0.975,x)))
  
  plot(range(x), c(0.75, 1.35), type="n", 
       xlab="",ylab="", xaxt="n", main=title[k], cex.lab=1.3, cex.axis=1.3, cex.main=1.5)
  axis(1, at=round(seq(1,T,2)), labels=round(seq(T.from,T.to,2)), las=0, cex.axis=1.3)
  X.Vec <- c(x, tail(x, 1), rev(x), x[1])
  Y.Vec <- c(q1, tail(q2, 1), rev(q2), q1[1])
  polygon(X.Vec, Y.Vec, col = "grey", border = NA)
  lines(temporal.mean)
  abline(h=1, lty=2)
  
  k <- k+1
}  
dev.off()


###########################################################################################
## Figure 4: Maps of posterior median estimates of LUNG cancer mortality risks $R_{itj}$ ##
###########################################################################################
risk <- matrix(Model$summary.fitted.values$`0.5quant`[Model$.args$data$ID.disease==1],S,T,byrow=F)
colnames(risk) <- paste("Year",seq(T.from,T.to),sep=".")

carto.aux <- cbind(carto,risk)
paleta <- brewer.pal(8,"YlOrRd")
values <- c(-Inf,0.77,0.83,0.91,1,1.1,1.20,1.30,Inf)

Map <- tm_shape(carto.aux) +
  tm_polygons(col=paste("Year",c(2006,2009,2012,2014,2017,2020),sep="."), id="name", palette=paleta,
              title="", legend.show=T, legend.reverse=T,
              style="fixed", breaks=values, interval.closure="left") +
  tm_layout(main.title="Lung cancer mortality data", main.title.position=0.25, panel.label.size=1.5,
            panel.labels=paste("Year",c(2006,2009,2012,2014,2017,2020)),
            legend.outside=T, legend.outside.position="right", legend.frame=F,
            legend.outside.size=0.2, outer.margins=c(0.02,0.01,0.02,0.01)) + 
  tm_facets(nrow=2, ncol=3)

tmap_mode("plot")
tmap_save(Map, file="Figure4.eps")


#################################################################################################
## Figure 5: Maps of posterior median estimates of COLORECTAL cancer mortality risks $R_{itj}$ ##
#################################################################################################
risk <- matrix(Model$summary.fitted.values$`0.5quant`[Model$.args$data$ID.disease==2],S,T,byrow=F)
colnames(risk) <- paste("Year",seq(T.from,T.to),sep=".")

carto.aux <- cbind(carto,risk)
paleta <- brewer.pal(8,"YlOrRd")
values <- c(-Inf,0.77,0.83,0.91,1,1.1,1.20,1.30,Inf)

Map <- tm_shape(carto.aux) +
  tm_polygons(col=paste("Year",c(2006,2009,2012,2014,2017,2020),sep="."), id="name", palette=paleta,
              title="", legend.show=T, legend.reverse=T,
              style="fixed", breaks=values, interval.closure="left") +
  tm_layout(main.title="Colorectal cancer mortality data", main.title.position=0.2, panel.label.size=1.5,
            panel.labels=paste("Year",c(2006,2009,2012,2014,2017,2020)),
            legend.outside=T, legend.outside.position="right", legend.frame=F,
            legend.outside.size=0.2, outer.margins=c(0.02,0.01,0.02,0.01)) + 
  tm_facets(nrow=2, ncol=3)

tmap_mode("plot")
tmap_save(Map, file="Figure5.eps")

##############################################################################################
## Figure 6: Maps of posterior median estimates of STOMACH cancer mortality risks $R_{itj}$ ##
##############################################################################################
risk <- matrix(Model$summary.fitted.values$`0.5quant`[Model$.args$data$ID.disease==3],S,T,byrow=F)
colnames(risk) <- paste("Year",seq(T.from,T.to),sep=".")

carto.aux <- cbind(carto,risk)
paleta <- brewer.pal(8,"YlOrRd")
values <- c(-Inf,0.77,0.83,0.91,1,1.1,1.20,1.30,Inf)

Map <- tm_shape(carto.aux) +
  tm_polygons(col=paste("Year",c(2006,2009,2012,2014,2017,2020),sep="."), id="name", palette=paleta,
              title="", legend.show=T, legend.reverse=T,
              style="fixed", breaks=values, interval.closure="left") +
  tm_layout(main.title="Stomach cancer mortality data", main.title.position=0.2, panel.label.size=1.5,
            panel.labels=paste("Year",c(2006,2009,2012,2014,2017,2020)),
            legend.outside=T, legend.outside.position="right", legend.frame=F,
            legend.outside.size=0.2, outer.margins=c(0.02,0.01,0.02,0.01)) + 
  tm_facets(nrow=2, ncol=3)

tmap_mode("plot")
tmap_save(Map, file="Figure6.eps")

###########################################################################################
## Figure 7: Maps of posterior median estimates of LOCP cancer mortality risks $R_{itj}$ ##
###########################################################################################
risk <- matrix(Model$summary.fitted.values$`0.5quant`[Model$.args$data$ID.disease==4],S,T,byrow=F)
colnames(risk) <- paste("Year",seq(T.from,T.to),sep=".")

carto.aux <- cbind(carto,risk)
paleta <- brewer.pal(8,"YlOrRd")
values <- c(-Inf,0.77,0.83,0.91,1,1.1,1.20,1.30,Inf)

Map <- tm_shape(carto.aux) +
  tm_polygons(col=paste("Year",c(2006,2009,2012,2014,2017,2020),sep="."), id="name", palette=paleta,
              title="", legend.show=T, legend.reverse=T,
              style="fixed", breaks=values, interval.closure="left") +
  tm_layout(main.title="LOCP cancer mortality data", main.title.position=0.25, panel.label.size=1.5,
            panel.labels=paste("Year",c(2006,2009,2012,2014,2017,2020)),
            legend.outside=T, legend.outside.position="right", legend.frame=F,
            legend.outside.size=0.2, outer.margins=c(0.02,0.01,0.02,0.01)) + 
  tm_facets(nrow=2, ncol=3)

tmap_mode("plot")
tmap_save(Map, file="Figure7.eps")
