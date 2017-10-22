


#Housekeeping

rm(list=ls(all=TRUE))
library(ggplot2)
#----------------------------------------------------
#Define the parameters of the graphs to be generated#
#----------------------------------------------------
setwd("/Users/rodrigoazuero/Dropbox/Didris/FINAL/Graficas")



#Read the file

PR<-read.csv("Private_institutions.csv")

#For quality of graph, put the value in between for NA:
for(ii in 1:14){
  for(jj in 2:8){
    a=PR[ii,jj]
    if(is.na(a) & ii<14){
      PR[ii,jj]=(PR[ii-1,jj]+PR[ii+1,jj])/2
    }
  }
}

#setting up basics of the graph
sizeline=2
#Setting size of lines

p<-ggplot(PR, aes(year,size=sizesT)) + 
  geom_line(aes(y = varBRA, colour = "Brazil"),size=sizeline) + 
  geom_line(aes(y = varCHI, colour = "Chile"),size=sizeline)+
  geom_line(aes(y = varCOL, colour = "Colombia"),size=sizeline)+
  geom_line(aes(y = varFRA, colour = "France"),size=sizeline)+
  geom_line(aes(y = varITA, colour = "Italy"),size=sizeline)+
  geom_line(aes(y = varSPA, colour = "Spain"),size=sizeline)+
  geom_line(aes(y = varUSA, colour = "USA"),size=sizeline)

p<-p+labs(title="",x="Year",y="%")
p<-p+scale_colour_brewer(palette="Set1",name="")
p<-p + theme(
  plot.title = element_text(size = rel(1.5)),
  axis.title=element_text(size = rel(1.5)),
  axis.text.x=element_text(size = rel(1.5)),
  axis.text.y=element_text(size = rel(1.5)),
  legend.text=element_text(size = rel(1.2)))


pline<-p
dev.set()
pdf(file='Private_institutionsR.pdf', width=7, height=4)
pline
dev.off()


pline<-p
dev.set()
pdf(file='/Users/rodrigoazuero/Dropbox/Didris/Document/Private_institutionsR.pdf', width=7, height=4)
pline
dev.off()

  
  
#==================  
#Quality evolution
#==================


Q<-read.csv("componentescalidad.csv")
sizeline=2

p <- ggplot(Q, aes(period, theta_dif))
p<-p+geom_smooth(span=1.3,se=FALSE,color="#3399CC",size=sizeline)
p<-p+labs(title="",x="Year",y="Difference in test scores - deciles")
p<-p+scale_colour_brewer(palette="Set1",name="")
p<-p + theme(
  plot.title = element_text(size = rel(1.5)),
  axis.title=element_text(size = rel(1.5)),
  axis.text.x=element_text(size = rel(1.5)),
  axis.text.y=element_text(size = rel(1.5)))

pline<-p
dev.set()
pdf(file='theta_difference.pdf', width=7, height=4)
pline
dev.off()


pline<-p
dev.set()
pdf(file='/Users/rodrigoazuero/Dropbox/Didris/Document/theta_difference.pdf', width=7, height=4)
pline
dev.off()




p <- ggplot(Q, aes(period, prof_dif))
p<-p+geom_smooth(span=1.3,se=FALSE,color="#3399CC",size=sizeline)
p<-p+labs(title="",x="Year",y="Difference in professors")
p<-p+scale_colour_brewer(palette="Set1",name="")
p<-p + theme(
  plot.title = element_text(size = rel(1.5)),
  axis.title=element_text(size = rel(1.5)),
  axis.text.x=element_text(size = rel(1.5)),
  axis.text.y=element_text(size = rel(1.5)))

pline<-p
dev.set()
pdf(file='prof_difference.pdf', width=7, height=4)
pline
dev.off()


pline<-p
dev.set()
pdf(file='/Users/rodrigoazuero/Dropbox/Didris/Document/prof_difference.pdf', width=7, height=4)
pline
dev.off()

#--------------------#
#Numero universidades#
#--------------------#


Q<-read.csv("numero_universidades.csv")
sizeline=2

p <- ggplot(Q, aes(anno, existe))
p<-p+geom_line(size=sizeline,color="#3399CC") 
p<-p+ylim(c(0,320))
p<-p+labs(title="",x="Year",y="Number of universities")
p<-p + theme(
  plot.title = element_text(size = rel(1.5)),
  axis.title=element_text(size = rel(1.5)),
  axis.text.x=element_text(size = rel(1.5)),
  axis.text.y=element_text(size = rel(1.5)))

pline<-p
dev.set()
pdf(file='numero_universidades.pdf', width=7, height=4)
pline
dev.off()


pline<-p
dev.set()
pdf(file='/Users/rodrigoazuero/Dropbox/Didris/Document/numero_universidades.pdf', width=7, height=4)
pline
dev.off()




p <- ggplot(Q, aes(anno, promedio))
p<-p+geom_line(size=sizeline,color="#3399CC") 
p<-p+labs(title="",x="Year",y="Average size (students)")
p<-p + theme(
  plot.title = element_text(size = rel(1.5)),
  axis.title=element_text(size = rel(1.5)),
  axis.text.x=element_text(size = rel(1.5)),
  axis.text.y=element_text(size = rel(1.5)))

pline<-p
dev.set()
pdf(file='numero_universidades2.pdf', width=7, height=4)
pline
dev.off()


pline<-p
dev.set()
pdf(file='/Users/rodrigoazuero/Dropbox/Didris/Document/numero_universidades2.pdf', width=7, height=4)
pline
dev.off()


#===========================
#Enteringcohorttotal_lpoly
#---------------------------
Q<-read.csv("enteringcohorttotal_lpoly.csv")

source('/Users/rodrigoazuero/Dropbox/Didris/FINAL/Graficas/ggplot_dual_axis.R', echo=TRUE)
library(grid)
library(gtable)

sizeline=1.6
# Create plots
sizeline=2
head(Q)
rel2<-Q$rel
Q<-cbind(Q,rel2)
p <- ggplot(Q, aes(x=period2))
p<-p+geom_line(aes(x=period2,y=rel,color=c("% Students with financial aid")))
p<-p+geom_line(aes(x=period2,y=rel2,color=c(" Students")))
p<-p+scale_color_manual(name="",values=c("#3399CC","red"))
p<-p+labs(x="Period", y="% Students with financial aid")
p<-p + theme(
  plot.title = element_text(size = rel(1.5)),
  axis.title=element_text(size = rel(1.5)),
  axis.text.y=element_text(size = rel(1.8)),
  axis.text.y=element_text(size = rel(1.8)),
  legend.text=element_text(size = rel(1.2)),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.position= c(0.6, 0.2),
  legend.background = element_rect(fill=alpha('white', 0.0)))

plot1.y=p

p <- ggplot(Q, aes(x=period2,y=v1,color="Size of entering cohort"))
p<-p+geom_smooth(span=0.3,se="FALSE")
p<-p+scale_color_manual(name="",values="red")
p<-p+labs(x="Period", y="Size of entering cohort")
p<-p + theme(
  plot.title = element_text(size = rel(1.5)),
  axis.title=element_text(size = rel(1.5)),
  axis.text.x=element_text(size = rel(1.8)),
  axis.text.y=element_text(size = rel(1.8)),
  legend.text=element_text(size = rel(1.2)),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.position= c(0.6, 0.1),
  legend.background = element_rect(fill=alpha('white', 0.0)))


plot2.y=p
# Run function
ggplot_dual_axis(plot1.y, plot2.y, "y")


pline<-p
dev.set()
pdf(file='enteringcohorttotal_lpoly.pdf', width=7, height=4)
ggplot_dual_axis(plot1.y, plot2.y, "y")
dev.off()


pline<-p
dev.set()
pdf(file='/Users/rodrigoazuero/Dropbox/Didris/Document/enteringcohorttotal_lpoly.pdf', width=7, height=4)
ggplot_dual_axis(plot1.y, plot2.y, "y")
dev.off()



#======
#Income credits
#======

Q<-read.csv("income_credits.csv")

sizeline=1.6
# Create plots
sizeline=2
head(Q)
sal_pr2<-Q$sal_pr
Q<-cbind(Q,sal_pr2)
p <- ggplot(Q, aes(x=periodo))
p<-p+geom_line(aes(x=periodo,y=sal_pr,color=c("% with aid")))
p<-p+geom_line(aes(x=periodo,y=sal_pr2,color=c("Average income")))
p<-p+scale_color_manual(name="",values=c("red","#3399CC"))
p<-p+labs(x="Period", y="Minimum wages (for income)")
p<-p + theme(
  plot.title = element_text(size = rel(1.5)),
  axis.title=element_text(size = rel(1.5)),
  axis.text.y=element_text(size = rel(1.8)),
  axis.text.x=element_text(size = rel(1.8)),
  legend.text=element_text(size = rel(1.4)),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.position= c(0.8, 0.4),
  legend.background = element_rect(fill=alpha('white', 0.0)))

plot1.y=p

p <- ggplot(Q, aes(x=periodo,y=acc_pr,color="%"))
p<-p+geom_line()
p<-p+scale_color_manual(name="",values="red")
p<-p+labs(x="Period", y="%")
p<-p + theme(
  plot.title = element_text(size = rel(1.5)),
  axis.title=element_text(size = rel(1.5)),
  axis.text.x=element_text(size = rel(1.8)),
  axis.text.y=element_text(size = rel(1.8)),
  legend.text=element_text(size = rel(1.4)),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.position= c(0.8, 0.2),
  legend.background = element_rect(fill=alpha('white', 0.0)))


plot2.y=p
# Run function
ggplot_dual_axis(plot1.y, plot2.y, "y")


pline<-p
dev.set()
pdf(file='income_credits.pdf', width=7, height=4)
#ggplot_dual_axis(plot1.y, plot2.y, "y")
plot1.y
dev.off()


pline<-p
dev.set()
pdf(file='/Users/rodrigoazuero/Dropbox/Didris/Document/income_credits.pdf', width=7, height=4)
ggplot_dual_axis(plot1.y, plot2.y, "y")
dev.off()



##############

##############


##############




Q<-read.csv("/Users/rodrigoazuero/Dropbox/Didris/FINAL/Graficas/categoriasicfescomparadas.csv")


p <- ggplot(Q, aes(x=periodo,y=icf_pr,color="%"))
p<-p+geom_line()
p<-p+scale_color_manual(name="",values="#3399CC")
p<-p+labs(x="Period", y="%")
p<-p + theme(
  plot.title = element_text(size = rel(1.5)),
  axis.title=element_text(size = rel(1.5)),
  axis.text.x=element_text(size = rel(1.8)),
  axis.text.y=element_text(size = rel(1.8)),
  legend.text=element_text(size = rel(1.4)),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.position=c(-10,-10),
  legend.background = element_rect(fill=alpha('white', 0.0)))


pline<-p
dev.set()
pdf(file='Categorias_ICFES_comparadas_lp.pdf', width=7, height=4)
pline
dev.off()


pline<-p
dev.set()
pdf(file='/Users/rodrigoazuero/Dropbox/Didris/Document/Categorias_ICFES_comparadas_lp.pdf', width=7, height=4)
pline
dev.off()



#=====================#
#Quality Evolution    #
#=====================#

Q<-read.csv("/Users/rodrigoazuero/Dropbox/Didris/FINAL/Graficas/Quality_evolution.csv")

sizeline=2

Q$tier2=as.factor(Q$tier2)
p <- ggplot(Q, aes(period,quality,color=tier2))
p<-p+geom_smooth(span=1.3)
p<-p+scale_fill_manual(name = "", values = c("brown","#3399CC"))
p<-p+labs(title="",x="Year",y="Quality")
p<-p+scale_colour_brewer(palette="Set1",name="")
p<-p + theme(
  plot.title = element_text(size = rel(1.5)),
  axis.title=element_text(size = rel(1.5)),
  axis.text.x=element_text(size = rel(1.5)),
  axis.text.y=element_text(size = rel(1.5)))



pline<-p
dev.set()
pdf(file='Quality_evolution.pdf', width=7, height=4)
pline
dev.off()


dev.set()
pdf(file='/Users/rodrigoazuero/Dropbox/Didris/Document/Quality_evolution.pdf', width=7, height=4)
pline
dev.off()


#======================#
#Quality Difference    #
#======================#

Q<-read.csv("/Users/rodrigoazuero/Dropbox/Didris/FINAL/Graficas/Quality_difference.csv")

sizeline=2

p <- ggplot(Q, aes(period,qual_dif2))
p<-p+geom_smooth(span=0.8,se="FALSE",colour="#3399CC")
p<-p+scale_fill_manual(name = "", values = c("brown","#3399CC"))
p<-p+labs(title="",x="Year",y="Ratio")
p<-p+scale_colour_brewer(palette="Set1",name="")
p<-p + theme(
  plot.title = element_text(size = rel(1.5)),
  axis.title=element_text(size = rel(1.5)),
  axis.text.x=element_text(size = rel(1.5)),
  axis.text.y=element_text(size = rel(1.5)))



pline<-p
dev.set()
pdf(file='Quality_difference.pdf', width=7, height=4)
pline
dev.off()


dev.set()
pdf(file='/Users/rodrigoazuero/Dropbox/Didris/Document/Quality_difference.pdf', width=7, height=4)
pline
dev.off()




#=============#
#Tuition evolution
#=============#



Q<-read.csv("/Users/rodrigoazuero/Dropbox/Didris/FINAL/Graficas/Tuition_evolution.csv")

sizeline=2


Q$tier2=as.factor(Q$tier2)
#BUSINESS
p <- ggplot(Q, aes(year,color=tier2))
p<-p+geom_line(aes(x=year,y=adm_))
p<-p+scale_colour_brewer(palette="Set1",name="")
p<-p+labs(title="Business",x="Year",y="Price")
p<-p + theme(
  plot.title = element_text(size = rel(1.5)),
  axis.title=element_text(size = rel(1.5)),
  axis.text.x=element_text(size = rel(1.5)),
  legend.text=element_text(size=rel(1.5)),
  axis.text.y=element_text(size = rel(1.5)),
  legend.position="bottom")
p1<-p


#LAW
p <- ggplot(Q, aes(year,color=tier2))
p<-p+geom_line(aes(x=year,y=der_))
p<-p+scale_colour_brewer(palette="Set1",name="")
p<-p+labs(title="Law",x="Year",y="Price")
p<-p + theme(
  plot.title = element_text(size = rel(1.5)),
  axis.title=element_text(size = rel(1.5)),
  axis.text.x=element_text(size = rel(1.5)),
  axis.text.y=element_text(size = rel(1.5)),
  legend.text=element_text(size=rel(1.5)),
  legend.position="bottom")
p2<-p

#Engineering
p <- ggplot(Q, aes(year,color=tier2))
p<-p+geom_line(aes(x=year,y=ing_))
p<-p+scale_colour_brewer(palette="Set1",name="")
p<-p+labs(title="Engineering",x="Year",y="Price")
p<-p + theme(
  plot.title = element_text(size = rel(1.5)),
  axis.title=element_text(size = rel(1.5)),
  axis.text.x=element_text(size = rel(1.5)),
  axis.text.y=element_text(size = rel(1.5)),
  legend.text=element_text(size=rel(1.5)),
  legend.position="bottom")
p3<-p

#Medicine
p <- ggplot(Q, aes(year,color=tier2))
p<-p+geom_line(aes(x=year,y=ing_))
p<-p+scale_colour_brewer(palette="Set1",name="")
p<-p+labs(title="Medicine",x="Year",y="Price")
p<-p + theme(
  plot.title = element_text(size = rel(1.5)),
  axis.title=element_text(size = rel(1.5)),
  axis.text.x=element_text(size = rel(1.5)),
  axis.text.y=element_text(size = rel(1.5)),
  legend.text=element_text(size=rel(1.5)),
  legend.position="bottom")
p4<-p


library(gridExtra)
library(grid)
grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}

p<-grid_arrange_shared_legend(p2, p1,p3,p4)



dev.set()
pdf(file='Tuition_evolution.pdf', width=7, height=4)
grid_arrange_shared_legend(p2, p1,p3,p4)
dev.off()


dev.set()
pdf(file='/Users/rodrigoazuero/Dropbox/Didris/Document/Tuition_evolution.pdf', width=7, height=4)
grid_arrange_shared_legend(p2, p1,p3,p4)
dev.off()



#=========================#
#====RESULTADOS ECAES=====#
#=========================#

#File using export delimited using "/Users/rodrigoazuero/Dropbox/Didris/FINAL/Graficas/ECAES.csv"
#Do file generating this dataset is /Users/rodrigoazuero/Dropbox/Didris/Saber Pro


Q<-read.csv("/Users/rodrigoazuero/Dropbox/Didris/FINAL/Graficas/ECAES.csv")

spaning=1

Q$tier2=as.factor(Q$tier2)
p <- ggplot(Q, aes(period2,color=tier2,group=tier2))
p<-p+geom_smooth(span=spaning,aes(x=period2,y=comescri),se="FALSE")
p<-p+scale_colour_brewer(palette="Set1",name="")
p<-p+labs(title="Written comprehension",x="Year",y="Score")
p<-p + theme(
  plot.title = element_text(size = rel(1.5)),
  axis.title=element_text(size = rel(1.5)),
  axis.text.x=element_text(size = rel(1.5)),
  axis.text.y=element_text(size = rel(1.5)),
  legend.text=element_text(size=rel(1.5)),
  legend.position="bottom")
p1<-p


dev.set()
pdf(file='ECAESwriting.pdf', width=7, height=4)
p
dev.off()


dev.set()
pdf(file='/Users/rodrigoazuero/Dropbox/Didris/Document/ECAESwriting.pdf', width=7, height=4)
p
dev.off()




p <- ggplot(Q, aes(period2,color=tier2,group=tier2))
p<-p+geom_smooth(span=spaning,aes(x=period2,y=complect),se="FALSE")
p<-p+scale_colour_brewer(palette="Set1",name="")
p<-p+labs(title="Reading comprehension",x="Year",y="Score")
p<-p + theme(
  plot.title = element_text(size = rel(1.5)),
  axis.title=element_text(size = rel(1.5)),
  axis.text.x=element_text(size = rel(1.5)),
  axis.text.y=element_text(size = rel(1.5)),
  legend.text=element_text(size=rel(1.5)),
  legend.position="bottom")
p2<-p


dev.set()
pdf(file='ECAESreading.pdf', width=7, height=4)
p
dev.off()


dev.set()
pdf(file='/Users/rodrigoazuero/Dropbox/Didris/Document/ECAESreading.pdf', width=7, height=4)
p
dev.off()



#================================#
#Academic production per students#
#================================#


Q<-read.csv("/Users/rodrigoazuero/Dropbox/Didris/FINAL/Graficas/Art_books_perstudent.csv")

spaning=1

Q$tier2=as.factor(Q$tier2)
p <- ggplot(Q, aes(year,color=tier2,group=tier2))
p<-p+geom_line(aes(x=year,y=artipref))
p<-p+scale_colour_brewer(palette="Set1",name="")
p<-p+labs(title="Articles published by professors/ student body",x="Year",y="Ratio")
p<-p + theme(
  plot.title = element_text(size = rel(1.5)),
  axis.title=element_text(size = rel(1.5)),
  axis.text.x=element_text(size = rel(1.5)),
  axis.text.y=element_text(size = rel(1.5)),
  legend.text=element_text(size=rel(1.5)),
  legend.position="bottom")
p1<-p


dev.set()
pdf(file='Art_perstudents.pdf', width=7, height=4)
p
dev.off()


dev.set()
pdf(file='/Users/rodrigoazuero/Dropbox/Didris/Document/Art_perstudents.pdf', width=7, height=4)
p
dev.off()


#Books

Q$tier2=as.factor(Q$tier2)
p <- ggplot(Q, aes(year,color=tier2,group=tier2))
p<-p+geom_line(aes(x=year,y=bookpref))
p<-p+scale_colour_brewer(palette="Set1",name="")
p<-p+labs(title="Books published by professors/ student body",x="Year",y="Ratio")
p<-p + theme(
  plot.title = element_text(size = rel(1.5)),
  axis.title=element_text(size = rel(1.5)),
  axis.text.x=element_text(size = rel(1.5)),
  axis.text.y=element_text(size = rel(1.5)),
  legend.text=element_text(size=rel(1.5)),
  legend.position="bottom")
p2<-p


dev.set()
pdf(file='Books_perstudents.pdf', width=7, height=4)
p
dev.off()


dev.set()
pdf(file='/Users/rodrigoazuero/Dropbox/Didris/Document/Books_perstudents.pdf', width=7, height=4)
p
dev.off()



#===========================#
#= Professors by categories=#
#===========================#


#Code where the function to put multiple plots is stored
source('/Users/rodrigoazuero/Dropbox/BACKUPRODRIGO/Research/Chile/RR/RfunctionsUsed/multiplot.R', echo=TRUE)


Q<-read.csv("/Users/rodrigoazuero/Dropbox/Didris/FINAL/Graficas/Planta_profesores.csv")

spaning=1



#PHD
p <- ggplot(Q, aes(year))
p<-p+geom_smooth(span=spaning,aes(x=period2,y=difCAT1_2),se="FALSE",color="#3399CC")
p<-p+labs(title="PhD",x="Year",y="Difference")
p<-p + theme(
  plot.title = element_text(size = rel(1.5)),
  axis.title=element_text(size = rel(1.5)),
  axis.text.x=element_text(size = rel(1.5)),
  axis.text.y=element_text(size = rel(1.5)),
  legend.text=element_text(size=rel(1.5)),
  legend.position="bottom")
p1<-p


#Masters
p <- ggplot(Q, aes(year))
p<-p+geom_smooth(span=spaning,aes(x=period2,y=difCAT2_2),se="FALSE",color="#3399CC")
p<-p+labs(title="Masters",x="Year",y="Difference")
p<-p + theme(
  plot.title = element_text(size = rel(1.5)),
  axis.title=element_text(size = rel(1.5)),
  axis.text.x=element_text(size = rel(1.5)),
  axis.text.y=element_text(size = rel(1.5)),
  legend.text=element_text(size=rel(1.5)),
  legend.position="bottom")
p2<-p



#Undergraduate studies
p <- ggplot(Q, aes(year))
p<-p+geom_smooth(span=spaning,aes(x=period2,y=difCAT3_2),se="FALSE",color="#3399CC")
p<-p+labs(title="Bachelor's degree",x="Year",y="Difference")
p<-p + theme(
  plot.title = element_text(size = rel(1.5)),
  axis.title=element_text(size = rel(1.5)),
  axis.text.x=element_text(size = rel(1.5)),
  axis.text.y=element_text(size = rel(1.5)),
  legend.text=element_text(size=rel(1.5)),
  legend.position="bottom")
p3<-p


#All professors
p <- ggplot(Q, aes(year))
p<-p+geom_smooth(span=spaning,aes(x=period2,y=dif_ratioTOTAL2),se="FALSE",color="#3399CC")
p<-p+labs(title="Total Professors",x="Year",y="Difference")
p<-p + theme(
  plot.title = element_text(size = rel(1.5)),
  axis.title=element_text(size = rel(1.5)),
  axis.text.x=element_text(size = rel(1.5)),
  axis.text.y=element_text(size = rel(1.5)),
  legend.text=element_text(size=rel(1.5)),
  legend.position="bottom")
p0<-p


dev.set()
pdf(file='Planta_profesores.pdf', width=7, height=4)
multiplot(p1,p2,p3,cols=2)
dev.off()


dev.set()
pdf(file='/Users/rodrigoazuero/Dropbox/Didris/Document/Planta_profesores.pdf', width=7, height=4)
multiplot(p1,p2,p3,cols=2)
dev.off()


