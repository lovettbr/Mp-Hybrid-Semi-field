#Transgenic Metarhizium rapidly kills mosquitoes in
#a malaria-endemic region of Burkina Faso
#21 April 2019
#Brian Lovett
#Department of Entomology
#University of Maryland
#lovettbr@umd.edu

####Load packages####
#install.packages(c("cowplot", "tidyverse", "plotly", "survival", "MASS", "transformr", "scales", "reshape2", "RColorBrewer"))
#install.packages("gganimate", dependencies=TRUE)
library(cowplot)
library(tidyverse)
library(plotly)
library(survival)
library(MASS)
library(transformr)
library(gganimate)
library(scales)
library(reshape2)
library(RColorBrewer)

####This script creates an R function to generate raincloud plots, then simulates ####
### data for plots. If using for your own data, you only need lines 1-80. 
### It relies largely on code previously written by David Robinson 
### (https://gist.github.com/dgrtwo/eb7750e74997891d7c20) and ggplot2 by H Wickham
### From https://github.com/RainCloudPlots/RainCloudPlots
### source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")
#Load packages
#packages <-c("ggplot2", "dplyr", "lavaan", "plyr", "cowplot", "rmarkdown", "readr", "caTools", "bitops")
#if(length(setdiff(packages, rownames(installed.packages())))>0){install.packages(setdiff(packages, rownames(installed.packages())))}
library(ggplot2)

# Defining the geom_flat_violin function. Note: the below code modifies the 
# existing github page by removing a parenthesis in line 50

"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomFlatViolin <-
  ggproto("GeomFlatViolin", Geom,
          setup_data = function(data, params) {
            data$width <- data$width %||%
              params$width %||% (resolution(data$x, FALSE) * 0.9)
            
            # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
            data %>%
              group_by(group) %>%
              mutate(ymin = min(y),
                     ymax = max(y),
                     xmin = x,
                     xmax = x + width / 2)
            
          },
          
          draw_group = function(data, panel_scales, coord) {
            # Find the points for the line to go all the way around
            data <- transform(data, xminv = x,
                              xmaxv = x + violinwidth * (xmax - x))
            
            # Make sure it's sorted properly to draw the outline
            newdata <- rbind(plyr::arrange(transform(data, x = xminv), y),
                             plyr::arrange(transform(data, x = xmaxv), -y))
            
            # Close the polygon: set first and last point the same
            # Needed for coord_polar and such
            newdata <- rbind(newdata, newdata[1,])
            
            ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
          },
          
          draw_key = draw_key_polygon,
          
          default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
                            alpha = NA, linetype = "solid"),
          
          required_aes = c("x", "y")
  )


####Release-recapture Survival####
dat = read.csv("Release_recapture_field.csv")
colnames(dat)[5:18]=1:14

dat2 = dat %>%
  gather(key="Day", value = "Died", -Replicate, -Treatment, -Feeding.Status, -Location) %>%
  dplyr::select(Replicate, Treatment, Day, Died) %>%
  dplyr::mutate(Died=as.numeric(Died), Day=as.factor(Day)) %>%
  dplyr::mutate(Day=factor(Day, levels(Day)[c(1, 7:14, 2:6, 15)])) %>%
  dplyr::group_by(Replicate, Treatment, Day) %>%
  dplyr::summarize(Died=sum(Died)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Replicate, Treatment) %>%
  dplyr::mutate(Percent=1-cumsum(Died)/sum(Died), n=sum(Died)) %>%
  dplyr::filter(Day!="Still.Alive") %>%
  dplyr::mutate(Day=as.numeric(Day)) %>%
  dplyr::ungroup() 

dat2.infected = dat %>%
  gather(key="Day", value = "Died", -Replicate, -Treatment, -Feeding.Status, -Location) %>%
  dplyr::select(Replicate, Treatment, Day, Died) %>%
  dplyr::filter(Day!="Still.Alive") %>%
  dplyr::mutate(Died=as.numeric(Died), Day=as.factor(Day)) %>%
  dplyr::mutate(Day=factor(Day, levels(Day)[c(1, 7:14, 2:6, 15)])) %>%
  dplyr::group_by(Replicate, Treatment, Day) %>%
  dplyr::summarize(Died=sum(Died)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Replicate, Treatment) %>%
  dplyr::mutate(Percent=1-cumsum(Died)/sum(Died), n=sum(Died)) %>%
  dplyr::mutate(Day=as.numeric(Day)) %>%
  dplyr::ungroup() 

dat3= dat2 %>%
  dplyr::group_by(Treatment, Day) %>%
  dplyr::summarize(per.mean=mean(Percent), se=sd(Percent)/sqrt(length(Percent)), Replicates=length(Percent))

limits=aes(ymax=per.mean+se, ymin=per.mean-se)
cbPalette <- c("#007acc", "#cc0000", "#29a329")
theme = theme_bw()+theme(text = element_text(size=25), axis.title.x = element_text(size=30), axis.title.y = element_text(size=30), title = element_text(size=35), legend.position="none")
plt=ggplot(dat3, aes(Day, per.mean, color=Treatment))+geom_line(size=3)+
  geom_errorbar(limits, width=.2, size=2)+theme+
  scale_colour_manual(values=cbPalette)+
  xlab("Days after capture from semi-field")+ylab("Percent survival")+
  scale_y_continuous(labels=percent, limits=c(0:1))+scale_x_continuous(breaks=0:max(dat2$Day))
plt

rere.dat.LT.infected=subset(dat2.infected, Day > 3 & Day <10) %>%
  dplyr::group_by(Replicate, Treatment) %>%
  dplyr::mutate(Alive=Percent*n, Dead=n-Alive) %>%
  dplyr::select(Replicate, Treatment, Day, Dead, Alive) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Replicate, Treatment) %>%
  dplyr::summarize(LT=as.numeric(dose.p(glm(cbind(Alive,Dead)~Day, binomial), p=surv.per))) %>%
  dplyr::mutate(LT = replace(LT, which(LT<0 | LT>30), NA))

rere.dat.LT.infected.error=rere.dat.LT.infected %>% ungroup() %>%
  dplyr::group_by(Treatment) %>%
  dplyr::summarize(mean= mean(LT, na.rm=T), se=sd(LT, na.rm=T)/sqrt(length(LT[!is.na(LT)])), Replicate=length(LT[!is.na(LT)]))

dat.loc=dat %>%
  gather(key="Day", value = "Died", -Replicate, -Treatment, -Feeding.Status, -Location) %>%
  group_by(Location, Replicate, Treatment) %>%
  dplyr::summarize(mosq=sum(Died, na.rm=T)) %>%
  ungroup() %>%
  group_by(Location) %>%
  dplyr::summarize(mean=mean(mosq), se=sd(mosq)/sqrt(length(mosq)), Replicates=length(mosq))

dat.loc.per=dat %>%
  gather(key="Day", value = "Died", -Replicate, -Treatment, -Feeding.Status, -Location) %>%
  group_by(Replicate, Treatment) %>%
  dplyr::mutate(n=sum(Died, na.rm=T)) %>%
  ungroup() %>%
  group_by(Location, Replicate, Treatment) %>%
  dplyr::summarize(mosq.per=sum(Died, na.rm=T)/unique(n)) %>%
  ungroup() %>%
  group_by(Location) %>%
  dplyr::summarize(per=mean(mosq.per), se=sd(mosq.per)/sqrt(length(mosq.per)), Replicates=length(mosq.per))

dat.fed.per=dat %>%
  gather(key="Day", value = "Died", -Replicate, -Treatment, -Feeding.Status, -Location) %>%
  group_by(Replicate, Treatment) %>%
  dplyr::mutate(n=sum(Died, na.rm=T)) %>%
  ungroup() %>%
  group_by(Feeding.Status, Replicate, Treatment) %>%
  dplyr::summarize(mosq.per=sum(Died, na.rm=T)/unique(n)) %>%
  ungroup() %>%
  group_by(Feeding.Status) %>%
  dplyr::summarize(per=mean(mosq.per), se=sd(mosq.per)/sqrt(length(mosq.per)), Replicates=length(mosq.per))

theme = theme_bw()+theme(text = element_text(size=25), axis.title.x = element_text(size=30), axis.title.y = element_text(size=30), title = element_text(size=35), legend.position="none")
limits=aes(ymax=per+se, ymin=per-se)
loc.per.plt=ggplot(dat.loc.per, aes(Location, per, fill=Location))+geom_bar(stat="identity", position="dodge")+geom_errorbar(limits, width=.2, size=2)+theme+scale_fill_brewer(palette = "Dark2")+
  xlab("Location captured")+ylab("Number of mosquitoes")+scale_y_continuous(labels = percent_format(accuracy=1))
loc.per.plt

####Mortality####
dat.m <- read.csv("Generational_mortality_field.csv")

dat.m2 = dat.m %>%
  dplyr::group_by(Day, Treatment, Replicate) %>%
  dplyr::summarize(sm=sum(Adults)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Treatment=factor(Treatment, levels = c("Control", "RFP", "Hybrid")))

####Generational Development Fecundity####
gen.dat <- read.csv("Generational_development_field.csv")

gen.dat.sum= gen.dat %>%
  dplyr::filter(Stage!="None") %>%
  dplyr::group_by(Replicate, Day, Treatment, Stage) %>%
  dplyr::summarize(sum=sum(Number, na.rm=T)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Treatment=factor(Treatment, levels=levels(Treatment)[c(1,3,2)]))

theme = theme_bw()+theme(text = element_text(size=30), axis.title.x = element_text(size=40), axis.title.y = element_text(size=40), title = element_text(size=40), plot.title = element_text(hjust = 0.5))
cbPalette <- c("#007acc","#cc0000", "#29a329")
gen.plt=ggplot()+ scale_fill_brewer(palette="Spectral")+ scale_color_manual(values=cbPalette)+
  geom_bar(data=gen.dat.sum, aes(Day, sum, fill=Stage), stat="identity")+ facet_wrap(~Replicate+Treatment) + theme +
  geom_line(data=dat.m2, aes(Day, sm, color=Treatment), size=2) +
  ylab("Individuals")+xlab("Time (Days)")
gen.plt

gen.plt=ggplot()+ scale_fill_brewer(palette="Spectral")+ scale_color_manual(values=cbPalette)+
  geom_bar(data=subset(gen.dat.sum, Replicate==1), aes(Day, sum, fill=Stage), stat="identity")+facet_wrap(~Treatment) + theme +
  geom_line(data=subset(dat.m2, Replicate==1), aes(Day, sm, color=Treatment), size=2) +
  ylab("Individuals")+xlab("Time (Days)")+
  ggtitle("Semi-Field Established Population\nF1 and F2 Offspring")
gen.plt

ggplotly(gen.plt)

####Persistence####
per.dat <- read.csv("Persistence_field.csv")

per.dat.mort = per.dat %>%
  group_by(Week, Replicate, Treatment, Species) %>%
  dplyr::mutate(Percent=1-cumsum(Mortality)/sum(Mortality), n=sum(Mortality)) %>%
  dplyr::filter(Day!="Alive") %>%
  dplyr::mutate_at("Day", as.character) %>%
  dplyr::mutate_at("Day", as.numeric)

per.dat.mn= per.dat.mort %>%
  dplyr::group_by(Species, Week, Treatment, Day) %>%
  dplyr::summarize(mean=mean(Percent), se=sd(Percent)/sqrt(length(Percent)))

limits=aes(ymax=mean+se, ymin=mean-se)
theme = theme_bw()+theme(text = element_text(size=20), axis.title.x = element_text(size=30), axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size=20), axis.text.y = element_text(size=25), title = element_text(size=35), plot.title = element_text(hjust = 0.5), legend.title = element_text(size=25), legend.text = element_text(size=20))
cbPalette <- c("#007acc", "#29a329", "#cc0000")

per.plt1=ggplot(per.dat.mn, aes(Day, mean, color=Treatment))+geom_line(size=3)+
  geom_errorbar(limits, width=.4, size=2)+theme+scale_colour_manual(values=cbPalette)+
  xlab("Days after exposure")+ylab("Percent survival")+
  scale_y_continuous(labels=percent)+scale_x_continuous(breaks=0:max(per.dat.mn$Day))+facet_wrap(~Species+Week)
per.plt1

per.plt2=ggplot(per.dat.mn, aes(Day, mean, color=Treatment, frame=Week))+geom_line(size=3)+
  geom_errorbar(limits, width=.4, size=2)+theme+scale_colour_manual(values=cbPalette)+
  xlab("Days after exposure")+ylab("Percent survival")+
  scale_y_continuous(labels=percent)+scale_x_continuous(breaks=0:max(per.dat.mn$Day))+facet_wrap(~Species)+
  transition_time(Week) +
  labs(title = 'Week: {frame_time}')+
  ease_aes('sine-in-out')

animate(per.plt2, fps=5, width=1000)
#anim_save("per.plt2.gif")

####Lab Fecundity####
f.lab.dat= read.csv("Development_lab.csv")

f.lab.dat4=f.lab.dat %>%
  dplyr::group_by(Replicate, Treatment) %>%
  dplyr::summarize(Adults=mean(Adults, na.rm=T), Pupae=mean(Pupae, na.rm=T),
                   Larvae=mean(Larvae, na.rm=T), Eggs=mean(Eggs, na.rm=T))

f.lab.dat5=f.lab.dat %>%
  dplyr::filter(Eggs!=0) %>%
  dplyr::group_by(Replicate, Treatment) %>%
  dplyr::summarize(Lay=mean(Time.to.lay, na.rm=T)) %>%
  dplyr::ungroup()

t.test(f.lab.dat5[f.lab.dat5$Treatment=="Control",]$Lay, f.lab.dat5[f.lab.dat5$Treatment=="hybrid",]$Lay, paired =T)
t.test(f.lab.dat5[f.lab.dat5$Treatment=="Control",]$Lay, f.lab.dat5[f.lab.dat5$Treatment=="rfp",]$Lay, paired =T)

f.lab.dat5.error=f.lab.dat5 %>%
  dplyr::select(Lay, Treatment) %>%
  dplyr::group_by(Treatment) %>%
  dplyr::summarize(Lay.mn=mean(Lay), Lay.se=sd(Lay, na.rm = T)/sqrt(length(Lay)))

f.lab.dat6=f.lab.dat %>%
  dplyr::filter(!is.na(Time.to.lay)) %>%
  dplyr::select(Replicate, Treatment, Time.to.lay) %>%
  dplyr::group_by(Treatment) %>%
  dplyr::count(Time.to.lay) %>%
  dplyr::right_join(data.frame(Treatment=rep(c("Control", "hybrid", "rfp"), 13), Time.to.lay=rep(0:12, each=3), n=0), by=c("Treatment", "Time.to.lay")) %>%
  dplyr::group_by(Treatment, Time.to.lay) %>%
  dplyr::mutate(n=sum(n.x, n.y, na.rm=T)) %>%
  dplyr::group_by(Treatment) %>%
  dplyr::select(Treatment, Time.to.lay, n) %>%
  dplyr::mutate(Lay.sum=cumsum(n))

f.lab.dat7= f.lab.dat %>%
  dplyr::mutate(Laid=Eggs>0, Laid=replace_na(Laid, "FALSE")) %>%
  dplyr::select(Treatment, Laid) %>%
  dplyr::group_by(Treatment) %>%
  dplyr::summarize(total=length(Laid)) %>%
  dplyr::left_join(f.lab.dat6, by=c("Treatment")) %>%
  dplyr::mutate(Percent.laid=Lay.sum/total)

f.lab.dat7$Treatment=recode(f.lab.dat7$Treatment, hybrid="Hybrid", rfp="RFP")
f.lab.dat7$Treatment=factor(f.lab.dat7$Treatment, levels=levels(f.lab.dat7$Treatment)[c(1,3,2)])

lay.plt2=ggplot(f.lab.dat7)+geom_step(aes(x=Time.to.lay, y=Percent.laid, color=Treatment), size=2)+theme+
  scale_color_manual(values=cbPalette)+scale_y_continuous(labels=percent, limits=c(0,1))+scale_x_continuous(breaks=seq(0,12, by=2))+
  ylab("Percent of egg laying females") + xlab("Days post infection") + theme(legend.title=element_text(size=25), legend.text=element_text(size=20))
lay.plt2

f.lab.dat7=f.lab.dat %>%
  group_by(Replicate, Treatment) %>%
  dplyr::select(Adults, Pupae, Larvae, Eggs, Replicate, Treatment) %>%
  gather(key=Stage, value=count, -Replicate, -Treatment) %>%
  dplyr::ungroup() %>%
  group_by(Treatment, Stage) %>%
  dplyr::summarize(mn=mean(count, na.rm=T), sd=sd(count, na.rm = T), reps=length(count))

f.lab.dat8=f.lab.dat %>%
  group_by(Replicate, Treatment) %>%
  dplyr::select(Adults, Pupae, Larvae, Eggs, Replicate, Treatment) %>%
  gather(key=Stage, value=count, -Replicate, -Treatment) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Stage=as.factor(Stage), Treatment=recode(Treatment, rfp="RFP", hybrid="Hybrid")) %>%
  dplyr::mutate(Stage=factor(Stage, levels(Stage)[c(2:4, 1)]), Treatment=factor(Treatment, levels(Treatment)[c(2,3,1)]))

raincloud_theme = theme(
  text = element_text(size = 20),
  axis.title.x = element_text(size = 25),
  axis.title.y = element_text(size = 25),
  axis.text = element_text(size = 20, face="bold"),
  axis.text.x = element_text(angle = 45, vjust = 0.5, face="bold"),
  legend.title=element_text(size=20),
  legend.text=element_text(size=20),
  legend.position = "right",
  plot.title = element_text(lineheight=.8, face="bold", size = 16),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

g2 <- ggplot(data = subset(f.lab.dat8, count!=0), aes(y = count, x = Treatment, fill = Stage)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .7, adjust=7) +
  geom_point(aes(y = count, color = Stage),position = position_jitterdodge(0.2, dodge=.3), size = 1, alpha = 0.8) +
  geom_boxplot(width = .3, outlier.shape = NA, alpha=0.7) +
  scale_fill_brewer(palette="Dark2")+ scale_color_brewer(palette="Dark2")+
  theme_bw() + xlab("Treatment")+ylab("Number of Mosquitoes")+ coord_flip()+
  raincloud_theme
g2

f.lab.dat2= f.lab.dat %>%
  dplyr::group_by(Replicate, Treatment) %>%
  dplyr::summarize(Adults=sum(Adults, na.rm=T), Pupae=sum(Pupae, na.rm=T)-Adults,
                   Larvae=sum(Larvae, na.rm=T)-Pupae-Adults, Eggs=sum(Eggs, na.rm=T)-Larvae-Pupae-Adults,
                   Total=sum(Adults, Pupae, Larvae, Eggs), Fecundity=sum(Adults, Pupae, Larvae, Eggs),
                   Fertility=Adults/Fecundity) %>%
  dplyr::select(Replicate, Treatment, Adults, Pupae, Larvae, Eggs) %>%
  gather(key=Stage, value=count, -Replicate, -Treatment) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Stage=as.factor(Stage), Treatment=recode(Treatment, rfp="RFP", hybrid="Hybrid")) %>%
  dplyr::mutate(Stage=factor(Stage, levels(Stage)[c(2:4, 1)]), Treatment=factor(Treatment, levels(Treatment)[c(1,3,2)])) %>%
  dplyr::group_by(Treatment, Stage) %>%
  dplyr::summarize(mean=mean(count, na.rm=T), se=sd(count, na.rm = T)/sqrt(length(count)), n=length(count))

f.lab.dat2.individual= f.lab.dat %>%
  dplyr::group_by(Treatment, Individual.Mosquito.Number) %>%
  dplyr::summarize(Adults=sum(Adults, na.rm=T), Pupae=sum(Pupae, na.rm=T),
                   Larvae=sum(Larvae, na.rm=T), Eggs=sum(Eggs, na.rm=T),
                   Fecundity=sum(Adults, Pupae, Larvae, Eggs),
                   Fertility=Adults/Fecundity) %>%
  dplyr::select(Treatment, Adults, Pupae, Larvae, Eggs) %>%
  gather(key=Stage, value=count, -Treatment) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Stage=as.factor(Stage), Treatment=recode(Treatment, rfp="RFP", hybrid="Hybrid")) %>%
  dplyr::mutate(Stage=factor(Stage, levels(Stage)[c(2:4, 1)]), Treatment=factor(Treatment, levels(Treatment)[c(1,3,2)])) %>%
  group_by(Treatment, Stage) %>%
  dplyr::summarize(mean=mean(count, na.rm=T), se=sd(count, na.rm = T)/sqrt(length(count)), n=length(count))

f.lab.dat8= f.lab.dat %>%
  dplyr::mutate(Laid=Eggs>0, Laid=replace_na(Laid, "FALSE")) %>%
  dplyr::select(Treatment, Replicate, Laid, Time.to.lay) %>%
  dplyr::group_by(Treatment, Replicate) %>%
  dplyr::summarize(Percent.Laid=length(Laid[Laid=="TRUE"])/length(Laid)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Treatment=factor(Treatment, levels(Treatment)[c(1,3,2)]), Treatment=recode(Treatment, rfp="RFP", hybrid="Hybrid")) %>%
  dplyr::group_by(Treatment) %>%
  dplyr::summarize(mean=mean(Percent.Laid, na.rm=T), se=sd(Percent.Laid, na.rm = T)/sqrt(length(Percent.Laid)), n=length(Percent.Laid))

limits=aes(ymin=mean-se, ymax=mean+se)
cbPalette <- c("#007acc", "#cc0000", "#29a329")

gen.plt1=ggplot(data=f.lab.dat8, aes(x=Treatment, y=mean, fill=Treatment))+geom_bar(stat="identity", position="dodge")+
  geom_errorbar(limits, position=position_dodge(.9), width=.2) + theme +
  scale_fill_manual(values=cbPalette)+ scale_y_continuous(labels=percent, limits=c(0,1))+
  ylab("Percent of egg-laying females")
gen.plt1

gen.plt2=ggplot(data=f.lab.dat2.individual, aes(x=Treatment, y=mean, fill=Stage))+geom_bar(stat="identity", position="dodge")+
  geom_errorbar(limits, position=position_dodge(.9), width=.2) + theme + scale_fill_brewer(palette="Spectral")+
  ylab("Mean number of mosquitoes")
gen.plt2
