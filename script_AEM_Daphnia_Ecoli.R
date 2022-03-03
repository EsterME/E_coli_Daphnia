setwd("C:/Users/Admin/Desktop/CNR/E. coli")
library("ggplot2")
library("reshape2")
library("cowplot")
ERIC<-read.csv("ERIC.csv")
row.names(ERIC)<-ERIC[,1]
ERIC<-ERIC[-1,]
tERIC<-t(ERIC)
dERIC<-dist(ERIC, method = "binary")
plot(hclust(dERIC, method="average"),hang=-1, main='binary', sub='', xlab='', cex=1.0) #plot cluster analysis of betapair
plot(hclust(dERIC, method="complete"),hang=-1, main='binary', sub='', xlab='', cex=1.0) #plot cluster analysis of betapair


###Pie charts for Sequence
setwd("C:/Users/Admin/Desktop/CNR/E.coli/E. coli/Daphnia/")
ST<-read.csv("ST_summary.csv")
head(ST)
ST1727<-
  ggplot(ST, aes(x="", y=ST1727, fill=X)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  guides(fill=guide_legend(nrow=3,byrow=TRUE, title = ""))+
  ggtitle("ST1727")+
  scale_fill_manual(values=c("maroon", "gold", "indianred3", "grey83", "tan", "aquamarine4", "darkgreen", "steelblue"))+
  theme_void()+
  theme(legend.position = 'top',  legend.spacing.x = unit(0.1, 'cm'), legend.spacing.y = unit(0.001, 'cm'))+
  theme(plot.title = element_text(hjust = 0.5))+
  annotate("text", x = 0, y = 0, label = "n=111")
ST38<-ggplot(ST, aes(x="", y=ST38, fill=X)) +
  geom_bar(stat="identity", width=1, color="white", show.legend = F) +
  ggtitle("ST38")+
  coord_polar("y", start=0) +
  scale_fill_manual(values=c("maroon", "gold", "indianred3", "grey83", "tan", "aquamarine4", "darkgreen", "steelblue"))+
  theme(legend.position = "none")+
   theme_void()+
  theme(plot.title = element_text(hjust = 0.5))+
  annotate("text", x = 0, y = 0, label = "n=3979")

ST3573<-ggplot(ST, aes(x="", y=ST3573, fill=X)) +
  geom_bar(stat="identity", width=1, color="white", show.legend = F) +
  ggtitle("ST3573")+
  coord_polar("y", start=0) +
  scale_fill_manual(values=c("maroon", "gold", "indianred3", "grey83", "tan", "aquamarine4", "darkgreen", "steelblue"))+
   theme_void()+
  theme(plot.title = element_text(hjust = 0.5))+
  annotate("text", x = 0, y = 0, label = "n=1")
ST4166<-
  ggplot(ST, aes(x="", y=ST4166, fill=X)) +
  geom_bar(stat="identity", width=1, color="white", show.legend = F) +
  coord_polar("y", start=0) +
  ggtitle("ST4166")+
  scale_fill_manual(values=c("maroon", "gold", "indianred3", "grey83", "tan", "aquamarine4", "darkgreen", "steelblue"))+
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5))+
  annotate("text", x = 0, y = 0, label = "n=3")
legend <- cowplot::get_legend(ST1727)
plot_grid(legend,ST1727+ theme(legend.position = "none"),ST38,ST3573,ST4166, ncol = 1, rel_heights = c(6,10,15,7,7))

rownames(ST)<-ST[,1]
ST<-ST[,-1]
prST<-prop.table(as.matrix(ST),2)


##Release
Re<-read.csv("release.csv")

Re$treatment[Re$treatment=="C"]<-"Control"
Re$treatment[Re$treatment=="D"]<-"Daphnia"

ggplot(Re, aes(x=strain, y=events, color=treatment, fill=treatment)) +
  geom_boxplot(size=1, alpha=0.5)  + 
  geom_point(position = position_jitterdodge(jitter.width = 0.1), color="black", size=1) +
  labs(y=expression(italic(E.coli) ~ ml^{-1}), )+
  theme_ipsum(plot_title_size = 12)+
    scale_color_manual(values=c("#8c8c8c", "#cccccc"))+
  scale_fill_manual(values=c("#8c8c8c", "#cccccc"))+
  theme(legend.position = "bottom",legend.text=element_text(size=14), legend.title = element_text(size=14), axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16), axis.title.x = element_text(size = 16))



library(MASS)
library(car)
library(emmeans)
library(ggplot2)
Re157<-Re[1:9,]
Re1<-Re[10:19,]
glm157<-lm(log(events)~treatment, data=Re157)
glm1<-lm(log(events)~treatment, data=Re1)
check_model(glm157)
Anova(glm157)

###uidA
uid<-read.csv("uida.csv")


library("ggplot2")
library("hrbrthemes")
ggplot(uid, aes(x=simple, y=uidA, color=simple)) +
  geom_boxplot(size=1, alpha=0.8)  + 
  ylab("uidA per animal or ml-1")+
  xlab("")+
  scale_color_manual(values=c("#8c8c8c", "#8c8c8c", "#8c8c8c"))+
  #scale_fill_manual(values=c("#8c8c8c", "#cccccc", "black"))+
  theme_ipsum(plot_title_size = 12)+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0), "cm"))



###Chemostat experiment


CFU<-read.csv("cfu_n.csv")
rownames(CFU)<-CFU[,1]
CFU<-as.matrix(CFU[,-1])
colnames(CFU)<-c("1","3","6","8","10","daph_end")



library("reshape2")
library("ggplot2")
library("hrbrthemes")
library("RColorBrewer")
library("performance")
CFU2<-as.matrix(CFU[,1:5])
mCFU<-melt(CFU2)

di<-ggplot(mCFU, aes(x=Var2, y=value))+
  geom_point(aes(color=Var1, size=3), show.legend = F)+
  geom_line(aes(group=Var1, color=Var1)) +
  labs(y = expression(paste(log ~  italic(   E.coli) ~ ml^{-1} )),x="days", color="# Daphnia")+
  scale_color_gradient(low="cyan3", high="black")+
  #theme(legend.text=element_text(size=14), legend.title = element_text(size=14), axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
        #axis.title.y = element_text(size = 16), axis.title.x = element_text(size = 16))+
  theme_ipsum(base_size = 14, axis_title_size = 14)
di + scale_y_log10()

CFU<-read.csv("cfu_n.csv")
rownames(CFU)<-CFU[,1]
CFU<-as.matrix(CFU[,-1])
CFU2<-as.data.frame(CFU[,1:5])


nrD<-as.integer(rownames(CFU))
nrE<-as.integer(CFU[,6])
glmB<-glm.nb(CFU2$X10~nrD)

check_model(glmB)
Anova(glmB)

glmE<-glm.nb(CFU2$X10~nrE)

check_model(glmE)
Anova(glmE)




cfu<-read.csv("cfu_new.csv")

plot(cfu$red,cfu$green)
plot(cfu$ED1, cfu$ED157)
library("reshape2")
library("ggplot2")
library("viridis")
library("hrbrthemes")
library("cowplot")

cfuRel<-cfu[,c(1,6:10)]
mcfu<-melt(cfuRel)

library(lme4)
library(performance)

#mcfu<-subset(mcfu, mcfu$date!="23-Mar")


mcfuA<-subset(mcfu, treatment=="alive")
mcfuC<-subset(mcfu, treatment=="carapax")
mcfuG<-subset(mcfu, treatment=="gut")
mcfuLB<-subset(mcfu, treatment=="LB")

glm2<-glmer(as.integer(mcfu$value)~treatment+(1|date), data=mcfu, family=poisson(link="log"))
Anova(glm2)

lm2<-lmer(log(as.integer(mcfu$value)+1)~variable+treatment+(1|date), data=mcfu)
check_model(lm2)




glmnb1<-glm.nb(as.integer(mcfu$value)~treatment, data=mcfu)

glmnb2<-glmer(as.integer(mcfuA$value)~variable+(1|date), family=neg.bin(theta =16973003770 ),data=mcfuA)


lm2<-lmer(as.integer(mcfu$value)~treatment+(1|date), data=mcfu)
check_model(lm2)



mcfu<-melt(cfuRel)
mcfu<-subset(mcfu, mcfu$date!="28-Mar")

mcfuA<-subset(mcfu, treatment=="alive")
mcfuC<-subset(mcfu, treatment=="carapax")
mcfuG<-subset(mcfu, treatment=="gut")
mcfuLB<-subset(mcfu, treatment=="LB")

PA<-ggplot(mcfuA, aes(x=date, y=value, color=variable, fill=variable)) +
  geom_boxplot(size=1, alpha=0.5)  + 
  geom_point(position = position_jitterdodge(), color="black", size=0.2) +
  ggtitle("Alive daphnids")+
  labs(y = expression(paste(italic(E.coli) ~ ml^{-1} ))) +
  theme_ipsum(plot_title_size = 16)+
  scale_y_log10(limits=c(10^4, 10^9))+
  scale_color_manual(values=c("#8c8c8c", "#cccccc"))+
  scale_fill_manual(values=c("#8c8c8c", "#cccccc"))+
  theme(legend.position = "none", legend.text=element_text(size=14), legend.title = element_text(size=14), axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 12))+
  theme(plot.margin = unit(c(0, 0.5, 0.5, 0), "cm"))

PC<-ggplot(mcfuC, aes(x=date, y=value, color=variable, fill=variable)) +
  geom_boxplot(size=1, alpha=0.5)  + 
  geom_point(position = position_jitterdodge(), color="black", size=0.2) +
  ggtitle("Carapax")+
  labs(y = expression(paste(italic(E.coli) ~ ml^{-1} )),fill="strain", color="strain") +
  theme_ipsum(plot_title_size = 16)+
  scale_y_log10(limits=c(10^4, 10^9))+
  scale_color_manual(values=c("#8c8c8c", "#cccccc"))+
  scale_fill_manual(values=c("#8c8c8c", "#cccccc"))+
  theme(legend.position = "bottom", legend.text=element_text(size=14), legend.title = element_text(size=14), axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 12))+
  theme(plot.margin = unit(c(0, 0.5, 0.5, 0), "cm"))

PG<-ggplot(mcfuG, aes(x=date, y=value, color=variable, fill=variable)) +
  geom_boxplot(size=1, alpha=0.5)  + 
  geom_point(position = position_jitterdodge(), color="black", size=0.2) +
  ggtitle("Gut")+
  ylab("")+
  labs(fill="strain", color="strain")+
  theme_ipsum(plot_title_size = 16)+
  scale_y_log10(limits=c(10^4, 10^9))+
  scale_color_manual(values=c("#8c8c8c", "#cccccc"))+
  scale_fill_manual(values=c("#8c8c8c", "#cccccc"))+
  theme(legend.position = "bottom", legend.text=element_text(size=14), legend.title = element_text(size=14), axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 12))+
  theme(plot.margin = unit(c(0, 0.5, 0.5, 0), "cm"))

PLB<-ggplot(mcfuLB, aes(x=date, y=value, color=variable, fill=variable)) +
  geom_boxplot(size=1, alpha=0.5)  + 
  geom_point(position = position_jitterdodge(), color="black", size=0.2) +
  ggtitle("Strains only")+
  ylab("")+
  scale_color_manual(values=c("#8c8c8c", "#cccccc"))+
  scale_fill_manual(values=c("#8c8c8c", "#cccccc"))+
  theme_ipsum(plot_title_size = 16)+
  scale_y_log10(limits=c(10^4, 10^9))+
  theme(legend.position = "none", legend.text=element_text(size=14), legend.title = element_text(size=14), axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 12))+
  theme(plot.margin = unit(c(0, 0.5, 0.5, 0), "cm"))

plot_grid(PA,PLB,PC, PG, rel_heights = c(1,1.1))





mcfu2<-subset(mcfu, mcfu$date!="23-Mar")


mcfuA<-subset(mcfu2, treatment=="alive")
mcfuC<-subset(mcfu2, treatment=="carapax")
mcfuG<-subset(mcfu2, treatment=="gut")
mcfuLB<-subset(mcfu2, treatment=="LB")

A<-ggplot(mcfuA, aes(x=variable, y=value)) +
  geom_boxplot(size=1, alpha=0.5)  + 
  geom_jitter(color="black", size=1, alpha=0.8, width=0.2) +
  ggtitle("Alive daphnids")+
  ylab("E.coli ml-1")+
  theme_ipsum(plot_title_size = 12)+
  scale_y_log10(limits=c(10^3, 10^9))+
  scale_color_manual(values=c("#8c8c8c", "#cccccc"))+
  scale_fill_manual(values=c("#8c8c8c", "#cccccc"))+
  theme(legend.position = "bottom")

library(MASS)
glm2<-glm.nb(as.integer(mcfuA$value)~variable+(1|date),data=mcfuA)

glm2<-glmer(as.integer(mcfuA$value)~variable+(1|date), family=neg.bin(theta =16973003770 ),data=mcfuA)
library(performance)
check_model(glm2)
Anova(glm2)

#theme(plot.margin = unit(c(0, 0.5, 0.5, 0), "cm"))
C<-ggplot(mcfuC, aes(x=variable, y=value)) +
  geom_boxplot(size=1, alpha=0.5)  + 
  geom_jitter(color="black", size=1, alpha=0.8) +
  ggtitle("carapax")+
  ylab("E.coli ml-1")+
  theme_ipsum(plot_title_size = 12)+
  scale_y_log10(limits=c(10^3, 10^9))+
  scale_color_manual(values=c("#8c8c8c", "#cccccc"))+
  scale_fill_manual(values=c("#8c8c8c", "#cccccc"))+
  theme(legend.position = "bottom")
G<-ggplot(mcfuG, aes(x=variable, y=value)) +
  geom_boxplot(size=1, alpha=0.5)  + 
  geom_jitter(color="black", size=1, alpha=0.8) +
  ggtitle("gut")+
  ylab("E.coli ml-1")+
  theme_ipsum(plot_title_size = 12)+
  scale_y_log10(limits=c(10^3, 10^9))+
  scale_color_manual(values=c("#8c8c8c", "#cccccc"))+
  scale_fill_manual(values=c("#8c8c8c", "#cccccc"))+
  theme(legend.position = "bottom")
LB<-ggplot(mcfuLB, aes(x=variable, y=value)) +
  geom_boxplot(size=1, alpha=0.5)  + 
  geom_jitter(color="black", size=1, alpha=0.8) +
  ggtitle("alone")+
  ylab("E.coli ml-1")+
  theme_ipsum(plot_title_size = 12)+
  scale_y_log10(limits=c(10^3, 10^9))+
  scale_color_manual(values=c("#8c8c8c", "#cccccc"))+
  scale_fill_manual(values=c("#8c8c8c", "#cccccc"))+
  theme(legend.position = "bottom")
grid.arrange(A,LB,C, G, nrow=2)
plot_grid(A,LB,C, G, rel_heights = c(1,1.1))

