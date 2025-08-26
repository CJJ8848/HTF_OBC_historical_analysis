setwd('/Users/jiajuncui/Desktop/2023-2024\ PhD\ ucl/2024_aMeta/wholepipeAt_Ps/2024_233_analysis/mapdamage/mapdamage2024_29/')
Atwd='/Users/jiajuncui/Desktop/2023-2024\ PhD\ ucl/2024_aMeta/wholepipeAt_Ps/2024_233_analysis/mapdamage/mapdamage2024_29/toAt/'

setwd(Atwd)
tableAt<-read.table(paste(Atwd,'allinoneCtoT.txt',sep=''),header = F)

colnames(tableAt)<-c("pos" , "value"      ,   "group"    ,"Pos"   ,  "total"  ,       "base"     ,    "proportion", "pos_from_5end" ,"sample"   )



#all46 

tableAt17<-read.table('/Users/jiajuncui/Desktop/2023-2024\ PhD\ ucl/2024_aMeta/wholepipeAt_Ps/2024_233_analysis/mapdamage/mapdamage2024_22/for17Hold/usedmaptoAt/allinoneCtoT.txt',header=F)
colnames(tableAt17)<-c("pos" , "value"      ,   "group"    ,"Pos"   ,  "total"  ,       "base"     ,    "proportion", "pos_from_5end" ,"sample"   )

tableallAt<-rbind(tableAt,tableAt17)

library(ggplot2)
ggplot(tableallAt,aes(pos_from_5end,proportion,group=sample,color='#BA8921'))+geom_line()+(ylab('C to T Frequency'))+ theme_bw() + theme(panel.grid=element_blank())+scale_color_manual(values = c("#BA8921"))+theme(axis.text=element_text(size=12)
)+  ylim(0, 0.04)
ggsave(paste('toAtallinone_46','.png',sep=''),width = 10, height = 6)


tableAtno137no27ESP<-tableallAt[tableallAt$sample !='PL0137_mapped_At_q20.sorted.mapDamage'&tableallAt$sample !='results_27.ESP_1975_mapped_At.q20',]

library(ggplot2)
ggplot(tableAtno137no27ESP,aes(pos_from_5end,proportion,group=sample,color='#BA8921'))+geom_line()+(ylab('C to T Frequency'))+ theme_bw() + theme(panel.grid=element_blank())+scale_color_manual(values = c("#BA8921"))+theme(axis.text=element_text(size=12)
)+  ylim(0, 0.04)
ggsave(paste('toAtallinone_no137no27ESP','.png',sep=''),width = 10, height = 6)









pswd='/Users/jiajuncui/Desktop/2023-2024\ PhD\ ucl/2024_aMeta/wholepipeAt_Ps/2024_233_analysis/mapdamage/mapdamage2024_29/tops/'
setwd(pswd)
tablePs<-read.table(paste(pswd,'allinoneCtoT.txt',sep=''),header = F)

colnames(tablePs)<-c("pos" , "value"      ,   "group"    ,"Pos"   ,  "total"  ,       "base"     ,    "proportion", "pos_from_5end" ,"sample"   )

tablePs108<-tablePs[tablePs$sample =='PL0108.mapped_to_Pseudomonas.dd.q20.mapDamage',]
tablePs_17<-read.table('/Users/jiajuncui/Desktop/2023-2024\ PhD\ ucl/2024_aMeta/wholepipeAt_Ps/2024_233_analysis/mapdamage/mapdamage2024_22/for17Hold/usedforCFML_ps/allinoneCtoT.txt',header=F)
#/Users/jiajuncui/Desktop/2023-2024\ PhD\ ucl/2024_aMeta/wholepipeAt_Ps/mapdamage/mapdamage2024_22/for17Hold/usedforCFML_ps/allinoneCtoT.txt
colnames(tablePs_17)<-c("pos" , "value"      ,   "group"    ,"Pos"   ,  "total"  ,       "base"     ,    "proportion", "pos_from_5end" ,"sample"   )

tablePsall<-rbind(tablePs,tablePs_17)

library(ggplot2)
ggplot(tablePsall,aes(pos_from_5end,proportion,group=sample,color='#BA8921'))+geom_line()+(ylab('C to T Frequency'))+ theme_bw() + theme(panel.grid=element_blank())+scale_color_manual(values = c("#BA8921"))+theme(axis.text=element_text(size=12)
) + ylim(0, 0.04)
ggsave(paste('topsallinone_46','.png',sep=''),width = 10, height = 6)
#dfall[which(dfall$value>=0.0075),]

