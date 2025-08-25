setwd('/Users/cuijiajun/Desktop/2023-2024\ PhD\ ucl/2024_aMeta/wholepipeAt_Ps/2024_233_analysis/mapdamage/mapdamage2024_29/toAt/')
args <- commandArgs(trailingOnly = TRUE)
print(args)
tableAt<-read.table(paste('/Users/cuijiajun/Desktop/2023-2024\ PhD\ ucl/2024_aMeta/wholepipeAt_Ps/2024_233_analysis/mapdamage/mapdamage2024_29/toAt/',args[1],'/misincorporation.txt',sep=''),header = TRUE)
samplename=args[2]
#samplename='PL0065_mapped_At_q20.sorted.mapDamage'
#tableAt<-read.table('./PL0065_mapped_At_q20.sorted.mapDamage/misincorporation.txt',header = T)
subtable<-tableAt[which(tableAt$End=='5p'),3:22]
#freq ref:#occurrences of mutations / #occurrences of the reference nucleotide
#https://ginolhac.github.io/mapDamage/#a4
#subtable[,9:29]<-subtable[,9:29]/subtable[,8]
totaldfall<-subtable[,c(3:7)]
#sum strand + and -
A<-aggregate(A ~ Pos, data = totaldfall, sum)
A$total<-A$A
A<-A[,c(1,3)]
A$base='A'
C<-aggregate(C ~ Pos, data = totaldfall, sum)
C$total<-C$C
C<-C[,c(1,3)]
C$base='C'

G<-aggregate(G ~ Pos, data = totaldfall, sum)
G$total<-G$G
G<-G[,c(1,3)]
G$base='G'

T<-aggregate(T ~ Pos, data = totaldfall, sum)
T$total<-T$T
T<-T[,c(1,3)]
T$base='T'

dfall<-NULL
#df1<-aggregate(C.T ~ Pos, data = subtable, sum)
#df2<-data.frame(pos=df1$Pos,value=df1$C.T,group='CtoT')
#rbind(df2,dfall)
for (i in 9:20){
  df1<-aggregate(subtable[,i] ~ Pos, data = subtable, sum)
  df2<-data.frame(pos=df1$Pos,value=df1[,2],group=colnames(subtable)[i])
  dfall<-rbind(df2,dfall)
}

#the order of total should be GGTTCCAATACG
unique(dfall$group)
#total site of ACGT list for 840 (70x12 substitution type)
total<-rbind(G,G,T,T,C,C,A,A,T,A,C,G)
library(ggplot2)
dfall<-cbind(dfall,total)
dfall$proportion<-dfall$value/dfall$total
dfall$pos_from_5end<-dfall$pos



dfone<-dfall[which(dfall$group=='C.T'),]
dfone$sample<-samplename
colnames(dfone)
write.table(dfone,'./allinoneCtoT.txt', append = TRUE,row.names = FALSE,col.names = FALSE,quote=FALSE)

