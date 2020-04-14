pch_list=c(16,1,3,4,8,15,17,18)

### Prepare input

# Read input data
blast <-read.table(file="BLAST/runtimes/runtimesBLAST.txt", header=TRUE)
plast <-read.table(file="PLAST/runtimes/runtimesPLAST.txt", header=TRUE)
mmseqs<-read.table(file="MMseqs2/runtimes/runtimesMMseqs2.txt", header=TRUE)
ublast<-read.table(file="UBLAST/runtimes/runtimesUBLAST.txt", header=TRUE)
blat  <-read.table(file="BLAT/runtimes/runtimesBLAT.txt", header=TRUE)

# Only take PLAST without quorum
plast<-subset(plast,plast$Quorum==0)


# Take mean of several runs
blast<-aggregate(blast,by=list(blast$Colors),FUN=mean)
mmseqs<-aggregate(mmseqs,by=list(mmseqs$Colors),FUN=mean)
ublast<-aggregate(ublast,by=list(ublast$Colors),FUN=mean)
blat<-aggregate(blat,by=list(blat$Colors),FUN=mean)
plast<-aggregate(plast,by=list(plast$Colors),FUN=mean)


# Get min and max number of colors
min_c<-min(blast$Colors,plast$Colors,mmseqs$Colors,ublast$Colors,blat$Colors)
max_c<-max(blast$Colors,plast$Colors,mmseqs$Colors,ublast$Colors,blat$Colors)

# Selection of samples for wide range plot, i.e., remove everything between min_c and 500 colors
blast_big<-subset(blast,blast$Colors>=500 | blast$Colors==min_c)
plast_big<-subset(plast,plast$Colors>=500 | plast$Colors==min_c)
mmseqs_big=subset(mmseqs,mmseqs$Colors>=500 | mmseqs$Colors==min_c)
ublast_big=subset(ublast,ublast$Colors>=500 | ublast$Colors==min_c)
blat_big=subset(blat,blat$Colors>=500 | blat$Colors==min_c)

# Selection of samples for small range plot, i.e., remove everything between larger than 500 colors
blast_small<-subset(blast,blast$Colors<=500)
plast_small<-subset(plast,plast$Colors<=500)
mmseqs_small=subset(mmseqs,mmseqs$Colors<=500)
ublast_small=subset(ublast,ublast$Colors<=500)
blat_small=subset(blat,blat$Colors<=500)

#get max runtime not considering outlier UBLAST and round a bit
max_t_small<-max(blast_small$UserTime,plast_small$UserTime,mmseqs_small$UserTime,blat_small$UserTime)
max_t_y_small<-ceiling(max_t_small/20)*20
max_t_y<-max_t_y_small*10

#get max memory not considering outlier MMseqs2 and round a bit
max_m_small<-max(blast_small$Memory,plast_small$Memory,ublast_small$Memory,blat_small$Memory)
max_m_y_small<-ceiling(max_m_small/1024/100)*100
max_m_y<-ceiling(max_m_y_small*10/1024)


### Plot runtimes for min_c, 500..max_c

pdf(file="runtime.pdf")

plot(blast_big$Colors,blast_big$UserTime,xlab="Number of genomes",ylab="Run time (seconds)",type="b",xlim=c(0,max_c),ylim=c(0,max_t_y),pch=pch_list[1])
points(plast_big$Colors,plast_big$UserTime,type="b",pch=pch_list[2])
points(mmseqs_big$Colors,mmseqs_big$UserTime,type="b",pch=pch_list[6])
points(ublast_big$Colors,ublast_big$UserTime,type="b",pch=pch_list[7])
points(blat_big$Colors,blat_big$UserTime,type="b",pch=pch_list[8])

# add rectangle
rect(xleft=-20, ybottom=-5, xright=520, ytop=max_t_y_small+5, col=NULL, border=par("fg"), lty="dashed", lwd=par("lwd"), xpd=FALSE)

# add legend
# legend("topleft",inset=0.02, c("PLAST","BLAST","BLAT","MMseqs2", "UBLAST"), pch=pch_list[c(2,1,8,6,7)])

dev.off()


### Plot runtimes for min_c..500

pdf(file="runtime_small.pdf")
plot(blast$Colors,blast$UserTime,xlab="Number of genomes",ylab="Run time (seconds)",type="b",xlim=c(0,500),ylim=c(0,max_t_y_small),frame.plot=FALSE,pch=pch_list[1])
points(plast$Colors,plast$UserTime,type="b",pch=pch_list[2])
points(mmseqs$Colors,mmseqs$UserTime,type="b",pch=pch_list[6])
points(ublast$Colors,ublast$UserTime,type="b",pch=pch_list[7])
points(blat$Colors,blat$UserTime,type="b",pch=pch_list[8])

# add dotted margin
box(which = "plot", lty = "dashed")

# add legend
legend("topleft",inset=0.02, c("PLAST","BLAST","BLAT","MMseqs2", "UBLAST"), pch=pch_list[c(2,1,8,6,7)])

dev.off()


### Plot memory for min_c, 500..max_c

pdf(file="memory.pdf")
plot(blast_big$Colors,blast_big$Memory/1024^2,xlab="Number of genomes",ylab="Memory peak (Gb)",type="b",xlim=c(0,max_c),ylim=c(0,max_m_y),pch=pch_list[1])
points(plast_big$Colors,plast_big$Memory/1024^2,type="b",pch=pch_list[2])
points(mmseqs_big$Colors,mmseqs_big$Memory/1024^2,type="b",pch=pch_list[6])
points(ublast_big$Colors,ublast_big$Memory/1024^2,type="b",pch=pch_list[7])
points(blat_big$Colors,blat_big$Memory/1024^2,type="b",pch=pch_list[8])

# add rectangle
rect(xleft=-20, ybottom=-0.1, xright=520, ytop=max_m_y_small/1024, col=NULL, border=par("fg"), lty="dashed", lwd=par("lwd"), xpd=FALSE)

# add legend
# legend("topleft",inset=0.02, c("PLAST","BLAST","BLAT","MMseqs2", "UBLAST"), pch=pch_list[c(2,1,8,6,7)])

dev.off()


### Plot memory for min_c..500

pdf(file="memory_small.pdf")
plot(blast$Colors,blast$Memory/1024,xlab="Number of genomes",ylab="Memory peak (Mb)",type="b",xlim=c(0,500),ylim=c(0,max_m_y_small),frame.plot=FALSE,pch=pch_list[1])
points(plast$Colors,plast$Memory/1024,type="b",pch=pch_list[2])
points(mmseqs$Colors,mmseqs$Memory/1024,type="b",pch=pch_list[6])
points(ublast$Colors,ublast$Memory/1024,type="b",pch=pch_list[7])
points(blat$Colors,blat$Memory/1024,type="b",pch=pch_list[8])

# add dotted margin
box(which = "plot", lty = "dashed")

# add legend
# legend("topleft",inset=0.02, c("PLAST","BLAST","BLAT","MMseqs2", "UBLAST"), pch=pch_list[c(2,1,8,6,7)])

dev.off()
