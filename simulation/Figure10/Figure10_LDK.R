#Figure10_LDK

SK <- as.matrix(read.table("SK.csv",header = TRUE,sep=","))
LDK4 <- as.matrix(read.table("LDK4.csv",header = TRUE,sep=","))
LDK9 <- as.matrix(read.table("LDK9.csv",header = TRUE,sep=","))
LDK16 <- as.matrix(read.table("LDK16.csv",header = TRUE,sep=","))
LDK25 <- as.matrix(read.table("LDK25.csv",header = TRUE,sep=","))
LDK36 <- as.matrix(read.table("LDK36.csv",header = TRUE,sep=","))
LDK49 <- as.matrix(read.table("LDK49.csv",header = TRUE,sep=","))
LDK64 <- as.matrix(read.table("LDK64.csv",header = TRUE,sep=","))
LDK81 <- as.matrix(read.table("LDK81.csv",header = TRUE,sep=","))
LDK100 <- as.matrix(read.table("LDK100.csv",header = TRUE,sep=","))

timeSK <- as.matrix(read.table("timeSK.csv",header = TRUE,sep=","))
timeLDK4 <- as.matrix(read.table("timeLDK4.csv",header = TRUE,sep=","))
timeLDK9 <- as.matrix(read.table("timeLDK9.csv",header = TRUE,sep=","))
timeLDK16 <- as.matrix(read.table("timeLDK16.csv",header = TRUE,sep=","))
timeLDK25 <- as.matrix(read.table("timeLDK25.csv",header = TRUE,sep=","))
timeLDK36 <- as.matrix(read.table("timeLDK36.csv",header = TRUE,sep=","))
timeLDK49 <- as.matrix(read.table("timeLDK49.csv",header = TRUE,sep=","))
timeLDK64 <- as.matrix(read.table("timeLDK64.csv",header = TRUE,sep=","))
timeLDK81 <- as.matrix(read.table("timeLDK81.csv",header = TRUE,sep=","))
timeLDK100 <- as.matrix(read.table("timeLDK100.csv",header = TRUE,sep=","))

x <- c(1,4,9,16,25,36,49,64,81,100)

df5 <- cbind(SK,LDK4,LDK9,LDK16, LDK25, LDK36, LDK49, LDK64, LDK81, LDK100)
df6 <- cbind(timeSK[,1],timeLDK4[,1],timeLDK9[,1],timeLDK16[,1], timeLDK25[,1], timeLDK36[,1], timeLDK49[,1], timeLDK64[,1], timeLDK81[,1], timeLDK100[,1])

filename <- "LDKbp.jpeg" 
jpeg(filename, width=1024, height=1024, res = 300)
par(oma=c(1,1,.1,.1), mar=c(2.5,1,1,.5))
boxplot(df5, xaxt="n", yaxt="n")
axis(2, at=c(0.105,0.115,0.125,0.135), labels=c(0.105,0.115,0.125,0.135))
axis(1, at=1:10, labels=c("SK",paste0("M=",x[-1])),cex.axis=0.9, las=2)
dev.off()

filename <- "LDKbp_time.jpeg" 
jpeg(filename, width=1024, height=1024, res = 300)
par(oma=c(1,1,.1,.1), mar=c(2.5,1,1,.5))
boxplot(df6, xaxt="n", yaxt="n")
axis(2, at=c(50,100,150,200,250), labels=c(50,100,150,200,250))
axis(1, at=1:10, labels=c("SK",paste0("M=",x[-1])),cex.axis=0.9, las=2)
dev.off()

