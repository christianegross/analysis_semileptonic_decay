library("hadron")
args <- commandArgs(trailingOnly = TRUE)

hlt <- try(read.table(sprintf("%s/outputY/Y.dat", args[1])))
Rcopyhlt <- try(read.table(sprintf("%s/outputR/Y.dat", args[1])))
Rfittedamplitude <- try(read.table(sprintf("%s/outputRfittedamplitude/Y.dat", args[1])))

if(inherits(hlt, "try-error") || inherits(Rcopyhlt, "try-error") || inherits(Rfittedamplitude, "try-error")) {
  stop("there was a problem reading in the data")
}
pdf(sprintf("%s/compareZ.pdf", args[1]), title="")

for(iy in 1:5) {
  plotwitherror(x=hlt$V6[hlt$V5==iy], y=abs(hlt$V7[hlt$V5==iy]), dy=hlt$V8[hlt$V5==iy], 
                col="black", pch=0, xlab="t/a", ylab="Y", main=paste("Y", iy), log="y", 
                ylim=range(abs(hlt$V7[hlt$V5==iy])))
  plotwitherror(x=Rcopyhlt$V6[Rcopyhlt$V5==iy], y=abs(Rcopyhlt$V7[Rcopyhlt$V5==iy]), dy=Rcopyhlt$V8[Rcopyhlt$V5==iy], 
                col="red", pch=1, rep=T)
  plotwitherror(x=Rfittedamplitude$V6[Rfittedamplitude$V5==iy], y=abs(Rfittedamplitude$V7[Rfittedamplitude$V5==iy]), dy=Rfittedamplitude$V8[Rfittedamplitude$V5==iy], 
                col="blue", pch=2, rep=T)
  legend(x="topright", legend=c("hlt", "R hlt", "R fitted amplitude"), col=c("black", "red", "blue"), pch=0:2)
}
for(iy in 1:5) {
  plotwitherror(x=hlt$V6[hlt$V5==iy], y=hlt$V8[hlt$V5==iy], 
                col="black", pch=0, xlab="t/a", ylab="Y", main=paste("Y", iy, "error"), log="y")
  plotwitherror(x=Rcopyhlt$V6[Rcopyhlt$V5==iy], y=Rcopyhlt$V8[Rcopyhlt$V5==iy], 
                col="red", pch=1, rep=T)
  plotwitherror(x=Rfittedamplitude$V6[Rfittedamplitude$V5==iy], y=Rfittedamplitude$V8[Rfittedamplitude$V5==iy], 
                col="blue", pch=2, rep=T)
  legend(x="topright", legend=c("hlt", "R hlt", "R fitted amplitude"), col=c("black", "red", "blue"), pch=0:2)
}
for(iy in 1:5) {
  plotwitherror(x=hlt$V6[hlt$V5==iy], y=Rcopyhlt$V7[Rcopyhlt$V5==iy]/hlt$V7[hlt$V5==iy], 
                col="black", pch=0, xlab="t/a", ylab="Y", main=paste("Y", iy, "ratio to hlt result"), ylim=c(0.99, 1.01))
  plotwitherror(x=hlt$V6[hlt$V5==iy], y=Rfittedamplitude$V7[Rfittedamplitude$V5==iy]/hlt$V7[hlt$V5==iy], 
                col="red", pch=1, rep=T)
  legend(x="topright", legend=c("R hlt", "R fitted amplitude"), col=c("black", "red", "blue"), pch=0:2)
}


hlt <- try(read.table(sprintf("%s/outputY/Y.log", args[1]), skip=36+128))
Rcopyhlt <- try(read.table(sprintf("%s/outputR/Y.log", args[1]), skip=20, nrows = 64))
# Rfittedamplitude <- try(read.table(sprintf("%s/outputRfittedamplitude/Y.log", args[1]), skip=20, nrows = 20))

# print(hlt)
# print(Rcopyhlt)
# print(Rfittedamplitude)

hlt <- hlt[2:65, ]

plotwitherror(x=hlt$V2, y=abs(hlt$V5), dy=hlt$V6, 
              col="black", pch=0, xlab="t/a", ylab="Y", main=paste("masses"), ylim=range(hlt$V5[hlt$V5>0])
)
plotwitherror(x=Rcopyhlt$V2, y=abs(Rcopyhlt$V3), dy=Rcopyhlt$V4,
              col="red", pch=1, rep=T)
# plotwitherror(x=Rfittedamplitude$V2, y=abs(Rfittedamplitude$V3), dy=Rfittedamplitude$V4, 
#               col="blue", pch=2, rep=T)
legend(x="topright", legend=c("hlt", "R hlt", "R fitted amplitude"), col=c("black", "red", "blue"), pch=0:2)


plotwitherror(x=hlt$V2, y=abs(Rcopyhlt$V3/hlt$V5-1), 
              col="black", pch=0, xlab="t/a", ylab="Y", main=paste("abs(ratio masses-1)"),
              log="y"
)
legend(x="topright", legend=c("R hlt", "R fitted amplitude"), col=c("black", "red"), pch=0:2)

