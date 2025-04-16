library("hadron")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args)==1)
mode <- args[1]
stopifnot(mode=="DG" || mode=="DM" || mode=="DM2")

V_cs <- 0.975
V_cd <- 0.221
V_su <- V_cd
G_F <- 1.1663788e-5
m_Ds <- 1.96835 # gev


if(mode=="DG") {
  zlist <- c(0, 1, 2, 3)
  zlistint <- c(0, 1, 2)
  conversionfactor <- 1/(48*pi^4/m_Ds^3) * G_F^2
  ylab <- "48 pi^4 / m_Ds^3 DG/Dq^2"
  normalise <- F
}

if(mode=="DM") {
  zlist <- c(0, 1, 2, 3, 4)
  zlistint <- c(0, 1, 2, 3)
  conversionfactor <- 1/(96*pi^4/m_Ds^3) * m_Ds^1 * G_F^2
  ylab <- "96 pi^4 / m_Ds^3 DM/Dq^2"
  normalise <- T
}

if(mode=="DM2") {
  zlist <- c(0, 1, 2, 3, 4, 5)
  zlistint <- c(0, 1, 2, 3, 4)
  conversionfactor <- 1/(960*pi^4/m_Ds^3) * m_Ds^2 * G_F^2
  ylab <- "96 pi^4 / m_Ds^3 DM^(2)/Dq^2"
  normalise <- T
}
savefolder <- "tables_fnfour_20"

onlypositive <- T


pdf(sprintf("plots/comparelimits_%s.pdf", mode), title="")
VAE <- read.table(sprintf("%s/%s_VAE_epslim_aic.csv", savefolder, mode), header=TRUE)
VAEcomb <- read.table(sprintf("%s/%s_VAE_epslim_aic_combined.csv", savefolder, mode), header=TRUE)

diff <- data.frame(channel=c(), kernel=c(), errtype=c(), iz=c(), reldiff=c(), errdiff=c())
for(errtype in c("stat", "sys", "vol", "tot")) {
  for(channel in c("cd", "cs")) {
    for(kernel in c("sigmoid", "erf")) {
      VEA <- read.table(sprintf("%s/%s_VEA_contlim_%s_%s_aic.csv", savefolder, mode, channel, kernel), header=TRUE)
      maskVAE <- VAE$errtype==errtype & VAE$channel==channel & VAE$kernel==kernel
      maskVAEcomb <- VAEcomb$errtype==errtype & VAEcomb$channel==channel & VAEcomb$kernel=="combined"
      maskVEA <- VEA$errtype==errtype
      plotwitherror(x=(VAE$theta[maskVAE]*0.0934516)^2, y=VAE$DGDq2[maskVAE]*conversionfactor, dy=VAE$dDGDq2[maskVAE]*conversionfactor, 
                    col=VAE$iz[maskVAE]+1, pch=1, xlab="q^2", ylab=ylab, 
                    main=paste("DG", channel, kernel, "errtype", errtype), 
                    xlim=c(0, 0.9))
      plotwitherror(x=VEA$q[maskVEA]^2+0.01, y=VEA$lim[maskVEA]*conversionfactor, dy=VEA$dlim[maskVEA]*conversionfactor, 
                    col=VEA$iz[maskVEA]+1, pch=2, rep=T)
      plotwitherror(x=(VAEcomb$theta[maskVAEcomb]*0.0934516)^2+0.02, y=VAEcomb$DGDq2[maskVAEcomb]*conversionfactor, dy=VAEcomb$dDGDq2[maskVAEcomb]*conversionfactor, 
                    col=VAEcomb$iz[maskVAEcomb]+1, pch=4, rep=T)
      legend(x="topright", legend=c("a, eps", "eps, a", "a, eps, comb", paste("Z", zlist)), col=c(1, 1, 1, zlist+1), pch=c(1, 2, 4, rep(3, length(zlist))))
      
      for(iz in zlist) {
        diff <- rbind(diff, data.frame(channel=channel, kernel=kernel, errtype=errtype, iz=iz, 
                                       reldiff=c((VAE$DGDq2[maskVAE & VAE$iz==iz]*conversionfactor - VEA$lim[maskVEA & VEA$iz==iz]) / sqrt(VAE$dDGDq2[maskVAE & VAE$iz==iz]^2*conversionfactor^2 + VEA$dlim[maskVEA & VEA$iz==iz]^2)), 
                                       errdiff=c(VAE$dDGDq2[maskVAE & VAE$iz==iz]*conversionfactor - VEA$dlim[maskVEA & VEA$iz==iz])))
      }
    }
  }
}

write.table(x=diff, file=sprintf("%s/%s_diff_VAE_VEA_afterlimits.csv", savefolder, mode), row.names=F, col.names=T)
dev.off()


errlist <- c("stat", "sys", "vol", "tot")
res <- data.frame(kernel=c(), errtype=c(), order=c(), total=c(), dtotal=c(), totalbs=c(), na=c())
reslist <- list()
if(normalise) normalisation <- readRDS(sprintf("%s/DG_final_result.RDS", savefolder))


for (errtype in errlist) {
  VEAerf         <- rep(0, 1001)
  VEAsigmoid     <- rep(0, 1001)
  VAEcomb        <- rep(0, 1001)
  suB64erf       <- rep(0, 1001)
  suB64sigmoid   <- rep(0, 1001)
  suD96erf       <- rep(0, 1001)
  suD96sigmoid   <- rep(0, 1001)
  
  ## cd
  VEAerfdata     <- readRDS(sprintf("%s/%s_VEA_%s_%s_int_%s.RDS", savefolder, mode, "cd", "erf", "aic"))
  VEAsigmoiddata <- readRDS(sprintf("%s/%s_VEA_%s_%s_int_%s.RDS", savefolder, mode, "cd", "sigmoid", "aic"))
  VAEcombdata    <- readRDS(sprintf("%s/%s_VAE_int_aic_combined.RDS", savefolder, mode))
  VEAerftab      <- read.table(sprintf("%s/%s_VEA_%s_%s_int_%s.csv", savefolder, mode, "cd", "erf", "aic"), header=T)
  VEAsigmoidtab  <- read.table(sprintf("%s/%s_VEA_%s_%s_int_%s.csv", savefolder, mode, "cd", "sigmoid", "aic"), header=T)
  VAEcombtab     <- read.table(sprintf("%s/%s_VAE_int_aic_combined.csv", savefolder, mode), header=T)
  for(iz in zlistint) {
    VEAerf     <- VEAerf     + c(VEAerftab$intspline[VEAerftab$name=="contlim" & VEAerftab$iz==iz & VEAerftab$errtype==errtype], VEAerfdata[[paste0("contlim", errtype, "iz", iz, "bsintspline")]]) * V_cd^2
    VEAsigmoid <- VEAsigmoid + c(VEAsigmoidtab$intspline[VEAsigmoidtab$name=="contlim" & VEAsigmoidtab$iz==iz & VEAsigmoidtab$errtype==errtype], VEAsigmoiddata[[paste0("contlim", errtype, "iz", iz, "bsintspline")]]) * V_cd^2
    VAEcomb    <- VAEcomb    + c(VAEcombtab$intspline[VAEcombtab$channel=="cd" & VAEcombtab$iz==iz & VAEcombtab$errtype==errtype], VAEcombdata[[paste0("cd", "combined", errtype, "iz", iz, "bsintspline")]]) * V_cd^2
  }
  
  ## cs
  VEAerfdata     <- readRDS(sprintf("%s/%s_VEA_%s_%s_int_%s.RDS", savefolder, mode, "cs", "erf", "aic"))
  VEAsigmoiddata <- readRDS(sprintf("%s/%s_VEA_%s_%s_int_%s.RDS", savefolder, mode, "cs", "sigmoid", "aic"))
  VAEcombdata    <- readRDS(sprintf("%s/%s_VAE_int_aic_combined.RDS", savefolder, mode))
  VEAerftab      <- read.table(sprintf("%s/%s_VEA_%s_%s_int_%s.csv", savefolder, mode, "cs", "erf", "aic"), header=T)
  VEAsigmoidtab  <- read.table(sprintf("%s/%s_VEA_%s_%s_int_%s.csv", savefolder, mode, "cs", "sigmoid", "aic"), header=T)
  VAEcombtab     <- read.table(sprintf("%s/%s_VAE_int_aic_combined.csv", savefolder, mode), header=T)
  for(iz in zlistint) {
    VEAerf     <- VEAerf     + c(VEAerftab$intspline[VEAerftab$name=="contlim" & VEAerftab$iz==iz & VEAerftab$errtype==errtype], VEAerfdata[[paste0("contlim", errtype, "iz", iz, "bsintspline")]]) * V_cs^2
    VEAsigmoid <- VEAsigmoid + c(VEAsigmoidtab$intspline[VEAsigmoidtab$name=="contlim" & VEAsigmoidtab$iz==iz & VEAsigmoidtab$errtype==errtype], VEAsigmoiddata[[paste0("contlim", errtype, "iz", iz, "bsintspline")]]) * V_cs^2
    VAEcomb    <- VAEcomb    + c(VAEcombtab$intspline[VAEcombtab$channel=="cs" & VAEcombtab$iz==iz & VAEcombtab$errtype==errtype], VAEcombdata[[paste0("cs", "combined", errtype, "iz", iz, "bsintspline")]]) * V_cs^2
  }
  
  ## su
  VEAsuerfdata     <- readRDS(sprintf("%s/%s_VEA_%s_%s_int_%s.RDS", savefolder, mode, "su", "erf", "su"))
  VEAsusigmoiddata <- readRDS(sprintf("%s/%s_VEA_%s_%s_int_%s.RDS", savefolder, mode, "su", "sigmoid", "su"))
  VEAsuerftab      <- read.table(sprintf("%s/%s_VEA_%s_%s_int_%s.csv", savefolder, mode, "su", "erf", "su"), header=T)
  VEAsusigmoidtab  <- read.table(sprintf("%s/%s_VEA_%s_%s_int_%s.csv", savefolder, mode, "su", "sigmoid", "su"), header=T)
  for(iz in zlistint) {
    if((errtype=="stat" || errtype=="sys") && onlypositive) {
      if(VEAsuerftab$intspline[VEAsuerftab$name=="B64" & VEAsuerftab$iz==iz & VEAsuerftab$errtype==errtype] > 0)                 suB64erf     <- suB64erf     + c(VEAsuerftab$intspline[VEAsuerftab$name=="B64" & VEAsuerftab$iz==iz & VEAsuerftab$errtype==errtype], VEAsuerfdata[[paste0("B64", errtype, "iz", iz, "bsintspline")]]) * V_su^2
      if(VEAsusigmoidtab$intspline[VEAsusigmoidtab$name=="B64" & VEAsusigmoidtab$iz==iz & VEAsusigmoidtab$errtype==errtype] > 0) suB64sigmoid <- suB64sigmoid + c(VEAsusigmoidtab$intspline[VEAsusigmoidtab$name=="B64" & VEAsusigmoidtab$iz==iz & VEAsusigmoidtab$errtype==errtype], VEAsusigmoiddata[[paste0("B64", errtype, "iz", iz, "bsintspline")]]) * V_su^2
      if(VEAsuerftab$intspline[VEAsuerftab$name=="D96" & VEAsuerftab$iz==iz & VEAsuerftab$errtype==errtype] > 0)                 suD96erf     <- suD96erf     + c(VEAsuerftab$intspline[VEAsuerftab$name=="D96" & VEAsuerftab$iz==iz & VEAsuerftab$errtype==errtype], VEAsuerfdata[[paste0("D96", errtype, "iz", iz, "bsintspline")]]) * V_su^2
      if(VEAsusigmoidtab$intspline[VEAsusigmoidtab$name=="D96" & VEAsusigmoidtab$iz==iz & VEAsusigmoidtab$errtype==errtype] > 0) suD96sigmoid <- suD96sigmoid + c(VEAsusigmoidtab$intspline[VEAsusigmoidtab$name=="D96" & VEAsusigmoidtab$iz==iz & VEAsusigmoidtab$errtype==errtype], VEAsusigmoiddata[[paste0("D96", errtype, "iz", iz, "bsintspline")]]) * V_su^2
    } else if((errtype=="stat" || errtype=="sys") && !onlypositive) {
      suB64erf     <- suB64erf     + c(VEAsuerftab$intspline[VEAsuerftab$name=="B64" & VEAsuerftab$iz==iz & VEAsuerftab$errtype==errtype], VEAsuerfdata[[paste0("B64", errtype, "iz", iz, "bsintspline")]]) * V_su^2
      suB64sigmoid <- suB64sigmoid + c(VEAsusigmoidtab$intspline[VEAsusigmoidtab$name=="B64" & VEAsusigmoidtab$iz==iz & VEAsusigmoidtab$errtype==errtype], VEAsusigmoiddata[[paste0("B64", errtype, "iz", iz, "bsintspline")]]) * V_su^2
      suD96erf     <- suD96erf     + c(VEAsuerftab$intspline[VEAsuerftab$name=="D96" & VEAsuerftab$iz==iz & VEAsuerftab$errtype==errtype], VEAsuerfdata[[paste0("D96", errtype, "iz", iz, "bsintspline")]]) * V_su^2
      suD96sigmoid <- suD96sigmoid + c(VEAsusigmoidtab$intspline[VEAsusigmoidtab$name=="D96" & VEAsusigmoidtab$iz==iz & VEAsusigmoidtab$errtype==errtype], VEAsusigmoiddata[[paste0("D96", errtype, "iz", iz, "bsintspline")]]) * V_su^2
    }
  }
  
  ## normalisation DM, DM2
  if(normalise) {
    VEAerf       <- VEAerf       / normalisation[[paste0(errtype, "VEAerf")]]
    VEAsigmoid   <- VEAsigmoid   / normalisation[[paste0(errtype, "VEAsigmoid")]]
    VAEcomb      <- VAEcomb      / normalisation[[paste0(errtype, "VAEcomb")]]
    suB64erf     <- suB64erf     / normalisation[[paste0(errtype, "VEAerf")]]
    suB64sigmoid <- suB64sigmoid / normalisation[[paste0(errtype, "VEAsigmoid")]]
    suD96erf     <- suD96erf     / normalisation[[paste0(errtype, "VEAerf")]]
    suD96sigmoid <- suD96sigmoid / normalisation[[paste0(errtype, "VEAsigmoid")]]
  }
  
  ## save results
  res <- rbind(res, data.frame(kernel="erf", errtype=errtype, order="VEA", total=VEAerf[1], dtotal=sd(VEAerf[2:1001], na.rm=T), totalbs=mean(VEAerf[2:1001], na.rm=T), na=length(which(is.na(VEAerf)))))
  res <- rbind(res, data.frame(kernel="sigmoid", errtype=errtype, order="VEA", total=VEAsigmoid[1], dtotal=sd(VEAsigmoid[2:1001], na.rm=T), totalbs=mean(VEAsigmoid[2:1001], na.rm=T), na=length(which(is.na(VEAsigmoid)))))
  res <- rbind(res, data.frame(kernel="comb", errtype=errtype, order="VAE", total=VAEcomb[1], dtotal=sd(VAEcomb[2:1001], na.rm=T), totalbs=mean(VAEcomb[2:1001], na.rm=T), na=length(which(is.na(VAEcomb)))))
  res <- rbind(res, data.frame(kernel="erf", errtype=errtype, order="VEA_su_B64", total=suB64erf[1], dtotal=sd(suB64erf[2:1001], na.rm=T), totalbs=mean(suB64erf[2:1001], na.rm=T), na=length(which(is.na(suB64erf)))))
  res <- rbind(res, data.frame(kernel="sigmoid", errtype=errtype, order="VEA_su_B64", total=suB64sigmoid[1], dtotal=sd(suB64sigmoid[2:1001], na.rm=T), totalbs=mean(suB64sigmoid[2:1001], na.rm=T), na=length(which(is.na(suB64sigmoid)))))
  res <- rbind(res, data.frame(kernel="erf", errtype=errtype, order="VEA_su_D96", total=suD96erf[1], dtotal=sd(suD96erf[2:1001], na.rm=T), totalbs=mean(suD96erf[2:1001], na.rm=T), na=length(which(is.na(suD96erf)))))
  res <- rbind(res, data.frame(kernel="sigmoid", errtype=errtype, order="VEA_su_D96", total=suD96sigmoid[1], dtotal=sd(suD96sigmoid[2:1001], na.rm=T), totalbs=mean(suD96sigmoid[2:1001], na.rm=T), na=length(which(is.na(suD96sigmoid)))))
  reslist[[paste0(errtype, "VEAerf")]] <- VEAerf  * conversionfactor
  reslist[[paste0(errtype, "VEAsigmoid")]] <- VEAsigmoid  * conversionfactor
  reslist[[paste0(errtype, "VAEcomb")]] <- VAEcomb  * conversionfactor
}

res$total <- res$total * conversionfactor
res$dtotal <- res$dtotal * conversionfactor
res$totalbs <- res$totalbs * conversionfactor
reslist <- reslist

res
finalcomment <- ""
if(normalise) finalcomment <- "_normalised"

write.table(res, sprintf("%s/%s_final_result%s.csv", savefolder, mode, finalcomment))
saveRDS(reslist, sprintf("%s/%s_final_result%s.RDS", savefolder, mode, finalcomment))
warnings()
