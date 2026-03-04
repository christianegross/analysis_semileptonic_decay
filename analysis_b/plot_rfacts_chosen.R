library("hadron")
library("optparse")
library("christianesfunctions")


if (TRUE) {
  # set option list
  option_list <- list(
    make_option(c("-i", "--inputname"), type = "character", default = "",
                help = "basename of the table of chosen rfacts [default %default]"),
    make_option(c("-b", "--basename"), type = "character", default = "",
                help = "base name for the output pdf [default %default]"),
    make_option(c("-m", "--mode"), type = "character", default = "DG",
                help = "mode: DG, DM, DM2 [default %default]"),
    make_option(c("-k", "--kernel"), type = "character", default = "sigmoid",
                help = "kernel [default %default]"),
    make_option(c("-f", "--folder"), type = "character", default = "-1",
                help = "folder with results [default %default]"),
    make_option(c("-n", "--normnumber"), type = "integer", default = 0,
                help = "index of norm to be analysed [default %default]")#,
    #~     make_option(c("-b", "--bmass"), type = "integer", default = "-1",
    #~     help = "index for bmass [default %default]")
  )
  parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
  args <- parse_args(parser, positional_arguments = 0)
  opt <- args$options
}

if(opt$folder=="-1") {
  opt$folder <- sprintf("~/Documents/heavymesons/data/%s/%s/%s", opt$channel, opt$ensemble, opt$momentum)
}
if(opt$basename=="") opt$basename <- opt$input


if(opt$mode=="DG") {
  modefolder <- "DGammaDq2"
  zmax <- 2
} else if(opt$mode=="DM") {
  modefolder <- "DMDq2"
  zmax <- 3
} else if(opt$mode=="DM2") {
  modefolder <- "DM2Dq2"
  zmax <- 4
}


erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

chosen <- read.table(opt$input, header=T, sep=",")
if(length(chosen$ieps)==0) chosen <- read.table(opt$input, header=T, sep=" ")
if(length(chosen$ieps)==0) stop("cannot read in resultfile")

res <- data.frame(ik=c(), spectreflag=c(), lambda_start=c(), lambdalambda_start=c(), Bnorm=c(), 
                  A0=c(), AA0_min=c(), AA0_ref=c(), AA0=c(), BnormB_ref=c(), BnormB=c(), 
                  C_ref=c(), C=c(), rho=c(), drho_stat=c(), drho_syst=c(), drho_tot=c(), resflag=c(), 
                  rfact=c(), iz=c(), ieps=c(), eps=c())
ressyserr <- data.frame(ik=c(), spectreflag=c(), lambda_start=c(), lambdalambda_start=c(), Bnorm=c(), 
                        A0=c(), AA0_min=c(), AA0_ref=c(), AA0=c(), BnormB_ref=c(), BnormB=c(), 
                        C_ref=c(), C=c(), rho=c(), drho_stat=c(), drho_syst=c(), drho_tot=c(), resflag=c(), 
                        rfact=c(), iz=c(), ieps=c(), eps=c())

for (iz in 0:zmax) {
  for(ieps in unique(chosen$ieps)) {
    
    rfact <- chosen$statrfact[chosen$iz==iz & chosen$ieps==ieps]
    filename <- sprintf("%s/%s_%.2e/output%s/%s_%d_iset_0_ieps_%d_icomb_0.dat", opt$folder, opt$kernel, 10^rfact, modefolder, opt$mode, iz, ieps)
    datastab <- try(read.table(filename, fill=TRUE, skip=3,
                               col.names=c("ik", "spectreflag", "lambda_start",
                                           "lambdalambda_start", "Bnorm", "A0", "AA0_min",
                                           "AA0_ref", "AA0", "BnormB_ref", "BnormB",
                                           "C_ref", "C", "rho", "drho_stat", "drho_syst",
                                           "drho_tot", "resflag")), silent=T)
    if(!inherits(datastab, "try-error")) {
      #~                     datastabaa0 <- datastab[datastab$spectreflag==1 & datastab$resflag==1 & datastab$ik, ]
      datastabaa0 <- datastab[datastab$spectreflag==1 & datastab$ik==opt$n, ]
      datastabaa0 <- datastabaa0[length(datastabaa0$ik), ]
      #~                     datastabaa0$rfact <- rfact
      datastabaa0 <- cbind(datastabaa0, data.frame(rfact=rfact, iz=iz, ieps=ieps, eps=chosen$eps[chosen$iz==iz & chosen$ieps==ieps]))
#~       print(datastabaa0)
      res <- rbind(res, datastabaa0)
      
    }
    filenamesyserr <- sprintf("%s/%s_%.2e/output%s/%s_%d_iset_0_ieps_%d_icomb_0.dat", opt$folder, opt$kernel, 10^(rfact-1), modefolder, opt$mode, iz, ieps)
    datastabsyserr <- try(read.table(filenamesyserr, fill=TRUE, skip=3,
                                     col.names=c("ik", "spectreflag", "lambda_start",
                                                 "lambdalambda_start", "Bnorm", "A0", "AA0_min",
                                                 "AA0_ref", "AA0", "BnormB_ref", "BnormB",
                                                 "C_ref", "C", "rho", "drho_stat", "drho_syst",
                                                 "drho_tot", "resflag")), silent=T)
    if(!inherits(datastabsyserr, "try-error")) {
      datastabaa0syserr <- datastabsyserr[datastabsyserr$spectreflag==1 & datastabsyserr$resflag==1, ]
      datastabaa0syserr <- datastabsyserr[datastabsyserr$spectreflag==1 & datastabsyserr$ik==opt$n, ]
      datastabaa0syserr <- datastabaa0syserr[length(datastabaa0syserr$ik), ]
      #~                     datastabaa0$rfact <- rfact
      datastabaa0syserr <- cbind(datastabaa0syserr, data.frame(rfact=rfact, iz=iz, ieps=ieps, eps=chosen$eps[chosen$iz==iz & chosen$ieps==ieps]))
      ressyserr <- rbind(ressyserr, datastabaa0syserr)
      
    }
  }
}


stopifnot(dim(ressyserr)==dim(res))

res$pull <- (res$rho-ressyserr$rho)/res$drho_stat
res$dsys <- abs(res$rho-ressyserr$rho)*erf(abs(res$pull)/sqrt(2))
seed <- paste0(charToRaw(paste0(as.character(opt$momentum), as.character(opt$bmass), "s", collapse="")), collapse="") ## convert characteristics of system to long string
  set.seed(as.integer(as.numeric(seed)%%131071)) ## convert string to numeric, reduce number of digits by taking modulo the 17th Mersenne Prime



#~ print(res)

pdf(paste0(opt$basename, ".pdf"), title="")
for(iz in 0:zmax) {
  plotwitherror(x=res$eps[res$iz==iz], y=res$rho[res$iz==iz], dy=res$drho_stat[res$iz==iz],
                ylab="rho", xlab="eps/m_H", main=paste("Z", iz))
}
dev.off()

print(paste0(opt$basename, ifelse(opt$basename==opt$input, "_chosen_data", ""), ".csv"))
write.table(res, paste0(opt$basename, ifelse(opt$basename==opt$input, "_chosen_data", ""), ".csv"), row.names=F)

boots <- chosenDGarray(basepath=opt$folder, tablechosen=paste0(opt$basename, ifelse(opt$basename==opt$input, "_chosen_data", ""), ".csv"), ik=opt$normnumber)[1, , , ]
saveRDS(boots, file=paste0(opt$basename, "_boots_err_stat.RDS"))
boots_sys <- addsyshltDGarray(boots, paste0(opt$basename, ifelse(opt$basename==opt$input, "_chosen_data", ""), ".csv"))
saveRDS(boots_sys, file=paste0(opt$basename, "_boots_err_stat_HLT.RDS"))



pdf(paste0(opt$basename, "_kernel.pdf"), title="")
for (iz in 0:zmax) {
  for(ieps in unique(chosen$ieps)) {
    
    rfact <- chosen$statrfact[chosen$iz==iz & chosen$ieps==ieps]
    eps <- chosen$eps[chosen$iz==iz & chosen$ieps==ieps]
    filename <- sprintf("%s/%s_%.2e/output%s/%s_%d_iset_0_ieps_%d_icomb_0.dat", opt$folder, opt$kernel, 10^rfact, modefolder, opt$mode, iz, ieps)
    data <- try(read.table(filename, fill=TRUE,
                           col.names=c("ik", "spectreflag", "omega",
                                       "kernel", "kernelbar", "kernelbarmkernel", "aM",
                                       paste0("column", 8:18))))
    if(!inherits(data, "try-error")) {
      data <- data[-seq(1, 9), ]
      data$ik <- as.integer(data$ik)
      data <- data[data$spectreflag == 0, ]
      plot(NA, xlim=c(min(data$omega), max(data$omega)),
           ylim=c(min(c(0, data$kernelbar)), max(data$kernelbar)),
           xlab="omega", ylab="kernel",
           main=paste("iz", iz, "ieps", ieps, "eps", eps, "rfact", rfact))
      
      lines(x=data$omega[data$ik == opt$normnumber],
            y=data$kernel[data$ik == opt$normnumber],
            col="black")
      
      lines(x=data$omega[data$ik == opt$normnumber],
            y=data$kernelbar[data$ik == opt$normnumber],
            col="red")
      legend(x="topright", col=c("black", "red"),
             lty=c(1, 1), legend=c("kernel", "reconstructed"))
    }
  }
}
