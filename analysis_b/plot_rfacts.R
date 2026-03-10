library("hadron")
library("optparse")


if (TRUE) {
  # set option list
  option_list <- list(
    make_option(c("-c", "--channel"), type = "character", default = "cs",
                help = "decay channel [default %default]"),
    make_option(c("-e", "--ensemble"), type = "character", default = "cB211.07.64",
                help = "ensemble [default %default]"),
    make_option(c("-t", "--momentum"), type = "character", default = "th2",
                help = "momentum [default %default]"),
    make_option(c("-m", "--mode"), type = "character", default = "DG",
                help = "mode: DG, DM, DM2 [default %default]"),
    make_option(c("-k", "--kernel"), type = "character", default = "sigmoid",
                help = "kernel [default %default]"),
    make_option(c("-f", "--folder"), type = "character", default = "-1",
                help = "folder with results [default %default]"),
    make_option(c("-p", "--plotfolder"), type = "character", default = "-1",
                help = "folder where plots and tables are stored [default %default]"),
    make_option(c("-b", "--bmass"), type = "integer", default = "-1",
                help = "index for bmass [default %default]"),
    make_option(c("-l", "--lowr"), type = "integer", default = "-2",
                help = "log10 of lowest analysed A/B factor [default %default]"),
    make_option(c("--plotlowr"), type = "integer", default = "3",
                help = "log10 of lowest analysed A/B factor to be shown in plotting [default %default]"),
    make_option(c("-u", "--highr"), type = "integer", default = "6",
                help = "log10 of highest analysed A/B factor [default %default]"),
    make_option(c("-n", "--normnumber"), type = "integer", default = 0,
                help = "index of norm to be analysed [default %default]"),
    make_option(c("--drawhorizontal"), action = "store_true", default = FALSE,
                help = "if true, draws horizontal lines for each endpoint [default %default]"),
    make_option(c("--epsmax"), type = "integer", default = 40,
                help = "index of highest smearing parameter to be considered [default %default]"),
    make_option(c("--numsuggest"), type = "integer", default = 2,
                help = "suggest the point with highest A/A0 that is compatible with the numsuggest points with smaller A/A0 [default %default]"),
    make_option(c("--compareerror"), type = "character", default = "sum",
                help = "Give the mode with which to calculate the compatibility matrix [default %default]")
    
  )
  parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
  args <- parse_args(parser, positional_arguments = 0)
  opt <- args$options
}

if(opt$folder=="-1") {
  opt$folder <- sprintf("~/Documents/heavymesons/data/%s/%s/%s", opt$channel, opt$ensemble, opt$momentum)
}
if(opt$plotfolder=="-1") {
  opt$plotfolder <- opt$folder
}


stopifnot(opt$highr > opt$lowr)

if(!(opt$compareerror %in% c("sum", "sqrtsum", "lower"))) stop("invalid mode for building the compatibility matrix given!")

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

rfacts <- 10^(opt$lowr:opt$highr)

bmassaddon <- ifelse(opt$bmass==-1, "", sprintf("_bmass_%d", opt$bmass))


pdf(sprintf("%s/rfact_%s_%s%s_ik%s_%s.pdf", opt$plotfolder, opt$mode, opt$kernel, bmassaddon, opt$normnumber, opt$momentum), title="", height=11.7, width=16.6)

filenamerfacts <- sprintf("%s/rfact_%s_%s%s_ik%d_%s.csv", opt$plotfolder, opt$mode, opt$kernel, bmassaddon, opt$normnumber, opt$momentum)

resultsselected <- file.exists(filenamerfacts)
if(resultsselected)  {
  rselection <- read.table(filenamerfacts, header=T, sep=",")
  if(length(rselection$ieps)==0) rselection <- read.table(filenamerfacts, header=T, sep=" ")
  if(length(rselection$ieps)==0) stop("cannot read in resultfile")
}

suggestedrfact <- c()
numrfactsuggested <- 0
numrfacttaken <- 0

## opaque gray

pcol <- col2rgb("gray", alpha=TRUE)/255 
pcol[4] <- 0.5
pcol <- rgb(red=pcol[1],green=pcol[2],blue=pcol[3],alpha=pcol[4])

for(iz in 0:zmax) {
  for(iepsilon in 0:min(opt$epsmax, 40)) {
    res <- data.frame(ik=c(), spectreflag=c(), lambda_start=c(), lambdalambda_start=c(), Bnorm=c(), A0=c(), AA0_min=c(), AA0_ref=c(), AA0=c(), BnormB_ref=c(), BnormB=c(), C_ref=c(), C=c(), rho=c(), drho_stat=c(), drho_syst=c(), drho_tot=c(), resflag=c(), rfact=c())
    additionalres <- data.frame(ik=c(), spectreflag=c(), lambda_start=c(), lambdalambda_start=c(), Bnorm=c(), A0=c(), AA0_min=c(), AA0_ref=c(), AA0=c(), BnormB_ref=c(), BnormB=c(), C_ref=c(), C=c(), rho=c(), drho_stat=c(), drho_syst=c(), drho_tot=c(), resflag=c())
    
    for(rfact in rfacts) {
      filename <- sprintf("%s/%s_%.2e_multi/output%s/%s_%d_iset_0_ieps_%d_icomb_0.dat", opt$folder, opt$kernel, rfact, modefolder, opt$mode, iz, iepsilon)
      datastab <- try(read.table(filename, fill=TRUE, skip=9,
                                 col.names=c("ik", "spectreflag", "lambda_start",
                                             "lambdalambda_start", "Bnorm", "A0", "AA0_min",
                                             "AA0_ref", "AA0", "BnormB_ref", "BnormB",
                                             "C_ref", "C", "rho", "drho_stat", "drho_syst",
                                             "drho_tot", "resflag")), silent=T)
      
      
      if(!inherits(datastab, "try-error")) {
        datastabaa0 <- datastab[datastab$spectreflag==1 & datastab$ik==opt$norm, ]
        additionalres <- rbind(additionalres, datastab[datastab$spectreflag==1 & datastab$ik!=opt$norm, ])
        datastabaa0$rfact <- rfact
        res <- rbind(res, datastabaa0)
        
      }
    }
    if(length(res$AA0_ref) > 1) {
      masklim <- res$resflag==1 & res$rfact>=10^opt$plotlowr & res$rfact<=10^opt$highr
      additionalres <- additionalres[order(additionalres$AA0_ref), ]
      
      ## plot additonal norms as grey points in background
      plot(NA, xlab="A/A0", ylab="rho", 
           main=paste(opt$channel, opt$ensemble, opt$momentum, opt$mode, "eps", iepsilon, "Z", iz, bmassaddon), 
           log="x", 
           ylim=range(res$rho[masklim] - 1.02*res$drho_stat[masklim], res$rho[masklim] + 1.02*res$drho_stat[masklim]),
           xlim=range(res$AA0_ref[res$resflag==0]))
      for(ik in unique(additionalres$ik)) {
      polygon(x=c(additionalres$AA0_ref[additionalres$resflag==0 & additionalres$ik==ik], 
                  rev(additionalres$AA0_ref[additionalres$resflag==0 & additionalres$ik==ik])),
              y=c(additionalres$rho[additionalres$resflag==0 & additionalres$ik==ik]+additionalres$drho_stat[additionalres$resflag==0 & additionalres$ik==ik], 
                  rev(additionalres$rho[additionalres$resflag==0 & additionalres$ik==ik]-additionalres$drho_stat[additionalres$resflag==0 & additionalres$ik==ik])), col=pcol, lty=0, lwd=0.001, border=pcol)
      }
      
      plotwitherror(x=res$AA0_ref[res$resflag==0], y=res$rho[res$resflag==0], dy=res$drho_stat[res$resflag==0], 
                    col=log10(res$rfact[res$resflag==0])+5, rep=T)
      grid()
      abline(v=outer(1:9, 10^((opt$lowr-1):(opt$highr+2))), col = "lightgray", lty = "dotted")
      abline(v=(res$AA0_ref)[res$resflag==1], col=log10(res$rfact[res$resflag==1])+5)
      abline(h=(res$rho)[res$resflag==1] + (res$drho_stat)[res$resflag==1], col=log10(res$rfact[res$resflag==1])+5, lty=4)
      abline(h=(res$rho)[res$resflag==1] - (res$drho_stat)[res$resflag==1], col=log10(res$rfact[res$resflag==1])+5, lty=4)
      
      legend(x="topright", legend=sprintf("%.2e", rfacts), col=log10(rfacts)+5, pch=c(21), pt.bg=log10(rfacts)+5)
      
      ## calculate suggested A/B fact
      ## select point that agrees with the next numsuggest smaller points at the level of the points with lower A/A0
      ## rfact ascends -> we have to consider the vector in reverse, higher rfact means higher A/A0
      respoints <- (res$rho)[res$resflag==1]
      drespoints <- (res$drho_stat)[res$resflag==1]
      
      if(opt$compareerror=="lower")   mat <- t(t(outer(respoints, respoints, '-'))/drespoints)
      if(opt$compareerror=="sqrtsum") mat <- t(t(outer(respoints, respoints, '-'))/outer(drespoints, drespoints, function(x, y) sqrt(x^2+y^2)))
      if(opt$compareerror=="sum")     mat <- t(t(outer(respoints, respoints, '-'))/outer(drespoints, drespoints, '+'))
      mat <- abs(mat)
      #~             print(mat)
      mat <- abs(mat) < 1
      #~             print(mat)
      foundsuggestion <- FALSE
      rfactindex <- opt$highr - opt$lowr + 2
      #~             print(rfactindex)
      while(!foundsuggestion & rfactindex > (opt$numsuggest+1)) {
        rfactindex <- rfactindex - 1
        foundsuggestion <- all(mat[rfactindex, (rfactindex - 1):(rfactindex - opt$numsuggest)])
      }
      if(!foundsuggestion) {
        print(paste(opt$ensemble, opt$momentum, opt$bmass, "iz", iz, "ieps", iepsilon, "no point can be suggested!"))
        suggestedrfact <- append(suggestedrfact, NA)
      }
      if(foundsuggestion) {
        abline(v=(res$AA0_ref)[res$resflag==1][rfactindex], col="red", lwd=5, lty=4)
        suggestedrfact <- append(suggestedrfact, log10(rfacts[rfactindex]))
      }
      
      if(resultsselected){
        rfactindexsolve <- rselection$statrfact[rselection$iz==iz & rselection$ieps==iepsilon] - opt$lowr + 1
        abline(v=(res$AA0_ref)[res$resflag==1][rfactindexsolve], lwd=10, lty=3, col="darkgreen")
        legend(x="topleft", legend=c("selected", "suggested"), col=c("darkgreen", "red"), lty=c(3, 4), lwd=5)
        if(foundsuggestion) {
          numrfacttaken <- numrfacttaken + as.integer(rselection$statrfact[rselection$iz==iz & rselection$ieps==iepsilon] == suggestedrfact[length(suggestedrfact)])
          numrfactsuggested <- numrfactsuggested + 1
          
        }
      } else {
        legend(x="topleft", legend=c("suggested"), col=c("red"), lty=c(4), lwd=5)
      }
    }
    
  }
}


#~ print(sprintf("%s/%s_%.2e/output%s/%s.dat", opt$folder, opt$kernel, rfacts[1], modefolder, modefolder))
tmp <- read.table(sprintf("%s/%s_%.2e/output%s/%s.dat", opt$folder, opt$kernel, rfacts[1], modefolder, modefolder), col.names=c("mH", "iset", "w", "icomb", "idg", "ieps", "eps", "rho", "stat", "sys", "tot"))
#~ print(tmp)
tmp <- tmp[tmp$icomb==0 & tmp$idg <= zmax & tmp$ieps <= opt$epsmax, ]
#~ print(tmp)

#~ print(tmp$mH)
#~ print(tmp$w)
#~ print(tmp$idg)
#~ print(tmp$ieps)
#~ print(tmp$eps)
#~ print(suggestedrfact)
setfinalrresult <- data.frame(m=tmp$mH, w=tmp$w, iz=tmp$idg, ieps=tmp$ieps, eps=tmp$eps, statrfact=suggestedrfact)
if(!file.exists(filenamerfacts)){
  write.table(x=setfinalrresult, file=filenamerfacts, col.names=T, row.names=F, sep=",")
  write(x = paste("##", opt$ensemble), file = filenamerfacts, append = T)
} 

if(numrfactsuggested!=0) {
  print(paste("took", numrfacttaken/numrfactsuggested, "of the suggested points"))
  cat(paste("took", numrfacttaken/numrfactsuggested, "of the suggested points", opt$mode, opt$kernel, bmassaddon, opt$normnumber, opt$momentum), file=sprintf("%s/takensuggestions.txt", opt$plotfolder), append=T)
}
