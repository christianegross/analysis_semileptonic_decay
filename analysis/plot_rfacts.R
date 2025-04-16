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
    make_option(c("-b", "--bmass"), type = "integer", default = "-1",
    help = "index for bmass [default %default]")
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
}

if(opt$folder=="-1") {
    opt$folder <- sprintf("~/Documents/heavymesons/data/%s/%s/%s", opt$channel, opt$ensemble, opt$momentum)
}


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

rfacts <- 10^(1:5)

bmassaddon <- ifelse(opt$bmass==-1, "", sprintf("_bmass_%d", opt$bmass))


pdf(sprintf("%s/rfact_%s_%s%s.pdf", opt$folder, opt$mode, opt$kernel, bmassaddon), title="")

filenamerfacts <- sprintf("%s/setrfacts_%s_%s%s.csv", opt$folder, opt$mode, opt$kernel, bmassaddon)

resultsselected <- file.exists(filenamerfacts)
if(resultsselected) rselection <- read.table(filenamerfacts, header=T, sep=",")

for(iz in 0:zmax) {
    for(iepsilon in 0:40) {
        res <- data.frame(ik=c(), spectreflag=c(), lambda_start=c(), lambdalambda_start=c(), Bnorm=c(), A0=c(), AA0_min=c(), AA0_ref=c(), AA0=c(), BnormB_ref=c(), BnormB=c(), C_ref=c(), C=c(), rho=c(), drho_stat=c(), drho_syst=c(), drho_tot=c(), resflag=c(), rfact=c())
        
        for(rfact in rfacts) {
            filename <- sprintf("%s/%s_%.2e/output%s/%s_%d_iset_0_ieps_%d_icomb_0.dat", opt$folder, opt$kernel, rfact, modefolder, opt$mode, iz, iepsilon)
                  datastab <- try(read.table(filename, fill=TRUE, skip=3,
                                             col.names=c("ik", "spectreflag", "lambda_start",
                                                         "lambdalambda_start", "Bnorm", "A0", "AA0_min",
                                                         "AA0_ref", "AA0", "BnormB_ref", "BnormB",
                                                         "C_ref", "C", "rho", "drho_stat", "drho_syst",
                                                         "drho_tot", "resflag")), silent=T)
                                                         
                  
                  if(!inherits(datastab, "try-error")) {
                    datastabaa0 <- datastab[datastab$spectreflag==1, ]
                    datastabaa0$rfact <- rfact
                    res <- rbind(res, datastabaa0)
                    
            }
        }
        if(length(res$AA0_ref) > 1) {
            plotwitherror(x=res$AA0_ref/res$BnormB*res$rfact, y=res$rho, dy=res$drho_stat, col=log10(res$rfact), xlab="A/B", ylab="rho", main=paste(opt$channel, opt$ensemble, opt$momentum, opt$mode, "eps", iepsilon, "Z", iz, bmassaddon), log="x")
            grid()
            abline(v=outer(1:9, 10^(1:7)), col = "lightgray", lty = "dotted")
            abline(v=(res$AA0_ref/res$BnormB*res$rfact)[res$resflag==1], col=log10(res$rfact[res$resflag==1]))
                        
            legend(x="topleft", legend=sprintf("%.2e", rfacts), col=log10(rfacts), pch=c(1))
            
            if(resultsselected){
                abline(v=rselection$statrfact[rselection$iz==iz & rselection$ieps==iepsilon], lwd=2)
                legend(x="topright", legend="selected", col="black", lty=1, lwd=2)
            } 
        }
    
    }
}


tmp <- read.table(sprintf("%s/%s_%.2e/output%s/%s.dat", opt$folder, opt$kernel, rfacts[1], modefolder, modefolder), col.names=c("mH", "iset", "w", "icomb", "idg", "ieps", "eps", "rho", "stat", "sys", "tot"))
tmp <- tmp[tmp$icomb==0 & tmp$idg <= zmax, ]

setfinalrresult <- data.frame(m=tmp$mH, w=tmp$w, iz=tmp$idg, ieps=tmp$ieps, eps=tmp$eps, statrfact=NA)
if(!file.exists(filenamerfacts)){
    write.table(x=setfinalrresult, file=filenamerfacts, col.names=T, row.names=F)
} 

