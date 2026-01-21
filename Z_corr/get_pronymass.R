library("hadron", lib.loc="/hadron/gross/hadron_prony/hadron_lib")
library(optparse)

if (TRUE) {
    # set option list
option_list <- list(
    make_option(c("--file"), type = "character", default = "",
    help = "file from which the mass is read [default %default]"),
    make_option(c("--name"), type = "character", default = "pronymass",
    help = "name of pdf [default %default]"),
    make_option(c("--boot.R"), type = "integer", default = 1000,
    help = "bootstrap samples [default %default]"),
    make_option(c("--doubleboot.R"), type = "integer", default = 50,
    help = "double bootstrap samples [default %default]"),
    make_option(c("--boot.l"), type = "integer", default = 5,
    help = "bootstrap samples [default %default]"),
    make_option(c("--t0"), type = "integer", default = 1,
    help = "t0 of thc method [default %default]"),
    make_option(c("--onlymass"), action = "store_true", default = FALSE,
    help = "store effective mass correlators without prony data [default %default]"),
    make_option(c("--spectruminfos"), action = "store_true", default = FALSE,
    help = "store additional file with only t0 and evs value without any bootstrap values for plotting spectrum [default %default]")

)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
}



readcffrompreprocessed <- function (realfile, path  = "", imfile="", im=F, symmetrise=T) {
  ## determine number of measurements
  nmeas <- read.table(paste0(path, realfile), nrows=1)
  Time <- nmeas$V2
  confno <- nmeas$V1
  dat <- read.table(paste0(path, realfile), skip=1)
  tmp <- matrix(data=dat$V2, nrow=Time, ncol=confno)
  Dm <- dim(tmp)
  ret <- cf_meta(nrObs = 1, Time =  Time, nrStypes = 1)
  if (im==F) ret <- cf_orig(ret, cf = t(tmp))
  if (im==T) {
    stopifnot(file.exists(paste0(path, imfile)))
    datim <- read.table(paste0(path, imfile), skip=1)
    tmpim <- matrix(data=datim$V2, nrow=Time, ncol=confno)
    ret <- cf_orig(ret, cf = t(tmp), icf=t(tmpim))
  }
  ret$conf.index <- confno
  if(symmetrise) ret <- symmetrise.cf(ret)
  
  return(invisible(ret))
}

mass <- readcffrompreprocessed(realfile=opt$file)
mass <- unsymmetrise.cf(mass)
mass <- bootstrap.cf(mass, boot.l=opt$boot.l, boot.R=opt$boot.R)
mass <- double_bootstrap.cf(mass, dbboot.R=opt$doubleboot.R)

pronyres <- bootstrap.truncated.pgevm(mass, t0=opt$t0)
pronymass <- pgevm2effectivemass(pronyres, errortype = "dbboot")

saveRDS(pronyres, file=sprintf("prony_%s.RDS", opt$name))
if(!opt$onlymass) saveRDS(pronymass, file=sprintf("pronymass_%s.RDS", opt$name))
if(opt$onlymass) {
    pronymass$pgevm <- NA
    saveRDS(pronymass, file=sprintf("pronyonlymass_%s.RDS", opt$name))
}

if(opt$spectruminfos) {
    tmp <- pronyres
    tmp$evs.tsboot <- NA
    tmp$evs.dbboot <- NA
    tmp$cf$cf.tsboot$t <- NA
    tmp$cf$cf.tsboot$doubleboot$cf <- NA
    saveRDS(pronymass, file=sprintf("pronyspectruminfos_%s.RDS", opt$name))
}

pdf(sprintf("THC_spectrum_%s.pdf", opt$name), title="")
plot(pronyres)
plot(pronymass)
plot(pronymass, ylim=mean(pronymass$effMass) + c(-1, 1)*sd(pronymass$effMass))
dev.off()
