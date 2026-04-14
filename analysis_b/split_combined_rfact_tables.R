library("hadron")
library("optparse")


if (TRUE) {
    # set option list
option_list <- list(
    make_option(c("-c", "--channel"), type = "character", default = "cs",
    help = "decay channel [default %default]"),
    make_option(c("-t", "--momentum"), type = "character", default = "th2",
    help = "momentum [default %default]"),
    make_option(c("-m", "--mode"), type = "character", default = "DG",
    help = "mode: DG, DM, DM2 [default %default]"),
    make_option(c("-k", "--kernel"), type = "character", default = "sigmoid",
    help = "kernel [default %default]"),
    make_option(c("-f", "--folder"), type = "character", default = "-1",
    help = "base folder where split results are written [default %default]"),
    make_option(c("-p", "--plotfolder"), type = "character", default = "-1",
    help = "folder where combined results are stored [default %default]"),
    make_option(c("-b", "--bmass"), type = "integer", default = "-1",
    help = "index for bmass [default %default]"),
    make_option(c("-n", "--normnumber"), type = "integer", default = "0",
    help = "index of norm to be analysed [default %default]"),
    make_option(c("--includeE112"), action = "store_true", default = FALSE,
                help = "if true, include E112 [default %default]")
    
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
}

if(opt$folder=="-1") {
    opt$folder <- sprintf("~/Documents/heavymesons/data/%s/%s", opt$channel, opt$momentum)
}
if(opt$plotfolder=="-1") {
    opt$plotfolder <- opt$folder
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


bmassaddon <- ifelse(opt$bmass==-1, "", sprintf("_bmass_%d", opt$bmass))



if(opt$includeE112) ensembles <- c("cB211.07.64", "cC211.06.80", "cD211.054.96", "cE211.044.112")
if(!opt$includeE112) ensembles <- c("cB211.07.64", "cC211.06.80", "cD211.054.96")

combined <- read.table(sprintf("%s/rfact_%s_%s%s_ik%d_%s_comb.csv", opt$plotfolder, opt$mode, opt$kernel, bmassaddon, opt$normnumber, opt$momentum), header=T)
if(length(combined$ieps)==0) combined <- read.table(sprintf("%s/rfact_%s_%s%s_ik%d_%s_comb.csv", opt$plotfolder, opt$mode, opt$kernel, bmassaddon, opt$normnumber, opt$momentum), header=T, sep=",")
if(length(combined$ieps)==0) stop("cannot read in resultfile")
#~ print(dim(combined))

for(ensemble in ensembles) {
    tmp <- combined[combined$ensemble==ensemble,]
#~     print(dim(tmp))
    tmp <- tmp[, -7]
#~     print(dim(tmp))
    tmp <- tmp[order(tmp$iz, tmp$ieps), ]
    write.table(tmp, sprintf("%s/%s/%s/plots_backcomparison_3_sqrtsum/rfact_%s_%s%s_ik%d_%s.csv", opt$folder, ensemble, opt$momentum, opt$mode, opt$kernel, bmassaddon, opt$normnumber, opt$momentum), row.names=F, quote=F)
    write(x = paste("##", ensemble), file = sprintf("%s/%s/%s/plots_backcomparison_3_sqrtsum/rfact_%s_%s%s_ik%d_%s.csv", opt$folder, ensemble, opt$momentum, opt$mode, opt$kernel, bmassaddon, opt$normnumber, opt$momentum), append = T)
}
