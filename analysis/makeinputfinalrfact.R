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
    help = "folder with results [default %default]")
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

rselection <- read.table(sprintf("%s/setrfacts_%s_%s.csv", opt$folder, opt$mode, opt$kernel), header=T, sep=",")

makeinputfile <- sprintf("%s/chosenrfactinput_%s_%s.txt", opt$folder, opt$mode, opt$kernel)
commandfile <- sprintf("%s/chosenrfactcommands_%s_%s.txt", opt$folder, opt$mode, opt$kernel)

cat("#! /bin/bash", sprintf("mkdir -p %s/inputchosenrfact", opt$folder) , append=F, sep="\n", file=makeinputfile)
cat("#! /bin/bash" , append=F, sep="\n", file=commandfile)

for(i in seq_along(rselection$m)) {
    statfolder <- sprintf("%s/%s_iz_%d_ieps_%d_stat", opt$folder, opt$kernel, rselection$iz[i], rselection$ieps[i])
    sysfolder <- sprintf("%s/%s_iz_%d_ieps_%d_sys", opt$folder, opt$kernel, rselection$iz[i], rselection$ieps[i])
    filename <- sprintf("%s/inputchosenrfact/%s_%s_iz_%d_ieps_%d", opt$folder, modefolder, opt$kernel, rselection$iz[i], rselection$ieps[i])
    jobname <- sprintf("%s_%s_%s_%s_%s_iz_%d_ieps_%d", opt$channel, opt$ensemble, opt$momentum, opt$mode, opt$kernel, rselection$iz[i], rselection$ieps[i])
    cat(sprintf("mkdir -p %s", statfolder), 
    sprintf("mkdir -p %s", sysfolder),
    sprintf("sed -e 's/finalrfact/%f/' -e 's/sigmafinal/%f/' %s/%s_%s_finalrfact.in > %s_stat.in", rselection$statrfact[i], rselection$eps[i], opt$folder, modefolder, opt$kernel, filename),
    sprintf("sed -e 's/finalrfact/%f/' -e 's/sigmafinal/%f/' %s/%s_%s_finalrfact.in > %s_sys.in", rselection$statrfact[i]/10, rselection$eps[i], opt$folder, modefolder, opt$kernel, filename),
    append=T, sep="\n", file=makeinputfile)
    cat(sprintf("sbatch --job-name=%s scriptdgdq2.sh %s_stat.in %s", jobname, filename, statfolder),
    sprintf("sbatch --job-name=%s scriptdgdq2.sh %s_sys.in %s", jobname, filename, sysfolder),
    append=T, sep="\n", file=commandfile)
    
}

print(sprintf("source %s", makeinputfile))
print(sprintf("source %s", commandfile))
