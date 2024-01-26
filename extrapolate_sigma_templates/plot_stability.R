library(hadron)
library(optparse)

if (TRUE) {
    # set option list
option_list <- list(
make_option(c("--scriptpath"), type="character", default="./",
            help = "path to where the script with the function definitions is stored"),
make_option(c("--path"), type="character", default="./",
            help = "path to where the data are stored"),
make_option(c("--plotpath"), type="character", default="./",
            help = "path to where the plots are stored"),
make_option(c("--neps"), type="integer", default=10,
            help = "number of epsilons that are considered"),
make_option(c("--nset"), type="integer", default=1,
            help = "number of sets that are considered"),
make_option(c("--nnorm"), type="integer", default=2,
            help = "number of norms that are considered"),
make_option(c("--comment"), type="character", default="",
            help = "comment to put on end of files"),
make_option(c("--calcZ"), action = "store_true", default = FALSE, 
            help = "if given, plot the Z-correlators"),
make_option(c("--inputfile"), type="character", default="",
            help = "inputfile to read in norms")
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options


source(paste0(opt$scriptpath, "/functions_stability_plots.R"))
}



pdf(sprintf("%s/stability_lambda_%s.pdf", opt$plotpath, opt$comment), title="")

plot_stability_lambda(path=opt$path, nset=opt$nset, neps=opt$neps, nnorm =opt$nnorm, inputfile=opt$inputfile)

pdf(sprintf("%s/stability_AA0_%s.pdf", opt$plotpath, opt$comment), title="")

plot_stability_AA0(path=opt$path, nset=opt$nset, neps=opt$neps, nnorm =opt$nnorm, inputfile=opt$inputfile)

pdf(sprintf("%s/reconstruct_kernel%s.pdf", opt$plotpath, opt$comment), title="")

plot_reconstruction(path=opt$path, nset=opt$nset, neps=opt$neps, nnorm=opt$nnorm, inputfile=opt$inputfile)


if (opt$calcZ) {
    pdf(sprintf("%s/Z.pdf", opt$plotpath), title="")
    plotZ(opt$path)
}
