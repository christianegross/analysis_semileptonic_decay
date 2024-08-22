#~ parent <- "/hiskp4/gross/heavymesons/data"
parent <- "~/Documents/heavymesons/data/newinput"
folders <- c("cB211.07.64", "cB211.07.96", "cC211.06.80_600", "cD211.054.96", "cE211.044.112", "cB211.07.48_300", "cB211.07.48_400", "cB211.07.64_48_36")
nameshort <- c("B64", "B96", "C80", "D96", "E112", "B48_300", "B48_400", "B64_48_36")
subfolders  <- c("th1", "th2", "th3", "th4", "th5", "th6", "th7", "th8", "th9", "th9.5")

infos <- read.table("parameters_input_files.csv", row.names=1, header=TRUE, sep=",", comment.char = "#")
#~ print(infos)
names(infos)

system(sprintf("mkdir -p %s/cd/", parent))
system(sprintf("mkdir -p %s/cs/", parent))

index <- 2

if(FALSE) {
for(channel in c("cd", "cs")) {
    for(kernel in c("sigmoid", "erf")){
        for (i in seq_along(folders)) {
            folder <- folders[i]
                ensinfos <- infos[folder]
    #~             print(ensinfos)
            fileDG <- sprintf("%s/%.2d-%s-%s-%s-DGDq2.Rmd", parent, index, channel, folder, kernel)
            fileDM <- sprintf("%s/%.2d-%s-%s-%s-DMDq2.Rmd", parent, index, channel, folder, kernel)
    #~         system(sprintf("touch %s %s", fileDG, fileDM))
            celltitleDG <- sprintf("%s%s%s", nameshort[i], channel, kernel)
            print(ensinfos["T", ])
            cat(    
            sprintf("# %s - %s - %s", folder, channel, kernel), 
            "", 
            sprintf("```{r setup%s, include=FALSE}", celltitleDG), 
            "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE, results = FALSE)", 
            "source(\"/home/gross/Documents/heavymesons/scripts_christiane/extrapolate_sigma_templates/functions_stability_plots.R\")", 
            "source(\"/home/gross/Documents/heavymesons/scripts_christiane/extrapolate_sigma_templates/calc_DGDq2.R\")", 
            "library(\"hadron\")", 
            "```", 
            "", 
            "", 
            "", 
            "## Infos", 
            "", 
            sprintf("```{r %s effmass}", celltitleDG),
            sprintf("mass <- read.table(\"/home/gross/Documents/heavymesons/data/%s/%s/th1/outputY/Y.log\", fill=TRUE, skip=%d+35, ", channel, folder, ensinfos["T", ]), 
            "                   col.names=c(\"iset\", \"t\", \"H\", \"DH\",", 
            "                        \"emH\", \"DemH\", \"mH\", \"DmH\",", 
            "                        paste0(\"col\", 9:20)))", 
            "", 
            "cutlow <- 10", 
            "cuthigh <- 15", 
            "", 
            "for(iset in seq(min(mass$iset), max(mass$iset))){", 
            "    timeextent <- length(mass$t[mass$iset==iset])", 
            "", 
            "plot(NA, xlim = c(min(mass$t[mass$iset==iset]), max(mass$t[mass$iset==iset])),", 
            "         ylim = c(min(mass$emH[mass$iset==iset][2:(timeextent-2)] - mass$DemH[mass$iset==iset][2:(timeextent-2)]),", 
            "            max(mass$emH[mass$iset==iset][2:(timeextent-2)] + mass$DemH[mass$iset==iset][2:(timeextent-2)])),", 
            "         xlab=\"t/a\", ylab=\"m_eff(t)\", main=paste(\"set\", iset))", 
            "", 
            "polygon(x=c(mass$t[mass$iset==iset & mass$mH > 0], rev(mass$t[mass$iset==iset & mass$mH > 0])),", 
            "    y=c(mass$mH[mass$iset==iset & mass$mH > 0] - mass$DmH[mass$iset==iset & mass$mH > 0],", 
            "        mass$mH[mass$iset==iset & mass$mH > 0] + mass$DmH[mass$iset==iset & mass$mH > 0]),", 
            "        col=\"red\")", 
            "", 
            "meffstring <- tex.catwitherror(x=max(mass$mH[mass$iset==iset & mass$mH > 0]), dx=max(mass$DmH[mass$iset==iset & mass$mH > 0]), ", 
            "                               digits = 2, with.dollar = FALSE)", 
            "", 
            "plotwitherror(x=mass$t[mass$iset==iset], y=mass$emH[mass$iset==iset],", 
            "    dy=mass$DemH[mass$iset==iset], rep=TRUE)", 
            "", 
            "", 
            "plot(NA, xlim = c(min(mass$t[mass$iset==iset]), max(mass$t[mass$iset==iset])),", 
            "         ylim = c(min(mass$emH[mass$iset==iset][cutlow:(timeextent-cuthigh)] - mass$DemH[mass$iset==iset][cutlow:(timeextent-cuthigh)]),", 
            "            max(mass$emH[mass$iset==iset][cutlow:(timeextent-cuthigh)] + mass$DemH[mass$iset==iset][cutlow:(timeextent-cuthigh)])),", 
            "         xlab=\"t/a\", ylab=\"m_eff(t)\", main=paste(\"set\", iset))", 
            "", 
            "polygon(x=c(mass$t[mass$iset==iset & mass$mH > 0], rev(mass$t[mass$iset==iset & mass$mH > 0])),", 
            "    y=c(mass$mH[mass$iset==iset & mass$mH > 0] - mass$DmH[mass$iset==iset & mass$mH > 0],", 
            "        mass$mH[mass$iset==iset & mass$mH > 0] + mass$DmH[mass$iset==iset & mass$mH > 0]),", 
            "        col=\"red\")", 
            "", 
            "plotwitherror(x=mass$t[mass$iset==iset], y=mass$emH[mass$iset==iset],", 
            "    dy=mass$DemH[mass$iset==iset], rep=TRUE)", 
            "}", 
            "```", 
            "", 
            "", 
            sprintf("```{r %sinfo, results=TRUE}", celltitleDG), 
            "infos <- data.frame(name=c(\"Volume\", \"afm\", \"mpigev\", \"ZA\", \"ZV\", \"nconf\", \"tsink\", \"tj\", ", 
            "                        \"tmeff1\", \"tmeff2\", \"meff\", \"nnorms\", \"nsigma\"),", 
            sprintf("                    value=c(\"%d^3x%d\", \"%f\", \"%f\", \"%f\", \"%f\", ", ensinfos["L", ], ensinfos["T", ], ensinfos["afm", ], ensinfos["mpi", ], ensinfos["za", ], ensinfos["zv", ]), 
            sprintf("                        \"300\", \"%d\", \"%d\", \"%d\", \"%d\", meffstring, \"3\", \"%d\"))", ensinfos["tsink", ], ensinfos["tj", ], ensinfos["tfitlower", ], ensinfos["tfitupper", ], ensinfos[paste0("nsigma_", kernel), ]), 
            "knitr::kable(infos, col.names=c(\"\", \"\"), caption = \"Information about the ensemble\")", 
            sprintf("sigmas <- c(\"%s\")", paste(ensinfos[paste0("sigma_", kernel, 1:ensinfos[paste0("nsigma_", kernel), ]), ], collapse="\", \"")), 
            "thetas <- data.frame(names=c(\"th1\", \"th2\", \"th3\", \"th4\", \"th5\", \"th6\", \"th7\", \"th8\", \"th9\", \"th9.5\"),", 
            sprintf("                     values=c(\"%s\"))", paste(ensinfos[subfolders, ], collapse="\", \"")), 
            "knitr::kable(thetas, col.names=c(\"\", \"\"), caption = \"Used thetas\")", 
            "```", 
            "", 
            "The sigma we used are `r pander::p(sigmas, wrap='')`.", 
            file=fileDG, sep="\n", append=FALSE)
            system(sprintf("cp %s %s", fileDG, fileDM))
            for (subfolder in subfolders) {
                cat(
                sprintf("## %s", subfolder),
                "",
                sprintf("```{r %s%spath}", celltitleDG, subfolder),
                sprintf("%spath <- \"/home/gross/Documents/heavymesons/data/%s/%s/%s/%s_new/outputDGammaDq2/\"", subfolder, channel, folder, subfolder, kernel),
                "```",
                "",
                "",
                "### Z",
                "",
                sprintf("```{r %s%sZ}", celltitleDG, subfolder),
                sprintf("plotZ(path=%spath)", subfolder),
                "```",
                "",
                "### stability plots A/A0",
                "",
                sprintf("```{r %s%sstabilaa0, eval=!(knitr::is_html_output())}", celltitleDG, subfolder),
                sprintf("plot_stability_AA0(path=%spath, inputfile=paste0(%spath, \"/../../DGammaDq2_%s_new.in\"), nnorm=3, neps=%d)", subfolder, subfolder, kernel, ensinfos[paste0("nsigma_", kernel), ]),
                "```",
                "",
                "",
                "### stability plots lambda",
                "",
                sprintf("```{r %s%sstabillambda, eval=!(knitr::is_html_output())}", celltitleDG, subfolder),
                sprintf("plot_stability_lambda(path=%spath, inputfile=paste0(%spath, \"/../../DGammaDq2_%s_new.in\"), nnorm=3, neps=%d)", subfolder, subfolder, kernel, ensinfos[paste0("nsigma_", kernel), ]),
                "```",
                "",
                "### Kernel reconstruction",
                "",
                sprintf("```{r %s%skernel, eval=!(knitr::is_html_output())}", celltitleDG, subfolder),
                sprintf("plot_reconstruction(path=%spath, inputfile=paste0(%spath, \"/../../DGammaDq2_%s_new.in\"), nnorm=3, neps=%d)", subfolder, subfolder, kernel, ensinfos[paste0("nsigma_", kernel), ]),
                "```\n",
                file=fileDG, sep="\n", append=TRUE)
                
                cat(
                sprintf("## %s", subfolder),
                "",
                sprintf("```{r %s%spath}", celltitleDG, subfolder),
                sprintf("%spath <- \"/home/gross/Documents/heavymesons/data/%s/%s/%s/%s_new/outputDMDq2/\"", subfolder, channel, folder, subfolder, kernel),
                "```",
                "",
                "",
                "### Z",
                "",
                sprintf("```{r %s%sZ}", celltitleDG, subfolder),
                sprintf("plotZ(path=%spath, mode=\"DM\")", subfolder),
                "```",
                "",
                "### stability plots A/A0",
                "",
                sprintf("```{r %s%sstabilaa0, eval=!(knitr::is_html_output())}", celltitleDG, subfolder),
                sprintf("plot_stability_AA0(path=%spath, inputfile=paste0(%spath, \"/../../DMDq2_%s_new.in\"), nnorm=3, neps=%d, mode=\"DM\")", subfolder, subfolder, kernel, ensinfos[paste0("nsigma_", kernel), ]),
                "```",
                "",
                "",
                "### stability plots lambda",
                "",
                sprintf("```{r %s%sstabillambda, eval=!(knitr::is_html_output())}", celltitleDG, subfolder),
                sprintf("plot_stability_lambda(path=%spath, inputfile=paste0(%spath, \"/../../DMDq2_%s_new.in\"), nnorm=3, neps=%d, mode=\"DM\")", subfolder, subfolder, kernel, ensinfos[paste0("nsigma_", kernel), ]),
                "```",
                "",
                "### Kernel reconstruction",
                "",
                sprintf("```{r %s%skernel, eval=!(knitr::is_html_output())}", celltitleDG, subfolder),
                sprintf("plot_reconstruction(path=%spath, inputfile=paste0(%spath, \"/../../DMDq2_%s_new.in\"), nnorm=3, neps=%d, mode=\"DM\")", subfolder, subfolder, kernel, ensinfos[paste0("nsigma_", kernel), ]),
                "```\n",
                
                file=fileDM, sep="\n", append=TRUE)
            }
            index <- index+1
        }
    }
}
}

if(TRUE) {
index <- 1
#~ subfolders  <- c("th1", "th2", "th3", "th4", "th5", "th6")
subfolders  <- c("th2_su", "th4_su", "th6_su", "th8_su", "th9.5_su", "th1", "th2", "th3", "th4", "th5", "th6")
folders <- c("cB211.07.64")
for(channel in c("su")) {
    for(kernel in c("sigmoid", "erf")){
        for (i in seq_along(folders)) {
            folder <- folders[i]
                ensinfos <- infos[folder]
    #~             print(ensinfos)
            fileDG <- sprintf("%s/%.2d-%s-%s-%s-DGDq2.Rmd", parent, 100+index, channel, folder, kernel)
            fileDM <- sprintf("%s/%.2d-%s-%s-%s-DMDq2.Rmd", parent, 100+index, channel, folder, kernel)
    #~         system(sprintf("touch %s %s", fileDG, fileDM))
            celltitleDG <- sprintf("%s%s%s", nameshort[i], channel, kernel)
            print(ensinfos["T", ])
            cat(    
            sprintf("# %s - %s - %s", folder, channel, kernel), 
            "", 
            sprintf("```{r setup%s, include=FALSE}", celltitleDG), 
            "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE, results = FALSE)", 
            "source(\"/home/gross/Documents/heavymesons/scripts_christiane/extrapolate_sigma_templates/functions_stability_plots.R\")", 
            "source(\"/home/gross/Documents/heavymesons/scripts_christiane/extrapolate_sigma_templates/calc_DGDq2.R\")", 
            "library(\"hadron\")", 
            "```", 
            "", 
            "", 
            "", 
            "## Infos", 
            "", 
            sprintf("```{r %s effmass, eval=T}", celltitleDG),
            sprintf("mass <- read.table(\"/home/gross/Documents/heavymesons/data/%s/%s/th1/outputY/Y.log\", fill=TRUE, skip=%d+35, ", channel, folder, ensinfos["T", ]), 
            "                   col.names=c(\"iset\", \"t\", \"H\", \"DH\",", 
            "                        \"emH\", \"DemH\", \"mH\", \"DmH\",", 
            "                        paste0(\"col\", 9:20)))", 
            "", 
            "cutlow <- 10", 
            "cuthigh <- 15", 
            "", 
            "for(iset in seq(min(mass$iset), max(mass$iset))){", 
            "    timeextent <- length(mass$t[mass$iset==iset])", 
            "", 
            "plot(NA, xlim = c(min(mass$t[mass$iset==iset]), max(mass$t[mass$iset==iset])),", 
            "         ylim = c(min(mass$emH[mass$iset==iset][2:(timeextent-2)] - mass$DemH[mass$iset==iset][2:(timeextent-2)]),", 
            "            max(mass$emH[mass$iset==iset][2:(timeextent-2)] + mass$DemH[mass$iset==iset][2:(timeextent-2)])),", 
            "         xlab=\"t/a\", ylab=\"m_eff(t)\", main=paste(\"set\", iset))", 
            "", 
            "polygon(x=c(mass$t[mass$iset==iset & mass$mH > 0], rev(mass$t[mass$iset==iset & mass$mH > 0])),", 
            "    y=c(mass$mH[mass$iset==iset & mass$mH > 0] - mass$DmH[mass$iset==iset & mass$mH > 0],", 
            "        mass$mH[mass$iset==iset & mass$mH > 0] + mass$DmH[mass$iset==iset & mass$mH > 0]),", 
            "        col=\"red\")", 
            "", 
            "meffstring <- tex.catwitherror(x=max(mass$mH[mass$iset==iset & mass$mH > 0]), dx=max(mass$DmH[mass$iset==iset & mass$mH > 0]), ", 
            "                               digits = 2, with.dollar = FALSE)", 
            "", 
            "plotwitherror(x=mass$t[mass$iset==iset], y=mass$emH[mass$iset==iset],", 
            "    dy=mass$DemH[mass$iset==iset], rep=TRUE)", 
            "", 
            "", 
            "plot(NA, xlim = c(min(mass$t[mass$iset==iset]), max(mass$t[mass$iset==iset])),", 
            "         ylim = c(min(mass$emH[mass$iset==iset][cutlow:(timeextent-cuthigh)] - mass$DemH[mass$iset==iset][cutlow:(timeextent-cuthigh)]),", 
            "            max(mass$emH[mass$iset==iset][cutlow:(timeextent-cuthigh)] + mass$DemH[mass$iset==iset][cutlow:(timeextent-cuthigh)])),", 
            "         xlab=\"t/a\", ylab=\"m_eff(t)\", main=paste(\"set\", iset))", 
            "", 
            "polygon(x=c(mass$t[mass$iset==iset & mass$mH > 0], rev(mass$t[mass$iset==iset & mass$mH > 0])),", 
            "    y=c(mass$mH[mass$iset==iset & mass$mH > 0] - mass$DmH[mass$iset==iset & mass$mH > 0],", 
            "        mass$mH[mass$iset==iset & mass$mH > 0] + mass$DmH[mass$iset==iset & mass$mH > 0]),", 
            "        col=\"red\")", 
            "", 
            "plotwitherror(x=mass$t[mass$iset==iset], y=mass$emH[mass$iset==iset],", 
            "    dy=mass$DemH[mass$iset==iset], rep=TRUE)", 
            "}", 
            "```", 
            "", 
            "", 
            sprintf("```{r %sinfo, results=TRUE, eval=T}", celltitleDG), 
            "infos <- data.frame(name=c(\"Volume\", \"afm\", \"mpigev\", \"ZA\", \"ZV\", \"nconf\", \"tsink\", \"tj\", ", 
            "                        \"tmeff1\", \"tmeff2\", \"meff\", \"nnorms\", \"nsigma\"),", 
            sprintf("                    value=c(\"%d^3x%d\", \"%f\", \"%f\", \"%f\", \"%f\", ", ensinfos["L", ], ensinfos["T", ], ensinfos["afm", ], ensinfos["mpi", ], ensinfos["za", ], ensinfos["zv", ]), 
            sprintf("                        \"300\", \"%d\", \"%d\", \"%d\", \"%d\", meffstring, \"3\", \"%d\"))", ensinfos["tsink", ], ensinfos["tj", ], ensinfos["tfitlower", ], ensinfos["tfitupper", ], ensinfos[paste0("nsigma_", kernel), ]), 
            "knitr::kable(infos, col.names=c(\"\", \"\"), caption = \"Information about the ensemble\")", 
            sprintf("sigmas <- c(\"%s\")", paste(ensinfos[paste0("sigma_", kernel, 1:ensinfos[paste0("nsigma_", kernel), ]), ], collapse="\", \"")), 
            "thetas <- data.frame(names=c(\"th2_su\", \"th4_su\", \"th6_su\", \"th8_su\", \"th9.5_su\", \"th1\", \"th2\", \"th3\", \"th4\", \"th5\", \"th6\"),", 
            sprintf("                     values=c(\"%s\"))", paste(ensinfos[subfolders, ], collapse="\", \"")), 
            "knitr::kable(thetas, col.names=c(\"\", \"\"), caption = \"Used thetas\")", 
            "```", 
            "", 
            "The sigma we used are `r pander::p(sigmas, wrap='')`.", 
            file=fileDG, sep="\n", append=FALSE)
            system(sprintf("cp %s %s", fileDG, fileDM))
            for (subfolder in subfolders) {
                cat(
                sprintf("\n## %s", subfolder),
                "",
                sprintf("```{r %s%spath}", celltitleDG, subfolder),
                sprintf("%spath <- \"/home/gross/Documents/heavymesons/data/%s/%s/%s/%s_new/outputDGammaDq2/\"", subfolder, channel, folder, subfolder, kernel),
                "```",
                "",
                "",
                "### Z",
                "",
                sprintf("```{r %s%sZ}", celltitleDG, subfolder),
                sprintf("plotZ(path=%spath)", subfolder),
                "```",
                "",
                "### stability plots A/A0",
                "",
                sprintf("```{r %s%sstabilaa0, eval=!(knitr::is_html_output())}", celltitleDG, subfolder),
                sprintf("plot_stability_AA0(path=%spath, inputfile=paste0(%spath, \"/../../DGammaDq2_%s_new.in\"), nnorm=3, neps=%d, comblist=c(0))", subfolder, subfolder, kernel, ensinfos[paste0("nsigma_", kernel), ]),
                "```",
                "",
                "",
                "### stability plots lambda",
                "",
                sprintf("```{r %s%sstabillambda, eval=(knitr::is_html_output())}", celltitleDG, subfolder),
                sprintf("plot_stability_lambda(path=%spath, inputfile=paste0(%spath, \"/../../DGammaDq2_%s_new.in\"), nnorm=3, neps=%d, comblist=c(0))", subfolder, subfolder, kernel, ensinfos[paste0("nsigma_", kernel), ]),
                "```",
                "",
                "### Kernel reconstruction",
                "",
                sprintf("```{r %s%skernel, eval=!(knitr::is_html_output())}", celltitleDG, subfolder),
                sprintf("plot_reconstruction(path=%spath, inputfile=paste0(%spath, \"/../../DGammaDq2_%s_new.in\"), nnorm=3, neps=%d, comblist=c(0))", subfolder, subfolder, kernel, ensinfos[paste0("nsigma_", kernel), ]),
                "```\n",
                file=fileDG, sep="\n", append=TRUE)
                
                cat(
                sprintf("## %s", subfolder),
                "",
                sprintf("```{r %s%spath}", celltitleDG, subfolder),
                sprintf("%spath <- \"/home/gross/Documents/heavymesons/data/%s/%s/%s/%s_new/outputDMDq2/\"", subfolder, channel, folder, subfolder, kernel),
                "```",
                "",
                "",
                "### Z",
                "",
                sprintf("```{r %s%sZ}", celltitleDG, subfolder),
                sprintf("plotZ(path=%spath, mode=\"DM\")", subfolder),
                "```",
                "",
                "### stability plots A/A0",
                "",
                sprintf("```{r %s%sstabilaa0, eval=!(knitr::is_html_output())}", celltitleDG, subfolder),
                sprintf("plot_stability_AA0(path=%spath, inputfile=paste0(%spath, \"/../../DMDq2_%s_new.in\"), nnorm=3, neps=%d, mode=\"DM\")", subfolder, subfolder, kernel, ensinfos[paste0("nsigma_", kernel), ]),
                "```",
                "",
                "",
                "### stability plots lambda",
                "",
                sprintf("```{r %s%sstabillambda, eval=!(knitr::is_html_output())}", celltitleDG, subfolder),
                sprintf("plot_stability_lambda(path=%spath, inputfile=paste0(%spath, \"/../../DMDq2_%s_new.in\"), nnorm=3, neps=%d, mode=\"DM\")", subfolder, subfolder, kernel, ensinfos[paste0("nsigma_", kernel), ]),
                "```",
                "",
                "### Kernel reconstruction",
                "",
                sprintf("```{r %s%skernel, eval=!(knitr::is_html_output())}", celltitleDG, subfolder),
                sprintf("plot_reconstruction(path=%spath, inputfile=paste0(%spath, \"/../../DMDq2_%s_new.in\"), nnorm=3, neps=%d, mode=\"DM\")", subfolder, subfolder, kernel, ensinfos[paste0("nsigma_", kernel), ]),
                "```\n",
                
                file=fileDM, sep="\n", append=TRUE)
            }
            index <- index+1
        }
    }
}
}
