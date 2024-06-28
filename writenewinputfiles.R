parent <- "/hiskp4/gross/heavymesons/data"
#~ parent <- "~/Documents/heavymesons/data/newinput"
folders <- c("cB211.07.64", "cB211.07.96", "cC211.06.80", "cD211.054.96", "cE211.044.112", "cB211.07.48_300", "cB211.07.48_400", "cB211.07.64_48_36")
folders <- c("cB211.07.48_charm", "cB211.07.48_charm_strange")
subfolders  <- c("th1", "th2", "th3", "th4", "th5", "th6", "th7", "th8", "th9", "th9.5")
subfolders  <- c("th2", "th4", "th6", "th8", "th9.5")

infos <- read.table("parameters_input_files.csv", row.names=1, header=TRUE, sep=",", comment.char = "#")
#~ print(infos)
names(infos)

system(sprintf("mkdir -p %s/cd/", parent))
system(sprintf("mkdir -p %s/cs/", parent))

for(channel in c("cd", "cs")) {
    for (folder in folders) {
        system(sprintf("mkdir -p %s/%s/%s", parent, channel, folder))
            ensinfos <- infos[folder]
#~             print(ensinfos)
        for (subfolder in subfolders) {
            path <- sprintf("%s/%s/%s/%s", parent, channel, folder, subfolder)
            system(sprintf("mkdir -p %s/", path))
            system(sprintf("mkdir -p %s/erf_new", path))
            system(sprintf("mkdir -p %s/sigmoid_new", path))
#~             print(sprintf("rm -rfv %s/outputY", path))
            Yin <- sprintf("%s/Y.in", path)
            sigmoidin <- sprintf("%s/DGammaDq2_sigmoid_new.in", path)
            momentsigmoidin <- sprintf("%s/DMDq2_sigmoid_new.in", path)
            erfin <- sprintf("%s/DGammaDq2_erf_new.in", path)
            momenterfin <- sprintf("%s/DMDq2_erf_new.in", path)
            system(sprintf("rm -f %s %s %s %s %s", Yin, sigmoidin, momentsigmoidin, erfin, momenterfin))
            cat("[Run Parameters]", 
                sprintf("L\t\t\t\t%d", ensinfos["L", ]), 
                sprintf("T\t\t\t\t%d", ensinfos["T", ]), 
                sprintf("afm\t\t\t\t%f", ensinfos["afm", ]), 
                sprintf("mpigev\t\t\t%f", ensinfos["mpi", ]), 
                sprintf("za\t\t\t\t%f", ensinfos["zv", ]), 
                sprintf("zv\t\t\t\t%f", ensinfos["za", ]), 
                file=Yin, sep="\n", append=FALSE)

            cat("\n[Analysis Parameters]",
                "\ndatatype\t\tRAW\t\t\t# RAW,JACK,BOOT",
                "statistics\t\tBOOT\t\t\t# JACK,BOOT",
                "nboot\t\t\t1000", 
                sprintf("bs\t\t\t\t%d", ensinfos["binlength", ]),
                "nsets\t\t\t1", file=Yin, sep="\n", append=TRUE)
                
            cat("\n[Set 1]",
                sprintf("theta			0.0 0.0 %f", ensinfos[subfolder, ]), 
                "n1\t\t\t\t1.0 0.0 0.0", 
                "n2\t\t\t\t0.0 1.0 0.0",
                sprintf("path\t\t\t%s/4Pts_%s_%s_data/", path, channel, subfolder),
                "Hfile\t\t\tDs.dat # we do not need the light meson, so it was not computed",
                "Lfile\t\t\tDs.dat",
                "L0file\t\t\tDs.dat",
                sprintf("tsink\t\t\t%d", ensinfos["tsink", ]),
                sprintf("tj\t\t\t\t%d", ensinfos["tj", ]),
                sprintf("tfith\t\t\t%d %d", ensinfos["tfitlower", ], ensinfos["tfitupper", ]),
                sprintf("tfitl\t\t\t%d %d", ensinfos["tfitlower", ], ensinfos["tfitupper", ]),
                sprintf("tfitl0\t\t\t%d %d", ensinfos["tfitlower", ], ensinfos["tfitupper", ]), 
                file=Yin, sep="\n", append=TRUE)
                
                
            cat("[MPFR]", 
                "precision\t\t400", 
                "\n[COMBOS]", 
                "type ZERO", 
                "\n[Y file]", 
                sprintf("filename\t\t%s/outputY/Y.bin", path), 
                "\n[Integration Space]", 
                "endpoint_type\t\tPINF\t\t# NUM,PINF", 
                sprintf("e0\t\t\t\t%.3f", ensinfos[paste0("e0_", channel), ]), 
                "am_H*e1\t\t\t2.0", 
                file=sigmoidin, sep="\n", append=FALSE)
                
            cat("\n[Norms]", 
                "nnorms\t\t\t\t3", 
                "norm1_power\t\t\t0.0", 
                "norm1_am_H*alpha\t0.0", 
                "norm2_power\t\t\t0.0", 
                "norm2_am_H*alpha\t1.99",
                "norm3_power\t\t\t0.0",  
                "norm3_am_H*alpha\t-1.99", 
                file=sigmoidin, sep="\n", append=TRUE)
                
                
            cat("\n[HLT]", 
                "tmin\t\t\t1", 
                sprintf("Nt\t\t\t\t%d", ensinfos["Nt", ]), 
                "covariance\t\tDIAG\t\t\t# FULL, DIAG", 
                "solve\t\t\tPLATEAUX\t\t# MULTI,PLATEAUX,AoB", 
                "nsolve\t\t\t6", 
                "nerr\t\t\t2", 
                "start\t\t\t0.3", 
                file=sigmoidin, sep="\n", append=TRUE)
            
            system(sprintf("cp %s %s", sigmoidin, momentsigmoidin))
            
            cat("normalization\t\t1\t\t# normalization=1 computes only dM, normalization=0 computes dM/dG", 
                file=momentsigmoidin, sep="\n", append=TRUE)
                
            cat("#solve\t\t\tAoB\t\t\t# MULTI,PLATEAUX,AoB\n#rfact\t\t\t1.0", 
                file=sigmoidin, sep="\n", append=TRUE)
            cat("#solve\t\t\tAoB\t\t\t# MULTI,PLATEAUX,AoB\n#rfact\t\t\t1.0", 
                file=momentsigmoidin, sep="\n", append=TRUE)

                
            
            system(sprintf("cp %s %s", sigmoidin, erfin))
            system(sprintf("cp %s %s", momentsigmoidin, momenterfin))
            
            namessigmoid <- paste0("sigma_sigmoid", 1:ensinfos["nsigma_sigmoid", ])
            nameserf <- paste0("sigma_erf", 1:ensinfos["nsigma_erf", ])
            cat("\n[Kernel]\n",
                "kernel\t\tSIGMOID # SIGMOID, ERF",
                "\n[Kernel Parameters]", 
                sprintf("nsigma\t%d", ensinfos["nsigma_sigmoid", ]),
                paste0("sigma", 1:ensinfos["nsigma_sigmoid", ], " ", ensinfos[namessigmoid, ]),
                file=sigmoidin, sep="\n", append=TRUE)
            
            cat("\n[Kernel]\n",
                "kernel\t\tSIGMOID # SIGMOID, ERF", 
                "\n[Kernel Parameters]", 
                sprintf("nsigma\t%d", ensinfos["nsigma_sigmoid", ]),
                paste0("sigma", 1:ensinfos["nsigma_sigmoid", ], " ", ensinfos[namessigmoid, ]),
                file=momentsigmoidin, sep="\n", append=TRUE)
            
            cat("\n[Kernel]\n",
                "kernel\t\tERF # SIGMOID, ERF", 
                "\n[Kernel Parameters]", 
                sprintf("nsigma\t%d", ensinfos["nsigma_erf", ]),
                paste0("sigma", 1:ensinfos["nsigma_erf", ], " ", ensinfos[nameserf, ]),
                file=erfin, sep="\n", append=TRUE)
            
            cat("\n[Kernel]\n",
                "kernel\t\tERF # SIGMOID, ERF", 
                "\n[Kernel Parameters]", 
                sprintf("nsigma\t%d", ensinfos["nsigma_erf_dmdq2", ]),
                paste0("sigma", 1:ensinfos["nsigma_erf_dmdq2", ], " ", ensinfos[nameserf[1:ensinfos["nsigma_erf_dmdq2", ]], ]),
                file=momenterfin, sep="\n", append=TRUE)



        }
    }
}
