library("hadron")


# read in: develop new function, probably based on readtextcf
# takes filenames, number of lines to read, starting line
# how to get starting line: either hardcode 16 values or write c++ function to read commented out lines
# prefer hardcode


##reads in only first contractions, e.g. only A0, not possible to select A1
read_in_single_spin <- function(filelist, correlator="mes_contr_H_C_S_H", 
    Time=128, symmetrise=FALSE, sym.vec=c(1),
                            nts = Time, index="P5P5"){
    list <- read.table(paste(filelist[1], "/", correlator, sep=""), comment.char="", 
            skip=2)
#~     print(head(list, n=150))
    starts <- seq(1, length(list$V1), by=Time+1)
    spin_combs <- list[starts, 2]
#~     print(spin_combs)
#~     print(starts)
#~     print(list[starts, ])
    blanklineoffset=seq(1, length(spin_combs))
    
    tmp <- array(seq(1, (Time)*2*length(filelist), by=1), dim=c(2*length(filelist), Time))
    ## how exactly are the data in cf arranged: How do they have to be arranged to correspond to this in tmp?
    #~ ar <- array(seq(1, 15), dim=c(3, 5))
    #~ print(ar)
    #~ ar[, 1:2] <- c(20, 21, 22, 23, 24, 25)
    #~ print(ar)
#~     for( index in spin_combs[1:4]){
#~         print(index)
#~     print(starts[spin_combs==index])
#~     print(spin_combs==index)
    for(filenumber in 1:length(filelist)){
        tmp2 <- read.table(paste(filelist[filenumber], "/", correlator, sep=""), comment.char="", 
            skip=3+starts[spin_combs==index]+blanklineoffset[spin_combs==index], nrows=Time)
        tmp[c(2*filenumber-1, 2*filenumber), ] <- t(as.matrix(tmp2))
#~     print(tmp2)
#~     }
    }
#~     print(dim(tmp))
#~     tmp <- t(tmp)
#~     print(dim(tmp))
    
    
    total_nts <- length(filelist) # nts*length(smear_combs_to_read)*nrow(combs_to_read)
    
#~     realcols <- seq(1,2*total_nts,2)
#~     imagcols <- seq(2,2*total_nts,2)
    
    realcols <- seq(1,2*total_nts,2)
    imagcols <- seq(2,2*total_nts,2)
    
    cf <- cf_meta(nrObs = 1, Time=Time, nrStypes = 1, symmetrised = FALSE)
    cf <- cf_orig(cf, cf = tmp[realcols,,drop=FALSE], icf = tmp[imagcols,,drop=FALSE])
    
#~     print(cf$cf)
    
    if(symmetrise){
        # in some cases it makes sense to store only a subset of the time slices of a
        # correlation function. In this case, symmetrisation is not possible unless
        # the missing time slices are reconstructed or added manually somehow.
        if( nts != Time ){
        stop("The time extent and the number of time slices in the correlator do not agree, cannot symmetrise!")
        }
#~         cf <- isym(cf, sym.vec)
        cf <- symmetrise.cf(cf, sym.vec)
    }
    return (invisible(cf))
}

######################################################
######################################################
######################################################
######################################################
######################################################
######################################################


filelist <- getorderedfilelist(path="/home/gross/Documents/heavymesons/out_th0/", 
        basename="b", ending="/mes_contr_H_C_S_H")
print(length(filelist))
filelist2 <- getorderedfilelist(path="/home/gross/Documents/heavymesons/out_th0/", 
        basename="a", ending="/mes_contr_H_C_S_H")
print(length(filelist2))
filelist <- c(filelist, filelist2)

filelist <- substr(filelist, 1, nchar(filelist[1])-17)
#~ filenumbers <- substr(filelist, nchar(filelist[1])-4, nchar(filelist[1])-1)
#~ print(filenumbers)

print(paste(filelist[1], "/mes_contr_H_C_S_H", sep=""))

## P5P5 = C(t)

Time <- 128
index <- "P5P5"

scalar <- read_in_single_spin(filelist, symmetrise=TRUE, index=index)
scalar <- bootstrap.cf(scalar, boot.R=20)
#~ print(class(wl))
#~ summary(wl)
#~ print(wl)

effmass <- bootstrap.effectivemass(scalar, type="solve")

pdf("tryreadin.pdf", title="")

plot(scalar, log="y")
plot(effmass)

## AA-VA-AV+VV for C_mumu=C_00
Time <- 128
index <- "A0P5"

#~ aa <- read_in_single_spin(filelist, symmetrise=TRUE, index="A0P5", correlator="mes_contr_H_C_Dth0_A_C_P_H_H_S_H")
#~ aa <- bootstrap.cf(aa, boot.R=20)

#~ va <- read_in_single_spin(filelist, symmetrise=TRUE, index="V0P5", correlator="mes_contr_H_C_Dth0_V_C_P_H_H_S_H")
#~ va <- bootstrap.cf(va, boot.R=20)

#~ av <- read_in_single_spin(filelist, symmetrise=TRUE, index="A0P5", correlator="mes_contr_H_C_Dth0_A_C_P_H_H_S_H")
#~ av <- bootstrap.cf(av, boot.R=20)

#~ vv <- read_in_single_spin(filelist, symmetrise=TRUE, index="V0P5", correlator="mes_contr_H_C_Dth0_V_C_P_H_H_S_H")
#~ vv <- bootstrap.cf(vv, boot.R=20)


#~ plusterms <- add.cf(aa, vv)
#~ minusterms <- add.cf(av, va)
#~ czerozero <- add.cf(plusterms, minusterms, a=1., b=-1.)
#~ effmass <- bootstrap.effectivemass(czerozero, type="log")

#~ plot(czerozero, log="y")
#~ plot(effmass)

aa <- readnissatextcf(filelist, smear_combs_to_read=c("mes_contr_H_C_Dth0_A_C_P_H_H_S_H"), Time=128, 
        corrtype="newcorr", combs_to_read=data.frame(op1_idx="C_H", op2_idx="Dth0_A0_C_P_H_H_S_H", spin_comb="A0P5"),
        symmetrise=TRUE, sym.vec=c(1))
aa <- bootstrap.cf(aa, boot.R=20)

av <- readnissatextcf(filelist, smear_combs_to_read=c("mes_contr_H_C_Dth0_V_C_P_H_H_S_H"), Time=128, 
        corrtype="newcorr", combs_to_read=data.frame(op1_idx="C_H", op2_idx="Dth0_V0_C_P_H_H_S_H", spin_comb="A0P5"),
        symmetrise=TRUE, sym.vec=c(1))
av <- bootstrap.cf(av, boot.R=20)

va <- readnissatextcf(filelist, smear_combs_to_read=c("mes_contr_H_C_Dth0_A_C_P_H_H_S_H"), Time=128, 
        corrtype="newcorr", combs_to_read=data.frame(op1_idx="C_H", op2_idx="Dth0_A0_C_P_H_H_S_H", spin_comb="V0P5"),
        symmetrise=TRUE, sym.vec=c(1))
va <- bootstrap.cf(va, boot.R=20)

vv <- readnissatextcf(filelist, smear_combs_to_read=c("mes_contr_H_C_Dth0_V_C_P_H_H_S_H"), Time=128, 
        corrtype="newcorr", combs_to_read=data.frame(op1_idx="C_H", op2_idx="Dth0_V0_C_P_H_H_S_H", spin_comb="V0P5"),
        symmetrise=TRUE, sym.vec=c(1))
vv <- bootstrap.cf(vv, boot.R=20)


plot(aa, log="y", main="AA", xlab="t/a", ylab="A_0A_0")
plot(av, log="y", main="AV", xlab="t/a", ylab="A_0V_0")
plot(va, log="y", main="VA", xlab="t/a", ylab="V_0A_0")
plot(vv, log="y", main="VV", xlab="t/a", ylab="V_0V_0")

plusterms <- add.cf(aa, vv)
minusterms <- add.cf(av, va)
czerozero <- add.cf(plusterms, minusterms, a=1., b=-1.)
effmass <- bootstrap.effectivemass(czerozero, type="log")

plot(czerozero, log="y")
plot(effmass)
