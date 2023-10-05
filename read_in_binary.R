## Read in binary files with results from DGammaDq2.bin

#~ writing functions:


#~ void jack_store(FILE *ofp,jack_t *thisjack)
#~ {
#~    fwrite(&(thisjack->boot),sizeof(int),1,ofp);
#~    fwrite(&(thisjack->n),sizeof(int),1,ofp);
#~    fwrite(thisjack->f,sizeof(double),thisjack->n+4,ofp);
#~ }

read_jack <- function(to.read, endian="little") {
    boot     <- readBin(to.read, "integer", 1, endian)
    n        <- readBin(to.read, "integer", 1, endian)
    bsamples <- readBin(to.read, "double", n+4, endian)
    return(list(boot=boot, n=n, bsamples = bsamples))
}

#~ void store_multi_solve(FILE* ofp,int nk,hlt_spectre_t *spectre,int nsteps)
#~ {
#~    int i,n,nmax,m,nm;

#~    if(nk<1)
#~       error("nk must be >=1");   
   
#~    for(i=0;i<nk;++i)
#~    {
#~       if(spectre+i==NULL)
#~ 	 error("invalid hlt_spectre");   
   
#~       n=spectre[i].n-nsteps;
#~       if(n<0 && nsteps>0)
#~ 	 error("nsteps appears to be too large");   
#~    }
#~    nm=spectre[nk-1].rho[n].n;

#~    fwrite(&(nk),sizeof(int),1,ofp);
#~    fwrite(&(nsteps),sizeof(int),1,ofp);
#~    fwrite(&(nm),sizeof(int),1,ofp);
   
#~    for(i=0;i<nk;++i)
#~    {
#~       fwrite(&(spectre[i].sys),sizeof(double),1,ofp);

#~       if(nsteps>0)
#~       {
#~ 	 nmax=nsteps;
#~ 	 n=spectre[i].n-nsteps;
#~       }
#~       else
#~       {
#~ 	 nmax=spectre[i].n;
#~ 	 n=0;
#~       }
#~       fwrite(&(nmax),sizeof(int),1,ofp);
      
#~       for(m=0;m<nmax;++m)
#~       {
#~ 	 fwrite(spectre[i].lambda+n+m,sizeof(double),1,ofp);
#~ 	 fwrite(spectre[i].Bnorm+n+m,sizeof(double),1,ofp);
#~ 	 fwrite(spectre[i].A0ABCW+NFUNC*(n+m),sizeof(double),NFUNC,ofp);
#~ 	 fwrite(spectre[i].A0ABCW_ref+NFUNC*(n+m),sizeof(double),NFUNC,ofp);

#~ 	 jack_store(ofp,spectre[i].rho+n+m);
#~       }
#~    }
#~ }

read_multi_solve <- function(to.read, endian="little", NFUNC = NFUNC) {
    tmp            <- list(sys=NA, nmax=NA, lambda=list(), Bnorm=list(), 
                            A0ABCW=list(), A0ABCW_ref=list(), boottype=list(), bootnumber=list(), 
                            DGDq2mean=list(), DGDq2meanwobias=list(), DGDq2sd=list(), 
                            DGDq2bias=list(), DGDq2boots=list())
    tmp$nk         <- readBin(to.read, "integer", 1, endian)
    tmp$nsteps     <- readBin(to.read, "integer", 1, endian)
    tmp$nm         <- readBin(to.read, "integer", 1, endian)
    for ( i in seq(1, tmp$nk)) {
        tmp$sys[i]               <- readBin(to.read, "double",  1, endian)
        tmp$nmax[i]              <- readBin(to.read, "integer", 1, endian)
        tmp$lambda[[i]]          <- rep(NA, tmp$nmax[i])
        tmp$Bnorm[[i]]           <- rep(NA, tmp$nmax[i])
        tmp$A0ABCW[[i]]          <- array(NA, dim=c(tmp$nmax[i], NFUNC))
        tmp$A0ABCW_ref[[i]]      <- array(NA, dim=c(tmp$nmax[i], NFUNC))
        tmp$boottype[[i]]        <- rep(NA, tmp$nmax[i])
        tmp$bootnumber[[i]]      <- rep(NA, tmp$nmax[i])
        tmp$DGDq2mean[[i]]       <- rep(NA, tmp$nmax[i])
        tmp$DGDq2meanwobias[[i]] <- rep(NA, tmp$nmax[i])
        tmp$DGDq2sd[[i]]         <- rep(NA, tmp$nmax[i])
        tmp$DGDq2bias[[i]]       <- rep(NA, tmp$nmax[i])
        tmp$DGDq2boots[[i]]      <- list()
        
        for(n in seq(1, tmp$nmax[i])) {
            tmp$lambda[[i]][n]          <- readBin(to.read, "double", 1,     endian)
            tmp$Bnorm[[i]][n]           <- readBin(to.read, "double", 1,     endian)
            tmp$A0ABCW[[i]][n, ]        <- readBin(to.read, "double", NFUNC, endian)
            tmp$A0ABCW_ref[[i]][n, ]    <- readBin(to.read, "double", NFUNC, endian)
            boottmp                     <- read_jack(to.read, endian)
            tmp$boottype[[i]][n]        <- boottmp$boot
            tmp$bootnumber[[i]][n]      <- boottmp$n
            tmp$DGDq2mean[[i]][n]       <- boottmp$bsamples[(boottmp$n-3)]
            tmp$DGDq2meanwobias[[i]][n] <- boottmp$bsamples[(boottmp$n-2)]
            tmp$DGDq2sd[[i]][n]         <- boottmp$bsamples[(boottmp$n-1)]
            tmp$DGDq2bias[[i]][n]       <- boottmp$bsamples[boottmp$n]
            tmp$DGDq2boots[[i]][[n]]    <- boottmp$bsamples[1:(boottmp$n-4)]
        }
    }
    return(tmp)
}


#~ input parameters:
endian <- "little"
NCOMBS <- 5
NDG <- 4
NFUNC <- 6
filename <- "/home/gross/Documents/heavymesons/data/cB211ab.07.64/t48_64/outputDGammaDq2/DGammaDq2.bin"
filename <- "/home/gross/Documents/heavymesons/scripts_alessandro/DGammaDq2.bin"
to.read <- file(filename, "rb")

#~ void store_DG(FILE *ofp,int nsets, DG_t *DG,int nsteps)
#~ {
#~    int iset,idg,n,idx;
   
#~    fwrite(&(nsets),sizeof(int),1,ofp);

#~    for(iset=0;iset<nsets;++iset)
#~    {
#~       fwrite(&(DG[iset].L),sizeof(int),1,ofp);
#~       fwrite(&(DG[iset].T),sizeof(int),1,ofp);
#~       fwrite(&(DG[iset].tsinktj2),sizeof(int),1,ofp);
#~       fwrite(&(DG[iset].tj1tsrc),sizeof(int),1,ofp);

#~       fwrite(&(DG[iset].afm),sizeof(double),1,ofp);
#~       fwrite(&(DG[iset].w),sizeof(double),1,ofp);
#~       fwrite(&(DG[iset].mpigev),sizeof(double),1,ofp);
#~       fwrite(&(DG[iset].mH),sizeof(double),1,ofp);
#~       fwrite(&(DG[iset].mL),sizeof(double),1,ofp);
#~       fwrite(&(DG[iset].eL),sizeof(double),1,ofp);

#~       fwrite(&(DG[iset].neps),sizeof(int),1,ofp);
#~       fwrite(&(DG[iset].nnorms),sizeof(int),1,ofp);

#~       n=NDG*NCOMBS*DG[iset].neps;
#~       for(idg=0;idg<n;++idg)
#~       {
#~ 	 idx=idg*DG[iset].nnorms;
#~ 	 store_multi_solve(ofp,
#~ 			   DG[iset].nnorms,
#~ 			   DG[iset].spectre+idx,
#~ 			   nsteps);
#~       }
#~    }
#~ }
nsets <- readBin(to.read, "integer", 1, endian)
lambdas <- c()
alldata <- list()

for ( iset in seq(1, nsets)) {
    metadata          <- list()
    
    metadata$L        <- readBin(to.read, "integer", 1, endian)
    metadata$NT       <- readBin(to.read, "integer", 1, endian)
    metadata$tsinktj2 <- readBin(to.read, "integer", 1, endian)
    metadata$tj1tsrc  <- readBin(to.read, "integer", 1, endian)
    metadata$afm      <- readBin(to.read, "double",  1, endian)
    metadata$w        <- readBin(to.read, "double",  1, endian)
    metadata$mpigev   <- readBin(to.read, "double",  1, endian)
    metadata$mH       <- readBin(to.read, "double",  1, endian)
    metadata$mL       <- readBin(to.read, "double",  1, endian)
    metadata$eL       <- readBin(to.read, "double",  1, endian)
    metadata$neps     <- readBin(to.read, "integer", 1, endian)
    metadata$nnorm    <- readBin(to.read, "integer", 1, endian)
    for(idg in seq(1, NDG*NCOMBS*metadata$neps)){
        print(paste("idg=", idg))
        alldata[[idg]] <- read_multi_solve(to.read, endian=endian, NFUNC=NFUNC)
#~         print(data)
        lambdas <- append(x=lambdas, values=c(alldata[[idg]]$lambda[[1]][1], alldata[[idg]]$lambda[[2]][1]))
    }
    
}
print(readBin(to.read, "integer", 1, endian))

sort( sapply(ls(),function(x){object.size(get(x))})) 
