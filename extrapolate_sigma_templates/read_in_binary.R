## Read in binary files with results from DGammaDq2.bin

#~ writing functions:


#~ void jack_store(FILE *ofp,jack_t *thisjack)
#~ {
#~    fwrite(&(thisjack->boot),sizeof(int),1,ofp);
#~    fwrite(&(thisjack->n),sizeof(int),1,ofp);
#~    fwrite(thisjack->f,sizeof(double),thisjack->n+4,ofp);
#~ }

## boottype: 0 jackknife, 1 bootstrap
## n=nboots
## n bootstrapsamples, mean, mean without bias, sd, bias
read_jack <- function(to.read, endian="little") {
    boot     <- readBin(to.read, "integer", 1, endian = endian)
    n        <- readBin(to.read, "integer", 1, endian = endian)
    bsamples <- readBin(to.read, "double", n+4, endian = endian)
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

## nk: number fo norms
## nmax: amount of times information is stored, could store partial results
## sys: systematic error, same for all norms
## lambda: tradeoff parameter
# NFUNC: data pointsfor the functional, for current norm and minimum of all norms. 
## ? Bnorm, A0, A/A0_min, A/A0, Bnorm*B, C ?
## functionals for the HLT
## DGdq2: result, differentail decay rate
read_multi_solve <- function(to.read, endian="little", NFUNC = 6) {
    res            <- list(A0ABCW=list(), A0ABCW_ref=list(), DGDq2boots=list())
    tmp            <- list(sys=NA, nmax=NA, lambda=list(), Bnorm=list(), 
                            A0ABCW=list(), A0ABCW_ref=list(), boottype=list(), bootnumber=list(), 
                            DGDq2mean=list(), DGDq2meanwobias=list(), DGDq2sd=list(), 
                            DGDq2bias=list(), DGDq2boots=list())
    resdf <- data.frame(nk=NA, nsteps=NA, nm=NA, sys=NA, nmax=NA, lambda=NA, 
                        Bnorm=NA, boottype=NA, bootnumber=NA, DGDq2mean=NA, 
                        DGDq2meanwobias=NA, DGDq2sd=NA, DGDq2bias=NA)
    tmp$nk         <- readBin(to.read, "integer", 1, endian = endian)
    tmp$nsteps     <- readBin(to.read, "integer", 1, endian = endian)
    tmp$nm         <- readBin(to.read, "integer", 1, endian = endian)
    for ( i in seq(1, tmp$nk)) {
        tmp$sys[i]               <- readBin(to.read, "double",  1, endian = endian)
        tmp$nmax[i]              <- readBin(to.read, "integer", 1, endian = endian)
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
            tmp$lambda[[i]][n]          <- readBin(to.read, "double", 1,     endian = endian)
            tmp$Bnorm[[i]][n]           <- readBin(to.read, "double", 1,     endian = endian)
            tmp$A0ABCW[[i]][n, ]        <- readBin(to.read, "double", NFUNC, endian = endian)
            tmp$A0ABCW_ref[[i]][n, ]    <- readBin(to.read, "double", NFUNC, endian = endian)
            boottmp                     <- read_jack(to.read, endian = endian)
            tmp$boottype[[i]][n]        <- boottmp$boot
            tmp$bootnumber[[i]][n]      <- boottmp$n
            tmp$DGDq2mean[[i]][n]       <- boottmp$bsamples[(boottmp$n+1)]
            tmp$DGDq2meanwobias[[i]][n] <- boottmp$bsamples[(boottmp$n+2)]
            tmp$DGDq2sd[[i]][n]         <- boottmp$bsamples[(boottmp$n+3)]
            tmp$DGDq2bias[[i]][n]       <- boottmp$bsamples[(boottmp$n+4)]
            tmp$DGDq2boots[[i]][[n]]    <- boottmp$bsamples[1:boottmp$n]
        }
    }
    return(tmp)
}

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

#~       fwrite(DG[iset].eps,sizeof(double),DG[iset].neps,ofp);
      
      
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

## NCOMBS: full, vpar, apar, vperp, aperp
## NDG: Z0, Z1, Z2, sum
## NFUNC: see store_multi_solve
read_in_DGDq2 <- function(filename, resultpath="./", endian=.Platform$endian, NCOMBS=5, NDG=4, NFUNC=6, write=FALSE, savename=-1) {
    to.read <- file(paste0(resultpath, "/", filename), "rb")
    
    nsets <- readBin(to.read, "integer", 1, endian = endian)
    
    result <- list(nsets=nsets, rep(NA, nsets))
    
    for ( iset in seq(1, nsets)) {
        metadata          <- data.frame(L=readBin(to.read, "integer", 1, endian = endian))
        metadata$T        <- readBin(to.read, "integer", 1, endian = endian)
        metadata$tsinktj2 <- readBin(to.read, "integer", 1, endian = endian)
        metadata$tj1tsrc  <- readBin(to.read, "integer", 1, endian = endian)
        metadata$afm      <- readBin(to.read, "double",  1, endian = endian)
        metadata$w        <- readBin(to.read, "double",  1, endian = endian)
        metadata$mpigev   <- readBin(to.read, "double",  1, endian = endian)
        metadata$mH       <- readBin(to.read, "double",  1, endian = endian)
        metadata$mL       <- readBin(to.read, "double",  1, endian = endian)
        metadata$eL       <- readBin(to.read, "double",  1, endian = endian)
        metadata$neps     <- readBin(to.read, "integer", 1, endian = endian)
        metadata$nnorm    <- readBin(to.read, "integer", 1, endian = endian)
#~         print(metadata)
        epsilons          <- readBin(to.read, "double",  metadata$neps, endian = endian)
        result[[iset+1]]  <- list()
        result[[iset+1]]$metadata <- metadata
        result[[iset+1]]$epsilons <- epsilons
        
        for(id in seq(0, NDG*NCOMBS*metadata$neps - 1)){
            ## ieps+DG->neps*(icomb+NCOMBS*idg)
            ieps <- id %% metadata$neps
            icomb <- ((id-ieps)/metadata$neps) %% NCOMBS
            iz    <- ((id-ieps)/metadata$neps-icomb)/NCOMBS
            name <- paste0("id", id, "ieps", ieps, "icomb", icomb, "iz", iz)
            tmpres <- read_multi_solve(to.read, endian=endian, NFUNC=NFUNC)
            result[[iset+1]][[name]] <- tmpres
        }
    }
    if(write) {
        filenamesave <- savename
        if (savename == -1) filenamesave <- paste0("DGammaDq2.RData")
        saveRDS(paste0(resultpath, "/", result), filename)
    }
    close(to.read)
    return (result)
}

#~ int indexDG(DG_t *DG,int inorm,int ieps,int icomb,int idg)
#~ {
#~    return inorm+DG->nnorms*(ieps+DG->neps*(icomb+NCOMBS*idg));
#~ }

## for reading in, nnorms should always be set to one

indexDG <- function(inorm, ieps, icomb=0, idg, nnorms=1, neps, NCOMBS=5) {
    return (inorm + nnorms*(ieps + neps * (icomb + NCOMBS * idg)))
}


#~ void store_partial_DG(FILE *ofp,int iset,int icomb,DG_t *DG,int nsteps)
#~ {
#~    int idg,ieps,idx;
   
#~    fwrite(&(iset),sizeof(int),1,ofp);
#~    fwrite(&(icomb),sizeof(int),1,ofp);

#~    fwrite(&(DG[iset].L),sizeof(int),1,ofp);
#~    fwrite(&(DG[iset].T),sizeof(int),1,ofp);
#~    fwrite(&(DG[iset].tsinktj2),sizeof(int),1,ofp);
#~    fwrite(&(DG[iset].tj1tsrc),sizeof(int),1,ofp);

#~    fwrite(&(DG[iset].afm),sizeof(double),1,ofp);
#~    fwrite(&(DG[iset].w),sizeof(double),1,ofp);
#~    fwrite(&(DG[iset].mpigev),sizeof(double),1,ofp);
#~    fwrite(&(DG[iset].mH),sizeof(double),1,ofp);
#~    fwrite(&(DG[iset].mL),sizeof(double),1,ofp);
#~    fwrite(&(DG[iset].eL),sizeof(double),1,ofp);

#~    fwrite(&(DG[iset].neps),sizeof(int),1,ofp);
#~    fwrite(&(DG[iset].nnorms),sizeof(int),1,ofp);

#~    fwrite(DG[iset].eps,sizeof(double),DG[iset].neps,ofp);
            
#~    for(idg=0;idg<NDG;++idg)
#~    {
#~       for(ieps=0;ieps<DG[iset].neps;++ieps)
#~       {
#~ 	 idx=indexDG(DG+iset,0,ieps,icomb,idg);
#~ 	 store_multi_solve(ofp,
#~ 			   DG[iset].nnorms,
#~ 			   DG[iset].spectre+idx,
#~ 			   nsteps);
#~       }
#~    }
#~ }

## Explanation see DGDq2, here only one combination of full, a, v, par, perp
read_in_DGDq2_partial <- function(filename, resultpath="./", endian=.Platform$endian, NCOMBS=5, NDG=4, NFUNC=6, write=FALSE, savename=-1) {
    to.read <- file(paste0(resultpath, "/", filename), "rb")
    

    metadata          <- data.frame(iset=readBin(to.read, "integer", 1, endian = endian))
    metadata$icomb    <- readBin(to.read, "integer", 1, endian = endian)
    metadata$L        <- readBin(to.read, "integer", 1, endian = endian)
    metadata$T        <- readBin(to.read, "integer", 1, endian = endian)
    metadata$tsinktj2 <- readBin(to.read, "integer", 1, endian = endian)
    metadata$tj1tsrc  <- readBin(to.read, "integer", 1, endian = endian)
    metadata$afm      <- readBin(to.read, "double",  1, endian = endian)
    metadata$w        <- readBin(to.read, "double",  1, endian = endian)
    metadata$mpigev   <- readBin(to.read, "double",  1, endian = endian)
    metadata$mH       <- readBin(to.read, "double",  1, endian = endian)
    metadata$mL       <- readBin(to.read, "double",  1, endian = endian)
    metadata$eL       <- readBin(to.read, "double",  1, endian = endian)
    metadata$neps     <- readBin(to.read, "integer", 1, endian = endian)
    metadata$nnorm    <- readBin(to.read, "integer", 1, endian = endian)
    print(metadata)
    epsilons          <- readBin(to.read, "double",  metadata$neps, endian = endian)
    result  <- list()
    result$metadata <- metadata
    result$epsilons <- epsilons
    
    for(idg in seq(0, NDG-1)){
        for (ieps in seq(0, neps-1)) {
        ## ieps+DG->neps*(icomb+NCOMBS*idg)
        idx <- indexDG(inorm=0, ieps=ieps, icomb=metadata$icomb, idg=idg, nnorms=metadata$nnorm, neps=metadata$neps, NCOMBS)
        name <- paste0("id", idx, "ieps", ieps, "icomb", metadata$icomb, "iz", idg)
        tmpres <- read_multi_solve(to.read, endian=endian, NFUNC=NFUNC)
        result[[name]] <- tmpres
    }
    }

    if(write) {
        filename <- savename
        if (savename == -1) filename <- paste0("DGammaDq2_partial_icomb", metadata$icomb, ".RData")
        saveRDS(result, filename)
    }
    close(to.read)
    return (result)
}

#~ int indexDGell(DGell_t *DGell,int inorm,int ieps,int iell,int idg)
#~ {
#~    return inorm+DGell->nnorms*(ieps+DGell->neps*(iell+DGell->nell*idg));
#~ }

indexDGELL <- function(inorm, ieps, iell, idg, nnorms, neps, nell) {
    inorm + nnorms*(ieps + neps*(iell + nell*idg))
}

#~ void store_DGell(FILE *ofp,int nsets, DGell_t *DGell,int nsteps)
#~ {
#~    int iset,idg,n,idx;
   
#~    fwrite(&(nsets),sizeof(int),1,ofp);

#~    for(iset=0;iset<nsets;++iset)
#~    {
#~       fwrite(&(DGell[iset].L),sizeof(int),1,ofp);
#~       fwrite(&(DGell[iset].T),sizeof(int),1,ofp);
#~       fwrite(&(DGell[iset].tsinktj2),sizeof(int),1,ofp);
#~       fwrite(&(DGell[iset].tj1tsrc),sizeof(int),1,ofp);

#~       fwrite(&(DGell[iset].afm),sizeof(double),1,ofp);
#~       fwrite(&(DGell[iset].w),sizeof(double),1,ofp);
#~       fwrite(&(DGell[iset].mpigev),sizeof(double),1,ofp);
#~       fwrite(&(DGell[iset].mH),sizeof(double),1,ofp);
#~       fwrite(&(DGell[iset].mL),sizeof(double),1,ofp);
#~       fwrite(&(DGell[iset].eL),sizeof(double),1,ofp);

#~       fwrite(&(DGell[iset].nell),sizeof(int),1,ofp);
#~       fwrite(&(DGell[iset].neps),sizeof(int),1,ofp);
#~       fwrite(&(DGell[iset].nnorms),sizeof(int),1,ofp);

#~       fwrite(&(DGell[iset].reps),sizeof(double),1,ofp);
#~       fwrite(DGell[iset].eps,sizeof(double),DGell[iset].neps,ofp);

#~       fwrite(DGell[iset].ell,sizeof(double),DGell[iset].nell,ofp);
      
      
#~       n=NDGELL*DGell[iset].nell*DGell[iset].neps;
#~       for(idg=0;idg<n;++idg)
#~       {
#~ 	 idx=idg*DGell[iset].nnorms;
#~ 	 store_multi_solve(ofp,
#~ 			   DGell[iset].nnorms,
#~ 			   DGell[iset].spectre+idx,
#~ 			   nsteps);
#~       }
#~    }
#~ }

## NDGELL: Y1 to Y5 plus sum
read_in_DGell <- function(filename, resultpath="./", endian=.Platform$endian, NDGELL=6, NFUNC=6, write=FALSE, savename=-1) {
    to.read <- file(paste0(resultpath, "/", filename), "rb")
    
    nsets <- readBin(to.read, "integer", 1, endian = endian)
    
    result <- list(nsets=nsets, rep(NA, nsets))
    
    for ( iset in seq(1, nsets)) {
        metadata          <- data.frame(L=readBin(to.read, "integer", 1, endian = endian))
        metadata$T        <- readBin(to.read, "integer", 1, endian = endian)
        metadata$tsinktj2 <- readBin(to.read, "integer", 1, endian = endian)
        metadata$tj1tsrc  <- readBin(to.read, "integer", 1, endian = endian)
        metadata$afm      <- readBin(to.read, "double",  1, endian = endian)
        metadata$w        <- readBin(to.read, "double",  1, endian = endian)
        metadata$mpigev   <- readBin(to.read, "double",  1, endian = endian)
        metadata$mH       <- readBin(to.read, "double",  1, endian = endian)
        metadata$mL       <- readBin(to.read, "double",  1, endian = endian)
        metadata$eL       <- readBin(to.read, "double",  1, endian = endian)
        metadata$nell     <- readBin(to.read, "integer", 1, endian = endian)
        metadata$neps     <- readBin(to.read, "integer", 1, endian = endian)
        metadata$nnorm    <- readBin(to.read, "integer", 1, endian = endian)
        metadata$reps     <- readBin(to.read, "double",  1, endian = endian)
        print(metadata)
        epsilons          <- readBin(to.read, "double",  metadata$neps, endian = endian)
        ells              <- readBin(to.read, "double",  metadata$nell, endian = endian)
        result[[iset+1]]  <- list()
        result[[iset+1]]$metadata <- metadata
        result[[iset+1]]$epsilons <- epsilons
        result[[iset+1]]$ells     <- ells
        
        for(id in seq(0, NDGELL*metadata$nell*metadata$neps - 1)){
            ## inorm+DGell->nnorms*(ieps+DGell->neps*(iell+DGell->nell*idg))
            ## norm = 0 
            ## ieps+DGell->neps*(iell+DGell->nell*idg))
            ieps <- id %% metadata$neps
            iell <- ((id-ieps)/metadata$neps) %% metadata$nell
            icomb    <- ((id-ieps)/metadata$neps-iell)/NDGELL
            name <- paste0("id", id, "ieps", ieps, "iell", iell, "icomb", icomb)
            tmpres <- read_multi_solve(to.read, endian=endian, NFUNC=NFUNC)
            result[[iset+1]][[name]] <- tmpres
        }
    }
    if(write) {
        filename <- savename
        if (savename == -1) filename <- paste0("DGammaDellDq2.RData")
        saveRDS(result, filename)
    }
    close(to.read)
    return (result)
}



#~ void store_partial_DGell(FILE *ofp,int iset,int iell,DGell_t *DGell,int nsteps)
#~ {
#~    int idg,ieps,idx;
   
#~    fwrite(&(iset),sizeof(int),1,ofp);
#~    fwrite(&(iell),sizeof(int),1,ofp);

#~    fwrite(&(DGell[iset].L),sizeof(int),1,ofp);
#~    fwrite(&(DGell[iset].T),sizeof(int),1,ofp);
#~    fwrite(&(DGell[iset].tsinktj2),sizeof(int),1,ofp);
#~    fwrite(&(DGell[iset].tj1tsrc),sizeof(int),1,ofp);

#~    fwrite(&(DGell[iset].afm),sizeof(double),1,ofp);
#~    fwrite(&(DGell[iset].w),sizeof(double),1,ofp);
#~    fwrite(&(DGell[iset].mpigev),sizeof(double),1,ofp);
#~    fwrite(&(DGell[iset].mH),sizeof(double),1,ofp);
#~    fwrite(&(DGell[iset].mL),sizeof(double),1,ofp);
#~    fwrite(&(DGell[iset].eL),sizeof(double),1,ofp);

#~    fwrite(&(DGell[iset].nell),sizeof(int),1,ofp);
#~    fwrite(&(DGell[iset].neps),sizeof(int),1,ofp);
#~    fwrite(&(DGell[iset].nnorms),sizeof(int),1,ofp);

#~    fwrite(&(DGell[iset].reps),sizeof(double),1,ofp);
#~    fwrite(DGell[iset].eps,sizeof(double),DGell[iset].neps,ofp);

#~    fwrite(DGell[iset].ell,sizeof(double),DGell[iset].nell,ofp);
            
#~    for(idg=0;idg<NDGELL;++idg)
#~    {
#~       for(ieps=0;ieps<DGell[iset].neps;++ieps)
#~       {
#~ 	 idx=indexDGell(DGell+iset,0,ieps,iell,idg);
#~ 	 store_multi_solve(ofp,
#~ 			   DGell[iset].nnorms,
#~ 			   DGell[iset].spectre+idx,
#~ 			   nsteps);
#~       }
#~    }
#~ }

read_in_DGDell_partial <- function(filename, resultpath="./", endian=.Platform$endian, NDGELL=6, NFUNC=6, write=FALSE, savename=-1) {
    to.read <- file(paste0(resultpath, "/", filename), "rb")
    

    metadata          <- data.frame(iset=readBin(to.read, "integer", 1, endian = endian))
    metadata$iell     <- readBin(to.read, "integer", 1, endian = endian)
    metadata$L        <- readBin(to.read, "integer", 1, endian = endian)
    metadata$T        <- readBin(to.read, "integer", 1, endian = endian)
    metadata$tsinktj2 <- readBin(to.read, "integer", 1, endian = endian)
    metadata$tj1tsrc  <- readBin(to.read, "integer", 1, endian = endian)
    metadata$afm      <- readBin(to.read, "double",  1, endian = endian)
    metadata$w        <- readBin(to.read, "double",  1, endian = endian)
    metadata$mpigev   <- readBin(to.read, "double",  1, endian = endian)
    metadata$mH       <- readBin(to.read, "double",  1, endian = endian)
    metadata$mL       <- readBin(to.read, "double",  1, endian = endian)
    metadata$eL       <- readBin(to.read, "double",  1, endian = endian)
    metadata$nell     <- readBin(to.read, "integer", 1, endian = endian)
    metadata$neps     <- readBin(to.read, "integer", 1, endian = endian)
    metadata$nnorm    <- readBin(to.read, "integer", 1, endian = endian)
    metadata$reps     <- readBin(to.read, "double",  1, endian = endian)
    print(metadata)
    epsilons          <- readBin(to.read, "double",  metadata$neps, endian = endian)
    ells              <- readBin(to.read, "double",  metadata$nell, endian = endian)
    result  <- list()
    result$metadata <- metadata
    result$epsilons <- epsilons
    
    for(idg in seq(0, NDGELL-1)){
        for (ieps in seq(0, neps-1)) {
#~     inorm + nnorms*(ieps + neps*(iell + nell*idg))
        idx <- indexDGELL(inorm=0, ieps=ieps, iell=metadata$iell, idg=idg, nnorms=metadata$nnorm, neps=metadata$neps, nell=metadata$nell)
        name <- paste0("id", idx, "ieps", ieps, "iell", metadata$iell, "icomb", idg)
        tmpres <- read_multi_solve(to.read, endian=endian, NFUNC=NFUNC)
        result[[name]] <- tmpres
    }
    }

    if(write) {
        filename <- savename
        if (savename == -1) filename <- paste0("DGammaDELL_partial_iell_", metadata$iell, ".RData")
        saveRDS(result, filename)
    }
    close(to.read)
    return (result)
}
