eva_metric <- function(Y,Y_pred){
    ## Y: true value
    ## Y_pred: continuous, predicted probability

    if(all(Y%in%c(0,1))) rtype = 'bin' else rtype = 'cont'
    if(rtype == 'bin'){
        Ybin = as.numeric(Y_pred>=0.5)
        out = xtabs(~Ybin+Y)
        if(nrow(out)==0) {
            out = rep(NA,4)
        } else {
            if(nrow(out)==1) out = rbind(c(0,0),out)
            ## TN FN
            ## FP TP
            PPV = out[2,2]/sum(out[2,]) # positive predictive value
            NPV = out[1,1]/sum(out[1,]) # negative predictive value
            out = c(PPV = PPV, NPV = NPV, clas_error = mean((Y - Ybin)^2),auc = verification::roc.area(Y,as.matrix(Y_pred))$A)
        }
    } else {
        out = c(MSE = mean((Y-Y_pred)^2))
    }

    return(out)
}
      


disseq <- function(seq,k){
    nseq = nchar(seq)
    out = NULL
    if(nseq>k){
        startp = 1:(nseq-(k-1))
        endp = k:nseq
        out = substring(seq,startp,endp)
    } 
    return(sort(out))
}


mk_Z <- function(sid,aaSeq,abundance,k,remove.ht = TRUE,keep.all=FALSE,normalization=TRUE){
    nAA = length(Biostrings::AA_STANDARD)
    usid = unique(sid)
    np = length(usid)    
    qc = apply(expand.grid(replicate(k, AA_STANDARD, simplify=FALSE)),1,paste,collapse="")
    if(k == 2) qc[which(qc=="NA")] = 'xNA'
    Zmat = matrix(0,nrow=np,ncol=nAA^k)
    rownames(Zmat) = usid
    colnames(Zmat) = qc
    if(remove.ht){
        ## remove head and tail 
        newSeq = substr(aaSeq,5,nchar(aaSeq)-2)
    }
    for(i in 1:np){
        pick = which(sid==usid[i])
        qout = unlist(rep(lapply(newSeq[pick],disseq,k=k),abundance[pick]))
        if(k==2|"NA"%in%qout)  qout[qout=='NA'] = 'xNA'
        tout = table(qout)
        Zmat[i,names(tout)] = tout
    }
    if(!keep.all) Zmat = Zmat[,colSums(Zmat==0)<np] # remove cols with all 0s
    if(normalization){
        sfactor <- apply(Zmat,2,function(x) quantile(x[x>0],0.75))
        Zmat = t(t(Zmat)/sfactor)
    }
    return(Zmat)
}


seqhom <- function(seq1, seq2,type){ # similarity between seqs
    s12 <- Biostrings::pairwiseAlignment(seq1, seq2, substitutionMatrix = type)
    return(s12@score)    
}


psd <- function(S){ # check psd
  if (isSymmetric(S)==FALSE){
    stop("S should be symmetric.")
  }
  es = eigen(S)
  S <- as.matrix(S)
  Q <- es$values
  m <- sum(Q>=0)
  P <- es$vectors
  n <- dim(S)[1]
  my.psd <- matrix(0, n, n)
  rownames(my.psd) <- rownames(S)
  colnames(my.psd) <- colnames(S)
  for (j in 1: m){
    my.psd <- my.psd + Q[j]*P[,j]%*%t(P[,j])
  }
  return(my.psd)
}

mk_K <- function(sid,aaSeq,abundance,refm = 'BLOSUM62'){
    ns = length(aaSeq)
    dvec = rep(NA,ns)
    for(i in 1:ns)
        dvec[i] = seqhom(aaSeq[i],aaSeq[i],refm)

    usid = unique(sid)
    nus = length(usid)
    smat = matrix(NA,nus,nus)
    rownames(smat) = colnames(smat) = usid
    for(i in 1:nus)
        for(j in i:nus){
        sdat1 = aaSeq[which(sid==usid[i])]
        sdat2 = aaSeq[which(sid==usid[j])]
        dvec1 = dvec[which(sid==usid[i])]
        dvec2 = dvec[which(sid==usid[j])]
        abund1 = abundance[which(sid==usid[i])]
        abund2 = abundance[which(sid==usid[j])]
        submat = matrix(NA,length(sdat1),length(sdat2))
        for(l in 1:length(sdat1))
            for(t in 1:length(sdat2))
                submat[l,t] = seqhom(sdat1[l],sdat2[t],refm)/sqrt(dvec1[l]*dvec2[t])
        s1 = sum(apply(submat,1,max)*abund1)
        s2 = sum(apply(submat,2,max)*abund2)
        smat[i,j] = smat[j,i] = (s1+s2)/sum(c(abund1,abund2))
        }
    K = psd(smat)
    return(K)
}
