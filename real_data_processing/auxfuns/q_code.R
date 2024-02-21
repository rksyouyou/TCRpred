

disseq <- function(seq,q){
    nseq = nchar(seq)
    out = NULL
    if(nseq>q){
        startp = 1:(nseq-(q-1))
        endp = q:nseq
        out = substring(seq,startp,endp)
    } 
    return(sort(out))
}


qcode <- function(sid,aaSeq,abundance,q=3,remove.ht = TRUE,keep.all=FALSE){
    nAA = length(AA_STANDARD)
    usid = unique(sid)
    np = length(usid)    
    qc = apply(expand.grid(replicate(q, AA_STANDARD, simplify=FALSE)),1,paste,collapse="")
    if(q == 2) qc[which(qc=="NA")] = 'xNA'
    qmat = matrix(0,nrow=np,ncol=nAA^q)
    rownames(qmat) = usid
    colnames(qmat) = qc
    if(remove.ht){
        ## remove head and tail 
        newSeq = substr(aaSeq,5,nchar(aaSeq)-2)
    }
    for(i in 1:np){
        pick = which(sid==usid[i])
        qout = unlist(rep(lapply(newSeq[pick],disseq,q=q),abundance[pick]))
        if(q==2|"NA"%in%qout)  qout[qout=='NA'] = 'xNA'
        tout = table(qout)
        qmat[i,names(tout)] = tout
    }
    if(!keep.all) qmat = qmat[,colSums(qmat==0)<np] # remove cols with all 0s
    return(qmat)
}

