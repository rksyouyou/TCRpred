## survival + covariates + VDJ info
### screen by TCR data 
## TCR
## 1. remove out_of_frame/stop codon/special character
## 2. remove seq with abundance 1 (as the author suggested)
## 3. remove patient with 1 seq
refm = 'PAM250'
## Due to the TCGA policy, we are not permitted to upload following data online. The TCR_raw can be obtained by processing the RNAseq data by using TRUST4. The homoMat can be constructed by the TCRhom method. More details can be found in the manuscript. 
TCR_raw = readRDS('TRUST4_V32_TRBC_filtered.rds')
homoMat = readRDS(paste0('similar_matrix_samples_',refm,'.rds'))
dset = c("LUSC","LUAD")
N = 8

for(dtype in dset){
    ## load data: survival 
    con = gzcon(url(paste0('https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-',dtype,'.survival.tsv')))
    txt = readLines(con)
    sdat = read.csv(textConnection(txt),fill=TRUE,sep='\t',check.names = FALSE,row.names=1)
    close(con)
    ## load data: covariantes
    con = gzcon(url(paste0('https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-',dtype,'.GDC_phenotype.tsv.gz')))
    txt = readLines(con)
    pdat = read.csv(textConnection(txt),fill=TRUE,sep='\t',check.names = FALSE,row.names=1)
    close(con)
    ## data pre
    pdat1 = pdat[,c('submitter_id','age_at_index.demographic','gender.demographic',"tumor_stage.diagnoses")]
    colnames(sdat) = c('OS','submitter_id','OStime')
    colnames(pdat1) = c('submitter_id','age','gender','stage')
    ## remove missing 
    pdat1 = pdat1[!(is.na(pdat1$age)|!(pdat1$gender%in%c('female','male'))),]
    ## combine stages
    pdat1$combine_stage = rep(NA,nrow(pdat1))
    tmp = substr(pdat1$stage,1,8)
    id4 = which(tmp=='stage iv')
    pdat1$combine_stage[which(tmp=='stage iv')] = 'stage iv'
    tmp = substr(pdat1$stage,1,9)
    id3 = which(tmp=='stage iii')
    pdat1$combine_stage[id3] = 'stage iii'
    tmp = substr(pdat1$stage,1,8)
    id2 = which(tmp=='stage ii')
    pdat1$combine_stage[setdiff(id2,id3)] = 'stage ii'
    tmp = substr(pdat1$stage,1,7)
    id1 = which(tmp=='stage i')
    pdat1$combine_stage[setdiff(setdiff(id1,id2),id4)] = 'stage i'
    pdat1$combine_stage[is.na(pdat1$combine_stage)] = 'not report'
    ## remove duplicates 
    sdat1 = sdat[which(!duplicated(sdat$submitter_id)),]
    pdat1 = pdat1[which(!duplicated(pdat1$submitter_id)),]
    ## merge 
    dat = merge(sdat1,pdat1,by='submitter_id')
    ##
    temp = TCR_raw[which(TCR_raw$pid == paste0('TCGA-',dtype)),c('Aliquote_id','ParticipantBarcode','X.count','CDR3aa','V','D','J')]
    ## beta chain 
    tmp = table(temp$Aliquote_id)
    temp = temp[temp$Aliquote_id%in%names(which(tmp>=N)),]
    ## 
    ## 01A primary solid tumor
    ## 06A/06B: metastatic
    ## 07A: additionalmetastatic
    ## 11A: solid tumor normal
    ##
    tutype = do.call(rbind,strsplit(temp[,1],split = '-'))[,4]
    tout = xtabs(~temp$ParticipantBarcode+tutype)
    checkid = which(rowSums(tout!=0)>1)
    rid = NULL
    if(length(checkid)>0){
        for(i in names(checkid)){
            aid = which(temp$ParticipantBarcode==i)
            tmp1 = tutype[temp$ParticipantBarcode==i]
            for(j in c("01A", "06A", "06B", "07A", "11A")) if(j%in%tmp1) break
            rid = c(rid,aid[which(tmp1!=j)])
        }
    }
    ## same tumor type, different samples.
    ## the similarity matrix is calculated by Aliquote_id not ParticipantBarcode
    ## thus duplicated Aliquote_id samples should be removed
    if(length(rid)>0) temp1 = temp[-rid,] else temp1 = temp
    rid = NULL
    if(length(unique(temp1$Aliquote_id))>length(unique(temp1$ParticipantBarcode))){
        subAid = substr(unique(temp1$Aliquote_id),1,12)
        checkid = subAid[which(duplicated(subAid))]
        for(i in checkid){
            tmp1 = temp1[which(temp1$ParticipantBarcode==i),]
            keepid = names(which.max(table(tmp1$Aliquote_id)))
            rid = c(rid,which(temp1$ParticipantBarcode == i&temp1$Aliquote_id != keepid))
        }
    }
    if(length(rid)>0) temp1 = temp1[-rid,]
    indv = intersect(unique(temp1$ParticipantBarcode),dat$submitter_id)
    dat1 = dat[dat$submitter_id%in%indv,] # response & covariate 
    TCR_beta = temp1[temp1$ParticipantBarcode%in%indv,] # TCR
    ##
    pid = unique(TCR_beta$Aliquote_id)
    mat = homoMat[pid,pid]
    rnmat = rownames(mat)
    new_rnmat = substr(rnmat,1,12)
    if(!all(new_rnmat%in%indv)) stop('Error')
    rownames(mat) = colnames(mat) = new_rnmat
    ## shannon entropy
    sum_abundance = aggregate(TCR_beta$X.count,by=list(TCR_beta$ParticipantBarcode),FUN="sum")
    colnames(sum_abundance) = c('ParticipantBarcode','sum_abund')
    TCR_beta= merge(sum_abundance,TCR_beta,by='ParticipantBarcode')
    TCR_beta$q = TCR_beta$X.count/TCR_beta$sum_abund
    TCR_beta$qln1 = -TCR_beta$q*log(TCR_beta$q)
    shannon<-aggregate(TCR_beta$qln1,by=list(TCR_beta$ParticipantBarcode),FUN="sum")
    names(shannon)<-c("submitter_id","shannon")
    pheno = merge(dat1,shannon,by='submitter_id')
    ## make VDJ matrix
    V_uniq = unique(TCR_beta$V)
    D_uniq = unique(TCR_beta$D)
    J_uniq = unique(TCR_beta$J)
    V_len = length(V_uniq)
    D_len = length(D_uniq)
    J_len = length(J_uniq)
    nsubj = length(unique(TCR_beta$ParticipantBarcode))
    Vmat = matrix(0,nsubj,V_len)
    Dmat = matrix(0,nsubj,D_len)
    Jmat = matrix(0,nsubj,J_len)
    colnames(Vmat) = V_uniq
    colnames(Dmat) = D_uniq
    colnames(Jmat) = J_uniq
    rownames(Vmat) = rownames(Dmat) = rownames(Jmat) = pheno$submitter_id
    ## 
    for(i in pheno$submitter_id){
        temp = TCR_beta[which(TCR_beta$ParticipantBarcode==i),]
        for(j in 1:nrow(temp)){
            vname = temp$V[j]
            Vmat[i,vname] = Vmat[i,vname] + temp[j,'X.count']
            ##
            dname = temp$D[j]
            Dmat[i,dname] = Dmat[i,dname] + temp[j,'X.count']
            ##
            jname = temp$J[j]
            Jmat[i,jname] = Jmat[i,jname] + temp[j,'X.count']
        }
    }
    if(any(colnames(Vmat)=='.')) Vmat = Vmat[,-which(colnames(Vmat)=='.')]
    if(any(colnames(Dmat)=='.')) Dmat = Dmat[,-which(colnames(Dmat)=='.')]
    if(any(colnames(Jmat)=='.')) Jmat = Jmat[,-which(colnames(Jmat)=='.')]
    VDJ = cbind(Vmat,Dmat,Jmat)
    ## 
    K = mat
    pheno = pheno[match(rownames(K),pheno[,1]),]
    VDJ = VDJ[match(rownames(K),rownames(VDJ)),]
    ## 
    out = list(pheno = pheno,TCRB = TCR_beta,K = mat,VDJ=VDJ)
    saveRDS(out,paste0(dtype,'_surv_covariate_TCR_N',N,'_bin_PAM250.rds'))
}

####

source('auxfuns/q_code.R')
source('auxfuns/psd.R')
library(Biostrings) ## AA standard



## dynamic cut 
## prepare the data
cutinfo_out = NULL
for(dtype in dset){
    dat = readRDS(paste0(path,'dat/',dtype,'_surv_covariate_TCR_N',N,'_bin_',refm,'.rds'))
    ##
    pid = which(dat$pheno$combine_stage=='stage i')
    pheno = dat$pheno[pid,]
    K = dat$K[pid,pid]
    VDJ = dat$VDJ[pid,] # VDJ matrix 
    TCRB = dat$TCRB[dat$TCRB$ParticipantBarcode%in%rownames(K),]
    ## mk binary outcome
    ## decide the cutoff
    cutvec = quantile(sort(pheno$OStime),seq(0.2,0.8,0.05))
    ## 
    rmat = NULL
    for(cutoff in cutvec){
        rid = which((pheno$OStime<cutoff&(pheno$OS==0)))
        if(length(rid)>0) pheno_tmp = pheno[-rid,]
        pheno_tmp$LSy = as.numeric(pheno_tmp$OStime>=cutoff)
        rmat = rbind(rmat,c(cutoff,mean(pheno_tmp$LSy),sum(pheno_tmp$LSy==0),sum(pheno_tmp$LSy==1),nrow(pheno_tmp)))
    }
    pick_cut = which.min(abs(rmat[,2]-0.5))
    ## 
    cutinfo = matrix(rmat[pick_cut,-2],nrow=1)
    colnames(cutinfo) = c('cutoff','short','long','sample.size')
    rownames(cutinfo) = paste0(dtype,'-N',N,'-v',v)
    cutinfo_out = rbind(cutinfo_out,cutinfo)
    cutoff = cutvec[pick_cut]
    if(abs(rmat[pick_cut,2]-0.5)<0.15){ #
        rid = which(pheno$OStime<cutoff&(pheno$OS==0))
        if(length(rid)>0) pheno = pheno[-rid,]
        pheno$LSy = as.numeric(pheno$OStime>=cutoff)
        ## similarity matrix
        K = K[pheno$submitter_id,pheno$submitter_id]
        ## mk Z matrix, match order
        Zmat = qcode(TCRB$ParticipantBarcode,TCRB$CDR3aa,TCRB$X.count,q=3)
        if(!identical(rownames(Zmat),rownames(K)))
            Zmat = Zmat[match(rownames(K),rownames(Zmat)),]
        Zmat = Zmat[,colSums(Zmat==0)<nrow(pheno)]
        ##
        VDJ = VDJ[pheno$submitter_id,]
        ##
        Y = matrix(pheno$LSy,ncol=1)
        vset = c('age','gender','shannon')
        X= model.matrix(~.,pheno[,vset]) # include intercept, model matrix 
        K[is.na(K)] = 0
        K = (K + t(K))
        diag(K)= 1
        K = psd(K,M=nrow(K))
        ## compute m
        es = eigen(K)
        Q <- es$values
        m <- sum(Q>=0)
        ##
        indat = list(Y=Y,X=X,K=K,Z=Zmat)
        indat_VDJ = list(Y=Y,X=X,K=K,Z=cbind(Zmat,VDJ))
        saveRDS(indat,paste0(dtype,'_N',N,'_s1v',v,'_',refm,'.rds'))
        saveRDS(indat_VDJ,paste0(dtype,'_N',N,'_VDJ_s1v',v,'_',refm,'.rds'))
    }
}



#### compute K2 #####
dtype1 = 'LUAD'
dtype2 = 'LUSC'
homomat = readRDS('similar_matrix_samples_PAM250.rds')
n = nrow(homomat)

for(i in 1:(n-1))
    for(j in (i+1):n)
        homomat[i,j] = homomat[j,i]


rawdat1 = readRDS(paste0(dtype1,'_surv_covariate_TCR_N8_bin_PAM250.rds'))
rawdat2 = readRDS(paste0(dtype2,'_surv_covariate_TCR_N8_bin_PAM250.rds'))

grp1 = unique(rawdat1$TCRB$ParticipantBarcode)
grp2 = unique(rawdat2$TCRB$ParticipantBarcode)

aid1 = rawdat1$TCRB$Aliquote_id
aid2 = rawdat2$TCRB$Aliquote_id
names(aid1) = rawdat1$TCRB$ParticipantBarcode
names(aid2) = rawdat2$TCRB$ParticipantBarcode

##
tmp = homomat[aid1[grp1],aid2[grp2]]
rownames(tmp) = grp1
colnames(tmp) = grp2
saveRDS(tmp,paste0('homo_',dtype1,'_',dtype2,'_PAM250.rds'))


