# heritability modeling
# write table of variance and covariance
# complement genes in twins

# housekeeping
workdir <- "/data/swe_gwas/ABZ/immunegenes/"
setwd(workdir)

libs <- c("dplyr", "psych", "ggplot2","nlme")
invisible(lapply(libs, require, character.only = TRUE))

################################ swedish twins
# load data
# add updated dx info
load("/data/swe_gwas/ABZ/RNA_GWAS/swedenclean.rdata")
dx = read.table("/data/swe_gwas/ABZ/RNA_GWAS/swedishNP_PCA.txt",header=T)
swedenclean$dx_new <- dx[match(swedenclean$StudyID,dx$IID),7] %>% as.factor()

# zygosity + tvab
zyg              <- read.table(paste0(workdir,"zygosity_alltwins.txt"),header=T)
swedenclean$tvab <- zyg[match(swedenclean$StudyID,zyg$IID),4]
swedenclean$zyg  <- zyg[match(swedenclean$StudyID,zyg$IID),3]

# phenotype (superior frontal cortical thickness)
pheno <- read.table("supfrontal_thickness.txt",header=T)
pheno$IID <- substr(pheno$ID,1,10)
swedenclean$SFthick <- pheno[match(swedenclean$StudyID,pheno$IID),"supfront"]

# input lists
# compgenes <- c("ILMN_1670305", "ILMN_1695588", "ILMN_2307740", "ILMN_2347789", 
#                "ILMN_1800540", "ILMN_1671928","ILMN_1692651", "ILMN_1746819")
# swedecompgenes <- swedenclean[,colnames(swedenclean) %in% compgenes] 
# df             <- cbind(swedenclean[,c(1,2,18567:18568)],swedecompgenes)
df             <- swedenclean %>% select(StudyID, Family, tvab, zyg, SFthick)

# MZ vs. DZ
MZ <- df[df$zyg=="MZ",]
DZ <- df[df$zyg=="DZ",]

# reshape to wide format
# 43 unique MZ FID; 34 complete pairs (pheno = 22 complete)
# 66 unique DZ FID; 50 compelte pairs (pheno = 26 complete)
MZ.wide <- reshape(MZ,idvar="Family",direction="wide",timevar="tvab")
DZ.wide <- reshape(DZ,idvar="Family",direction="wide",timevar="tvab")

# no NAs (complete pairs) 
# (will throw out people without pheno info if included in df)
MZcomp <- na.omit(MZ.wide)
DZcomp <- na.omit(DZ.wide)

# write var and covar table
# 10 columns of info
# 8 rows of genes (or 1 row for phenotype)
# vartable <- matrix(nrow=length(compgenes),ncol=10)
# for (i in 1:length(compgenes)){
#   t1gene <- paste0(compgenes[i],".1")
#   t2gene <- paste0(compgenes[i],".2")
#   vartable[i,1] <- var(MZcomp[,t1gene])
#   vartable[i,2] <- cov(MZcomp[,t1gene],MZcomp[,t2gene],
#                        use="pairwise.complete.obs")
#   vartable[i,3] <- var(MZcomp[,t2gene])
#   vartable[i,4] <- var(DZcomp[,t1gene])
#   vartable[i,5] <- cov(DZcomp[,t1gene],DZcomp[,t2gene],
#                        use="pairwise.complete.obs")
#   vartable[i,6] <- var(DZcomp[,t2gene])
#   vartable[i,7] <- length(rownames(MZcomp))
#   vartable[i,8] <- length(rownames(DZcomp))
#   vartable[i,9] <- paste0("X1",compgenes[i])
#   vartable[i,10] <- paste0("X2",compgenes[i])
# }

vartable <- matrix(nrow=1,ncol=10)
for (i in 1:1){
  t1gene <- "SFthick.1"
  t2gene <- "SFthick.2"
  vartable[i,1] <- var(MZcomp[,t1gene])
  vartable[i,2] <- cov(MZcomp[,t1gene],MZcomp[,t2gene],
                       use="pairwise.complete.obs")
  vartable[i,3] <- var(MZcomp[,t2gene])
  vartable[i,4] <- var(DZcomp[,t1gene])
  vartable[i,5] <- cov(DZcomp[,t1gene],DZcomp[,t2gene],
                       use="pairwise.complete.obs")
  vartable[i,6] <- var(DZcomp[,t2gene])
  vartable[i,7] <- length(rownames(MZcomp))
  vartable[i,8] <- length(rownames(DZcomp))
  vartable[i,9] <- paste0("X1","SFthickness")
  vartable[i,10] <- paste0("X2","SFthickness")
}

vartable <- as.data.frame(vartable)
vartable[,1:6] <- lapply(vartable[,1:6],as.character)
vartable[,1:6] <- lapply(vartable[,1:6],as.numeric)
names(vartable) <- c("MZtwin1var","MZcovar","MZtwin2var",
                     "DZtwin1var","DZcovar","DZtwin2var",
                     "NUM_PAIRS_MZ","NUM_PAIRS_DZ","T1_VAR","T2_VAR")

# zero out negative covariance
# vartable$MZcovar <- ifelse(vartable$MZcovar < 0,0,vartable$MZcovar)
# vartable$DZcovar <- ifelse(vartable$DZcovar < 0,0,vartable$DZcovar)

# write csv file
# write.table(vartable,paste0(workdir,"comp_h2/complement_covar.csv"),
#             sep=",",col.names=T,row.names=F,quote=F)
# write.table(vartable, paste0(workdir,"comp_h2/complement_covar2.csv"),
#             sep=",", col.names=T, row.names=F, quote=F)
# write.table(vartable, paste0(workdir,"comp_h2/complement_covar3.csv"),
#             sep=",", col.names=T, row.names=F, quote=F)
write.table(vartable, paste0(workdir,"comp_h2/complement_covar4.csv"),
            sep=",", col.names=T, row.names=F, quote=F)

################################ finnish twins
# load data
# add updated dx info
load("/data/swe_gwas/ABZ/NDE1/finlandclean.rdata")
colnames(finlandclean)[2:18567] = colnames(finlandclean)[1:18566]
finlandclean <- finlandclean[,2:18567]

FID <- read.table(paste0(workdir,"comp_h2/missingFID_finns.txt"),header=T)
FID <- FID[FID$IID %in% finlandclean$Persnb,]
finlandclean[57:73,2] <- FID[,2]

# tvab
famID <- unique(finlandclean$Pairnb)
twin1 <- finlandclean[!duplicated(finlandclean$Pairnb),]
twin2 <- finlandclean[!(finlandclean$Persnb %in% twin1$Persnb),]
twin1$tvab <- 1
twin2$tvab <- 2
df2 <- rbind(twin1,twin2)

# input lists
compgenes <- c('ILMN_1662523', "ILMN_1668996", "ILMN_1670305", "ILMN_1695588", 
               "ILMN_1737242", "ILMN_2307740", "ILMN_2360307")
temp <- df2[,colnames(df2) %in% compgenes] 
df2 <- cbind(df2[,c(1:3,7,18567)],temp)

# MZ vs. DZ
df2 <- df2[order(df2$Zygosity),]
MZ <- df2[1:32,]
MZ <- MZ[order(MZ$tvab),]
DZ <- df2[33:73,]
DZ <- DZ[order(DZ$tvab),]

# reshape to wide format
# 16 unique MZ FID; 16 complete pairs
# 21 unique DZ FID; 20 compelte pairs
MZ.wide <- reshape(MZ,idvar="Pairnb",direction="wide",timevar="tvab")
DZ.wide <- reshape(DZ,idvar="Pairnb",direction="wide",timevar="tvab")

# no NAs (complete pairs) - also with complete pheno info
MZcomp <- na.omit(MZ.wide)
DZcomp <- na.omit(DZ.wide)

# write var and covar table
# 5 columns of info
# 7 rows of genes
vartable <- matrix(nrow=length(compgenes),ncol=10)
for (i in 1:length(compgenes)){
  t1gene <- paste0(compgenes[i],".1")
  t2gene <- paste0(compgenes[i],".2")
  vartable[i,1] <- var(MZcomp[,t1gene])
  vartable[i,2] <- cov(MZcomp[,t1gene],MZcomp[,t2gene],
                       use="pairwise.complete.obs")
  vartable[i,3] <- var(MZcomp[,t2gene])
  vartable[i,4] <- var(DZcomp[,t1gene])
  vartable[i,5] <- cov(DZcomp[,t1gene],DZcomp[,t2gene],
                       use="pairwise.complete.obs")
  vartable[i,6] <- var(DZcomp[,t2gene])
  vartable[i,7] <- length(rownames(MZcomp))
  vartable[i,8] <- length(rownames(DZcomp))
  vartable[i,9] <- paste0("X1",compgenes[i])
  vartable[i,10] <- paste0("X2",compgenes[i])
}

vartable <- as.data.frame(vartable)
vartable[,1:6] <- lapply(vartable[,1:6],as.character)
vartable[,1:6] <- lapply(vartable[,1:6],as.numeric)
names(vartable) <- c("MZtwin1var","MZcovar","MZtwin2var",
                     "DZtwin1var","DZcovar","DZtwin2var",
                     "NUM_PAIRS_MZ","NUM_PAIRS_DZ","T1_VAR","T2_VAR")

# zero out negative covariance
vartable$MZcovar <- ifelse(vartable$MZcovar < 0,0,vartable$MZcovar)
vartable$DZcovar <- ifelse(vartable$DZcovar < 0,0,vartable$DZcovar)

# write csv file
write.table(vartable,paste0(workdir,"/comp_h2/complement_covar_finns.csv"),
            sep=",",col.names=T,row.names=F,quote=F)

################################ compare mean expression swe v finn
df2 <- df2[,c(1:2,5,3,6:12)]
df$study <- "swe"
df2$study <- "finn"
names(df2) <- names(df)
df2$zyg <- as.factor(df2$zyg) 
levels(df2$zyg) <- c("MZ","DZ",NA)
all <- rbind(df,df2)
all$Family <- as.factor(all$Family)

describeBy(all[,5:11],all$sample)
study <- matrix(nrow=length(compgenes),ncol=4)
for (i in 1:length(compgenes)){
  f <- formula(paste(compgenes[i], "study", sep="~"))
  m1 <- lme(fixed=f, random = ~1|Family, data=all, na.action=na.omit)
  tval <- summary(m1)$tTable[2,4]
  study[i,1] <- compgenes[i]
  study[i,2] <- tval
  study[i,3] <- sqrt(tval^2/(tval^2+summary(m1)$tTable[2,3]))
  study[i,4] <- summary(m1)$tTable[2,5]
}

study <- as.data.frame(study)
names(study) <- c("gene","t","corr","p")

################################ combined h2 - swe + fin
# load data - swedish
# add updated dx info
load("/data/swe_gwas/ABZ/RNA_GWAS/swedenclean.rdata")
dx = read.table("/data/swe_gwas/ABZ/RNA_GWAS/swedishNP_PCA.txt",header=T)
swedenclean$dx_new <- dx[match(swedenclean$StudyID,dx$IID),7] %>% as.factor()
controls <- swedenclean[swedenclean$dx_new==7,1]

# load data - finnish
# add updated dx info
load("/data/swe_gwas/ABZ/NDE1/finlandclean.rdata")
colnames(finlandclean)[2:18567] = colnames(finlandclean)[1:18566]
finlandclean <- finlandclean[,2:18567]

FID <- read.table(paste0(workdir,"comp_h2/missingFID_finns.txt"),header=T)
FID <- FID[FID$IID %in% finlandclean$Persnb,]
finlandclean[57:73,2] <- FID[,2]

# zygosity + tvab
zyg <- read.table(paste0(workdir,"zygosity_alltwins.txt"),header=T)
swedenclean$tvab <- zyg[match(swedenclean$StudyID,zyg$IID),4]
swedenclean$zyg <- zyg[match(swedenclean$StudyID,zyg$IID),3]

famID <- unique(finlandclean$Pairnb)
twin1 <- finlandclean[!duplicated(finlandclean$Pairnb),]
twin2 <- finlandclean[!(finlandclean$Persnb %in% twin1$Persnb),]
twin1$tvab <- 1
twin2$tvab <- 2
df2 <- rbind(twin1,twin2)

# input lists
compgenes <- c('ILMN_1662523', "ILMN_1668996", "ILMN_1670305", "ILMN_1695588", 
               "ILMN_1737242", "ILMN_2307740", "ILMN_2360307")

swedecompgenes <- swedenclean[,colnames(swedenclean) %in% compgenes] 
df <- cbind(swedenclean[,c(1,2,18567:18568)],swedecompgenes)
names(df)[1:4] <- c("IID","FID","tvab","zyg")

temp <- df2[,colnames(df2) %in% compgenes] 
df2 <- cbind(df2[,c(1,2,18567,3)],temp)
names(df2)[1:4] <- c("IID","FID","tvab","zyg")
df2$zyg <- factor(df2$zyg, levels = c(2,1,9)) 
levels(df2$zyg) <- c("DZ","MZ",NA) # 1 = MZ; 2 = DZ

dfall <- rbind(df,df2)
dfall <- dfall[!is.na(dfall$zyg),]

# z-score to swedish controls
means <- apply(dfall[dfall$IID %in% controls,5:11],2,mean)
sds <- apply(dfall[dfall$IID %in% controls,5:11],2,sd)
zs <- (dfall[,5:11] - means)/sds 
names(zs) <- for (i in 1:7){paste(compgenes[i],i,sep=".")}
dfall <- cbind(dfall,zs)

# MZ vs. DZ
MZ    <- dfall[dfall$zyg=="MZ",]
DZ    <- dfall[dfall$zyg=="DZ",]

# reshape to wide format
# 59 unique MZ FID; 50 complete pairs
# 85 unique DZ FID; 69 compelte pairs
MZ.wide <- reshape(MZ,idvar="FID",direction="wide",timevar="tvab")
DZ.wide <- reshape(DZ,idvar="FID",direction="wide",timevar="tvab")

# no NAs (complete pairs) - also with complete pheno info
MZcomp <- na.omit(MZ.wide)
DZcomp <- na.omit(DZ.wide)

# write var and covar table
# 5 columns of info
# 7 rows of genes
vartable <- matrix(nrow=length(compgenes),ncol=10)
for (i in 1:length(compgenes)){
  t1gene <- paste0(compgenes[i],".1")
  t2gene <- paste0(compgenes[i],".2")
  vartable[i,1] <- var(MZcomp[,t1gene])
  vartable[i,2] <- cov(MZcomp[,t1gene],MZcomp[,t2gene],
                       use="pairwise.complete.obs")
  vartable[i,3] <- var(MZcomp[,t2gene])
  vartable[i,4] <- var(DZcomp[,t1gene])
  vartable[i,5] <- cov(DZcomp[,t1gene],DZcomp[,t2gene],
                       use="pairwise.complete.obs")
  vartable[i,6] <- var(DZcomp[,t2gene])
  vartable[i,7] <- length(rownames(MZcomp))
  vartable[i,8] <- length(rownames(DZcomp))
  vartable[i,9] <- paste0("X1",compgenes[i])
  vartable[i,10] <- paste0("X2",compgenes[i])
}

vartable <- as.data.frame(vartable)
vartable[,1:6] <- lapply(vartable[,1:6],as.character)
vartable[,1:6] <- lapply(vartable[,1:6],as.numeric)
names(vartable) <- c("MZtwin1var","MZcovar","MZtwin2var",
                     "DZtwin1var","DZcovar","DZtwin2var",
                     "NUM_PAIRS_MZ","NUM_PAIRS_DZ","T1_VAR","T2_VAR")

# zero out negative covariance
# vartable$MZcovar <- ifelse(vartable$MZcovar < 0,0,vartable$MZcovar)
# vartable$DZcovar <- ifelse(vartable$DZcovar < 0,0,vartable$DZcovar)

# write csv file
write.table(vartable,paste0(workdir,"comp_h2/complement_covar_swefin.csv"),
            sep=",",col.names=T,row.names=F,quote=F)

