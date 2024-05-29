### version-1 codes----
load("1.LUAD.mRNA.data.related.Rdata")
load("2.LUAD.Cli.related.Rdata")
load("3.LUAD.dat1.Rdata")

svdata1 <- cbind(meta1[,c("time","event")],t(mRNA.data.tumor.wt2.FerrDb1))

svdata2 <- svdata1
for(i in colnames(svdata2)[3:ncol(svdata2)]){svdata2[,i] <- ifelse(svdata2[,i]> median(svdata2[,i]),1,0)}

#二分变量差异表达基因的单因素Cox
library("survival")
library("survminer")


covariates2 <- colnames(svdata2[,3:ncol(svdata2)])
univ_formulas2 <- sapply(covariates2,
                         function(x) as.formula(paste('Surv(time,event)~', "\`",x,"\`",sep="")))
univ_models2 <- lapply(univ_formulas2, function(x){coxph(x, data = svdata2)})

univ_results2 <- lapply(univ_models2,
                        function(x){ 
                          x <- summary(x)
                          p.value<-signif(x$wald["pvalue"], digits=2)
                          wald.test<-signif(x$wald["test"], digits=2)
                          beta<-signif(x$coef[1], digits=2);#coeficient beta
                          HR <-signif(x$coef[2], digits=2);#exp(beta)
                          HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                          HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                          HR <- paste0(HR, " (", 
                                       HR.confint.lower, "-", HR.confint.upper, ")")
                          res2<-c(beta, HR, wald.test, p.value)
                          names(res2)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                         "p.value")
                          return(res2)
                          #return(exp(cbind(coef(x),confint(x))))
                        })
res2 <- t(as.data.frame(univ_results2, check.names = FALSE))
OS.res2 <- as.data.frame(res2)

single_cox2 <- as.data.frame(OS.res2)
sig_gene_multi_cox2 <- single_cox2[single_cox2$p.value<=0.05,]
sig_gene_multi_cox2$gene <- rownames(sig_gene_multi_cox2)
unicox.res.sig2 <- sig_gene_multi_cox2

multi.ferr.deg2 <- svdata1[,rownames(unicox.res.sig2)]
multi.ferr.deg2 <- cbind(svdata2[,1:2],multi.ferr.deg2) 


mySurv2 <-with(multi.ferr.deg2[,1:2],Surv(time,event))
multi_COX2 <-coxph(mySurv2 ~ ., data=multi.ferr.deg2[,3:ncol(multi.ferr.deg2)])  
x2<-summary(multi_COX2)
x2
ggforest(multi_COX2,data=multi.ferr.deg2)

step.multi_Cox2 <- step(multi_COX2,direction = "both")

y2 <- summary(step.multi_Cox2)
y2
y3 <- data.frame(y2)
ggforest(step.multi_Cox2,data=multi.ferr.deg2)

library(survivalROC)

riskScore2=predict(step.multi_Cox2, type="risk", newdata=multi.ferr.deg2)



CoxGenes2=rownames(as.data.frame(y2$coefficients))
CoxGene2=gsub("`", "", CoxGenes2)


save(step.multi_Cox2,riskScore2,CoxGenes2,file="./step.multi_Cox.Rdata")



outCol2=c("time", "event", CoxGenes2)
riskOut2=cbind(multi.ferr.deg2[,outCol2], riskScore2)
riskOut2=cbind(id=rownames(riskOut2), riskOut2)
#write.table(riskOut2, file="step.riskScore2.txt", sep="\t", quote=F, row.names=F)


outTab2=data.frame()
outTab2=cbind(
  coef=y2$coefficients[,"coef"],
  HR=y2$conf.int[,"exp(coef)"],
  HR.95L=y2$conf.int[,"lower .95"],
  HR.95H=y2$conf.int[,"upper .95"],
  pvalue=as.numeric(y2$coefficients[,"Pr(>|z|)"]))
outTab2=cbind(id=row.names(outTab2),outTab2)
outTab2=gsub("`","",outTab2)
#write.table(outTab2,file="step.multi.Cox2.txt",sep="\t",row.names=F,quote=F)


new_data2 <- svdata2[,1:2]
new_data2$fp <- riskScore2
library(timeROC)

result2_1 <-with(new_data2, timeROC(T=time,
                                    delta=event,
                                    marker=fp,
                                    cause=1,
                                    times=c(12,36,60),
                                    iid = TRUE))
dat2 = data.frame(fpr = as.numeric(result2_1$FP),
                  tpr = as.numeric(result2_1$TP),
                  time = rep(as.factor(c(12,36,60)),each = nrow(result2_1$TP)))

##c("#92C5DE", "#F4A582", "#66C2A5")

library(ggplot2)
ggplot() + 
  geom_line(data = dat2,aes(x = fpr, y = tpr,color = time),size = 1) + 
  scale_color_manual(name = NULL,values = c("#CC6633", "#9966CC", "#66CC00"),
                     labels = paste0("AUC of ",c(1,3,5),"-y survival: ",
                                     format(round(result2_1$AUC,2),nsmall = 2)))+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
        legend.position = c(0.765,0.125))+
  scale_x_continuous(expand = c(0.005,0.005))+
  scale_y_continuous(expand = c(0.005,0.005))+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  coord_fixed()

#用survivalROC寻找最佳截断点

library(survivalROC)
result2_2 <-with(new_data2, 
                 survivalROC(Stime=time,
                             status=event,
                             marker=fp,
                             predict.time=12,
                             method="KM"))
res.cut2<- result2_2$cut.values[which.max(result2_2$TP-result2_2$FP)]
new_data2$risk = ifelse(new_data2$fp< res.cut2,"low","high")
fit2 <- survfit(Surv(time, event)~risk, data=new_data2)
ggsurvplot(fit2, data=new_data2, pval=TRUE,palette=c('red','green'),
           risk.table = TRUE, conf.int = F) 
#ggsurvplot(fit2, data=new_data2, pval=TRUE,palette=c('#CC6633','#66CC00'),risk.table = TRUE, conf.int = TRUE) 


#开始绘制风险模型的生存点图
fp <- riskScore2
phe <- new_data2
fp_dat=data.frame(patientid=1:length(fp),fp=as.numeric(sort(fp)))
##添加风险分组，以风险评分的cutpoint值将患者分为两组，大于cutpoint的患者为高风险组，小于或等于cutpoint的患者为低风险组
fp_dat$riskgroup= ifelse(fp_dat$fp < res.cut2,"low","high")


sur_dat=data.frame(patientid=1:length(fp),time=phe[names(sort(fp)),'time'],event=phe[names(sort(fp)),'event']) 
sur_dat$event=ifelse(sur_dat$event==0,'alive','death')
sur_dat$event=factor(sur_dat$event,levels = c("death","alive"))
exp_dat=svdata1[names(sort(fp)),CoxGene2]
#fp_dat用来绘制第一幅图
#sur_dat用来绘制第二幅图
#exp_dat用来绘制第三幅图

###第一个图
library(ggplot2)
p1=ggplot(fp_dat,aes(x=patientid,y=fp))+geom_point(aes(color=riskgroup))+
  scale_colour_manual(values = c("red","green"))+
  theme_bw()+labs(x="Patient ID(increasing risk score)",y="Risk score")+
  geom_hline(yintercept=median(fp_dat$fp),colour="black", linetype="dotted",size=0.8)+
  geom_vline(xintercept=sum(fp_dat$riskgroup=="low"),colour="black", linetype="dotted",size=0.8)
p1

#第二个图
p2=ggplot(sur_dat,aes(x=patientid,y=time))+geom_point(aes(col=event))+theme_bw()+
  scale_colour_manual(values = c("red","green"))+
  labs(x="Patient ID(increasing risk score)",y="Survival time(year)")+
  geom_vline(xintercept=sum(fp_dat$riskgroup=="low"),colour="black", linetype="dotted",size=0.8)
p2

#第三个图
library(pheatmap)
mycolors <- colorRampPalette(c("white", "green", "red"), bias = 1.2)(100)
tmp=t(scale(exp_dat))
tmp[tmp > 1] = 1
tmp[tmp < -1] = -1


annotation_col = data.frame(fp_dat$fp)
colnames(annotation_col) <- "riskScore"
rownames(annotation_col) <- rownames(exp_dat)
p3=pheatmap(tmp,annotation_col = annotation_col,col= mycolors,show_colnames = F,cluster_cols = F)
p3

#拼图实现三图联动
library(ggplotify)
plots = list(p1,p2,as.ggplot(as.grob(p3)))
library(gridExtra)
lay1 = rbind(c(rep(1,7)),c(rep(2,7)),c(rep(3,7))) #布局矩阵
grid.arrange(grobs = plots, layout_matrix = lay1, heigths = c(2,3,2),weights=c(10,10,10))


library(ggplotify)
plots2 = list(p1,p2)
library(gridExtra)
lay2 = rbind(c(rep(1,7)),c(rep(2,7))) #布局矩阵
grid.arrange(grobs = plots2, layout_matrix = lay2, heigths = c(2,3),weights=c(10,10))

library(survival)
library(rms)
dd <- datadist(meta1)
options(datadist="dd")
f <- cph(Surv(time, event)~ age1 + gender_group+riskScore, data = meta1, 
         x = T, y = T, surv = T)
surv <- Survival(f)
s1 <- function(x)surv(1*12, x)
s3 <- function(x)surv(3*12, x)
s5 <- function(x)surv(5*12, x)
s10 <- function(x)surv(120, x)
s20 <- function(x)surv(241.6, x)

nom <- nomogram(f, fun = list(s1, s3, s5), fun.at = c(0.05, seq(0.1, 0.9, by=0.05), 0.95), 
                funlabel = c("1-year survival", '3-year survival', '5-year survival'))
plot(nom)


f2 <- psm(Surv(time, event)~ riskScore, data = meta1, x=T, y=T, dist='lognormal')

cal1 <- calibrate(f2, 
                  cmethod='KM', 
                  method="boot", 
                  u=24, # u需要与之前模型中定义好的time.inc一致，即365或730；
                  m=60, #每次抽样的样本量，
                  B=1000)
plot(cal1)

meta1$patient <- rownames(meta1)
annotation_col$patient <- rownames(annotation_col)
meta1 <- merge.data.frame(meta1, annotation_col, by = "patient")



library(stringi)
library(GOplot)

GO <- read.table(file = "./GO.txt", header = T, sep = "\t")



normal.sample <- substr(colnames(mRNA.data.normal), 1, 12)
tumor.sample <- substr(colnames(mRNA.data.tumor.wt2.FerrDb1), 1, 12)

paird.tumor <- mRNA.data.tumor.wt2.FerrDb1[-120, which(tumor.sample %in% normal.sample)]
paird.normal <- mRNA.data.normal[rownames(paird.tumor), which(normal.sample %in% tumor.sample)]

ddit4 <- data.frame(normal = as.numeric(paird.normal['DDIT4', ]), tumor = as.numeric(paird.tumor['DDIT4', ]))
rownames(ddit4) <- colnames()

ddit4.plot <- gather(ddit4, key = "key", value = "value")
ggplot(ddit4.plot, aes(x= key, y=value, color = key))+
  geom_boxplot()+
  geom_jitter(width = 0.25)+
  stat_compare_means(method = "t.test",)+
  theme_bw()+
  labs(x = "", y = "log2(FPKM+1)")+
  theme(axis.title.y = element_text(size = 20), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20))


mRNA.data.tumor.wt2.FerrDb1 <- mRNA.data.tumor.wt2.FerrDb1[-120, ]
mRNA.data.normal.FerrDb1 <- mRNA.data.normal[rownames(mRNA.data.tumor.wt2.FerrDb1), ]

forDEG <- cbind.data.frame(mRNA.data.normal.FerrDb1, mRNA.data.tumor.wt2.FerrDb1)
forDEG <- 2^forDEG-1
forDEG <- round(forDEG)

condition <- factor(c(rep("normal",58),rep("tumor",286)))
colData <- data.frame(row.names=colnames(forDEG), condition)
dds <- DESeqDataSetFromMatrix(countData = forDEG, colData = colData, design = ~ condition)
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE) 
res <- results(dds1)
# res格式转化：用data.frame转化为表格形式
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
# 依次按照pvalue值log2FoldChange值进行排序
res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
res1$gene <- rownames(res1)
res1$change <- ifelse(res1$log2FoldChange > 1 & res1$padj < 0.05, "up", 
                      ifelse(res1$log2FoldChange < -1 & res1$padj < 0.05, "down", "none"))
res1$label <- ifelse(res1$gene %in% c('DDIT4', "NOX4", "GCLC"), res1$gene, "")

library(tidyverse)
UP.gene <- res1 %>% filter(change == "up") %>% pull(gene)
risk.gene <- c('DDIT4', "NOX4", "DUOX1", "PRKCA", "GCLC", "PARK7")
intersect(UP.gene, risk.gene) #"GCLC"  "NOX4"  "DDIT4"

ggplot(res1, aes(x=log2FoldChange, y=-log10(padj), color=change))+
  geom_point(alpha=0.4, size=3.5)+
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_hline(yintercept = 5, lty=2)+
  geom_vline(xintercept = c(-1, 1), lty=2)


p+geom_label_repel(data = res1, aes(x=log2FoldChange, y=-log10(padj), label=label), 
                   size = 4, color = "black", 
                   box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.8, "lines"), 
                   segment.color = "black",
                   segment.size = 0.4)

library(GSVA)
risklist <- list(risk = c("DDIT4", "NOX4", "GCLC"))

riskscore.ssgsea1 <- gsva(as.matrix(mRNA.data.tumor.wt2.FerrDb1), risklist, method="gsva", verbose=FALSE, kcdf = "Gaussian")
riskscore.ssgsea1 <- as.data.frame(t(riskscore.ssgsea1))

riskscore_surv <- cbind.data.frame(meta1, riskscore.ssgsea1)
riskscore_surv$patient <- rownames(riskscore_surv)  
riskscore_surv$TCGA.Participant.Barcode <- substr(riskscore_surv$patient, 1, 12)


library(timeROC)
library(survival)

time_roc_res <- timeROC(
  T = riskscore_surv$time,
  delta = riskscore_surv$event,
  marker = riskscore_surv$risk,
  cause = 1,
  weighting="marginal",
  times = c(1 * 12, 3 * 12, 5 * 12, 10*12),
  ROC = TRUE,
  iid = TRUE
)
time_roc_res$AUC

plot(time_roc_res, time=1 * 12, col = "red", title = FALSE)  
plot(time_roc_res, time=3 * 12, add=TRUE, col="blue") 
plot(time_roc_res, time=5 * 12, add=TRUE, col="green") 
plot(time_roc_res, time=10 * 12, add=TRUE, col="orange") 

legend("bottomright",c("1 Years" ,"3 Years", "5 Years", "10 Years"),
       col=c("red", "blue", "green", "orange"), lty=1, lwd=2)

dat2 = data.frame(fpr = as.numeric(time_roc_res$FP),
                  tpr = as.numeric(time_roc_res$TP),
                  time = rep(as.factor(c(1 * 12, 3 * 12, 5 * 12, 10*12)),each = nrow(time_roc_res$TP)))
ggplot() + 
  geom_line(data = dat2,aes(x = fpr, y = tpr,color = time),size = 1) + 
  scale_color_manual(name = NULL,values = c("#CC6633", "#9966CC", "#66CC00", "orange"),
                     labels = paste0("AUC of ",c(1,3,5,10),"-y survival: ",
                                     format(round(time_roc_res$AUC,2),nsmall = 2)))+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
        legend.position = c(0.765,0.125))+
  scale_x_continuous(expand = c(0.005,0.005))+
  scale_y_continuous(expand = c(0.005,0.005))+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  coord_fixed()


library(rms)
dd=datadist(riskscore_surv)
options(datadist="dd")

f <- cph(Surv(time, event)~ age1+gender_group+risk+stage_group, data = riskscore_surv, 
         x = T, y = T, surv = T)
surv <- Survival(f)
s1 <- function(x)surv(1*12, x)
s3 <- function(x)surv(3*12, x)
s5 <- function(x)surv(5*12, x)
s10 <- function(x)surv(10*12, x)

nom <- nomogram(f, fun = list(s1, s3, s5, s10), fun.at = c(0.05, seq(0.1, 0.9, by=0.2), 0.95), 
                funlabel = c("1-year survival", "3-year survival", '5-year survival', '10-year survival'))
plot(nom, xfrac = 0.45)

f2 <- psm(Surv(time, event)~ age1+gender_group+risk+stage_group, data = riskscore_surv, x=T, y=T, dist='lognormal')

cal1 <- calibrate(f2, 
                  cmethod='KM', 
                  method="boot", 
                  u=1*12, # u需要与之前模型中定义好的time.inc一致，即365或730；
                  m=57, #每次抽样的样本量，
                  B=1000)
plot(cal1)

res.cut <- surv_cutpoint(riskscore_surv, #数据集
                         time = "time", #生存状态
                         event = "event", #生存时间
                         variables = "risk" #需要计算的数据列名
)
res.cat <- surv_categorize(res.cut)
fit1 <- survfit(Surv(time, event)~risk, data=res.cat)
ggsurvplot(fit1, data=res.cat, pval=TRUE,palette = "jco",
           risk.table = TRUE, conf.int = F) 


HR <- sig_gene_multi_cox2
HR.split <- strsplit(HR$`HR (95% CI for HR)`, split = "\\s")

for (i in 1:nrow(HR)) {
  HR$HR.value[i] <-HR.split[[i]][1] 
  HR$CI[i] <- HR.split[[i]][2]
}

HR$CI <- gsub(pattern = "\\(", replacement = "", x = HR$CI)
HR$CI <- gsub(pattern = "\\)", replacement = "", x = HR$CI)

CI <- strsplit(HR$CI, split = "-")

for (i in 1:nrow(HR)) {
  HR$CI.L[i] <- CI[[i]][1] 
  HR$CI.H[i] <- CI[[i]][2]
}

HR$HR.value <- as.numeric(HR$HR.value)
HR$CI.L <- as.numeric(HR$CI.L)
HR$CI.H <- as.numeric(HR$CI.H)

HR$pvalue <- as.numeric(HR$p.value)
HR$index <- c(1:nrow(HR))
HR1 <- HR[order(HR$HR.value), ]
HR1$index <- rownames(HR1)
HR1$index <- factor(HR1$index, levels = c(HR1$index))

library(ggplot2)
ggplot(HR1, aes(x=HR.value, y=index))+
  geom_point()+
  annotate("segment", x = 1.20, xend = 2.50, y=1, yend = 1, size=1,
           arrow=arrow(ends = "both", angle = 90, length = unit(.2, "cm")))

p <- ggplot(HR1, aes(x=HR.value, y=index))+
  geom_point(shape=22, size=4.2, fill="black")
p
for (i in 1:31) {
  p <- p+annotate("segment", x = HR1$CI.L[i], xend = HR1$CI.H[i], y=i, yend = i, size=1,
                  arrow=arrow(ends = "both", angle = 90, length = unit(.1, "cm")))
}
p

p1 <- p+geom_vline(aes(xintercept=1), linetype="dashed")
p1

for (i in 1:31) {
  p1 <- p1 + annotate("text", x=-0.3, y=i, label=HR1$pvalue[i], color="black")
}
p1

p2 <- p1 + theme_bw()+xlim(c(-1.8, 2.7))
p2
for (i in 1:31) {
  p2 <- p2 + annotate("text", x=-1.5, y=i, label=HR1$HR.value[i], color="black")
}
p2

for (i in 1:31) {
  p2 <- p2 + annotate("text", x=-0.8, y=i, label=HR1$CI[i], color = "black")
}
p2

p2+theme_bw()+theme(panel.grid = element_blank())+
  labs(x = "HR", y = "")+
  theme(axis.title.y = element_text(size = 10, color = "black"), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"))


HR.multi <- as.data.frame(cbind(y2$coefficients, y2$conf.int))
colnames(HR.multi) <- c("coef", "HR", "HR.se", "z", "p.value", "HR1", "HR2", "CI.L", "CI.H")

HR.multi1 <- round(HR.multi, 2)
HR.multi1$CI <- paste(HR.multi1$CI.L, HR.multi1$CI.H, sep = "-")
HR.multi1$index <- c(1:nrow(HR.multi1))
HR.multi1 <- HR.multi1[order(HR.multi1$HR), ]
HR.multi1$index <- rownames(HR.multi1)
HR.multi1$index <- factor(HR.multi1$index, levels = c(HR.multi1$index))

ggplot(HR.multi1, aes(x=HR, y=index))+
  geom_point()+
  annotate("segment", x = 1.20, xend = 2.50, y=1, yend = 1, size=1,
           arrow=arrow(ends = "both", angle = 90, length = unit(.2, "cm")))

p <- ggplot(HR.multi1, aes(x=HR, y=index))+
  geom_point(shape=22, size=4.2, fill="black")
p
for (i in 1:13) {
  p <- p+annotate("segment", x = HR.multi1$CI.L[i], xend = HR.multi1$CI.H[i], y=i, yend = i, size=1,
                  arrow=arrow(ends = "both", angle = 90, length = unit(.1, "cm")))
}
p

p1 <- p+geom_vline(aes(xintercept=1), linetype="dashed")
p1

for (i in 1:13) {
  p1 <- p1 + annotate("text", x=-0.3, y=i, label=HR.multi1$p.value[i], color="black")
}
p1

p2 <- p1 + theme_bw()+xlim(c(-2, 5))
p2
for (i in 1:13) {
  p2 <- p2 + annotate("text", x=-2, y=i, label=HR.multi1$HR[i], color="black")
}
p2

for (i in 1:13) {
  p2 <- p2 + annotate("text", x=-1.3, y=i, label=HR.multi1$CI[i], color = "black")
}
p2

p2+theme_bw()+theme(panel.grid = element_blank())+
  labs(x = "HR", y = "")+
  theme(axis.title.y = element_text(size = 10, color = "black"), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"))

genes <- as.character(rev(c(HR.multi1$index)))

res.genes <- res1[genes, ]
res.genes$x <- "1"
res.genes$gene <- factor(res.genes$gene, levels = rev(c("PARK7","PRKCA","DUOX1","NOX4","DDIT4","GCLC","ARRDC3","ACSF2","SNCA","GDF15","PEBP1","HERPUD1","MFN2")))
ggplot(data = res.genes, aes(x=x, y=gene, color=log2FoldChange, size=-log10(padj)))+
  geom_point()+
  scale_color_gradient2(low = "cyan", high = "red", midpoint = 0)+
  theme_bw()+
  theme(panel.grid.major.y = element_blank())+
  theme(axis.title.y = element_text(size = 10, color = "black"), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"))

LUAD.FPKM <- read_tsv(file = "TCGA-LUAD.htseq_fpkm.tsv.gz")
gene.table <- read.table(file = "gencode.v22.annotation.gene.probeMap", header = T, sep = "\t")
gene.table <- gene.table[, c(1, 2)]
colnames(gene.table)[1] <- c("Ensembl_ID")

LUAD.FPKM <- merge.data.frame(LUAD.FPKM, gene.table, by = "Ensembl_ID")
LUAD.FPKM <- as.data.frame(avereps(LUAD.FPKM[, -c(1, 587)], ID = LUAD.FPKM$gene))

normal.wt.plot <- LUAD.FPKM[genes, colnames(mRNA.data.normal.FerrDb1)]
tumor.wt.plot <- LUAD.FPKM[genes, colnames(mRNA.data.tumor.wt2.FerrDb1)]

ntplot <- cbind.data.frame(normal.wt.plot, tumor.wt.plot)
library(pheatmap)
anno <- data.frame(anno=c(rep("normal", 35), rep("tumor", 35)))
rownames(anno) <- colnames(paired.plot)

ntplot1 <- as.matrix(scale(ntplot))
ntplot2 <- ifelse(ntplot1 < -2, -2, 
                  ifelse(ntplot1 > 0.2, 2, ntplot1))
paired.plot <- cbind.data.frame(paird.normal, paird.tumor)
paired.plot <- paired.plot[genes, ]

pheatmap(scale(paired.plot), annotation_col = anno, cluster_rows = F, cluster_cols = F, show_colnames = F, scale = "row", border = F)



ddit4.n <- data.frame(DDIT4 = as.numeric(paird.normal['DDIT4', ]), sample = colnames(paird.normal), group = "Normal") %>% arrange(sample) %>% transform(index = 1:nrow(ddit4))
ddit4.t <- data.frame(DDIT4 = as.numeric(paird.tumor['DDIT4', ]), sample = colnames(paird.tumor), group = "Tumor")%>% arrange(sample) %>% transform(index = 1:nrow(ddit4))
ddit4 <- rbind.data.frame(ddit4.n, ddit4.t)
ddit4$group <- factor(ddit4$group, levels = c('Tumor', "Normal"))

ggplot(ddit4, aes(x= group, y=DDIT4, color = group))+
  geom_boxplot(aes(color = group))+
  geom_point(size = 3,shape = 1)+
  geom_line(aes(group = index), color = 'grey', lwd = 0.14)+
  stat_compare_means(method = "t.test",paired = TRUE, comparisons=list(c("Normal", "Tumor")))+
  theme_bw()+
  theme(panel.grid = element_line(linetype = "dashed"))+
  labs(x="", y="Log2(DDIT4+1)")+
  theme(axis.title.y = element_text(size = 14, color = "black"), 
        axis.text.x = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(size = 14, color = "black"),
        legend.position = "none")

NOX4.n <- data.frame(NOX4 = as.numeric(paird.normal['NOX4', ]), sample = colnames(paird.normal), group = "Normal") %>% arrange(sample) %>% transform(index = 1:35)
NOX4.t <- data.frame(NOX4 = as.numeric(paird.tumor['NOX4', ]), sample = colnames(paird.tumor), group = "Tumor")%>% arrange(sample) %>% transform(index = 1:35)
NOX4 <- rbind.data.frame(NOX4.n, NOX4.t)
NOX4$group <- factor(NOX4$group, levels = c('Tumor', "Normal"))

ggplot(NOX4, aes(x= group, y=NOX4, color = group))+
  geom_boxplot(aes(color = group))+
  geom_point(size = 3,shape = 1)+
  geom_line(aes(group = index), color = 'grey', lwd = 0.14)+
  stat_compare_means(method = "t.test",paired = TRUE, comparisons=list(c("Normal", "Tumor")))+
  theme_bw()+
  theme(panel.grid = element_line(linetype = "dashed"))+
  labs(x="", y="Log2(NOX4+1)")+
  theme(axis.title.y = element_text(size = 14, color = "black"), 
        axis.text.x = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(size = 14, color = "black"),
        legend.position = "none")

GCLC.n <- data.frame(GCLC = as.numeric(paird.normal['GCLC', ]), sample = colnames(paird.normal), group = "Normal") %>% arrange(sample) %>% transform(index = 1:35)
GCLC.t <- data.frame(GCLC = as.numeric(paird.tumor['GCLC', ]), sample = colnames(paird.tumor), group = "Tumor")%>% arrange(sample) %>% transform(index = 1:35)
GCLC <- rbind.data.frame(GCLC.n, GCLC.t)
GCLC$group <- factor(GCLC$group, levels = c('Tumor', "Normal"))

ggplot(GCLC, aes(x= group, y=GCLC, color = group))+
  geom_boxplot(aes(color = group))+
  geom_point(size = 3,shape = 1)+
  geom_line(aes(group = index), color = 'grey', lwd = 0.14)+
  stat_compare_means(method = "t.test",paired = TRUE, comparisons=list(c("Normal", "Tumor")))+
  theme_bw()+
  theme(panel.grid = element_line(linetype = "dashed"))+
  labs(x="", y="Log2(GCLC+1)")+
  theme(axis.title.y = element_text(size = 14, color = "black"), 
        axis.text.x = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(size = 14, color = "black"),
        legend.position = "none")


GENE_3_CLI <- t(mRNA.data.tumor.wt2.FerrDb1) %>% as.data.frame() %>% select(DDIT4, NOX4, GCLC)
GENE_3_CLI <- cbind.data.frame(meta1, GENE_3_CLI)


res.cut <- surv_cutpoint(GENE_3_CLI, #数据集
                         time = "time", #生存状态
                         event = "event", #生存时间
                         variables = c("DDIT4", "NOX4", "GCLC") #需要计算的数据列名
)

res.cat <- surv_categorize(res.cut)
fit1 <- survfit(Surv(time, event)~GCLC, data=res.cat)
ggsurvplot(fit1, data=res.cat, pval=TRUE,palette = "jco",
           risk.table = TRUE, conf.int = F, title ="GCLC") 


res.cut <- surv_cutpoint(riskscore_surv, #数据集
                         time = "time", #生存状态
                         event = "event", #生存时间
                         variables = "risk" #需要计算的数据列名
)



riskscore_surv$riskgroup2 <- ifelse(riskscore_surv$risk <= median(riskscore_surv$risk), "low", "high")


res.cat <- surv_categorize(res.cut)
fit1 <- survfit(Surv(time, event)~risk, data=res.cat)
ggsurvplot(fit1, data=res.cat, pval=TRUE,palette = "jco",
           risk.table = TRUE, conf.int = F) 

fit12 <- survfit(Surv(time, event)~riskgroup2, data=riskscore_surv)
ggsurvplot(fit12, data=riskscore_surv, pval=TRUE,palette = "jco",
           risk.table = TRUE, conf.int = F) 




riskscore_surv$riskgroup <- res.cat$risk

III.surv <- riskscore_surv %>% filter(stage.short == "III")
res.cutIII <- surv_cutpoint(III.surv, #数据集
                          time = "time", #生存状态
                          event = "event", #生存时间
                          variables = "risk" #需要计算的数据列名
)
res.catIII <- surv_categorize(res.cutIII)

fitIII <- survfit(Surv(time, event)~risk, data=res.catIII)
ggsurvplot(fitIII, data=res.catIII, pval=TRUE,palette = "jco",
           risk.table = TRUE, conf.int = F) 

riskscore_surv$stage.short <- ifelse(riskscore_surv$stage_group %in% c("IB"), "I", 
                                     ifelse(riskscore_surv$stage_group %in% c("II", 'IIA', 'IIB'), "II", "III"))

riskscore_surv$patient <- rownames(riskscore_surv)
###triplots----
#开始绘制风险模型的生存点图
#fp <- score.clinical$score
#phe <- score.clinical
#fp_dat=data.frame(patientid=1:length(fp),fp=as.numeric(sort(fp)))
#添加风险分组，以风险评分的cutpoint值将患者分为两组，大于cutpoint的患者为高风险组，小于或等于cutpoint的患者为低风险组
#fp_dat$riskgroup= ifelse(fp_dat$fp < res.cut1,"low","high")
fp_dat <- riskscore_surv[, c(1,20, 24, 23)]
fp_dat <- fp_dat[order(fp_dat$risk, decreasing = F), ]
fp_dat$patientid <- c(1:nrow(fp_dat))

sur_dat <- riskscore_surv[, c(20,24,1, 13)]
sur_dat=sur_dat[order(sur_dat$risk, decreasing = F), ]
sur_dat$event=factor(sur_dat$vital_status, levels = c("Dead","Alive"))
sur_dat$patientid <- c(1:nrow(sur_dat))

exp_dat=mRNA.data.tumor.wt2.FerrDb1[c("DDIT4", "NOX4", "GCLC"), fp_dat$patient]
#exp_dat <- log10(exp_dat+1)

library(ggplot2)
p1=ggplot(fp_dat,aes(x=patientid,y=risk))+
  geom_point(aes(color=riskgroup2))+
  scale_colour_manual(values = c("#0073C2FF", "#EFC000FF"))+
  theme_bw()+labs(x="Patient ID(increasing risk score)",y="Risk score")+
  geom_hline(yintercept=median(riskscore_surv$risk), colour="black", linetype="dotted",size=0.8)+
  geom_vline(xintercept= 143,colour="black", linetype="dotted",size=0.8)
p1

p2=ggplot(sur_dat,aes(x=patientid,y=time))+
  geom_point(aes(col=event))+theme_bw()+
  scale_colour_manual(values = c("#0073C2FF", "#EFC000FF"))+
  labs(x="Patient ID(increasing risk score)",y="Survival time(Month)")+
  geom_vline(xintercept=143,colour="black", linetype="dotted",size=0.8)
p2

library(pheatmap)
mycolors <- colorRampPalette(c("white", "#EFC000FF", "#0073C2FF"), bias = 1.2)(100)
tmp=log2(exp_dat+1)
#tmp[tmp > 4] = 4
#tmp[tmp < 2] = 2
annotation_col = data.frame(fp_dat$risk)
colnames(annotation_col) <- "riskScore"
rownames(annotation_col) <- fp_dat$patient
tmp1 <- tmp[, rownames(annotation_col)]
pheatmap(tmp1,annotation_col = annotation_col,col= mycolors,show_colnames = F,cluster_cols = F)


pheatmap(tmp,annotation_col = annotation_col,col= mycolors,show_colnames = F,cluster_cols = F, cluster_rows = F)
p3 <- pheatmap(tmp1,annotation_col = annotation_col,col= mycolors,show_colnames = F,cluster_cols = F)
p3
library(ggplotify)
plots = list(p1,p2,as.ggplot(as.grob(p3)))
library(gridExtra)
lay = rbind(c(rep(1,7)),c(rep(2,7)),c(rep(3,7))) #布局矩阵
grid.arrange(grobs = plots, layout_matrix = lay, heigths = c(2,3,2),weights=c(10,10,10))

f2 <- psm(Surv(time, event)~ age1+gender_group+risk+stage_group, data = riskscore_surv, x=T, y=T, dist='lognormal')

cal1 <- calibrate(f2, 
                  cmethod='KM', 
                  method="boot", 
                  u=1*12, # u需要与之前模型中定义好的time.inc一致，即365或730；
                  m=70, #每次抽样的样本量，
                  B=2000)

cal3 <- calibrate(f2, 
                  cmethod='KM', 
                  method="boot", 
                  u=3*12, # u需要与之前模型中定义好的time.inc一致，即365或730；
                  m=70, #每次抽样的样本量，
                  B=2000)

cal5 <- calibrate(f2, 
                  cmethod='KM', 
                  method="boot", 
                  u=5*12, # u需要与之前模型中定义好的time.inc一致，即365或730；
                  m=70, #每次抽样的样本量，
                  B=2000)


plot(cal1,lwd = 2,lty = 0,errbar.col = c("#2166AC"),
     bty = "l", #只画左边和下边框
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#2166AC"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#2166AC"), pch = 16)
mtext("")

plot(cal3,lwd = 2,lty = 0,errbar.col = c("#B2182B"),
     xlim = c(0,1),ylim= c(0,1),col = c("#B2182B"),add = T)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#B2182B"), pch = 16)

plot(cal5,lwd = 2,lty = 0,errbar.col = c("#EFC000FF"),
     xlim = c(0,1),ylim= c(0,1),col = c("#EFC000FF"),add = T)
lines(cal5[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#EFC000FF"), pch = 16)

abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

legend("topleft", #图例的位置
       legend = c("1-year","3-year", "5-year"), #图例文字
       col =c("#2166AC","#B2182B", "#EFC000FF"), #图例线的颜色，与文字对应
       lwd = 1,#图例中线的粗细
       cex = 1,#图例字体大小
       bty = "n")#不显示图例边框


low.sample <- riskscore_surv %>% filter(riskgroup2 == "low") %>% pull(patient)
high.sample <- riskscore_surv %>% filter(riskgroup2 == "high") %>% pull(patient)
low.exp <- mRNA.data[, low.sample]
high.exp <- mRNA.data[, high.sample]

lh.exp <- cbind.data.frame(low.exp, high.exp)
lh.exp.ct <- 2^lh.exp - 1
lh.exp.ct <- round(lh.exp.ct)

library(DESeq2)
condition <- factor(c(rep("low",143),rep("high",143)))
colData <- data.frame(row.names=colnames(lh.exp.ct), condition)
dds <- DESeqDataSetFromMatrix(countData = lh.exp.ct, colData = colData, design = ~ condition)
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE) 
res <- results(dds1)
# res格式转化：用data.frame转化为表格形式
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
# 依次按照pvalue值log2FoldChange值进行排序
res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
res1$gene <- rownames(res1)
res1$change <- ifelse(res1$log2FoldChange > 1 & res1$padj < 0.05, "down", 
                      ifelse(res1$log2FoldChange < -1 & res1$padj < 0.05, "up", "none"))

x <- as.data.frame(t(lh.exp.ct[c("DDIT4", "NOX4", "GCLC"), ]))

write.table(res1, file = "DEGs_lh.txt", quote = F, sep = "\t")

table(res1$change)

norm_counts <- counts(dds1, normalized = T)
write.table(norm_counts, file = 'DESEQ2_normed_matrix.txt', quote = F, sep = "\t")


GSEAplot.high <- read_tsv(file = "./gsea_report_for_HIGH_1690373068171.tsv")
GSEAplot.low <- read_tsv(file = "./gsea_report_for_LOW_1690373068171.tsv")

GSEAPLOT <- rbind.data.frame(GSEAplot.high, GSEAplot.low)
GSEAPLOT$x <- "1"

GSEAPLOT$name3 <- str_remove(string = GSEAPLOT$NAME, pattern = "HALLMARK_")
GSEAPLOT$name3 <- factor(GSEAPLOT$name3, levels = GSEAPLOT$name3)
GSEAPLOT$SIG <- ifelse(abs(GSEAPLOT$NES) > 1 & GSEAPLOT$`NOM p-val` < 0.05 & GSEAPLOT$`FDR q-val` <0.25, "SIG", "None")


ggplot(GSEAPLOT, aes(x = x, y = name3))+
  geom_point(aes(fill = NES, size = -log10(`FDR q-val`+ 1)), shape = 22)+
  scale_fill_gradient2(low = "cyan", high = "red", midpoint = 0)+
  theme_bw()+
  labs(x="", y="Hallmarks")+
  theme(axis.title.y = element_text(size = 10, color = "black"), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 10, color = "black"))+
  scale_y_discrete(position = "right")


GO <- read_csv(file = "./gProfiler_GO.csv")
GO.BP <- GO %>% filter(source == "GO:BP")
go.plt.bp = data.frame(Category = GO.BP$source,Term = GO.BP$term_name, Genes = GO.BP$intersections, adj_pval = GO.BP$adjusted_p_value)
go.plt.bp <- go.plt.bp[c(14,15,40,48), ]
genelist=data.frame(ID = res1$gene, logFC = res1$log2FoldChange)

row.names(genelist)=genelist[,1]
circ <- circle_dat(go.plt.bp, genelist) 
termNum = 4
termNum=ifelse(nrow(go.plt.bp)<termNum,nrow(go.plt.bp),termNum)
geneNum = nrow(genelist) 

chord <- chord_dat(circ, genelist[1:geneNum,], go.plt.bp$Term[1:termNum])
GOChord(chord,
        space = 0.02,           # 基因之间的间距
        gene.order = 'logFC',    # 排序基因
        gene.space = 0.2,       # 基因离圆圈距离
        gene.size = 3,           # 基因字体大小
        border.size = 0.1,               # 线的大小
        process.label = 10, 
        lfc.col = c("cyan", "red"),
        ribbon.col=brewer.pal(4, "Set3"))


sig.GO <- GO %>% filter(adjusted_p_value <= 0.05)

sig.GO$Term1 <- factor(sig.GO$term_name, levels = sig.GO$term_name)
sig.GO <- sig.GO[1:127, ]
sig.GO10 <- sig.GO[c(1:10, 19:28, 101:110), ]
ggplot(sig.GO10, aes(x=Term1, y=-log10(adjusted_p_value), fill=source))+
  geom_bar(stat="identity",position = "dodge")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.2, hjust = 1))+
  facet_grid(.~source, scales = 'free')+
  theme(legend.position = "none")+
  theme(axis.title.y = element_text(size = 13, color = "black"), 
        axis.text.x = element_text(size = 13, color = "black"), 
        axis.text.y = element_text(size = 13, color = "black"))


LM22 <- read_rds(file = "./TCGA_LM22Data.rds.gz")
LM22.LUAD <- LM22$LM22$LUAD

wt.sample <- substr(colnames(lh.exp), 1, 12)
lm22.wt <- LM22.LUAD[LM22.LUAD$TCGA.Participant.Barcode %in% wt.sample, ]
lm22.wt1 <- lm22.wt[!(is.na(lm22.wt$B.Cells.Memory)), ]
colData$TCGA.Participant.Barcode <- substr(rownames(colData),1, 12)

patient_risk_group <- riskscore_surv[, c(20,21,22,23)]
patient_risk_group$TCGA.Participant.Barcode <- substr(patient_risk_group$patient,1, 12)

lm22.wt2 <- left_join(lm22.wt1, colData, by = "TCGA.Participant.Barcode")
lm22.wt3 <- left_join(lm22.wt1, patient_risk_group, by = "TCGA.Participant.Barcode")
lm22.wt4 <- lm22.wt3[, -c(24, 25)]

lm22.long <- gather(lm22.wt4, key = "cell", value = "proportion", -"TCGA.Participant.Barcode", -"riskgroup2")

ggplot(data = lm22.long, aes(x = riskgroup2, y = proportion, fill = riskgroup2)) +
  geom_boxplot()+
  stat_compare_means(aes(label = paste0("p = ", ..p.format..)), label.x = 1.35, vjust = 0.8, size = 4.5) +
  guides(fill = FALSE)+
  labs(x = "",
       y = "Cell composition fraction")+
  facet_wrap("cell", scales = "free", nrow = 3)+
  theme(axis.title=element_text(face = 'bold',size=12),
        axis.text = element_text(size = 12, colour = 'black'),
        axis.ticks.length=unit(.2, "cm"),
        axis.ticks = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black",size = 0, fill=NA),
        plot.margin = margin(1, 1, 1, 1, "cm"))



ESTIMATE <- read.table(file = "./pics new/LUAD_ESTIMATE.txt", header = T, sep = "\t")
ESTIMATE$patient <- paste0(ESTIMATE$ID, "A")
ESTIMATE$TCGA.Participant.Barcode <- substr(ESTIMATE$ID,1, 12)
ESTIMATE <- merge.data.frame(patient_risk_group, ESTIMATE, by = "TCGA.Participant.Barcode")


ggplot(ESTIMATE, aes(x= riskgroup2, y=Immune_score, color = riskgroup2))+
  geom_boxplot(aes(color = riskgroup2))+
  geom_point(size = 3,shape = 1)+
  stat_compare_means()+
  theme_bw()+
  theme(panel.grid = element_line(linetype = "dashed"))+
  labs(x="", y="Stromal_score")+
  theme(axis.title.y = element_text(size = 14, color = "black"), 
        axis.text.x = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(size = 14, color = "black"),
        legend.position = "none")

library(GSEABase)
negainf <- getGmt("./GOBP_PURINE_CONTAINING_COMPOUND_METABOLIC_PROCESS.v2023.1.Hs.gmt")
negainfpp <- list(negainf = negainf[["GOBP_PURINE_CONTAINING_COMPOUND_METABOLIC_PROCESS"]]@geneIds)
negainfpp.score <- gsva(as.matrix(mRNA.data.tumor.wt1), negainfpp, method="gsva", verbose=FALSE, kcdf = "Gaussian")
negainfpp.score <- as.data.frame(t(negainfpp.score))
negainfpp.score <- cbind.data.frame(negainfpp.score, patient_risk_group)


cor.test(negainfpp.score$negainf, negainfpp.score$risk)
ggplot(negainfpp.score, aes(x= riskgroup2, y=negainf, color = riskgroup2))+
  geom_boxplot(aes(color = riskgroup2))+
  geom_point(size = 3,shape = 1)+
  stat_compare_means(method = "t.test")



negainf.exp <- mRNA.data.tumor.wt1[negainf, ] %>% na.omit() %>% t() %>% as.data.frame()
negainf.exp$risk <- riskscore_surv2$risk


ggplot(negainfpp.score, aes(x= riskgroup, y=negainf, color = riskgroup))+
  geom_boxplot(aes(color = riskgroup))+
  geom_point(size = 3,shape = 1)+
  stat_compare_means()+
  theme_bw()+
  theme(panel.grid = element_line(linetype = "dashed"))+
  labs(x="", y="negainf_score")+
  theme(axis.title.y = element_text(size = 14, color = "black"), 
        axis.text.x = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(size = 14, color = "black"),
        legend.position = "none")

riskscore_surv2 <- riskscore_surv[, c(22,23,20)]
riskscore_surv2$TCGA.Participant.Barcode <- substr(riskscore_surv2$patient, 1, 12)
lm22_riskscore <- merge.data.frame(riskscore_surv2, lm22.wt2, by = 'TCGA.Participant.Barcode')  


cor.res <- cor.test(negainf.exp$ACE, negainf.exp$risk)
cor.table <- data.frame(gene = colnames(negainf.exp)[1:109], corr = rep(1, 109), p.val = rep(1, 109))

cor.res[["p.value"]]; cor.res[["estimate"]]

for (i in c(1:109)) {
  cor.res <- cor.test(negainf.exp$risk, negainf.exp[, i])
  cor.table$corr[i] <- cor.res[["estimate"]][["cor"]]
  cor.table$p.val[i] <- cor.res[["p.value"]]
}

cor.table$relation <- ifelse(cor.table$corr > 0.2 & cor.table$p.val < 0.05, "Positive related", 
                             ifelse(cor.table$corr < -0.2 & cor.table$p.val < 0.05, "Negative related", "Not related"))  

cor.table <- cor.table[-102, ]
cor.table$label <- ifelse(cor.table$relation == "Not related", "", cor.table$gene)

ggplot(cor.table, aes(x=corr, y=-log10(p.val), color = relation))+
  geom_point()+
  scale_color_manual(values = c("#00BFFF", "#999999","#F8766D"))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", color = "#999999")+
  geom_text_repel(aes(label = label))+
  theme(axis.title.y = element_text(size = 14, color = "black"), 
        axis.title.x = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(size = 14, color = "black"),
        legend.position = "none")


FPKM.TIDE <- LUAD.FPKM[, colnames(mRNA.data.tumor.wt2.FerrDb1)]
FPKM.TIDE <- 2^FPKM.TIDE - 1
Expr <- t(apply(FPKM.TIDE, 1, function(x)x-(mean(x))))

IMMSUP.exp <- LUAD.FPKM[c("TIGIT", "PDCD1", "LAG3", "CTLA4", "CD274"), colnames(mRNA.data.tumor.wt2.FerrDb1)] %>% na.omit() %>% t() %>% as.data.frame()
IMMSUP.exp$risk <- patient_risk_group$risk
IMMSUP.exp$riskgroup <- patient_risk_group$riskgroup2

ggplot(IMMSUP.exp, aes(x=riskgroup, y = TIGIT, color = riskgroup))+
  geom_boxplot()+
  stat_compare_means(method = "t.test")+
  theme_bw()+
  labs(x="")+
  theme(panel.grid = element_line(linetype = "dashed"))+
  theme(axis.title.y = element_text(size = 14, color = "black"), 
        axis.title.x = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(size = 14, color = "black"),
        legend.position = "none")

TIDE <- read_csv("./TIDE_LUAD_WT.csv")

TIDE_risk <- merge.data.frame(patient_risk_group, TIDE, by = "patient")
TIDE_risk$Responder1 <- ifelse(TIDE_risk$Responder == T, "Responder", "Nonresponder")
TIDE_risk1 <- TIDE_risk[order(TIDE_risk$risk), ]



ggplot(TIDE_risk, aes(x=Responder1, y = risk, color = Responder1))+
  geom_boxplot()+
  stat_compare_means(method = "t.test")+
  theme_bw()+
  labs(x="")+
  theme(panel.grid = element_line(linetype = "dashed"))+
  theme(axis.title.y = element_text(size = 14, color = "black"), 
        axis.title.x = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(size = 14, color = "black"),
        legend.position = "none")


ggplot(TIDE_risk, aes(x=Responder1, y = MDSC, color = Responder1))+
  geom_boxplot()+
  stat_compare_means(method = "t.test")+
  theme_bw()+
  labs(x="")+
  theme(panel.grid = element_line(linetype = "dashed"))+
  theme(axis.title.y = element_text(size = 14, color = "black"), 
        axis.title.x = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(size = 14, color = "black"),
        legend.position = "none")


stack <- as.data.frame(table(TIDE_risk$Responder1, TIDE_risk$riskgroup2))
stack1 <- plyr::ddply(stack, "Var2", transform, percent=Freq/sum(Freq)*100)
ggplot(stack1, aes(x=Var2, y=percent, fill=Var1))+
  geom_bar(stat = "identity")+
  theme_bw()+
  scale_fill_manual(values = c("#156077", "#f46f20"))+
  theme(panel.grid = element_line(linetype = "dashed"))+
  theme(axis.title.y = element_text(size = 14, color = "black"), 
        axis.title.x = element_text(size = 0, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(size = 14, color = "black"))

forChi <- matrix(c(91,52,47,96), nrow = 2, ncol = 2)
chisq.test(forChi)

ggplot(TIDE_risk, aes(x=riskgroup2, y = TIDE, color = riskgroup2))+
  geom_boxplot()+
  stat_compare_means(method = "t.test")+
  theme_bw()+
  labs(x="")+
  theme(panel.grid = element_line(linetype = "dashed"))+
  theme(axis.title.y = element_text(size = 14, color = "black"), 
        axis.title.x = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(size = 14, color = "black"),
        legend.position = "none")



CAFgene <- c("COL1A1","COL1A2","COL5A1","ACTA2","FGF2","FAP","LRP1","CD248","COL6A1","COL6A2","COL6A3","CXCL12","FBLN1","LUM","MFAP5","MMP3","MMP2","PDGFRB","PDGFRA")

CAFexp <- LUAD.FPKM[CAFgene, riskscore_surv2$patient]
CAFexp1 <- CAFexp[, rownames(annocol.caf)]
annocol.caf <- data.frame(group = riskscore_surv2$riskgroup2, risk = riskscore_surv2$risk)
rownames(annocol.caf) <- riskscore_surv2$patient
annocol.caf <- annocol.caf[order(annocol.caf$risk), ]
annocol.caf$response <- TIDE_risk1$Responder1

groupcolor <- c("blue", "red")
names(groupcolor) <- c("low", "high")
pheatmap(scales::rescale(as.matrix(CAFexp1), to = c(-1, 1)),
         show_colnames = F,
         cluster_cols = F, 
         cluster_rows = T, 
         annotation_col = annocol.caf, 
         color = colorRampPalette(colors = c("#2EA9DF","white","#DB4D6D"))(100),
         annotation_colors = groupcolor
         )

MDSCgene <- unique(c("IDO1","ARG1","IL10","CYBB","PTGS2","IL4I1","IL6","CSF2","CSF3","CXCL12","CCL26","IL6","CXCL8","CXCL5","CSF1R","CSF2RA","CSF3R","CXCR4","IL6R","CXCR2","CCL15","CSF1"))
MDSCexp <- LUAD.FPKM[MDSCgene, riskscore_surv2$patient]
MDSCexp1 <- MDSCexp[, rownames(annocol.caf)]
annocol.MDSC <- data.frame(group = riskscore_surv2$riskgroup2, risk = riskscore_surv2$risk)
rownames(annocol.MDSC) <- riskscore_surv2$patient
annocol.MDSC <- annocol.MDSC[order(annocol.MDSC$risk), ]
annocol.MDSC$response <- TIDE_risk1$Responder1

pheatmap(scales::rescale(as.matrix(MDSCexp1), to = c(-1, 1)),
         show_colnames = F,
         cluster_cols = F, 
         cluster_rows = T, 
         annotation_col = annocol.MDSC, 
         color = colorRampPalette(colors = c("#2EA9DF","white","#DB4D6D"))(100))


lm.test.tide <- TIDE_risk[, c(3,15)]
colnames(lm.test.tide) <- c("risk", "ex")
lm.test.tide1 <- scales::rescale(as.matrix(lm.test.tide), to = c(-1, 1)) %>% as.data.frame()

cor.test(riskscore.ssgsea.test$risk, TIDE_risk$Exclusion)

ggplot(TIDE_risk, aes(x=risk, y=CAF))+
  geom_point(color = "#757575")+
  geom_smooth(method = "lm", formula = y~x)+
  theme_bw()+
  labs(x = "Risk score", y = "CAF")+
  theme(panel.grid = element_line(linetype = "dashed"))+
  theme(axis.title.y = element_text(size = 14, color = "black"), 
        axis.title.x = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(size = 14, color = "black"),
        legend.position = "none")

TME <- read_rds(file = "../ImmType_Mut_20220801.rds.gz")

LUAD.TME <- TME$TMB$LUAD
LUAD.subtype <- TME$Subtypes$LUAD
colnames(LUAD.subtype)[1] <- c('TCGA.Participant.Barcode')
LUAD.subtype.RISK <- merge.data.frame(patient_risk_group, LUAD.subtype, by = "TCGA.Participant.Barcode")

stack.suntype <- as.data.frame(table(LUAD.subtype.RISK$MFP, LUAD.subtype.RISK$riskgroup2))
stack.suntype1 <- plyr::ddply(stack.suntype, "Var2", transform, percent=Freq/sum(Freq)*100)
ggplot(stack.suntype1, aes(x=Var1, y=percent, fill=Var2))+
  geom_bar(stat = "identity")+
  theme_bw()+
  scale_fill_manual(values = c("#156077", "#f46f20", "blue", "red"))+
  theme(panel.grid = element_line(linetype = "dashed"))+
  theme(axis.title.y = element_text(size = 14, color = "black"), 
        axis.title.x = element_text(size = 0, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(size = 14, color = "black"))

ggplot(LUAD.subtype.RISK, aes(x=riskgroup2, y = TMB, color = riskgroup2))+
  geom_boxplot()+
  stat_compare_means()+
  theme_bw()+
  labs(x="")+
  theme(panel.grid = element_line(linetype = "dashed"))+
  theme(axis.title.y = element_text(size = 14, color = "black"), 
        axis.title.x = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(size = 14, color = "black"),
        legend.position = "none")

colnames(TIL) <- paste0(colnames(TIL), "A")
colnames(TIL) <- gsub(pattern = "\\.", replacement = "-", colnames(TIL))
riskscore_surv2 <- riskscore_surv[order(riskscore_surv$riskgroup2), ]
TIL.RISK <- TIL[, which(colnames(TIL) %in% riskscore_surv2$patient)]

TIL.RISK.sample <- data.frame(patient = colnames(TIL.RISK))
TIL.RISK.GROUP <- merge.data.frame(TIL.RISK.sample, patient_risk_group, by = "patient")
TIL.RISK.GROUP <- TIL.RISK.GROUP[order(TIL.RISK.GROUP$riskgroup2), ]
TIL.RISK <- TIL.RISK[, TIL.RISK.GROUP$patient]

library(pheatmap)
mycolors <- colorRampPalette(c("white", "#EFC000FF", "#0073C2FF"), bias = 1.2)(100)

annotation_col_til = data.frame(TIL.RISK.GROUP$riskgroup2)
colnames(annotation_col_til) <- "Riskgroup"
rownames(annotation_col_til) <- TIL.RISK.GROUP$patient

pheatmap(TIL.RISK,annotation_col = annotation_col_til,show_colnames = F,cluster_cols = F, scale = "row", cluster_rows = F)

TIL.RISK1 <- as.data.frame(t(TIL.RISK))
TIL.RISK1$patient <- rownames(TIL.RISK1)
TIL.RISK1$group <- TIL.RISK.GROUP$riskgroup2
TIL.RISK.long <- gather(TIL.RISK1, key = "cell", value = "abundance", -"patient", -"group")

ggplot(data = TIL.RISK.long, aes(x = group, y = abundance, fill = group)) +
  geom_boxplot()+
  stat_compare_means(aes(label = paste0("p = ", ..p.format..)), label.x = 1.35, vjust = 0.8, size = 4.5) +
  guides(fill = FALSE)+
  labs(x = "",
       y = "Cell composition fraction")+
  facet_wrap("cell", scales = "free", nrow = 3)+
  theme(axis.title=element_text(face = 'bold',size=12),
        axis.text = element_text(size = 12, colour = 'black'),
        axis.ticks.length=unit(.2, "cm"),
        axis.ticks = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black",size = 0, fill=NA),
        plot.margin = margin(1, 1, 1, 1, "cm"))



NCI60.exp <- read.table(file = "./NCI60_exp.txt", header = T, sep = "\t")
rownames(NCI60.exp) <- NCI60.exp$Gene.name
NCI60.exp <- NCI60.exp[, -1]

NCI60.RISK <- gsva(as.matrix(NCI60.exp), risklist, method="gsva", verbose=FALSE, kcdf = "Gaussian")
NCI60.RISK1 <- as.data.frame(t(NCI60.RISK))

NCI60.IC50 <- readxl::read_xlsx(path = "./DTP_NCI60_ZSCORE.xlsx")
NCI60.IC50.FDA <- NCI60.IC50[NCI60.IC50$`FDA status` != "-", ]
library(limma)
NCI60.IC50.FDA <- avereps(NCI60.IC50.FDA[, -c(1,3,4,5,6, 67,68)], ID = NCI60.IC50.FDA$`Drug name`)
NCI60.IC50.FDA <- as.data.frame(NCI60.IC50.FDA)
rownames(NCI60.IC50.FDA) <- NCI60.IC50.FDA$`Drug name`
NCI60.IC50.FDA <- NCI60.IC50.FDA[, -1]


library(impute)
mat <- as.matrix(NCI60.IC50.FDA)
mat <- matrix(as.numeric(mat), ncol = 60, nrow = 767, dimnames = list(rownames = rownames(NCI60.IC50.FDA), colnames = colnames(NCI60.IC50.FDA)))
mat1 <- impute.knn(mat)

NCI60.RISK2 <- cbind(NCI60.RISK1, as.data.frame(t(mat)))

NCI60.RISK.COR <- data.frame(Drug = colnames(NCI60.RISK2)[2:768], corr = rep(1, 767), p.val = rep(1, 767))
for (i in c(1:767)) {
  cor.res <- cor.test(NCI60.RISK2$risk, NCI60.RISK2[, i+1])
  NCI60.RISK.COR$corr[i] <- cor.res[["estimate"]][["cor"]]
  NCI60.RISK.COR$p.val[i] <- cor.res[["p.value"]]
}

NCI60.RISK.COR$Sig <- ifelse(abs(NCI60.RISK.COR$corr) > 0.2 & NCI60.RISK.COR$p.val < 0.05, "Sig", "None")
NCI60.RISK.COR.sig <- NCI60.RISK.COR %>% filter(Sig == "Sig")
NCI60.RISK.COR.sig1 <- NCI60.RISK.COR.sig[order(NCI60.RISK.COR.sig$corr, decreasing = F), ]
NCI60.RISK.COR.sig1$targets <- Targets
NCI60.RISK.COR.sig2 <- NCI60.RISK.COR.sig1[-c(3,4,7), ]

ggplot(NCI60.RISK.COR.sig2[-3, ], aes(x=Drug, y=targets, color = corr, size = -log10(p.val + 1)))+
  geom_point()+
  scale_color_gradient2(low = "#74759b", high = "red", midpoint = 0)+
  theme_bw()+
  theme(panel.grid = element_line(linetype = "dashed"))+
  theme(axis.text.x = element_text(angle = 90,vjust = 1,hjust = 1))+
  theme(axis.title.y = element_text(size = 14, color = "black"), 
        axis.title.x = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(size = 14, color = "black"))+
  labs(x="", y="")
  

cmap <- read_table(file = "./top150_cmap.txt")

cmap_50 <- cmap[cmap$Score < -50, ]



MDSCgene <- unique(c("HPRT1", "IMPDH1", "IMPDH2", "PPAT"))
MDSCexp <- LUAD.FPKM[MDSCgene, riskscore_surv2$patient]
MDSCexp1 <- MDSCexp[, rownames(annocol.caf)]
annocol.MDSC <- data.frame(group = riskscore_surv2$riskgroup2, risk = riskscore_surv2$risk)
rownames(annocol.MDSC) <- riskscore_surv2$patient
annocol.MDSC <- annocol.MDSC[order(annocol.MDSC$risk), ]
annocol.MDSC$response <- TIDE_risk1$Responder1

pheatmap(scales::rescale(as.matrix(MDSCexp1), to = c(-1, 1)),
         show_colnames = F,
         cluster_cols = F, 
         cluster_rows = T, 
         annotation_col = annocol.MDSC, 
         color = colorRampPalette(colors = c("#2EA9DF","white","#DB4D6D"))(100))


library(cBioPortalData)
cbio <- cBioPortal()
studies = getStudies(cbio)

clinical = clinicalData(cbio, 'luad_tcga_pan_can_atlas_2018')
LUAD.MSI <- data.frame(patient = clinical$patientId, MSI = clinical$MSI_SCORE_MANTIS, TMB = clinical$TMB_NONSYNONYMOUS)
LUAD.MSI$patient <- paste0(LUAD.MSI$patient, "-01A")
LUAD.MSI <- merge.data.frame(LUAD.MSI, patient_risk_group, by = "patient")
LUAD.MSI$MSI <- as.numeric(LUAD.MSI$MSI)
LUAD.MSI$TMB <- as.numeric(LUAD.MSI$TMB)


ggplot(LUAD.MSI, aes(x=riskgroup2, y = MSI, color = riskgroup2))+
  geom_boxplot()+
  stat_compare_means(method = "t.test")+
  theme_bw()+
  labs(x="")+
  theme(panel.grid = element_line(linetype = "dashed"))+
  theme(axis.title.y = element_text(size = 14, color = "black"), 
        axis.title.x = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(size = 14, color = "black"),
        legend.position = "none")

ggplot(LUAD.MSI, aes(x=riskgroup2, y = TMB, color = riskgroup2))+
  geom_boxplot()+
  stat_compare_means(method = "t.test")+
  theme_bw()+
  labs(x="")+
  theme(panel.grid = element_line(linetype = "dashed"))+
  theme(axis.title.y = element_text(size = 14, color = "black"), 
        axis.title.x = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(size = 14, color = "black"),
        legend.position = "none")

celline <- read.table(file = "./Cell_line_exp_array.txt", header = T, sep = "\t")
GDSC_AUC <- read.table(file = "./GDSC_AUC.txt", header = T, sep = "\t")

celline1 <- as.data.frame(avereps(celline[, -1], ID = celline$GENE_SYMBOLS))

celline_risk <- gsva(as.matrix(celline1), risklist, method="gsva", verbose=FALSE, kcdf = "Gaussian")
celline_risk <- as.data.frame(t(celline_risk))

my_comparisons <- list( c("D", "F"), c("D", "IE"), c("D", "IE/F"),  c("F", "IE"), c("F", "IE/F"), c("IE", "IE/F"))
ggplot(LUAD.subtype.RISK, aes(x=MFP, y = risk, color = MFP))+
  geom_boxplot()+
  stat_compare_means(comparisons = my_comparisons)+
  theme_bw()+
  labs(x="")+
  theme(panel.grid = element_line(linetype = "dashed"))+
  theme(axis.title.y = element_text(size = 14, color = "black"), 
        axis.title.x = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(size = 14, color = "black"),
        legend.position = "none")

TME.list <- list(MHCI = c("HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2","TAPBP"),
                 MHCII = c("HLA-DRA","HLA-DRB1","HLA-DMA","HLA-DPA1","HLA-DPB1","HLA-DMB","HLA-DQB1","HLA-DQA1","CIITA"),
                 Co_activation_molecules=c("CD28","CD40","TNFRSF4","ICOS","TNFRSF9","CD27","CD80","CD86","CD40LG","CD83","TNFSF4","ICOSLG","TNFSF9","CD70"),
                 Effector_cells = c("IFNG","GZMA","GZMB","PRF1","GZMK","ZAP70","GNLY","FASLG","TBX21","EOMES","CD8A","CD8B"),
                 Effector_cell_traffic = c("CXCL9","CXCL10","CXCL11","CX3CL1","CCL3","CCL4","CX3CR1","CCL5","CXCR3"),
                 NK_cells = c("NKG7","CD160","CD244","NCR1","KLRC2","KLRK1","CD226","GZMH","GNLY","IFNG","KIR2DL4","EOMES","GZMB","FGFBP2","KLRF1","SH2D1B","NCR3"),
                 T_cells = c("TBX21","ITK","CD3D","CD3E","CD3G","TRAC","TRBC1","TRBC2","CD28","CD5","TRAT1"),
                 B_cells = c("CD19","MS4A1","TNFRSF13C","CR2","TNFRSF17","TNFRSF13B","CD22","CD79A","CD79B","BLK","FCRL5","PAX5","STAP1"),
                 M1_signature=c("NOS2","TNF","IL1B","SOCS3","CMKLR1","IRF5","IL12A","IL12B","IL23A"),
                 Th1_signature=c("IFNG","IL2","CD40LG","IL21","TBX21","STAT4","IL12RB2"),
                 Antitumor_cytokines=c("TNF","IFNB1","IFNA2","CCL3","TNFSF10","IL21"),
                 Checkpoint_molecules=c("PDCD1","CD274","CTLA4","LAG3","PDCD1LG2","BTLA","HAVCR2","TIGIT","VSIR","C10orf54"),
                 Treg=c("FOXP3","CTLA4","IL10","TNFRSF18","CCR8","IKZF4","IKZF2"),
                 Treg_and_Th2_traffic=c("CCL17","CCL22","CCL1","CCL28","CCR4","CCR8","CCR10"),
                 Neutrophil_signature=c("MPO","ELANE","PRTN3","CTSG","CXCR1","CXCR2","FCGR3B","CD177","FFAR2","PGLYRP1"),
                 Granulocyte_traffic=c("CXCL8","CXCL1","CXCL2","CXCL5","CCL11","KITLG","CXCR1","CXCR2","CCR3"),
                 Immune_Suppression_by_Myeloid_Cells=c("IDO1","ARG1","IL10","CYBB","PTGS2","IL4I1","IL6"),
                 Myeloid_cells_traffic=c("CSF2","CSF3","CXCL12","CCL26","IL6","CXCL8","CXCL5","CSF1R","CSF2RA","CSF3R","CXCR4","IL6R","CXCR2","CCL15","CSF1"),
                 Tumor_associated_Macrophages=c("IL10","MRC1","MSR1","CD163","CSF1R","IL4I1","SIGLEC1","CD68"),
                 Macrophage_and_DC_traffic=c("CCL2","CCL7","CCL8","XCL1","CCR2","XCR1","CSF1R","CSF1"),
                 Th2_signature=c("IL4","IL5","IL13","IL10","GATA3","CCR4"),
                 Protumor_cytokines=c("IL10","TGFB1","TGFB2","TGFB3","IL22","MIF","IL6"),
                 Cancer_associated_fibroblasts=c("COL1A1","COL1A2","COL5A1","ACTA2","FGF2","FAP","LRP1","CD248","COL6A1","COL6A2","COL6A3","CXCL12","FBLN1","LUM","MFAP5","MMP3","MMP2","PDGFRB","PDGFRA"),
                 Matrix=c("FN1","COL1A1","COL1A2","COL4A1","COL3A1","VTN","LGALS7","LGALS9","LAMA3","LAMB3","LAMC2","TNC","ELN","COL5A1","COL11A1"),
                 Matrix_remodeling=c("CA9","MMP9","MMP2","MMP1","MMP3","MMP12","MMP7","MMP11","PLOD2","ADAMTS4","ADAMTS5","LOX"),
                 Angiogenesis=c("VEGFA","VEGFB","VEGFC","PDGFC","CXCL8","CXCR2","FLT1","PGF","CXCL5","KDR","ANGPT1","ANGPT2","TEK","VWF","CDH5"),
                 Endothelium=c("NOS3","KDR","FLT1","VCAM1","VWF","CDH5","MMRN1","ENG","CLEC14A","MMRN2"),
                 Tumor_proliferation_rate=c("MKI67","ESCO2","CETN3","CDK2","CCND1","CCNE1","AURKA","AURKB","E2F1","MYBL2","BUB1","PLK1","CCNB1","MCM2","MCM6"),
                 EMT_signature=c("SNAI1","SNAI2","TWIST1","TWIST2","ZEB1","ZEB2","CDH2"))

LUAD.myTME <- gsva(as.matrix(LUAD.FPKM), TME.list, method="ssgsea", verbose=FALSE)
LUAD.myTME <- as.data.frame(t(LUAD.myTME))
LUAD.myTME.risk <- LUAD.myTME[patient_risk_group$patient, ]
LUAD.myTME.risk <- cbind(LUAD.myTME.risk, patient_risk_group[, c(1,2)])

TME.RISK.COR <- data.frame(Feature = colnames(LUAD.myTME.risk)[1:29], corr = rep(1, 29), p.val = rep(1, 29))
for (i in c(1:29)) {
  cor.res <- cor.test(LUAD.myTME.risk$risk, LUAD.myTME.risk[, i])
  TME.RISK.COR$corr[i] <- cor.res[["estimate"]][["cor"]]
  TME.RISK.COR$p.val[i] <- cor.res[["p.value"]]
}



TME.res.nonres <- merge.data.frame(LUAD.subtype.RISK, TIDE_risk1, by = "TCGA.Participant.Barcode")

library(rmda)
library(caret)
library(rms)


riskscore_surv$riskgroup_3 <- ifelse(riskscore_surv$riskgroup2 == "high", 1, 0)
riskscore_surv$stage_group_2 <- ifelse(riskscore_surv$stage_group %in% c("IB"), 1, 
                                     ifelse(riskscore_surv$stage_group %in% c("II", "IIA", "IIB"), 2, 3))

full.model <- decision_curve(riskgroup_3 ~ age1 + gender_group + event + time + stage_group_2,
                             data = riskscore_surv,
                             thresholds = seq(0,1, by = 0.05),
                             bootstraps = 10)


plot_decision_curve(full.model,
                    curve.names ="Ful1 Model",
                    cost.benefit.axis =FALSE,
                    #col= c('red'blue'),
                    confidence.intervals=FALSE,
                    standardize = T)


###stage ~ group stackplot ----
stage_stack <- riskscore_surv[, c(20,21,26)]

stage_stack_plt <- as.data.frame(table(stage_stack$stage.short, stage_stack$riskgroup2))
stage_stack_plt1 <- plyr::ddply(stage_stack_plt, "Var2", transform, percent=Freq/sum(Freq)*100)
ggplot(stage_stack_plt1, aes(x=Var2, y=percent, fill=Var1))+
  geom_bar(stat = "identity")+
  theme_bw()+
  scale_fill_manual(values = c("#a8d8ea","#fcbad3", "#aa96da"))+
  theme(panel.grid = element_line(linetype = "dashed"))+
  theme(axis.title.y = element_text(size = 14, color = "black"), 
        axis.title.x = element_text(size = 0, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(size = 14, color = "black"))

stage_risk_chi <- matrix(stage_stack_plt1$Freq, nrow = 2,byrow = T)
chisq.test(stage_risk_chi, correct = F)

ggplot(stage_stack, aes(x=stage.short, y=risk, fill=stage.short))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values = c("#a8d8ea","#fcbad3", "#aa96da"))+
  theme(panel.grid = element_line(linetype = "dashed"))+
  theme(axis.title.y = element_text(size = 14, color = "black"), 
        axis.title.x = element_text(size = 0, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(size = 14, color = "black"))+
  stat_compare_means(comparisons = list(c("I", "II"), c("I", "III")))

GSEA_OR <- GSEAPLOT
GSEA_OR$x <- "TCGA Cohort"
GSEA_OR$name3 <- str_to_sentence(gsub(pattern = "_", replacement = " ", x = GSEA_OR$name3))
GSEA_OR <- GSEA_OR[, c(13,14,6,7,15)]
colnames(GSEA_OR)[c(3,4)] <- c("NES_TCGA", "NOM_q_TCGA")

GSEA_V1 <- read_tsv("./GSEA_V1.tsv")
GSEA_V1$x <- "GSE31210"
GSEA_V1$name3 <- str_remove(string = GSEA_V1$NAME, pattern = "HALLMARK_")
GSEA_V1$name3 <- str_to_sentence(gsub(pattern = "_", replacement = " ", x = GSEA_V1$name3))
GSEA_V1$SIG <- ifelse(abs(GSEA_V1$NES) > 1 & GSEA_V1$`NOM p-val` < 0.05 & GSEA_V1$`FDR q-val` <0.25, "SIG", "None")
GSEA_V1 <- GSEA_V1[, c(13,14,6,7,15)]
colnames(GSEA_V1)[c(3,4)] <- c("NES_V1", "NOM_q_V1")

GSEA_V2 <- read_tsv("./GSEA_V2.tsv")
GSEA_V2$x <- "oncoSG"
GSEA_V2$name3 <- str_remove(string = GSEA_V2$NAME, pattern = "HALLMARK_")
GSEA_V2$name3 <- str_to_sentence(gsub(pattern = "_", replacement = " ", x = GSEA_V2$name3))
GSEA_V2$SIG <- ifelse(abs(GSEA_V2$NES) > 1 & GSEA_V2$`NOM p-val` < 0.05 & GSEA_V2$`FDR q-val` <0.25, "SIG", "None")
GSEA_V2 <- GSEA_V2[, c(13,14,6,7,15)]
colnames(GSEA_V2)[c(3,4)] <- c("NES_V2", "NOM_q_V2")

GSEA_V4 <- read_tsv("./GSEA_V4.tsv")
GSEA_V4$x <- "GSE11969"
GSEA_V4$name3 <- str_remove(string = GSEA_V4$NAME, pattern = "HALLMARK_")
GSEA_V4$name3 <- str_to_sentence(gsub(pattern = "_", replacement = " ", x = GSEA_V4$name3))
GSEA_V4$SIG <- ifelse(abs(GSEA_V4$NES) > 1 & GSEA_V4$`NOM p-val` < 0.05 & GSEA_V4$`FDR q-val` <0.25, "SIG", "None")
GSEA_V4 <- GSEA_V4[, c(13,14,6,7,15)]
colnames(GSEA_V4)[c(3,4)] <- c("NES_V4", "NOM_q_V4")

NES.rank <- merge.data.frame(GSEA_OR[, c(2:4)], GSEA_V1[, c(2:4)])
NES.rank <- merge.data.frame(NES.rank, GSEA_V2[,c(2:4)])
NES.rank <- merge.data.frame(NES.rank, GSEA_V4[,c(2:4)])
NES.rank$mean.NES <- rowSums(NES.rank[, c(2,4,6,8)])/4
NES.rank <- NES.rank[order(NES.rank$mean.NES, decreasing = T), ]

NES.plot <- NES.rank[, c(1:9)]
NES.plot_nes <- gather(NES.plot[,c(1,2,4,6,8)], key="key", value = "value", -"name3")
NES.plot_q <- gather(NES.plot[,c(1,3,5,7,9)], key="qkey", value = "qvalue", -"name3")
NES.dotplot <- cbind(NES.plot_nes, NES.plot_q[, c(2,3)])
NES.dotplot <- NES.dotplot[, -4]
NES.dotplot$name3 <- factor(NES.dotplot$name3, levels = rev(NES.rank$name3))
ggplot(NES.dotplot, aes(x = key, y = name3))+
  geom_point(aes(color = value, size = -log10(qvalue + 1)), shape = 16)+
  scale_color_gradient2(low = "#7ac7c4", high = "#f73859", midpoint = 0)+
  theme_bw()+
  labs(x="", y="Hallmarks")+
  theme(axis.title.y = element_text(size = 10, color = "black"), 
        axis.text.x = element_text(size = 10, color = "black", angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(size = 10, color = "black"))+
  scale_y_discrete(position = "left")+
  theme(legend.position = "bottom")


NES.rank$name3<- factor(NES.rank$name3, levels = rev(NES.rank$name3))
ggplot(NES.rank, aes(y=name3, x=mean.NES+1.5))+
  geom_bar(stat="identity",position="dodge", fill = "#c9d6df")+
  geom_text(aes(x=mean.NES+1.6,y=name3,label= round(mean.NES, digits = 2)))+
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(), 
        panel.grid = element_blank())
p1+p2


### DDIT4 ~ EMT genes----
library(clusterProfiler)
Hallmark50 <- read.gmt("./h.all.v2023.1.Hs.symbols.gmt") 
Hallmark <- lapply(unique(Hallmark50$term), function(x){print(x);Hallmark50$gene[Hallmark50$term == x]})
names(Hallmark) <- unique(Hallmark50$term)

EMT.exp <- na.omit(mRNA.data.tumor.wt1[c("DDIT4", Hallmark$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION), ]) %>% t() %>% as.data.frame
EMT.cor.table.DDIT4 <- data.frame(gene = colnames(EMT.exp), corr = rep(1, ncol(EMT.exp)), p.val = rep(1, ncol(EMT.exp)))

for (i in c(1:ncol(EMT.exp))) {
  cor.res <- cor.test(EMT.exp$DDIT4, EMT.exp[, i], method = "spearman")
  EMT.cor.table.DDIT4$corr[i] <- cor.res[["estimate"]][["rho"]]
  EMT.cor.table.DDIT4$p.val[i] <- cor.res[["p.value"]]
}
EMT.cor.table.DDIT4$relation <- ifelse(EMT.cor.table.DDIT4$corr > 0.3 & EMT.cor.table.DDIT4$p.val < 0.01, "Positive", 
                                       ifelse(EMT.cor.table.DDIT4$corr < -0.3 & EMT.cor.table.DDIT4$p.val < 0.01, "Negative", "None"))  

pie.emt <- data.frame(table(EMT.cor.table.DDIT4$relation))
pie.emt$pct <- round((pie.emt$Freq/sum(pie.emt$Freq) * 100),1)
ggplot(pie.emt, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 1) +    
  coord_polar(theta = "y") + 
  labs(x = "", y = "", title = "") + 
  theme(axis.ticks = element_blank()) + 
  theme(legend.title = element_blank(), legend.position = "top") + 
  scale_fill_manual(values = c("#d3d4d8","#ffaaa5")) + 
  theme(axis.text.x = element_blank(),
        panel.background = element_blank(), 
        panel.grid = element_blank())

### DDIT4 ~ HYPOXIA ----
HPX.exp <- na.omit(mRNA.data.tumor.wt1[c("DDIT4", Hallmark$HALLMARK_HYPOXIA), ]) %>% t() %>% as.data.frame
HPX.cor.table.DDIT4 <- data.frame(gene = colnames(HPX.exp), corr = rep(1, ncol(HPX.exp)), p.val = rep(1, ncol(HPX.exp)))

for (i in c(1:ncol(HPX.exp))) {
  cor.res <- cor.test(HPX.exp$DDIT4, HPX.exp[, i])
  HPX.cor.table.DDIT4$corr[i] <- cor.res[["estimate"]][["cor"]]
  HPX.cor.table.DDIT4$p.val[i] <- cor.res[["p.value"]]
}
HPX.cor.table.DDIT4$relation <- ifelse(HPX.cor.table.DDIT4$corr > 0.3 & HPX.cor.table.DDIT4$p.val < 0.01, "Positive", 
                                       ifelse(HPX.cor.table.DDIT4$corr < -0.3 & HPX.cor.table.DDIT4$p.val < 0.01, "Negative", "None"))  

pie.HPX <- data.frame(table(HPX.cor.table.DDIT4$relation))
pie.HPX$pct <- round((pie.HPX$Freq/sum(pie.HPX$Freq) * 100),1)
ggplot(pie.HPX, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 1) +    
  coord_polar(theta = "y") + 
  labs(x = "", y = "", title = "") + 
  theme(axis.ticks = element_blank()) + 
  theme(legend.title = element_blank(), legend.position = "top") + 
  scale_fill_manual(values = c("#d3d4d8","#ffaaa5")) + 
  theme(axis.text.x = element_blank(),
        panel.background = element_blank(), 
        panel.grid = element_blank())

### DDIT4 ~ GLYCOLYSIS ----
GLY.exp <- na.omit(mRNA.data.tumor.wt1[c("DDIT4", Hallmark$HALLMARK_GLYCOLYSIS), ]) %>% t() %>% as.data.frame
GLY.cor.table.DDIT4 <- data.frame(gene = colnames(GLY.exp), corr = rep(1, ncol(GLY.exp)), p.val = rep(1, ncol(GLY.exp)))

for (i in c(1:ncol(GLY.exp))) {
  cor.res <- cor.test(GLY.exp$DDIT4, GLY.exp[, i])
  GLY.cor.table.DDIT4$corr[i] <- cor.res[["estimate"]][["cor"]]
  GLY.cor.table.DDIT4$p.val[i] <- cor.res[["p.value"]]
}
GLY.cor.table.DDIT4$relation <- ifelse(GLY.cor.table.DDIT4$corr > 0.3 & GLY.cor.table.DDIT4$p.val < 0.01, "Positive", 
                                       ifelse(GLY.cor.table.DDIT4$corr < -0.3 & GLY.cor.table.DDIT4$p.val < 0.01, "Negative", "None"))  

pie.GLY <- data.frame(table(GLY.cor.table.DDIT4$relation))
pie.GLY$pct <- round((pie.GLY$Freq/sum(pie.GLY$Freq) * 100),1)
ggplot(pie.GLY, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 1) +    
  coord_polar(theta = "y") + 
  labs(x = "", y = "", title = "") + 
  theme(axis.ticks = element_blank()) + 
  theme(legend.title = element_blank(), legend.position = "top") + 
  scale_fill_manual(values = c("#d3d4d8","#ffaaa5")) + 
  theme(axis.text.x = element_blank(),
        panel.background = element_blank(), 
        panel.grid = element_blank())

### DDIT4 ~ MTOC1 ----
MTOC.exp <- na.omit(mRNA.data.tumor.wt1[c("DDIT4", Hallmark$HALLMARK_MTORC1_SIGNALING), ]) %>% t() %>% as.data.frame
MTOC.cor.table.DDIT4 <- data.frame(gene = colnames(MTOC.exp), corr = rep(1, ncol(MTOC.exp)), p.val = rep(1, ncol(MTOC.exp)))

for (i in c(1:ncol(MTOC.exp))) {
  cor.res <- cor.test(MTOC.exp$DDIT4, MTOC.exp[, i])
  MTOC.cor.table.DDIT4$corr[i] <- cor.res[["estimate"]][["cor"]]
  MTOC.cor.table.DDIT4$p.val[i] <- cor.res[["p.value"]]
}
MTOC.cor.table.DDIT4$relation <- ifelse(MTOC.cor.table.DDIT4$corr > 0.3 & MTOC.cor.table.DDIT4$p.val < 0.01, "Positive", 
                                        ifelse(MTOC.cor.table.DDIT4$corr < -0.3 & MTOC.cor.table.DDIT4$p.val < 0.01, "Negative", "None"))  

pie.MTOC <- data.frame(table(MTOC.cor.table.DDIT4$relation))
pie.MTOC$pct <- round((pie.MTOC$Freq/sum(pie.MTOC$Freq) * 100),1)
ggplot(pie.MTOC, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 1) +    
  coord_polar(theta = "y") + 
  labs(x = "", y = "", title = "") + 
  theme(axis.ticks = element_blank()) + 
  theme(legend.title = element_blank(), legend.position = "top") + 
  scale_fill_manual(values = c("#d3d4d8","#ffaaa5")) + 
  theme(axis.text.x = element_blank(),
        panel.background = element_blank(), 
        panel.grid = element_blank())

### DDIT4 ~ ANGIO ----
AGG.exp <- na.omit(mRNA.data.tumor.wt1[c("DDIT4", Hallmark$HALLMARK_ANGIOGENESIS), ]) %>% t() %>% as.data.frame
AGG.cor.table.DDIT4 <- data.frame(gene = colnames(AGG.exp), corr = rep(1, ncol(AGG.exp)), p.val = rep(1, ncol(AGG.exp)))

for (i in c(1:ncol(AGG.exp))) {
  cor.res <- cor.test(AGG.exp$DDIT4, AGG.exp[, i])
  AGG.cor.table.DDIT4$corr[i] <- cor.res[["estimate"]][["cor"]]
  AGG.cor.table.DDIT4$p.val[i] <- cor.res[["p.value"]]
}
AGG.cor.table.DDIT4$relation <- ifelse(AGG.cor.table.DDIT4$corr > 0.3 & AGG.cor.table.DDIT4$p.val < 0.01, "Positive", 
                                       ifelse(AGG.cor.table.DDIT4$corr < -0.3 & AGG.cor.table.DDIT4$p.val < 0.01, "Negative", "None"))  

pie.AGG <- data.frame(table(AGG.cor.table.DDIT4$relation))
pie.AGG$pct <- round((pie.AGG$Freq/sum(pie.AGG$Freq) * 100),1)
ggplot(pie.AGG, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 1) +    
  coord_polar(theta = "y") + 
  labs(x = "", y = "", title = "") + 
  theme(axis.ticks = element_blank()) + 
  theme(legend.title = element_blank(), legend.position = "top") + 
  scale_fill_manual(values = c("#d3d4d8","#ffaaa5")) + 
  theme(axis.text.x = element_blank(),
        panel.background = element_blank(), 
        panel.grid = element_blank())



### NOX4 ~ EMT genes----
EMT.exp <- na.omit(mRNA.data.tumor.wt1[c("NOX4", Hallmark$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION), ]) %>% t() %>% as.data.frame
EMT.cor.table.NOX4 <- data.frame(gene = colnames(EMT.exp), corr = rep(1, ncol(EMT.exp)), p.val = rep(1, ncol(EMT.exp)))

for (i in c(1:ncol(EMT.exp))) {
  cor.res <- cor.test(EMT.exp$NOX4, EMT.exp[, i], method = "spearman")
  EMT.cor.table.NOX4$corr[i] <- cor.res[["estimate"]][["rho"]]
  EMT.cor.table.NOX4$p.val[i] <- cor.res[["p.value"]]
}
EMT.cor.table.NOX4$relation <- ifelse(EMT.cor.table.NOX4$corr > 0.3 & EMT.cor.table.NOX4$p.val < 0.01, "Positive", 
                                      ifelse(EMT.cor.table.NOX4$corr < -0.3 & EMT.cor.table.NOX4$p.val < 0.01, "Negative", "None"))  

NOX4.pie.emt <- data.frame(table(EMT.cor.table.NOX4$relation))
NOX4.pie.emt$pct <- round((NOX4.pie.emt$Freq/sum(NOX4.pie.emt$Freq) * 100),1)
ggplot(NOX4.pie.emt, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 1) +    
  coord_polar(theta = "y") + 
  labs(x = "", y = "", title = "") + 
  theme(axis.ticks = element_blank()) + 
  theme(legend.title = element_blank(), legend.position = "top") + 
  scale_fill_manual(values = c("#d3d4d8","#ffaaa5")) + 
  theme(axis.text.x = element_blank(),
        panel.background = element_blank(), 
        panel.grid = element_blank())

### NOX4 ~ HYPOXIA ----
HPX.exp <- na.omit(mRNA.data.tumor.wt1[c("NOX4", Hallmark$HALLMARK_HYPOXIA), ]) %>% t() %>% as.data.frame
HPX.cor.table.NOX4 <- data.frame(gene = colnames(HPX.exp), corr = rep(1, ncol(HPX.exp)), p.val = rep(1, ncol(HPX.exp)))

for (i in c(1:ncol(HPX.exp))) {
  cor.res <- cor.test(HPX.exp$NOX4, HPX.exp[, i])
  HPX.cor.table.NOX4$corr[i] <- cor.res[["estimate"]][["cor"]]
  HPX.cor.table.NOX4$p.val[i] <- cor.res[["p.value"]]
}
HPX.cor.table.NOX4$relation <- ifelse(HPX.cor.table.NOX4$corr > 0.3 & HPX.cor.table.NOX4$p.val < 0.01, "Positive", 
                                      ifelse(HPX.cor.table.NOX4$corr < -0.3 & HPX.cor.table.NOX4$p.val < 0.01, "Negative", "None"))  

NOX4.pie.HPX <- data.frame(table(HPX.cor.table.NOX4$relation))
NOX4.pie.HPX$pct <- round((NOX4.pie.HPX$Freq/sum(NOX4.pie.HPX$Freq) * 100),1)
ggplot(NOX4.pie.HPX, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 1) +    
  coord_polar(theta = "y") + 
  labs(x = "", y = "", title = "") + 
  theme(axis.ticks = element_blank()) + 
  theme(legend.title = element_blank(), legend.position = "top") + 
  scale_fill_manual(values = c("#d3d4d8","#ffaaa5")) + 
  theme(axis.text.x = element_blank(),
        panel.background = element_blank(), 
        panel.grid = element_blank())

### NOX4 ~ GLYCOLYSIS ----
GLY.exp <- na.omit(mRNA.data.tumor.wt1[c("NOX4", Hallmark$HALLMARK_GLYCOLYSIS), ]) %>% t() %>% as.data.frame
GLY.cor.table.NOX4 <- data.frame(gene = colnames(GLY.exp), corr = rep(1, ncol(GLY.exp)), p.val = rep(1, ncol(GLY.exp)))

for (i in c(1:ncol(GLY.exp))) {
  cor.res <- cor.test(GLY.exp$NOX4, GLY.exp[, i])
  GLY.cor.table.NOX4$corr[i] <- cor.res[["estimate"]][["cor"]]
  GLY.cor.table.NOX4$p.val[i] <- cor.res[["p.value"]]
}
GLY.cor.table.NOX4$relation <- ifelse(GLY.cor.table.NOX4$corr > 0.3 & GLY.cor.table.NOX4$p.val < 0.01, "Positive", 
                                      ifelse(GLY.cor.table.NOX4$corr < -0.3 & GLY.cor.table.NOX4$p.val < 0.01, "Negative", "None"))  

NOX4.pie.GLY <- data.frame(table(GLY.cor.table.NOX4$relation))
NOX4.pie.GLY$pct <- round((NOX4.pie.GLY$Freq/sum(NOX4.pie.GLY$Freq) * 100),1)
ggplot(NOX4.pie.GLY, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 1) +    
  coord_polar(theta = "y") + 
  labs(x = "", y = "", title = "") + 
  theme(axis.ticks = element_blank()) + 
  theme(legend.title = element_blank(), legend.position = "top") + 
  scale_fill_manual(values = c("#d3d4d8","#ffaaa5")) + 
  theme(axis.text.x = element_blank(),
        panel.background = element_blank(), 
        panel.grid = element_blank())

### NOX4 ~ MTOC1 ----
MTOC.exp <- na.omit(mRNA.data.tumor.wt1[c("NOX4", Hallmark$HALLMARK_MTORC1_SIGNALING), ]) %>% t() %>% as.data.frame
MTOC.cor.table.NOX4 <- data.frame(gene = colnames(MTOC.exp), corr = rep(1, ncol(MTOC.exp)), p.val = rep(1, ncol(MTOC.exp)))

for (i in c(1:ncol(MTOC.exp))) {
  cor.res <- cor.test(MTOC.exp$NOX4, MTOC.exp[, i])
  MTOC.cor.table.NOX4$corr[i] <- cor.res[["estimate"]][["cor"]]
  MTOC.cor.table.NOX4$p.val[i] <- cor.res[["p.value"]]
}
MTOC.cor.table.NOX4$relation <- ifelse(MTOC.cor.table.NOX4$corr > 0.3 & MTOC.cor.table.NOX4$p.val < 0.01, "Positive", 
                                       ifelse(MTOC.cor.table.NOX4$corr < -0.3 & MTOC.cor.table.NOX4$p.val < 0.01, "Negative", "None"))  

NOX4.pie.MTOC <- data.frame(table(MTOC.cor.table.NOX4$relation))
NOX4.pie.MTOC$pct <- round((NOX4.pie.MTOC$Freq/sum(NOX4.pie.MTOC$Freq) * 100),1)
ggplot(NOX4.pie.MTOC, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 1) +    
  coord_polar(theta = "y") + 
  labs(x = "", y = "", title = "") + 
  theme(axis.ticks = element_blank()) + 
  theme(legend.title = element_blank(), legend.position = "top") + 
  scale_fill_manual(values = c("#d3d4d8","#ffaaa5")) + 
  theme(axis.text.x = element_blank(),
        panel.background = element_blank(), 
        panel.grid = element_blank())

### NOX4 ~ ANGIO ----
AGG.exp <- na.omit(mRNA.data.tumor.wt1[c("NOX4", Hallmark$HALLMARK_ANGIOGENESIS), ]) %>% t() %>% as.data.frame
AGG.cor.table.NOX4 <- data.frame(gene = colnames(AGG.exp), corr = rep(1, ncol(AGG.exp)), p.val = rep(1, ncol(AGG.exp)))

for (i in c(1:ncol(AGG.exp))) {
  cor.res <- cor.test(AGG.exp$NOX4, AGG.exp[, i])
  AGG.cor.table.NOX4$corr[i] <- cor.res[["estimate"]][["cor"]]
  AGG.cor.table.NOX4$p.val[i] <- cor.res[["p.value"]]
}
AGG.cor.table.NOX4$relation <- ifelse(AGG.cor.table.NOX4$corr > 0.3 & AGG.cor.table.NOX4$p.val < 0.01, "Positive", 
                                      ifelse(AGG.cor.table.NOX4$corr < -0.3 & AGG.cor.table.NOX4$p.val < 0.01, "Negative", "None"))  

NOX4.pie.AGG <- data.frame(table(AGG.cor.table.NOX4$relation))
NOX4.pie.AGG$pct <- round((NOX4.pie.AGG$Freq/sum(NOX4.pie.AGG$Freq) * 100),1)
ggplot(NOX4.pie.AGG, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 1) +    
  coord_polar(theta = "y") + 
  labs(x = "", y = "", title = "") + 
  theme(axis.ticks = element_blank()) + 
  theme(legend.title = element_blank(), legend.position = "top") + 
  scale_fill_manual(values = c("#d3d4d8","#ffaaa5")) + 
  theme(axis.text.x = element_blank(),
        panel.background = element_blank(), 
        panel.grid = element_blank())




### GCLC ~ EMT genes----
EMT.exp <- na.omit(mRNA.data.tumor.wt1[c("GCLC", Hallmark$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION), ]) %>% t() %>% as.data.frame
EMT.cor.table.GCLC <- data.frame(gene = colnames(EMT.exp), corr = rep(1, ncol(EMT.exp)), p.val = rep(1, ncol(EMT.exp)))

for (i in c(1:ncol(EMT.exp))) {
  cor.res <- cor.test(EMT.exp$GCLC, EMT.exp[, i], method="spearman")
  EMT.cor.table.GCLC$corr[i] <- cor.res[["estimate"]][["rho"]]
  EMT.cor.table.GCLC$p.val[i] <- cor.res[["p.value"]]
}
EMT.cor.table.GCLC$relation <- ifelse(EMT.cor.table.GCLC$corr > 0.3 & EMT.cor.table.GCLC$p.val < 0.01, "Positive", 
                                      ifelse(EMT.cor.table.GCLC$corr < -0.3 & EMT.cor.table.GCLC$p.val < 0.01, "Negative", "None"))  

GCLC.pie.emt <- data.frame(table(EMT.cor.table.GCLC$relation))
GCLC.pie.emt$pct <- round((GCLC.pie.emt$Freq/sum(GCLC.pie.emt$Freq) * 100),1)
ggplot(GCLC.pie.emt, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 1) +    
  coord_polar(theta = "y") + 
  labs(x = "", y = "", title = "") + 
  theme(axis.ticks = element_blank()) + 
  theme(legend.title = element_blank(), legend.position = "top") + 
  scale_fill_manual(values = c("#d3d4d8","#ffaaa5")) + 
  theme(axis.text.x = element_blank(),
        panel.background = element_blank(), 
        panel.grid = element_blank())

### GCLC ~ HYPOXIA ----
HPX.exp <- na.omit(mRNA.data.tumor.wt1[c("GCLC", Hallmark$HALLMARK_HYPOXIA), ]) %>% t() %>% as.data.frame
HPX.cor.table.GCLC <- data.frame(gene = colnames(HPX.exp), corr = rep(1, ncol(HPX.exp)), p.val = rep(1, ncol(HPX.exp)))

for (i in c(1:ncol(HPX.exp))) {
  cor.res <- cor.test(HPX.exp$GCLC, HPX.exp[, i])
  HPX.cor.table.GCLC$corr[i] <- cor.res[["estimate"]][["cor"]]
  HPX.cor.table.GCLC$p.val[i] <- cor.res[["p.value"]]
}
HPX.cor.table.GCLC$relation <- ifelse(HPX.cor.table.GCLC$corr > 0.3 & HPX.cor.table.GCLC$p.val < 0.01, "Positive", 
                                      ifelse(HPX.cor.table.GCLC$corr < -0.3 & HPX.cor.table.GCLC$p.val < 0.01, "Negative", "None"))  

GCLC.pie.HPX <- data.frame(table(HPX.cor.table.GCLC$relation))
GCLC.pie.HPX$pct <- round((GCLC.pie.HPX$Freq/sum(GCLC.pie.HPX$Freq) * 100),1)
ggplot(GCLC.pie.HPX, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 1) +    
  coord_polar(theta = "y") + 
  labs(x = "", y = "", title = "") + 
  theme(axis.ticks = element_blank()) + 
  theme(legend.title = element_blank(), legend.position = "top") + 
  scale_fill_manual(values = c("#d3d4d8","#ffaaa5")) + 
  theme(axis.text.x = element_blank(),
        panel.background = element_blank(), 
        panel.grid = element_blank())

### GCLC ~ GLYCOLYSIS ----
GLY.exp <- na.omit(mRNA.data.tumor.wt1[c("GCLC", Hallmark$HALLMARK_GLYCOLYSIS), ]) %>% t() %>% as.data.frame
GLY.cor.table.GCLC <- data.frame(gene = colnames(GLY.exp), corr = rep(1, ncol(GLY.exp)), p.val = rep(1, ncol(GLY.exp)))

for (i in c(1:ncol(GLY.exp))) {
  cor.res <- cor.test(GLY.exp$GCLC, GLY.exp[, i])
  GLY.cor.table.GCLC$corr[i] <- cor.res[["estimate"]][["cor"]]
  GLY.cor.table.GCLC$p.val[i] <- cor.res[["p.value"]]
}
GLY.cor.table.GCLC$relation <- ifelse(GLY.cor.table.GCLC$corr > 0.3 & GLY.cor.table.GCLC$p.val < 0.01, "Positive", 
                                      ifelse(GLY.cor.table.GCLC$corr < -0.3 & GLY.cor.table.GCLC$p.val < 0.01, "Negative", "None"))  

GCLC.pie.GLY <- data.frame(table(GLY.cor.table.GCLC$relation))
GCLC.pie.GLY$pct <- round((GCLC.pie.GLY$Freq/sum(GCLC.pie.GLY$Freq) * 100),1)
ggplot(GCLC.pie.GLY, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 1) +    
  coord_polar(theta = "y") + 
  labs(x = "", y = "", title = "") + 
  theme(axis.ticks = element_blank()) + 
  theme(legend.title = element_blank(), legend.position = "top") + 
  scale_fill_manual(values = c("#d3d4d8","#ffaaa5")) + 
  theme(axis.text.x = element_blank(),
        panel.background = element_blank(), 
        panel.grid = element_blank())

### GCLC ~ MTOC1 ----
MTOC.exp <- na.omit(mRNA.data.tumor.wt1[c("GCLC", Hallmark$HALLMARK_MTORC1_SIGNALING), ]) %>% t() %>% as.data.frame
MTOC.cor.table.GCLC <- data.frame(gene = colnames(MTOC.exp), corr = rep(1, ncol(MTOC.exp)), p.val = rep(1, ncol(MTOC.exp)))

for (i in c(1:ncol(MTOC.exp))) {
  cor.res <- cor.test(MTOC.exp$GCLC, MTOC.exp[, i])
  MTOC.cor.table.GCLC$corr[i] <- cor.res[["estimate"]][["cor"]]
  MTOC.cor.table.GCLC$p.val[i] <- cor.res[["p.value"]]
}
MTOC.cor.table.GCLC$relation <- ifelse(MTOC.cor.table.GCLC$corr > 0.3 & MTOC.cor.table.GCLC$p.val < 0.01, "Positive", 
                                       ifelse(MTOC.cor.table.GCLC$corr < -0.3 & MTOC.cor.table.GCLC$p.val < 0.01, "Negative", "None"))  

GCLC.pie.MTOC <- data.frame(table(MTOC.cor.table.GCLC$relation))
GCLC.pie.MTOC$pct <- round((GCLC.pie.MTOC$Freq/sum(GCLC.pie.MTOC$Freq) * 100),1)
ggplot(GCLC.pie.MTOC, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 1) +    
  coord_polar(theta = "y") + 
  labs(x = "", y = "", title = "") + 
  theme(axis.ticks = element_blank()) + 
  theme(legend.title = element_blank(), legend.position = "top") + 
  scale_fill_manual(values = c("#d3d4d8","#ffaaa5")) + 
  theme(axis.text.x = element_blank(),
        panel.background = element_blank(), 
        panel.grid = element_blank())

### GCLC ~ ANGIO ----
AGG.exp <- na.omit(mRNA.data.tumor.wt1[c("GCLC", Hallmark$HALLMARK_ANGIOGENESIS), ]) %>% t() %>% as.data.frame
AGG.cor.table.GCLC <- data.frame(gene = colnames(AGG.exp), corr = rep(1, ncol(AGG.exp)), p.val = rep(1, ncol(AGG.exp)))

for (i in c(1:ncol(AGG.exp))) {
  cor.res <- cor.test(AGG.exp$GCLC, AGG.exp[, i])
  AGG.cor.table.GCLC$corr[i] <- cor.res[["estimate"]][["cor"]]
  AGG.cor.table.GCLC$p.val[i] <- cor.res[["p.value"]]
}
AGG.cor.table.GCLC$relation <- ifelse(AGG.cor.table.GCLC$corr > 0.3 & AGG.cor.table.GCLC$p.val < 0.01, "Positive", 
                                      ifelse(AGG.cor.table.GCLC$corr < -0.3 & AGG.cor.table.GCLC$p.val < 0.01, "Negative", "None"))  

GCLC.pie.AGG <- data.frame(table(AGG.cor.table.GCLC$relation))
GCLC.pie.AGG$pct <- round((GCLC.pie.AGG$Freq/sum(GCLC.pie.AGG$Freq) * 100),1)
ggplot(GCLC.pie.AGG, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 1) +    
  coord_polar(theta = "y") + 
  labs(x = "", y = "", title = "") + 
  theme(axis.ticks = element_blank()) + 
  theme(legend.title = element_blank(), legend.position = "top") + 
  scale_fill_manual(values = c("#d3d4d8","#ffaaa5")) + 
  theme(axis.text.x = element_blank(),
        panel.background = element_blank(), 
        panel.grid = element_blank())


### preprocess GDSC_AUC and GDSC_cell_line_exp ----
LUAD.GDSC.AUC <- read.table(file = "./GDSC_AUC.txt", header = T, sep = "\t")
LUAD.GDSC.AUC <- GDSC_AUC[GDSC_AUC$Cell.line.cosmic.identifiers %in% colnames(cell.line.exp.luad), ]
rownames(LUAD.GDSC.AUC) <- c(1:nrow(LUAD.GDSC.AUC))
LUAD.GDSC.AUC <- LUAD.GDSC.AUC[order(LUAD.GDSC.AUC$Cell.line.cosmic.identifiers), ]
LUAD.GDSC.AUC$Cell.line.cosmic.identifiers <- as.character(LUAD.GDSC.AUC$Cell.line.cosmic.identifiers)
LUAD.GDSC.AUC.ANA <- LUAD.GDSC.AUC[, -2]
rownames(LUAD.GDSC.AUC.ANA) <- LUAD.GDSC.AUC.ANA$Cell.line.cosmic.identifiers
LUAD.GDSC.AUC.ANA <- LUAD.GDSC.AUC.ANA[, -1]


cell.line.exp <- read.table(file = "./GDSC/Cell_line_RMA_proc_basalExp.txt", header = T, sep = "\t", quote = "")
colnames(cell.line.exp) <- gsub(pattern = "DATA.", replacement = "", colnames(cell.line.exp))
cell.line.exp <- as.data.frame(avereps(cell.line.exp[, -1], ID = cell.line.exp$GENE_SYMBOLS))
COSMIC_cellline_LUAD <- read.xlsx("./COSMIC_CELL_LINE_NAME.xlsx") %>% filter(TCGA == "LUAD")
cell.line.exp.luad <- cell.line.exp[, colnames(cell.line.exp) %in% COSMIC_cellline_LUAD$COSMIC_ID]
cell.line.exp.luad <- as.data.frame(t(cell.line.exp.luad))
cell.line.exp.luad1 <- cell.line.exp.luad[LUAD.GDSC.AUC$Cell.line.cosmic.identifiers, ]

### EMT ~ AUC ----
EMT.pos.gene <- data.frame(Hallmark = "EMT", gene = EMT.gene)
HPX.pos.gene <- data.frame(Hallmark = "HPX", gene = HPX.gene)
GLY.pos.gene <- data.frame(Hallmark = "GLY", gene = GLY.gene)
MTOC.pos.gene <- data.frame(Hallmark = "MTOC", gene = MTOC.gene)
AGG.pos.gene <- data.frame(Hallmark = "AGG", gene = AGG.gene)

hallmark.pos.gene <- rbind.data.frame(EMT.pos.gene, HPX.pos.gene, GLY.pos.gene, MTOC.pos.gene, AGG.pos.gene)
HM.cl <- cell.line.exp.luad1[, which(colnames(cell.line.exp.luad1) %in% hallmark.pos.gene$gene)]
HM.drug.cor <- data.frame(HM.gene = c(rep(colnames(HM.cl), each = ncol(LUAD.GDSC.AUC.ANA))), 
                          Drug = c(rep(colnames(LUAD.GDSC.AUC.ANA), time = ncol(HM.cl))), 
                          p.value = c(rep(0, time = ncol(HM.cl))), 
                          cor = c(rep(0, time = ncol(HM.cl))))


for (i in 1:nrow(HM.drug.cor)) {
  cor.res <- cor.test(HM.cl[, HM.drug.cor$HM.gene[i]], LUAD.GDSC.AUC.ANA[, HM.drug.cor$Drug[i]], method = "spearman")
  HM.drug.cor$cor[i] <- cor.res[["estimate"]][["rho"]]
  HM.drug.cor$p.value[i] <- cor.res[["p.value"]]
}
HM.drug.cor$Response <- ifelse(HM.drug.cor$cor > 0.3 & HM.drug.cor$p.value < 0.05, "Sensitive", 
                               ifelse(HM.drug.cor$cor < 0.3 & HM.drug.cor$p.value < 0.05, "Resistance", "None"))
HM.drug.cor <- HM.drug.cor[HM.drug.cor$Response != "None", ]
table(HM.drug.cor$Response)
colnames(hallmark.pos.gene)[2] <- "HM.gene"
HM.drug.cor <- left_join(HM.drug.cor, hallmark.pos.gene, by="HM.gene")

drug_path <- read.xlsx("./drug_path.xlsx")
drug_path$Name <- gsub(pattern = "-", replacement = "\\.", drug_path$Name)
colnames(drug_path)[1] <- "Drug"

HM.drug.cor <- left_join(HM.drug.cor, drug_path, by="Drug")
HM.drug.cor$x <- c("1")
HM.drug.cor$gene <- factor(HM.drug.cor$HM.gene, levels = unique(HM.drug.cor$HM.gene))
ggplot(HM.drug.cor, aes(x=x,y=gene,fill=Hallmark))+
  geom_tile()
  theme_void()

AGG.plot <- HM.drug.cor %>% filter(Hallmark == "AGG")
ggplot(AGG.plot, aes(x=x,y=gene,fill=Hallmark))+
  geom_tile()

mmm <- as.data.frame(table(AGG.plot$HM.gene, AGG.plot$Response))
mmm <- mmm[order(mmm$Var2, mmm$Freq, decreasing = T), ]
mmm$Freq <- ifelse(mmm$Var2 == "Resistance", -mmm$Freq, mmm$Freq)
mmm$Var1 <- factor(mmm$Var1, levels = unique(mmm$Var1))

ggplot(mmm, aes( x = Var1, weight = Freq, fill = Var2))+
  geom_bar( position = "stack")+
  scale_fill_manual(values = c("#fcbad3", "#a8d8ea"))+
  theme_bw()+
  theme(element_blank())+
  labs(x="", y="Compound counts")+
  theme(axis.title.y = element_text(size = 10, color = "black"), 
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))+
  coord_flip()
  
  
     
nnn <- as.data.frame(table(AGG.plot$HM.gene, AGG.plot$`Targeted.process/pathway`, AGG.plot$Response))
gene_path_sum <- aggregate(Freq ~ Var1, nnn, sum)
colnames(gene_path_sum)[2] <- "Sum"
nnn <- left_join(nnn, gene_path_sum, by="Var1")
nnn$pct <- round((nnn$Freq/nnn$Sum * 100),1)
ggplot(nnn[nnn$Var1 == "APP", ], aes(x = "", y = pct, fill = Var2)) +
  geom_bar(stat = "identity", width = 1) +    
  coord_polar(theta = "y") + 
  labs(x = "", y = "", title = "") + 
  theme(axis.ticks = element_blank()) + 
  scale_fill_manual(values = MyColorCode) + 
  theme(legend.title = element_blank(), legend.position = "top") + 
  theme(axis.text.x = element_blank(),
        panel.background = element_blank(), 
        panel.grid = element_blank())


deconres  <- deconvolute(mRNA.data.tumor.wt1, method="quantiseq")

deconres.h <- deconres[, c("cell_type", riskscore_surv %>% filter(riskgroup2 == "high") %>% pull(patient))]
deconres.l <- deconres[, c(riskscore_surv %>% filter(riskgroup2 == "low") %>% pull(patient))]
deconres.hl <- cbind.data.frame(deconres.h, deconres.l)
rownames(deconres.hl) <- deconres.hl$cell_type
deconres.hl <- deconres.hl[, -1]
deconres.hl <- as.data.frame(t(deconres.hl))
deconres.hl$risk <- c(rep("high","143"), rep("low", "143"))
deconres.hl.long <- gather(deconres.hl, key = "key", value = "value", -risk)
#library(ggpubr)
ggboxplot(
  deconres.hl.long,
  x = "key",
  y = "value",
  color = "black",
  fill = "risk",
  xlab = "",
  ylab = "Cell composition",
  main = "TME Cell composition group by COX2 expression"
) +
  stat_compare_means(
    aes(group = risk),
    label = "p.signif", 
    method = "t.test",
    hide.ns = T,
    size = 4.5
  ) +
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  ))

ggplot(deconres.hl.long, aes(x=key, y=value, fill=risk))+
  geom_boxplot()


xcell.res <- xCellAnalysis(LUAD.FPKM.IMDECON)
xcell.res <- as.data.frame(t(xcell.res))
riskscore_surv <- cbind.data.frame(riskscore_surv, xcell.res)

score_imm_cor <- cor.test(riskscore_surv$risk, riskscore_surv$Adipocytes, method = "spearman")

riskscore_surv <- riskscore_surv[, -c(27:58)]
imm_str_sigs <- getGmt("../Ferro_3gene/immue_stroma.gmt")
imm_str_score <- GSVA::gsva(as.matrix(LUAD.FPKM.IMDECON), imm_str_sigs, method="ssgsea", verbose=FALSE, kcdf = "Gaussian")
imm_str_score <- as.data.frame(t(imm_str_score))
riskscore_surv <- cbind.data.frame(riskscore_surv, imm_str_score)
risk_imm <- data.frame(immune = colnames(riskscore_surv)[c(27:47)], cor = 0, p=0)
for (i in 1:nrow(risk_imm)) {
  score_imm_cor <- cor.test(riskscore_surv$risk, riskscore_surv[, 26+i], method = "spearman")
  risk_imm$cor[i] <- score_imm_cor[["estimate"]][["rho"]]
  risk_imm$p[i] <- score_imm_cor[["p.value"]]
}
write.csv(risk_imm, file = "TCGA_imm_str_cor.csv")

  
LUAD.array <- openxlsx::read.xlsx("./EGFR.N.clinical.xlsx")
LUAD.array$REDD1.mean <- ifelse(LUAD.array$REDD1 < mean(LUAD.array$REDD1), "low", "high")
LUAD.array$GCLC.mean <- ifelse(LUAD.array$GCLC < mean(LUAD.array$GCLC), "low", "high")

fit2 <- survfit(Surv(time, state)~median.group, data=LUAD.array)
ggsurvplot(fit2, data=LUAD.array, pval=TRUE,palette=c('red','green'),
           risk.table = TRUE, conf.int = F) 

res.cut <- surv_cutpoint(LUAD.array, #数据集
                         time = "time", #生存状态
                         event = "state", #生存时间
                         variables = "GCLC" #需要计算的数据列名
)
res.cat <- surv_categorize(res.cut)
fit1 <- survfit(Surv(time, state)~GCLC, data=res.cat)
ggsurvplot(fit1, data=res.cat, pval=TRUE,palette = "jco",
           risk.table = TRUE, conf.int = F) 

LUAD.array.plot <- LUAD.array[, c(1,2,5,6,8,9,10,12)]
LUAD.array.plot$x <- "x"
LUAD.array.plot <- LUAD.array.plot[order(LUAD.array.plot$score, decreasing = F), ]
LUAD.array.plot$patient <- factor(LUAD.array.plot$patient, levels = LUAD.array.plot$patient)


ggplot(LUAD.array.plot,aes(x=patient, y=x, fill=median.group))+
  geom_tile(color = "white", linewidth = 0.5)+
  coord_equal()+
  theme_void()+
  labs(fill = "FeSig")

ggplot(LUAD.array.plot,aes(x=patient, y=x, fill=as.character(state)))+
  geom_tile(color = "white", linewidth = 0.5)+
  coord_equal()+
  theme_void()+
  scale_fill_manual(values = c("#9CCF65", "#F05627"))+
  labs(fill = "Survial\nstate")

ggplot(LUAD.array.plot,aes(x=patient, y=x, fill=time))+
  geom_tile(color = "white", linewidth = 0.5)+
  coord_equal()+
  theme_void()+
  scale_fill_gradient2(low = "#BABEBF", mid = "#9EB4D3", high = "#FF9687")+
  labs(fill = "Survial\ntime (Month)")

ggplot(LUAD.array.plot,aes(x=patient, y=x, fill=score))+
  geom_tile(color = "white", linewidth = 0.5)+
  coord_equal()+
  theme_void()+
  scale_fill_gradient(low = "#4B81BF",high = "red")+
  labs(fill = "FeSig score")

ggplot(LUAD.array.plot,aes(x=patient, y=x, fill=GCLC))+
  geom_tile(color = "white", linewidth = 0.5)+
  coord_equal()+
  theme_void()+
  scale_fill_gradient(low = "#F2E9DB",high = "#6667AB")+
  labs(fill = "GCLC")

ggplot(LUAD.array.plot,aes(x=patient, y=x, fill=REDD1))+
  geom_tile(color = "white", linewidth = 0.5)+
  coord_equal()+
  theme_void()+
  scale_fill_gradient(low = "#F2E9DB",high = "#F67866")+
  labs(fill = "DDIT4")


ggplot(riskscore_surv, aes(x=riskgroup2, y=CD8T, fill=riskgroup2))+
  geom_boxplot()+
  stat_compare_means(method = "t.test")+
  theme_classic()+
  labs(x="", y="T cell infiltration score")+
  theme(axis.text.x = element_blank(),
        axis.text = element_text(size = 12, colour = "black"),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size =15), 
        legend.position = "none")

fibro.surv <- riskscore_surv[, c(13,14,20,21,44)]
fibro.surv$fib.group <- ifelse(fibro.surv$Fibroblasts < median(fibro.surv$Fibroblasts), "low", "high")
fit.fib <- survfit(Surv(time, event)~fib.group, data=fibro.surv)
ggsurvplot(fit.fib, data=fibro.surv, pval=TRUE,palette=c('red','green'),
           risk.table = TRUE, conf.int = F) 

res.cut <- surv_cutpoint(fibro.surv, #数据集
                         time = "time", #生存状态
                         event = "event", #生存时间
                         variables = "Fibroblasts" #需要计算的数据列名
)
res.cat <- surv_categorize(res.cut)
fit1 <- survfit(Surv(time, event)~Fibroblasts, data=res.cat)
ggsurvplot(fit1, data=res.cat, pval=TRUE,palette = "jco",
           risk.table = TRUE, conf.int = F) 

ggplot(riskscore_surv, aes(x=risk, y=Fibroblasts))+
  geom_point()+
  geom_smooth(method = "lm", formula = y~x)

### INLGF score ----
inlgf.set <- list(inlgf = c("CCN2","IGFBPL1","IGF1R","IGF2R","IGFALS","IGFBP1","IGFBP2","IGFBP3","IGFBP4","IGFBP5","IGFBP6","IGFBP7","INSR","ITGA6","ITGAV","ITGB3","ITGB4","LRP2","KAZALD1"))
inlgf.ssgsea <- gsva(as.matrix(mRNA.data.tumor.wt1), inlgf.set, method="ssgsea", verbose=FALSE, kcdf = "Gaussian")
inlgf.ssgsea <- as.data.frame(t(inlgf.ssgsea))
riskscore_surv$inlgf <- inlgf.ssgsea$inlgf

ggplot(riskscore_surv, aes(x=riskgroup2, y=inlgf, fill=riskgroup2))+
  geom_boxplot()+
  theme_classic()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"))+
  labs(x="",y="IGF binding score")+
  ggpubr::stat_compare_means()


ins <- LUAD.FPKM.IMDECON[inlgf.set$inlgf,] %>% na.omit()
ins.score <- colMeans(ins)
riskscore_surv$inlgf <- ins.score


LN <- riskscore_surv[, c(8,21)]
LN <- LN[LN$N != "NX", ]
LN$LN2 <- ifelse(LN$N == "N0", "Neg", "Pos")

LN.plot <- as.data.frame(table(LN$riskgroup2, LN$LN2))
colnames(LN.plot) <- c("Group", "LN", "Num")
ggplot(LN.plot, aes(x=LN, y=Num, fill=Group))+
  geom_bar(stat = "identity", position = "fill")

Ferro_Niche.ssgsea <- GSVA::gsva(as.matrix(mRNA.data.tumor.wt1), Ferro.niche.list, method="ssgsea", verbose=FALSE, kcdf = "Gaussian")
Ferro_Niche.ssgsea <- as.data.frame(t(Ferro_Niche.ssgsea))
cor.test(Ferro_Niche.ssgsea$Ferro_gsea, Ferro_Niche.ssgsea$niche4_L, method = "spearman")
cor.test(Ferro_Niche.ssgsea$Ferro_gsea, Ferro_Niche.ssgsea$niche4_T, method = "spearman")
cor.test(Ferro_Niche.ssgsea$Ferro_gsea, Ferro_Niche.ssgsea$niche4_LR, method = "spearman")
cor.test(Ferro_Niche.ssgsea$Ferro_gsea, Ferro_Niche.ssgsea$niche4_LRT, method = "spearman")



saveRDS(Ferro_Niche, file = "Ferro.niche.list.rds")

HL <- cbind.data.frame(high.exp, low.exp)
write.table(HL, file = "./Ferro_3gene/TCGA_LUAD_HL_EXP.txt", quote = F, sep = "\t")

cor.test(Ferro_Niche.ssgsea$Ferro_gsea, Ferro_Niche.ssgsea$niche2_L, method = "pearson")


ggplot(Ferro_Niche.ssgsea, aes(x=Ferro_gsea, y=niche4_LR))+
  geom_point(color="grey")+
  geom_smooth(method = "lm", se=T, color="red", formula = y ~ x)+
  theme_bw()+
  xlab("Ferroptosis score")+
  ylab("Nichenet ligand score")+
  theme(axis.text = element_text(color = "black", size=13),
        axis.title = element_text(color = "black", size=16))

LUAD.surv.all <- read.table("./LUAD_surv.txt", header = T, sep = "\t", quote = "")
LUAD.surv.all <- LUAD.surv.all[,c(1:4)]
LUAD.surv.all <- LUAD.surv.all[grep("-01", LUAD.surv.all$sample), ]
LUAD.surv.all$sample <- paste(LUAD.surv.all$sample, "A", sep = "")
LUAD.surv.mut <- LUAD.surv.all[LUAD.surv.all$sample %in% colnames(mRNA.data.tumor.mut), ]

riskscore.ssgsea.mut <- gsva(as.matrix(mRNA.data.tumor.mut), risklist, method="gsva", verbose=FALSE, kcdf = "Gaussian")
riskscore.ssgsea.mut <- as.data.frame(t(riskscore.ssgsea.mut))

riskscore.ssgsea.mut <- riskscore.ssgsea.mut %>% rownames_to_column(var = "sample")
riskscore_surv.mut <- merge.data.frame(LUAD.surv.mut, riskscore.ssgsea.mut, by = "sample")

riskscore_surv.mut$riskgroup <- ifelse(riskscore_surv.mut$risk <= median(riskscore_surv.mut$risk), "low", "high")

fit.mut <- survfit(Surv(OS.time, OS)~riskgroup, data=riskscore_surv.mut)
ggsurvplot(fit.mut, data=riskscore_surv.mut, pval=TRUE,palette = "jco",
           risk.table = TRUE, conf.int = F) 


res.cut.mut <- surv_cutpoint(riskscore_surv.mut, #数据集
                             time = "OS.time", #生存状态
                             event = "OS", #生存时间
                             variables = "risk" #需要计算的数据列名
)
res.cat.mut <- surv_categorize(res.cut.mut)
fit1.mut <- survfit(Surv(OS.time, OS)~risk, data=res.cat.mut)
ggsurvplot(fit1.mut, data=res.cat.mut, pval=TRUE,palette = "jco",
           risk.table = TRUE, conf.int = F) 

normal.3.gene <- mRNA.data.normal[c("DDIT4","NOX4","GCLC"), ]
normal.3.gene <- as.data.frame(t(normal.3.gene)) %>% mutate(group = "Normal")
MUT.3.gene <- mRNA.data.tumor.mut[c("DDIT4","NOX4","GCLC"), ]
MUT.3.gene <- as.data.frame(t(MUT.3.gene)) %>% mutate(group = "Mutant")

Norm_Mut <- rbind.data.frame(normal.3.gene, MUT.3.gene)
ggplot(Norm_Mut, aes(x=group, y=GCLC))+
  geom_boxplot()+
  stat_compare_means()

MUT.3.gene <- MUT.3.gene %>% rownames_to_column(var = "sample")
MUT.3.gene.surv <- merge.data.frame(MUT.3.gene, LUAD.surv.mut, by = "sample")

library(timeROC)
library(survival)

MUT.3.gene.surv <- merge.data.frame(MUT.3.gene.surv, riskscore.ssgsea.mut, by = "sample")
time_roc_res <- timeROC(
  T = MUT.3.gene.surv$OS.time,
  delta = MUT.3.gene.surv$OS,
  marker = MUT.3.gene.surv$risk,
  cause = 1,
  weighting="marginal",
  times = c(1 * 365, 3 * 365, 5 * 365),
  ROC = TRUE,
  iid = TRUE
)
time_roc_res$AUC

plot(time_roc_res, time=1 * 12, col = "red", title = FALSE)  
plot(time_roc_res, time=3 * 12, add=TRUE, col="blue") 
plot(time_roc_res, time=5 * 12, add=TRUE, col="green") 
plot(time_roc_res, time=10 * 12, add=TRUE, col="orange") 

legend("bottomright",c("1 Years" ,"3 Years", "5 Years", "10 Years"),
       col=c("red", "blue", "green", "orange"), lty=1, lwd=2)

dat2 = data.frame(fpr = as.numeric(time_roc_res$FP),
                  tpr = as.numeric(time_roc_res$TP),
                  time = rep(as.factor(c(1 * 365, 3 * 365, 5 * 365)),each = nrow(time_roc_res$TP)))
ggplot() + 
  geom_line(data = dat2,aes(x = fpr, y = tpr,color = time),size = 1) + 
  scale_color_manual(name = NULL,values = c("#CC6633", "#9966CC", "#66CC00", "orange"),
                     labels = paste0("AUC of ",c(1,3,5,10),"-y survival: ",
                                     format(round(time_roc_res$AUC,2),nsmall = 2)))+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
        legend.position = c(0.765,0.125))+
  scale_x_continuous(expand = c(0.005,0.005))+
  scale_y_continuous(expand = c(0.005,0.005))+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  coord_fixed()
