library("autoReg")
library(survminer)
library(survival)
library(forestplot)

setwd("E:/")#设置工作目录

#读入
risk=read.table("risk.txt", header=T, sep="\t", check.names=F, row.names=1)#读取风险文件
cli=read.table("clinical.txt", header=T, sep="\t", check.names=F, row.names=1)#读取临床文件

#合并
sameSample=intersect(row.names(cli),row.names(risk))
risk=risk[sameSample,]
cli=cli[sameSample,]
rt=cbind(time=risk[,1], state=risk[,2], cli, riskScore=risk[,(ncol(risk)-1)])
colnames(rt)

#设置为数值或分类变量,注意，这里不用列名的原因，
#是因为，T，不仅仅是字符的T分期，也是逻辑的true
rt[,3] = as.numeric(rt[,3])
rt[,4] = factor(rt[,4])
rt[,5] = factor(rt[,5])
rt[,6] = factor(rt[,6])
rt[,7] = factor(rt[,7])
rt[,8] = as.numeric(rt[,8])

#删除空值及NA
rt=rt[apply(rt,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
#paste(colnames(rt), collapse="+")

#cox回归模型构建
coxmod<-coxph(Surv(time, state)~.,
              data=rt)
summary(coxmod)

#可将threshold设置为0.05，既排除单因素大于0.05的
ft3 = autoReg(coxmod,uni=TRUE,threshold=0.05) 
myft(ft3)
#导出
write.table(ft3, file="result.txt", sep="\t", row.names=F, quote=F)

#paste(colnames(rt), collapse="+")
#纳入p小于0.05临床特征
multiCox<-coxph(Surv(time, state)~ Age+T+N+Stage+riskScore, #此处修改
              data=rt)
multiCoxSum = summary(multiCox)
multiTab=data.frame()
multiTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
multiTab=cbind(id=row.names(multiTab),multiTab)
write.table(multiTab, file="mul.cox.result.txt", sep="\t", row.names=F, quote=F)

#绘制森林图
pdf(file="forest.pdf", width=10, height=8, onefile = FALSE)
ggforest(multiCox, data = rt,
         main = "Hazard ratio",
         cpositions = c(0.02, 0.22, 0.4), 
         fontsize = 0.8, 
         refLabel = "reference", 
         noDigits = 3)
dev.off()
