library(Hmisc)

#加载包
library(survival)
library(regplot)
library(rms)

#设置工作目录
setwd("E:/")

#读取风险输入文件
risk=read.table("risk.txt", header=T, sep="\t", check.names=F, row.names=1)

#读取临床数据文件
cli=read.table("clinical.txt", header=T, sep="\t", check.names=F, row.names=1)
#删除包含unknow的样本
cli=cli[apply(cli,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
#将年龄这一列转化为数值型
cli$Age=as.numeric(cli$Age)

#合并数据
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
colnames(risk1)
#rt=cbind(risk1[,c("time", "state", "risk")], cli)
rt=cbind(risk1[,c("time", "state", "riskScore")], cli)
#rt=cbind(risk1, cli)

#cox回归
res.cox=coxph(Surv(time, state) ~ . , data = rt)
#绘制列线图

#单项得分，即图中的Points，表示每个变量在不同取值下所对应的单项分数，
#总得分，即Total Point，表示所有变量取值后对应的单项分数加起来合计的总得分。
#预测概率

nom1=regplot(res.cox,
             clickable=F,
             title="",
             points=TRUE,
             droplines=TRUE,
             observation=rt[50,],
             rank="sd",
             failtime = c(1,3,5),
             prfail = F)

#列线图打分
nomoRisk=predict(res.cox, data=rt, type="risk")
rt=cbind(risk1, Nomogram=nomoRisk)
outTab=rbind(ID=colnames(rt), rt)
write.table(outTab, file="nomoRisk.txt", sep="\t", col.names=F, quote=F)

#校准曲线
pdf(file="calibration.pdf", width=6, height=6)
#1年
f <- cph(Surv(time, state) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=1)
cal <- calibrate(f, cmethod="KM", method="boot", u=1, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-predicted OS (%)", ylab="Observed OS (%)", lwd=3, col="Firebrick2", sub=F)
#3年
f <- cph(Surv(time, state) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=3)
cal <- calibrate(f, cmethod="KM", method="boot", u=3, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", lwd=3, col="MediumSeaGreen", sub=F, add=T)
#5年
f <- cph(Surv(time, state) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=5)
cal <- calibrate(f, cmethod="KM", method="boot", u=5, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",  lwd=3, col="NavyBlue", sub=F, add=T)
legend('bottomright', c('1-year', '3-year', '5-year'),
       col=c("Firebrick3","MediumSeaGreen","NavyBlue"), lwd=3, bty = 'n')
dev.off()

