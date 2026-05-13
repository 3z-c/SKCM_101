library(Mime1)

##################################################################################独立验证
# 1. 读取训练数据
cat("1. 读取训练数据...\n")
input_data <- read.table("input(log2).txt", header = TRUE, sep = "\t")
cat("训练数据原始维度:", dim(input_data), "\n")

library(data.table)
# 2. 读取验证数据
cat("2. 读取验证数据...\n")
geo_data <- fread("GSE65904_with_clinical.txt", header = TRUE, sep = "\t", nThread = 4)
cat("验证数据原始维度:", dim(geo_data), "\n")
geo_data <- as.data.frame(geo_data)

# 3. 统一列名
cat("3. 统一列名格式...\n")
colnames(input_data)[1:3] <- c("ID", "OS.time", "OS")
colnames(geo_data)[1:3] <- c("ID", "OS.time", "OS")

# 4. 检查并处理生存时间零值
cat("4. 检查并处理生存时间零值...\n")

# 训练集
if (any(input_data$OS.time == 0, na.rm = TRUE)) {
  warning("训练集生存时间中存在零值，请检查数据。")
  zero_time_rows <- which(input_data$OS.time == 0)
  cat("训练集零值的行数:", length(zero_time_rows), "\n")
  cat("前10个零值行索引:", head(zero_time_rows, 10), "\n")
  # 将零值替换为一个非常小的正数
  input_data$OS.time[input_data$OS.time == 0] <- 0.5
  cat("已将所有零值替换为0.5\n")
} else {
  cat("训练集生存时间无零值\n")
}

# 验证集
if (any(geo_data$OS.time == 0, na.rm = TRUE)) {
  warning("验证集生存时间中存在零值，请检查数据。")
  zero_time_rows <- which(geo_data$OS.time == 0)
  cat("验证集零值的行数:", length(zero_time_rows), "\n")
  cat("前10个零值行索引:", head(zero_time_rows, 10), "\n")
  # 将零值替换为一个非常小的正数
  geo_data$OS.time[geo_data$OS.time == 0] <- 0.5
  cat("已将所有零值替换为0.5\n")
} else {
  cat("验证集生存时间无零值\n")
}

# 5. 获取基因列表
cat("5. 获取基因列表...\n")
train_genes <- colnames(input_data)[4:ncol(input_data)]  # 训练集基因
val_genes <- colnames(geo_data)[4:ncol(geo_data)]        # 验证集基因

cat("训练集基因数:", length(train_genes), "\n")
cat("验证集基因数:", length(val_genes), "\n")

# 6. 找到共同基因
cat("6. 寻找共同基因...\n")
common_genes <- intersect(train_genes, val_genes)
cat("共同基因数:", length(common_genes), "\n")
cat("基因重叠比例:", round(length(common_genes)/length(train_genes)*100, 1), "%\n")

# 7. 检查共同基因数量
if(length(common_genes) < 10) {
  stop("错误：共同基因太少（<10），无法进行有效分析！")
}

# 8. 只保留共同基因
cat("7. 筛选共同基因...\n")
train_final <- input_data[, c("ID", "OS.time", "OS", common_genes)]
val_final <- geo_data[, c("ID", "OS.time", "OS", common_genes)]

cat("训练集（筛选后）维度:", dim(train_final), "\n")
cat("验证集（筛选后）维度:", dim(val_final), "\n")

# 9. 创建数据列表
cat("8. 创建数据列表...\n")
list_train_vali_Data <- list(
  Dataset1 = train_final,  # 训练集
  Dataset2 = val_final     # 验证集
)

# 10. 读取并筛选候选基因
cat("9. 读取候选基因列表...\n")
genelist <- read.table("genelist.txt", header = FALSE, stringsAsFactors = FALSE)
genelist <- as.character(genelist$V1)
cat("原始候选基因数:", length(genelist), "\n")

# 筛选在共同基因中的候选基因
genelist_in_common <- intersect(genelist, common_genes)
cat("在共同基因中的候选基因数:", length(genelist_in_common), "\n")

if(length(genelist_in_common) < length(genelist)) {
  missing_genes <- setdiff(genelist, common_genes)
  cat("缺失的基因数:", length(missing_genes), "\n")
  cat("前10个缺失基因:", head(missing_genes, 10), "\n")
  # 使用共同基因中的候选基因
  genelist <- genelist_in_common
}

# 检查候选基因数量
if(length(genelist) < 5) {
  warning(paste("候选基因太少（", length(genelist), "个），模型性能可能受影响"))
} else {
  cat("最终使用的候选基因数:", length(genelist), "\n")
}

# 11. 查看数据
cat("\n=== 数据预览 ===\n")
cat("训练集前5行前6列:\n")
print(list_train_vali_Data[["Dataset1"]][1:5, 1:6])
cat("\n验证集前5行前6列:\n")
print(list_train_vali_Data[["Dataset2"]][1:5, 1:6])

#============ 第二部分：模型构建 ============
######################################构建预后模型
library(Matrix)
library(tidyr)
library(randomForestSRC)
library(purrr)
library(prodlim)
library(nlme)
library(dplyr)
library(MASS)
library(lattice)
library(mixOmics)
library(Mime1)
library(mime)

library(snow)
library(snowfall)
summary(train_data)
colSums(is.na(train_data))
sfInit(parallel = TRUE, cpus = 14, type = "SOCK")
#模型比较
res <- ML.Dev.Prog.Sig(train_data = list_train_vali_Data$Dataset1,
                       list_train_vali_Data = list_train_vali_Data,
                       unicox.filter.for.candi = T,
                       unicox_p_cutoff = 0.05,
                       candidate_genes = genelist,
                       mode = 'all', nodesize = 5, seed = 42)
# 清理资源
sfStop()

#=========================================================================================

# 只修改C-index的格式化，保留平均值格式化不变
func_text <- deparse(body(cindex_dis_all))

# 找到包含cindex$Cindex的格式化行
# 将cindex$Cindex的格式化从%.2f改为%.3f
func_text <- gsub('cindex\\$Cindex <- as\\.numeric\\(sprintf\\("%.2f", cindex\\$Cindex\\)\\)',
                  'cindex$Cindex <- as.numeric(sprintf("%.3f", cindex$Cindex))', 
                  func_text)

body(cindex_dis_all) <- parse(text = func_text)

# 运行函数
cindex_dis_all(res, 
               validate_set = names(list_train_vali_Data)[-1], 
               order = names(list_train_vali_Data), 
               width = 0.35)

#=========================================================================================
selected_model <- "CoxBoost + survival-SVM"  

#最优模型的C-index指数
cindex_dis_select(res,
                  model = selected_model,
                  order = names(list_train_vali_Data))

##############################################################################计算风险评分
str(res$riskscore[[selected_model]], max.level = 3)

if (!selected_model %in% names(res$riskscore)) {
  stop(paste("模型 '", selected_model, "' 不存在。可用模型：\n", 
             paste(names(res$riskscore), collapse = ", ")))
}

model_results <- res$riskscore[[selected_model]]

# 保存训练集风险评分
write.csv(
  data.frame(
    SampleID = model_results$Dataset1$ID,
    RiskScore = model_results$Dataset1$RS,
    Model = selected_model
  ),
  paste0("RiskScores_", gsub("[^A-Za-z0-9]+", "_", selected_model), "_Training.csv"),
  row.names = FALSE,
  quote = FALSE
)

# 保存验证集风险评分
write.csv(
  data.frame(
    SampleID = model_results$Dataset2$ID,
    RiskScore = model_results$Dataset2$RS,
    Model = selected_model
  ),
  paste0("RiskScores_", gsub("[^A-Za-z0-9]+", "_", selected_model), "_Validation.csv"),
  row.names = FALSE,
  quote = FALSE
)


#最优模型的生存曲线
survplot <- vector("list", 2)
for (i in 1:2) {
  print(survplot[[i]] <- rs_sur(res, 
                                model_name = selected_model,
                                dataset = names(list_train_vali_Data)[i],
                                median.line = "hv",
                                cutoff = 0.5,
                                conf.int = T,
                                xlab = "Day",
                                pval.coord = c(1000, 0.9)))
}
aplot::plot_list(gglist = survplot, ncol = 2)


library(survivalROC)
#
all.auc.1y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["Dataset1"]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 1,
                             auc_cal_method="KM")
all.auc.3y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["Dataset1"]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 3,
                             auc_cal_method="KM")
all.auc.5y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["Dataset1"]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 5,
                             auc_cal_method="KM")
#绘制所有模型预测的 1 年 AUC
auc_dis_all(all.auc.1y,
            dataset = names(list_train_vali_Data),
            validate_set=names(list_train_vali_Data)[-1],
            order= names(list_train_vali_Data),
            width = 0.35,
            year=1)
#出图，特定模型在不同数据集中的 ROC
roc_vis(all.auc.1y,
        model_name = selected_model,
        dataset = names(list_train_vali_Data),
        order= names(list_train_vali_Data),
        anno_position=c(0.65,0.55),
        year=1)
#出图，绘制特定模型在不同数据集中的 1 、 3 和 5 年 AUC：
auc_dis_select(list(all.auc.1y,all.auc.3y,all.auc.5y),
               model_name= selected_model,
               dataset = names(list_train_vali_Data),
               order= names(list_train_vali_Data),
               year=c(1,3,5))



####引入对比
# 读取签名数据
signature_data <- read.table("signature.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# 查看签名数据的前几行
head(signature_data)

rs.glioma.lgg.gbm <- cal_RS_pre.prog.sig(use_your_own_collected_sig = T,collected_sig_table = signature_data,
                                         list_input_data = list_train_vali_Data)


cc.glioma.lgg.gbm <- cal_cindex_pre.prog.sig(use_your_own_collected_sig = T,collected_sig_table = signature_data,
                                             list_input_data = list_train_vali_Data)

cindex_comp(cc.glioma.lgg.gbm,
            res,
            model_name= selected_model,
            dataset=names(list_train_vali_Data))

auc.glioma.lgg.gbm.1 <- cal_auc_pre.prog.sig(use_your_own_collected_sig = T,
                                             collected_sig_table = signature_data,
                                             list_input_data = list_train_vali_Data,AUC_time = 1,
                                             auc_cal_method = 'KM')

auc_comp(auc.glioma.lgg.gbm.1,
         all.auc.1y,
         model_name= selected_model,
         dataset=names(list_train_vali_Data))




library(DOSE)
res.feature.all <- ML.Corefeature.Prog.Screen(InputMatrix = list_train_vali_Data$Dataset1,
                                              candidate_genes = genelist,
                                              mode = "all_without_SVM",
                                              nodesize = 5,
                                              seed = 42)

core_feature_select(res.feature.all)


core_feature_rank(res.feature.all, top=20)

dataset_col<-c("#3182BDFF","#E6550DFF")
corplot <- list()
for (i in c(1:2)) {
  print(corplot[[i]]<-cor_plot(list_train_vali_Data[[i]],
                               dataset=names(list_train_vali_Data)[i],
                               color = dataset_col[i],
                               feature1="GSDMD",
                               feature2="CD274",
                               method="pearson"))
}
aplot::plot_list(gglist=corplot,ncol=2)
survplot <- vector("list",2) 
for (i in c(1:2)) {
  print(survplot[[i]]<-core_feature_sur("GSDMD", 
                                        InputMatrix=list_train_vali_Data[[i]],
                                        dataset = names(list_train_vali_Data)[i],
                                        #color=c("blue","green"),
                                        median.line = "hv",
                                        cutoff = 0.5,
                                        conf.int = T,
                                        xlab="Day",pval.coord=c(1000,0.9)))
}
aplot::plot_list(gglist=survplot,ncol=2)

