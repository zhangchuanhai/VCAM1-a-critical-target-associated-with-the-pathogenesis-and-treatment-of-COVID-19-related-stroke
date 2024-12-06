# ----0、条件设置----
# 系统报错改为英文
Sys.setenv(LANGUAGE = "en")
# 禁止转化为因子
options(stringsAsFactors = FALSE)
# 清空环境
rm(list=ls())
# 设置工作目录
setwd("/home/data/t180324/R/R_project/MR/")

library(TwoSampleMR)
library(openxlsx)
library(data.table)

# ---1、导入循环蛋白的数据----
# ebi-a-GCST011075 新冠重症
# ebi-a-GCST011081 新冠

# ebi-a-GCST90038613 中风
# ebi-a-GCST90018864 缺血性脑卒
# ebi-a-GCST90014123 腔隙性中风

# #-----读取exposure
repeat{try({
    ## 提取暴露变量'
    exposure_data <-extract_instruments(outcomes=c('ebi-a-GCST011075', # 新冠重症
                                                   'ebi-a-GCST011081' # 新冠
                                                   ),p1 = 1e-5,
                                 clump=TRUE, r2=0.01,kb=5000,
                                 access_token = NULL) #获取暴露数据
  })
  if(exists("exposure_data")){break}
  Sys.sleep(2)
};head(exposure_data) #查看暴露数据

#-----读取outcome_data
outcome_data <-extract_outcome_data(snps=exposure_data$SNP, outcomes=c("ebi-a-GCST90038613", # 中风
                                                                       "ebi-a-GCST90018864", # 缺血性脑卒
                                                                       "ebi-a-GCST90014123")) # 腔隙性中风"
save(exposure_data,outcome_data,file="./AIS.rdata")

#-----预处理数据
dat <- harmonise_data(exposure_data,outcome_data)
# write.csv(dat,"dat_harmonised_protein.csv")
#-----自选方法进行MR分析
res <- mr(dat,method_list = c("mr_ivw","mr_ivw_mre")) # 查询方法：mr_method_list()
res
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)
library(MRPRESSO)
mr_presso(BetaOutcome = "eaf.outcome", BetaExposure = "eaf.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  SignifThreshold = 0.05)
res_loo <- mr_leaveoneout(subset(dat,id.outcome=="ebi-a-GCST90018864"&id.exposure=="ebi-a-GCST011075"))
mr_leaveoneout_plot(res_loo)
p1 <- mr_scatter_plot(subset(res,id.outcome=="ebi-a-GCST90018864"&id.exposure=="ebi-a-GCST011075"), 
                      subset(dat,id.outcome=="ebi-a-GCST90018864"&id.exposure=="ebi-a-GCST011075"))
p1
res_single <- mr_singlesnp(subset(dat,id.outcome=="ebi-a-GCST90018864"&id.exposure=="ebi-a-GCST011075"))
mr_forest_plot(res_single)
#绘制漏斗图
mr_funnel_plot(res_single)
