#------------------------------------------------------------------------------#                                 
 
                    ###    Clinical Data Extraction    ###
 
#------------------------------------------------------------------------------#

#--------------------------------------------#
####        0. Install libraries          #### 
#--------------------------------------------#

library(dplyr)
library(useful)
library(janitor)
library(utils)
library(GEOquery)
library(affxparser)
library(stringr)
library(affyio)

#--------------------------------------------#
####         1. Source Functions          #### 
#--------------------------------------------#

source("/home/user/Documents/Files/Projects/Colorectal_Cancer_Superset/scripts/Metacohort_Generation/clinical_data_extraction_misc.R")


#--------------------------------------------#
####     2. Dataset Data Extraction       ####
#--------------------------------------------#


##### GSE14333. Jorissen #####

gse <- getGEO('GSE14333',GSEMatrix=TRUE)
GSE14333_RAW <- as.data.frame(t(pData(phenoData(gse[[1]]))))

GSE14333_CD <- modify_GSE14333_RAW(GSE_RAW = GSE14333_RAW)


##### GSE143985. Shinto E. #####

gse <- getGEO('GSE143985',GSEMatrix=TRUE)
GSE143985_RAW <- as.data.frame(t(pData(phenoData(gse[[1]]))))

GSE143985_CD <- modify_GSE143985_RAW(GSE_RAW = GSE143985_RAW)


##### GSE17536. Smith JJ #####

gse <- getGEO('GSE17536',GSEMatrix=TRUE)
GSE17536_RAW <- as.data.frame(t(pData(phenoData(gse[[1]]))))

GSE17536_CD <- modify_GSE17536_RAW(GSE_RAW = GSE17536_RAW)


##### GSE17537. Smith JJ. #####

gse <- getGEO('GSE17537',GSEMatrix=TRUE)
GSE17537_RAW <- as.data.frame(t(pData(phenoData(gse[[1]]))))

GSE17537_CD <- modify_GSE17537_RAW(GSE_RAW = GSE17537_RAW)


##### GSE33114. de Sousa E Melo F. #####

gse <- getGEO('GSE33114',GSEMatrix=TRUE)
GSE33114_RAW <- as.data.frame(t(pData(phenoData(gse[[1]]))))

GSE33114_CD <- modify_GSE33114_RAW(GSE_RAW = GSE33114_RAW)


##### GSE38832. Tripathi MK #####

gse <- getGEO('GSE38832',GSEMatrix=TRUE)
GSE38832_RAW <- as.data.frame(t(pData(phenoData(gse[[1]]))))

GSE38832_CD <- modify_GSE38832_RAW(GSE_RAW = GSE38832_RAW)


##### GSE39582. Marisa #####

gse <- getGEO('GSE39582',GSEMatrix=TRUE)
GSE39582_RAW <- as.data.frame(t(pData(phenoData(gse[[1]]))))

GSE39582_CD <- modify_GSE39582_RAW(GSE_RAW = GSE39582_RAW)


########### GSE37892. Laibe S ############################

gse <- getGEO('GSE37892',GSEMatrix=TRUE)
GSE37892_RAW <- as.data.frame(t(pData(phenoData(gse[[1]]))))

GSE37892_CD <- get_GSE7_CD(GSE7_RAW = GSE37892_RAW)

########### GSE92921. Gotoh K ############################

gse <- getGEO('GSE92921',GSEMatrix=TRUE)
GSE92921_RAW <- as.data.frame(t(pData(phenoData(gse[[1]]))))

GSE92921_CD <- get_GSE11_CD(GSE11_RAW = GSE92921_RAW)



#--------------------------------------------#
####       3. Pheno.data creation         ####
#--------------------------------------------#


df_list <- list("GSE14333_CD" = GSE14333_CD, 
                "GSE143985_CD" = GSE143985_CD, 
                "GSE17536_CD"=GSE17536_CD, 
                "GSE17537_CD"=GSE17537_CD,
                "GSE33114_CD"=GSE33114_CD, 
                "GSE38832_CD"=GSE38832_CD, 
                "GSE39582_CD"=GSE39582_CD)

df_list2 <- lapply(names(df_list), function(x) {
  a <- as.data.frame(t(as.matrix(df_list[[x]])))
})

names(df_list2) <- c("GSE14333_CD","GSE143985_CD","GSE17536_CD","GSE17537_CD", 
                     "GSE33114_CD","GSE38832_CD","GSE39582_CD")

saveRDS(df_list2, file = "list_of_samples.rds")

df_list2 <- readRDS("list_of_samples.rds")

df <- do.call("bind_rows", df_list2)


##### Variable Homogeneization ##### 

df <- homogenize_variables(df)


##### TXT files Creation ##### 

###### Extract dates from cel files ######

celfiles <-list.files("/home/user/Downloads/CEL_Files")
dates <- NULL

for (x in celfiles) {
  a <- get.celfile.dates(paste0("/home/user/Downloads/CEL_Files/", x))
  b <- as.character(a)
  dates <- c(dates, b)
}

d <- data.frame("date" = dates)
rownames(d) <- celfiles
d$GSM <- gsub("_.*", "", rownames(d))
d$GSM <- gsub("\\..*", "", d$GSM)

dataset <- as.data.frame(df$dataset)
rownames(dataset) <- rownames(df)
colnames(dataset) <- "dataset"

dd <- merge(d, dataset, by.x = "GSM", by.y = 0)
dd.necessary <- data.frame("GSM" = dd$GSM, "CEL_Date" = dd$date)

saveRDS(dd.necessary, file ="/home/user/Documents/Files/Projects/Colorectal_Cancer_Superset/date_GSM.rds")

###### File preparation ######

df$dfs_time <- as.numeric(df$dfs_time)
df$stage <- as.numeric(df$stage)

t.df <- as.data.frame(t(df))

elim <- df[is.na(df$dfs_time),]

del.samples <- is.na(elim$dfs_time) & is.na(elim$dss_time) & is.na(elim$os_time)
del.samples <- rownames(elim[del.samples,])

df_filt.samples <- which(colnames(t.df) %in% del.samples)
df.filt <- t.df[,-df_filt.samples]

t.df.filt <- as.data.frame(t(df.filt))
t.df.filt$dfs_time <- as.numeric(t.df.filt$dfs_time)
t.df.filt$stage <- as.integer(t.df.filt$stage)

t.df.filt <- merge(df, dd.necessary, by.x = 0, by.y = "GSM")
t.df.filt$age <- as.numeric(t.df.filt$age)

usable <- which(!is.na(t.df.filt$dfs_time) & !is.na(t.df.filt$dfs_event))
usable.t.df.filt <- t.df.filt[usable,]
rownames(usable.t.df.filt) <- usable.t.df.filt$Row.names
usable.t.df.filt$Row.names <- NULL

#Reorder columns to aid comprehension
usable.t.df.filt <- usable.t.df.filt[,c(1,20,21,2:19)]


###### Duplicated Samples ###### 

pheno.data <- usable.t.df.filt

pheno.data$Row.names <- rownames(pheno.data)

usable.t.df.filt.nodups <- pheno.data[!pheno.data$first_author == "Gotoh K",]

dups <- usable.t.df.filt.nodups %>% 
  group_by(age, gender, stage, dfs_event, dfs_time, CEL_Date) %>% 
  filter(n() >= 2) %>% 
  ungroup()

nodups <- usable.t.df.filt.nodups %>% 
  group_by(dataset, age, gender, stage, dfs_event, dfs_time, CEL_Date) %>% 
  filter(n() >= 2) %>% 
  ungroup()


dups.elim <- dups[duplicated(dups[,c(3,4,5,6,8,9)]),]
dups.noelim <- dups[!duplicated(dups[,c(3,4,5,6,8,9)]),]

usable.t.df.filt.nodups <- usable.t.df.filt.nodups[!usable.t.df.filt.nodups$Row.names %in% dups.elim$Row.names,]
usable.t.df.filt.nodups$Row.names <- NULL

saveRDS(usable.t.df.filt.nodups, file ="/home/user/Documents/Files/Projects/Colorectal_Cancer_Superset/Results_Files/pheno_data_july2024.rds")

###### SUMMARY ######

summary_table <- usable.t.df.filt %>%                                        
  group_by(dataset)%>%                                             
  summarise(                                                         
    samples = n(),
    Gender_M = paste0(sum(gender == "M"), " (", round((sum(gender == "M", na.rm=TRUE)/sum(!is.na(gender), na.rm=TRUE))*100, 2), "%)"),
    braf_info = sum(!is.na(braf_mut)),
    kras_info = sum(!is.na(kras_mut)),
    tp53_info = sum(!is.na(tp53_mut)),
    mean_dfs  = round(mean(dfs_time, na.rm=TRUE), 3),
    stage_1   = sum(stage == 1,na.rm=TRUE),
    stage_2   = sum(stage == 2,na.rm=TRUE),
    stage_3   = sum(stage == 3,na.rm=TRUE),
    stage_4   = sum(stage == 4,na.rm=TRUE),
    stage_NA  = sum(is.na(stage)),
    time_event = sum(!is.na(dfs_time) & !is.na(dfs_event)))

write.xlsx(summary_table, "summary_draft.xlsx")

summary_table_stage <- usable.t.df.filt %>%                                        
  group_by(stage)%>%                                             
  summarise(                                                         
    samples = n(),
    braf_Yes = paste0(sum(braf_mut == "1", na.rm=TRUE), " (", round((sum(braf_mut == "1", na.rm=TRUE)/sum(!is.na(braf_mut), na.rm=TRUE))*100, 2), "%)"),
    braf_No = paste0(sum(braf_mut == "0", na.rm=TRUE), " (", round((sum(braf_mut == "0", na.rm=TRUE)/sum(!is.na(braf_mut), na.rm=TRUE))*100, 2), "%)"),
    braf_Na = sum(is.na(braf_mut)),
    kras_Yes = paste0(sum(kras_mut == "1", na.rm=TRUE), " (", round((sum(kras_mut == "1", na.rm=TRUE)/sum(!is.na(kras_mut), na.rm=TRUE))*100, 2), "%)"),
    kras_No = paste0(sum(kras_mut == "0", na.rm=TRUE), " (", round((sum(kras_mut == "0", na.rm=TRUE)/sum(!is.na(kras_mut), na.rm=TRUE))*100, 2), "%)"),
    kras_Na = sum(is.na(kras_mut)),
    tp53_Yes = paste0(sum(tp53_mut == "1", na.rm=TRUE), " (", round((sum(tp53_mut == "1", na.rm=TRUE)/sum(!is.na(tp53_mut), na.rm=TRUE))*100, 2), "%)"),
    tp53_No = paste0(sum(tp53_mut == "0", na.rm=TRUE), " (", round((sum(tp53_mut == "0", na.rm=TRUE)/sum(!is.na(tp53_mut), na.rm=TRUE))*100, 2), "%)"),
    tp53_Na = sum(is.na(tp53_mut)))

write.xlsx(summary_table_stage, "summary_stage_draft.xlsx")





