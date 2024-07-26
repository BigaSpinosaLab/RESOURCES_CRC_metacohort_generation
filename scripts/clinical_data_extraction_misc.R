#------------------------------------------------------------------------------#                                  
 
                    ###    Clinical Data Extraction    ###
                    ###         Misc Functions         ###
 
#------------------------------------------------------------------------------#

# All of these functions are used in the "clinical_data_extraction.R" script.

#------------------------------------------------------------------------------#

####    Dataset modifications    ####

#' These functions are used to select the needed variables and 
#' make all of the variables have the same name between datasets
#' 
#' @param GSE_RAW           Data frame with phenoData of a dataset
#' @returns                 Data frame with selected variables

  ##### GSE14333. Jorissen #####
  modify_GSE14333_RAW <- function(GSE_RAW) {
  
  GSE14333_RAW <- GSE14333_RAW[-1,]
  
  Samples <- as.vector(GSE14333_RAW[1,])
  
  dataset <- rep("GSE14333", ncol(GSE14333_RAW))
  
  location <- as.vector(GSE14333_RAW[9,])
  location <- sapply(location, function(x) unlist(strsplit(x, split= ";"))[1])
  location <- unlist(gsub(".*: ", "", location))
  
  stage <- as.vector(GSE14333_RAW[9,])
  stage <- sapply(stage, function(x) unlist(strsplit(x, split= ";"))[2])
  stage <- unlist(gsub(".*: ", "", stage))
  
  age <- as.vector(GSE14333_RAW[9,])
  age <- sapply(age, function(x) unlist(strsplit(x, split= ";"))[3])
  age <- unlist(gsub(".*: ", "", age))
  
  gender <- as.vector(GSE14333_RAW[9,])
  gender <- sapply(gender, function(x) unlist(strsplit(x, split= ";"))[4])
  gender <- unlist(gsub(".*: ", "", gender))
  
  dfs_event <- GSE14333_RAW[9,]
  dfs_event <- sapply(dfs_event, function(x) unlist(strsplit(x, split= ";"))[6])
  dfs_event <- unlist(gsub(".*: ", "", dfs_event))
  
  dfs_time <- GSE14333_RAW[9,]
  dfs_time <- sapply(dfs_time, function(x) unlist(strsplit(x, split= ";"))[5])
  dfs_time <- unlist(gsub(".*: ", "", dfs_time))
  
  radio <- GSE14333_RAW[9,]
  radio <- sapply(radio, function(x) unlist(strsplit(x, split= ";"))[7])
  radio <- unlist(gsub(".*: ", "", radio))
  
  chemo <- GSE14333_RAW[9,]
  chemo <- sapply(chemo, function(x) unlist(strsplit(x, split= ";"))[8])
  chemo <- unlist(gsub(".*: ", "", chemo))
  
  GSE14333_CD <- data.frame("dataset"=dataset,
                            "age" = age,
                            "gender"=gender,
                            "stage"=stage, 
                            "location"=location,
                            "dfs_event"=dfs_event,
                            "dfs_time"=dfs_time,
                            "radio"=radio,
                            "chemo"=chemo)
  GSE14333_CD <- t(GSE14333_CD)
  
  colnames(GSE14333_CD) <- Samples
  
  t.GSE14333_CD <- t(GSE14333_CD)
  
  t.GSE14333_CD <- as.data.frame(t.GSE14333_CD) %>%
    mutate(dfs_event = case_when(dfs_event == "0"  ~ "1",
                                 dfs_event == "1"  ~ "0"))
  
  GSE14333_CD <- as.data.frame(t(t.GSE14333_CD))
  
  return(as.data.frame(GSE14333_CD))
}


  ##### GSE143985. Shinto E.##### 
  modify_GSE143985_RAW <- function(GSE_RAW) {
  
  GSE143985_RAW <- GSE143985_RAW[-1,]
  
  Samples <- as.vector(GSE143985_RAW[1,])
  
  dataset <- rep("GSE143985", ncol(GSE143985_RAW))
  
  stage <- GSE143985_RAW[9,]
  stage <- gsub("Stage: ", "", stage)
  
  dfs_event <- GSE143985_RAW[10,]
  dfs_event <- gsub("dfs_event: ", "", dfs_event)
  
  dfs_time_days <- GSE143985_RAW[11,]
  dfs_time_days <- gsub("dfs_time: ", "", dfs_time_days)
  
  dfs_time <- as.numeric(dfs_time_days) / 30
  
  chemo <- GSE143985_RAW[12,]
  chemo <- gsub("chemotherapy: ", "", chemo)
  
  braf_mut <- GSE143985_RAW[15,]
  braf_mut <- gsub("braf_mutation: ", "", braf_mut)
  
  kras_mut <- GSE143985_RAW[16,]
  kras_mut <- gsub("kras_mutation: ", "", kras_mut)
  
  tp53_mut <- GSE143985_RAW[17,]
  tp53_mut <- gsub("tp53_mutation: ", "", tp53_mut)
  
  GSE143985_CD <- data.frame("dataset"=dataset,
                             "stage"=stage, 
                             "dfs_event"=dfs_event, 
                             "dfs_time"=dfs_time, 
                             "dfs_time_days"=dfs_time_days,
                             "chemo"=chemo,
                             "braf_mut"=braf_mut,
                             "kras_mut"=kras_mut,
                             "tp53_mut"=tp53_mut)
  
  GSE143985_CD <- t(GSE143985_CD)
  
  
  colnames(GSE143985_CD) <- Samples
  
  return(as.data.frame(GSE143985_CD))
}


  ##### GSE17536. Smith JJ. ##### 
  modify_GSE17536_RAW <- function(GSE_RAW) {
  
  GSE17536_RAW <- GSE17536_RAW[-1,]
  
  Samples <- as.vector(GSE17536_RAW[1,])
  
  dataset <- rep("GSE17536", ncol(GSE17536_RAW))
  
  age <- GSE17536_RAW[9,]
  age <- gsub("age: ", "", age)
  
  gender <- GSE17536_RAW[10,]
  gender <- gsub("gender: ", "", gender)
  
  ethnicity <- GSE17536_RAW[11,]
  ethnicity <- gsub("ethnicity: ", "", ethnicity)
  
  stage <- GSE17536_RAW[12,]
  stage <- gsub("ajcc_stage: ", "", stage)
  
  grade <- GSE17536_RAW[13,]
  grade <- gsub("grade: ", "", grade)
  
  dss_event <- GSE17536_RAW[15,]
  dss_event <- gsub(".*: ", "", dss_event)
  
  dfs_event <- GSE17536_RAW[16,]
  dfs_event <- gsub(".*: ", "", dfs_event)
  
  dss_time <- GSE17536_RAW[18,]
  dss_time <- gsub(".*: ", "", dss_time)
  
  dfs_time <- GSE17536_RAW[19,]
  dfs_time <- gsub(".*: ", "", dfs_time)
  
  
  GSE17536_CD <- data.frame("dataset"=dataset,
                            "age"=age, 
                            "gender"=gender, 
                            "ethnicity"=ethnicity, 
                            "stage"=stage,
                            "grade"=grade,
                            "dss_event"=dss_event,
                            "dfs_event"=dfs_event,
                            "dss_time"=dss_time,
                            "dfs_time"=dfs_time)
  
  GSE17536_CD <- t(GSE17536_CD)
  
  colnames(GSE17536_CD) <- Samples
  
  
  return(as.data.frame(GSE17536_CD))
}


  ##### GSE17537. Smith JJ. ##### 
  modify_GSE17537_RAW <- function(GSE_RAW) {
  
  GSE17537_RAW <- GSE17537_RAW[-1,]
  
  Samples <- as.vector(GSE17537_RAW[1,])
  
  dataset <- rep("GSE17537", ncol(GSE17537_RAW))
  
  age <- GSE17537_RAW[9,]
  age <- gsub("age: ", "", age)
  
  gender <- GSE17537_RAW[10,]
  gender <- gsub("gender: ", "", gender)
  
  ethnicity <- GSE17537_RAW[11,]
  ethnicity <- gsub("ethnicity: ", "", ethnicity)
  
  stage <- GSE17537_RAW[12,]
  stage <- gsub("ajcc_stage: ", "", stage)
  
  grade <- as.character(GSE17537_RAW[45,])
  
  os_event <- as.character(GSE17537_RAW[47,])
  
  dfs_event <- as.character(GSE17537_RAW[41,])
  
  os_time <- as.character(GSE17537_RAW[46,])
  
  dfs_time <- as.character(GSE17537_RAW[42,])
  
  
  
  
  GSE17537_CD <- data.frame("dataset"=dataset,
                            "age"=age, 
                            "gender"=gender, 
                            "ethnicity"=ethnicity, 
                            "stage"=stage,
                            "grade"=grade,
                            "os_event"=os_event,
                            "dfs_event"=dfs_event,
                            "os_time"=os_time,
                            "dfs_time"=dfs_time)
  
  GSE17537_CD <- t(GSE17537_CD)
  
  colnames(GSE17537_CD) <- Samples
  
  
  return(as.data.frame(GSE17537_CD))
}


  ##### GSE33114. de Sousa E Melo F. ##### 
  modify_GSE33114_RAW <- function(GSE_RAW) {
  
  GSE33114_RAW2 <- GSE33114_RAW[-1,]
  
  GSE33114_RAW2 <- GSE33114_RAW2[,-c(1:12)]
  
  GSE33114_RAW2 <- GSE33114_RAW2[,-c((ncol(GSE33114_RAW2)-5):ncol(GSE33114_RAW2))]
  
  Samples <- as.vector(GSE33114_RAW2[1,])
  
  dataset <- rep("GSE33114", ncol(GSE33114_RAW2))
  
  gender <- GSE33114_RAW2[12,]
  gender <- gsub(".*: ", "", gender)
  
  dfs_event <- GSE33114_RAW2[13,]
  dfs_event <- gsub(".*: ", "", dfs_event)
  
  dfs_time_days <- GSE33114_RAW2[14,]
  dfs_time_days <- gsub(".*: ", "", dfs_time_days)
  
  dfs_time <- as.numeric(dfs_time_days) / 30
  
  GSE33114_CD <- data.frame("dataset"=dataset,
                            "gender"=gender, 
                            "dfs_event"=dfs_event,
                            "dfs_time"=dfs_time,
                            "dfs_time_days"=dfs_time_days)
  
  GSE33114_CD <- GSE33114_CD %>%
    mutate(dfs_event = case_when(dfs_event == "yes" | dfs_event == "Yes" | dfs_event == "1" ~ "Y",
                                 dfs_event == "no" | dfs_event == "No" | dfs_event == "0" ~ "N"))
  
  GSE33114_CD <- t(GSE33114_CD)
  
  colnames(GSE33114_CD) <- Samples
  
  return(as.data.frame(GSE33114_CD))
}


  ##### GSE38832. Tripathi MK #####  
  modify_GSE38832_RAW <- function(GSE_RAW) {
  
  GSE38832_RAW <- GSE38832_RAW[-1,]
  
  Samples <- as.vector(GSE38832_RAW[1,])
  
  dataset <- rep("GSE38832", ncol(GSE38832_RAW))
  
  stage <- GSE38832_RAW[9,]
  stage <- gsub(".*: ", "", stage)
  
  dfs_event <- GSE38832_RAW[10,]
  dfs_event <- gsub(".*: ", "", dfs_event)
  
  dfs_time <- GSE38832_RAW[11,]
  dfs_time <- gsub(".*: ", "", dfs_time)
  
  dss_event <- GSE38832_RAW[12,]
  dss_event <- gsub(".*: ", "", dss_event)
  
  dss_time <- GSE38832_RAW[13,]
  dss_time <- gsub(".*: ", "", dss_time)
  
  
  GSE38832_CD <- data.frame("dataset"=dataset,
                            "stage"=stage, 
                            "dfs_event"=dfs_event,
                            "dfs_time"=dfs_time,
                            "dss_event"=dss_event,
                            "dss_time"=dss_time)
  # 
  # GSE7_CD <- GSE7_CD %>%
  #   mutate(dfs_event = case_when(dfs_event == "yes" | dfs_event == "Yes" ~ "recurrence",
  #                                dfs_event == "no" | dfs_event == "No" ~ "no recurrence"))
  
  GSE38832_CD <- t(GSE38832_CD)
  
  colnames(GSE38832_CD) <- Samples
  
  return(as.data.frame(GSE38832_CD))
}


  ##### GSE39582. Marisa ##### 
  modify_GSE39582_RAW <- function(GSE_RAW) {
  
  GSE39582_RAW <- GSE39582_RAW[-1,]
  
  Samples <- as.vector(GSE39582_RAW[1,])
  
  dataset <- rep("GSE39582", ncol(GSE39582_RAW))
  
  age <- GSE39582_RAW[12,]
  age <- gsub(".*: ", "", age)
  
  gender <- GSE39582_RAW[11,]
  gender <- gsub(".*: ", "", gender)
  
  stage <- GSE39582_RAW[13,]
  stage <- gsub(".*: ", "", stage)
  
  location <- GSE39582_RAW[17,]
  location <- gsub(".*: ", "", location)
  
  dfs_event <- GSE39582_RAW[20,]
  dfs_event <- gsub(".*: ", "", dfs_event)
  
  dfs_time <- GSE39582_RAW[21,]
  dfs_time <- gsub(".*: ", "", dfs_time)
  
  os_event <- GSE39582_RAW[22,]
  os_event <- gsub(".*: ", "", os_event)
  
  os_time <- GSE39582_RAW[23,]
  os_time <- gsub(".*: ", "", os_time)
  
  chemo <- GSE39582_RAW[18,]
  chemo <- gsub(".*: ", "", chemo)
  
  braf_mut <- GSE39582_RAW[35,]
  braf_mut <- gsub(".*: ", "", braf_mut)
  
  kras_mut <- GSE39582_RAW[31,]
  kras_mut <- gsub(".*: ", "", kras_mut)
  
  tp53_mut <- GSE39582_RAW[27,]
  tp53_mut <- gsub(".*: ", "", tp53_mut)
  
  GSE39582_CD <- data.frame("dataset"=dataset,
                            "age" = age,
                            "gender"=gender, 
                            "stage"=stage, 
                            "location"=location,
                            "dfs_event"=dfs_event,
                            "dfs_time"=dfs_time,
                            "os_event"=os_event,
                            "os_time"=os_time,
                            "chemo"=chemo,
                            "braf_mut"=braf_mut,
                            "kras_mut"=kras_mut,
                            "tp53_mut"=tp53_mut)
  GSE39582_CD <- t(GSE39582_CD)
  
  colnames(GSE39582_CD) <- Samples
  
  # delete corrupt sample (GSM972425 --> 469)
  GSE39582_CD <- GSE39582_CD[,-469]
  
  # delete non tumoral samples (GSM1681353-71 --> Last 19 samples)
  GSE39582_CD <- GSE39582_CD[ , -((ncol(GSE39582_CD) - 18):ncol(GSE39582_CD))]
  
  
  return(as.data.frame(GSE39582_CD))
  }
  

#------------------------------------------------------------------------------#

####    Variable Homogeneization    ####

#' This function is used to homogenize the results from different datasets 
#' (dfs_event = "Yes/No", "Y/N", "1/0" ---> "1/0" for all datasets) 
#' 
#' @param df                Data frame with combined pheno.data of all datasets
#' @returns                 Data frame with homogenized variables

homogenize_variables <- function(df) {
df <- df %>%
  mutate(gender = case_when(gender == "m" | gender == "male" | gender == "Male" | gender == "M" ~ "M",
                            gender == "f" | gender == "female" | gender == "Female" | gender == "F" ~ "F")) %>%
  mutate(gender = factor(gender, levels = c("M", "F"))) %>%
  
  mutate(dfs_event = case_when(dfs_event == "yes" | dfs_event == "Yes" | dfs_event == "1" | dfs_event == "1 (cancer recurrence)" | dfs_event == "recurrence" | dfs_event == "Y" ~ "1",
                               dfs_event == "no" | dfs_event == "No" | dfs_event == "0" | dfs_event == "0 (no recurrence)" | dfs_event == "no recurrence" | dfs_event == "N" ~ "0")) %>%
  mutate(dfs_event = factor(dfs_event, levels = c("1", "0")))%>%  
  
  mutate(dss_event = case_when(dss_event == "1 (death from cancer)" | dss_event == "death" ~ "1",
                               dss_event == "0 (no death)" | dss_event == "no death" ~ "0")) %>%
  mutate(dss_event = factor(dss_event, levels = c("1", "0")))%>%
  
  mutate(braf_mut = case_when(braf_mut == "yes" | braf_mut == "Yes" | braf_mut == "1" | braf_mut == "Y" | braf_mut == "M" ~ "1",
                              braf_mut == "no" | braf_mut == "No" | braf_mut == "0" | braf_mut == "N"| braf_mut == "WT" ~ "0")) %>%
  mutate(braf_mut = factor(braf_mut, levels = c("1", "0")))%>%
  
  mutate(kras_mut = case_when(kras_mut == "yes" | kras_mut == "Yes" | kras_mut == "1" | kras_mut == "Y" | kras_mut == "M" ~ "1",
                              kras_mut == "no" | kras_mut == "No" | kras_mut == "0" | kras_mut == "N"| kras_mut == "WT" ~ "0")) %>%
  mutate(kras_mut = factor(kras_mut, levels = c("1", "0")))%>%
  
  mutate(tp53_mut = case_when(tp53_mut == "yes" | tp53_mut == "Yes" | tp53_mut == "1" | tp53_mut == "Y" | tp53_mut == "M" ~ "1",
                              tp53_mut == "no" | tp53_mut == "No" | tp53_mut == "0" | tp53_mut == "N"| tp53_mut == "WT" ~ "0")) %>%
  mutate(tp53_mut = factor(tp53_mut, levels = c("1", "0")))%>%
  
  mutate(chemo = case_when(chemo == "yes" | chemo == "Yes" | chemo == "1" | chemo == "Y" | chemo == "M" ~ "1",
                           chemo == "no" | chemo == "No" | chemo == "0" | chemo == "N"| chemo == "WT" ~ "0")) %>%
  mutate(chemo = factor(chemo, levels = c("1", "0")))%>%
  
  mutate(radio = case_when(radio == "yes" | radio == "Yes" | radio == "1" | radio == "Y" | radio == "M" ~ "1",
                           radio == "no" | radio == "No" | radio == "0" | radio == "N"| radio == "WT" ~ "0")) %>%
  mutate(radio = factor(radio, levels = c("1", "0")))%>%
  
  mutate(stage = case_when(stage == "1" | stage == "A" ~ "1",
                           stage == "2" | stage == "B" ~ "2",
                           stage == "3" | stage == "C" ~ "3",
                           stage == "4" | stage == "D" ~ "4")) %>%
  mutate(stage = factor(stage, levels = c("1", "2", "3", "4")))%>%
  
  mutate(location = case_when(location == "right" | location == "Right"  | location == "distal" | location == "Distal" ~ "Distal",
                              location == "left" | location == "Left" | location == "Rectum" | location == "proximal" | location == "Proximal" ~ "Proximal",
                              location == " " | location == "NA" | location == "Colon" ~ NA)) %>%
  
  mutate(first_author = case_when(dataset == "GSE14333" ~ "Jorissen RN",
                                  dataset == "GSE143985" ~ "Shinto E",
                                  dataset == "GSE17536" ~ "Smith JJ",
                                  dataset == "GSE17537" ~ "Smith JJ",
                                  dataset == "GSE33114" ~ "de Sousa E Melo F",
                                  dataset == "GSE38832" ~ "Tripathi MK",
                                  dataset == "GSE39582" ~ "Marisa L")) %>%
  
  
  mutate(dfs_time = na_if(dfs_time, "NA")) %>%
  mutate(dss_time = na_if(dss_time, "NA"))

}



