#.rs.restartR()
library(rcdk)
library(rJava)
library(dplyr)
library(extraTrees)
library(pROC)
library(stringr)  
options(java.parameters = "-Xmx10g")
#####main function######

ET_MT <- function(feature_type){

  #read file#
  path <-str_c('/home1/rhlin/Zebrafish_final202310_withMorphologyALLdata_ECFP6.xlsx')
  dat <- readxl::read_excel(path)
  
  if("SMILES" %in% colnames(dat))
  {dat<- subset( dat, select = -SMILES )
  }

  #data partition#
  n_train <- 7
  n_valid <- 1
  n_test <- 2
  
  
  ifelse (feature_type=='ECFP',ECFP <- T,ECFP <-F)    
  ifelse (feature_type=='mol2vec',mol2vec <- T,mol2vec <- F)
  ifelse (feature_type=='MMBART(BID)',MMBART_BID<- T,MMBART_BID<- F)
  ifelse (feature_type=='MMBART(B)',MMBART_B <- T,MMBART_B <- F)
  ifelse (feature_type=='MMBART(F)',MMBART_F <- T,MMBART_F <- F)    
  ifelse (feature_type=='smiles2vec',smiles2vec <- T,smiles2vec <- F)

  #setting#
  doTraining <- T
  modelType <- "et"
  #testSet <- "train"
  #testSet <- "valid"
  testSet <- "test"
  
  
  #dat <- getFeatures(dat,fp=fp,padel=padel,fpType = fpType,fp_output = T)
  selectCol <- c()
  
  if(ECFP){
    selectCol <- grep("^[0-9]+$",colnames(dat),value = T)
  }
  if(MMBART_F){
    selectCol <- grep("^[0-9]+$",colnames(dat),value = T)
  }
  if(MMBART_B){
    selectCol <- grep("^[0-9]+$",colnames(dat),value = T)
  }
  if(MMBART_BID){
    selectCol <- grep("^[0-9]+$",colnames(dat),value = T)
  }
  if(smiles2vec){
    selectCol <- grep("^[0-9]+$",colnames(dat),value = T)
  }
  if(mol2vec){
    selectCol <- grep("^mol2vec-[0-9]+$",colnames(dat),value = T)
  }
  
  
  model_dat <- do.call("rbind",lapply(2:49,function(t){
    task_dat <- select(dat,No,selectCol,colnames(dat)[t]) %>% rename(target=colnames(dat)[t])
    task_dat$task <- t-1
    task_dat$task_name <- colnames(dat)[t]
    task_dat<-task_dat[complete.cases(task_dat), ]
    task_dat <- task_dat[order(task_dat$target),]
    task_dat$train_test <- rep(c(rep("train",n_train),rep("valid",n_valid),rep("test",n_test)),length.out=nrow(task_dat))
    task_dat <- task_dat[order(task_dat$No),]
    return(task_dat)
  }))
  
  #model_dat<-model_dat%>% drop_na()  


  
  mtry <- floor(log2(1024))
  #train
  model_dat_train <- model_dat %>% filter(model_dat$train_test=="train") %>% select(-train_test)
  if(doTraining){
    if(modelType == "et"){
      set.seed(8)
      et <- extraTrees(x = as.matrix(model_dat_train %>% select(-task,-task_name,-No,-target)),
                       y = as.factor(model_dat_train$target),
                       numThreads = 16,
                       tasks = model_dat_train$task,
                       ntree = 500,
                       mtry=mtry,                
                       na.action='fuse')
    

      
    }else{
      sts <- lapply(unique(model_dat_train$task),function(t){
        print(t)
        model_dat_train <- model_dat_train %>% filter(model_dat_train$task == t)
        set.seed(8)
        st <- extraTrees(x = as.matrix(model_dat_train %>% select(-task,-task_name,-No,-target)),
                         y = as.factor(model_dat_train$target),
                         numThreads = 16,
                         tasks = NULL,                      
                         ntree = 500,
                         na.action='fuse'
        )
        
      
        
        
        return(st)
      })
      names(sts) <- unique(model_dat_train$task)
    }
    
    
  }else{
  }
  
  ###train####
  if(modelType == "et"){
    p <- predict(et, model_dat %>% filter(model_dat$train_test == testSet) %>% select(
      selectCol
    ),
    newtasks = filter(model_dat,model_dat$train_test == testSet)$task )
    prob <- predict(et, model_dat %>% filter(model_dat$train_test == testSet) %>% select(
      selectCol
    ),
    newtasks = filter(model_dat,model_dat$train_test == testSet)$task, probability=TRUE)    
    target <- filter(model_dat, model_dat$train_test == testSet)$target
    tp <- sum((target=='1') & (target==p))
    tn <- sum((target=='0') & (target==p))
    fp <- sum((target=='1') & (target!=p))
    fn <- sum((target=='0') & (target!=p))
    
    roc_obj <- pROC::roc(as.numeric(target), 
                         as.numeric(prob[,2]), 
                         levels = c(0, 1), direction = "<")      
    auc <- pROC::auc(roc_obj)      
    acc <- (tp+tn)/(tp+fp+fn+tn)
    precision <- tp/(tp+fp)
    recall <- tp/(tp+fn)
    speci <- tn/(fp+tn)
    auc<-round(auc, digits = 3)
    acc<-round(acc, digits = 3)
    precision<-round(precision, digits = 3)
    recall<-round(recall, digits = 3)   
    speci <-round(speci, digits = 3) 
    
    print(paste0("/",auc,"/",acc,"/",precision,"/",recall,
                 "/",speci))
    ##################################################################################    
    lapply(1:48,function(t){
      p <- predict(et, model_dat %>% filter(model_dat$train_test == testSet & model_dat$task == t) %>% select(
        selectCol
      ),
      newtasks = filter(model_dat,model_dat$train_test == testSet & model_dat$task == t)$task )
      prob <- predict(et, model_dat %>% filter(model_dat$train_test == testSet & model_dat$task == t) %>% select(
        selectCol
      ),
      newtasks = filter(model_dat,model_dat$train_test == testSet & model_dat$task == t)$task, probability=TRUE)   
      
      target <- filter(model_dat, model_dat$train_test == testSet & model_dat$task == t)$target
      tp <- sum((target=='1') & (target==p))
      tn <- sum((target=='0') & (target==p))
      fp <- sum((target=='1') & (target!=p))
      fn <- sum((target=='0') & (target!=p))
      roc_obj <- pROC::roc(as.numeric(target), 
                           as.numeric(prob[,2]), 
                           levels = c(0, 1), direction = "<")      
      auc <- pROC::auc(roc_obj)      
      acc <- (tp+tn)/(tp+fp+fn+tn)
      precision <- tp/(tp+fp)
      recall <- tp/(tp+fn)
      speci <- tn/(fp+tn)
      auc<-round(auc, digits = 3)
      acc<-round(acc, digits = 3)
      precision<-round(precision, digits = 3)
      recall<-round(recall, digits = 3)   
      speci <-round(speci, digits = 3) 
      
      print(paste0("/",auc,"/",acc,"/",precision,"/",recall,
                   "/",speci))
      
    })
    
    
  }else if(modelType =="st"){
    lapply(1:48,function(t){
      if(!doTraining){
        load(paste0("/home1/rhlin/fish/MLMTLPaper_File/R_model/-",t,ECFP,".rdata"))
      }else{
        st <- sts[[t]]
      }
      p <- predict(st, model_dat %>% filter(model_dat$train_test == testSet & model_dat$task == t) %>% select(
        selectCol
      ),
      newtasks = NULL )
      prob <- predict(st, model_dat %>% filter(model_dat$train_test == testSet & model_dat$task == t) %>% select(
        selectCol
      ),
      newtasks = NULL, probability=TRUE )      
      target <- filter(model_dat, model_dat$train_test == testSet & model_dat$task == t)$target
      tp <- sum((target=='1') & (target==p))
      tn <- sum((target=='0') & (target==p))
      fp <- sum((target=='1') & (target!=p))
      fn <- sum((target=='0') & (target!=p))
      
      roc_obj <- pROC::roc(as.numeric(target), 
                           as.numeric(prob[,2]), 
                           levels = c(0, 1), direction = "<")      
      auc <- pROC::auc(roc_obj)      
      acc <- (tp+tn)/(tp+fp+fn+tn)
      precision <- tp/(tp+fp)
      recall <- tp/(tp+fn)
      speci <- tn/(fp+tn)
      auc<-round(auc, digits = 3)
      acc<-round(acc, digits = 3)
      precision<-round(precision, digits = 3)
      recall<-round(recall, digits = 3)   
      speci <-round(speci, digits = 3) 
      
      print(paste0("/",auc,"/",acc,"/",precision,"/",recall,
                   "/",speci))
    })
  }
}



######################################
Operation
######################################
x.list <- c('ECFP')

for (i in x.list){
  print(i)
  ET_MT(i)
  
}
