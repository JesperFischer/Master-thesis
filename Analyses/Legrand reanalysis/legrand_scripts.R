
get_data_nort = function(rt_threshold_s = 0.1){
  
  select <- dplyr::select
  
  
  df = data.frame()
  for(files in list.files(here::here("..","CardioceptionPaper","data","raw","HRD"), recursive = T)){
    
    df1 = read.csv(paste0(here::here("..","CardioceptionPaper","data","raw","HRD"),"/",files))
    df1$id = sub("^(.*?)/.*", "\\1", files)
    
    
    df = rbind(df,df1)
  }
  df$session = 1
  
  df1 = df %>% filter(DecisionRT > rt_threshold_s)
  
  print(paste0("removed ", nrow(df)-nrow(df1), " trials"))
  
  df_ses1 = df1
  
  
  df = data.frame()
  for(files in list.files(here::here("..","CardioceptionPaper","data","raw","HRD2"), recursive = T)){
    
    df1 = read.csv(paste0(here::here("..","CardioceptionPaper","data","raw","HRD2"),"/",files))
    df1$id = sub("^(.*?)/.*", "\\1", files)
    
    
    df = rbind(df,df1)
  }
  df$session = 2
  
  df1 = df %>% filter(DecisionRT > rt_threshold_s)
  
  print(paste0("removed ", nrow(df)-nrow(df1), " trials"))
  
  
  df_ses2 = df1
  
  df_ses1 = df_ses1 %>% select(Condition,Modality,Decision,DecisionRT,Alpha,nTrials,id,session, EstimatedThreshold, EstimatedSlope)
  
  df_ses2 = df_ses2 %>% select(Condition,Modality,Decision,DecisionRT,Alpha,nTrials,id,session, EstimatedThreshold, EstimatedSlope)
  
  ids2 = unique(df_ses2$id)
  ids1 = unique(df_ses1$id)
  
  df_ses1 = df_ses1 %>% filter(id %in% ids2)
  df_ses2 = df_ses2 %>% filter(id %in% ids1)
  
  df = rbind(df_ses1, df_ses2)
  
  
  hrd<-df%>%
    #select conditions of interest  
    filter(Modality=="Intero")%>% 
    #condensate the table so that there iss one line per intensity
    group_by(Alpha,id,session)%>%
    summarise(x=mean(Alpha),
              n=sum(nTrials<1000),
              y=sum(Decision=="More"),
              s=0,
              ID=mean(as.numeric(str_replace(id,'sub_',''))))%>%
    ungroup()%>%
    select(x,n,y,ID,s,session)
  
  hrd2 = hrd
  
  #change id so that it start at 1 and goes up to N_participant in increment of 1  
  old_ID<-unique(hrd$ID)
  
  for(idx in 1:length(old_ID)){
    hrd$s[hrd$ID==old_ID[idx]]<-idx
  }
  hrd<-select(hrd,x,n,y,s,session)
  
  
  
  return(hrd)
  
}


get_data = function(rt_threshold_s = 0.1){
  
  select <- dplyr::select
  df = data.frame()
  for(files in list.files(here::here("Analyses","Legrand reanalysis","raw data","HRD"), recursive = T)){
    
    df1 = read.csv(paste0(here::here("Analyses","Legrand reanalysis","raw data","HRD"),"/",files))
    df1$id = sub("^(.*?)/.*", "\\1", files)
    
    
    df = rbind(df,df1)
  }
  
  
  df$session = 1
  
  df1 = df %>% filter(DecisionRT > rt_threshold_s)
  
  print(paste0("removed ", nrow(df)-nrow(df1), " trials"))
  
  
  df_ses1 = df1
  
  
  df = data.frame()
  for(files in list.files(here::here("Analyses","Legrand reanalysis","raw data","HRD2"), recursive = T)){
    
    df1 = read.csv(paste0(here::here("Analyses","Legrand reanalysis","raw data","HRD2"),"/",files))
    df1$id = sub("^(.*?)/.*", "\\1", files)
    
    
    df = rbind(df,df1)
  }
  df$session = 2
  
  
  df1 = df %>% filter(DecisionRT > rt_threshold_s)
  
  print(paste0("removed ", nrow(df)-nrow(df1), " trials"))
  
  
  df_ses2 = df1
  
  df_ses1 = df_ses1 %>% select(Condition,Modality,Decision,DecisionRT,Confidence,Alpha,nTrials,id,session, EstimatedThreshold, EstimatedSlope,DecisionRT)
  
  df_ses2 = df_ses2 %>% select(Condition,Modality,Decision,DecisionRT,Confidence,Alpha,nTrials,id,session, EstimatedThreshold, EstimatedSlope,DecisionRT)
  
  ids2 = unique(df_ses2$id)
  ids1 = unique(df_ses1$id)
  
  df_ses1 = df_ses1 %>% filter(id %in% ids2)
  df_ses2 = df_ses2 %>% filter(id %in% ids1)
  
  df = rbind(df_ses1, df_ses2)
  
  hrd<-df%>%
    #select conditions of interest  
    filter(Modality=="Intero")%>% 
    #condensate the table so that there iss one line per intensity
    group_by(Alpha,id,session,nTrials, Confidence)%>%
    summarise(x=mean(Alpha),
              n=sum(nTrials<1000),
              y=sum(Decision=="More"),
              RT = DecisionRT,
              Confidence = Confidence,
              s=0,
              ID=mean(as.numeric(str_replace(id,'sub_',''))))%>%
    ungroup()%>%
    select(x,n,y,ID,s,session,RT,Confidence)
  
  hrd2 = hrd
  
  #change id so that it start at 1 and goes up to N_participant in increment of 1  
  old_ID<-unique(hrd$ID)
  
  for(idx in 1:length(old_ID)){
    hrd$s[hrd$ID==old_ID[idx]]<-idx
  }
  hrd<-select(hrd,x,n,y,s,session, RT, Confidence)
  
  return(hrd)
  
}

