######################################################################################################################################################33
######################################################################################################################################################33
##                                                                                                                                                    ##
##  R-code: all analyses for Germain, Feng, et al. (2023) -  Species-specific traits mediate avian demographic responses under past climate change                                                                             ##
##        
##                                                                                                                                                    ##
##  
##                                                                                                                                                    ##
#######################################################################################################################################################33
######################################################################################################################################################33


library(lme4)
library(lmerTest)
library(Hmisc)
library (tseries)
library(car)
library(RColorBrewer)
library(rgdal)
library(plotrix)
library(rcompanion)
library(data.table)
library(dplyr)
library(jpeg)
library(png)
library(wesanderson)


library(RColorBrewer)
library(pheatmap)
library(phytools)

library(ape)
library(geiger)
library(ggplot2)
library(phylopath)
library(phylolm)

library(MuMIn)                          
library(predictmeans)
library(DHARMa)


#####################################################
#####################################################3
### code to extract raw Ne data from .txt files for each species and summarize overall trends


args <- commandArgs(trailingOnly = TRUE)
NET.file=args[1]
out.file=args[2]


filelist = list.files(pattern = ".*.txt")

## Create a dataframe which will have the summary output for each species
BIG_OUTPUT<-as.data.frame(filelist)
## add a column with just the B10k id minus the '.NE-T.txt'
BIG_OUTPUT$B10k_ID<-substr(BIG_OUTPUT$filelist,1,nchar(as.character(BIG_OUTPUT$filelist))-13)
BIG_OUTPUT$Time_Ne_max<-as.numeric(NA)
BIG_OUTPUT$Time_Ne_min<-as.numeric(NA)
BIG_OUTPUT$Range_Ne<-as.numeric(NA)
BIG_OUTPUT$Var_Ne<-as.numeric(NA)
BIG_OUTPUT$Mean_Ne<-as.numeric(NA)
BIG_OUTPUT$SD_Ne<-as.numeric(NA)
BIG_OUTPUT$count_Ne<-as.integer(NA)

    
All_PSMC<-data.frame(time = numeric(), Ne = numeric(), tag = character(), i = numeric(),time.l10 = numeric(), B10k_ID = character(), time_point = numeric(), normalized_Ne = numeric(), Mean_Ne = numeric(), SD_Ne = numeric(), count_Ne=integer(), stringsAsFactors = F)


## Run for each file in the file list 
for (aaa in 1:length(filelist))
{
  
  NET<-read.delim(filelist[aaa], header = T)      
  NET$time<-row.names(NET)
  names(NET)<-c("Ne","Time")
  row.names(NET)<-c(1:length(row.names(NET)))
  
  
  NET.df<-as.data.frame(NET)
  NET.df$time<-as.numeric(NET.df$Time)
  
  NET.df$time.l10<-log10(NET.df$time)
  NET.df$B10k_ID<-rep(filelist[aaa])
  NET.df$time_point<-as.numeric(1:nrow(NET.df))
  NET.df$normalized_Ne<-NET.df$Ne/max(NET.df$Ne)
  

  
  ## First cut off data at the 30kya point (and max 1mill), then locate max and min
  cut_recent_time<-subset(NET.df,time >= 30000) 
  cut_recent_time<-subset(cut_recent_time,time <= 1000000)
  
  if (length(row.names(cut_recent_time)) > 0)  ### to skip birds with no Ne values within the time period of interest
  {
    
     All_PSMC<-rbind(All_PSMC,cut_recent_time)
    
     temp_name<-substr(filelist[aaa],1,nchar(as.character(filelist[aaa]))-4)
   
     pdf(paste(temp_name,"_plot",".pdf", sep = ""),width= 7.02, height= 3.92)
          ggplot(NET.df)+geom_step(aes(x=time.l10,y=Ne),size=1.2)+
          geom_vline(xintercept=4.4)
     invisible(dev.off())
  
     
  find_max_time<-which(cut_recent_time$Ne == max(cut_recent_time$Ne))
  BIG_OUTPUT[aaa,"Time_Ne_max"]<-cut_recent_time[max(find_max_time),"time" ]
  
  find_min_time<-which(cut_recent_time$Ne == min(cut_recent_time$Ne))
  BIG_OUTPUT[aaa,"Time_Ne_min"]<-cut_recent_time[max(find_min_time),"time"]
  
  }
  
  if (is.na(BIG_OUTPUT[aaa,"Time_Ne_max"]) == FALSE) #two entries with NA for ne max, creates infinity value for range
  {
    BIG_OUTPUT[aaa,"Range_Ne"]<-(max(cut_recent_time$Ne) - min(cut_recent_time$Ne))
    BIG_OUTPUT[aaa,"Var_Ne"]<-var(cut_recent_time$Ne)
  }
  

  BIG_OUTPUT[aaa,"Mean_Ne"]<-mean(cut_recent_time$Ne)
  BIG_OUTPUT[aaa,"SD_Ne"]<-sd(cut_recent_time$Ne)
  BIG_OUTPUT[aaa,"count_Ne"]<-length(cut_recent_time$Ne)

}      


### Quick visual summaries
summary(BIG_OUTPUT$Range_Ne)
summary(BIG_OUTPUT$Var_Ne)
   

All_PSMC$B10k_ID<-substr(All_PSMC$B10k_ID,1,nchar(as.character(All_PSMC$B10k_ID))-13)


                      # # ###  Side-Bar, to generate csv files for figure 2B (curve examples, plug into loop above)
                      # #B10K-Casuarius_casuarius
                        #aaa<-43
                       # write.csv(NET.df, "curve_B10K-Casuarius_casuarius_UPDATED_NE_MAY_2022.csv")
                      #
                      # # #B10K-DU-030-22
                         #aaa<-249
                         #write.csv(NET.df, "curve_B10K-DU-030-22_UPDATED_NE_MAY_2022.csv")
                      # #
                      # #B10K-DU-001-46
                      # aaa<-97
                      # write.csv(NET.df, "curve_B10K-DU-001-46_UPDATED_NE_MAY_2022.csv")
                      #
                      # #APP-013
                        #aaa<-12
                        #write.csv(NET.df, "curve_APP-013_UPDATED_NE_MAY_2022.csv")
                      #
                      # #B10K-Nothocercus_nigrocapillus
                       # aaa<-258
                        # write.csv(NET.df, "curve_B10K-Nothocercus_nigrocapillus_UPDATED_NE_MAY_2022.csv")

                    # update
                      # B10K-DU-001-08
                     # aaa<-70
                     #  write.csv(NET.df, "curve_B10K-DU-001-08_UPDATED_NE_MAY_2022.csv")







#### load example curves for graphing later
Curve_B10K_Casuarius_casuarius<-read.csv('curve_B10K-Casuarius_casuarius_UPDATED_NE_MAY_2022.csv', header=T)
Curve_B10K_DU_001_46<-read.csv('curve_B10K-DU-001-46_UPDATED_NE_MAY_2022.csv', header=T)
Curve_APP_013<-read.csv('curve_APP-013_UPDATED_NE_MAY_2022.csv', header=T)
Curve_B10K_DU_030_22<-read.csv('curve_B10K-DU-030-22_UPDATED_NE_MAY_2022.csv', header=T)
curve_B10K_DU_001_08_UPDATED_NE<-read.csv('curve_B10K-DU-001-08_UPDATED_NE_MAY_2022.csv', header=T)


#####################################################
#####################################################3
### End





#####################################################
#####################################################3
###  Add in trait data and realm for each species


trait_data<-read.csv('B10K_Trait_Data_031319.csv', header=T)
trait_data$species_name_char<-as.character(trait_data$species)
levels(trait_data$threat_status)<-sapply(levels(trait_data$threat_status), FUN=function(x){ifelse(x %in% c('NA',' NA',''), NA,x)})
trait_data$threat_status_char<-as.character(trait_data$threat_status)

## disagreement in ending of species name for a few, fix in trait data
trait_data[which(trait_data[,"species_name_char"]=="Eolophus_roseicapilla"), 'species_name_char']<-"Eolophus_roseicapillus"
trait_data[which(trait_data[,"species_name_char"]=="Cercotrichas_coryphoeus"), 'species_name_char']<-"Cercotrichas_coryphaeus"


realmdata<-read.csv('Realms.csv', header=T)
realmdata$B10k_ID_char<-as.character(realmdata$B10k_ID)
realmdata$species_name_char<-as.character(realmdata$Species_Name_Standard)
realmdata$mass_m<-as.numeric(NA)
realmdata$mass_f<-as.numeric(NA)
realmdata$mass_unsexed<-as.numeric(NA)
realmdata$clutch_size<-as.numeric(NA)
realmdata$inc_duration<-as.numeric(NA)
realmdata$IUCN<-as.character(NA)
realmdata$migration<-as.character(NA)
realmdata$territory<-as.character(NA)
realmdata$classify<-as.character(NA)
realmdata$maximum_longevity_y<-as.numeric(NA)
realmdata$max_elevation<-as.numeric(NA)
realmdata$range_size_km2<-as.numeric(NA)
realmdata$kipps_distance<-as.numeric(NA)
realmdata$tarsus_length<-as.numeric(NA)
realmdata$Order<-as.character(NA)

realmdata$egg_mass<-as.numeric(NA)
realmdata$mortality_both<-as.numeric(NA)
realmdata$min_elevation<-as.numeric(NA)
realmdata$max_elevation<-as.numeric(NA)
realmdata$wing_chord<-as.numeric(NA)
realmdata$hand_wing_index<-as.numeric(NA)
realmdata$bill_total_culmen<-as.numeric(NA)
realmdata$bill_nares<-as.numeric(NA)
realmdata$bill_width<-as.numeric(NA)
realmdata$bill_depth<-as.numeric(NA)
realmdata$b_brain_size<-as.numeric(NA)





for (ccc in 1:nrow(realmdata))
{
  grab1<-realmdata[ccc,"species_name_char"]
  find2<-which(trait_data[,'species_name_char'] == grab1)
  realmdata[ccc,"mass_m"]<-trait_data[find2, "m_mass"]
  realmdata[ccc,"mass_f"]<-trait_data[find2, "f_mass"]
  realmdata[ccc,"mass_unsexed"]<-trait_data[find2, "unsexed_mass"]
  realmdata[ccc,"clutch_size"]<-trait_data[find2, "clutch_size"]
  realmdata[ccc,"inc_duration"]<-trait_data[find2, "inc_duration"]
  realmdata[ccc,"IUCN"]<-trait_data[find2, "threat_status_char"]
  realmdata[ccc,"migration"]<-as.character(trait_data[find2, "migration"])
  realmdata[ccc,"territory"]<-as.character(trait_data[find2, "territory"])
  realmdata[ccc,"classify"]<-as.character(trait_data[find2, "Ryan_Classify"])
  realmdata[ccc,"maximum_longevity_y"]<-trait_data[find2, "maximum_longevity_y"]
  realmdata[ccc,"range_size_km2"]<-trait_data[find2, "range_size_km2"]
  realmdata[ccc,"kipps_distance"]<-trait_data[find2, "kipps_distance"]
  realmdata[ccc,"tarsus_length"]<-trait_data[find2, "tarsus_length"]
  realmdata[ccc,"Order"]<-as.character(trait_data[find2, "order"])
  
  
  realmdata[ccc,"egg_mass"]<-trait_data[find2, "egg_mass"]
  realmdata[ccc,"mortality_both"]<-trait_data[find2, "mortality_both"]
  realmdata[ccc,"min_elevation"]<-trait_data[find2, "min_elevation"]
  realmdata[ccc,"max_elevation"]<-trait_data[find2, "max_elevation"]
  realmdata[ccc,"wing_chord"]<-trait_data[find2, "wing_chord"]
  realmdata[ccc,"hand_wing_index"]<-trait_data[find2, "hand_wing_index"]
  realmdata[ccc,"bill_total_culmen"]<-trait_data[find2, "bill_total_culmen"]
  realmdata[ccc,"bill_nares"]<-trait_data[find2, "bill_nares"]
  realmdata[ccc,"bill_width"]<-trait_data[find2, "bill_width"]
  realmdata[ccc,"bill_depth"]<-trait_data[find2, "bill_depth"]
  realmdata[ccc,"b_brain_size"]<-trait_data[find2, "b_brain_size"]

}




########## put it all together in the BIG OUTPUT file
BIG_OUTPUT$spcies_name<-as.character(NA)
BIG_OUTPUT$Realm<-as.character(NA)
BIG_OUTPUT$Order<-as.character(NA)
BIG_OUTPUT$Classify<-as.character(NA)
BIG_OUTPUT$mass_m<-as.numeric(NA)
BIG_OUTPUT$mass_f<-as.numeric(NA)
BIG_OUTPUT$mass_unsexed<-as.numeric(NA)
BIG_OUTPUT$clutch_size<-as.numeric(NA)
BIG_OUTPUT$inc_duration<-as.numeric(NA)
BIG_OUTPUT$egg_mass<-as.numeric(NA)
BIG_OUTPUT$maximum_longevity_y<-as.numeric(NA)
BIG_OUTPUT$mortality_both<-as.numeric(NA)

BIG_OUTPUT$IUCN<-as.character(NA)
BIG_OUTPUT$migration<-as.character(NA)
BIG_OUTPUT$territory<-as.character(NA)
BIG_OUTPUT$min_elevation<-as.numeric(NA)
BIG_OUTPUT$max_elevation<-as.numeric(NA)
BIG_OUTPUT$range_size_km2<-as.numeric(NA)

BIG_OUTPUT$kipps_distance<-as.numeric(NA)
BIG_OUTPUT$tarsus_length<-as.numeric(NA)
BIG_OUTPUT$wing_chord<-as.numeric(NA)
BIG_OUTPUT$hand_wing_index<-as.numeric(NA)
BIG_OUTPUT$bill_total_culmen<-as.numeric(NA)
BIG_OUTPUT$bill_nares<-as.numeric(NA)
BIG_OUTPUT$bill_width<-as.numeric(NA)
BIG_OUTPUT$bill_depth<-as.numeric(NA)
BIG_OUTPUT$b_brain_size<-as.numeric(NA)


for (bbb in 1:nrow(BIG_OUTPUT))
{
  focalb10_ID<-BIG_OUTPUT[bbb,'B10k_ID']
  find_focal<-which(realmdata[,'B10k_ID_char'] == focalb10_ID)
  BIG_OUTPUT[bbb,'spcies_name']<-as.character(realmdata[find_focal,'species_name_follow_UCE_tree_char'])
  BIG_OUTPUT[bbb,'Realm']<-as.character(realmdata[find_focal,'Realm_fixed'])
  BIG_OUTPUT[bbb,'Order']<-as.character(realmdata[find_focal,'Order'])
  BIG_OUTPUT[bbb, 'Classify']<-as.character(realmdata[find_focal, "classify"])
  BIG_OUTPUT[bbb,'mass_m']<-as.numeric(realmdata[find_focal,'mass_m'])
  BIG_OUTPUT[bbb,'mass_f']<-as.numeric(realmdata[find_focal,'mass_f'])
  BIG_OUTPUT[bbb,'mass_unsexed']<-as.numeric(realmdata[find_focal,'mass_unsexed'])
  BIG_OUTPUT[bbb,'clutch_size']<-as.numeric(realmdata[find_focal,'clutch_size'])
  BIG_OUTPUT[bbb,'inc_duration']<-as.numeric(realmdata[find_focal,'inc_duration'])
  BIG_OUTPUT[bbb,'IUCN']<-as.character(realmdata[find_focal,'IUCN'])
  BIG_OUTPUT[bbb,'migration']<-as.numeric(realmdata[find_focal,'migration'])
  BIG_OUTPUT[bbb,'territory']<-as.character(realmdata[find_focal,'territory'])
  BIG_OUTPUT[bbb,'maximum_longevity_y']<-as.numeric(realmdata[find_focal,'maximum_longevity_y'])
  BIG_OUTPUT[bbb,'range_size_km2']<-as.numeric(realmdata[find_focal,'range_size_km2'])
  BIG_OUTPUT[bbb,'kipps_distance']<-as.numeric(realmdata[find_focal,'kipps_distance'])
  BIG_OUTPUT[bbb,'tarsus_length']<-as.numeric(realmdata[find_focal,'tarsus_length'])
  BIG_OUTPUT[bbb,'egg_mass']<-as.numeric(realmdata[find_focal,'egg_mass'])
  BIG_OUTPUT[bbb,'mortality_both']<-as.numeric(realmdata[find_focal,'mortality_both'])
  BIG_OUTPUT[bbb,'min_elevation']<-as.numeric(realmdata[find_focal,'min_elevation'])
  BIG_OUTPUT[bbb,'max_elevation']<-as.numeric(realmdata[find_focal,'max_elevation'])
  BIG_OUTPUT[bbb,'wing_chord']<-as.numeric(realmdata[find_focal,'wing_chord'])
  BIG_OUTPUT[bbb,'hand_wing_index']<-as.numeric(realmdata[find_focal,'hand_wing_index'])
  BIG_OUTPUT[bbb,'bill_total_culmen']<-as.numeric(realmdata[find_focal,'bill_total_culmen'])
  BIG_OUTPUT[bbb,'bill_nares']<-as.numeric(realmdata[find_focal,'bill_nares'])
  BIG_OUTPUT[bbb,'bill_width']<-as.numeric(realmdata[find_focal,'bill_width'])
  BIG_OUTPUT[bbb,'bill_depth']<-as.numeric(realmdata[find_focal,'bill_depth'])
  BIG_OUTPUT[bbb,'b_brain_size']<-as.numeric(realmdata[find_focal,'b_brain_size'])
  
  }


            #### update Oct 2021 - Some manual fixes to IUCN status - determined via IUCN website (all species checked)
            
            
            BIG_OUTPUT[which(BIG_OUTPUT[,"spcies_name"]=="Antrostomus carolinensis"), 'IUCN']<-"NT"
            BIG_OUTPUT[which(BIG_OUTPUT[,"spcies_name"]=="Apteryx rowi"), 'IUCN']<-"VU"
            BIG_OUTPUT[which(BIG_OUTPUT[,"spcies_name"]=="Buceros rhinoceros"), 'IUCN']<-"VU"
            BIG_OUTPUT[which(BIG_OUTPUT[,"spcies_name"]=="Bucorvus abyssinicus"), 'IUCN']<-"VU"
            BIG_OUTPUT[which(BIG_OUTPUT[,"spcies_name"]=="Calcarius ornatus"), 'IUCN']<-"VU"
            BIG_OUTPUT[which(BIG_OUTPUT[,"spcies_name"]=="Callaeas wilsoni"), 'IUCN']<-"NT"
            BIG_OUTPUT[which(BIG_OUTPUT[,"spcies_name"]=="Casuarius casuarius"), 'IUCN']<-"LC"
            BIG_OUTPUT[which(BIG_OUTPUT[,"spcies_name"]=="Chaetura pelagica"), 'IUCN']<-"VU"
            BIG_OUTPUT[which(BIG_OUTPUT[,"spcies_name"]=="Edolisoma coerulescens"), 'IUCN']<-"LC"
            BIG_OUTPUT[which(BIG_OUTPUT[,"spcies_name"]=="Hemignathus wilsoni"), 'IUCN']<-"EN"
            BIG_OUTPUT[which(BIG_OUTPUT[,"spcies_name"]=="Lanius ludovicianus"), 'IUCN']<-"NT"
            BIG_OUTPUT[which(BIG_OUTPUT[,"spcies_name"]=="Nestor notabilis"), 'IUCN']<-"EN"
            BIG_OUTPUT[which(BIG_OUTPUT[,"spcies_name"]=="Pedionomus torquatus"), 'IUCN']<-"CR"
            BIG_OUTPUT[which(BIG_OUTPUT[,"spcies_name"]=="Pelecanus crispus"), 'IUCN']<-"NT"
            BIG_OUTPUT[which(BIG_OUTPUT[,"spcies_name"]=="Rissa tridactyla"), 'IUCN']<-"VU"
            BIG_OUTPUT[which(BIG_OUTPUT[,"spcies_name"]=="Setophaga coronata"), 'IUCN']<-"LC"



## do the same for the all_PSMC data

All_PSMC$realm<-as.character(NA)

for (ddd in 1:nrow(All_PSMC))
{
  focalb10_ID_PSMC<-All_PSMC[ddd,'B10k_ID']
  find_focal<-which(BIG_OUTPUT[,'B10k_ID'] == focalb10_ID_PSMC)
  All_PSMC[ddd,"realm"]<-BIG_OUTPUT[find_focal,"Realm"]
}



## also classify species as 'at risk' if they are critically endangered, endangered, or vulnerable

BIG_OUTPUT$Conservation<-as.character(NA)
BIG_OUTPUT$Pass_NonPass<-as.character(NA)
for (eee in 1:nrow(BIG_OUTPUT))
{
  
  focal_species_Cons<-BIG_OUTPUT[eee,"IUCN"]
  if (is.na(focal_species_Cons) == FALSE)
  {
          if ((focal_species_Cons == "CR") || (focal_species_Cons == "EN") || (focal_species_Cons == "VU")|| (focal_species_Cons == "NT"))
            {
            BIG_OUTPUT[eee, "Conservation"]<-"At_Risk"
          }else {
            BIG_OUTPUT[eee, "Conservation"]<-"Not_At_Risk"
          }
  }
  
  focal_species_Order<-BIG_OUTPUT[eee,"Order"]
  if (is.na(focal_species_Order) == FALSE)
  {
    if ((focal_species_Order == "Passeriformes"))
    {
      BIG_OUTPUT[eee, "Pass_NonPass"]<-"Passerine"
    }else {
      BIG_OUTPUT[eee, "Pass_NonPass"]<-"Non_Passerine"
    }
    
  }
  
}



BIG_OUTPUT$Conservation_f<-as.factor(BIG_OUTPUT$Conservation)
BIG_OUTPUT$Pass_NonPass_f<-as.factor(BIG_OUTPUT$Pass_NonPass)


#####################################################3
#####################################################
### End





#####################################################3
#####################################################
###  Cluster analysis for figure 1

NorNe.data<-read.table("30_1000.P121.normalized.txt",sep=",")
names(NorNe.data)[1]<-"id"
spe.sort<-read.table("263.spe.name.txt",sep="\t")

NorNe.data$id<-factor(NorNe.data$id,levels=spe.sort$V1)
NorNe.data.sort<-data.frame(NorNe.data %>% arrange(id))
row.names(NorNe.data.sort)<-NorNe.data.sort$id
NorNe.data.sort.m<-as.matrix(NorNe.data.sort[,2:length(names(NorNe.data.sort))])
spe.sort$V3<-ifelse(spe.sort$V3=="Passeriformes","Passerine","Non-Passerine")

annotation_row = data.frame(
  Order=spe.sort$V3
)
rownames(annotation_row)=rownames(NorNe.data.sort.m)
annotation_color = list(
  Order = c("Passerine" = "#ed8ec2","Non-Passerine" = "#009e73")
)

temp.nor<-
  pheatmap(NorNe.data.sort.m,display_numbers=F,
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
           angle_col=90,show_rownames = T,fontsize_row=3,show_colnames = F,legend = T,
           cluster_row = T, cluster_col = F,
           annotation_row=annotation_row,
           annotation_colors=annotation_color,
           main="30-1000(kya) normalized 121 sampling points",
           filename = "Ne.pdf", width=6.5, height=6.5
  )

cluster.k=7
rownames(annotation_row)=rownames(NorNe.data.sort.m)
Cluster.tmp<-factor(annotation_row$Cluster)
annotation_row$Cluster<-factor(ifelse(Cluster.tmp==1,7,
                                      ifelse(Cluster.tmp==2,5,
                                             ifelse(Cluster.tmp==3,3,
                                                    ifelse(Cluster.tmp==4,6,
                                                           ifelse(Cluster.tmp==5,2,
                                                                  ifelse(Cluster.tmp==6,4,1)))))),levels = c("1","2","3","4","5","6","7"))

annotation_color = list(
  Clade = c("Passerine" = "#ed8ec2","Non-Passerine" = "#009e73"),
  Cluster=c("#BF5B17","#386CB0","#FDC086","#F0027F","#BEAED4","#FFFF99","#7FC97F")
)

names(annotation_color$Cluster)=c("7","5","3","6","2","4","1")[1:cluster.k]
pheatmap(NorNe.data.sort.m,display_numbers=F,
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
         angle_col=90,show_rownames = F,fontsize_row=3,show_colnames = F,legend = T,
         cluster_row = T, cluster_col = F,
         annotation_row=annotation_row,
         annotation_colors=annotation_color,
         main="30-1000(kya) normalized 121 sampling points",
         filename = paste("cluster_",cluster.k,".pdf",sep=""), width=5.3, height=6.5
)

#####################################################3
#####################################################
### End




#####################################################3
#####################################################
### Adding this info to 'big output' dataset


New_Cluster_Groups<-read.csv("clusters_different_k.csv", header = T)
New_Cluster_Groups$B10k_ID_char<-as.character(New_Cluster_Groups$ID)
New_Cluster_Groups$New_Cluster_f<-as.character(New_Cluster_Groups$K.7)  #### if you change number of clusters, do it here.


BIG_OUTPUT$New_Clusters_MAY_2022<-as.character(NA)

for (xxx in 1:nrow(BIG_OUTPUT))
    {
         focal_species_group<-BIG_OUTPUT[xxx,"B10k_ID"]
         lookup_clusters<-which(New_Cluster_Groups$B10k_ID_char == focal_species_group)
       
         if (length(lookup_clusters) > 0)
         {
           BIG_OUTPUT[xxx, "New_Clusters_MAY_2022"]<-New_Cluster_Groups[lookup_clusters, "New_Cluster_f"]
         }
           
        
    }

BIG_OUTPUT$New_Clusters_MAY_2022_f<-as.factor(BIG_OUTPUT$New_Clusters_MAY_2022)



#### chi-square test looking at cons status and pass/non-pass across the four cluster groups

Cluster_info<-subset(BIG_OUTPUT, is.na(BIG_OUTPUT$New_Clusters_MAY_2022_f) == FALSE)


chi_clus<-table(Cluster_info$Conservation_f, Cluster_info$New_Clusters_MAY_2022_f)
chisq.test(chi_clus)
tapply(Cluster_info$New_Clusters_MAY_2022_f, Cluster_info$Conservation_f, summary)

chi_pass<-table(Cluster_info$Pass_NonPass_f, Cluster_info$New_Clusters_MAY_2022_f)
chisq.test(chi_pass)
tapply(Cluster_info$New_Clusters_MAY_2022_f, Cluster_info$Pass_NonPass_f, summary)

#####################################################3
#####################################################
### End








#####################################################3
#####################################################
### Analyzing Patterns over time


### FIG 2  - Species curves 
#### Passerine vs. Nonpasserine Ne
    
  #### Bring in CSV with log Ne values for all 120 time points. This uses the 263  species for which we have full values of Ne across the whole time period. 
 
 LogNe_all_species_over_time<-read.csv("Log_Ne_Full_Time_Period_equidistant_time_points.csv",header=TRUE)   
 LogNe_all_species_over_time$B10k_id_char<-as.character(LogNe_all_species_over_time$ID)
 LogNe_all_species_over_time$Pass_Nonpass<-as.character(NA)   

 
 ### add in whether species is passerine/nonpasserine    
 for (jjj in 1:nrow(LogNe_all_species_over_time))  
 {
   focal_species_pass_nonpass<-LogNe_all_species_over_time[jjj,"B10k_id_char"]
   
   lookup_ID<-which(BIG_OUTPUT$B10k_ID == focal_species_pass_nonpass)
   
   if (BIG_OUTPUT[lookup_ID,"Order"] == "Passeriformes")
         {
            LogNe_all_species_over_time[jjj, "Pass_Nonpass"]<-"Passerine"
         } else
         {
           LogNe_all_species_over_time[jjj, "Pass_Nonpass"]<-"Nonpasserine"
         } 
   
 }
    
 LogNe_all_species_over_time$Pass_Nonpass_f<-as.factor(LogNe_all_species_over_time$Pass_Nonpass)
 summary(LogNe_all_species_over_time$Pass_Nonpass_f)
 
 #### Plot average Ne over time (with error) by Passerine/nonPasserine   
 Passerines_Log_Ne_time<-subset(LogNe_all_species_over_time, Pass_Nonpass == "Passerine")   
     Passerines_mean_Log_Ne<-as.data.frame(c("Passerine_mean",(apply(Passerines_Log_Ne_time[,2:121],2,mean)))) #calculates the mean value of each column in the data frame
     Passerines_SD_Log_Ne<-as.data.frame(c("Passerine_SD",(apply(Passerines_Log_Ne_time[,2:121],2,sd))))
   
 Non_Passerines_Log_Ne_time<-subset(LogNe_all_species_over_time, Pass_Nonpass == "Nonpasserine")   
     Non_Passerines_mean_Log_Ne<-as.data.frame(c("Non_passerine_mean",(apply(Non_Passerines_Log_Ne_time[,2:121],2,mean)))) #calculates the mean value of each column in the data frame 
     Non_Passerines_SD_Log_Ne<-as.data.frame(c("Non_passerine_SD",(apply(Non_Passerines_Log_Ne_time[,2:121],2,sd))))  
    
Mean_Log_Ne_Over_Time<-as.data.frame(c(Passerines_mean_Log_Ne,Passerines_SD_Log_Ne,Non_Passerines_mean_Log_Ne,Non_Passerines_SD_Log_Ne))
colnames(Mean_Log_Ne_Over_Time)<-c("Passerine_mean","Passerine_SD","Non_passerine_mean","Non_passerine_SD")
Mean_Log_Ne_Over_Time<-Mean_Log_Ne_Over_Time[2:nrow(Mean_Log_Ne_Over_Time),]   
Mean_Log_Ne_Over_Time$timePoint<-1:120
## to put on log scale
Mean_Log_Ne_Over_Time$timePoint_equal<-seq(30000,1000000, by = ((1000000 - 30000)/119))
Mean_Log_Ne_Over_Time$timePoint_equal_LOG<-log10(seq(30000,1000000, by = ((1000000 - 30000)/119)))

Mean_Log_Ne_Over_Time$Passerine_mean<-as.numeric(as.character(Mean_Log_Ne_Over_Time$Passerine_mean))
Mean_Log_Ne_Over_Time$Passerine_SD<-as.numeric(as.character(Mean_Log_Ne_Over_Time$Passerine_SD))
Mean_Log_Ne_Over_Time$Non_passerine_mean<-as.numeric(as.character(Mean_Log_Ne_Over_Time$Non_passerine_mean))
Mean_Log_Ne_Over_Time$Non_passerine_SD<-as.numeric(as.character(Mean_Log_Ne_Over_Time$Non_passerine_SD))

#### T-test of differences among Passerine/non-passerine means
t.test(Mean_Log_Ne_Over_Time$Passerine_mean,Mean_Log_Ne_Over_Time$Non_passerine_mean, paired = TRUE)

        pass_pnd<-Mean_Log_Ne_Over_Time$Passerine_mean+Mean_Log_Ne_Over_Time$Passerine_SD
        pass_nsd<-Mean_Log_Ne_Over_Time$Passerine_mean-Mean_Log_Ne_Over_Time$Passerine_SD
        non_pass_pnd<-Mean_Log_Ne_Over_Time$Non_passerine_mean+Mean_Log_Ne_Over_Time$Non_passerine_SD
        non_pass_nsd<-Mean_Log_Ne_Over_Time$Non_passerine_mean-Mean_Log_Ne_Over_Time$Non_passerine_SD
    

## Plot Pass/Non-Pass figure        
        ry_green<-adjustcolor("#009e73", alpha.f=0.35)
        labelz=c(expression(paste("3 x 10"^"4")),rep("",6),expression("10"^"5"),rep("",8),expression("10"^"6"))
        at<-c(log10(seq(30000,90000, by = ((90000 - 30000)/6))), log10(seq(100000, 1000000, by = ((1000000 - 100000)/9))))
        
        
        plot(Passerine_mean~timePoint_equal_LOG, data = Mean_Log_Ne_Over_Time, 
             ylim = c(-30,160), pch = ".", col = "#ed8ec2", ylab = "Log_Ne", ann = F,axes = F)
        axis(side=2, at=c(-20,0,20,40,60,80,100,120), las=1,cex.axis = 1.3,cex.lab = 1.3, line = -1)
        axis(1, at=at,label=labelz, cex.axis = 1.3,cex.lab = 1.3)
        
        
        lines(Mean_Log_Ne_Over_Time$timePoint_equal_LOG, pass_pnd,col="#CC3399", lwd = 1)
        lines(Mean_Log_Ne_Over_Time$timePoint_equal_LOG, pass_nsd, col="#CC3399", lwd = 1)
        polygon(x=c(Mean_Log_Ne_Over_Time$timePoint_equal_LOG, rev(Mean_Log_Ne_Over_Time$timePoint_equal_LOG)), y=c(pass_pnd, rev(pass_nsd)), col="#ed8ec2", density = 200, angle=90)
        lines(Mean_Log_Ne_Over_Time$timePoint_equal_LOG, Mean_Log_Ne_Over_Time$Passerine_mean, col="#CC3399",lwd=4)
        
        
        polygon(x=c(Mean_Log_Ne_Over_Time$timePoint_equal_LOG, rev(Mean_Log_Ne_Over_Time$timePoint_equal_LOG)), y=c(non_pass_pnd, rev(non_pass_nsd)), col=ry_green, density = 200, angle=90)
        points(Non_passerine_mean~timePoint_equal_LOG, data = Mean_Log_Ne_Over_Time,pch = ".", col = "#009e73")
        lines(Mean_Log_Ne_Over_Time$timePoint_equal_LOG, non_pass_pnd,col="#006600",lwd=1)
        lines(Mean_Log_Ne_Over_Time$timePoint_equal_LOG, non_pass_nsd,col="#006600",lwd=1)
        lines(Mean_Log_Ne_Over_Time$timePoint_equal_LOG, Mean_Log_Ne_Over_Time$Non_passerine_mean, col="#006600",lwd=4)
        
        
        mtext("Time Since Present (Years)", side = 1, line = 3, adj = 0.5, cex = 1.4, outer = F)    
        mtext(expression(paste("Log Effective Population Size ("~italic("N")[italic("e")]~")")), side = 2, line = 2.5, adj = 0.2, cex = 1.4, outer = F)
        
        
        legend(x = 5.54, y = 115, bty = 'n',
               c("Non-Passerine","", "Passerine"),
               cex=1.2,
               title="",
               lty=1,
               pch=NA,
               lwd=4,
               col=c("#009e73", "white", "#ed8ec2"))
            
        
      
##### example curves from threatened species  (Fig 2B)
Curve_B10K_Casuarius_casuarius<-subset(Curve_B10K_Casuarius_casuarius, time >= 30000)
Curve_B10K_Casuarius_casuarius<-subset(Curve_B10K_Casuarius_casuarius, time <= 1050000)

Curve_B10K_DU_001_46<-subset(Curve_B10K_DU_001_46, time >= 30000)
Curve_B10K_DU_001_46<-subset(Curve_B10K_DU_001_46, time <= 1050000)

Curve_APP_013<-subset(Curve_APP_013, time >= 30000)
Curve_APP_013<-subset(Curve_APP_013, time <= 1050000)

#Curve_B10K_Nothocercus_nigrocapillus<-subset(Curve_B10K_Nothocercus_nigrocapillus, time >= 30000)
#Curve_B10K_Nothocercus_nigrocapillus<-subset(Curve_B10K_Nothocercus_nigrocapillus, time <= 1050000)

Curve_B10K_DU_030_22<-subset(Curve_B10K_DU_030_22, time >= 30000)
Curve_B10K_DU_030_22<-subset(Curve_B10K_DU_030_22, time <= 1050000)

curve_B10K_DU_001_08_UPDATED_NE<-subset(curve_B10K_DU_001_08_UPDATED_NE, time >= 30000)
curve_B10K_DU_001_08_UPDATED_NE<-subset(curve_B10K_DU_001_08_UPDATED_NE, time <= 1050000)


       #tiff('Fig_Species_Curves_IUCN_RED_UPDATED_Post_Sci_Rej_May_2022.tiff',width=8000,height=6000,units = 'px', res = 800)
        
                 label=c(expression(paste("3 x 10"^"4")),rep("",6),expression("10"^"5"),rep("",8),expression("10"^"6"))
                    at<-c(seq(3,9,1)*10^4,seq(1,9,1)*10^5,10^6)
                      
                      par(mar = c(5, 6, 0.5, 0.5))
                      plot(Curve_B10K_Casuarius_casuarius$time,Curve_B10K_Casuarius_casuarius$Ne,
                           ylim=c(0,25), las = 1, bty = "n",
                           type="l",col="#FF0000",lwd=3,
                           xlim=c(30000,1000000),xaxt="n",log = "x",
                           xlab="Time Since Present (Years)",
                           cex.axis = 1.3,cex.lab = 1.3,
                           #ylab=expression(paste("Effective Population Size (x10"^"4",")")))
                           ylab=expression(paste("Effective Population Size ("~italic("N")[italic("e")]~"x 10"^"4",")")))
                      
                      
                      lines(Curve_B10K_DU_030_22$time,Curve_B10K_DU_030_22$Ne, col = "#1F78B4", lwd = 3)
                      lines(Curve_APP_013$time,Curve_APP_013$Ne, col = "#FB9A99", lwd = 3)
                      lines(curve_B10K_DU_001_08_UPDATED_NE$time,curve_B10K_DU_001_08_UPDATED_NE$Ne, col = "#FDBF6F", lwd = 3)
                      
                      axis(1,at=at,label=label, cex.axis = 1.3,cex.lab = 1.3)
                      legend(x = 50000, y = 27.5, bty = 'n',
                             c("Southern cassowary","", "Yellowhead", "","MacQueen's bustard", "","White-crested guan"),
                             cex=1.68,
                             title="",
                             lty=1,
                             pch=NA,
                             lwd=4.2,
                             col=c("#FF0000","white","#1F78B4","white","#FB9A99","white","#FDBF6F"))
                      
        
        #dev.off()

#####################################################3
#####################################################
### End



     
                      
                      
                                       
                      
                      

#####################################################3
#####################################################
### Analyzing Patterns across realms

#### curve of all species (normalized)    
All_PSMC_cut<-subset(All_PSMC,time <= 1000000) 
All_PSMC_cut<-subset(All_PSMC_cut, time >= 30000)


All_PSMC_cut$Round_Time<-round(All_PSMC_cut$time, digits = -4) 
Mean_PSMC_All_Species_All_Times<-All_PSMC_cut %>%
  group_by (Round_Time) %>%
  summarize (Mean_ne = mean (normalized_Ne), SD_ne = sd(normalized_Ne), var_ne = var(normalized_Ne))  

Mean_PSMC_All_Species_All_Times$Time_kya<-Mean_PSMC_All_Species_All_Times$Round_Time/1000
Normalized_NE_Smooth<-smooth.spline(x = Mean_PSMC_All_Species_All_Times$Time_kya, y = Mean_PSMC_All_Species_All_Times$Mean_ne)

Mean_CI_All_Species_All_Times<-groupwiseMean(normalized_Ne ~ Round_Time,
                    data   = All_PSMC_cut,
                    conf   = 0.95,
                    digits = 3)

Mean_CI_All_Species_All_Times$Time_kya<-Mean_CI_All_Species_All_Times$Round_Time/1000

### with Confidence Intervals
Normalized_NE_Smooth_CI<-smooth.spline(x = Mean_CI_All_Species_All_Times$Time_kya, y = Mean_CI_All_Species_All_Times$Mean)
plot(Mean_CI_All_Species_All_Times$Time_kya,Mean_CI_All_Species_All_Times$Mean, ylim = c(0.41,0.7),ann = F, axes = F, pch = 20)
lines(Normalized_NE_Smooth_CI, lwd = 3, col = 'red')
axis(side=2, at=c(0.3,0.4,0.5,0.6,0.7), las = 1,cex.axis = 1.15,cex.lab = 1.15)
axis(side=1, at=c(0,200,400,600,800,1000),cex.axis = 1.15,cex.lab = 1.15)
minor.tick(nx=2,ny=2, tick.ratio=0.5)
mtext("Time Since Present (kya)", side = 1, line = 3, adj = 0.5, cex = 1.4, outer = F)
mtext(expression(paste(Normalized~italic("N")[italic("e")])), side = 2, line = 3, adj = 0.5, cex = 1.4, outer = F)
text(900,0.65,"C",cex=1.4)    

    
BIG_OUTPUT$Round_Time_Ne_max_kya<-round(BIG_OUTPUT$Time_Ne_max_kya, digits = -1)    ### to help visualize time breaks (1000 year intervals) on the same scale as normalized data
BIG_OUTPUT$Round_Time_Ne_min_kya<-round(BIG_OUTPUT$Time_Ne_min_kya, digits = -1)    


length(unique(All_PSMC_cut$B10k_ID))

######
######   Bring in Zoogeographic Realms
######

### turn Big0utput realm and IUCN data into a factor, basic summary stats
BIG_OUTPUT$Realm_f<-as.factor(BIG_OUTPUT$Realm)
BIG_OUTPUT$IUCN_f<-as.factor(BIG_OUTPUT$IUCN)
BIG_OUTPUT$migration_f<-as.factor(BIG_OUTPUT$migration)
BIG_OUTPUT$territory_f<-as.factor(BIG_OUTPUT$territory)

summary(BIG_OUTPUT$Realm_f)

tapply(BIG_OUTPUT$Range_Ne, BIG_OUTPUT$Realm_f, summary)
tapply(BIG_OUTPUT$Var_Ne, BIG_OUTPUT$Realm_f, summary)


# remove the two species with no Ne information (APP-025 and B10K-DU-013-51) since keeping them will mess up sample sizes downstream
BIG_OUTPUT<-subset(BIG_OUTPUT, is.na(BIG_OUTPUT$Time_Ne_max) == FALSE)




sample_size_double_check<-setDT(All_PSMC_cut)[, count := uniqueN(B10k_ID), by = realm]

All_PSMC_cut<-subset(All_PSMC_cut, is.na(All_PSMC_cut$realm) == FALSE)
All_PSMC_cut$realm_f<-as.factor(All_PSMC_cut$realm)
setDT(All_PSMC_cut)[, count := uniqueN(B10k_ID), by = realm]

##Some time periods for Madagascar only have one value, so UCL and LCL can't be calculated
Mean_CI_All_Realms<-groupwiseMean(normalized_Ne ~ Round_Time + realm,
                                  data = All_PSMC_cut,
                                  conf = 0.95,
                                  digits = 3)


Mean_CI_All_Realms$Time_kya<-Mean_CI_All_Realms$Round_Time/1000

Mean_CI_Afro<-subset(Mean_CI_All_Realms, realm == "Afrotropical")
Normalized_NE_Smooth_Afro<-smooth.spline(x = Mean_CI_Afro$Time_kya, y = Mean_CI_Afro$Mean)
Mean_CI_Aus<-subset(Mean_CI_All_Realms, realm == "Australian")
Normalized_NE_Smooth_Aus<-smooth.spline(x = Mean_CI_Aus$Time_kya, y = Mean_CI_Aus$Mean)
Mean_CI_Mad<-subset(Mean_CI_All_Realms, realm == "Madagascan")
Normalized_NE_Smooth_Mad<-smooth.spline(x = Mean_CI_Mad$Time_kya, y = Mean_CI_Mad$Mean)
Mean_CI_Near<-subset(Mean_CI_All_Realms, realm == "Nearctic")
Normalized_NE_Smooth_Near<-smooth.spline(x = Mean_CI_Near$Time_kya, y = Mean_CI_Near$Mean)
Mean_CI_Neotro<-subset(Mean_CI_All_Realms, realm == "Neotropical")
Normalized_NE_Smooth_Neotro<-smooth.spline(x = Mean_CI_Neotro$Time_kya, y = Mean_CI_Neotro$Mean)
Mean_CI_Ocean<-subset(Mean_CI_All_Realms, realm == "Oceanina")
Normalized_NE_Smooth_Ocean<-smooth.spline(x = Mean_CI_Ocean$Time_kya, y = Mean_CI_Ocean$Mean)
Mean_CI_Orient<-subset(Mean_CI_All_Realms, realm == "Oriental")
Normalized_NE_Smooth_Orient<-smooth.spline(x = Mean_CI_Orient$Time_kya, y = Mean_CI_Orient$Mean)
Mean_CI_Pale<-subset(Mean_CI_All_Realms, realm == "Palearctic")
Normalized_NE_Smooth_Pale<-smooth.spline(x = Mean_CI_Pale$Time_kya, y = Mean_CI_Pale$Mean)
Mean_CI_Pan<-subset(Mean_CI_All_Realms, realm == "Panamanian")
Normalized_NE_Smooth_Pan<-smooth.spline(x = Mean_CI_Pan$Time_kya, y = Mean_CI_Pan$Mean)
Mean_CI_Saha<-subset(Mean_CI_All_Realms, realm == "Saharo-Arabian")
Normalized_NE_Smooth_Saha<-smooth.spline(x = Mean_CI_Saha$Time_kya, y = Mean_CI_Saha$Mean)
Mean_CI_Sino<-subset(Mean_CI_All_Realms, realm == "Sino-Japanese")
Normalized_NE_Smooth_Sino<-smooth.spline(x = Mean_CI_Sino$Time_kya, y = Mean_CI_Sino$Mean)





############################################
## put it all together 


### For Normalized Figure
Normalized_values_All_Realms_263_spp<-read.csv("Normalized_Ne_Full_Time_Period_Updated_PSMC_2022_05_09_equidistant_time_points.csv", header = TRUE)
Normalized_values_All_Realms_263_spp$B10k_ID_char<-as.character(Normalized_values_All_Realms_263_spp$ID)
Normalized_values_All_Realms_263_spp$Realm_char<-as.character(NA)


### add in realm to separate
          for (kkk in 1:nrow(Normalized_values_All_Realms_263_spp))  
          {
            focal_species_realm_fig<-Normalized_values_All_Realms_263_spp[kkk,"B10k_ID_char"]
            
            lookup_ID_realm<-which(BIG_OUTPUT$B10k_ID == focal_species_realm_fig)
            
            Normalized_values_All_Realms_263_spp[kkk, "Realm_char"]<-BIG_OUTPUT[lookup_ID_realm,"Realm"]
          
          }



########  function to get both the mean/SD per realm, and curves for each species in the realm

            Transpose_get_mean_per_realm_Normalize<-function(Focal_Realm)
            {
              
              Normalized_values_Focal_Realm<-subset(Normalized_values_All_Realms_263_spp, Realm_char == Focal_Realm)
              t_Normalized_values_Focal_Realm<-transpose(Normalized_values_Focal_Realm[,2:121])
              Num_Columns_realm<-ncol(t_Normalized_values_Focal_Realm)
              t_Normalized_values_Focal_Realm$mean<-apply(t_Normalized_values_Focal_Realm[1:Num_Columns_realm], 1, mean)
              t_Normalized_values_Focal_Realm$SD<-apply(t_Normalized_values_Focal_Realm[1:Num_Columns_realm], 1, sd)
              t_Normalized_values_Focal_Realm$time<-as.numeric(NA)
              t_Normalized_values_Focal_Realm[1:120,"time"]<-c(1:120)
              ## new addition Oct 2021
              t_Normalized_values_Focal_Realm$Log10_time<-as.numeric(NA)
              t_Normalized_values_Focal_Realm[1:120,"Log10_time"]<-log10(seq(30000,1000000, by = ((1000000 - 30000)/119)))
              
              return(t_Normalized_values_Focal_Realm)
              
            }
            
            t_Normalized_values_All_Realms_263_spp_Afro<-Transpose_get_mean_per_realm_Normalize("Afrotropical")
            t_Normalized_values_All_Realms_263_spp_Aust<-Transpose_get_mean_per_realm_Normalize("Australian")
            t_Normalized_values_All_Realms_263_spp_Mad<-Transpose_get_mean_per_realm_Normalize("Madagascan")
            t_Normalized_values_All_Realms_263_spp_Near<-Transpose_get_mean_per_realm_Normalize("Nearctic")
            t_Normalized_values_All_Realms_263_spp_Neo<-Transpose_get_mean_per_realm_Normalize("Neotropical")
            t_Normalized_values_All_Realms_263_spp_Ocean<-Transpose_get_mean_per_realm_Normalize("Oceanina")
            t_Normalized_values_All_Realms_263_spp_Orient<-Transpose_get_mean_per_realm_Normalize("Oriental")
            t_Normalized_values_All_Realms_263_spp_Pale<-Transpose_get_mean_per_realm_Normalize("Palearctic")
            t_Normalized_values_All_Realms_263_spp_Panam<-Transpose_get_mean_per_realm_Normalize("Panamanian")
            t_Normalized_values_All_Realms_263_spp_Saharo<-Transpose_get_mean_per_realm_Normalize("Saharo-Arabian")
            t_Normalized_values_All_Realms_263_spp_Sino<-Transpose_get_mean_per_realm_Normalize("Sino-Japanese")
            
            
            #For Global average
            t_Normalized_values_All_Realms_263_spp<-transpose(Normalized_values_All_Realms_263_spp[,2:121])
            Num_Columns_realm<-ncol(t_Normalized_values_All_Realms_263_spp)
            t_Normalized_values_All_Realms_263_spp$mean<-apply(t_Normalized_values_All_Realms_263_spp[1:Num_Columns_realm], 1, mean)
            t_Normalized_values_All_Realms_263_spp$SD<-apply(t_Normalized_values_All_Realms_263_spp[1:Num_Columns_realm], 1, sd)
            t_Normalized_values_All_Realms_263_spp$time<-as.numeric(NA)
            t_Normalized_values_All_Realms_263_spp[1:120,"time"]<-c(1:120)
            ## new addition Oct 2021
            t_Normalized_values_All_Realms_263_spp$Log10_time<-as.numeric(NA)
            t_Normalized_values_All_Realms_263_spp[1:120,"Log10_time"]<-log10(seq(30000,1000000, by = ((1000000 - 30000)/119)))
            

### Figure
labelz=c(expression(paste("3 x 10"^"4")),rep("",6),expression("10"^"5"),rep("",8),expression("10"^"6"))
at<-c(log10(seq(30000,90000, by = ((90000 - 30000)/6))), log10(seq(100000, 1000000, by = ((1000000 - 100000)/9))))

        par(mar=c(3,4.5,0.5,0.5))  
        plot(t_Normalized_values_All_Realms_263_spp_Afro$Log10_time,t_Normalized_values_All_Realms_263_spp_Afro$mean, type = "n", ylim = c(0,1),
             ann = F, axes = F)
        axis(side=2, at=c(0.0,0.2,0.4,0.6,0.8,1), las = 1,cex.axis = 1.15,cex.lab = 1.15)
        axis(1, at=at,label=labelz, cex.axis = 1.15,cex.lab = 1.15)
        
        num_individuals<-(ncol(t_Normalized_values_All_Realms_263_spp_Afro)-4)   ## don't forget last three columns are means etc
        
        for (mmm in 1:num_individuals)                 
        {
          lines(x=t_Normalized_values_All_Realms_263_spp_Afro[,"Log10_time"], y=t_Normalized_values_All_Realms_263_spp_Afro[,mmm], lwd = 3, col = 'light grey')
        }
        points(t_Normalized_values_All_Realms_263_spp_Afro$Log10_time,t_Normalized_values_All_Realms_263_spp_Afro$mean, pch = 20, col = 'dodgerblue')
        arrows(x0=t_Normalized_values_All_Realms_263_spp_Afro$Log10_time, y0=t_Normalized_values_All_Realms_263_spp_Afro$mean-t_Normalized_values_All_Realms_263_spp_Afro$SD,
               x1=t_Normalized_values_All_Realms_263_spp_Afro$Log10_time, y1=t_Normalized_values_All_Realms_263_spp_Afro$mean+t_Normalized_values_All_Realms_263_spp_Afro$SD,
               angle = 0, col = "dodgerblue")
        mtext(text = "Afrotropical (38)",side = 1, line = -12.3, adj = 0.1, cex = 1.01, outer = F)
        
        
        
        par(mar=c(3,4.5,0.5,0.5))  
        plot(t_Normalized_values_All_Realms_263_spp_Aust$Log10_time,t_Normalized_values_All_Realms_263_spp_Aust$mean, type = "n", ylim = c(0,1),
             ann = F, axes = F)
        axis(side=2, at=c(0.0,0.2,0.4,0.6,0.8,1), las = 1,cex.axis = 1.15,cex.lab = 1.15)
        axis(1, at=at,label=labelz, cex.axis = 1.15,cex.lab = 1.15)
        
        num_individuals<-(ncol(t_Normalized_values_All_Realms_263_spp_Aust)-4)   ## don't forget last three columns are means etc
        
        for (mmm in 1:num_individuals)                 
        {
          lines(x=t_Normalized_values_All_Realms_263_spp_Aust[,"Log10_time"], y=t_Normalized_values_All_Realms_263_spp_Aust[,mmm], lwd = 3, col = 'light grey')
        }
        points(t_Normalized_values_All_Realms_263_spp_Aust$Log10_time,t_Normalized_values_All_Realms_263_spp_Aust$mean, pch = 20, col = '#1F78B4')
        arrows(x0=t_Normalized_values_All_Realms_263_spp_Aust$Log10_time, y0=t_Normalized_values_All_Realms_263_spp_Aust$mean-t_Normalized_values_All_Realms_263_spp_Aust$SD,
               x1=t_Normalized_values_All_Realms_263_spp_Aust$Log10_time, y1=t_Normalized_values_All_Realms_263_spp_Aust$mean+t_Normalized_values_All_Realms_263_spp_Aust$SD,
               angle = 0, col = "#1F78B4")
        mtext(text = "Australian (28)",side = 1, line = -12.3, adj = 0.1, cex = 1.01, outer = F)
        
        
        
        
        par(mar=c(3,4.5,0.5,0.5)) 
        plot(t_Normalized_values_All_Realms_263_spp_Mad$Log10_time,t_Normalized_values_All_Realms_263_spp_Mad$mean, type="n", ylim = c(0,1),
             ann = F, axes = F)
        axis(side=2, at=c(0.0,0.2,0.4,0.6,0.8,1), las = 1,cex.axis = 1.15,cex.lab = 1.15)
        axis(1, at=at,label=labelz, cex.axis = 1.15,cex.lab = 1.15)
        
        
        num_individuals<-(ncol(t_Normalized_values_All_Realms_263_spp_Mad)-4)   ## don't forget last three columns are means etc
        
        for (mmm in 1:num_individuals)                 
        {
          lines(x=t_Normalized_values_All_Realms_263_spp_Mad[,"Log10_time"], y=t_Normalized_values_All_Realms_263_spp_Mad[,mmm], lwd = 3, col = 'light grey')
        }
        points(t_Normalized_values_All_Realms_263_spp_Mad$Log10_time,t_Normalized_values_All_Realms_263_spp_Mad$mean, pch = 20, col = 'forestgreen')
        arrows(x0=t_Normalized_values_All_Realms_263_spp_Mad$Log10_time, y0=t_Normalized_values_All_Realms_263_spp_Mad$mean-t_Normalized_values_All_Realms_263_spp_Mad$SD,
               x1=t_Normalized_values_All_Realms_263_spp_Mad$Log10_time, y1=t_Normalized_values_All_Realms_263_spp_Mad$mean+t_Normalized_values_All_Realms_263_spp_Mad$SD,
               angle = 0, col = "darkolivegreen4")
        mtext(text = "Madagascan (4)",side = 1, line = -12.3, adj = 0.1, cex = 1.01, outer = F)
        
        
        
        ######### Next row
        par(mar=c(3,4.5,0.5,0.5)) 
        plot(t_Normalized_values_All_Realms_263_spp_Near$Log10_time,t_Normalized_values_All_Realms_263_spp_Near$mean, type="n", ylim = c(0,1),
             ann = F, axes = F)
        axis(side=2, at=c(0.0,0.2,0.4,0.6,0.8,1), las = 1,cex.axis = 1.15,cex.lab = 1.15)
        axis(1, at=at,label=labelz, cex.axis = 1.15,cex.lab = 1.15)
        
        
        num_individuals<-(ncol(t_Normalized_values_All_Realms_263_spp_Near)-4)   ## don't forget last three columns are means etc
        
        for (mmm in 1:num_individuals)                 
        {
          lines(x=t_Normalized_values_All_Realms_263_spp_Near[,"Log10_time"], y=t_Normalized_values_All_Realms_263_spp_Near[,mmm], lwd = 3, col = 'light grey')
        }
        points(t_Normalized_values_All_Realms_263_spp_Near$Log10_time,t_Normalized_values_All_Realms_263_spp_Near$mean, pch = 20, col = 'forestgreen')
        arrows(x0=t_Normalized_values_All_Realms_263_spp_Near$Log10_time, y0=t_Normalized_values_All_Realms_263_spp_Near$mean-t_Normalized_values_All_Realms_263_spp_Near$SD,
               x1=t_Normalized_values_All_Realms_263_spp_Near$Log10_time, y1=t_Normalized_values_All_Realms_263_spp_Near$mean+t_Normalized_values_All_Realms_263_spp_Near$SD,
               angle = 0, col = "forestgreen")
        mtext(text = "Nearctic (40)",side = 1, line = -12.3, adj = 0.1, cex = 1.01, outer = F)
        
        
        
        par(mar=c(3,4.5,0.5,0.5))  
        plot(t_Normalized_values_All_Realms_263_spp_Neo$Log10_time,t_Normalized_values_All_Realms_263_spp_Neo$mean, type="n", ylim = c(0,1),
             ann = F, axes = F)
        axis(side=2, at=c(0.0,0.2,0.4,0.6,0.8,1), las = 1,cex.axis = 1.15,cex.lab = 1.15)
        axis(1, at=at,label=labelz, cex.axis = 1.15,cex.lab = 1.15)
        
        num_individuals<-(ncol(t_Normalized_values_All_Realms_263_spp_Neo)-4)   ## don't forget last three columns are means etc
        
        for (mmm in 1:num_individuals)                 
        {
          lines(x=t_Normalized_values_All_Realms_263_spp_Neo[,"Log10_time"], y=t_Normalized_values_All_Realms_263_spp_Neo[,mmm], lwd = 3, col = 'light grey')
        }
        points(t_Normalized_values_All_Realms_263_spp_Neo$Log10_time,t_Normalized_values_All_Realms_263_spp_Neo$mean, pch = 20, col = 'palevioletred3')
        arrows(x0=t_Normalized_values_All_Realms_263_spp_Neo$Log10_time, y0=t_Normalized_values_All_Realms_263_spp_Neo$mean-t_Normalized_values_All_Realms_263_spp_Neo$SD,
               x1=t_Normalized_values_All_Realms_263_spp_Neo$Log10_time, y1=t_Normalized_values_All_Realms_263_spp_Neo$mean+t_Normalized_values_All_Realms_263_spp_Neo$SD,
               angle = 0, col = "palevioletred3")
        mtext(text = "Neotropical (60)",side = 1, line = -12.3, adj = 0.1, cex = 1.01, outer = F)
        
        
        
        
        par(mar=c(3,4.5,0.5,0.5))  
        plot(t_Normalized_values_All_Realms_263_spp_Ocean$Log10_time,t_Normalized_values_All_Realms_263_spp_Ocean$mean, type="n", ylim = c(0,1),
             ann = F, axes = F)
        axis(side=2, at=c(0.0,0.2,0.4,0.6,0.8,1), las = 1,cex.axis = 1.15,cex.lab = 1.15)
        axis(1, at=at,label=labelz, cex.axis = 1.15,cex.lab = 1.15)
        
        num_individuals<-(ncol(t_Normalized_values_All_Realms_263_spp_Ocean)-4)   ## don't forget last three columns are means etc
        
        for (mmm in 1:num_individuals)                 
        {
          lines(x=t_Normalized_values_All_Realms_263_spp_Ocean[,"Log10_time"], y=t_Normalized_values_All_Realms_263_spp_Ocean[,mmm], lwd = 3, col = 'light grey')
        }
        points(t_Normalized_values_All_Realms_263_spp_Ocean$Log10_time,t_Normalized_values_All_Realms_263_spp_Ocean$mean, pch = 20, col = '#E31A1C')
        arrows(x0=t_Normalized_values_All_Realms_263_spp_Ocean$Log10_time, y0=t_Normalized_values_All_Realms_263_spp_Ocean$mean-t_Normalized_values_All_Realms_263_spp_Ocean$SD,
               x1=t_Normalized_values_All_Realms_263_spp_Ocean$Log10_time, y1=t_Normalized_values_All_Realms_263_spp_Ocean$mean+t_Normalized_values_All_Realms_263_spp_Ocean$SD,
               angle = 0, col = "#E31A1C")
        mtext(text = "Oceanian (19)",side = 1, line = -12.3, adj = 0.1, cex = 1.01, outer = F)
        
        
        #### next row, y axis here
        par(mar=c(3,4.5,0.5,0.5))  
        plot(t_Normalized_values_All_Realms_263_spp_Orient$Log10_time,t_Normalized_values_All_Realms_263_spp_Orient$mean, type="n", ylim = c(0,1),
             ann = F, axes = F)
        axis(side=2, at=c(0.0,0.2,0.4,0.6,0.8,1), las = 1,cex.axis = 1.15,cex.lab = 1.15)
        axis(1, at=at,label=labelz, cex.axis = 1.15,cex.lab = 1.15)
        
        num_individuals<-(ncol(t_Normalized_values_All_Realms_263_spp_Orient)-4)   ## don't forget last three columns are means etc
        
        for (mmm in 1:num_individuals)                 
        {
          lines(x=t_Normalized_values_All_Realms_263_spp_Orient[,"Log10_time"], y=t_Normalized_values_All_Realms_263_spp_Orient[,mmm], lwd = 3, col = 'light grey')
        }
        points(t_Normalized_values_All_Realms_263_spp_Orient$Log10_time,t_Normalized_values_All_Realms_263_spp_Orient$mean, pch = 20, col = 'lightgoldenrod3')
        arrows(x0=t_Normalized_values_All_Realms_263_spp_Orient$Log10_time, y0=t_Normalized_values_All_Realms_263_spp_Orient$mean-t_Normalized_values_All_Realms_263_spp_Orient$SD,
               x1=t_Normalized_values_All_Realms_263_spp_Orient$Log10_time, y1=t_Normalized_values_All_Realms_263_spp_Orient$mean+t_Normalized_values_All_Realms_263_spp_Orient$SD,
               angle = 0, col = "lightgoldenrod3")
        mtext(text = "Oriental (12)",side = 1, line = -12.3, adj = 0.1, cex = 1.01, outer = F)
        mtext(expression(paste(Normalized~italic("N")[italic("e")])), side = 2, line = 2.8, adj = 10, cex = 1.3, outer = F)
         
        
        
        
        par(mar=c(3,4.5,0.5,0.5))  
        plot(t_Normalized_values_All_Realms_263_spp_Pale$Log10_time,t_Normalized_values_All_Realms_263_spp_Pale$mean, type="n", ylim = c(0,1),
             ann = F, axes = F)
        axis(side=2, at=c(0.0,0.2,0.4,0.6,0.8,1), las = 1,cex.axis = 1.15,cex.lab = 1.15)
        axis(1, at=at,label=labelz, cex.axis = 1.15,cex.lab = 1.15)
        
        
        num_individuals<-(ncol(t_Normalized_values_All_Realms_263_spp_Pale)-4)   ## don't forget last three columns are means etc
        
        for (mmm in 1:num_individuals)                 
        {
          lines(x=t_Normalized_values_All_Realms_263_spp_Pale[,"Log10_time"], y=t_Normalized_values_All_Realms_263_spp_Pale[,mmm], lwd = 3, col = 'light grey')
        }
        points(t_Normalized_values_All_Realms_263_spp_Pale$Log10_time,t_Normalized_values_All_Realms_263_spp_Pale$mean, pch = 20, col = '#FF7F00')
        arrows(x0=t_Normalized_values_All_Realms_263_spp_Pale$Log10_time, y0=t_Normalized_values_All_Realms_263_spp_Pale$mean-t_Normalized_values_All_Realms_263_spp_Pale$SD,
               x1=t_Normalized_values_All_Realms_263_spp_Pale$Log10_time, y1=t_Normalized_values_All_Realms_263_spp_Pale$mean+t_Normalized_values_All_Realms_263_spp_Pale$SD,
               angle = 0, col = "#FF7F00")
        mtext(text = "Palearctic (28)",side = 1, line = -12.3, adj = 0.1, cex = 1.01, outer = F)
        
        
        
        
        par(mar=c(3,4.5,0.5,0.5))  
        plot(t_Normalized_values_All_Realms_263_spp_Panam$Log10_time,t_Normalized_values_All_Realms_263_spp_Panam$mean, type="n", ylim = c(0,1),
             ann = F, axes = F)
        axis(side=2, at=c(0.0,0.2,0.4,0.6,0.8,1), las = 1,cex.axis = 1.15,cex.lab = 1.15)
        axis(1, at=at,label=labelz, cex.axis = 1.15,cex.lab = 1.15)
        
        num_individuals<-(ncol(t_Normalized_values_All_Realms_263_spp_Panam)-4)   ## don't forget last three columns are means etc
        
        for (mmm in 1:num_individuals)                 
        {
          lines(x=t_Normalized_values_All_Realms_263_spp_Panam[,"Log10_time"], y=t_Normalized_values_All_Realms_263_spp_Panam[,mmm], lwd = 3, col = 'light grey')
        }
        points(t_Normalized_values_All_Realms_263_spp_Panam$Log10_time,t_Normalized_values_All_Realms_263_spp_Panam$mean, pch = 20, col = 'mediumorchid2')
        arrows(x0=t_Normalized_values_All_Realms_263_spp_Panam$Log10_time, y0=t_Normalized_values_All_Realms_263_spp_Panam$mean-t_Normalized_values_All_Realms_263_spp_Panam$SD,
               x1=t_Normalized_values_All_Realms_263_spp_Panam$Log10_time, y1=t_Normalized_values_All_Realms_263_spp_Panam$mean+t_Normalized_values_All_Realms_263_spp_Panam$SD,
               angle = 0, col = "mediumorchid2")
        mtext(text = "Panamanian (12)",side = 1, line = -12.3, adj = 0.1, cex = 1.01, outer = F)
        
        
        
        ### next row
        par(mar=c(3,4.5,0.5,0.5))  
        plot(t_Normalized_values_All_Realms_263_spp_Saharo$Log10_time,t_Normalized_values_All_Realms_263_spp_Saharo$mean, type="n", ylim = c(0,1),
             ann = F, axes = F)
        axis(side=2, at=c(0.0,0.2,0.4,0.6,0.8,1), las = 1,cex.axis = 1.15,cex.lab = 1.15)
        axis(1, at=at,label=labelz, cex.axis = 1.15,cex.lab = 1.15)
        
        
        num_individuals<-(ncol(t_Normalized_values_All_Realms_263_spp_Saharo)-4)   ## don't forget last three columns are means etc
        
        for (mmm in 1:num_individuals)                 
        {
          lines(x=t_Normalized_values_All_Realms_263_spp_Saharo[,"Log10_time"], y=t_Normalized_values_All_Realms_263_spp_Saharo[,mmm], lwd = 3, col = 'light grey')
        }
        points(t_Normalized_values_All_Realms_263_spp_Saharo$Log10_time,t_Normalized_values_All_Realms_263_spp_Saharo$mean, pch = 20, col = 'purple4')
        arrows(x0=t_Normalized_values_All_Realms_263_spp_Saharo$Log10_time, y0=t_Normalized_values_All_Realms_263_spp_Saharo$mean-t_Normalized_values_All_Realms_263_spp_Saharo$SD,
               x1=t_Normalized_values_All_Realms_263_spp_Saharo$Log10_time, y1=t_Normalized_values_All_Realms_263_spp_Saharo$mean+t_Normalized_values_All_Realms_263_spp_Saharo$SD,
               angle = 0, col = "purple4")
        mtext(text = "Saharo-Arabian (3)",side = 1, line = -12.3, adj = 0.1, cex = 1.01, outer = F)
        
        
        
        
        par(mar=c(3,4.5,0.5,0.5))  
        plot(t_Normalized_values_All_Realms_263_spp_Sino$Log10_time,t_Normalized_values_All_Realms_263_spp_Sino$mean, type="n", ylim = c(0,1),
             ann = F, axes = F)
        axis(side=2, at=c(0.0,0.2,0.4,0.6,0.8,1), las = 1,cex.axis = 1.15,cex.lab = 1.15)
        axis(1, at=at,label=labelz, cex.axis = 1.15,cex.lab = 1.15)
        
        
        num_individuals<-(ncol(t_Normalized_values_All_Realms_263_spp_Sino)-4)   ## don't forget last three columns are means etc
        
        for (mmm in 1:num_individuals)                 
        {
          lines(x=t_Normalized_values_All_Realms_263_spp_Sino[,"Log10_time"], y=t_Normalized_values_All_Realms_263_spp_Sino[,mmm], lwd = 3, col = 'light grey')
        }
        points(t_Normalized_values_All_Realms_263_spp_Sino$Log10_time,t_Normalized_values_All_Realms_263_spp_Sino$mean, pch = 20, col = '#B15928')
        arrows(x0=t_Normalized_values_All_Realms_263_spp_Sino$Log10_time, y0=t_Normalized_values_All_Realms_263_spp_Sino$mean-t_Normalized_values_All_Realms_263_spp_Sino$SD,
               x1=t_Normalized_values_All_Realms_263_spp_Sino$Log10_time, y1=t_Normalized_values_All_Realms_263_spp_Sino$mean+t_Normalized_values_All_Realms_263_spp_Sino$SD,
               angle = 0, col = "#B15928")
        mtext(text = "Sino-Japanese (10)",side = 1, line = -12.3, adj = 0.1, cex = 1.01, outer = F)
        mtext("Time Since Present (kya)", side = 1, line = 2.8, adj = 0.5, cex = 1.2, outer = F)
        
        
        
        par(mar=c(3,4.5,0.5,0.5))  
        plot(t_Normalized_values_All_Realms_263_spp$Log10_time,t_Normalized_values_All_Realms_263_spp$mean, pch = 20, col = 'grey36', ylim = c(0,1),
             ann = F, axes = F)
        axis(side=2, at=c(0.0,0.2,0.4,0.6,0.8,1), las = 1,cex.axis = 1.15,cex.lab = 1.15)
        axis(1, at=at,label=labelz, cex.axis = 1.15,cex.lab = 1.15)
        
        num_individuals<-(ncol(t_Normalized_values_All_Realms_263_spp)-4)   ## don't forget last three columns are means etc
        
        points(t_Normalized_values_All_Realms_263_spp$Log10_time,t_Normalized_values_All_Realms_263_spp$mean, pch = 20, col = 'grey36')
        arrows(x0=t_Normalized_values_All_Realms_263_spp$Log10_time, y0=t_Normalized_values_All_Realms_263_spp$mean-t_Normalized_values_All_Realms_263_spp$SD,
               x1=t_Normalized_values_All_Realms_263_spp$Log10_time, y1=t_Normalized_values_All_Realms_263_spp$mean+t_Normalized_values_All_Realms_263_spp$SD,
               angle = 0, col = "grey36")
        mtext(text = "Global Average (263)",side = 1, line = -12.3, adj = 0.1, cex = 1.01, outer = F)



#####################################################3
#####################################################
### End






        
        
        

#####################################################3
#####################################################
### Trait analysis    
    
library("factoextra")
library("missForest")

   ### There are missing data for a number of fields - impute missing values via random forest 
  trait_data_random_forest <- BIG_OUTPUT[, c("B10k_ID","mass_unsexed","clutch_size","egg_mass","maximum_longevity_y",
                                             "mortality_both","inc_duration","min_elevation","max_elevation","range_size_km2","tarsus_length","wing_chord",
                                             "kipps_distance","hand_wing_index", "bill_total_culmen", "bill_nares","bill_width", "bill_depth", "b_brain_size", 
                                             "territory_f","migration_f",
                                             "Conservation_f", "Realm_f","Order","Pass_NonPass_f","Classify","New_Clusters_MAY_2022_f")]
              
          ## count missing values per column (trait)
          sapply(trait_data_random_forest, function(x) sum(is.na(x)))
  
  
          traits_imp<-missForest(trait_data_random_forest[,2:19]) # make sure to leave the B10k Id out of the process or you get an error, also leave out the things that can't/shouldn't be imputed (realm, conservation status)               
          traits_imp$ximp                      
          traits_imp$OOBerror    # looks good                  
          
          
          Trait_Dataset<-as.data.frame(traits_imp$ximp)
          Trait_Dataset_THIS_ONE<-cbind(BIG_OUTPUT$B10k_ID,Trait_Dataset,BIG_OUTPUT$Conservation_f,BIG_OUTPUT$Realm_f,BIG_OUTPUT$Order,BIG_OUTPUT$Pass_NonPass_f,BIG_OUTPUT$Classify,BIG_OUTPUT$New_Clusters_MAY_2022_f)
          ### reset column names for simplicity
          colnames(Trait_Dataset_THIS_ONE)<-c("B10k_ID","mass_unsexed","clutch_size","egg_mass","maximum_longevity_y",
                                              "mortality_both","inc_duration","min_elevation","max_elevation","range_size_km2","tarsus_length","wing_chord",
                                              "kipps_distance","hand_wing_index", "bill_total_culmen", "bill_nares","bill_width", "bill_depth", "b_brain_size", 
                                              "territory_f","migration_f",
                                              "Conservation_f", "Realm_f","Order","Pass_NonPass_f","Classify","New_Clusters_MAY_2022_f")
          Trait_Dataset_THIS_ONE$B10k_ID_char<-as.character(Trait_Dataset_THIS_ONE$B10k_ID)
          Trait_Dataset_THIS_ONE$Classify_f<-as.factor(Trait_Dataset_THIS_ONE$Classify)
    
          
          
   
### relating traits to demographic responses to warming/cooling          

        #### summarized pearson correlation coefficients between temperature change and demographic change during cooling and warming 
       Correlations_Warm_Cool<-read.csv("Pearson_correlation_Summary_Stats_No_Filter_May_2022.csv", header = T)
       Correlations_Warm_Cool$B10k_ID_char<-as.character(Correlations_Warm_Cool$ID)  
       Correlations_Warm_Cool$Overall_response_type_char<-as.character(Correlations_Warm_Cool$Overall_response_type)  
       
              
       Trait_Dataset_THIS_ONE$corr_Warm_1<-as.numeric(NA)
       Trait_Dataset_THIS_ONE$corr_Cool_1<-as.numeric(NA)
       Trait_Dataset_THIS_ONE$overall_response_type<-as.character(NA)
     
      
       for (ggg in 1:nrow(Trait_Dataset_THIS_ONE))
       {
         focal_species_B10k<-Trait_Dataset_THIS_ONE[ggg,"B10k_ID_char"]
         
         lookup_species_warm_cool<-which(Correlations_Warm_Cool$B10k_ID_char == focal_species_B10k)
         if (length(lookup_species_warm_cool) == 1)    
         {
           Trait_Dataset_THIS_ONE[ggg,"corr_Warm_1"]<-Correlations_Warm_Cool[lookup_species_warm_cool,"IN1_coef"]
           Trait_Dataset_THIS_ONE[ggg,"corr_Cool_1"]<-Correlations_Warm_Cool[lookup_species_warm_cool,"DC1_coef"]
           Trait_Dataset_THIS_ONE[ggg,"overall_response_type"]<-Correlations_Warm_Cool[lookup_species_warm_cool,"Overall_response_type_char"]
           
          }
       }
          
       
    ### Remove species where we don't have accurate demo data across both warming and cooling periods
       
       completeFun <- function(data, desiredCols) {
         completeVec <- complete.cases(data[, desiredCols])
         return(data[completeVec, ])
       }
       
       
       Trait_Dataset_THIS_ONE_Complete_Cases<-completeFun(Trait_Dataset_THIS_ONE, "corr_Warm_1")
              Trait_Dataset_THIS_ONE_Complete_Cases$ratio_brain_body<-as.numeric(Trait_Dataset_THIS_ONE_Complete_Cases$b_brain_size/Trait_Dataset_THIS_ONE_Complete_Cases$mass_unsexed)
       
           
       ### First thing is remove all the NN  entries. These are the cases where we don't have high enough confidence to assign a positive/negative correlation
           Trait_Dataset_THIS_ONE_Complete_Cases<-subset(Trait_Dataset_THIS_ONE_Complete_Cases, overall_response_type != "NN")    
            ## reduces dataset to 215
       
       
       
 
    
###  add in a column for Passerine/NonPasserine
          
Trait_Dataset_THIS_ONE_Complete_Cases$Passerine<-as.character(NA)

for (hhh in 1:nrow(Trait_Dataset_THIS_ONE_Complete_Cases))
{
  if (Trait_Dataset_THIS_ONE_Complete_Cases[hhh, "Order"] == "Passeriformes")
  {
    Trait_Dataset_THIS_ONE_Complete_Cases[hhh, "Passerine"]<-"Passerine"
  } else
  {
    Trait_Dataset_THIS_ONE_Complete_Cases[hhh, "Passerine"]<-"Nonpasserine"
  }

}

Trait_Dataset_THIS_ONE_Complete_Cases$Passerine_f<-as.factor(Trait_Dataset_THIS_ONE_Complete_Cases$Passerine)
summary(Trait_Dataset_THIS_ONE_Complete_Cases$Passerine_f)
                          



                
                
                ########################################################                             
                ## Mixed models - warming
                 
                scale_warm<-scale(Trait_Dataset_THIS_ONE_Complete_Cases$corr_Warm_1)
                
                warm_1_Model<-lmer(scale(corr_Warm_1) ~ (1|Passerine_f) 
                                            + scale(mass_unsexed) + scale(ratio_brain_body) 
                                            + scale(tarsus_length) + scale (bill_total_culmen)
                                            + scale(egg_mass) + scale(inc_duration) + scale(clutch_size)
                                            + scale(hand_wing_index),
                                            na.action = na.fail, 
                                            data = Trait_Dataset_THIS_ONE_Complete_Cases)
                
                
                summary(warm_1_Model)     
                r.squaredGLMM(warm_1_Model)
                
                allmodels_warm1_Model<-dredge(warm_1_Model)
                sub_allmodels_warm1_Model<-subset(allmodels_warm1_Model, delta<5)   
                
                ##### Model Average
                Model_Average_warm1<- model.avg(sub_allmodels_warm1_Model, fit = TRUE) # givesestimates, etc of averaged full and subset of all model combos
                summary(Model_Average_warm1)
                confint(Model_Average_warm1)
                names(Model_Average_warm1)
                
                
                
                
                
                
                
                ########################################################                             
                ## Mixed models - Cooling
                scale_cool<-scale(Trait_Dataset_THIS_ONE_Complete_Cases$corr_Cool_1)
                
                Cool1_Model<-lmer(scale(corr_Cool_1) ~ (1|Passerine_f) 
                                  + scale(mass_unsexed) + scale(ratio_brain_body) 
                                  + scale(tarsus_length) 
                                  + scale (bill_total_culmen)
                                  + scale(egg_mass) + scale(inc_duration) + scale(clutch_size)
                                  + scale(hand_wing_index),
                                   na.action = na.fail, 
                                   data = Trait_Dataset_THIS_ONE_Complete_Cases)
                
                summary(Cool1_Model)     
                r.squaredGLMM(Cool1_Model)
                
                allmodels_Cool1_Model<-dredge(Cool1_Model)
                sub_allmodels_Cool1_Model<-subset(allmodels_Cool1_Model, delta<5)  
                
                ##### Model Average
                Model_Average_cool1<- model.avg(sub_allmodels_Cool1_Model, fit = TRUE) # givesestimates, etc of averaged full and subset of all model combos
                summary(Model_Average_cool1)
                confint(Model_Average_cool1)
                names(Model_Average_cool1)





####Phylogenetic Path Analysis####
                
                models <- define_model_set(
                  aonea = c(RS ~ HWI, RS ~ BM, RS ~ EM),
                  aoneb = c(RS ~ HWI, RS ~ BILL, RS ~ EM),
                  aonec = c(RS ~ HWI, RS ~ BM, RS ~ BILL, RS ~ EM),
                  atwoa = c(RS ~ HWI, RS ~ BM, RS ~ ID),
                  atwob = c(RS ~ HWI, RS ~ BILL, RS ~ ID),
                  atwoc = c(RS ~ HWI, RS ~ BM, RS ~ BILL, RS ~ ID),
                  athreea = c(RS ~ HWI, RS ~ BM, RS ~ CS),
                  athreeb = c(RS ~ HWI, RS ~ BILL, RS ~ CS),
                  athreec = c(RS ~ HWI, RS ~ BM, RS ~ BILL, RS ~ CS),
                  afoura = c(RS ~ HWI, RS ~ BM, RS ~ EM, RS ~ ID),
                  afourb = c(RS ~ HWI, RS ~ BILL, RS ~ EM, RS ~ ID),
                  afourc = c(RS ~ HWI, RS ~ BM, RS ~ BILL, RS ~ EM, RS ~ ID),
                  afivea = c(RS ~ HWI, RS ~ BM, RS ~ EM, RS ~ CS),
                  afiveb = c(RS ~ HWI, RS ~ BILL, RS ~ EM, RS ~ CS),
                  afivec = c(RS ~ HWI, RS ~ BM, RS ~ BILL,RS ~ EM, RS ~ CS),
                  asixa = c(RS ~ HWI, RS ~ BM, RS ~ ID, RS ~ CS),
                  asixb = c(RS ~ HWI, RS ~ BILL, RS ~ ID, RS ~ CS),
                  asixc = c(RS ~ HWI, RS ~ BM, RS ~ BILL, RS ~ ID, RS ~ CS),
                  asevena = c(RS ~ HWI, RS ~ BM, RS ~ EM, RS ~ ID, RS ~ CS),
                  asevenb = c(RS ~ HWI, RS ~ BILL, RS ~ EM, RS ~ ID, RS ~ CS),
                  asevenc = c(RS ~ HWI, RS ~ BM, RS ~ BILL, RS ~ EM, RS ~ ID, RS ~ CS),
                  bonea = c(RS ~ HWI, RS ~ BM, RS ~ EM, EM ~ BM),
                  boneb = c(RS ~ HWI, RS ~ BILL, RS ~ EM, EM ~ BM),
                  bonec = c(RS ~ HWI, RS ~ BM, RS ~ BILL, RS ~ EM, EM ~ BM),
                  btwoa = c(RS ~ HWI, RS ~ BM, RS ~ ID, EM ~ BM),
                  btwob = c(RS ~ HWI, RS ~ BILL, RS ~ ID, EM ~ BM),
                  btwoc = c(RS ~ HWI, RS ~ BM, RS ~ BILL, RS ~ ID, EM ~ BM),
                  bthreea = c(RS ~ HWI, RS ~ BM, RS ~ CS, EM ~ BM),
                  bthreeb = c(RS ~ HWI, RS ~ BILL, RS ~ CS, EM ~ BM),
                  bthreec = c(RS ~ HWI, RS ~ BM, RS ~ BILL, RS ~ CS, EM ~ BM),
                  bfoura = c(RS ~ HWI, RS ~ BM, RS ~ EM, RS ~ ID, EM ~ BM),
                  bfourb = c(RS ~ HWI, RS ~ BILL, RS ~ EM, RS ~ ID, EM ~ BM),
                  bfourc = c(RS ~ HWI, RS ~ BM, RS ~ BILL, RS ~ EM, RS ~ ID, EM ~ BM),
                  bfivea = c(RS ~ HWI, RS ~ BM, RS ~ EM, RS ~ CS, EM ~ BM),
                  bfiveb = c(RS ~ HWI, RS ~ BILL, RS ~ EM, RS ~ CS, EM ~ BM),
                  bfivec = c(RS ~ HWI, RS ~ BM, RS ~ BILL, RS ~ EM, RS ~ CS, EM ~ BM),
                  bsixa = c(RS ~ HWI, RS ~ BM, RS ~ ID, RS ~ CS, EM ~ BM),
                  bsixb = c(RS ~ HWI, RS ~ BILL, RS ~ ID, RS ~ CS, EM ~ BM),
                  bsixc = c(RS ~ HWI, RS ~ BM, RS ~ BILL, RS ~ ID, RS ~ CS, EM ~ BM),
                  bsevena = c(RS ~ HWI, RS ~ BM, RS ~ EM, RS ~ ID, RS ~ CS, EM ~ BM),
                  bsevenb = c(RS ~ HWI, RS ~ BILL, RS ~ EM, RS ~ ID, RS ~ CS, EM ~ BM),
                  bsevenc = c(RS ~ HWI, RS ~ BM, RS ~ BILL, RS ~ EM, RS ~ ID, RS ~ CS, EM ~ BM),
                  conea = c(RS ~ HWI, RS ~ BM, RS ~ EM, HWI ~ BM, HWI ~ BILL),
                  coneb = c(RS ~ HWI, RS ~ BILL, RS ~ EM, HWI ~ BM, HWI ~ BILL),
                  conec = c(RS ~ HWI, RS ~ BM, RS ~ BILL, RS ~ EM, HWI ~ BM, HWI ~ BILL),
                  ctwoa = c(RS ~ HWI, RS ~ BM, RS ~ ID, HWI ~ BM, HWI ~ BILL),
                  ctwob = c(RS ~ HWI, RS ~ BILL, RS ~ ID, HWI ~ BM, HWI ~ BILL),
                  ctwoc = c(RS ~ HWI, RS ~ BM, RS ~ BILL, RS ~ ID, HWI ~ BM, HWI ~ BILL),
                  cthreea = c(RS ~ HWI, RS ~ BM, RS ~ CS, HWI ~ BM, HWI ~ BILL),
                  cthreeb = c(RS ~ HWI, RS ~ BILL, RS ~ CS, HWI ~ BM, HWI ~ BILL),
                  cthreec = c(RS ~ HWI, RS ~ BM, RS ~ BILL, RS ~ CS, HWI ~ BM, HWI ~ BILL),
                  cfoura = c(RS ~ HWI, RS ~ BM, RS ~ EM, RS ~ ID, HWI ~ BM, HWI ~ BILL),
                  cfourb = c(RS ~ HWI, RS ~ BILL, RS ~ EM, RS ~ ID, HWI ~ BM, HWI ~ BILL),
                  cfourc = c(RS ~ HWI, RS ~ BM, RS ~ BILL, RS ~ EM, RS ~ ID, HWI ~ BM, HWI ~ BILL),
                  cfivea = c(RS ~ HWI, RS ~ BM, RS ~ EM, RS ~ CS, HWI ~ BM, HWI ~ BILL),
                  cfiveb = c(RS ~ HWI, RS ~ BILL, RS ~ EM, RS ~ CS, HWI ~ BM, HWI ~ BILL),
                  cfivec = c(RS ~ HWI, RS ~ BM, RS ~ BILL, RS ~ EM, RS ~ CS, HWI ~ BM, HWI ~ BILL),
                  csixa = c(RS ~ HWI, RS ~ BM, RS ~ ID, RS ~ CS, HWI ~ BM, HWI ~ BILL),
                  csixb = c(RS ~ HWI, RS ~ BILL, RS ~ ID, RS ~ CS, HWI ~ BM, HWI ~ BILL),
                  csixc = c(RS ~ HWI, RS ~ BM, RS ~ BILL, RS ~ ID, RS ~ CS, HWI ~ BM, HWI ~ BILL),
                  csevena = c(RS ~ HWI, RS ~ BM, RS ~ EM, RS ~ ID, RS ~ CS, HWI ~ BM, HWI ~ BILL),
                  csevenb = c(RS ~ HWI, RS ~ BILL, RS ~ EM, RS ~ ID, RS ~ CS, HWI ~ BM, HWI ~ BILL),
                  csevenc = c(RS ~ HWI, RS ~ BM, RS ~ BILL, RS ~ EM, RS ~ ID, RS ~ CS, HWI ~ BM, HWI ~ BILL),
                  donea = c(RS ~ HWI, RS ~ BM, RS ~ EM, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  doneb = c(RS ~ HWI, RS ~ BILL, RS ~ EM, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  donec = c(RS ~ HWI, RS ~ BM, RS ~ BILL, RS ~ EM, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  dtwoa = c(RS ~ HWI, RS ~ BM, RS ~ ID, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  dtwob = c(RS ~ HWI, RS ~ BILL, RS ~ ID, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  dtwoc = c(RS ~ HWI, RS ~ BM, RS ~ BILL, RS ~ ID, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  dthreea = c(RS ~ HWI, RS ~ BM, RS ~ CS, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  dthreeb = c(RS ~ HWI, RS ~ BILL, RS ~ CS, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  dthreec = c(RS ~ HWI, RS ~ BM, RS ~ BILL, RS ~ CS, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  dfoura = c(RS ~ HWI, RS ~ BM, RS ~ EM, RS ~ ID, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  dfourb = c(RS ~ HWI, RS ~ BILL, RS ~ EM, RS ~ ID, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  dfourc = c(RS ~ HWI, RS ~ BM, RS ~ BILL, RS ~ EM, RS ~ ID, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  dfivea = c(RS ~ HWI, RS ~ BM, RS ~ EM, RS ~ CS, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  dfiveb = c(RS ~ HWI, RS ~ BILL, RS ~ EM, RS ~ CS, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  dfivec = c(RS ~ HWI, RS ~ BM, RS ~ BILL, RS ~ EM, RS ~ CS, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  dsixa = c(RS ~ HWI, RS ~ BM, RS ~ ID, RS ~ CS, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  dsixb = c(RS ~ HWI, RS ~ BILL, RS ~ ID, RS ~ CS, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  dsixc = c(RS ~ HWI, RS ~ BM, RS ~ BILL, RS ~ ID, RS ~ CS, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  dsevena = c(RS ~ HWI, RS ~ BM, RS ~ EM, RS ~ ID, RS ~ CS, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  dsevenb = c(RS ~ HWI, RS ~ BILL, RS ~ EM, RS ~ ID, RS ~ CS, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  dsevenc = c(RS ~ HWI, RS ~ BM, RS ~ BILL, RS ~ EM, RS ~ ID, RS ~ CS, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  eonea = c(RS ~ BM, RS ~ EM, HWI ~ BM, HWI ~ BILL),
                  eoneb = c(RS ~ BILL, RS ~ EM, HWI ~ BM, HWI ~ BILL),
                  eonec = c(RS ~ BM, RS ~ BILL, RS ~ EM, HWI ~ BM, HWI ~ BILL),
                  etwoa = c(RS ~ BM, RS ~ ID, HWI ~ BM, HWI ~ BILL),
                  etwob = c(RS ~ BILL, RS ~ ID, HWI ~ BM, HWI ~ BILL),
                  etwoc = c(RS ~ BM, RS ~ BILL, RS ~ ID, HWI ~ BM, HWI ~ BILL),
                  ethreea = c(RS ~ BM, RS ~ CS, HWI ~ BM, HWI ~ BILL),
                  ethreeb = c(RS ~ BILL, RS ~ CS, HWI ~ BM, HWI ~ BILL),
                  ethreec = c(RS ~ BM, RS ~ BILL, RS ~ CS, HWI ~ BM, HWI ~ BILL),
                  efoura = c(RS ~ BM, RS ~ EM, RS ~ ID, HWI ~ BM, HWI ~ BILL),
                  efourb = c(RS ~ BILL, RS ~ EM, RS ~ ID, HWI ~ BM, HWI ~ BILL),
                  efourc = c(RS ~ BM, RS ~ BILL, RS ~ EM, RS ~ ID, HWI ~ BM, HWI ~ BILL),
                  efivea = c(RS ~ BM, RS ~ EM, RS ~ CS, HWI ~ BM, HWI ~ BILL),
                  efiveb = c(RS ~ BILL, RS ~ EM, RS ~ CS, HWI ~ BM, HWI ~ BILL),
                  efivec = c(RS ~ BM, RS ~ BILL, RS ~ EM, RS ~ CS, HWI ~ BM, HWI ~ BILL),
                  esixa = c(RS ~ BM, RS ~ ID, RS ~ CS, HWI ~ BM, HWI ~ BILL),
                  esixb = c(RS ~ BILL, RS ~ ID, RS ~ CS, HWI ~ BM, HWI ~ BILL),
                  esixc = c(RS ~ BM, RS ~ BILL, RS ~ ID, RS ~ CS, HWI ~ BM, HWI ~ BILL),
                  esevena = c(RS ~ BM, RS ~ EM, RS ~ ID, RS ~ CS, HWI ~ BM, HWI ~ BILL),
                  esevenb = c(RS ~ BILL, RS ~ EM, RS ~ ID, RS ~ CS, HWI ~ BM, HWI ~ BILL),
                  esevenc = c(RS ~ BM, RS ~ BILL, RS ~ EM, RS ~ ID, RS ~ CS, HWI ~ BM, HWI ~ BILL),
                  fonea = c(RS ~ BM, RS ~ EM, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  foneb = c(RS ~ BILL, RS ~ EM, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  fonec = c(RS ~ BM, RS ~ BILL, RS ~ EM, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  ftwoa = c(RS ~ BM, RS ~ ID, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  ftwob = c(RS ~ BILL, RS ~ ID, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  ftwoc = c(RS ~ BM, RS ~ BILL, RS ~ ID, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  fthreea = c(RS ~ BM, RS ~ CS, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  fthreeb = c(RS ~ BILL, RS ~ CS, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  fthreec = c(RS ~ BM, RS ~ BILL, RS ~ CS, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  ffoura = c(RS ~ BM, RS ~ EM, RS ~ ID, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  ffourb = c(RS ~ BILL, RS ~ EM, RS ~ ID, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  ffourc = c(RS ~ BM, RS ~ BILL, RS ~ EM, RS ~ ID, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  ffivea = c(RS ~ BM, RS ~ EM, RS ~ CS, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  ffiveb = c(RS ~ BILL, RS ~ EM, RS ~ CS, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  ffivec = c(RS ~ BM, RS ~ BILL, RS ~ EM, RS ~ CS, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  fsixa = c(RS ~ BM, RS ~ ID, RS ~ CS, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  fsixb = c(RS ~ BILL, RS ~ ID, RS ~ CS, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  fsixc = c(RS ~ BM, RS ~ BILL, RS ~ ID, RS ~ CS, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  fsevena = c(RS ~ BM, RS ~ EM, RS ~ ID, RS ~ CS, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  fsevenb = c(RS ~ BILL, RS ~ EM, RS ~ ID, RS ~ CS, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  fsevenc = c(RS ~ BM, RS ~ BILL, RS ~ EM, RS ~ ID, RS ~ CS, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  gone = c(RS ~ HWI, RS ~ EM, EM ~ BM),
                  gtwo = c(RS ~ HWI, RS ~ ID, EM ~ BM),
                  gthree = c(RS ~ HWI, RS ~ CS, EM ~ BM),
                  gfour = c(RS ~ HWI, RS ~ EM, RS ~ ID, EM ~ BM),
                  gfive = c(RS ~ HWI, RS ~ EM, RS ~ CS, EM ~ BM),
                  gsix = c(RS ~ HWI, RS ~ ID, RS ~ CS, EM ~ BM),
                  gseven = c(RS ~ HWI, RS ~ EM, RS ~ ID, RS ~ CS, EM ~ BM),
                  hone = c(RS ~ HWI, RS ~ EM, HWI ~ BM, HWI ~ BILL),
                  htwo = c(RS ~ HWI, RS ~ ID, HWI ~ BM, HWI ~ BILL),
                  hthree = c(RS ~ HWI, RS ~ CS, HWI ~ BM, HWI ~ BILL),
                  hfour = c(RS ~ HWI, RS ~ EM, RS ~ ID, HWI ~ BM, HWI ~ BILL),
                  hfive = c(RS ~ HWI, RS ~ EM, RS ~ CS, HWI ~ BM, HWI ~ BILL),
                  hsix = c(RS ~ HWI, RS ~ ID, RS ~ CS, HWI ~ BM, HWI ~ BILL),
                  hseven = c(RS ~ HWI, RS ~ EM, RS ~ ID, RS ~ CS, HWI ~ BM, HWI ~ BILL),
                  ione = c(RS ~ HWI, RS ~ EM, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  itwo = c(RS ~ HWI, RS ~ ID, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  ithree = c(RS ~ HWI, RS ~ CS, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  ifour = c(RS ~ HWI, RS ~ EM, RS ~ ID, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  ifive = c(RS ~ HWI, RS ~ EM, RS ~ CS, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  isix = c(RS ~ HWI, RS ~ ID, RS ~ CS, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  iseven = c(RS ~ HWI, RS ~ EM, RS ~ ID, RS ~ CS, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  jzeroa = c(RS ~ HWI, RS ~ BM, EM ~ BM),
                  jzerob = c(RS ~ HWI, RS ~ BILL, EM ~ BM),
                  jzeroc = c(RS ~ HWI, RS ~ BM, RS ~ BILL, EM ~ BM),
                  kzeroa = c(RS ~ HWI, RS ~ BM, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  kzerob = c(RS ~ HWI, RS ~ BILL, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  kzeroc = c(RS ~ HWI, RS ~ BM, RS ~ BILL, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  lone = c(RS ~ EM, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  ltwo = c(RS ~ ID, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  lthree = c(RS ~ CS, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  lfour = c(RS ~ EM, RS ~ ID, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  lfive = c(RS ~ EM, RS ~ CS, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  lsix = c(RS ~ ID, RS ~ CS, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  lseven = c(RS ~ EM, RS ~ ID, RS ~ CS, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  mzeroa = c(RS ~ BM, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  mzerob = c(RS ~ BILL, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  mzeroc = c(RS ~ BM, RS ~ BILL, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  nzero = c(RS ~ HWI, EM ~ BM, HWI ~ BM, HWI ~ BILL),
                  .common = c(ID ~ EM, CS ~ EM, ID ~ CS, BILL ~ BM)
                )
                
                phy<-read.newick("tree.nwk")
                setwd("PPA/")
                
                life_traits<-read.csv(file=paste(label,"/Warming_Positive.PPA.txt",sep=""), sep="",dec=".",header=TRUE)
                life_traits$SP<-as.character(life_traits$species)
                rownames(life_traits)<-life_traits[,1]
                list2WP<-data.frame(RS=life_traits$Response,BILL=life_traits$Bill_length,CS=life_traits$Clutch_size,EM=life_traits$Egg_mass,ID=life_traits$Inc_duration,HWI=life_traits$HWI,BM=life_traits$Body_mass,row.names=life_traits$SP)
                list2WP$RS <- ifelse(life_traits$RS==1, "CASE", "OTHERS")
                
                life_traits<-read.csv(file=paste(label,"/Warimg_Negative.PPA.txt",sep=""), sep="",dec=".",header=TRUE)
                life_traits$SP<-as.character(life_traits$species)
                rownames(life_traits)<-life_traits[,1]
                list2WN<-data.frame(RS=life_traits$Response,BILL=life_traits$Bill_length,CS=life_traits$Clutch_size,EM=life_traits$Egg_mass,ID=life_traits$Inc_duration,HWI=life_traits$HWI,BM=life_traits$Body_mass,row.names=life_traits$SP)
                list2WN$RS <- ifelse(life_traits$RS==2, "CASE", "OTHERS")
                
                life_traits<-read.csv(file=paste(label,"/Climate_Sensitive.PPA.txt",sep=""), sep="",dec=".",header=TRUE)
                life_traits$SP<-as.character(life_traits$species)
                rownames(life_traits)<-life_traits[,1]
                list2CS<-data.frame(RS=life_traits$Response,BILL=life_traits$Bill_length,CS=life_traits$Clutch_size,EM=life_traits$Egg_mass,ID=life_traits$Inc_duration,HWI=life_traits$HWI,BM=life_traits$Body_mass,row.names=life_traits$SP)
                list2CS$RS <- ifelse(life_traits$RS==2, "CASE_ONE", "CASE_TWO")
                
                life_traits<-read.csv(file=paste(label,"/Consistent_Decrease.PPA.txt",sep=""), sep="",dec=".",header=TRUE)
                life_traits$SP<-as.character(life_traits$species)
                rownames(life_traits)<-life_traits[,1]
                list2CD<-data.frame(RS=life_traits$Response,BILL=life_traits$Bill_length,CS=life_traits$Clutch_size,EM=life_traits$Egg_mass,ID=life_traits$Inc_duration,HWI=life_traits$HWI,BM=life_traits$Body_mass,row.names=life_traits$SP)
                list2CD$RS <- ifelse(life_traits$RS==1, "CASE", "OTHERS")
                
                life_traits<-read.csv(file=paste(label,"/Climate_Warming.PPA.txt",sep=""), sep="",dec=".",header=TRUE)
                life_traits$SP<-as.character(life_traits$species)
                rownames(life_traits)<-life_traits[,1]
                list2Warm<-data.frame(RS=life_traits$Response,BILL=life_traits$Bill_length,CS=life_traits$Clutch_size,EM=life_traits$Egg_mass,ID=life_traits$Inc_duration,HWI=life_traits$HWI,BM=life_traits$Body_mass,row.names=life_traits$SP)
                list2Warm$RS <- ifelse(life_traits$RS==1, "CASE_ONE", "CASE_TWO")
                
                #Fig4A#
                result<-phylo_path(models,list2WP,phy,btol=45,lower.bound=0,log.alpha.bound=6)
                output <- summary(result)
                output
                best_model <- best(result)
                best_model
                average_model <- average(result)
                average_model
                
                #Fig4B#
                result<-phylo_path(models,list2WN,phy,btol=45,lower.bound=0)
                output <- summary(result)
                output
                best_model <- best(result)
                best_model
                average_model <- average(result)
                average_model
                
                #Fig4C#
                result<-phylo_path(models,list2CS,phy,btol=45,lower.bound=0)
                output <- summary(result)
                output
                best_model <- best(result)
                best_model
                average_model <- average(result)
                average_model
                
                #Fig4D#
                result<-phylo_path(models,list2CD,phy,btol=45,lower.bound=0)
                output <- summary(result)
                output
                best_model <- best(result)
                best_model
                average_model <- average(result)
                average_model
                
                #SFig12E#
                result<-phylo_path(models,list2Warm,phy,btol=45,lower.bound=0,log.alpha.bound=5)
                output <- summary(result)
                output
                best_model <- best(result)
                best_model
                average_model <- average(result)
                average_model







#####################################################3
#####################################################
### End







#####################################################3
#####################################################
### contemporary data 

contemporary_data<-read.csv("Contemporary data for mix model.csv", header = T, stringsAsFactors = F)
  
contemporary_data$Latitude_num<-as.numeric(contemporary_data$Latitude)  
contemporary_data$Abs_Lat_num<-as.numeric(contemporary_data$abs_lat)   
contemporary_data$clutch_size_num<-as.numeric(contemporary_data$clutch_size)  
contemporary_data$clutch_size_log_num<-as.numeric(contemporary_data$clutch_size_log)   
contemporary_data$inc_bodymass_num<-as.numeric(contemporary_data$inc_bodymass)  
contemporary_data$inc_bodymass_num_log<-as.numeric(contemporary_data$inc_bodymass_log)   
contemporary_data$body_mass.log_num<-as.numeric(contemporary_data$body_mass.log.)   
  

contemporary_data_Complete_Cases<-completeFun(contemporary_data, "clutch_size_num")
contemporary_data_Complete_Cases<-completeFun(contemporary_data_Complete_Cases, "Egg.mass..g.")
contemporary_data_Complete_Cases<-completeFun(contemporary_data_Complete_Cases, "incubation_d")
contemporary_data_Complete_Cases<-completeFun(contemporary_data_Complete_Cases, "body_mass")
contemporary_data_Complete_Cases<-completeFun(contemporary_data_Complete_Cases, "Beak_Length..culmen.")
contemporary_data_Complete_Cases<-completeFun(contemporary_data_Complete_Cases, "HWI")



Abs_Latitude_MODEL<-lm(scale(Abs_Lat_num) ~ 
                     + scale(clutch_size_num) + scale(Egg.mass..g.) 
                   + scale(incubation_d) +  scale(body_mass.log_num)
                   + scale(Beak_Length..culmen.) + scale(HWI), 
                   na.action = na.fail, 
                   data = contemporary_data_Complete_Cases)


summary(Abs_Latitude_MODEL)     
r.squaredGLMM(Abs_Latitude_MODEL) 
confint(Abs_Latitude_MODEL)



allmodels_Abs_Latitude_MODEL<-dredge(Abs_Latitude_MODEL)


sub_allmodels_Abs_Latitude_MODEL<-subset(allmodels_Abs_Latitude_MODEL, delta<5) 



#tiff('Effects_output_Abs_Latitude_MODEL.tiff',width=6400,height=6000,units = 'px', res = 800)

## set CI lines to black
t1 = trellis.par.get("plot.line")
t1$col <- "black"
trellis.par.set("plot.line",t1)

## turn off grid lines
d1 <- trellis.par.get("dot.line")
d1$lwd <- 0  ## hack -- set line width to 0
d1$lty <- "blank"  ## better hack, set line type to blank
trellis.par.set("dot.line",d1)

Abs_Latitude_PLOT<-Dotplot(Effects_output_Abs_Latitude_MODEL$Component ~ Cbind(Effects_output_Abs_Latitude_MODEL$Estimate,Effects_output_Abs_Latitude_MODEL$LCI,Effects_output_Abs_Latitude_MODEL$UCI),
                           xlim=c(-0.4,0.4), xlab="",ylab="",col="black",
                           panel=function(x,y){
                             panel.Dotplot(x, y, col="black",cex = 1.3, pch = 19,  # adjust point size and shape here
                                           panel.abline(v=0.0, lty=2))},  ## to add the verticle line
                           strip=strip.custom(style=1, bg="lightgrey"))


## adjust location of tick marks manually, then add in the y axis labels, 
xpos <- seq(from=(-0.4), by=0.2, to=0.4)
tlabs<-c("Incubation duration", "HWI", "Egg mass", "Clutch size","Body mass","Bill length") ## has to be in opposite alphabetical order

update(Abs_Latitude_PLOT, auto.key=list(columns=2), xlab=list("Estimate",cex = 1.3), ylab=list("Trait",cex = 1.3),
       scales=list(tck= c(1,0), x=list(at=xpos), y=list(labels=tlabs), cex =1.3))  ## tck sets the tick marks to only the left and botto   

#dev.off()









#######################################################
##########  FIN  
######################################################




