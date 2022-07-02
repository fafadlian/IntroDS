#Libraries
library(tidyverse)
library(dplyr)
library(ggplot2)
library(WDI)
library(corrplot)
library(GGally)
library(factoextra)
library(cluster)
library(fpc)
library(ggpubr)
library(MASS)
library(lubridate)
library(ggstatsplot)

#WDI Loading
new_wdi_cache<-WDIcache()

WDI<-WDI(country = "all",
         indicator = c("NY.GDP.MKTP.KD.ZG",
                       "NV.IND.MANF.ZS", 
                       "NV.SRV.TOTL.ZS", 
                       "NV.IND.TOTL.ZS", 
                       "NV.AGR.TOTL.ZS"),
         start = 2020,
         end = 2020,
         extra = TRUE,
         cache = new_wdi_cache)

#COVID19 Data Loading
covidData<-read.csv("covidData.csv")

covidDF<-data.frame(covidData[,(colnames(covidData) %in% c('iso_code',
                                                            'continent',
                                                            'location',
                                                            'date',
                                                            'total_cases_per_million',
                                                            'total_deaths_per_million',
                                                            'total_vaccinated_per_hundred',
                                                            'gdp_per_capita',
                                                            'life_expectancy',
                                                            'human_development_index',
                                                            'population_density',
                                                            'hospital_beds_per_thousand',
                                                            'total_vaccinations_per_hundred'))])

#MFeature Engineering on Stringency Index Value
covidDF_stringency<-data.frame(covidData[, (colnames(covidData) %in% c('iso_code',
                                                                       'location',
                                                                       'date',
                                                                       "stringency_index"))])

stringency_2020<-covidDF_stringency %>% filter(date <= as.Date("2020-12-31"))

#Value for maximum and mean stringency index
stringency_mean<-stringency_2020 %>%
  group_by(location) %>%
  summarise(stringency_index = mean(stringency_index, na.rm = TRUE))

colnames(stringency_mean)<-c("location", "mean_stringency_idx")

stringency_max<-stringency_2020 %>%
  group_by(location) %>%
  summarise(stringency_index = max(stringency_index, na.rm = TRUE))

colnames(stringency_max)<-c("location", "max_stringency_idx")


#Creating new column names for WDI Data
WDI_DF<-data.frame(WDI$country,
                   WDI$NY.GDP.MKTP.KD.ZG,
                   WDI$NV.SRV.TOTL.ZS,
                   WDI$NV.IND.TOTL.ZS,
                   WDI$NV.AGR.TOTL.ZS)

colnames(WDI_DF)<-c("location",
                    "GDP_growth_rate_2020",
                    "service_sect",
                    "industry_sect",
                    "agriculture_sect")



#Date Filtering and Datasets Merging
covidDF<-filter(covidDF, date=="2020-12-31")
covidDF<-merge(covidDF, WDI_DF,by = "location")
covidDF<-merge(covidDF, stringency_mean, by = "location")
covidDF<-merge(covidDF, stringency_max, by = "location")


#Dealing with data quality
#NA Value Replacement
covidDF<-mutate_at(covidDF, c("total_vaccinations_per_hundred"), ~replace(., is.na(.),0.0))

#NA Value Deletion
covidDF<-na.omit(covidDF)


#Outliers Deletion with IQR Method
Q<-quantile(covidDF$GDP_growth_rate_2020, probs = c(.25, .75), na.rm = FALSE)
iqr<-IQR(covidDF$GDP_growth_rate_2020)
up<-Q[2]+1.5*iqr #upper interquartile 
low<-Q[2]-1.5*iqr #lower interquartile

covidDF_no_outlier<-subset(covidDF, 
                           covidDF$GDP_growth_rate_2020>low &
                             covidDF$GDP_growth_rate_2020<up) #Deleting datapoints outside interquartile range


#Comparing dataset before and after cleaning
covidDF_no_outlier_check<-covidDF_no_outlier
covidDF_check<-covidDF

covidDF_no_outlier_check$remark<-"Cleaned"
covidDF_check$remark<-"Original"

covidDF_check_ori_and_no_outliers<-rbind(covidDF_no_outlier_check, covidDF_check)

ggplot(covidDF_check_ori_and_no_outliers, aes(remark, GDP_growth_rate_2020))+
  geom_boxplot()+
  geom_jitter()+
  ggtitle("Original and Cleaned Data Comparisson")+
  xlab("Remark")+
  ylab("GDP Growth Rate 2020")+
  theme(plot.title = element_text(hjust = 0.5))

#Correlation Matric for All Variables
cor_data1<-select_if(covidDF_no_outlier, is.numeric)
cor = cor(cor_data1)
corrplot(cor, method = 'number')

#Data Scaling
covidDF_clean<-na.omit(covidDF_no_outlier)
covidDF_sclaed<-scale(covidDF_clean[,c(4:17)])

#Elbow Analysis
fviz_nbclust(
  covidDF_sclaed,
  FUNcluster = kmeans,
  method = "wss",
  diss = NULL,
  k.max = 10,
  nboot = 100,
  verbose = interactive(),
  barfill = "steelblue",
  barcolor = "steelblue",
  linecolor = "steelblue",
  print.summary = TRUE)

#KMEANS Algorithm with 3 clusters
pamk_result<-pamk(covidDF_sclaed,krange=3,criterion="asw", usepam=TRUE,
                  scaling=FALSE, alpha=0.001, diss=inherits(data, "dist"),
                  critout=FALSE, ns=10, seed=NULL)

clusters<-data.frame(pamk_result$pamobject$clustering)

covidDF_clean<-cbind(covidDF_clean, clusters)
colnames(covidDF_clean)[18]<-"clusters"

covidDF_clean_scaled<-data.frame(cbind(covidDF_clean[,1:3], covidDF_sclaed, clusters))

#correlation matrix check before and after scaling (result : same correlation value between variables)
cor_data_clean<-select_if(covidDF_clean, is.numeric)
cor = cor(cor_data_clean)
corrplot(cor, method = 'number')

cor_data_clean_scaled<-select_if(covidDF_clean_scaled, is.numeric)
cor = cor(cor_data_clean_scaled)
corrplot(cor, method = 'number')




#Principal COmponent Analysis
covidDF_clean_scaled_for_PCA<-subset(covidDF_clean_scaled, select=-c(GDP_growth_rate_2020))
pca<-prcomp(covidDF_clean_scaled_for_PCA[,4:16])
covidPCA<-data.frame(covidDF_clean_scaled,
                     PC1=pca$x[,1],
                     PC2=pca$x[,2])

colnames(covidPCA)[18]<-"clusters"

#Visualise Datapoints in the PCA Axes
ggplot(complete_data,aes(PC1, PC2))+
  geom_point(aes(col=as.character(clusters), shape = income), size = 3)+
  ggtitle("Visualisation with PCA Dimension")+
  xlab("PC1")+
  ylab("PC2")+
  labs(col = "Clusters", shape = "Income Level")

WDI_select<-data.frame(WDI$country, WDI$income)
colnames(WDI_select)<-c("location", "income")

complete_data<-merge(covidPCA, WDI_select, by = )
complete_data<-subset(complete_data, location != "World")

complete_data_clean<-merge(covidDF_clean, WDI_select, by = )
complete_data_clean<-subset(complete_data_clean, location != "World")

cor_data_clean_scaled_PCA<-select_if(complete_data, is.numeric)
cor_data_clean_scaled_PCA<-subset(cor_data_clean_scaled_PCA, select=-c(clusters))
cor = cor(cor_data_clean_scaled_PCA)
corrplot(cor, method = 'number')




#Scatter Plot Service Sector vs GDP Growth Rate
ggplot(complete_data_clean,aes(GDP_growth_rate_2020, service_sect))+
  geom_point(aes(col=as.character(clusters), shape = income), size = 3)+
  geom_smooth()+
  ggtitle("Relationship of GDP Growth Rate 2020 and Service Sector")+
  xlab("GDP Growth Rate 2020")+
  ylab("Service Sector (% of GDP)")+
  labs(col = "Clusters", shape = "Income")

#Scatter Plot Total Deaths vs GDP Growth Rate
ggplot(complete_data_clean,aes(GDP_growth_rate_2020, total_deaths_per_million))+
  geom_point(aes(col=as.character(clusters), shape = income), size = 3)+
  geom_smooth(method = lm)+
  ggtitle("Relationship of GDP Growth Rate 2020 and Total Deaths per Millionr")+
  xlab("GDP Growth Rate 2020")+
  ylab("Total Deaths Per Million")+
  labs(col = "Clusters", shape = "Income")

#Scatter Plot Total Cases vs GDP Growth Rate
ggplot(complete_data_clean,aes(GDP_growth_rate_2020, total_cases_per_million))+
  geom_point(aes(col=as.character(clusters), shape = income), size = 3)+
  geom_smooth(method = lm)+
  ggtitle("Relationship of GDP Growth Rate 2020 and Total Cases per Millionr")+
  xlab("GDP Growth Rate 2020")+
  ylab("Total Cases Per Million")+
  labs(col = "Clusters", shape = "Income")

#Scatter Plot Agriculture Sector vs GDP Growth Rate
ggplot(complete_data_clean,aes(GDP_growth_rate_2020, agriculture_sect))+
  geom_point(aes(col=as.character(clusters), shape = income), size = 3)+
  geom_smooth(method = lm)+
  ggtitle("Relationship of GDP Growth Rate 2020 and Agriculture Sector")+
  xlab("GDP Growth Rate 2020")+
  ylab("Agriculture Sector (% of GDP)")+
  labs(col = "Clusters", shape = "Income")

#Scatter Plot Human Development Index vs GDP Growth Rate
ggplot(complete_data_clean,aes(GDP_growth_rate_2020, human_development_index))+
  geom_point(aes(col=as.character(clusters), shape = income), size = 3)+
  geom_smooth(method = lm)+
  ggtitle("Relationship of GDP Growth Rate 2020 and Human Development Index")+
  xlab("GDP Growth Rate 2020")+
  ylab("Human Development Index")+
  labs(col = "Clusters", shape = "Income")

#Scatter Plot Life Expectancy vs GDP Growth Rate
ggplot(complete_data_clean,aes(GDP_growth_rate_2020, life_expectancy))+
  geom_point(aes(col=as.character(clusters), shape = income), size = 3)+
  geom_smooth(method = lm)+
  ggtitle("Relationship of GDP Growth Rate 2020 and Life Expectancy")+
  xlab("GDP Growth Rate 2020")+
  ylab("Life Expectancy")+
  labs(col = "Clusters", shape = "Income")


#Boxplot Clusters and GDP Growth Rate
ggplot(complete_data_clean, aes(as.character(clusters), GDP_growth_rate_2020))+
  geom_boxplot()+
  ggtitle("Cluster vs GDP Growth Rate 2020")+
  xlab("clusters")+ylab("GDP Growth Rate 2020 (%)")



#Create data frames for each clusters
cluster1<-filter(complete_data_clean, clusters == 1)
cluster2<-filter(complete_data_clean, clusters == 2)
cluster3<-filter(complete_data_clean, clusters == 3)

#correlation matrix for cluster 1
cor_data_m1<-select_if(cluster1, is.numeric)
cor = cor(cor_data_m1)
corrplot(cor, method = 'number')

#correlation matrix for cluster 2
cor_data_m2<-select_if(cluster2, is.numeric)
cor = cor(cor_data_m2)
corrplot(cor, method = 'number')

#correlation matrix for cluster 3
cor_data_m3<-select_if(cluster3, is.numeric)
cor = cor(cor_data_m3)
corrplot(cor, method = 'number')

#Linear Regression Model for Each Cluster
#Multiple Linear Regression for Attempt 3 (Cluster 1)
model1<-lm(formula = GDP_growth_rate_2020~
             total_deaths_per_million+
             human_development_index+
             service_sect+
             max_stringency_idx, 
           data = cluster1)

summary(model1)

#Multiple Linear Regression for Attempt 3 (Cluster 2)
model2<-lm(formula = GDP_growth_rate_2020~
             total_deaths_per_million+
             human_development_index+
             service_sect+
             max_stringency_idx,  
           data = cluster2)

summary(model2)

#Multiple Linear Regression for Attempt 3 (Cluster 3)
model3<-lm(formula = GDP_growth_rate_2020~
             total_deaths_per_million+
             human_development_index+
             service_sect+
             max_stringency_idx, 
           data = cluster3)

summary(model3)



##Multiple Linear Regression for Attempt 5 (PC1 and PC2 Variables)
model_all_PCA<-lm(formula = GDP_growth_rate_2020~
                    PC1+
                    PC2, 
                  data = complete_data)

summary(model_all_PCA)

#Multiple Linear Regression for Attempt 4 (All clusters with original variables)
model_all_clean<-lm(formula = GDP_growth_rate_2020~total_deaths_per_million+
                      human_development_index+
                      service_sect+
                      max_stringency_idx, data = covidDF_clean)


summary(model_all_clean)



#Plotting Performance from each attempts


gdp_prediction_cluster1<-data.frame(location = cluster1$location, 
                                actual = cluster1$GDP_growth_rate_2020, 
                                pred = model1$fitted.values)

gdp_prediction_cluster2<-data.frame(location = cluster2$location, 
                                    actual = cluster2$GDP_growth_rate_2020, 
                                    pred = model2$fitted.values)

gdp_prediction_cluster3<-data.frame(location = cluster3$location, 
                                    actual = cluster3$GDP_growth_rate_2020, 
                                    pred = model3$fitted.values)

#Plotting for Attempt 1 (Cluster 1)
cluster1_plot<-ggplot(data = gdp_prediction_cluster1, aes(actual, pred))+
  geom_point()+
  geom_abline(mapping=aes(slope=1, intercept = 0), color='red')+
  geom_text(aes(label = location), vjust = 0, nudge_y=-0.1, check_overlap = TRUE)+
  ggtitle("Cluster 1 Actual vs Predicted GDP Growth Rate 2020")+
  xlab("Actual GDP Growth Rate 2020 (%)")+
  ylab("Predicted GDP Growth Rate 2020 (%)")
cluster1_plot

#Plotting for Attempt 2 (Cluster 2)
cluster2_plot<-ggplot(data = gdp_prediction_cluster2, aes(actual, pred))+
  geom_point()+geom_abline(mapping=aes(slope=1, intercept = 0), color='red')+
  geom_text(aes(label = location), vjust = 0, nudge_y=-0.2, check_overlap = TRUE)+
  ggtitle("Cluster 2 Actual vs Predicted GDP Growth Rate 2020")+
  xlab("Actual GDP Growth Rate 2020 (%)")+
  ylab("Predicted GDP Growth Rate 2020 (%)")
cluster2_plot

#Plotting for Attempt 3 (Cluster 3)
cluster3_plot<-ggplot(data = gdp_prediction_cluster3, aes(actual, pred))+
  geom_point()+geom_abline(mapping=aes(slope=1, intercept = 0), color='red')+
  geom_text(aes(label = location), vjust = 0, nudge_y=-0.2, check_overlap = TRUE)+
  ggtitle("Cluster 3 Actual vs Predicted GDP Growth Rate 2020")+
  xlab("Actual GDP Growth Rate 2020 (%)")+
  ylab("Predicted GDP Growth Rate 2020 (%)")
cluster3_plot

#Plotting for Attempt 5 (PCA Regression)
gdp_prediction_all_pca<-data.frame(location = complete_data$location,
                                   actual = complete_data$GDP_growth_rate_2020,
                                   pred = modelPCA$fitted.values)

gdp_prediction_all_pca<-merge(complete_data, gdp_prediction_all_pca, by = "location")

all_pca<-ggplot(data = gdp_prediction_all_pca, aes(actual, pred))+
  geom_point(aes(col=as.character(clusters), shape = income), size = 3)+
  geom_abline(mapping=aes(slope=1, intercept = 0), color='red')+
  ggtitle("Actual vs Predicted GDP Growth Rate 2020 (PCA)")+
  xlab("Actual GDP Growth Rate 2020 (%) - scaled")+
  ylab("Predicted GDP Growth Rate 2020 (%) - scaled")+
  geom_text(aes(label = location), vjust = 0, nudge_y=-0.05, check_overlap = TRUE)+
  labs(col = "Clusters", shape = "Income")

all_pca

#Plotting for Attempt 4 (Unclustered with Original Variable)
gdp_prediction_all_clean<-data.frame(location = complete_data_clean$location, 
                                     actual = complete_data_clean$GDP_growth_rate_2020, 
                                     pred = model_all_clean$fitted.values)

gdp_prediction_all_clean<-merge(gdp_prediction_all_clean, WDI_select,by = "location")
gdp_prediction_all_clean<-merge(gdp_prediction_all_clean, complete_data_clean, by = "location")

all_clean<-ggplot(data = gdp_prediction_all_clean, aes(actual, pred))+
  geom_point(aes(col=as.character(clusters), shape = income.y), size = 3)+
  geom_abline(mapping=aes(slope=1, intercept = 0), color='red')+
  ggtitle("Actual vs Predicted GDP Growth Rate 2020")+
  xlab("Actual GDP Growth Rate 2020 (%)")+
  ylab("Predicted GDP Growth Rate 2020 (%)")+
  geom_text(aes(label = location), vjust = 0, nudge_y=-0.3, check_overlap = TRUE)+
  labs(col = "Clusters", shape = "Income")
geom_text(aes(label = location))
all_clean



