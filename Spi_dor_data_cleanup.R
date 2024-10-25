####We have to separate data clean up for each species to properly include
# zero capture dates. 

#This is for Spi.dors
pkgs <-
  c("tibble",
    "readr",
    "here",
    "tidyverse",
    "mgcv",
    "janitor",
    "writexl",
    "readxl")

vapply(
  pkgs,
  library,
  FUN.VALUE = logical(1L),
  character.only = TRUE,
  logical.return = TRUE
)


df_muscidae <- readxl::read_xlsx("Data/Zackenberg_Muscidae_SL.xlsx", sheet = "Sheet1")

unique(df_muscidae$Species_name)

#We need to include zero capture dates. We can get this information from the 
# full dataset in the GEM database
df7 <- read.csv("Data/EMdata_final_AbundancePTD.csv",sep=",",stringsAsFactors = FALSE, header = TRUE)
#df8<-subset(df7,select=-c(X))


#df7$Year<- as.factor(df7$Year)#den var en integer, men det er nemmere at arbejde med den som en faktor.

class(df7$Abundance)
df7$Abundance<-as.numeric(df7$Abundance)

#Subset Muscidae and relevant plots

df_musc_all <- subset(df7, SpeciesID == "ANMU")

df_musc_all <- subset(df_musc_all, Plot == "Art2" | Plot == "Art3" | Plot == "Art5")

df_musc_all <- subset(df_musc_all, Year < 2015)

#Reorder rows based on year and DOY
df_muscidae %>%
  arrange(Year, DOY) -> df_muscidae

#Rename rows in plot column
df_muscidae$Plot[df_muscidae$Plot == "2"] <- "Art2" 
df_muscidae$Plot[df_muscidae$Plot == "3"] <- "Art3" 
df_muscidae$Plot[df_muscidae$Plot == "5"] <- "Art5" 



#Calculating abundance in each trap
#Note: In 1996 traps are called GUL. Indicates that cups have been pooled?
df_muscidae%>%
  group_by(Year,Plot,Month,DOY,Species_name)%>%
  summarise(A=sum(Trap == "A"),B=sum(Trap == "B"), GUL=sum(Trap == "GUL"),GUL_unknown=sum(Trap == "GUL?")) -> df_new

df_new%>%
  group_by(Year,Plot,Month,DOY,Species_name)%>%
  summarise(Abundance=sum(A+B+GUL+GUL_unknown)) -> df1a

#Because all four cups were pooled in 1996, the abundance is divided by two
df1a$Abundancecorr <- ifelse(df1a$Year<1997, df1a$Abundance/2, df1a$Abundance)
class(df1a$Abundancecorr)

df1 <- df1a %>% 
  select(Year, Plot, DOY, Species_name, Abundancecorr)

df1 = df1 %>% rename("Abundance"="Abundancecorr")

#We need to merge/bind rows based on Plot, Year and DOY

#Subset for each species before merging????

df_Spidor <- subset(df1, Species_name == "Spi.dor")


df1b <- merge(df_musc_all, df_Spidor, by = c("Plot", "Year", "DOY"), all = T)

df2 <- df1b %>% 
  select(Year, Plot, DOY, Species_name, Abundance.y)

df2 = df2 %>% rename("Abundance"="Abundance.y")


df2$Abundance[is.na(df2$Abundance)] <- 0 #Replace 'na' with zero
df2$Species_name[is.na(df2$Species_name)] <- "Spi.dor"

#Implementing criteria for abundance and events
#Abundance = At least 3 individuals
#Events = At least 3 cases with abundance criteria fulfilled

df2$Event<-ifelse(df2$Abundance>0,1,0)

df2%>%
  group_by(Plot,Year)%>%
  summarise(TotalAbundance=sum(Abundance),TotalEvents=sum(Event))->df2a

df2a$TotalAbunAndEventCriteria<-ifelse(df2a$TotalAbundance>3&df2a$TotalEvents>2,1,0)

df2a%>%
  group_by(Plot)%>%
  summarise(TotalYear=sum(TotalAbunAndEventCriteria))->df2c

df2$TotalYear <- (df2c$TotalYear[match(paste0(df2$Plot),paste0(df2c$Plot))])

df2a$Include<-ifelse(df2a$TotalAbundance>3&df2a$TotalEvents>2,1,0)#Need at least 50 individuals in a season and 3 capture events

#Filter original data for sampling criterias
df2$Include <- (df2a$Include[match(paste0(df2$Year,df2$Plot),paste0(df2a$Year,df2a$Plot))])

write.csv(df2, file = "Data/Species_data_clean\\Spi_dor.csv", row.names=FALSE)

