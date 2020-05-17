library(wesanderson)
library(ggplot2)
library(tidyverse)
library(ggthemes)
library(reshape2) #change the data format from long to wide and vice versa

#Let's make a dataset here
A = c(528, 621, 10, 0, 307, 7091, 0)
B = c(548, 1735, 23, 0, 1807, 2401, 0)
C = c(2241, 3076, 358, 47, 2948, 3593, 131)
D = c(2944, 3963, 898, 792, 1740, 2345, 1713)
E = c(3094, 6221, 2242, 2918, 2667, 2063, 2187)
O = c(84, 199, 0, 0, 3, 33525, 0)
Assembler = c("Canu", "Flye", "Miniasm", "Raven", "Redbean", "Shasta", "SmartDN")

#Let's combine these data
df3 = data.frame(O, A, B, C, D, E, Assembler)

#Now let's change the data frame to long format
df3 <- melt(df3, id.vars = 'Assembler')

#We can change the variable names if needed, may be there are more efficient methods but this is my own methode but will work for small data.
df4 <- within(df3, levels(variable)[levels(variable) == "O"] <- "0-1000bp")
df3 <- within(df4, levels(variable)[levels(variable) == "A"] <- "1000-4999bp")
df4 <- within(df3, levels(variable)[levels(variable) == "B"] <- "5000-9999bp")
df3 <- within(df4, levels(variable)[levels(variable) == "C"] <- "10000-24999bp")
df4 <- within(df3, levels(variable)[levels(variable) == "D"] <- "25000-49999bp")
df3 <- within(df4, levels(variable)[levels(variable) == "E"] <- "50000+bp")

#Check the dataframe now
df3

#Let's plot a line graph now (see it for yourself, you will understand what each line is doing)
ggplot(data = df3, mapping = aes(x=variable, y=value, color = Assembler, group=1)) +
  geom_line() +
  facet_grid(rows = vars(Assembler)) +
  labs(title = "Contig's length distribution",
       x = 'bp',
       y= "No. of Contigs") +
  theme_bw() +
  theme(text = element_text(size = 16))+
  theme(plot.title = element_text(hjust = 0.5))

# we need ggforce to use zoom function I think
install.packages("ggforce")
library(ggforce)

#Let's plot the same data in bargraph with zoom function (from 0-11000 y coordinates)
ggplot(data = df3, aes(fill= variable, y=value, x=Assembler))+
  geom_bar(position = "stack", stat = "identity")+
  theme_bw(base_size = 15)+
  labs(y = "No. of Contigs") +
  ggtitle("Contigs distribution") +
  theme(plot.title = element_text(hjust = 0.5))+
  facet_zoom(ylim = c(0, 11000))


#Lets make another data here
Greater_0bp = c(9439, 15815, 3531, 3757, 9472, 51018, 4031)
Greater_999bp = c(9355, 15616, 3531, 3757, 9469, 17493, 4031)
Greater_4999bp = c(8827, 14995, 3521, 3757, 9162, 10402, 4031)
Greater_9999bp = c(8279, 13260, 3498, 3757, 7355, 8001, 4031)
Greater_24999bp = c(6038, 10184, 3140, 3710, 4407, 4408, 3900)
Greater_49999bp = c(3094, 6221, 2242, 2918, 2667, 2063, 2187)

#Let's combine these data together with Assembler from last data set
df = data.frame(Greater_0bp, Greater_999bp, Greater_4999bp, Greater_9999bp, Greater_24999bp, Greater_49999bp, Assembler)
#Check out how the data looks like
df

#when we melt data, we basically create only one column with the values
#and one column with the variable i.e. your x,y,z
df <- melt(df, id.vars = 'Assembler')

#check how data looks like
df

#Now let's plot another bar graph
ggplot(df, aes(x=Assembler, y=value, fill=variable)) +
  geom_bar(stat='identity', position='dodge') +
  labs(y = "No._of_Contigs") +
  ggtitle("Contigs distribution") +
  theme_bw(base_size = 20) +
  theme(plot.title = element_text(hjust = 0.5))

#Another data here (Total length of different assemblers)
Total_length = c(442.148406, 787.657953, 256.365835, 330.277576, 455.550962, 372.483022, 259.363014)
#Combine the data
df1 = data.frame(Total_length, Assembler)
#melt the data to long format
df1 <- melt(df1, id.vars = 'Assembler')

#check data
df1

#plot a simple bargraph
ggplot(df1, aes(x=Assembler, y=value, fill=Assembler)) +  #df1 is the data,, we ussed assemblers to fill
  geom_bar(stat='identity', width = 0.5, position='dodge') +
  labs(y = "Mb") +
  ggtitle("Total length") +
  theme_bw(base_size = 20) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))

#I just learned how to setup working direcotory, YaY
setwd("erg")
getwd()
list.files()

#lets make another data set and  melt as before #For busco
Complete_BUSCO = c(31.68, 50.5, 0, 19.47, 38.61, 18.81, 9.57)
Partial_BUSCO = c(23.43, 18.48, 0.33, 15.51, 18.81, 9.24, 5.94)
dfn = data.frame(Complete_BUSCO, Partial_BUSCO, Assembler)
dfn <- melt(dfn, id.vars = 'Assembler')
dfn

#LEts plot a bar graph
ggplot(dfn, aes(x=Assembler, y=value, fill=variable)) +
  geom_bar(stat='identity', position='dodge') +
  labs(y = "%") +
  ylim(0, 75) +
  ggtitle("RAW assembly BUSCO score") +
  theme(legend.title = element_blank())+
  theme_bw(base_size = 20) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label=round(value, 2), fontface="bold", vjust=2), position=position_dodge(width=0.9), size=3.5)

#YAY, Now I learned how to import CSV file from the computer
df2 <- read.csv("buscodata.csv", header = TRUE, sep = ",")

#change the variable names
dfx <- within(df2, levels(variable_1)[levels(variable_1) == "No_trimming"] <- "No trim")
df2 <- within(dfx, levels(variable_1)[levels(variable_1) == "After_trimming"] <- "Trimmed")
df2

#Lets plot a bar graph
ggplot(df2, aes(x=variable_1, y=value, fill=factor(variable_2, c("Complete_BUSCO", "Partial_BUSCO")))) +
  geom_bar(position=position_dodge(width=0.9), stat="identity") +
  facet_grid(~Assemblers) +
  geom_text(aes(label=round(value, 2), fontface="bold", vjust=2), position=position_dodge(width=0.9), size=3.5) +
  theme_bw(base_size = 20) +
  labs(y = "%")+
  theme(legend.title=element_blank(), legend.position="bottom",
        axis.text.y=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_text(face="bold"), axis.ticks.x=element_blank(), 	axis.ticks.y=element_blank(),
        panel.background=element_blank(), panel.grid.minor=element_blank(), panel.grid.major=element_blank())

#another data set imported
finalstat <- read.csv("finalstat.csv",header = TRUE, sep = ",")
finalstat

#Plotting another bargraph
ggplot(finalstat, aes(x=Assemblers, y=X., fill=BUSCO)) +
  geom_bar(stat='identity', position='dodge') +
  labs(y = "%") +
  ggtitle("BUSCO score") +
  theme_bw(base_size = 20) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label=round(X., 2), fontface="bold", vjust=2), position=position_dodge(width=0.9), size=3.5) +
  geom_text(aes(y=5, label=Polishing, fontface="bold"), color="white", position=position_dodge(width=0.9), angle=90, size=4)
