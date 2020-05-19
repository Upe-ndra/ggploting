library(wesanderson)
library(ggplot2)
library(tidyverse)
library(ggthemes)
library(reshape2) #change the data format from long to wide and vice versa

# Let's start with making a dataset:
# This data is the grouped contig length distribution from quast
A = c(528, 621, 10, 0, 307, 7091, 0)
B = c(548, 1735, 23, 0, 1807, 2401, 0)
C = c(2241, 3076, 358, 47, 2948, 3593, 131)
D = c(2944, 3963, 898, 792, 1740, 2345, 1713)
E = c(3094, 6221, 2242, 2918, 2667, 2063, 2187)
O = c(84, 199, 0, 0, 3, 33525, 0)
Assembler = c("Canu", "Flye", "Miniasm", "Raven", "Redbean", "Shasta", "SmartDN")

# Let's combine these data
data1 = data.frame(O, A, B, C, D, E, Assembler)

# Now let's change the data frame to long format
data1 <- melt(data1, id.vars = 'Assembler')

# We can change the variable names if needed, may be there are more efficient methods but this is my own methode but will work for small data.
data2 <- within(data1, levels(variable)[levels(variable) == "O"] <- "0-1000bp")
data1 <- within(data2, levels(variable)[levels(variable) == "A"] <- "1000-4999bp")
data2 <- within(data1, levels(variable)[levels(variable) == "B"] <- "5000-9999bp")
data1 <- within(data2, levels(variable)[levels(variable) == "C"] <- "10000-24999bp")
data2 <- within(data1, levels(variable)[levels(variable) == "D"] <- "25000-49999bp")
data1 <- within(data2, levels(variable)[levels(variable) == "E"] <- "50000+bp")

# Check the dataframe now
data1

# Let's plot a line graph now (see it for yourself, you will understand what each line is doing)
ggplot(data = data1, mapping = aes(x=variable, y=value, color = Assembler, group=1)) +
  geom_line() +
  facet_grid(rows = vars(Assembler)) +
  labs(title = "Contig's length distribution",
       x = 'bp',
       y= "No. of Contigs") +
  theme_bw() +
  theme(text = element_text(size = 16))+
  theme(plot.title = element_text(hjust = 0.5))

# we need ggforce to use zoom function
install.packages("ggforce")
library(ggforce)

# Let's plot the same data in bargraph with zoom function (from 0-11000 y coordinates)
ggplot(data = data1, aes(fill= variable, y=value, x=Assembler))+
  geom_bar(position = "stack", stat = "identity")+
  theme_bw(base_size = 15)+
  labs(y = "No. of Contigs") +
  ggtitle("Contigs distribution") +
  theme(plot.title = element_text(hjust = 0.5))+
  facet_zoom(ylim = c(0, 11000))


# Let's make another data here
Greater_0bp = c(9439, 15815, 3531, 3757, 9472, 51018, 4031)
Greater_999bp = c(9355, 15616, 3531, 3757, 9469, 17493, 4031)
Greater_4999bp = c(8827, 14995, 3521, 3757, 9162, 10402, 4031)
Greater_9999bp = c(8279, 13260, 3498, 3757, 7355, 8001, 4031)
Greater_24999bp = c(6038, 10184, 3140, 3710, 4407, 4408, 3900)
Greater_49999bp = c(3094, 6221, 2242, 2918, 2667, 2063, 2187)

# Combine these data together with Assembler from last data set
data2 = data.frame(Greater_0bp, Greater_999bp, Greater_4999bp, Greater_9999bp, Greater_24999bp, Greater_49999bp, Assembler)
# Check out how the data looks like
data2

# when we melt data, we basically create only one column with the values and one column with the variable i.e. your x,y,z
data2 <- melt(data2, id.vars = 'Assembler')

# check how the data looks like
data2

# Now let's plot another bar graph
ggplot(data2, aes(x=Assembler, y=value, fill=variable)) +
  geom_bar(stat='identity', position='dodge') +
  labs(y = "No._of_Contigs") +
  ggtitle("Contigs distribution") +
  theme_bw(base_size = 20) +
  theme(plot.title = element_text(hjust = 0.5))

# Another data here (Total length of different assemblers)
Total_length = c(442.148406, 787.657953, 256.365835, 330.277576, 455.550962, 372.483022, 259.363014)
#Combine the data with Assembler
data3 = data.frame(Total_length, Assembler)
# melt the data to long format
data3 <- melt(data3, id.vars = 'Assembler')

# check the data
data3

# plot a simple bargraph
ggplot(data3, aes(x=Assembler, y=value, fill=Assembler)) + 
  geom_bar(stat='identity', width = 0.5, position='dodge') +
  labs(y = "Mb") +
  ggtitle("Total length") +
  theme_bw(base_size = 20) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))

# lets make another data set and  melt as before #For busco
Complete_BUSCO = c(31.68, 50.5, 0, 19.47, 38.61, 18.81, 9.57)
Partial_BUSCO = c(23.43, 18.48, 0.33, 15.51, 18.81, 9.24, 5.94)
data4 = data.frame(Complete_BUSCO, Partial_BUSCO, Assembler)        #merging the data
data4 <- melt(data4, id.vars = 'Assembler')                         #changing the data format
data4

# let's plot a bar graph
ggplot(data4, aes(x=Assembler, y=value, fill=variable)) +
  geom_bar(stat='identity', position='dodge') +
  labs(y = "%") +
  ylim(0, 75) +
  ggtitle("RAW assembly BUSCO score") +
  theme(legend.title = element_blank())+
  theme_bw(base_size = 20) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label=round(value, 2), fontface="bold", vjust=2), position=position_dodge(width=0.9), size=3.5)

# Now let's import another data in .csv file from computer
setwd("/path/to/working/directory")                               #setup working directory
getwd()                                                           #print the working directory
list.files()                                                      #list all the files in working directory

data5 <- read.csv("buscodata.csv", header = TRUE, sep = ",")      #import the data

# change the variable names
data5 <- within(data5, levels(variable_1)[levels(variable_1) == "No_trimming"] <- "No trim")
data5 <- within(data5, levels(variable_1)[levels(variable_1) == "After_trimming"] <- "Trimmed")
data5

# Lets plot a bar graph
ggplot(data5, aes(x=variable_1, y=value, fill=factor(variable_2, c("Complete_BUSCO", "Partial_BUSCO")))) +
  geom_bar(position=position_dodge(width=0.9), stat="identity") +
  facet_grid(~Assemblers) +
  geom_text(aes(label=round(value, 2), fontface="bold", vjust=2), position=position_dodge(width=0.9), size=3.5) +
  theme_bw(base_size = 20) +
  labs(y = "%")+
  theme(legend.title=element_blank(), legend.position="bottom",
        axis.text.y=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_text(face="bold"), axis.ticks.x=element_blank(), 	axis.ticks.y=element_blank(),
        panel.background=element_blank(), panel.grid.minor=element_blank(), panel.grid.major=element_blank())

# let's import another data set
data6 <- read.csv("finalstat.csv",header = TRUE, sep = ",")
data6

# let's plot another bargraph
ggplot(data6, aes(x=Assemblers, y=X., fill=BUSCO)) +
  geom_bar(stat='identity', position='dodge') +
  labs(y = "%") +
  ggtitle("BUSCO score") +
  theme_bw(base_size = 20) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label=round(X., 2), fontface="bold", vjust=2), position=position_dodge(width=0.9), size=3.5) +
  geom_text(aes(y=5, label=Polishing, fontface="bold"), color="white", position=position_dodge(width=0.9), angle=90, size=4)
