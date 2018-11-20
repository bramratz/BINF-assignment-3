#Going to be investigating whether closely related species live in the smae or geographically close locations, or whether species are spread randomly geographically regardlss of relateness. As a secondary question I'm going to look at what proportion of Daphnia species were sampled from Canada.

#install necessary packages
install.packages("tidyverse")
library(tidyverse)
install.packages("readr")
library(readr)

#Download the Daphnia data file from the BOLD database
DaphniaC <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Daphnia&format=tsv")

#Write the file to hard disk
write_tsv(DaphniaC, "DaphniaC_BOLD_data.tsv")

#Read file from hard disk 
DaphniaC <- read_tsv("DaphniaC_BOLD_data.tsv")

#Summary of the data
summary(DaphniaC)

#More information.
str(DaphniaC)

#column names. What types of data is available 
names(DaphniaC)

#geographic data available
#how many unique countries are in the data set 
unique(DaphniaC$country)

#what do the fist 100 entries look like, are there NA's?
head(DaphniaC$country, 100)

#how many unique species names are present in the data set
length(unique(DaphniaC$species_name))

#how many unique accession numbers are present in the data 
length(unique(DaphniaC$genbank_accession))

#make df
DaphniaC_df <- DaphniaC[,c("species_name", "genbank_accession", "country", "lat", "lon", "nucleotides")]
View(DaphniaC_df)

#remove NA's from the dataframe
DaphniaC_df2 <- DaphniaC_df[complete.cases(DaphniaC_df), ]
View(DaphniaC_df2)

#compare accession numbers from the df generated from BOLD with the list of species from part B once the sequences were clustered into OTU. looking for overlap in accession numbers in order to look at which species present i can make a phylogeny for in order to look at relatedness between species

#install necessary packages
source("http://bioconductor.org/biocLite.R")
biocLite("Biostrings")
library(Biostrings)
biocLite("muscle")
library(muscle)
biocLite("DECIPHER")
library(DECIPHER)
biocLite("ape")
library(ape)

#first have to run alignment on the data with muscle 
#preform a multiple sequence alignment
#first convert df to DNAstringset
DaphniaC_convert <- DNAStringSet(DaphniaC_df2$nucleotides)

#checking that the class converted
class(DaphniaC_convert)

#Now, align sequences.
DaphniaC_alignment <- DNAStringSet(muscle(DaphniaC_convert, maxiters = 3, diags = TRUE))

#Writing to file
writeXStringSet(DaphniaC_alignment, file="DaphniaC_alignment.fas", format = "fasta", width = 1000)

#Adding acession numbers as names for the sequences.
names(DaphniaC_alignment) <- DaphniaC_df2$genbank_accession

#Viewing the alignmnet
DaphniaC_alignment

#converting data
DaphniaC_DNAbin <- as.DNAbin(DaphniaC_alignment)

#Creating a distance matrix to make dendrogram using the default N model
DaphniaC_distancematrix <- dist.dna(DaphniaC_DNAbin, model = "N", as.matrix = TRUE, 
                                       pairwise.deletion = TRUE)

#Clustering into OTUs
DaphniaC_cluster <- IdClusters(DaphniaC_distancematrix,
                                  method = "single",
                                  cutoff= 0.03,
                                  showPlot = FALSE,
                                  type = "clusters",
                                  processors = 2,
                                  verbose = TRUE)

#create dendrogram
DaphniaC_dendrogram <- IdClusters(DaphniaC_distancematrix,
                                     method = "single",
                                     cutoff= 0.03,
                                     showPlot = TRUE,
                                     type = "both",
                                     processors = 2,
                                     verbose = TRUE)


#going to uses species present in phylogeny to see location in the world and plot this on a map using ggplot and maps. After going to look at the cloest related species that are furtherest apart in giegraphical distance vs. closest species geographically while being furthest apart in terms of relateness. look for dissimilarity analysis to see why this is the case or if anything interesting arises from this.
#create maps
#load libraries
install.packages("raster")
install.packages("rgdal")
install.packages("dismo")
install.packages("XML")
install.packages("maps")
install.packages("maptools")
install.packages("ggplot2")
library(raster)
library(rgdal)
library(dismo)
library(XML)
library(maps)
library(maptools)
library(ggplot2)

#show world map
map()

#using lat and lon locate where individuals are found on the map
DaphniaC_map <- points(DaphniaC_df2$lon, DaphniaC_df2$lat, col = "red", pch = 19, xlab = "Longitude", ylab = "Latitude")


#zoom in
DaphniaC_map_canada <- plot(DaphniaC_df2$lon, DaphniaC_df2$lat, xlim = c(-160, -25), ylim = c(20, 90), col = "red", pch = 19, xlab = "Longitude", ylab = "Latitude")
map(add = T)

#try again to create base world map
mp <- NULL
mapworld <- borders("world", colour = "grey50", fill = "grey50")
mp <- ggplot() + mapworld

#layer species on top
mp <- mp + geom_point(aes(x = DaphniaC_df2$lon, y = DaphniaC_df2$lat), color = "blue", size = 3)
mp

#load the ggmap package
devtools::install_github("dkahle/ggmap")
#load libraries 
library(stringr)
install.packages("viridis")
library(viridis)

#add data points
gg1 +
  geom_point(data = DaphniaC_df2, aes(lon, lat), color = "red", size = 5) +
  ggtitle("Species Distribution") +
  geom_text(data = DaphniaC_df2, aes(lon, lat, label = species_name))


#The map shows the distrubution of Daphnia over the world. based on this we can see that Daphnia is located on every continent with the exception of south america. It is also clear that the majority live north of the equator. 
#I wasn't able to assign data labels to the map in a way that allows me to pin point which species lives where which would alllow me to determine the species distribution, showing if more close related species live geographically closer to each other or not. The next step in this would be to preform a species accuulation curve by country to determine species richness. The end goal of this would be to determine which country has the greatest diversity of daphnia. further research could look at whether this diversity is effected by habitat diversity.
