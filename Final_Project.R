######################################
# BIOL 465 Research Project
# Eric von Amsberg, Katelyn Cline
# Fall 2020
# Temperature and NDVI as Explanatory Variables of Beetle, Spider, and Caterpillar Abundance in the Continental US
######################################

library(dplyr)
library(lubridate)
library(raster)
library(readr)
library(MODISTools)

# read in data
# source: https://caterpillarscount.unc.edu/dataDownload/
# The Caterpillars Count! data set has arthropod observations from all over the US over several years.
bugs = read.csv("CaterpillarsCountData.csv", header = TRUE)

# exploring the data, finding out how many sites there are, how many orders of bugs are found, and which dates have the most observations to make decisions about what to analyze
bugs_site = arrange(count(bugs, SiteName), desc(n))
bugs_orders = arrange(count(bugs, UpdatedArthropodGroup), desc(n))
bugs_dates = arrange(count(bugs, LocalDate), desc(n))

# change date from type character to type Date so we can filter by date
bugs$LocalDate = as.Date(bugs$LocalDate, format = "%Y-%m-%d")

# after meeting with Dr. Hurlbert and looking through the data, we decided to analyze the abundance data of spiders, caterpillars, and beetles in June 2019. These three orders are the most commonly observed orders in the whole dataset. June 2019 has lots of observations and was recommended by Hurlbert.

# filter Caterpillars Count! data to only records in June 2019
bugs_june19 = bugs %>%
  filter(month(bugs$LocalDate, TRUE) == "Jun", 
         year(bugs$LocalDate) == "2019") 

# count and arrange the list of sites in June 2019 to see how many observations each site has. n is the total number of observations for the given site in June 2019.
sites_june19 = arrange(count(bugs_june19, SiteName), n)
# names(sites_june19)[2] = 'totalObservationsinJune19' ---> tried renaming the column but too much depends on it later on

# determine the number of distinct branches at each site so we can ultimately calculate abundance per branch at each site
branches = bugs_june19 %>%
  filter(UpdatedArthropodGroup %in% c("caterpillar", "spider", "beetle")) %>%
  group_by(SiteName) %>%
  count(SurveyLocationCode) %>%
  group_by(SiteName) %>%
  summarise(numBranches = sum(n))

# join list of sites with their lat/long
bugs_lat_long = dplyr::select(bugs, c('SiteName', 'Latitude', 'Longitude'))
sites_june19 = left_join(sites_june19, bugs_lat_long, by = 'SiteName')
sites_june19 = distinct(sites_june19)
# so sites_june19 is the list of the sites where observations were made in June 2019 and their latitude and longitude

# filter bugs to observations of caterpillars, beetles, and spiders in June 2019
caterpillars_june19 = bugs %>%
  filter(month(bugs$LocalDate, TRUE) == "Jun", 
         year(bugs$LocalDate) == "2019",
         UpdatedArthropodGroup == "caterpillar") 

beetles_june19 = bugs %>%
  filter(month(bugs$LocalDate, TRUE) == "Jun", 
         year(bugs$LocalDate) == "2019",
         UpdatedArthropodGroup == "beetle")

spiders_june19 = bugs %>%
  filter(month(bugs$LocalDate, TRUE) == "Jun", 
         year(bugs$LocalDate) == "2019",
         UpdatedArthropodGroup == "spider")

#*****************AVERAGE TEMP IN JUNE 2019*****************
#defines spatial resolution
library(raster)
tavg= getData("worldclim", var="tmean", res=10)

#match site corrds with intersection of climate raster
coords = sites_june19[, c('Longitude', 'Latitude')]
points = SpatialPoints(coords, proj4string = tavg@crs)

#extract pulls climate data at each coordinate. stack of 12 layers one for each month
tavg_data = extract(tavg, points)
head(tavg_data)

#selects only the values for june. div by 10 bc decimal temp
june19_tavg = tavg_data[,6]/10

#creates df with temp data joined
sites_june19= cbind.data.frame(sites_june19,june19_tavg = june19_tavg)


#*****************NDVI DATA*****************
# working with MODIS
modisSites = sites_june19
modisSites = modisSites %>% 
  dplyr::rename(site_name = SiteName,
         lat = Latitude,
         lon = Longitude) %>%
  dplyr::select(-c(n))

# retrieve MODIS data
ndvi = mt_batch_subset(df = modisSites, product = "MOD13Q1", band = "250m_16_days_NDVI", start = "2019-06-01", end = "2019-06-16")
# refine big modis df to just site and ndvi value
ndvi_small = ndvi %>% 
  dplyr::select(c('site', 'value')) %>%
  dplyr::rename(NDVI = value)

#rescale NDVI -1 < x < 1
ndvi_small$NDVI= ndvi_small$NDVI*.0001
# join to site data so sites_june19 has all the site specific data
sites_june19 = left_join(sites_june19, ndvi_small, by = c('SiteName' = 'site'))

# 2 points are negative which usually indicates water instead of vegetation. Using Google Maps, I found the sites via their lat/long and they are very close to the coast. I'm going to remove them from the data set because they create huge outliers and are incorrect ndvi values for the location.
sites_june19 = sites_june19[-c(2, 28),]

#*****************ABUNDANCE*****************
# calculate abundance data for each order in each location for the whole month of june
cat_abundance = caterpillars_june19 %>%
  group_by(SiteName) %>%
  summarise(caterpillarAbundance = sum(ArthropodQuantity))

beetle_abundance = beetles_june19 %>%
  group_by(SiteName) %>%
  summarise(beetleAbundance = sum(ArthropodQuantity))

spider_abundance = spiders_june19 %>%
  group_by(SiteName) %>%
  summarise(spiderAbundance = sum(ArthropodQuantity))

# join these abundance values to sites dataframe
sites_june19 = left_join(sites_june19, cat_abundance, by = 'SiteName')
sites_june19 = left_join(sites_june19, beetle_abundance, by = 'SiteName')
sites_june19 = left_join(sites_june19, spider_abundance, by = 'SiteName')
# turn all the NAs into zero
sites_june19[is.na(sites_june19)] = 0

# join branches to sites_june19 & then divide abundance values by branch values
sites_june19 = left_join(sites_june19, branches, by = 'SiteName')
sites_june19$corrCatAbundance = sites_june19$caterpillarAbundance / sites_june19$numBranches
sites_june19$corrBeetleAbundance = sites_june19$beetleAbundance / sites_june19$numBranches
sites_june19$corrSpiderAbundance = sites_june19$spiderAbundance / sites_june19$numBranches

#**************Variance Partitioning***************
#caterpillar
lm.cat.ndvi = lm(corrCatAbundance ~ NDVI, data = sites_june19)
lm.cat.temp = lm(corrCatAbundance ~ june19_tavg, data = sites_june19)
lm.cat.ndvi_temp = lm(corrCatAbundance ~ NDVI + june19_tavg, data = sites_june19)
cat.R1 = summary(lm.cat.ndvi)$r.squared
cat.R2 = summary(lm.cat.temp)$r.squared
cat.R12 = summary(lm.cat.ndvi_temp)$r.squared

#beetle
lm.beet.ndvi = lm(corrBeetleAbundance ~ NDVI, data = sites_june19)
lm.beet.temp = lm(corrBeetleAbundance ~ june19_tavg, data = sites_june19)
lm.beet.ndvi_temp = lm(corrBeetleAbundance ~ NDVI + june19_tavg, data = sites_june19)
beet.R1 = summary(lm.beet.ndvi)$r.squared
beet.R2 = summary(lm.beet.temp)$r.squared
beet.R12 = summary(lm.beet.ndvi_temp)$r.squared

#spider
lm.spid.ndvi = lm(corrSpiderAbundance ~ NDVI, data = sites_june19)
lm.spid.temp = lm(corrSpiderAbundance ~ june19_tavg, data = sites_june19)
lm.spid.ndvi_temp = lm(corrSpiderAbundance ~ NDVI + june19_tavg, data = sites_june19)
spid.R1 = summary(lm.spid.ndvi)$r.squared
spid.R2 = summary(lm.spid.temp)$r.squared
spid.R12 = summary(lm.spid.ndvi_temp)$r.squared


#***************** PLOTS ********************
# set up plot window to make 2 x 3 view of graphs
par(mfrow = c(2, 3))
plot(sites_june19$NDVI, sites_june19$corrCatAbundance, xlab = "NDVI", ylab = "Abundance/Branches per Site", main = "Caterpillar NDVI vs Abundance", col = "springgreen4")
abline(lm.cat.ndvi, col = "springgreen4")
plot(sites_june19$NDVI, sites_june19$corrBeetleAbundance, xlab = "NDVI", ylab = "Abundance/Branches per Site", main = "Beetle NDVI vs Abundance", col = "sienna2")
abline(lm.beet.ndvi, col = "sienna2")
plot(sites_june19$NDVI, sites_june19$corrSpiderAbundance, xlab = "NDVI", ylab = "Abundance/Branches per Site", main = "Spider NDVI vs Abundance", col = "royalblue1")
abline(lm.spid.ndvi, col = "royalblue1")
plot(sites_june19$june19_tavg, sites_june19$corrCatAbundance, xlab = "Avg Temp June 2019 (*C)", ylab = "Abundance/Branches per Site", main = "Caterpillar Temperature vs Abundance", col = "springgreen4")
abline(lm.cat.temp, col = "springgreen4")
plot(sites_june19$june19_tavg, sites_june19$corrBeetleAbundance, xlab = "Avg Temp June 2019 (*C)", ylab = "Abundance/Branches per Site", main = "Beetle Temperature vs Abundance", col = "sienna2")
abline(lm.beet.temp, col = "sienna2")
plot(sites_june19$june19_tavg, sites_june19$corrSpiderAbundance, xlab = "Avg Temp June 2019 (*C)", ylab = "Abundance/Branches per Site", main = "Spider Temperature vs Abundance", col = "royalblue1")
abline(lm.spid.temp, col = "royalblue1")
# reset plot window
par(mfrow = c(1,1)) 

#************************ANCOVA***************************

# rearrange data to fit ancova format from class
# Create 3 df one comparing cat to beetle, one beet to spider, one spider to cat

# caterpillar compared to beetle df
CatBee_sites_june19= cbind(rep(sites_june19[,1],2), rep(sites_june19[,5:6],1 ), stack(sites_june19[,11:12]), c(rep(0,nrow(sites_june19)), rep(1,nrow(sites_june19))))
names(CatBee_sites_june19) = c('site','temp','NDVI','abundance','taxon','group')

# temp abundance relationship 
ancova.catbeeT= lm(abundance ~ group + temp + temp * group, data=CatBee_sites_june19)
# NDVI abundance relationship 
ancova.catbeeN = lm(abundance ~ group + NDVI + NDVI * group, data=CatBee_sites_june19)

# beetle compared to spider 
BeeSpi_sites_june19= cbind(rep(sites_june19[,1],2), rep(sites_june19[,5:6],1 ), stack(sites_june19[,12:13]), c(rep(0,nrow(sites_june19)), rep(1,nrow(sites_june19))))
names(BeeSpi_sites_june19) = c('site','temp','NDVI','abundance','taxon','group')

# temp abundance relationship 
ancova.beespiT= lm(abundance ~ group + temp + temp * group, data=BeeSpi_sites_june19)
# NDVI abundance relationship 
ancova.beespiN = lm(abundance ~ group + NDVI + NDVI * group, data=BeeSpi_sites_june19)

# spider compared to caterpillar
SpiCat_sites_june19= cbind(rep(sites_june19[,1],2), rep(sites_june19[,5:6],1 ), stack(sites_june19[,c(13,11)]), c(rep(0,nrow(sites_june19)), rep(1,nrow(sites_june19))))
names(SpiCat_sites_june19) = c('site','temp','NDVI','abundance','taxon','group')

# temp abundance relationship
ancova.spicatT= lm(abundance ~ group + temp + temp * group, data=SpiCat_sites_june19)
# NDVI abundance relationship 
ancova.spicatN = lm(abundance ~ group + NDVI + NDVI * group, data=SpiCat_sites_june19)

# Results of ANCOVA 

# Temp Abundance Relationship 
summary(ancova.catbeeT)
summary(ancova.beespiT)
summary(ancova.spicatT)

#NDVI Abundance Relationship
summary(ancova.catbeeN)
summary(ancova.beespiN)
summary(ancova.spicatN)

#B0 = Intercept = Intercept of first bug order
#B1 = group = How similar is intercept for other bug order
#B2 = temp  = Slope of first bug order 
#B3 = group:temp = difference in slope compared to first order

#*********NDVI and temp relationship plot*************
plot(sites_june19$june19_tavg, sites_june19$NDVI, xlab = "Avg Temp June 2019 (*C)", ylab = "NDVI", main = "Temperature vs NDVI", col = "darkmagenta")
enviro = lm(sites_june19$NDVI~ sites_june19$june19_tavg)
abline(enviro, col = "darkmagenta")

plot(sites_june19$Latitude, sites_june19$Longitude, main = "Lat vs Long", xlab = "Latitude", ylab = "Longitude", col = "chartreuse4")
