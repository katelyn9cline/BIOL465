library(dplyr)

# read in data
bbs = read.csv('HW9_NC_warblers.csv', stringsAsFactors = FALSE)
mass = read.csv('HW9_avian_body_masses.csv', stringsAsFactors = FALSE)

# Step 1: number of sites occupied by each species
occupancy = count(bbs, SpeciesName)
names(occupancy)[2] = 'occ'  

  # choose species
  chestnut_warb = bbs[bbs$SpeciesName == 'Chestnut-sided Warbler', ]
  chestnut_warb = filter(bbs, SpeciesName == 'Chestnut-sided Warbler')
  amer_red = filter(bbs, SpeciesName == 'American Redstart')
  
  #count species richness for route
  route = count(bbs, Route)
  
# Step 2: avg pop density across all sites
bbs_by_species = group_by(bbs, SpeciesName)
pop = summarize(bbs_by_species, meanN = mean(Abundance))

max_route = group_by(bbs, Route)
max_abun = summarise(max_route, max(Abundance))

# Step 3: join
birds1 = left_join(occupancy, pop, 
                   by = c('SpeciesName' = 'SpeciesName'))
birds1 = left_join(occupancy, pop, by = 'SpeciesName')
birds2 = left_join(birds1, mass, by = c('SpeciesName' = 'CommonName'))

pop = bbs %>%
  group_by(SpeciesName) %>%
  summarize(meanN = mean(Abundance)) %>%
  left_join(mass, by = c('SpeciesName' = 'CommonName'))

max_abun_pipe = bbs %>%
  group_by(Route) %>%
  summarise(maxAbund = max(Abundance), 
            maxSpecies = SpeciesName[Abundance == maxAbund])

# Step 4: plotting

graph_data = left_join(pop, occupancy, by = 'SpeciesName')
par(mfrow = c(1, 2), mar = c(5, 5, 1, 1))
plot(graph_data$Mass_g, graph_data$occ, 
     xlab = 'Bird Mass (g)', 
     ylab = 'Site Occupancy',
     main = 'Mass vs Occupancy for NC Warblers')
plot(graph_data$Mass_g, graph_data$meanN,
     xlab = 'Bird Mass (g)',
     ylab = 'Abundance',
     main = 'Mass vs Abundance for NC Warblers')
lm.occ = lm(graph_data$occ ~ graph_data$Mass_g)
summary(lm.occ)
lm.abund = lm(graph_data$meanN ~ graph_data$Mass_g)
summary(lm.abund)

par(mfrow = c(1, 1))
plot(graph_data$occ, graph_data$meanN, 
     xlab = 'Site Occupancy',
     ylab = 'Abundance',
     main = 'Occupancy vs Abundance for NC Warblers')
lm.occ_abund = lm(graph_data$meanN ~ graph_data$occ)
summary(lm.occ_abund)
