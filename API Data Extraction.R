#Setting up environment
require(ggplot2)
require(dplyr)
require(rebird)
require(geosphere)
require(pbapply)

#setwd('path')
rm = list(ls())

#Group all hotspots by site and get aggregated co-ordinates
locsAll <- ebirdhotspotlist("GB-ENG-NFK") |>
  select(locId, locName, locLat = lat, locLng = lng) |>
  
  #Hotspots at the same site are often named in format [site name]--[specific
  #subsection of site] or [site name]([specific stipulation]), so grouping
  #hotspots by site is done by removing strings after -- or (; where distant 
  #sites share names, this is handled in nonDupe.
  mutate(nonDupe = duplicated(locName) | duplicated(locName, fromLast = T),
         siteName = trimws(sub("--.*|\\(.*","",locName))) |>
  
  #Removing historical stakeouts
  filter(!grepl("stakeout|Stakeout",siteName)) |>
  
  #Manual corrections
  mutate(siteName = ifelse(
    siteName == "Th Broads NP",
    "The Broads NP",
    siteName
  )) |>
  
  #Calculation of aggregated co-ordinates
  group_by(siteName) |>
  mutate(
    siteLat = ifelse(n()>1 & !nonDupe, mean(locLat), locLat),
    siteLng = ifelse(n()>1 & !nonDupe, mean(locLng), locLng)
  ) |>
  ungroup()

#Get English (UK) names for taxa instead of English (US)
ukNames <- ebirdtaxonomy(locale = "en_UK") |>
  select(speciesCode, comName)

#Generating a list of all species recorded in Norfolk in the last 30 days
spAll <- ebirdregion("GB-ENG-NFK", provisional = T, back=30) |>
  
  #Rolling observations at subspecies level into parent species
  left_join(ebirdtaxonomy() |> select(speciesCode, reportAs, category), by = "speciesCode") |>
  mutate(speciesCode = ifelse(
    category == "subspecies" & !is.na(reportAs),
    reportAs,
    speciesCode
  )) |>
  
  #Changing names from USA output to UK one
  select(-comName) |>
  left_join(ukNames, by = "speciesCode") |>
  
  #Filtering world list to remove hybrids and aggregates
  filter(speciesCode %in% (ebirdtaxonomy() |>
           filter(category == "species") |>
           pull(speciesCode))) |>
  pull(speciesCode)

#Compute data frame of most recent observations of all species across all sites
#in Norfolk.
obsAll <- pblapply(spAll, function(sp){
  ebirdregion("GB-ENG-NFK",sp, provisional = T, back=30) |>
    
    #Adjust to English (UK) names
    select(-comName) |>
    left_join(ukNames, by = "speciesCode") |>
    
    #Group all non-hotspot checklists to nearest hotspot within 2km
    select(comName, sciName, locId, locName, obsDt, obsValid, lat, lng) |>
    
    rowwise() |>
    mutate(siteName = ifelse(
      locId %in% locsAll$locId,
      locsAll$locName[match(locId, locsAll$locId)],
      ifelse(
        min(distHaversine(c(lng, lat), cbind(locsAll$locLng, locsAll$locLat))) <= 2000,
        locsAll$locName[which.min(distHaversine(c(lng, lat), cbind(locsAll$locLng, locsAll$locLat)))],
        NA_character_
      )
    )) |>
    ungroup() |>
    
    #Get site co-ordinates from locsAll, else leave unchanged
    mutate(
      siteLat = ifelse(
        locId %in% locsAll$locId,
        locsAll$siteLat[match(locId, locsAll$locId)],
        lat
      ),
      
      siteLng = ifelse(
        locId %in% locsAll$locId,
        locsAll$siteLng[match(locId, locsAll$locId)],
        lng
      )
    ) |>
    select(-c(locId, locName, lat, lng)) |>
    filter(!is.na(siteName))
}) |>
  bind_rows()

write.csv(obsAll, paste0("All Obs ", format(Sys.Date(), "%d%m"), ".csv"), row.names = F)
