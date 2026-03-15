# Birds of Norfolk: Web Development and Data Visualisation
This repo contains the code and resources used in developing the new Birds of Norfolk website.

Here follows information about the contents of the repo.

## Recent Observations Map
The `Observation Visualisations.html` file contains the code used to generate
a sightings map for recent observations of birds in Norfolk. The map is 
searchable by species and displays the locations of all eBird checklists 
reporting the species which have been submitted in Norfolk in the last 30 days.

### Data Collection and Wrangling
The data is sourced directly from [eBird](https://ebird.org/home), a citizen science project which
enables users to submit birds they see whilst out birding. The `API Data Extraction.R` script
handles the data extraction process by:

- Calling the [eBird API](https://documenter.getpostman.com/view/664302/S1ENwy59) through the `rebird`
wrapper package to extract observation data
- Adjusting species names to the UK taxonomy
- Grouping all observations to the nearest eBird Hotspot within 1km
- Clustering remaining observations using `dbscan` and assigning names using
reverse geocoding via the `tidygeocoder` package.

Observations not assigned to a cluster in the last step are removed; this is to 
preserve observer privacy as many remaining observations are people's private
home addresses.

The data collection script runs every day at 00:00 UTC via GitHub Actions.

