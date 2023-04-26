# Estancia Suitability Analysis
This repository contains code used to identify geographic trends in the remains of Andean *estancias*. This project was active from January-May 2023 and was completed towards the fulfillment of the Vanderbilt Masters in Data Science capstone requirement.

# Repository Organization

This repository contains Jupyter Python notebooks used to perform analysis, as well as one pure python file hs_geofuncs.py which provides utility functions that are used in the notebooks. Code was developed for use on the Microsoft Planetary Computer and will not work outside this platform.

- estancias_clim.ipynb is a notebook exploring attaching climate data to estancia locations. This direction was not ultimately pursued in the final analysis.
- estancias_dem.ipynb is a notebook exploring atteaching elevation data and derived terrain data to estancia locations - the functions and approach developed in this notebook are used in the latter half of point_analysis.ipynb
-estancias_ndvi.ipynb is a notebook exploring attaching vegetation data to estancia locations - the function and approach developed in this notebook are used in the "NDVI" section of point_analysis.ipynb
-point_analysis.ipynb represents the ultimate analysis, pulling together approaches from the previous three notebooks into a cohesive overview of estancia distributions
-logistic_regression.ipynb is a secondary approach that uses a logistic regression model to discriminate between real and fake estancia locations. 
-hs_geofuncs.py is a collection of python code used to facilitate the use of Microsoft Planetary Computer's data catalog to extract geospatial data to points.

# Background

The remains of historic ranches, or “estancias” as they are called in the local Spanish, dot the landscape of the highland Andes of modern-day Peru. Spanning a temporal extent of over 400 years, these estancias represent a longstanding history of pastoralism in post-Inca society. By studying these estancias, researchers can gain insights into the complex interplay between agricultural and pastoralist lifestyles in colonial times, providing a better understanding of the social and economic dynamics of the past. The Geospatial Platform for Andean Culture, History, and Archaeology (GeoPACHA) facilitates these pursuits as a cross-institutional campaign to provide a comprehensive documentation of archaeological sites in the Andean region.

# Data

The data for this project was manually collected via the GeoPACHA platform. The data itself is not public and has not been made available in this repository. It is referred to in the code as "Estancias_corrals_from_sat_imagery.csv","Southwest_survey_polygon.shp", and "wernke_logres.csv". The data consists of the longitude-latitude coordinates of 9,500 estancias in the survey region.

<img width="214" alt="image" src="https://user-images.githubusercontent.com/46304188/234697706-6a266109-ba0c-40a7-9eee-f4dbccabb823.png">




