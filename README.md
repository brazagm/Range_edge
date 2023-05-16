# Morphological variability throughout the geographic range of the Gray slender opossum  
 
This project aims to investigate the spatial patterns of the intraspecific variability of skull traits throughout the geographic range of a marsupial from the Brazilian Atlantic Forest, the Gray slender opossum (<i>Marmosops incanus</i>). We aimed to answer whether the distance from the range edge and the environmental suitability explain the geographic variation of morphological variability of the species.

![Upper and lower skull view and lateral view of the mandible from a specimen of Gray slender opossum, <i>Marmosops incanus</i>.](https://github.com/brazagm/Range_edge/blob/master/README_skull.png?raw=true)

### Project structure
<pre>
Range_edge/  
├── README.md             # overview of the project  
├── data/                 # data files used in the project  
│   └── shapes/           # shape files used in the analyses
├── processed_data/       # intermediate files from the analysis  
├── results/              # results of the analysis (data, tables, figures)
│   ├── figures/          # figures for results  
│   └── niche_models/     # niche model results in raster files  
└── src/                  # contains all code in the project  
</pre>

### Methods
We analysed specimens deposited in the main biological collections in southeastern Brazil. The morphological traits were characterised by linear measures of the skull and environmental suitability was estimated by Ecological Niche Models.
The analyses were structure in three R scripts:
<pre>
#### 01_morpho_analyses.R
The morphological data preparation and analyses. Import occurrence records, populations and trait measures of the species. Export measurement errors for each trait, test sexual dimorphism and calculate trait measures corrected by size.   

#### 02_niche_modeling.R
Ecological niche modeling procedures to estimate environmental suitability for the speices. Same script version from <i>Braz et al. (2020) Interspecific competition constrains local abundance in highly suitable areas. Ecography 43: 1560–1570</i>.

#### 03_statistical_analyses.R
Statistical analyses to test the relationship between morphological variability, environmental suitability and distance from the range edge throughout the geographic range of the Gray slender opossum, <i>Marmosops incanus</i>. Models describing the relationship between dependent-independent variables were selected using the Akaike criteria.
</pre>
