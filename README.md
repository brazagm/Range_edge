Dataset: Morphological variability decreases in populations living in less suitable environments and close to the range edges
-----

*Contextualization*: This dataset contains data and R codes to investigate the spatial patterns of the intraspecific variability of skull traits throughout the geographic range of a marsupial from the Brazilian Atlantic Forest, the Gray slender opossum (<i>Marmosops incanus</i>). We aimed to answer whether the distance from the range edge and the environmental suitability explain the geographic variation of morphological variability of the species.

*Methods*: We analysed specimens deposited in the main biological collections in southeastern Brazil. The morphological traits were characterised by linear measures of the skull and environmental suitability was estimated by Ecological Niche Models. Distance from the range edge was calculated based on the minimum linear distance between populations and the closest range limit defined by a Minimum Convex Polygon including all the occurrence records of the species. Our hypotheses were tested by fitting linear regressions and selecting the best models that explained the spatial variation of morphological variability based on the Akaike information criterion.

*Results*: Morphological variability of <i>M. incanus</i> increases with local environmental suitability and distance from the range edges. However, the relationship between morphological variability and environmental suitability depends on the niche modelling method.

![<font size="2">*Upper and lower skull view and lateral view of the mandible from a specimen of Gray slender opossum, <i>Marmosops incanus</i>.</font>*](https://github.com/brazagm/Range_edge/blob/master/README_skull.png?raw=true)


### 1. Description of the data and file structure

The project follows a simple organization structure in which the R codes are found in sources directory `./src` and use the other directories in which input files `./data`, intermediate output files `./processed_data` and final output files `./results` are located or created.  

Only the input data is provided here! All outputs are produced by scripts located in `./src` and they are not provided here. Directory/subdirectory structures are provided for the appropriate code behavior.

<pre>
Range_edge/  
├── README.md             # overview of the project  
├── data/                 # data files used in the project  
│   └── shapes/           # shape files used in the analyses
│   └── layers/           # raster files used in the analyses
├── processed_data/       # intermediate files from the analysis  
├── results/              # results of the analysis (data, tables, figures)
│   ├── figures/          # figures for results  
│   └── niche_models/     # niche model results in raster files  
└── src/                  # contains all codes in the project  
</pre>

<br>

#### 1.1. `./data` directory and files

This directory is where all the input data needed for the project's analyses are located -- i.e., files located here come from, for example, spreadsheets with data collected in fieldwork, measurements of study objects, data from literature, and so on...

Occurrence records information `./data/data_records.csv` and morphological data `./data/data_morphology.csv` for *M. incanus* are the main input files for this project. Occurrence records of the species were obtained from the specimens deposited in the main biological collections from southeastern Brazil and from the 'Abundance of small mammals in the Atlantic Forest' dataset, a.k.a. ASMAF ([Figueiredo et al. 2017](http://doi.wiley.com/10.1002/ecy.2005)). Morphological data were obtained by visiting the main biological collections from southeastern Brazil to collect linear measures of the skull of *M. incanus*'s specimens.

**Important**: *Detailed information about collecting and organizing the information contained in spreadsheets and other input files can be found in the metadata file within this directory* `./data/metadata.txt`. *Empty cells must be considered as NA values. Empty cells were used in cases of missing information or measures that could not be taken. Scripts in this project are already prepared to work with this empty cells as NAs.*

Shapefiles for the world borders polygons are located within shapes directory `./data/shapes` and were downloaded from the World Border Dataset in [Thematic Mapping API](https://thematicmapping.org/downloads/world_borders.php).

**Important**: *World borders shapefile is not included in this project. Shape files must be downloaded [here](https://thematicmapping.org/downloads/world_borders.php) and placed in the project directory for them *`./data/shapes`.

Raster files for ecological niche modelling are located within layers directory `./data/layers`. Raster are organized in two types of environmental layers: (1) `./data/layers/WorldClim25_Global_v2` contains the 19 Bioclimatic variables in a 2.5 arc-minutes resolution downloaded from the [WorldClim version 2](https://www.worldclim.org/data/worldclim21.html); and (2) `./data/layers/MOD13A13`
contains the productivity variable, the sum of positive values of the integrative normalised difference vegetation index (INDVI), derived from the [MOD13A3 v006](https://lpdaac.usgs.gov/products/mod13a3v006/) data product.

**Important**: *WorldClim's layers are not included in this project. Raster files must be downloaded [here](https://www.worldclim.org/data/worldclim21.html) and placed in the project directory for them *`./data/layers/WorldClim25_Global_v2`*. The sum of positive values of the integrative normalised difference vegetation index (INDVI) raster *`./data/layers/MOD13A3/MOD13A3_INDVI.tif` *was calculated based on the average of NDVI monthly values available in the MOD13A3 v006 dataset at a 30 arc-seconds resolution. INDVI was then resampled through bilinear interpolation to a 2.5 arc-minutes resolution. Please, for more details see the [published article](http://doi.wiley.com/10.1111/jbi.14687) regarding this project.*

<br>

#### 1.2. `./processed_data` directory and files
This directory is where all the intermediate outputs created by the project's analyses are located.

File names located in this folder are identified by a number that refers to their parent code. For example: files whose name begins with '01' are created by the '01_morpho_analyses' script `./src/01_morpho_analyses.R`, and those beginning with '03' are created by the '03_statistical_analyses' script `./src/03_statistical_analyses.R`. 

**Important**: *All intermediate outputs are produced by scripts located in* `./src` *and they are not provided here. Only directory/subdirectory structures are provided for appropriate code behavior. The analytical steps for producing these outputs are described in detailed comments found in the codes.*

<br>

#### 1.3. `./results` directory and files
This directory is where all the outputs for the project's results are located.

File names located in this folder are identified by a number that refers to their parent code. For example: files whose name begins with '01' are created by the '01_morpho_analyses' script `./src/01_morpho_analyses.R`, and those beginning with '03' are created by the '03_statistical_analyses' script `./src/03_statistical_analyses.R`. 

`./results` directory contains two folders for outputs: (1) `./results/figures` contains all the figures produced by the analyses; and (2) `./results/niche_models`
contains all the Ecological Niche Modelling results, including raster files with the environmental suitability for each algorithm (.tif files) and performance metrics for them (.csv files).

**Important**: *All final outputs are produced by scripts located in* `./src` *and they are not provided here. Only directory/subdirectory structures are provided for appropriate code behavior. The analytical steps for producing these outputs are described in detailed comments found in the codes.*

<br>

#### 1.4. `./src` directory and files

This directory contains all the codes needed for the project. The analyses are structure in three R scripts:

  - `01_morpho_analyses.R` -- morphological data preparation and analyses. Import occurrence records, populations and trait measures of the species. Export measurement errors for each trait, test sexual dimorphism and calculate trait measures corrected by size.
  
  - `02_niche_modeling.R` -- Ecological Niche Modeling procedures to estimate environmental suitability for the speices. Same script version from *Braz et al. (2020) Interspecific competition constrains local abundance in highly suitable areas. Ecography 43: 1560–1570*.
  
  - `03_statistical_analyses.R` -- statistical analyses to test the relationship between morphological variability, environmental suitability and distance from the range edge throughout the geographic range of the Gray slender opossum, *Marmosops incanus*. Models describing the relationship between dependent-independent variables were selected using the Akaike criteria.
  
  - `function_corner_text.R` -- function created by Eran Raviv 2015 for facilitate addition of textboxes located by the corners in R plots.

<br>

## Sharing/Access information
Please, in case of any doubt feel free to contact the corresponding author at brazagm@gmail.com. We kindly request that the use of the data and/or codes available here, modified or not, be communicated to the authors. 

Links to other publicly accessible locations of the data:  

  - [Git repository](https://github.com/brazagm/Range_edge)  
  - [Figshare](https://figshare.com/articles/dataset/Marmosops_incanus_occurrence_records_and_morphological_data/22790000)
  
Related research:  

  - Braz, A.G.; Figueiredo, M.S.L.; Weber, M.M.; Grelle, C.E.V. 2023 Morphological variability decreases in populations living in less suitable environments and close to the range edges. *Journal of Biogeography*, *in prep.*. doi: [10.1111/jbi.14687](http://doi.wiley.com/10.1111/jbi.14687)

