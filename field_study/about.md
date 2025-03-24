## Field Study
This section contains the field study conducted along the coastline of New Jersey. It has the data aggregated and processed as well as the code (in Python) utilized to plot our results.

### Field Data
The data used in this analysis is located in the excel and csv file. The excel was included to incorporate the excel functions used as well as for enhanced readability. The csv file was created from the excel file. Underlying data in the excel were obtained from:
* Landsat 7 & 8 satellite imagery collected from USGS Earth Explorer (more info below)
* Zillow ZHVI data from Zillow Housing Research
* S&P CoreLogic Case-Shiller NY-New York Home Price Index from the Federal Reserve Bank of St. Louis
* Beach nourishment data from The Program for the Study of Developed Shorelines at Western Carolina University

### Couplets Study
The python code in couplets_study_final.py is used to analyze the field data in the csv file field_study_couplets.csv. 

### Geospatial Data
The geopackages containing digitized beach widths used in this field data are included in this repository. Beach widths were digitized in QGIS 3.28 2-Firenze using Landsat 7 & 8 imagery from USGS Earth Explorer. The satellite images utilized were:
* LC08_L1TP_013032_20200613_20200823_02_T1
* LC08_L1TP_014033_20200722_20200911_02_T1
* LE07_L1TP_014032_20000808_20200917_02_T1
* LE07_L1TP_014033_20000707_20200918_02_T1
* LE07_L1TP_014032_19990806_20200918_02_T1
* LE07_L1TP_014033_19990705_20200918_02_T1
