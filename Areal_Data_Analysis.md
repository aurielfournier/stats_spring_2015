# Areal Data Analysis
Jackson Cothren  
February 17, 2015  

# Introduction to areal data

Areal data differs from point pattern and continuous data in the form of the data itself. While continuous data involves. Point patterns can be thought of as samples from a continuous space (but where the lcoation of the event is of important) and continuos spatial distributions (such as temperature readings at various point locations) are estimated from samples. Areal, or lattice data, on the other hand, involves aggregated quantities within some relevant spatial partition of a given region (such as census tracts within a city, counties within a statem, watersheds in ecoregion, etc.). Often, the boundaries have little or nothing to do with the aggregated data, having been form for other reasons. Even though more and more non-aggregated data is collected (mobile phone locations and tweets, for example) much of the social science data available to us has been aggregated for anonymity reasons.   

## Issues with areal data

There are several issues you must considered when dealing with data that has been aggregated to areal units.

- __Modifiable areal unit problem (MAUP)__: results of statistical anlaysis may and often do depend on the specific geographic unit used in the study. Both size and shape of the units can have an effect. We'll see several examples when we discuss spatial autocorrelation measures.

- __Ecological fallacy__: results obtained from aggregated data cannot be assumed to apply to individual people. For instance, assume that you measured the math scores of a particular classroom and found that they had the highest average score in the district. Later you run into one of the kids from that class and you think to yourself "she must be a math whiz." That is a fallacy. Just because she comes from the class with the highest average doesn't mean that she is automatically a high-scorer in math. She could be the lowest math scorer in a class that otherwise consists of math geniuses.

- __Non-uniformity of space__: phenomena  are not distributed evenly in space.

- __Edge issues__: edges of the map, beyond which there is no data, can significantly affect results (we've seen this in point pattern analysis). 

- __Spatial autocorrelation__: data from locations near to each other are usually more similar than data from locations far away from each other. This essentially leads to observations which are not independent. 

Virtually all disciplines which deal with areal entities have dealt with these problems. Geographer's typically refer to it as _Tobler's Law_ (which only partially describes) while others will refer to _Galton's problem_. The problem is to establish how many effectively indpedendent observations are present, when arbitrary boundaries have been used to tesselate the study area. In 1889, Galton question questioned whether observations of marriage laws across areal entities constituted independent observations, since they could just reflect a general pattern from which they had all descended. So positive spatial dependence tends to reduce the amount of information contained in the observations, because the proximate observations can in part be used predict each other. 

## Spatial concepts

- __Distance__: the magnitude of spatial separation. Euclidean (straight line) distances often only an approximation and don't represent actual distance well (travel time on road networks, delays in connecting flights, more). 

- __Adjacency or neighborhood__: nominal or sometimes binary (0,1) equivalent of distance. Levels of adjacency exist: 1st, 2nd, 3rd  nearest neighbor and so on. 2nd order adjacency is connection to another unit through one other unit. 3rd order is connection to another unit through two other units.

- __Interaction__: the strength of the relationship between entities. This is typically an inverse function of distance but it doesn't have to be. It can also represent the physical ease with which people or wild animals move between units, the length of a shared boundary, or cooperation between governmental agencies (free trade for example), for example. 

## Load necessary R packages

There are several packages that we can use to work with areal data. We'll start with these six and show you a quick way to load them all at once.


```r
x <- c("spdep", "ggplot2", "ggmap", "rgdal", "rgeos", "maptools", "dplyr", "tidyr") # the packages
```


```r
install.packages(x) # warning: this may take a number of minutes
```


```r
lapply(x, library, character.only = TRUE) # load the required packages
```

```
## Loading required package: sp
## Loading required package: Matrix
## rgdal: version: 0.9-1, (SVN revision 518)
## Geospatial Data Abstraction Library extensions to R successfully loaded
## Loaded GDAL runtime: GDAL 1.11.1, released 2014/09/24
## Path to GDAL shared files: C:/Users/jcothren/Documents/R/win-library/3.1/rgdal/gdal
## GDAL does not use iconv for recoding strings.
## Loaded PROJ.4 runtime: Rel. 4.8.0, 6 March 2012, [PJ_VERSION: 480]
## Path to PROJ.4 shared files: C:/Users/jcothren/Documents/R/win-library/3.1/rgdal/proj
## rgeos version: 0.3-8, (SVN revision 460)
##  GEOS runtime version: 3.4.2-CAPI-1.8.2 r3921 
##  Polygon checking: TRUE 
## 
## Checking rgeos availability: TRUE
## 
## Attaching package: 'dplyr'
## 
## The following objects are masked from 'package:rgeos':
## 
##     intersect, setdiff, union
## 
## The following object is masked from 'package:stats':
## 
##     filter
## 
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
## 
## 
## Attaching package: 'tidyr'
## 
## The following object is masked from 'package:Matrix':
## 
##     expand
```

# Working with spatial data in R

Before we jump into defining neighborhoods, we need to spend a little time working with some spatial data in more common GIS formats. 

First, lets look at some data which reports the locations (in X and Y coordinates, unkown system at the moment) in London of cycle hires. 


```r
cycle<- read.csv("London_cycle_hire_locs.csv", header=T)
class(cycle)
```

```
## [1] "data.frame"
```

```r
# Inspect the column headings
head(cycle)
```

```
##                         Name      Village      X      Y Capacity
## 1 Kensington Olympia Station      Olympia 524384 179210       25
## 2            Ilchester Place   Kensington 524844 179509       24
## 3            Chepstow Villas Notting Hill 524897 180779       17
## 4           Turquoise Island Notting Hill 524940 181022       21
## 5         West Cromwell Road Earl's Court 525174 178737       24
## 6           Pembridge Villas Notting Hill 525179 180668       16
```

```r
# Summarize the columns
summary(cycle)
```

```
##                    Name            Village          X         
##  Abbey Orchard Street:  1   Marylebone : 17   Min.   :524384  
##  Abingdon Villas     :  1   Kensington : 11   1st Qu.:527956  
##  Alderney Street     :  1   Bloomsbury : 10   Median :529945  
##  Altab Ali Park      :  1   Holborn    : 10   Mean   :529833  
##  Ampton Street       :  1   Clerkenwell:  9   3rd Qu.:531887  
##  Appold Street       :  1   The Borough:  9   Max.   :534723  
##  (Other)             :286   (Other)    :226                   
##        Y             Capacity   
##  Min.   :177853   Min.   :13.0  
##  1st Qu.:179386   1st Qu.:16.0  
##  Median :180943   Median :18.0  
##  Mean   :180760   Mean   :21.8  
##  3rd Qu.:182010   3rd Qu.:26.0  
##  Max.   :183624   Max.   :54.0  
## 
```

A simple plot of the XY coordinates shows the spatial data.


```r
plot(cycle$X, cycle$Y)
```

![](Areal_Data_Analysis_files/figure-html/unnamed-chunk-5-1.png) 

However, just like in a GIS, this is not really a spatial data object in R (package `sp`). We can create a true spatial data object (and enable all the analytical functionality that goes with it) with the following command. By the way, most packages that work with spatial data in R have been written or updated to work some or all of the `sp` data objects.


```r
coordinates(cycle)<- c("X", "Y")
```

Now, look at the created object..


```r
class(cycle)
```

```
## [1] "SpatialPointsDataFrame"
## attr(,"package")
## [1] "sp"
```

```r
str(cycle)
```

```
## Formal class 'SpatialPointsDataFrame' [package "sp"] with 5 slots
##   ..@ data       :'data.frame':	292 obs. of  3 variables:
##   .. ..$ Name    : Factor w/ 292 levels "Abbey Orchard Street",..: 144 135 51 261 282 199 202 2 183 260 ...
##   .. ..$ Village : Factor w/ 64 levels "Aldgate","Angel",..: 37 26 36 36 14 36 26 26 36 14 ...
##   .. ..$ Capacity: int [1:292] 25 24 17 21 24 16 37 17 16 18 ...
##   ..@ coords.nrs : int [1:2] 3 4
##   ..@ coords     : num [1:292, 1:2] 524384 524844 524897 524940 525174 ...
##   .. ..- attr(*, "dimnames")=List of 2
##   .. .. ..$ : NULL
##   .. .. ..$ : chr [1:2] "X" "Y"
##   ..@ bbox       : num [1:2, 1:2] 524384 177853 534723 183624
##   .. ..- attr(*, "dimnames")=List of 2
##   .. .. ..$ : chr [1:2] "X" "Y"
##   .. .. ..$ : chr [1:2] "min" "max"
##   ..@ proj4string:Formal class 'CRS' [package "sp"] with 1 slot
##   .. .. ..@ projargs: chr NA
```

The class has become a `SpatialPointsDataFrame` which is a type of S4 object that requires handling slightly differently. The `str()` output contains lots of `@` symbols which denote a different slot (or collection of data types). Typing `cycle@data` will extract the attribute data. The X and Y locational information can now be found in the `coords` slot, while the `bbox` slot contains the bounding box coordinates and the `pro4string` slot contains the projection, or CRS (coordinate reference system) information. 


```r
# the attribute slot
head(cycle@data)
```

```
##                         Name      Village Capacity
## 1 Kensington Olympia Station      Olympia       25
## 2            Ilchester Place   Kensington       24
## 3            Chepstow Villas Notting Hill       17
## 4           Turquoise Island Notting Hill       21
## 5         West Cromwell Road Earl's Court       24
## 6           Pembridge Villas Notting Hill       16
```

```r
# the coordinate slot
head(cycle@coords)
```

```
##           X      Y
## [1,] 524384 179210
## [2,] 524844 179509
## [3,] 524897 180779
## [4,] 524940 181022
## [5,] 525174 178737
## [6,] 525179 180668
```

```r
# the bounding box slot
cycle@bbox
```

```
##      min    max
## X 524384 534723
## Y 177853 183624
```

```r
# the projection string (empty at this point)
cycle@proj4string
```

```
## CRS arguments: NA
```

We have not specified yet spectified the data projection that slot is empty at the moment. We therefore need to refer to the correct Proj4 string information. These are loaded with the rgdal package and can simply be referred to with an ID. To see the available Coordiante Reference Systems (CRS's) you use the following code. Note that this all comes from the Proj4 libary used in most open-source GIS packages. ESRI uses the same codes (for the most part).


```r
EPSG <- make_EPSG()
head(EPSG)
```

```
##   code                                               note
## 1 3819                                           # HD1909
## 2 3821                                            # TWD67
## 3 3824                                            # TWD97
## 4 3889                                             # IGRS
## 5 3906                                         # MGI 1901
## 6 4001 # Unknown datum based upon the Airy 1830 ellipsoid
##                                                                                            prj4
## 1 +proj=longlat +ellps=bessel +towgs84=595.48,121.69,515.35,4.115,-2.9383,0.853,-3.408 +no_defs
## 2                                                         +proj=longlat +ellps=aust_SA +no_defs
## 3                                    +proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs
## 4                                    +proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs
## 5                            +proj=longlat +ellps=bessel +towgs84=682,-203,480,0,0,0,0 +no_defs
## 6                                                            +proj=longlat +ellps=airy +no_defs
```

Our cycle data is in British National Grid. We can search for this within the EPSG object as follows...


```r
with(EPSG, EPSG[grep("British National", note),])
```

```
##       code                                note
## 3424 27700 # OSGB 1936 / British National Grid
##                                                                                                          prj4
## 3424 +proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +datum=OSGB36 +units=m +no_defs
```

and see that the code we are after is 27700. You can try this other systems. For example, to find the UTM Zone 15 code, you can type...


```r
# search for the string "zone=15" in the prj4 field
utm_codes <- with(EPSG, EPSG[grep("+zone=15", prj4),])
utm_codes[,c("code","note")]
```

```
##       code                                   note
## 490   2027             # NAD27(76) / UTM zone 15N
## 1613  3159           # NAD83(CSRS) / UTM zone 15N
## 2176  3722       # NAD83(NSRS2007) / UTM zone 15N
## 2199  3745           # NAD83(HARN) / UTM zone 15N
## 2416  4488 # Mexican Datum of 1993 / UTM zone 15N
## 3131 26715                 # NAD27 / UTM zone 15N
## 3279 26915                 # NAD83 / UTM zone 15N
## 3658 31969           # SIRGAS 2000 / UTM zone 15N
## 3855 32215                # WGS 72 / UTM zone 15N
## 3915 32315                # WGS 72 / UTM zone 15S
## 3975 32415              # WGS 72BE / UTM zone 15N
## 4035 32515              # WGS 72BE / UTM zone 15S
## 4096 32615                # WGS 84 / UTM zone 15N
## 4164 32715                # WGS 84 / UTM zone 15S
```

There are 14, one for each datum and northern and southern hemisphere. We'll do this again later with the Arkansas State Plane North system. 

We can assign the British National Grid CRS to our cycle data easily enough.


```r
# Save the string to BNG
BNG <- CRS("+init=epsg:27700")

# assign string to cycle spatial data object
proj4string(cycle) <- BNG

# look at the result to be sure
cycle@proj4string
```

```
## CRS arguments:
##  +init=epsg:27700 +proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717
## +x_0=400000 +y_0=-100000 +datum=OSGB36 +units=m +no_defs
## +ellps=airy
## +towgs84=446.448,-125.157,542.060,0.1502,0.2470,0.8421,-20.4894
```

From this point we can combine the data with other spatial information and also perform transformations on the data. This is especially useful when exporting to other software. Shapefiles are extremely simple to import/export to/from an R session and are handled as spatial objects in the same way as above. In this case we are going to load in a SpatialPolygonsDataframe. We can specify the CRS when the data at this stage (it is _BNG_ as above).


```r
sport <- readShapePoly("london_sport.shp", proj4string= BNG)
class(sport)
```

```
## [1] "SpatialPolygonsDataFrame"
## attr(,"package")
## [1] "sp"
```

Note the shapefile was imported as a `SpatialPolygonsDataFrame`. Take a look at the attribute table headings (these are the values stored in the `data` slot).


```r
names(sport)
```

```
## [1] "ons_label"  "name"       "Partic_Per" "Pop_2001"
```

Also verify the correct CRS, 


```r
# should be the same as the cycle data
sport@proj4string
```

```
## CRS arguments:
##  +init=epsg:27700 +proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717
## +x_0=400000 +y_0=-100000 +datum=OSGB36 +units=m +no_defs
## +ellps=airy
## +towgs84=446.448,-125.157,542.060,0.1502,0.2470,0.8421,-20.4894
```

Let's plot both objects to see what we have spatially (`plot` is overloaded for spatial objects).


```r
plot(sport, col="lightgrey")
plot(cycle, add=T, col= "white", pch=21)
```

![](Areal_Data_Analysis_files/figure-html/unnamed-chunk-16-1.png) 

Refer to http://spatial.ly/2013/12/introduction-spatial-data-ggplot2/ to learn more about plotting and analysis of these data.

To export the data as a shapefile use the following syntax for the point data:

```r
# writePointsShape(cycle, "cycle.shp")
# writePolyShape(sport, "london_sport.shp")
```

Much of this content was adapted from a worksheet written by Dr Gavin Simpson (UCL Geography).

# Defining neighborhoods and weights

Before we can discuss methods to deal with the issues outlined above, we must establish a framework in which to think about the spatial concepts outlined above. Up to know, we have thought of spatial relationships in terms of distance - distance point events, distance of events from a covariate, etc. When we have aggregated _marked events_ to areal units, our concept of spatial relationship has to include other forms of relationship besides distance. Measuring distances between, say, centroids of polygons is not always (is rarely) a good measure of spatial proximity due to variable sized units and the unknown distribution of events within the unit. Instead, we tend to think in terms of _neighborhoods_ of _connected_ units. Connectivity aomong areal units is typically defined as a shared polygon point or points (i.e. a shared vertex or edge). However, connectivity could also be defined shared environmental, political or infrastructure entities (e.g. respectively, a stream reach or animal habitat, a trade agreeement, or a railroad). Connectivity must be based on the goals of the study itself and as such its definition represents significant _a priori_ knowledge. It is a critical component of model design and its impact on the outcome should not be underestimated. 

After a _neigbhorhood_ is defined, the relative strength of the connected units within it may be considered. For example, longer shared boundaries between units may warrant a stronger connection. More or larger highways between units may warrent stronger relative connections. The strength of the connections is considered a _weight_ and stored in a _weight matrix_. The following section illustrates these concepts with a simple example. The section after that formalizes the concepts and introduces the data structures and commands used by the `spded` package.  

## Simple example

Consider the following tesselation of a region (selected tracts from Washington County). ![Washington County Tracts](C:\Users\jcothren\Dropbox\teaching\GEOS4863\Areal_Data_Notes\Simple_Example.png)

We have looked at that application of distance matrices in point pattern analysis. In areal data analysis a more common and useful concept is that of a neigbhorhood. A neighborhood consists of a polygon (an areal unit) and all polygons that are connected to it in a way defined by the researcher and, hopefully, based on some meaningful geometric attributes. We'll start simply and assume a neighborhood is composed of a tract and all tracts that share a border (i.e. at least two sequential vertices) with it. Neighborhoods can be represented mathematically as an _adjacency matrix_.


```r
# create a list of polygon names (from the image)
poly.names <- c('A','B','C','D','E','F','G','H')

# manually create the adjacency matrix
adjacency.matrix <- matrix(c(0,0,0,0,1,0,0,0,  0,0,1,0,0,0,0,0,  0,1,0,1,0,0,1,1, 0,0,1,0,1,1,1,0, 1,0,0,1,0,1,0,0, 0,0,0,1,1,0,1,1,  0,0,1,1,0,1,0,1, 0,0,1,0,0,1,1,0),nrow=8,ncol=8,byrow=TRUE,dimnames = list(poly.names,poly.names)) 

adjacency.matrix
```

```
##   A B C D E F G H
## A 0 0 0 0 1 0 0 0
## B 0 0 1 0 0 0 0 0
## C 0 1 0 1 0 0 1 1
## D 0 0 1 0 1 1 1 0
## E 1 0 0 1 0 1 0 0
## F 0 0 0 1 1 0 1 1
## G 0 0 1 1 0 1 0 1
## H 0 0 1 0 0 1 1 0
```

The adjacency matrix is a compact structure which tells us the connectivity between areal units. It has a one-to-one relationship with an undirected graph.  

Second order adjacency can be easily computed with this data structure...


```r
second.order <- adjacency.matrix %*% adjacency.matrix
second.order
```

```
##   A B C D E F G H
## A 1 0 0 1 0 1 0 0
## B 0 1 0 1 0 0 1 1
## C 0 0 4 1 1 3 2 1
## D 1 1 1 4 1 2 2 3
## E 0 0 1 1 3 1 2 1
## F 1 0 3 2 1 4 2 1
## G 0 1 2 2 2 2 4 2
## H 0 1 1 3 1 1 2 3
```

Here, the entries tell you how many second order connections there are between pairs of polygons. While this is no longer binary, it does tell you that, for example, __A__ has only one second order connection to __F__ (through __E__). Likewise, __C__ has 3 second order connections to __F__ (through __G__, __H__, and __D__).


```r
third.order <- second.order %*% adjacency.matrix
third.order
```

```
##   A B  C  D E  F  G H
## A 0 0  1  1 3  1  2 1
## B 0 0  4  1 1  3  2 1
## C 1 4  4 10 4  5  9 9
## D 1 1 10  6 7 10 10 5
## E 3 1  4  7 2  7  4 4
## F 1 3  5 10 7  6 10 9
## G 2 2  9 10 4 10  8 8
## H 1 1  9  5 4  9  8 4
```

Now, because this dataset is so small, almost every polygon can be connected to another in "three steps". Only __A__ and __B__ remain unconnected.

Often, you will work with _row normalized_ adjacency matrices.


```r
#row normalized
adjacency.matrix.rn <- matrix(c(0,0,0,0,0,1,0,0,  0,0,1,0,0,0,0,0,  0,1/3,0,1/3,0,0,1/3,0, 0,0,1/4,0,1/4,1/4,1/4,0, 1/3,0,0,1/3,0,1/3,0,0, 0,0,0,1/4,1/4,0,1/4,1/4,  0,0,1/4,1/4,0,1/4,0,1/4, 0,0,1/3,0,0,1/3,1/3,0),nrow=8,ncol=8,byrow=TRUE,dimnames = list(poly.names,poly.names))

adjacency.matrix.rn
```

```
##           A         B         C         D    E         F         G    H
## A 0.0000000 0.0000000 0.0000000 0.0000000 0.00 1.0000000 0.0000000 0.00
## B 0.0000000 0.0000000 1.0000000 0.0000000 0.00 0.0000000 0.0000000 0.00
## C 0.0000000 0.3333333 0.0000000 0.3333333 0.00 0.0000000 0.3333333 0.00
## D 0.0000000 0.0000000 0.2500000 0.0000000 0.25 0.2500000 0.2500000 0.00
## E 0.3333333 0.0000000 0.0000000 0.3333333 0.00 0.3333333 0.0000000 0.00
## F 0.0000000 0.0000000 0.0000000 0.2500000 0.25 0.0000000 0.2500000 0.25
## G 0.0000000 0.0000000 0.2500000 0.2500000 0.00 0.2500000 0.0000000 0.25
## H 0.0000000 0.0000000 0.3333333 0.0000000 0.00 0.3333333 0.3333333 0.00
```

We'll examine the effect of this later. So the _adjacency matrix_ represents our neighborhood. We can create an `spdep` `weights list` object from the matrix using the `mat2listw` command. 


```r
# skip nb object, go straight to listw from matrix
W <- mat2listw(adjacency.matrix)
W
```

```
## Characteristics of weights list object:
## Neighbour list object:
## Number of regions: 8 
## Number of nonzero links: 24 
## Percentage nonzero weights: 37.5 
## Average number of links: 3 
## 
## Weights style: M 
## Weights constants summary:
##   n nn S0 S1  S2
## M 8 64 24 48 336
```

```r
# but a neighborhood object is part of the listw
W$neighbours
```

```
## Neighbour list object:
## Number of regions: 8 
## Number of nonzero links: 24 
## Percentage nonzero weights: 37.5 
## Average number of links: 3
```

```r
W.rn <- mat2listw(adjacency.matrix.rn)
```

Note that this object contains a number of important attributes which describe the overall connectivity of the units (all the neighborhoods, not just one of them). In this case we see that the 8 regions (units) have 24 connections out of a possible 64 or 37.5%. The average region has three connections. By querying the `W$neighbours` field (this is an `nb` object we'll see later) we can see our adjancency matrix stored differently.


```r
W$neighbours[2]  # unit 2 (B) connections
```

```
## [[1]]
## C 
## 3
```

```r
W$neighbours[3]  # unit 3 (C) connections
```

```
## [[1]]
## B D G H 
## 2 4 7 8
```

The corresponding weights are stored in `W$weights`. In this case, they are all 1's.

![Simple Example with Attributes](C:\Users\jcothren\Dropbox\teaching\GEOS4863\Areal_Data_Notes\Simple_Example_Attributes.png)


```r
# add attributes for polygons A-H
y = c(27,51,45,47,54,18,7,52)
```

## Creating neighbors in `spdep`

The R package `spdep` uses a neighborhood object, or `nb`, rather than an adjacency matrix. In this section we'll work through several different ways of defining neighbors. In the previous section we defined an adjacency matrix manually and then used it to create a `nb`. This helps you understand what's happening a little under the hood but it's also good to know in case you have odd connectivity criteria. Typically the creation of the object is handled internally and is based on a number of well-known criteria.

### Contiguity Neighbors

Polygon contiguity defined by shared boundaries or vertices is the most common connectivity method used. It is straighforward but you have to be careful, though, to make sure that your polygon dataset is topogologicaly intact. If there are overlaps among polygons or slivers between edges (i.e, the dataset must be _planar_), then your `nb` will not be correct. We'll not worry about that here. There are no R packages which allow you to this but QGIS and ArcGIS all have topology enforcing functions. This violates our "reproducibility" requirement in terms of scientific software but in this case it is more of a data issue than an analysis issue. Lattice structure are just assumed to be planar. However, the `poly2nb` function in `spdep` has a `snap` argument you can set to compensate for a non-planar dataset. The distance you specify for this argument is the distance at which two vertices are considered equivalent. 

In `polynb`, the default connectivity condition is that two polygons share at least one vertex. This is the so-called _queen-contiguity_ condition. A more restrictive condition is _rook-continguity_ in which neighbors must share an edege. Let's create queen and rook contiguity `nb` objects and examine them using tract polygons from a Syracuse NY. 


```r
# load and plot the polygon layer
syracuse <- readOGR(dsn=".", layer="Syracuse")  #note that there is no projection information (bad!)
```

```
## OGR data source with driver: ESRI Shapefile 
## Source: ".", layer: "Syracuse"
## with 63 features and 17 fields
## Feature type: wkbPolygon with 2 dimensions
```

```r
plot(syracuse,col='grey', main="Queen contiguity in blue, Rook in white")

# create continuity nb
cuse_nb_queen <- poly2nb(syracuse, queen=TRUE, row.names <- syracuse$AREAKEY)
cuse_nb_rook <- poly2nb(syracuse, queen=FALSE, row.names <- syracuse$AREAKEY) 

# compare using the plot command and the nb objects
coords <- coordinates(syracuse) # extract the centroids of each tract
plot(cuse_nb_queen, coords, col='blue', add=TRUE, lwd=5)
plot(cuse_nb_rook, coords, col='white', add=TRUE)
```

![](Areal_Data_Analysis_files/figure-html/unnamed-chunk-25-1.png) 

Note that rook contiguity is a strict sub-set of queen contiguity.


```r
summary(cuse_nb_queen)
```

```
## Neighbour list object:
## Number of regions: 63 
## Number of nonzero links: 346 
## Percentage nonzero weights: 8.717561 
## Average number of links: 5.492063 
## Link number distribution:
## 
##  1  2  3  4  5  6  7  8  9 
##  1  1  5  9 14 17  9  6  1 
## 1 least connected region:
## 36067005602 with 1 link
## 1 most connected region:
## 36067002900 with 9 links
```

```r
summary(cuse_nb_rook)
```

```
## Neighbour list object:
## Number of regions: 63 
## Number of nonzero links: 308 
## Percentage nonzero weights: 7.760141 
## Average number of links: 4.888889 
## Link number distribution:
## 
##  1  2  3  4  5  6  7  8 
##  1  1  7 18 15 11  9  1 
## 1 least connected region:
## 36067005602 with 1 link
## 1 most connected region:
## 36067005500 with 8 links
```

### Graph-based Neighbors

In the case of contingity neighbors we used the entire polygon to define connectivity. It is also possible to use the polygon centroid (or another point, perhaps a weighted centroid) as a representative point. Once representative points are available, then other measures can be used to define neighborhoods including graph-based neighbors, distance thresholds, and k-nearest neighbors.

In general, we use a class called a proximity graph. There is a simply a graph in which two vertices are connected
by an edge if and only if the vertices satisfy particular geometric requirements. The most direct graph representation of neighbors is to make a Delaunay triangulation of the points. The neighbor relationships are defined by the triangulation, which extends outwards to the convex hull of the points and which is planar. Note that graph-based representations construct the interpoint relationships based on Euclidean distance, with no option to use Great Circle distances
for geographical coordinates. Because it joins distant points around the convex hull, it may be worthwhile to thin the triangulation as a Sphere of Influence (SOI) graph, removing links that are relatively long. Points are SOI neighbours if circles centred on the points, of radius equal to the points’ nearest neighbour distances, intersect in two places. The following code computes 1) the Delaunay triangulation of the centroids; 2) Sphere of Influence (SOI) neighbors; 3) Gabriel graph neighbors; and 4) relative graph neighbors. 


```r
IDs <- row.names(as(syracuse, "data.frame"))
coords <- coordinates(syracuse)
cuse_delaunay <- tri2nb(coords, row.names = IDs)
```

```
## 
##      PLEASE NOTE:  The components "delsgs" and "summary" of the
##  object returned by deldir() are now DATA FRAMES rather than
##  matrices (as they were prior to release 0.0-18).
##  See help("deldir").
##  
##      PLEASE NOTE: The process that deldir() uses for determining
##  duplicated points has changed from that used in version
##  0.0-9 of this package (and previously). See help("deldir").
```

```r
cuse_soi <- graph2nb(soi.graph(cuse_delaunay, coords), row.names = IDs)
cuse_gabriel <- graph2nb(gabrielneigh(coords), row.names = IDs)
cuse_relative <- graph2nb(relativeneigh(coords), row.names = IDs)
```


```r
par(mfrow = c(2,2))

plot(syracuse,col='grey', main="delaunay")
plot(cuse_delaunay, coords, col='blue', add=TRUE)

plot(syracuse,col='grey', main="soi") #removing long edges
plot(cuse_soi, coords, col='blue', add=TRUE)

plot(syracuse,col='grey', main="gabriel") #removing long edges
plot(cuse_gabriel, coords, col='blue', add=TRUE)

plot(syracuse,col='grey', main="relative") #removing long edges
plot(cuse_relative, coords, col='blue', add=TRUE)
```

![](Areal_Data_Analysis_files/figure-html/unnamed-chunk-28-1.png) 

```r
par(mfrow=c(1,1))
```

Delaunay triangulation neighbors and SOI neighbors are symmetric by design – if $i$ is a neighbour of $j$, then $j$ is a neighbour of $i$. The Gabriel graph is also a subgraph of the Delaunay triangulation, retaining a different set of neighbours (Matula and Sokal, 1980). It does not, however, guarantee symmetry; the same applies to Relative graph neighbours (Toussaint, 1980). The `graph2nb` function takes a `sym` argument to insert links to restore symmetry, but the graphs then no longer exactly fulfill their neighbor criteria. All the graph-based neighbour schemes always ensure that all the points will have at least one neighbour. Subgraphs of the full triangulation may also have more than one graph after trimming. The function `is.symmetric.nb` can be used to check for symmetry, with argument `force=TRUE` if the symmetry attribute is to be overridden, and `n.comp.nb` reports the number of graph components and the components to which points belong (after enforcing symmetry, because the algorithm assumes that the graph is not directed). When there are more than one graph component, the matrix representation of the spatial weights can become block-diagonal if observations are appropriately sorted. This is especially important in large datasets.

### Distance-based neighbors

Another method is to choose the $k$ nearest points as neighbors. This method is adaptive in that distance between neighbors varies across the study area, taking account of differences in the densities of areal entities. Naturally, in the majority of cases, it leads to asymmetric neighbours, but it will ensure that all areas have $k$ neighbours. The `knearneigh` function returns an intermediate form converted to an `nb` object by `knn2nb`. Note that `knearneigh` can also take a `longlat` argument to handle geographical coordinates and uses short geodesics as distances rather than euclidean distance.


```r
cuse_knn1 <- knn2nb(knearneigh(coords, k = 1), row.names = IDs) # closest neighbor
cuse_knn2 <- knn2nb(knearneigh(coords, k = 2), row.names = IDs) # two closest neighbors
cuse_knn4 <- knn2nb(knearneigh(coords, k = 4), row.names = IDs) # four closest neighbors

par(mfrow = c(2,2))

plot(syracuse,col='grey', main="nearest neighbor")
plot(cuse_knn1, coords, col='blue', add=TRUE)

plot(syracuse,col='grey', main="two nearest neighbors") 
plot(cuse_knn2, coords, col='blue', add=TRUE)

plot(syracuse,col='grey', main="four nearest neighbors") 
plot(cuse_knn4, coords, col='blue', add=TRUE)

par(mfrow=c(1,1))
```

![](Areal_Data_Analysis_files/figure-html/unnamed-chunk-29-1.png) 

The figures above show the neighbour relationships for $k = 1,2,4$, with many components for $k = 1$. If necessary, $k$-nearest neighbour objects can be made symmetrical using the `make.sym.nb` function. The large number of disjoint connected subgraphs in `cuse_knn1` while the other two graphs are completely connected.

And we can check for symmetry and connectivity again.


```r
nb_l <- list(k1 = cuse_knn1, k2 = cuse_knn2, k4 = cuse_knn4)
sapply(nb_l, function(x) is.symmetric.nb(x, verbose = FALSE, force = TRUE))
```

```
##    k1    k2    k4 
## FALSE FALSE FALSE
```

```r
sapply(nb_l, function(x) n.comp.nb(x)$nc)
```

```
## k1 k2 k4 
## 15  1  1
```

The $k = 1$ object is also useful in finding the minimum distance at which all areas have a distance-based neighbour. Using the `nbdists` function, we can calculate a list of vectors of distances corresponding to the neighbor object, here for first nearest neighbors. The greatest value will be the minimum distance needed to make sure that all the areas are linked to at least one neighbour. The `dnearneigh` function is used to find neighbors with an interpoint distance, with arguments `d1` and `d2` setting the lower and upper distance bounds; it can also take a longlat argument to handle geographical coordinates.


```r
cuse_nn_dists <- unlist(nbdists(cuse_knn1, coords)) # a little tedious to get the distances from this list
summary(cuse_nn_dists)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   395.7   587.3   700.1   760.4   906.1  1545.0
```

Now we can find the maximum distance between two neighbors and use that distance to compute neighbhors.


```r
cuse_max_1nn <- max(cuse_nn_dists) # which we saw in the summary above

cuse_075_nb <- dnearneigh(coords, d1 = 0, d2 = 0.75 * cuse_max_1nn, row.names = IDs)
cuse_100_nb <- dnearneigh(coords, d1 = 0, d2 = 1 * cuse_max_1nn, row.names = IDs)
cuse_150_nb <- dnearneigh(coords, d1 = 0, d2 = 1.5 * cuse_max_1nn, row.names = IDs)

par(mfrow = c(2,2))

plot(syracuse,col='grey', main="0 to 0.75*max distance")
plot(cuse_075_nb, coords, col='blue', add=TRUE)

plot(syracuse,col='grey', main="0 to max distance") 
plot(cuse_100_nb, coords, col='blue', add=TRUE)

plot(syracuse,col='grey', main="0 to 1.5 * max distance") 
plot(cuse_150_nb, coords, col='blue', add=TRUE)

par(mfrow=c(1,1))
```

![](Areal_Data_Analysis_files/figure-html/unnamed-chunk-32-1.png) 

And check for symmetric and disjoint subgraphs.


```r
nb_l <- list(d1 = cuse_075_nb, d2 = cuse_100_nb , d3 = cuse_150_nb)
sapply(nb_l, function(x) is.symmetric.nb(x, verbose = FALSE, force = TRUE))  # symmetric?
```

```
##   d1   d2   d3 
## TRUE TRUE TRUE
```

```r
sapply(nb_l, function(x) n.comp.nb(x)$nc) # number of disjoint subgraphs
```

```
## d1 d2 d3 
##  4  1  1
```

The figure above shows how the numbers of distance-based neighbors increase with moderate increases in distance. Moving from 0:75 times the minimum all-included distance (1158 m), to the all-included distance (1545 m), and 1.5 times the minimum all-included distance (2317 m), the numbers of links grow rapidly. This is a major problem when some of the first nearest neighbor distances in a study area are much larger than others, since to avoid no-neighbour areal entities, the distance criterion will need to be set such that many areas have many neighbours. In Syracuse, the census tracts are of similar areas, but were we to try to use the distance-based neighbour criterion on the eight-county study area, the smallest distance securing at least one neighbour for every areal entity is over 38 km.


```r
# load the entire 8-county study area
NY <- readOGR(dsn=".", layer="NewYorkTracts")
```

```
## OGR data source with driver: ESRI Shapefile 
## Source: ".", layer: "NewYorkTracts"
## with 281 features and 17 fields
## Feature type: wkbPolygon with 2 dimensions
```

```r
NY_nb <- poly2nb(NY, queen=FALSE)

plot(NY, col='lightgrey')
plot(NY_nb, coordinates(NY), add=TRUE, col='blue')
```

![](Areal_Data_Analysis_files/figure-html/unnamed-chunk-34-1.png) 

```r
ny_dists <- unlist(nbdists(NY_nb, coordinates(NY)))
summary(ny_dists)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    82.7  1543.0  3449.0  5871.0  9147.0 38440.0
```

If the areal entities are approximately regularly spaced, using distance-based neighbours is not necessarily a problem. Provided that care is taken to handle the side effects of “weighting” areas out of the analysis, using lists of neighbours with no-neighbour areas is not necessarily a problem either, but certainly ought to raise questions. Different disciplines handle the definition of neighbours in their own ways by convention; in particular, it seems that ecologists frequently use distance bands. If many distance bands are used, then the results begin to approach the variogram, although the underlying understanding of spatial autocorrelation seems to be by contagion rather than continuous. 

### Higher-order neighbors

Distance bands can be generated by using a sequence of `d1` and `d2` argument values for the `dnearneigh` function. In this way we can construct a _spatial autocorrelogram_ as understood in ecology. In other conventions, correlograms are constructed by taking an input list of neighbours as the first-order sets, and stepping out across the graph to second-, third-, and higher-order neighbours (as we saw using the adjacency matrix above) based on the number of links traversed, but not permitting cycles, which could risk making $i$ a neighbour of $i$ itself. The `nblag` function takes an existing neighbour list and returns a list of lists, from first to `maxlag` order neighbours. It's easier to show this than explain.


```r
# we'll need to work with adjacency matrices for this so require the igraph package
require(igraph)
```

```
## Loading required package: igraph
```

```r
# create the list of nb lists
cuse_nb_lags <- nblag(cuse_nb_queen, maxlag = 9)

# create a table for easier viewing
Table <- matrix(data=NA, nrow=63, ncol=9)
```

Table 1 shows how the wave of connectedness in the graph spreads to the third order, receding to the eighth order, and dying away at the ninth order – there are no tracts nine steps from each other in this graph. Both the distance bands and the graph step order approaches to spreading neighbourhoods can be used to examine the shape of relationship intensities in space, like the variogram, and can be used in attempting to look at the effects of scale.

### Grid neighbors

When the data are known to be arranged in a regular, rectangular grid, the `cell2nb` function can be used to construct neighbour lists, including those on a torus. These are useful for simulations, because, since all areal entities have equal numbers of neighbours, and there are no edges, the structure of the graph is as neutral as can be achieved. Neighbours can either be of type rook or queen (in image processing and GIS this is typically called 4-connected and 8-connected, respectively).


```r
cell2nb(7, 7, type = "rook", torus = TRUE)
```

```
## Neighbour list object:
## Number of regions: 49 
## Number of nonzero links: 196 
## Percentage nonzero weights: 8.163265 
## Average number of links: 4
```

When a regular, rectangular grid is not complete, then we can use knowledge of the cell size stored in the grid topology to create an appropriate list of neighbours, using a tightly bounded distance criterion. Neighbour lists of this kind are commonly found in ecological assays, such as studies of species richness at a national or continental scale. It is also in these settings, with moderately large $n$, here $n = 3,103$, that the use of a sparse, list based representation shows its strength. Handling a 281x281 matrix for the eight-county census tracts is feasible, easy for a 63x63 matrix for Syracuse census tracts, but demanding for a 3103 x 3103 matrix.


```r
data(meuse.grid)
head(meuse.grid) # note the x and y columns with coordinate
```

```
##        x      y part.a part.b      dist soil ffreq
## 1 181180 333740      1      0 0.0000000    1     1
## 2 181140 333700      1      0 0.0000000    1     1
## 3 181180 333700      1      0 0.0122243    1     1
## 4 181220 333700      1      0 0.0434678    1     1
## 5 181100 333660      1      0 0.0000000    1     1
## 6 181140 333660      1      0 0.0122243    1     1
```

```r
coordinates(meuse.grid) <- c("x", "y") # use the x and y columns to convert data.frame into SpatialPointsDataFrame
gridded(meuse.grid) <- TRUE # confirm it is a grid type
summary(meuse.grid)
```

```
## Object of class SpatialPixelsDataFrame
## Coordinates:
##      min    max
## x 178440 181560
## y 329600 333760
## Is projected: NA 
## proj4string : [NA]
## Number of points: 3103
## Grid attributes:
##   cellcentre.offset cellsize cells.dim
## x            178460       40        78
## y            329620       40       104
## Data attributes:
##      part.a           part.b            dist        soil     ffreq   
##  Min.   :0.0000   Min.   :0.0000   Min.   :0.0000   1:1665   1: 779  
##  1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:0.1193   2:1084   2:1335  
##  Median :0.0000   Median :1.0000   Median :0.2715   3: 354   3: 989  
##  Mean   :0.3986   Mean   :0.6014   Mean   :0.2971                    
##  3rd Qu.:1.0000   3rd Qu.:1.0000   3rd Qu.:0.4402                    
##  Max.   :1.0000   Max.   :1.0000   Max.   :0.9926
```

```r
dst <- max(slot(slot(meuse.grid, "grid"), "cellsize")) # find the distance between grid points 

# 4-connnected (rook)
mg_nb <- dnearneigh(coordinates(meuse.grid), 0, dst) # create neigbhorhood contiguity by distance (from 0 to dst)
mg_nb
```

```
## Neighbour list object:
## Number of regions: 3103 
## Number of nonzero links: 12022 
## Percentage nonzero weights: 0.1248571 
## Average number of links: 3.874315
```

```r
table(card(mg_nb))
```

```
## 
##    1    2    3    4 
##    1  133  121 2848
```

```r
# 8-connected (queen)
mg_nb_8 <- dnearneigh(coordinates(meuse.grid), 0, dst*sqrt(2)) # create neigbhorhood contiguity by distance (from 0 to dist*sqrt(2))
mg_nb_8
```

```
## Neighbour list object:
## Number of regions: 3103 
## Number of nonzero links: 23920 
## Percentage nonzero weights: 0.2484263 
## Average number of links: 7.708669
```

```r
table(card(mg_nb_8))
```

```
## 
##    3    4    5    6    7    8 
##    2   82   93   79  129 2718
```

There is a simpler way to do this but in the case of 4 and 8 connectivity. However, the distance method is more general and allows you to create distance rings as before. This is an important concept we'll encounter later. To see the simpler version, load this reasonably large Land Use Land Cover (LULC) dataset from Arkansas.


```r
# read lulc layer and convert numbers to factors (lulc classes found in the associated xml file)
lulc <- readGDAL('IMAGE_DBO_LULC_FALL_CAST2006.tif')
```

```
## IMAGE_DBO_LULC_FALL_CAST2006.tif has GDAL driver GTiff 
## and has 2824 rows and 4550 columns
```

```r
summary(lulc)
```

```
## Object of class SpatialGridDataFrame
## Coordinates:
##         min       max
## x  513085.5  642760.5
## y 3893613.0 3974097.0
## Is projected: TRUE 
## proj4string :
## [+proj=utm +zone=15 +datum=NAD83 +units=m +no_defs +ellps=GRS80
## +towgs84=0,0,0]
## Grid attributes:
##   cellcentre.offset cellsize cells.dim
## x          513099.7     28.5      4550
## y         3893627.2     28.5      2824
## Data attributes:
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    11.0   100.0   100.0   119.5   100.0   210.0
```

Note that this dataset is relatively small by todays standards and while R doesn't have trouble working with it, larger datasets can pose a problem. We need to convert the numbers to factors (classes).


```r
lulc$band1 <- as.factor(lulc$band1)
summary(lulc)
```

```
## Object of class SpatialGridDataFrame
## Coordinates:
##         min       max
## x  513085.5  642760.5
## y 3893613.0 3974097.0
## Is projected: TRUE 
## proj4string :
## [+proj=utm +zone=15 +datum=NAD83 +units=m +no_defs +ellps=GRS80
## +towgs84=0,0,0]
## Grid attributes:
##   cellcentre.offset cellsize cells.dim
## x          513099.7     28.5      4550
## y         3893627.2     28.5      2824
## Data attributes:
##      11      13      31      41      51     100     201     202     203 
##  367095   88271   14017  256361  752300 8210548  194657   24436   15240 
##     205     206     208     209     210 
##   36466    5784   66703 1523360 1293962
```

Now the summary tells us how many cells of each class are in the data set. We can easily create a 4-connected `nb` object.


```r
# create weight matrix, rook, and note that we don't need the dataset itself (just it's size)
# lulc_nb <- cell2nb(nrow(lulc),ncol(lulc), type='rook')
# W <- nb2listw(w.grid)
```

This is not an efficient operation but we'll use the `lulc_nb` later when we discuss join-count methods.

## Defining spatial weights

Once you've defined a neighborhood you can begin to think about spatial weights as they apply to understanding spatial correlation and modeling. Spatial weights can be seen as a list of weights indexed by a list of neighbors, where the weight of the link between $i$ and $j$ is the $k$ th element of the $i$ th weights list component, and $k$ tells youwhich of the $i$ th neighbor list component values is equal to $j$. If you look at the simple example we started with this will make sense...


```r
# list of neighbors of C
W$neighbours[3]
```

```
## [[1]]
## B D G H 
## 2 4 7 8
```

```r
# and their weights
W$weights[3]
```

```
## [[1]]
## [1] 1 1 1 1
```

It is simple to create a weights object, `listw`, from an `nb` object. Let's work with the 8-county dataset from New York.


```r
NY_listw <- nb2listw(NY_nb, style = "B")

# tenth tract in the list neighbors
NY_listw$neighbours[10]
```

```
## [[1]]
## [1]  5  6  9 11 12
```

```r
# and their weights
NY_listw$weights[10]
```

```
## [[1]]
## [1] 1 1 1 1 1
```

The `nb2listw` command creates the new object with binary, `style = "B"`, weights. We have other options for `style`: `W` is row standarized, `C` is globaly standardized (all weights sum to the number of units), `S` is variance-stablising (see `?nb2listw` for the reference that explains), and 'U' is equal to `C` divided by the number of neighbors (sums over all links is 1). The weights for `S` vary less than for `style="W"`  


```r
# row normalized weights
NY_listw <- nb2listw(NY_nb, style = "W")

# tenth tract in the list neighbors
NY_listw$neighbours[10]
```

```
## [[1]]
## [1]  5  6  9 11 12
```

```r
NY_listw$neighbours[15]
```

```
## [[1]]
## [1]  1 14 16 50
```

```r
# and their weights
NY_listw$weights[10]
```

```
## [[1]]
## [1] 0.2 0.2 0.2 0.2 0.2
```

```r
NY_listw$weights[15]
```

```
## [[1]]
## [1] 0.25 0.25 0.25 0.25
```

If you think about row-normalization for a moment you can see how it might be useful. While the `B` form sums the attributes of neighbors, the `W` form averages the attributes of neighbors. Because units at the edge of the study area tend to have fewer neighbors, their weights tend to be exagerated compared to highly connected units. Note also that when row-normalization is chosen, the weight matrix is no longer symmetric. 

We can use the `unlist` function again to get some statistics on weights.


```r
NY_listw_B <- nb2listw(NY_nb, style = "B")
NY_listw_C <- nb2listw(NY_nb, style = "C")
NY_listw_W <- nb2listw(NY_nb, style = "W")
NY_listw_U <- nb2listw(NY_nb, style = "U")

# summaries of each
summary(unlist(NY_listw_B$weights)) 
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##       1       1       1       1       1       1
```

```r
summary(unlist(NY_listw_C$weights))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.1839  0.1839  0.1839  0.1839  0.1839  0.1839
```

```r
summary(unlist(NY_listw_W$weights))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 0.09091 0.14290 0.16670 0.18390 0.20000 1.00000
```

```r
summary(unlist(NY_listw_U$weights))
```

```
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
## 0.0006545 0.0006545 0.0006545 0.0006545 0.0006545 0.0006545
```

```r
# sums of each
sum(unlist(NY_listw_B$weights))
```

```
## [1] 1528
```

```r
sum(unlist(NY_listw_C$weights))
```

```
## [1] 281
```

```r
sum(unlist(NY_listw_W$weights))
```

```
## [1] 281
```

```r
sum(unlist(NY_listw_U$weights))
```

```
## [1] 1
```

You can also use the `glist` parameter to define your own list of weights. It must take the form of a list of weight vectors corresponding to the `neighbours` list. This is perhaps the most useful. For example, suppose you think that the strength of neighbor relationships attenuates with distance (as is often the case). The weights then would be inversely proportional to some power of the distance. You can use `nbdists` to calculate distances for an `nb` object (it will only calculate distances between neighbors). 


```r
# supply both nb object and centroids
dists <- nbdists(NY_nb, coordinates(NY))

# this returns a list object of type nbdist, here the distances between the neighbors of unit 3
dists[3]
```

```
## [[1]]
## [1] 1160.921 1804.868 1685.356 1543.600 1677.131
```

```r
# now, calculate inverse distance weights, using lapply to calculate for all members of the dists lists
idw <- lapply(dists, function(x,p) 1/(x/1000)^p, 1)
NY_linkw <- nb2listw(NY_nb, glist = idw, style = 'B')

#this returns a list of weights
idw[3]
```

```
## [[1]]
## [1] 0.8613852 0.5540572 0.5933465 0.6478361 0.5962564
```

```r
# that are included in the NY_weights object
NY_linkw$weights[3]
```

```
## [[1]]
## [1] 0.8613852 0.5540572 0.5933465 0.6478361 0.5962564
```

```r
summary(unlist(NY_linkw$weights))
```

```
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
##  0.02602  0.10930  0.28990  0.45470  0.64820 12.09000
```

```r
# if the connection strength falls off more rapidly, we can use the inverse of the squared distance
idw <- lapply(dists, function(x,p) 1/(x/1000)^p, 2)
NY_linkw <- nb2listw(NY_nb, glist = idw, style = 'B')
summary(unlist(NY_linkw$weights))
```

```
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
##   0.00068   0.01195   0.08404   0.55800   0.42020 146.20000
```


# Measures of Spatial Autocorrelation

There are many measures of spatial autocorrelation. In general, we place them in to categories - global and local. To illustrate the difference, we will derive Moran's I statistic for which there is both a local and global version.

## Moran's _I_

__This section is under development...__

Moran's I is simply the ratio of the covariance between a unit and the weighted sum of it's neighbors (it's spatial lag) and the variance of the entire study area. The ratio is adjusted using the weights so that it falls between -1 and 1.

For unit I, the local Moran statistic is defined (by Luc Anselin) as...
$$I_{i} = \frac{n}{\sum_{j=1}^{n} w_{ij}}\frac{(y_{i} - \bar{y}){\sum_{j=1}^{n} (y_{j} - \bar{y})w_{ij}}}{\sum_{i=1}^{n}(y_i-\bar{y})^{2}}$$.

The first term is the weight scaling term while the second term is the covariance/variance ratio. It is important to note that $\sigma_{y}^2$ is the biased versus (using $n$ rather than $n-1$).

If the covariance term is positive then the unit and it's spatial lag are positively correlated to some degree. If negative then the are negatively correlated. Spatial lags are an important concept.


```r
# using the simple example from about
y = c(27,51,45,47,54,18,7,52)

# create the difference from mean vector
yb = y -  mean(y)

# the spatial lags can be created using the weight objects we created above. First the "B" (binary) version
lag.listw(W,yb)
```

```
## [1]  16.375   7.375   6.500 -26.500 -20.875   9.500  11.500 -42.875
```

Note that the lag vector is the sum of the values of it's neighbors (difference from the overal mean). If we use the row normalized `listw` object, we get the averge of the values.  This illustrates one of the differences between row normalization...


```r
lag.listw(W.rn, yb)
```

```
## [1] -19.625000   7.375000  -2.625000  -6.625000  -6.958333   2.375000
## [7]   2.875000 -14.291667
```

The Moran plot is used to display these...


```r
moran.plot(yb,listw = W)
```

![](Areal_Data_Analysis_files/figure-html/unnamed-chunk-48-1.png) 

The quadrants tell you the relative differences between a unit and it's neighbors. The upper left quadrant show units with values below the mean but with neighbors above. Upper right - both are above the mean. Lower left, both are below the mean. Lower right, unit is above but neighbors are below.

Now, let's manually construct the local Moran I statistic using the simple example.


```r
## Moran's I (local -> global)
n <- length(yb)

I <- vector(length=8)
for (local in seq(1:8)) {
  first.term <- n/sum(unlist(W.rn$weights[local])) 
  second.term <-  (yb[local] * sum(yb[unlist(W.rn$neighbours[local])] %*% unlist(W.rn$weights[local])))/sum(yb^2)
  
  I[local] <- first.term * second.term
}

# local Moran (LISA) for each unit
I
```

```
## [1]  0.74740969  0.35357043 -0.06939233 -0.22262672 -0.40841968 -0.16706805
## [7] -0.31559787 -0.73639503
```

The global Moran's I is the average of the locals...


```r
global.I <- sum(I)/n
global.I
```

```
## [1] -0.1023149
```

Compare this to the `moran.test` output

```r
moran.test(y,W.rn)
```

```
## 
## 	Moran's I test under randomisation
## 
## data:  y  
## weights: W.rn  
## 
## Moran I statistic standard deviate = 0.1761, p-value = 0.4301
## alternative hypothesis: greater
## sample estimates:
## Moran I statistic       Expectation          Variance 
##       -0.10231494       -0.14285714        0.05297952
```



```r
# usually call moran.test and moran.mc

# to test significance use inferential techniques based on probability distriubtions
# I.test <- moran.test(y,W,alternative='two.sided')
# I.test

# try with row normalization (effectively down-weighting the influence of D, F and G, the most connected polygons)
# I.test.rn <- moran.test(y,W.rn,alternative='two.sided')
```

Not much difference here, same results basically although if alpha=5% then one passes, one fails on strict adherence to inference. But here our "sample" is so small this is an insignificant difference. We'll look a more realistic result later.

## Geary's _C_

$$C = \frac{(n-1)\sum_{i}^{n}\sum_{j}^{n} w_{ij}(y_i-y_j)^2}{2(\sum_{i}^{n}\sum_{j}^{n}w_{ij},\sum_{i}^{n}(y_i-\bar{y})^2)}$$

Values range from 0 (perfect correlation) to 2 (perfect dispersion).

## Getis-Ord 

This is actually a _hot-spot_ tool as such identifies units which are different from their neightbors. It is, as such, a local static (no global equivalent). For a unit $i$

$$G^{*}_{i} = \frac{\sum_{j=1}^{n}w_{ij}x{j} - \bar{X}\sum_{j=1}^{n}w_{ij}}{S\sqrt{\frac{[n\sum_{j=1}^{n}w^{2}_{ij} - (\sum_{j=1}^{n}w_{ij})^2]}{n-1}}}$$

This is standardized score. Higher positive scores imply clustering of high positive values. Higher negative scores indicate clustering of high negative values.


```r
localG(y, W, zero.policy=NULL, spChk=NULL)
```

```
## [1]  0.85719462  0.54562160  0.43032440 -0.85676616 -0.59080479 -0.07576748
## [7] -0.32957653 -1.53751699
## attr(,"gstari")
## [1] FALSE
## attr(,"call")
## localG(x = y, listw = W, zero.policy = NULL, spChk = NULL)
## attr(,"class")
## [1] "localG"
```

## Join count statistics

Moran's I and Geary's C can only be applied to continuous data.  When dealing with categorical data, a measure called Joins-Count statistic is used (see Unwin and O'Sullivan, page 211) to assess clustering or dispersion. This approach is similar to many fragmentation statistics and measures the occurrence of similar neighbors versus dissimilar neibhgbors. For example, it quantifies the frequency of a forest/urban vs. forest/forest vs. forest/water vs. water/forest vs. water/urban contiguity. 

So, for example, in the simplist case a binary variable is mapped into two categories (_Black_ and _White_), such that a join or edge, is classified as either $WW (0-0)$, $BB (1-1)$, $BW (0-1)$. Join count statistics can tell you if

* the number of $BW$ joins is significantly lower than what we would expect by chance (_clustering_)
* the number of $BW$ joins is signficanlty higher than what we would expect by chance (_dispersion_)
* the number of $BW$ joins is approximately the same as what we would expect by chance

the respective probabilities of observing the two types of units are $$P_B = \frac{n_B}{n}    P_W = \frac{n-n_W}{n} = 1 - P_B$$ 

The probability of $BB$ and $WW$ in two adjacent cells are $$P_{BB} = P_{B}P_{B}  P_{WW} = (1 - P_{B})(1 - P_{B}) = (1 - P_{B})^2$$

The probability of $BW$ in two adjacent cells is $$P_{BW} = P_{B}(1 - P_{B}) + (1 - P_{B})P_{B} = 2P_P{B}(1 - P_{B})$$

We can count the number of joins (binary or multiple levels) and compare it to our expections $$E[BB] = \frac{1}{2}\sum_i\sum_j w_{ij}P_{B}^2$$ and
$$E[BW] = \frac{1}{2}\sum_i\sum_j w_{ij}2P_{B}(1 - P_{B})$$

The variance for each of these can be computed as well but is quite complex. Given the expection and variance of the join counts, we can compute a test statistic $$Z(BW) = \frac{BW - E[BW]}{\sqrt{\sigma^{2}_{BW}}}$$

This is what `joincount.test` does. We can also use `jointcount.mc` in much the same way we used `moran.mc` to randomly redistribute the values and compute a distribute of outcomes. `jointcount.multi` allows a non-binary variable of interest. The LULC landcover example shows a practical example of its use.


```r
# read lulc layer and convert numbers to factors (lulc classes found in xml file)
lulc<-readGDAL('IMAGE_DBO_LULC_FALL_CAST2006.tif')
```

```
## IMAGE_DBO_LULC_FALL_CAST2006.tif has GDAL driver GTiff 
## and has 2824 rows and 4550 columns
```

```r
lulc$band1 <- as.factor(lulc$band1)

# reduce size for this document (this is a case where AHPCC might be useful)
lulc <- lulc[1500:1700,1500:1700]
 
# create weight matrix (simple rook)
w.grid <- cell2nb(nrow(lulc),ncol(lulc), type='rook')
W <- nb2listw(w.grid)

lulc.jc <- joincount.multi(lulc$band1,W)
head(lulc.jc)
```

```
##         Joincount     Expected     Variance   z-value
## 11:11       76.00 4.673391e+00 2.266293e+00  47.37983
## 13:13        1.00 2.252475e-03 1.125555e-03  29.73975
## 31:31        5.50 1.227723e-02 6.129424e-03  70.09424
## 41:41      113.75 9.807178e-01 4.835793e-01 162.16495
## 51:51      749.00 1.633951e+02 6.766728e+01  71.18940
## 100:100  12594.75 8.982556e+03 4.985554e+02 161.77608
```


## A case study - North Carolina SIDS
A classic case-study in epidemiology examines the spatial characteristics of SIDs in North Carolina. 


```r
nc <- readOGR(dsn=".", layer="sids")
```

```
## OGR data source with driver: ESRI Shapefile 
## Source: ".", layer: "sids"
## with 100 features and 20 fields
## Feature type: wkbPolygon with 2 dimensions
```

```r
nc@proj4string <- CRS("+proj=longlat +datum=NAD27")

# transform to state plane system
nc_SP <- spTransform(nc, CRS("+init=epsg:3358"))

# first, QUEEN contiguity
IDs <- row.names(as(nc_SP, 'data.frame'))
nc_nbq <- poly2nb(nc,row.names=IDs)
 
# see what information you can get about the spatial connections...
nc_nbq
```

```
## Neighbour list object:
## Number of regions: 100 
## Number of nonzero links: 490 
## Percentage nonzero weights: 4.9 
## Average number of links: 4.9
```

```r
summary(nc_nbq)
```

```
## Neighbour list object:
## Number of regions: 100 
## Number of nonzero links: 490 
## Percentage nonzero weights: 4.9 
## Average number of links: 4.9 
## Link number distribution:
## 
##  2  3  4  5  6  7  8  9 
##  8 15 17 23 19 14  2  2 
## 8 least connected regions:
## 20 21 26 27 64 68 74 88 with 2 links
## 2 most connected regions:
## 48 62 with 9 links
```

```r
# cardinality (i.e. the number of connections for each row)
card(nc_nbq)
```

```
##   [1] 6 4 3 4 3 5 6 5 5 3 7 7 5 6 3 3 4 6 8 3 2 2 5 4 6 6 2 2 7 5 6 5 5 7 7
##  [36] 3 5 3 5 4 6 7 7 6 5 3 5 4 9 4 7 5 3 6 6 6 5 3 6 5 3 7 9 7 2 4 4 5 2 3
##  [71] 7 3 4 7 2 6 5 5 5 6 6 7 4 6 4 5 4 4 2 4 3 7 5 5 4 6 8 6 5 4
```

```r
# and plot...
plot(nc_SP,border='darkgrey', col='lightgrey')
plot(nc_nbq,coordinates(nc),col='blue', add=T)
title(main='Queen Contiguity Connectivity')
text(coordinates(nc), label=nc$FIPSNO, cex=0.5)
```

![](Areal_Data_Analysis_files/figure-html/unnamed-chunk-55-1.png) 

```r
# now, ROOK contiguity
nc_nbr <- poly2nb(nc,queen=F)
card(nc_nbr)
```

```
##   [1] 6 4 3 4 3 5 6 5 5 3 6 6 5 6 3 3 4 5 8 3 2 2 4 4 6 6 2 2 6 5 6 5 5 6 5
##  [36] 3 5 3 5 4 5 6 7 5 4 3 4 4 9 4 6 5 3 6 5 6 5 3 6 5 3 6 8 5 2 4 4 5 2 3
##  [71] 7 3 4 7 2 6 4 5 4 5 6 7 3 5 3 5 4 3 2 4 3 6 4 5 4 6 8 6 5 4
```

```r
plot(nc,border='gray')
plot(nc_nbr,coordinates(nc),col='green',add=T)
title(main='Rook Contiguity Connectivity')
text(coordinates(nc), label=nc$FIPSNO, cex=0.5)
```

![](Areal_Data_Analysis_files/figure-html/unnamed-chunk-55-2.png) 

```r
# first, use nbdists to compute the distances between the neighbors, e.g. k=1 again, using nbdists as we did a few lines up

# "W" = row-standardized, weights sum to 1 for all entities, weights represent a percentage influence of each neighbor on the entity.
nc_nbq_wr <- nb2listw(nc_nbq, style = "W")

# "B" = binary (ie. no normalization, 1=linked, 0=not-linked), sums of weights differ according to number of neighbors
nc_nbq_wb <- nb2listw(nc_nbq, style="B")
 
# "C" = complete set of weights for all links sum to the number of entities
nc_nbq_wc <- nb2listw(nc_nbq, style="C")
# 
# "U" = complete set of weights sum to unity 
nc_nbq_wu <- nb2listw(nc_nbq, style = "U")
 
# Spatial Autocorrelation - compute Moran's I using the different listw's you created for the NC dataset. First the SID74 variable (number of SID's cases by county in 1974)
moran.test(nc_SP$SID74,nc_nbq_wr)
```

```
## 
## 	Moran's I test under randomisation
## 
## data:  nc_SP$SID74  
## weights: nc_nbq_wr  
## 
## Moran I statistic standard deviate = 2.5192, p-value = 0.00588
## alternative hypothesis: greater
## sample estimates:
## Moran I statistic       Expectation          Variance 
##       0.147740529      -0.010101010       0.003925567
```

```r
moran.plot(nc_SP$SID74,nc_nbq_wr)
```

![](Areal_Data_Analysis_files/figure-html/unnamed-chunk-55-3.png) 

The slope of the line in the `moran.plot` result is simply `lm(wx~x)` (the least squares line fit to the polygon value versus the weighted sum of its neighbors). 

What effect does the rown normalization have on these results.

```r
# Style B (binary)
moran.test(nc_SP$SID74,nc_nbq_wb)
```

```
## 
## 	Moran's I test under randomisation
## 
## data:  nc_SP$SID74  
## weights: nc_nbq_wb  
## 
## Moran I statistic standard deviate = 2.1707, p-value = 0.01498
## alternative hypothesis: greater
## sample estimates:
## Moran I statistic       Expectation          Variance 
##       0.119089049      -0.010101010       0.003542176
```

This has a l arge effect - __3__ times more likely to be random. The normalized W creates greater influence of lesser connected counties. This always be investigated.
 
Sometimes, it's better to try a monte carlo permutation approach...


```r
nc.moran.perms <- moran.mc(nc_SP$SID74,nc_nbq_wr,nsim=1000)

# this creates an object that can be plotted
plot(nc.moran.perms)
```

![](Areal_Data_Analysis_files/figure-html/unnamed-chunk-57-1.png) 

We can also use join counts in the case of areal data. Note that, in this case we have to break, or classify, the continuous data into factors.


```r
# easy example designed by classifying NC SIDS data as above median and below median
# summary(nc_SP$SID74)
# sids.rank <- cut(nc_SP$SID74, breaks=c(-1,4,45),labels=c('below-median','above-median'))
# names(sids.rank) <- rownames(nc_SP$names)
# nc.jc.factor <- joincount.mc(sids.rank, nc_nbq_wr,nsim=1000)
# plot(nc.jc.factor)
```

## Use of caution in autocorrelation tests

When using any of these autocorrelation tests you should always consider these items.

1. Autocorrelation tests are highly sensitive to spatial patterning in the variable of interest from any source. But by assuming that the regression model removes such systematic spatial patterning, spatial autocorrelation tests do not always produce insights into the model.

2. All of these tests are also highly sensitive to one's choice of spatial weights. Where the weights do not reflect the true structure of spatial interaction, estimated autocorrelation (or lack thereof) may actually stem from model misspecification.

3. As originally designed, spatial autocorrelation tests assumed there are __no__ neighborless units in the study area. When this assumption is violated, the size of $n$ may be adjusted (reduced) to reflect the fact that some units are effectively being ignored. Not doing so will generally bias the absolute value for the autocorrelation statistic upward and the variance downward.

# Modeling spatial effects

## The Spatial AutoRegressive (SAR) model

At the beginning of the course, we summarized the basic linear model, $$y_i = X_i\beta + \epsilon_i$$ where $\epsilon_i ~ N(0,\Sigma), i = 1,...,n$. If the observed values are independent or have known correlation, $\Sigma$, then the Best Linear Uniformly Unbiased Estimate (BLUUE) of $\beta$, $$\hat{\beta} = (X^T\Sigma^{-1}X)^{-1} (X^T\Sigma^{-1}y)$$ and is also the solution that minimizes the weighted sum of residuals squared (the Weighted Least Squares Solution (WLESS). The variance of the parameter vector is $$Var[\hat{\beta}] = \sigma^2(X^T\Sigma^{-1}X)^{-1}$$ where $\sigma^2$ is the unknown error variance, estimated as $\tilde{\epsilon}^T\Sigma^{-1}\tilde{\epsilon}/(n-k)$, with $\tilde{\epsilon} = y - X\hat{\beta}$. This functionality is implemented in base R using the function `lm` (for linear model).

Let's look at an example using the Columbus Crime dataset. We are assuming that the `CRIME` observations are not correlated in any way. 


```r
columbus <- readOGR(dsn='.',layer = 'columbus')
```

```
## OGR data source with driver: ESRI Shapefile 
## Source: ".", layer: "columbus"
## with 49 features and 20 fields
## Feature type: wkbPolygon with 2 dimensions
```

```r
col.ols <- lm(CRIME ~ INC + HOVAL, data=columbus)
summary(col.ols)
```

```
## 
## Call:
## lm(formula = CRIME ~ INC + HOVAL, data = columbus)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -34.418  -6.388  -1.580   9.052  28.649 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  68.6190     4.7355  14.490  < 2e-16 ***
## INC          -1.5973     0.3341  -4.780 1.83e-05 ***
## HOVAL        -0.2739     0.1032  -2.654   0.0109 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 11.43 on 46 degrees of freedom
## Multiple R-squared:  0.5524,	Adjusted R-squared:  0.5329 
## F-statistic: 28.39 on 2 and 46 DF,  p-value: 9.341e-09
```

Here, we are attempting to explain `CRIME` as a linear function of `INC` and `HOVAL` (income and housing values). Furthermore, we are assuming that covariance of `CRIME` is $\Sigma = \sigma^{2}I_{n}$. That is, the independent variables are uncorrelated among themselves. Under these assumptions, the `lm` function provides a solution in which the intercept is significant as is income (as income declines, crime tends to rise). Housing values do not seem to affect crime level. The residuals, $\tilde{\epsilon}$, are stored in the object. 


```r
head(col.ols$residuals)
```

```
##           0           1           2           3           4           5 
##   0.3465419  -3.6947990  -5.2873940 -19.9855151   6.4475490  -9.0734793
```

```r
stem(col.ols$residuals)
```

```
## 
##   The decimal point is 1 digit(s) to the right of the |
## 
##   -3 | 4
##   -2 | 0
##   -1 | 6543
##   -0 | 9999886665544432222
##    0 | 011134456779
##    1 | 01112333456
##    2 | 9
```

The stem plots shows a normal distribution of residuals. However, we know that residuals also have spatial locations (the census block groups in Columbus). Let's plot them that way. `sp` provides an easy function for doing this...


```r
# check that the residuals are in the same order as the original data.frame (they are)

# add the residuals to the data frame in columbus
columbus@data[,'residuals'] <- col.ols$residuals
spplot(columbus, zcol='residuals')
```

![](Areal_Data_Analysis_files/figure-html/unnamed-chunk-61-1.png) 

It is somewhat easy to notice here that the residuals seem to show some spatial clustering. We can, of course, test this using `moran.plot`. We need to define a `listw` object defining the neighborhoods. Let's assume `QUEEN` contiguity and use a row normalized weight list. Therefore, we are comparing the average of the neighborhood residuals to the units residuals.


```r
col.wq <- nb2listw(poly2nb(columbus, queen =TRUE), style='W')
moran.plot(col.ols$residuals,listw = col.wq)
```

![](Areal_Data_Analysis_files/figure-html/unnamed-chunk-62-1.png) 

Clearly this is positive global autocorrelation and the global `moran.test` shows it.


```r
col.moran.mc <- moran.mc(col.ols$residuals, listw = col.wq, nsim = 100)
plot(col.moran.mc)
```

![](Areal_Data_Analysis_files/figure-html/unnamed-chunk-63-1.png) 

Finally, we can use the following function applied directly to the `lm` object. 


```r
lm.morantest(col.ols, col.wq, zero.policy=NULL, alternative = "greater", spChk=NULL, resfun=weighted.residuals)
```

```
## 
## 	Global Moran's I for regression residuals
## 
## data:  
## model: lm(formula = CRIME ~ INC + HOVAL, data = columbus)
## weights: col.wq
## 
## Moran I statistic standard deviate = 2.8393, p-value = 0.00226
## alternative hypothesis: greater
## sample estimates:
## Observed Moran's I        Expectation           Variance 
##        0.222109407       -0.033418335        0.008099305
```

So there seems to be some level (though not necessarily a significant amount - depends on your threshold) of autocorrelation remaining in the residuals. We have at least two choices:

1. We could try to re-specify the model (with different dependent variables or transformation of dependent variables) and see if the autocorrelation was actually due to spatial autocorrelation of a dependent variable or variables. This would be true if you think the demographic values can completely predict crime rates.

2. We could try a spatial autoregessive model (which we'll develop below) which assumes the the crime rates themselves are spatially autocorrelated. This is the path you would take after exhausting the limits of the model (choice 1) or if you make the assumption that some interaction across units (census blocks) is cause crime to be spatially autocorrelated. 

The choices you make (and likely it will be an iteration of these choices) may depend on your knowledge of the domain. For example, your understanding of how crime propagates through a city will affect your choice.

When you begin working on choice 2, your task is to build an appropriate covariance matrix that models spatial interactions of crime rates.

Beginning with the model detailed above $$y_i = X_i\beta + \epsilon_i$$ and explore what happens when we assume that two neighbors $i$ and $j$ interact...

$$y_{i} = \alpha_{j}y_{j} + X_{i}\beta + \epsilon_{i} \\
y_{j} = \alpha_{i}y_{i} + X_{j}\beta + \epsilon_{j} \\
\epsilon_{i} \sim N(0,\sigma^2), i = 1 \\
\epsilon_{j} \sim N(0,\sigma^2), j = 2 $$

essentially meaning that observed values at location $i$ depend on those at location $j$, and vice versa. We also assume here that the data generating process is "simultaneous". Time dependent models are harder.

With $n$ observations, we generalize the interactions as

$$y_{i} = \rho \sum^{n}_{j=1}W_{ij}y_{j} + X_{i}\beta + \epsilon_{i} \\
\epsilon_{i} \sim N(0,\sigma^2), i = 1,...,n$$

This equation is for each $i$ observation and it's interactions with all the other observations at different locations. For all the observations we can write this a matrix equation (as before)

$$y = \rho Wy + X\beta + \epsilon \\
\epsilon \sim N(0, \sigma^2I_{n}) $$

where $W$ is the spatial weights matrix, $\rho$ is the _spatial autoregressive scalar parameter_ and $I_{n}$ is the $n x n$ identify matrix. The absence of subscripts imply a vector or matrix. 

You will often see this written as $$(y - \rho W y) = X \beta + \epsilon$$. The new term introduced to the model $\rho W y$ is the spatial lag we have seen before.

In order to get an idea of the effect of $\rho$ consider the following :

* When $\rho = 0$, the variable is not spatially autocorrelated. Information about a measurement in one location gives us no additional information about the value in neighboring locations (this is spatial independence and the model reverts to ordinary least squares).
* When $\rho > 0$, the variable is positively spatially spatially autocorrelated. Neighboring values tend to be simiar to each other (clustering).
* When $\rho < 0$, the variable is negatively negatively spatially autocorrelated. Neighboring values tend to be different from each other (segregation).

Consider the following simulation using the columbus `QUEEN` conitiguity and `style='W'` row normalized weights. We are bypassing the `linkw` object here in favor of the $W$ form.  


```r
# # create the nb object
# col.nb <-poly2nb(columbus, queen =TRUE) 
# 
# # define n and create a vector randomly distributed independent variables
# y <- rnorm(n,mean=0,sd=1)  # column vector
# 
# # create three different spatial lags by changing rho
# lag.zero <- invIrM(col.nb, 0.01, style='W',method='solve') %*% y
# lag.pos <- invIrM(col.nb, 0.9, style='W',method='solve') %*% y
# lag.neg <- invIrM(col.nb, -0.9, style='W',method='solve') %*% y
# 
# # now check with a Moran spatial autocorrelation test
# moran.test(lag.zero, col.wq,alternative='two.sided')
# moran.test(lag.pos, col.wq,alternative='two.sided')
# moran.test(lag.neg, col.wq,alternative='two.sided')
```

So it turns out that spatial interaction goes a bit beyond just a spatially structured covariance matrix. The $\beta$ parameters in our model reflect the short-run direct impact of $X_{i}$ on $y_{i}$. However, we must account for the indirect impact of $X_{i}$ on $y_{i}$, from the influence $y_{i}$ exerts on it's neighbors $y_{i}$ which in turn feeds back into $y_{i}$. We've created a highly non-linear model and weighted least squares no longer gives an unbiased estimate. The solution to the spatial autoregressive (SAR) models is very complex. Luckily, `R` handles the complexities for us. 


```r
# SAR model for Columbus crime data (using QUEEN contiguity and row-normalized weights)
col.sar <- lmSLX(CRIME ~ INC + HOVAL, data=columbus, listw = col.wq)

# compare to OLS
summary(col.sar)
```

```
## 
## Call:
## lm(formula = nfo, data = data, na.action = na.action)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -36.342  -7.662  -0.013   7.978  25.572 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  74.5534     6.7156  11.101 2.42e-14 ***
## INC          -1.0974     0.3738  -2.936  0.00528 ** 
## HOVAL        -0.2944     0.1017  -2.896  0.00587 ** 
## WX.INC       -1.3987     0.5601  -2.497  0.01633 *  
## WX.HOVAL      0.2148     0.2079   1.033  0.30712    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 10.91 on 44 degrees of freedom
## Multiple R-squared:  0.6105,	Adjusted R-squared:  0.5751 
## F-statistic: 17.24 on 4 and 44 DF,  p-value: 1.413e-08
```

```r
summary(col.ols)
```

```
## 
## Call:
## lm(formula = CRIME ~ INC + HOVAL, data = columbus)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -34.418  -6.388  -1.580   9.052  28.649 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  68.6190     4.7355  14.490  < 2e-16 ***
## INC          -1.5973     0.3341  -4.780 1.83e-05 ***
## HOVAL        -0.2739     0.1032  -2.654   0.0109 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 11.43 on 46 degrees of freedom
## Multiple R-squared:  0.5524,	Adjusted R-squared:  0.5329 
## F-statistic: 28.39 on 2 and 46 DF,  p-value: 9.341e-09
```


## The Spatial Error Model (SEM)

Suppose that $y$ is explained entirely by two variables $x$ and $z$, where $x,z \sim N(0, I_{n})$ and are indepdendent $y = x\beta + z\theta$.

If $z$ is not observed the vector $z\theta$ is nested in teh error term $\epsilon$ and $y = x\beta + \epsilon$. $z$ is any variable that is unaccounted for in the model and is often called a _latent_ variable. Often, it is a variable that is hard to quantify.

In addition the bias induced by latent variable (omitted variable), a motiviation for this Spatial Error Model or SEM, is spatial heterogneity. Suppose we have a panel data set with multiple observations in each unit and we want our model to incorporate individual effects, we can include an $n x 1$ vector $a$ of individual intercepts for each unit $y = a + X\beta$. We can treat $a$ as a vector of spatial random effects and assume that it follows a spatial autoregressive process $$a = \lambda W a + \epsilon \\ a = (I_{n} - \lambda W)^{-1}\epsilon$$ where $\epsilon \sim N(0, \sigma^2I_{n})$ is a vector of disturbances. Substituting this into the model we get 

$$y = X\beta + a = X\beta + (I_{n} - \lambda W)^{-1}\epsilon$$



## Continuation of the SIDS case-study

### Probability Mapping
We will focus on probability mapping for disease rates data. Typically, we have counts of the incidence of some disease by spatial unit, associated with counts of populations at risk. The task is then to try to establish whether any spatial units seem to be characterised by higher or lower counts of cases than might have been expected in general terms (Bailey and Gatrell, 1995).
An early approach by Choynowski (1959), described by Cressie and Read (1985)
and Bailey and Gatrell (1995), assumes, given that the true rate for the spatial units is small, that as the population at risk increases to infinity, the spatial unit case counts are Poisson with mean value equal to the population at risk times the rate for the study area as a whole. Choynowski’s approach folds the two tails of the measured probabilities together, so that small values, for a chosen $\alpha$, occur for spatial units with either unusually high or low rates. For this reason, the high and low counties are plotted separately


```r
ch <- choynowski(nc$SID74, nc$BIR74)
nc$ch_pmap_low <- ifelse(ch$type, ch$pmap, NA)
nc$ch_pmap_high <- ifelse(!ch$type, ch$pmap, NA)
prbs <- c(0, 0.001, 0.01, 0.05, 0.1, 1)
nc$high = cut(nc$ch_pmap_high, prbs)
nc$low = cut(nc$ch_pmap_low, prbs)

#Probability map of North Carolina counties, SIDS cases 1974–78, alpha = 0.05, reproducing Cressie and Read (1985)
spplot(nc, c("low", "high"), col.regions = grey.colors(5))
```

![](Areal_Data_Analysis_files/figure-html/unnamed-chunk-67-1.png) 

```r
pmap <- probmap(nc$SID74, nc$BIR74)
nc$pmap <- pmap$pmap
brks <- c(0, 0.001, 0.01, 0.025, 0.05, 0.95, 0.975, 0.99,0.999, 1)
library(RColorBrewer)
spplot(nc, "pmap", at = brks, col.regions = rev(brewer.pal(9,"RdBu")))
```

![](Areal_Data_Analysis_files/figure-html/unnamed-chunk-67-2.png) 

```r
 hist(nc$pmap, main = "")
```

![](Areal_Data_Analysis_files/figure-html/unnamed-chunk-67-3.png) 

One ad-hoc way to assess the impact of the possible failure of our assumption
that the counts follow the Poisson distribution is to estimate the dispersion by fitting a generalized linear model of the observed counts including only the intercept(null model) and offset by the observed population at risk (suggested by Marilia Carvalho and associates):


```r
res <- glm(SID74 ~ offset(log(BIR74)), data = nc, family = "quasipoisson")
nc$stdres <- rstandard(res)
brks <- c(-4, -3, -2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2, 3, 4)
spplot(nc, "stdres", at = brks, col.regions = rev(brewer.pal(11,"RdBu")))
```

![](Areal_Data_Analysis_files/figure-html/unnamed-chunk-68-1.png) 

### Exploration of the SIDs Data (following Cressie and Read, 1985)

One of the first steps taken by Cressie and Read (1985) is to try to bring out spatial
trends by dividing North Carolina up into 4×4 rough rectangles. Just to see how this
works, let us map these rough rectangles before proceeding further.


```r
nc$both <- factor(paste(nc$L_id, nc$M_id, sep = ":"))
nboth <- length(table(unclass(nc$both)))
spplot(nc, "both", col.regions = sample(rainbow(nboth)))
```

![](Areal_Data_Analysis_files/figure-html/unnamed-chunk-69-1.png) 

Cressie constructs a transformed SIDS rates variable, 1974–78, for his analyses. We can replicate his stem-and-leaf figure on p. 396 in the book, taken from Cressie and Read (1989):


```r
nc$ft.SID74 <- sqrt(1000) * (sqrt(nc$SID74/nc$BIR74) + sqrt((nc$SID74 + 1)/nc$BIR74))
stem(round(nc$ft.SID74, 1), scale = 2)
```

```
## 
##   The decimal point is at the |
## 
##   0 | 9
##   1 | 111244
##   1 | 567789999
##   2 | 0011111222334444
##   2 | 55555666677778999999999
##   3 | 000111122333333344444444
##   3 | 5568999
##   4 | 013344
##   4 | 555557
##   5 | 2
##   5 | 
##   6 | 3
```

### Median polish smoothing

Cressie (1991, pp. 46–48, 393–400) discusses in some detail how smoothing may be
used to partition the variation in the data into smooth and rough. In order to try it out on the North Carolina SIDS data set, we will use a coarse gridding into four columns and four rows given by Cressie (1991, pp. 553–554), where four grid cells are empty; these are given by variables $L_{id}$ and $M_{id}$ in object `nc`. Next we aggregate the number of live births and the number of SIDS cases 1974–1978 for the grid cells:

```r
mBIR74 <- tapply(nc$BIR74, nc$both, sum)
mSID74 <- tapply(nc$SID74, nc$both, sum)
```

Using the same Freeman-Tukey transformation as is used for the county data, we
coerce the data into a correctly configured matrix, some of the cells of which are empty. The medpolish function is applied to the matrix, being told to remove empty cells; the function iterates over the rows and columns of the matrix using median to extract an overall effect, row and column effects, and residuals:



