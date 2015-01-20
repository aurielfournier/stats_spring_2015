# review and download the data file BayouMetoQA.csv

#first look at the csv file in a text editor, note the header, the commas, missing data, etc.

#							PREPARE AND INSPECT THE DATA

# set working directory (to wherever you put your data file...)
setwd("~/SourceTree/Stats_Spring_2015")

# verify working directory
getwd()

# load data from a comma-separated-value formatted file
check = read.csv('BayouMetoQA.csv')

# or you could use the slightly more complicated read.table function
check = read.table('BayouMetoQA.csv',header=TRUE,sep=',')

# examine the resulting data.frame
head(check)

# another, more detailed digest of the data
str(check)

# extract the number of rows and cols to a variable (maybe we'll need it later)
n.obs = nrow(check)
n.attributes = ncol(check)

# add a column to the data.frame containing the delta between the GPS and LiDAR Z's.  NOTE: Negative values indicate the lidar return is higher than the GPS return (we'll maintain these values because there is information here)
check$delta = check$Z-check$LidarZ
head(check)
str(check)

# one of the most effective ways to explore the distribution, called the stem and leaf plot, is now taught in 6th grade
stem(check$delta)		# is this left skewed or right skewed?

# note the one distinctly different value, what about in absolute value terms?
stem(abs(check$delta))		# even more pronounced, sRtill skewed?


# discuss observation or data types: continuous (ordinal,interval,ratio)and categorical (ordinal,nominal/discrete)


#							CREATE TABLES CLASSIFYING OBS BY TYPE

# look at the frequency distribution of the three categorical variables (AKA Factors) - Region, LC, Class

# first by Region

# apply the table functionbar
Region.freq = table(check$Region)
Region.freq

# convert to column format for easier reading
cbind(Region.freq)

# do the same for LC (land cover) type
LC.freq = table(check$LC)

# display as a column, but don't save
cbind(LC.freq)

# we may also look at a joint frequency distribution
joint = table(check$LC, check$Region)

# what kind of data type is "joint"?
is.data.frame(joint)
is.character(joint)
is.vector(joint)
is.matrix(joint)   #yes, with colnames and rownames defined

# or simply query it...
class(joint)

# we can extract data from the table or the data.frame in a convenient way using indexing
check[,'LC']	# land cover classes by observation
check[1,]		# first GPS point and all attributes
check[check$LC=='Asphalt',]		# all attributes of GPS points collected on Asphalt
check[abs(check$delta)>0.18,] 	# all attribures of GPS points differing from Lidar by >18cm

joint['Asphalt',]	# number of Asphalt point
joint[,'North']		# nubmer of points in the North by type, etc.
 
# you can also generate proportion tables, 1=sum by row, 2=by column, '' by total
prop.table(joint,1)
prop.table(joint,2)
prop.table(joint)

# to save the relative frequencey distribution of the LC types to variable
LC.relfreq = LC.freq / n.obs


#						DISPLAY TABLES GRAPHICALLY

# list the current variables in the workspace
ls()

# enough tables and numbers, create some graphs

# simple barplot
barplot(LC.freq)
colors = c('black','brown','green','pink','cyan')
barplot(LC.freq,col=colors)

# note that the plot command is over-ridden for data.frames depending on inputs if we use the formula data type which is creating using the ~ symbol.  ~ can be read as "desribed by"
class(check$delta ~ check$LC)

# using the formula datatype for continuous data:
plot(check$delta ~ check$Z)		# creates a scatterplot

# while using the formual with an factor as the independent variable:
plot(check$delta ~ check$LC)	# creates a multi boxplot
plot(check$delta ~ check$Region)

# pie charts are, in general, better for showing relative frequency graphically
pie(LC.freq,col=colors)


#						DESCRIBE THE DATA

# the presence of categories implies an expected behavior difference (and we saw that above) so we may summarize continuous variables by category

# to dig a little deeper to better understand what R is doing, create a logical vector for the presence of a particular category
Asphalt.index = check$LC=="Asphalt"
Asphalt.index
check[Asphalt.index,]

Forest.index = check$LC == 'Forest'
Forest.index
check[Forest.index,]

# compare the means elevation delta's. Here we use Asphalt.index to index the rows and the column name "delta" to select the column
mean(check[Asphalt.index,'delta'])
mean(check[Forest.index,'delta'])

# this is called slicing, notice the substantially different means (but is the difference significant given the variation we see? - this is what ANOVA tells us)

# shouldn't there be an easier way of getting summaries by category? of course there is...
tapply(check$delta,check$LC,mean)
cbind(tapply(check$delta,check$LC,mean))

# also try lapply (returns list) and sapply (returns simple string)
#what about variance within the LC type
cbind(tapply(check$delta,check$LC,var))

# or, more intuitively, standard deviation
cbind(sqrt(tapply(check$delta,check$LC,var)))
cbind(tapply(check$delta,check$LC,sd))

# or range
cbind(tapply(check$delta,check$LC,range))
cbind(tapply(check$delta,check$LC,min))
cbind(tapply(check$delta,check$LC,max))
cbind(tapply(check$delta,check$LC,quantile))

# notice that the large variation in scrub, grass and forest is due to a large -/+ range

# lets look only at absolute differences
cbind(tapply(abs(check$delta),check$LC,sd))

# is there a large difference in mean and variation across flights?
cbind(tapply(check$delta,check$Region,mean))
cbind(tapply(check$delta,check$Region,var))
cbind(tapply(abs(check$delta),check$Region,mean))
cbind(tapply(abs(check$delta),check$Region,var))

# evidently, yes. but perhaps the LC distribution is different. this is a job for ANOVA (later in the semester)

# finally, before we move to analyzing the distributions of the continuous variables, note the very useful summary function
tapply(check$delta,check$Region,summary)

# remember to look by at these summary results when discussing skewness and kurtosis
 
# frequency distribution of quantitative data

# we may reduce typing by attaching "check" to the workspace. we may refer to the columns without the $ notation
attach(check)

# verify the range of the data
range(delta)

# in order to look at frequency distributions, we need to group the continuous variables
breaks = seq(-0.5,0.15,by=0.125)
breaks

# bin, or collect, the delta values according to these sub-intervals, the cut function will define intervals based on breaks and the classify each element in "delta" by the classifications
delta.cut = cut(delta,breaks,right=FALSE)

# note that cut creates factors with length(breaks)-1 levels
delta.cut

# now we can tabulate the delta values by the interval in which they fall using the table function which counts by levels
delta.freq = table(delta.cut)
cbind(delta.freq)

# it is more efficient to use the hist (histogram) function:
# produces a very useful histogram object and can compute breaks using the rules described in BBR ('Scott','Sturges',others)
delta.hist = hist(delta,breaks='Scott',main='GPS - LiDAR')

# look at delta.hist closely...
summary(delta.hist)
str(delta.hist)

# ...and extract midpoints of the intervals
delta.hist$mids

# QUESTION: why are both intensities and densities and why are they the same in our dataset?

# now compute relative frequency distribution
delta.relfreq = delta.freq / n.obs		# use n.obs from earlier computation
cbind(delta.freq,delta.relfreq)			# why is this different from hist densities/intensities? - breaks are different.

# lets compute an empirical quantile plot, using the formula notation
p=seq(0,1,0.01);
r=quantile(delta,probs=p);
plot(r~p,type='l',xlab='percent',ylab='value')

# What happens as breaks are redefined?

# how many breaks should we use? try Scott's rule (page 49 of BBR)
bin.width = 3.5*sd(delta)*nrow(check)^(-1/3)
bin.width

brk.a = seq(-0.5,0.50,by=0.2)
hist(delta,breaks=brk.a,col='blue',main='-0.5 to 0.5 by 0.125 (brk.a)')

# the bin width is slightly smaller than the separation between the two most distance "adjacent" values...
diff(sort(delta))

# if one of these two data points happen to like near a breakpoint, then there will be an empty bin making observation look like an outlier.
brk.b = seq(-0.55,0.8,by=0.125)
hist(delta,breaks=brk.b,col='red',main='-0.55 to 0.8 by 0.125 (brk.b)')

# another seemingly small change has a significant effect on the distribution - whereas brk.a looks symmetric, this distribution looks left skewed
brk.c = seq(-0.426,0.8,by=0.125)
hist(delta,breaks=brk.c,col='green',main='-0.426 to 0.8 by 0.125 (brk.c)')

# it is a simple matter to look at cumulative frequency distributions using the cumsum function
delta.cumfreq = cumsum(delta.freq)
cbind(delta.cumfreq)

# plot the cumulative distribution
cumfreq0 = c(0,cumsum(delta.freq))
plot(breaks,cumfreq0,type='l',main='Cumulative Distribution of Delta')
lines(brk.a,cumfreq0)

# now, it is perhaps more apparent that the cumulative frequency of absolute values would make more sense. try it with absolute values

#we've seen that delta does depend on flight and LC.  does is also depend on position (you would expect yes since the flights are N and S) 

#simple scatter plots will help
plot(check$Z,delta,xlab='Elevation',ylab='Delta')
plot(check$X,delta,xlab='Easting',ylab='Delta')
plot(check$Y,delta,xlab='Northing',ylab='Delta')
 
# finally, let's look at the delta's summarized graphically
boxplot(delta~check$LC)
boxplot(abs(delta)~check$LC)

#quantile calculations in small samples such as this can be deceiving so the IQR results are likewise misleading.  In this case a stripchart is useful for plotting the raw data
stripchart(abs(check$delta) ~ check$LC,method='jitter')
 
#note the distributions by land cover - the veg classes show more negative values (why?)
