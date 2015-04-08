## Field Experiment Modeling Day
## April 8 2015
## Population Ecology Spring 2015
## Auriel Fournier

## this will be different because we won't have the same matrix all the time
## before we've had for loops looping over a matrix varying a certain parameter

## lambda in a density dependent model won't work because lambda is going to
## vary a lot, wiht dramatic growth and then near 0 at equillibrium

## adult population size is the primary end point
## equilliburium population size

## we cannnot rely on eiganvalues and eiganvectors

## for each value of a parameter we are going to run a time series and 
## take the parameter at teh end of the sensativity analysis

## matrix will now be IN the for loop, instead of above it

## four element matrix model
##   
## Juv  (juv-juv)     (adult-juv) 
## Adl  (juv - adult) (adult-adult)

## J	0			Fecundity(Adults) - Fa
## A 	Survival(Juv) St	Survival (Adults) - Sa

# needs

## num of tadpoles produced Tt = sex ratio * Adults(t) * Fecundity * Survival(eggs)

## Survival of tadpoles St = Stmax /(1-D*Tadpoles(t))^gamma
		## gamma is what tells you if it's over or under compensating 

## Fecundity (Fa) clutch size * survival(eggs) * survival(tadpoles)(function of time) * survival(metamorphs) * sex ratio
phi = 5000
se = .6
sm = .2
sj = .2
sa = .75
stmax = .8
p = .5
d = .05
gamma = 1
tmax = 200

mat <- matrix(c(0,sj,ft,sa),ncol=2, nrow=2)

N <- matrix(nrow=2,ncol=tmax)
N[,1] <- 1


for(t in 2:tmax){
tt = N[2,t-1] * p * phi * se
st = (stmax)/((1+d*tt)^gamma)
ft = p * phi * se * st * sm
mat <- matrix(c(0,sj,ft,sa),ncol=2, nrow=2)
newN <- mat %*% N[,t-1]
N[,t] <- newN
}


plot(seq(1:tmax), N[2,])

max(N[2,])



####

## Sensativity

####


phi = 5000
se = .6
sm = .2
sj = .2
sa = seq(0,1,by=.01)
stmax = .8
p = .5
d = .05
gamma = 1
tmax = 1000

mat <- matrix(c(0,sj,ft,sa),ncol=2, nrow=2)

N <- matrix(nrow=2,ncol=tmax)
N[,1] <- 1

savad <- matrix(ncol=2, nrow=length(sa))
savad[,1] <- sa

for(i in 1:length(sa)){
	for(t in 2:tmax){
		tt = N[2,t-1] * p * phi * se
		st = (stmax)/((1+d*tt)^gamma)
		ft = p * phi * se * st * sm
		mat <- matrix(c(0,sj,ft,sa[i]),ncol=2, nrow=2)
		newN <- mat %*% N[,t-1]
		N[,t] <- newN
	}
savad[i,2] <- newN[2,]
}




###

## egg survival

###


phi = 5000
se = seq(0,1,by=.01)
sm = .2
sj = .2
sa = .75
stmax = .8
p = .5
d = .05
gamma = 1
tmax = 1000

mat <- matrix(c(0,sj,ft,sa),ncol=2, nrow=2)

N <- matrix(nrow=2,ncol=tmax)
N[,1] <- 1

savegg <- matrix(ncol=2, nrow=length(se))
savegg[,1] <- se

for(i in 1:length(se)){
	for(t in 2:tmax){
		tt = N[2,t-1] * p * phi * se[i]
		st = (stmax)/((1+d*tt)^gamma)
		ft = p * phi * se[i] * st * sm
		mat <- matrix(c(0,sj,ft,sa),ncol=2, nrow=2)
		newN <- mat %*% N[,t-1]
		N[,t] <- newN
	}
savegg[i,2] <- newN[2,]
}



###

## plots

###

par(mfrow=c(3,1))
plot(seq(1:tmax), N[2,], main="initial model")
plot(savegg, main="egg survival")
plot(savad, main="adult survival")




