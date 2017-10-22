# (c) Liisa Loog 2017
#
# Code for Sliding Window Analyses of Mobility through time.
#
# Ref. Loog et. al. (2017) Estimating mobility using sparse data: Application to  human genetic variation. PNAS DOI: XXX


library(vegan)
GEN.DIST.MAT <- read.csv("GEN.DIST.MAT.csv", head = T) # Loading the Genetic distance matrix
GEO.DIST.MAT <- read.csv("GEO.DIST.MAT.csv", head = T) # Loading the Geographic distance matrix
data <- read.csv("DATE_DATA.csv", head = T) # Loading sample date information
version <- 1.0 # Version nr for the output file
nindividuals <- length(data[,1])

############################## START ##############################
nr.run <- 100 # Number of resampling runs
window.moves = 125 # Number of windows

results.array = array(data=NA, dim=c(window.moves,14,nr.run))
corvec.array = array(data=NA, dim=c(window.moves,200,nr.run))

for (m in seq(1,nr.run)){
	print(paste("run nr ", m))
############################## Draw Dates From Sample Date Ranges ##############################

	datefrom <- as.numeric(data[,3])
	dateto <- as.numeric(data[,2])
	random.date <- vector(mode="numeric", length= nindividuals)

	for(i in seq(1, nindividuals)){
		rdate = runif(1, datefrom[i], dateto[i])
		random.date[i] <- rdate
	}

	data <- cbind(data[,1:3], as.numeric(random.date))

##############################  CALCULATE TEMPORAL DISTANCES ##############################

	TIME.DIST.MAT <- matrix(nrow= nindividuals,ncol= nindividuals)

	for(i in seq(1,nindividuals)){
#	print(i)
		for(j in seq(1,nindividuals)){
			TIME.DIST.MAT[i,j] <- abs(data[i,4] - data[j,4])
		}
	}
############################## SLIDING WINDOW PARAMETERS ##############################

	start.windows.at = 16000 # Starting time of the Sliding Window analyses (Years Ago)
	gap.window = 4000 # Width of the sliding window (in Years)
	move.window.by = 100 # Stride (in Years)

	for(w in seq(1: window.moves)) {
############################## Define the Window ##############################

		time.from = start.windows.at
		time.to = start.windows.at - gap.window
		start.windows.at = time.from - move.window.by

		data.positions <- which((data[,17] <= time.from) & (data[,17] >= time.to))
		nsamples <- length(data.positions)

#		data.positions <- sample(data.positions,1, replace = F, ) # <- Use this line for jackknifing

############################## Genetic Distances ##############################
		
		gen.dist.mat <- matrix(NA, nsamples, nsamples)
		gen.dist.mat <- GEN.DIST.MAT[data.positions,data.positions]

############################## Geographic Distances ##############################

		geo.dist.mat <- matrix(NA, nsamples, nsamples)
		geo.dist.mat <- GEO.DIST.MAT[data.positions,data.positions]

###################################### Temporal Distances ###################################
		
		time.dist.mat <- matrix(NA, nsamples, nsamples)
		time.dist.mat <- TIME.DIST.MAT[data.positions,data.positions]

############################## MANTEL TESTING GEOG/TIME EUCLUD DERRIV ##############################

        nsim=200 # number of angles where correlation is calculated
		mcor_resvec = rep(0,nsim)

		scalevec <- tan(seq(0, ((pi/2) * (1-1e-8)), length.out= nsim))
		ang <- seq(0, (90 * (1-1e-8)), length.out= nsim)

		for(i in seq(1,nsim)){
			geotimedist <- (((geo.dist.mat )^2) + ((time.dist.mat * scalevec[i])^2))^0.5
			mcor_resvec[i] <- mantel(gen.dist.mat, geotimedist, permutations=0, method="pearson")$statistic
		}
		mcor_resvec[is.na(mcor_resvec)] <- 0
		maxratio = ang[mcor_resvec==max(mcor_resvec)]
		maxcor = max(mcor_resvec)
		mincor = min(mcor_resvec)
		ends = c(mcor_resvec[1],  mcor_resvec[nsim])
		peak = maxcor - max(ends)
		if(length(maxratio) > 1){
			maxratio = 0
		}

############################## RECORD RESULTS  ##############################
        results.array[w,1,m] <- w # window number
        results.array[w,2,m] <- time.from # Starting time of the window (Years Ago)
		results.array[w,3,m] <- time.to # End time of the window (Years Ago)
        results.array[w,4,m] <- nsamples # Number of samples in a window
        results.array[w,5,m] <- maxratio # Angle value at the point of highest correlation
        results.array[w,6,m] <- maxcor # Highest correlation
        results.array[w,7,m] <- peak # Difference between the highest correlation and the correlation between the highest of the ends (correlation eith either geographic or temporal distance matrix)
        results.array[w,8,m] <- mcor_resvec[1] # Correlation with the geographic distance matrix
		results.array[w,9,m] <- mcor_resvec[nsim] # Correlation with temporal distance matrix

        corvec.array[w,,m] <-mcor_resvec # Vector of correlation values given each angle (n = nsim)


############################## PERMUTE MATRICES FOR P-VALUE CALCULATION ##############################

############################## Geographic Distances ##############################

		data.positions.geo <- sample(data.positions, nsamples, replace = F, )
		p.geo.dist.mat <- GEO.DIST.MAT[data.positions.geo,data.positions.geo]

###################################### Temporal Distances ###################################
		
		data.positions.time <- sample(data.positions, nsamples, replace = F, )
		p.time.dist.mat <- TIME.DIST.MAT[data.positions.time,data.positions.time]

############################## MANTEL TESTING GEOG/TIME EUCLUD DERRIV ##############################

		nsim=200 # number of angles where correlation is calculated
		p.mcor_resvec = rep(0,nsim)

		scalevec <- tan(seq(0, ((pi/2) * (1-1e-8)), length.out= nsim))
		ang <- seq(0, (90 * (1-1e-8)), length.out= nsim)

		for(i in seq(1,nsim)){
			geotimedist <- (((p.geo.dist.mat )^2) + ((p.time.dist.mat * scalevec[i])^2))^0.5
			p.mcor_resvec[i] <- mantel(gen.dist.mat, geotimedist, permutations=0, method="pearson")$statistic
		}
		p.mcor_resvec[is.na(p.mcor_resvec)] <- 0
		p.maxratio = ang[p.mcor_resvec==max(p.mcor_resvec)]
		p.maxcor = max(p.mcor_resvec)
		p.mincor = min(p.mcor_resvec)
		p.ends = c(p.mcor_resvec[1],  p.mcor_resvec[nsim])
		p.peak = p.maxcor - max(p.ends)
		if(length(p.maxratio) > 1){
			p.maxratio = 0
		}
		

		
############################## RECORD RESULTS  ##############################

		results.array[w,10,m] <- p.maxratio
		results.array[w,11,m] <- p.maxcor
		results.array[w,12,m] <- p.peak
		results.array[w,13,m] <- p.mcor_resvec[1]
		results.array[w,14,m] <- p.mcor_resvec[nsim]
		
	}
}

file.name <- paste("Output", version, "R.data", sep = ".")

save.image(file = file.name)


