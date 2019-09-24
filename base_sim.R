## 3D encounter simulation
## Author: Juan S. Vargas
## Date: September 2019

# This is a simulation of particles moving in 3D space, and being detected by conical detectors. 
# This is used to test the 3D REM method of calculating density from encounter data.

# Clear workspace
rm(list = ls())

# Load packages
library(tidyverse)

# The first step is to set up the volume where individuals will be distributed. This is a cube.
# Define min and max coordinates. Since it is a cube all sides are of equal length
xmin = -5; xmax = 5;
# Calculate the total volume
totalVol = (xmax-xmin)^3

# The detectors are distributed within a smaller subsection of the total area
dxmin = -2; dxmax = 2;
sampVol = (dxmax-dxmin)^3

# Now we define the detection region of the detectors. This is defined by a detection distance and an opening angle
detRadius = 0.5
detAngle = pi/4
# The REM method is based around the determination of a profile area, which is the area being sampled by a detector instantaneously.
profArea = pi*detRadius^2*(2-2*cos(detAngle/2)+sin(detAngle/2))/4

# Movement parameters. For the simple ideal gas model this is only the total distance that will 
# be covered by each individual during the simulation
movDist = 10*(xmax-xmin)
stepSize = detRadius/10
nSteps = movDist/stepSize

# Simulation parameters. 
nIter = 20 # The number of iterations for each scenario

# Set the total number of detectors
detectorCases = c(5, 10, 20)

# Set the real density that we are trying to estimate
realDens = 1

# Determine the number of individuals by multiplying the density by the total volume
nInd = realDens*totalVol

# Create a data frame to be used for plotting ggplot. This df should have a column for the time, 
# one for the effort (time times number of detectors) and one for the density estimation. Rows will be added after every iteration
dfrows = nSteps*nDetectors*nIter
df <- data.frame(iteration = integer(), 
                 step = integer(), 
                 estdens = numeric(),
                 realdens = numeric(), 
                 detectors = integer())

for (nDetectors in detectorCases) {
  # Open iterations loop
  for (iter in 1:nIter) {
    # Distribute detectors randomly inside the small sampling volume
    detX = runif(nDetectors, dxmin, dxmax)
    detY = runif(nDetectors, dxmin, dxmax)
    detZ = runif(nDetectors, dxmin, dxmax)
    
    # Set the direction that each detector is facing randomly. This direction is defined by two angles for each detector.
    detDir1 = runif(nDetectors, 0, 2*pi) # azimuth angle drawn between 0 and 2 Pi
    detDir2 = acos(runif(nDetectors, -1, 1)) # colatitude angle is determined by drawing the cos of the angle from [-1:1]
    
    # We use the detection distance and the direction that each detector is facing to calculate a vector that defines the axis of the detection cone.
    # This vector is defined by three coordinates, x, y and z
    axX = detRadius*sqrt(1-detDir2^2)*cos(detDir1)
    axY = detRadius*sqrt(1-detDir2^2)*sin(detDir1)
    axZ = detRadius*detDir2
    
    # Distribute all individuals randomly in the larger cube
    indXi = runif(nInd, xmin, xmax)
    indYi = runif(nInd, xmin, xmax)
    indZi = runif(nInd, xmin, xmax)
    
    # Define a direction for each individual. This direction is constant throughout the simulation for the ideal gas model. The direction is defined by two angles, which are drawn at random from a uniform distribution
    indDir1 = runif(nInd, 0, 2*pi)
    indDir2 = runif(nInd, 0, pi)
    
    # Find the final position of each individual. From an initial position and with a given direction and movement distance.
    indXf = indXi + movDist*cos(indDir2)*cos(indDir1)
    indYf = indYi + movDist*cos(indDir2)*sin(indDir1)
    indZf = indZi + movDist*sin(indDir2)
    
    # Create an object to store the number of detections for each detector at each step.
    detectionCount = matrix(0, nDetectors, nSteps)
    # Create an object to store the number of individuals who are actually inside the sampling volume at each step
    iterCounter <- integer(nSteps)
    # Open loop for individuals
    for (ind in 1:nInd) {
      # To find if an individual has crossed the detection volume, we divide each individual trajectory into many points, and determine if each point is within the detection volume. This is determined based on the distance to the detector and the angle with respect to the axis
      # Divide each individual trajectory into many points. The number of points should be high enough that even small distances within the detection volume will be captured.
      trajX = seq(from = indXi[ind], to = indXf[ind], length.out = nSteps)
      trajY = seq(indYi[ind], indYf[ind], length.out = nSteps)
      trajZ = seq(indZi[ind], indZf[ind], length.out = nSteps)
      
      # Determine if the invidual is within the sampling volume
      inX <- trajX > dxmin & trajX < dxmax
      inY <- trajY > dxmin & trajY < dxmax
      inZ <- trajZ > dxmin & trajZ < dxmax
      inReg <- inX & inY & inZ
      # Report it in a single object per iteration for all inviduals
      iterCounter <- iterCounter+inReg
      
      # We need to calculate the difference betwen coordinates of each point and the detector. This will be used to calculate the dot product as well as the Euclidean distance . THis will create a matrix of dimensions (number of points on the trajectory, number of detectors)
      vecDiffX = outer(detX, trajX, '-')
      vecDiffY = outer(detY, trajY, '-')
      vecDiffZ = outer(detZ, trajZ, '-')
      
      # First we calculate the Euclidean distance between each point on the trajectory and the detector
      trajDist = sqrt(vecDiffX^2+vecDiffY^2+vecDiffZ^2)
      
      # From this we can determine if the point is close enought to be detected
      isClose = trajDist <= detRadius
      
      # Now we calculate the dot product between two vectors: the detection zone axis and the segment from the detector to points on the trajectory of the individual. THis gives a single matrix with dims n points on trajectory, n detectors
      dotProd = vecDiffX*axX + vecDiffY*axY + vecDiffZ*axZ
      # Alternative, using subindexing inside a different loop
      # dotProd = (trajX-detX(detector))*axX(detector)+
      #   (trajY-detY[detector])*axY[detector]+
      #   (trajZ-detZ[detector])*axZY[detector]
      
      # The next step is to calculate the product of the lengths of both vectors. For every case the axis measures the detection radius
      nProd = trajDist*detRadius
      
      # We determine if any point is visible from the detector.
      # THe point is visible if the angle between the vectors is less than half of the opening angle of the detector. This angle is calculated as the arccos of the quotient of the dot product and the product of the vectors
      vecAngles = acos(dotProd/nProd)
      isVisible = vecAngles<=detAngle/2
      
      # Now we combine the two logical matrices to determine, for each detector if a point in the trajectory of the individual is contained within the detection region
      detected = isVisible & isClose
      # We transform the detected matrix so that for every trajectory, every entry after the first detection is equal to 1 (already detected)
      for (i in 1:nrow(detected)) {
        detected[i,] = cummax(detected[i,])
      }
      # Add the transformed matrix to the detectionCount output. This object stores the result of a single iteration. It shows for every detector the number of detections at every time step, and therefore can be used to plot the density/accuracy as a funciton of effort (time*nDetectors)
      detectionCount = detectionCount + detected
    } # Close individuals loop
    
    # Calculate real density at every step. This gives a vector of length nSteps
    stepDens <- iterCounter/sampVol
    
    # Calculate density based on count frequency
    estDens = detectionCount/(profArea*stepSize)
    # Recalculate by including the number of steps up to that point
    stepCount = matrix(rep(1:nSteps, each = nDetectors),nrow = nDetectors)
    estDens = estDens/stepCount
    
    # Reorganize data into columns for data frame
    dfi <- as.data.frame(x = estDens)
    
    # Use tidyr's gather function to split the values into columns
    dfi <- gather(data = dfi,
                  key = "step",
                  value = "estdens")
    dfi$step <- rep(1:nSteps, each = nDetectors)
    dfi$realdens <- rep(stepDens, each = nDetectors)
    dfi$iteration <- iter
    dfi$detectors <- nDetectors
    dfi$effort <- dfi$step*nDetectors
    
    # Add observations from this iteration to the larger data frame
    df <- bind_rows(df, dfi)
    
    
  } # Close iterations loop
  
  df$error = (df$estdens-df$realdens)/df$realdens # This is the relative error for every point. Dividing this by the number of observations gives the scaled mean error, a measure of bias (per Walther)
  
  
} # Close detectors loop

finaldf <- df[df$step==nSteps,]
finalmean = finaldf %>% group_by(detectors, iteration) %>% summarise(mean = mean(estdens), error = mean(error), cv = sd(estdens)/mean(estdens))

meanDensdf = df %>% group_by(detectors, iteration, step, effort) %>% summarise(mean = mean(estdens), error = mean(error), cv = sd(estdens)/mean(estdens))

ggplot(meanDensdf, aes(x = effort, y = error)) +  
  # geom_point(aes(color = detectors), alpha=0.5) +
  # geom_point(alpha = 0.5) +
  # geom_line(aes(group = iteration, color = iteration)) +
  geom_smooth(aes(group = detectors, color = as.factor(detectors))) +
  labs(x = 'Effort', y = 'Relative error', color='Detectors') + 
  theme_bw(base_size = 18) +
  theme(legend.position = 'right') 
# facet_grid(.~detectors, scales = 'free')

ggplot(meanDensdf, aes(x = effort, y = error)) +  
  geom_point(aes(color = detectors), alpha=0.5) +
  labs(x = 'Effort', y = 'Relative error') + 
  theme_bw(base_size = 18)
# facet_grid(.~detectors)

ggplot(meanDensdf, aes(x = effort, y = cv)) +  
  geom_smooth(aes(color = as.factor(detectors))) + 
  scale_y_continuous(limits = c(0,2)) +
  labs(x = 'Effort', y = 'Variability (C.V.)', color = "Detectors") + 
  theme_bw(base_size = 18)

a=finalmean %>% group_by(detectors) %>% summarise(meanErr = mean(error), se = sd(error)/sqrt(nIter))
ggplot(a, aes(as.factor(detectors), meanErr, ymin = meanErr - se, ymax = meanErr+se))+
  geom_pointrange(size = 1)+
  labs(x = 'Number of detectors', y = 'Relative error') +
  theme_bw(base_size = 18)


