# Collinearity-equations-in-photogrametry
Collinearity equations to map a known 3D location to the image plane using internal and external orientation

###  The base R code for the photogrammetry lab exercise.
###
rm( list = ls() )

setwd("C:/Users/Olamidayo Fasalejo/OneDrive - University of Eastern Finland/Documents/Remote sensing/Lab 3")
?matrix
###
### External orientation of Lvl02-01385-Col.jpg.
###

# Camera's location (KKJ3 PCS, N60 elevation).
X0 <- 3597426.93     # Camera location X (meters)
Y0 <- 7022404.66     # Camera location Y (meters)
Z0 <- 3215.73        # Camera location Z (meters)

# Camera's heading (rotations).
O0 <- -0.0184150     # Rotation Omega w (radians)
P0 <- -0.0052709     # Rotation Phi   (radians)
K0 <-  0.6869564     # Rotation Kappa (radians)

###
### Internal orientation of the used digital camera.
###
cc <- 0.1014         # Camera constant (meters)
ps <- 0.000028125    # Pixel size in CCD (meters)
nc <- 3680           # Number of colums in CCD
nr <- 2400           # Number of rows in CCD
x0 <- -0.00018       # X offset of PP (meters; left from the image center)
y0 <- -0.00036       # Y offset of PP (meters; down from the image center)

###
### 3D location to project onto an image plane, in this case a lidar point.
###
X <- 3597544.50      # X location of a point (KKJ3)
Y <- 7022052.80      # Y location of a point (KKJ3)
Z <- 221.64          # Z location of a point (N60)

###
### Rotation matrix M.
###

#------------------------------------------------------------------#
#  DIY: Construct the rotation matrix M.                           #
#------------------------------------------------------------------#

M <- matrix( NA, 3, 3 )  # Rotation matrix
M
###
### Collinearity equations.

M[1,1]<-cos(P0)*cos(K0)

M[1,2]<-cos(O0)*sin(K0)+sin(O0)*sin(P0)*cos(K0)

M[1,3]<-sin(O0)*sin(K0)-cos(O0)*sin(P0)*cos(K0)
  
M[2,1]<--cos(P0)*sin(K0)

M[2,2]<-cos(O0)*cos(K0)-sin(O0)*sin(P0)*sin(K0)

M[2,3]<-sin(O0)*cos(K0)+(cos(O0)*sin(P0)*sin(K0))

M[3,1]<-sin(P0)

M[3,2]<- -sin(O0)*cos(P0)

M[3,3]<- cos(O0)*cos(P0)

M
t(M)-solve(M)



###

#------------------------------------------------------------------#
#  DIY: Compute the image space coordinates (in meters) with       #
#       respect to the image center (IC).                          #
#------------------------------------------------------------------#

xx0 <- -cc*(M[1,1]*(X-X0)+M[1,2]*(Y-Y0)+M[1,3]*(Z-Z0))/(M[3,1]*(X-X0)+M[3,2]*(Y-Y0)+M[3,3]*(Z-Z0))        # x coordinate with respect to the IC

                  
yy0 <- -cc*(M[2,1]*(X-X0)+M[2,2]*(Y-Y0)+M[2,3]*(Z-Z0))/(M[3,1]*(X-X0)+M[3,2]*(Y-Y0)+M[3,3]*(Z-Z0))     # y coordinate with respect to the IC

#------------------------------------------------------------------#
#  DIY: Then compute the image space coordinates (in meters) with  #
#       respect to the principal point (PP).                       #
#------------------------------------------------------------------#

x <- xx0+x0              # x coordinate with respect to the PP
y <- yy0+y0              # y coordinate with respect to the PP

###
### Pixel space coordinates.
###

#------------------------------------------------------------------#
#  DIY: Finally compute pixel space coordinates with respect to    #
#       the top left corner of the top left pixel.                 #
#------------------------------------------------------------------#

xpx <- nc/2+x/ps-0.5            # x coordinate in pixel space
ypx <- -(nr/2-y/ps+0.5)            # y coordinate in pixel space

###
### Print the result. Where is it located in the Lvl02-01385-Col.jpg?
###

cat( "\nPixel space coordinates:\nX = ", xpx, "\nY = ", ypx, "\n" )
