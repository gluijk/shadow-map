# Calculation of shadows projection from a digital elevation map (DEM)
# www.overfitting.net
# https://www.overfitting.net/2024/04/proyeccion-de-sombras-sobre-un-dem-con-r.html

library(data.table)  # fread()
library(terra)  # read GeoTIFF, plot, reprojection, crop and resample
library(tiff)  # save 16-bit TIFF's
library(png)  # save 8-bit PNG's


hillshademap=function(DEM, dx=25, dlight=c(0, 2, 3), gamma=1) {
    # hillshademap() inputs DEM data and outputs a hillshade matrix
    #
    # DEM: digital elevation map matrix
    # dx: DEM resolution/cell size (same units as elevation values)
    # dlight: lighting direction (3D vector defined from observer to light source):
    #   dlight=c(0, 2, 3)  # sunrise
    #   dlight=c(0, 0, 1)  # midday
    #   dlight=c(0,-2, 3)  # sunset
    # gamma: optional output gamma lift
    
    DIMY=nrow(DEM)
    DIMX=ncol(DEM)
    # If array turn its first dimension into a genuine matrix
    if (!is.matrix(DEM)) {
        print("WARNING: input DEM is not a matrix but an array. First dimension is used")
        DEM=matrix(DEM[,,1], nrow=DIMY, ncol=DIMX)
    }
    
    dlightM=sum(dlight^2)^0.5
    
    # Vectorial product to calculate n (orthogonal vector)
    nx = 2*dx*(DEM[1:(DIMY-2), 2:(DIMX-1)] - DEM[3:DIMY,     2:(DIMX-1)])
    ny = 2*dx*(DEM[2:(DIMY-1), 1:(DIMX-2)] - DEM[2:(DIMY-1), 3:DIMX])
    nz = 4*dx^2
    nM = (nx^2 + ny^2 + nz^2)^0.5
    
    # Dot product to calculate cos(theta)
    dn = dlight[1]*nx + dlight[2]*ny + dlight[3]*nz  # (DIMY-2)x(DIMX-2) matrix
    
    # Reflectance (=cos(theta))
    hillshadepre=dn/(dlightM*nM)
    hillshadepre[hillshadepre<0]=0  # clip negative values
    
    # Add 1-pix 'lost' borders
    hillshademap=matrix(0, nrow=DIMY, ncol=DIMX)
    hillshademap[2:(DIMY-1), 2:(DIMX-1)]=hillshadepre
    rm(hillshadepre)
    hillshademap[c(1,DIMY),]=hillshademap[c(2,DIMY-1),]
    hillshademap[,c(1,DIMX)]=hillshademap[,c(2,DIMX-1)]
    
    return(hillshademap^(1/gamma))
}


shadowsmap=function(DEM, dx=25, dlight=c(0, 2, 3)) {
    # shadowsmap() inputs DEM data and outputs shadows projection
    #
    # DEM: digital elevation map matrix
    # dx: DEM resolution/cell size (same units as elevation values)
    # dlight: lighting direction (3D vector defined from observer to light source):
    #   dlight=c(0, 2, 3)  # sunrise
    #   dlight=c(0, 0, 1)  # midday
    #   dlight=c(0,-2, 3)  # sunset
    
    # Loop version: 2.22min (FASTER than the vectorized version)
    
    # Check for invalid input parameters
    if (dlight[3]<=0) {
        print(paste0("ERROR: Z light source must be positive (Z=", dlight[3], ")"))
        return (-1)
    } else if (dlight[1] & dlight[2]) {
        print(paste0("ERROR: X or Y light source must be 0 (X=", dlight[1], ", Y=", dlight[2],  ")"))
        return (-1)
    }
    
    DIMY=nrow(DEM)
    DIMX=ncol(DEM)
    # If array turn its first dimension into a genuine matrix
    if (!is.matrix(DEM)) {
        print("WARNING: input DEM is not a matrix but an array. First dimension is used")
        DEM=matrix(DEM[,,1], nrow=DIMY, ncol=DIMX)
    }
    
    # Turn any light source into East light source setting dlightY
    if (!dlight[1] & !dlight[2]) {
        print(paste0("WARNING: zenithal Z light source, no shadows (X=", dlight[1], ", Y=", dlight[2], ")"))
        return (DEM*0+1)  # return white matrix
    } else if (dlight[2]<0) {  # West light source
        dlightY=-dlight[2]
        DEM=DEM[,ncol(DEM):1]  # transpose cols
    } else if (dlight[1]>0) {  # South light source
        dlightY=dlight[1]
        DEM=t(DEM)  # transpose all
    } else if (dlight[1]<0) {  # North light source
        dlightY=-dlight[1]
        DEM=t(DEM)  # transpose all
        DEM=DEM[,ncol(DEM):1]  # transpose cols
    } else {  # East light source (standard case)
        dlightY=dlight[2]
    }
    dlightZ=dlight[3]
    
    # DIMY and DIMX change if South/North light source
    DIMY=nrow(DEM)
    DIMX=ncol(DEM)
    
    shadows=DEM*0  # 0=no shadow, 1=shadow (no shadow at the beginning)
    for (y in 1:(DIMX-1)) {  # col DIMX excluded since it cannot get shadow from anyone
        LEN=DIMX-y+1
        dh=(LEN-1)*dx*dlightZ/dlightY
        for (x in 1:DIMY) {
            LIGHTPATH=seq(from=DEM[x,y], to=DEM[x,y]+dh, length.out=LEN)
            DELTA=DEM[x, y:DIMX]-LIGHTPATH
            if (max(DELTA)>0) shadows[x,y]=1  # shadow if some elevation protrudes above light path
        }
    }
    
    # Undo transformations if applied
    if (dlight[2]<0) {  # West light source
        shadows=shadows[,ncol(shadows):1]  # transpose cols
    } else if (dlight[1]>0) {  # South light source
        shadows=t(shadows)  # transpose all
    } else if (dlight[1]<0) {  # North light source
        shadows=shadows[,ncol(shadows):1]  # transpose cols
        shadows=t(shadows)  # transpose all
    }
    
    return(1-shadows)  # shadow=black
    
    
    # Vectorized version: 3.36min (SLOWER than the loop version)
    
    # First we vectorize the needed linear interpolation
    # seq.vectorized=Vectorize(seq.default, vectorize.args=c("from", "to"))
    # 
    # shadows=DEM*0  # 0=no shadow, 1=shadow (no shadow at the beginning)
    # for (y in 1:(DIMX-1)) {  # col DIMX excluded since it cannot get shadow from anyone
    #     LEN=DIMX-y+1
    #     dh=(LEN-1)*dx*dlightZ/dlightY
    #     LIGHTPATH=t(seq.vectorized(from=DEM[,y], to=DEM[,y]+dh, length.out=LEN))
    #     DELTA=DEM[, y:DIMX]-LIGHTPATH
    #     shadows[,y]=apply(DELTA, MARGIN=1, FUN=max)  # max of each row
    # }
    # shadows[shadows<=0]=0
    # shadows[shadows>0]=1
}


#################################################

# 1. READ RASTER DATA IN 4 TXT FILES OF RESOLUTION 25m

sierra_11=data.matrix(
    fread("PNOA_MDT25_ETRS89_HU30_0483_LID.txt", sep=" ", dec="."))
sierra_12=data.matrix(
    fread("PNOA_MDT25_ETRS89_HU30_0484_LID.txt", sep=" ", dec="."))
sierra_21=data.matrix(
    fread("PNOA_MDT25_ETRS89_HU30_0508_LID.txt", sep=" ", dec="."))
sierra_22=data.matrix(
    fread("PNOA_MDT25_ETRS89_HU30_0509_LID.txt", sep=" ", dec="."))


# Final crop (manually detected overlappings)
sierra=matrix(0, nrow=1508, ncol=2269)
sierra[1:759, 11:1148]=sierra_11
sierra[14:768, 1136:2269]=sierra_12
sierra[741:1499, 1:1141]=sierra_21
sierra[754:1508, 1129:2265]=sierra_22
sierra=sierra[14:1499, 11:2265]
rm(sierra_11, sierra_12, sierra_21, sierra_22)
dim(sierra)  # -> 1486 x 2255

RESOLUTION=25  # 25m/pixel
DIMY=nrow(sierra)
DIMX=ncol(sierra)
ALTMIN=min(sierra)
ALTMAX=max(sierra)

# Guardamos raster en TIFF
writeTIFF(sierra/ALTMAX, "sierra.tif", bits.per.sample=16, compression="LZW")


#################################################

# 2. GENERATE STANDARD HILLSHADE

hillshade=hillshademap(sierra, dx=RESOLUTION, dlight=c(0, 10, 5))

# Save hillshade
writeTIFF(hillshade, "hillshade.tif", bits.per.sample=16, compression="LZW")

# Display hillshade
image(t(hillshade[nrow(hillshade):1,]), useRaster=TRUE,
      col=c(gray.colors(256, start=0, end=1, gamma=1)),
      asp=nrow(hillshade)/ncol(hillshade), axes=FALSE)


#################################################

# 3. CALCULATE SHADOWS

# Test all light source directions -> OK
shadows=shadowsmap(DEM=sierra, dx=RESOLUTION, dlight=c(0, 40, 5))
writePNG(shadows, "shadowsEAST.png")

shadows=shadowsmap(DEM=sierra, dx=RESOLUTION, dlight=c(0, -40, 5))
writePNG(shadows, "shadowsWEST.png")

shadows=shadowsmap(DEM=sierra, dx=RESOLUTION, dlight=c(40, 0, 5))
writePNG(shadows, "shadowsSOUTH.png")

shadows=shadowsmap(DEM=sierra, dx=RESOLUTION, dlight=c(-40, 0, 5))
writePNG(shadows, "shadowsNORTH.png")


# Example c(0, 30, 5)
shadows=shadowsmap(DEM=sierra, dx=RESOLUTION, dlight=c(0, 30, 5))
writePNG(shadows, "shadows_5_30.png")

# Arbitrary angle rotationn with (0, 40, 5)
sierra_rotated=readTIFF("sierra_rotated30.tif")
sierra_rotated=sierra_rotated*ALTMAX  # restore height values

shadows=shadowsmap(DEM=sierra_rotated, dx=RESOLUTION, dlight=c(0, 40, 5))
writePNG(shadows, "shadows_rotated.png")

hillshade=hillshademap(sierra_rotated, dx=RESOLUTION, dlight=c(0, 10, 5))
writeTIFF(hillshade, "hillshade_rotated.tif", bits.per.sample=16, compression="LZW")
