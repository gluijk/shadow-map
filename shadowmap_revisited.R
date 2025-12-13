# C++ fast shadowmap implementation
# www.overfitting.net
# https://www.overfitting.net/2024/10/radiografia-de-tenerife-con-r.html

library(terra)  # build blur and resample functions
library(tiff)  # save 16-bit TIFF's
library(Rcpp)


# Generic array resample function
# works both for matrix (grayscale images) or 3-channel arrays (colour images)
arrayresample=function(img, DIMX, DIMY, method='bilinear') {
    require(terra)
    
    raster=rast(img)
    rasterrs=rast(nrows=DIMY, ncols=DIMX, extent=ext(raster))
    rasterrs=resample(raster, rasterrs, method=method)
    
    if (is.matrix(img)) return (matrix(as.array(rasterrs), nrow=nrow(rasterrs)))
    else return (as.array(rasterrs))  # convert back to matrix/array
}


# Hillshade calculation
hillshademap=function(DEM, dx=25, dlight=c(0, 2, 3), gamma=1) {
    # hillshademap() inputs DEM data and outputs a hillshade matrix
    #
    # DEM: digital elevation map matrix
    # dx: DEM resolution/cell size (same units as elevation values)
    # dlight: lighting direction. It can be defined in three ways:
    #   a) 1D value indicating light source Azimuth in degrees (0-360)
    #      (0=North, 90=East, 180=South, 270=West)
    #   b) 2D vector indicating light source (X,Y) coordinates
    #   c) 3D vector indicating light source (X,Y,Z) coordinates:
    #      (X=South, Y=East, Z=Up)
    #      dlight=c(0, 2, 3)  # sunrise
    #      dlight=c(0, 0, 1)  # midday
    #      dlight=c(0,-2, 3)  # sunset
    #   NOTE: both in a) and b) a 45º Elevation angle is applied
    # gamma: optional output gamma lift
    
    DIMY=nrow(DEM)
    DIMX=ncol(DEM)
    # If array turn its first dimension into a genuine matrix
    if (!is.matrix(DEM)) {
        print("WARNING: input DEM is not a matrix but an array. First dimension is used")
        DEM=matrix(DEM[,,1], nrow=DIMY, ncol=DIMX)
    }
    
    # Deal with lighting direction
    if (length(dlight)==1) dlight=c(-cos(dlight*pi/180), sin(dlight*pi/180))
    if (length(dlight)==2) dlight=c(dlight, (dlight[1]^2+dlight[2]^2)^0.5)
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


# Shadow projection calculation
shadowmap=function(DEM, dx=25, dlight=c(0, 2, 3)) {
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
    
    shadows=DEM*0  # no shadows at the beginning
    for (y in 1:(DIMX-1)) {  # col DIMX excluded since it cannot get shadow from anyone
        LEN=DIMX-y+1
        dh=(LEN-1)*dx*dlightZ/dlightY
        for (x in 1:DIMY) {
            LIGHTPATH=seq(from=DEM[x,y], to=DEM[x,y]+dh, length.out=LEN)
            DELTA=DEM[x, y:DIMX]-LIGHTPATH
            
            # 4 shadow styles:
            if (max(DELTA)>0) shadows[x,y]=1  # style='hard': shadow=1 if some elevation protrudes above light path
            # shadows[x,y]=shadows[x,y] + length(DELTA[DELTA>0])  # style='width': thickness of protrusion
            # shadows[x,y]=shadows[x,y] + max(0,DELTA[DELTA>0])  # style='max': highest protrusion
            # shadows[x,y]=shadows[x,y] + sum(DELTA[DELTA>0])  # style='mixed': width + height of protrusion
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
    
    return(1-shadows/max(shadows))  # normalize 0..1, shadow=black
}


cppFunction('
    NumericMatrix shadowmap_cpp(NumericMatrix DEM, double dx,
                                NumericVector dlight) {
    
        // ---- VALIDACIÓN ----
        if (dlight[2] <= 0) {
            stop("ERROR: Z light source must be positive");
        }
        if (dlight[0] != 0 && dlight[1] != 0) {
            stop("ERROR: X or Y light source must be 0");
        }
    
        int DIMY = DEM.nrow();
        int DIMX = DEM.ncol();
    
        double dX = dlight[0];
        double dY = dlight[1];
        double dZ = dlight[2];
    
        // ---- TRANSFORMACIÓN A CASO ESTÁNDAR (SOL DESDE EL ESTE) ----
        bool flip_cols = false;
        bool transpose  = false;
    
        double dlightY;
        double dlightZ = dZ;
    
        if (dX == 0 && dY == 0) {  
            // Luz cenital → sin sombras
            NumericMatrix out(DIMY, DIMX);
            for (int i = 0; i < DIMY; i++)
                for (int j = 0; j < DIMX; j++)
                    out(i,j) = 1.0;
            return out;
        } 
        else if (dY < 0) {          // oeste
            dlightY = -dY;
            flip_cols = true;
        } 
        else if (dX > 0) {          // sur
            dlightY = dX;
            transpose = true;
        } 
        else if (dX < 0) {          // norte
            dlightY = -dX;
            transpose = true;
            flip_cols = true;
        } 
        else {                      // este
            dlightY = dY;
        }
    
        // aplicar transformaciones al DEM
        if (transpose) {
            NumericMatrix tmp(DIMX, DIMY);
            for (int i = 0; i < DIMY; i++)
                for (int j = 0; j < DIMX; j++)
                    tmp(j,i) = DEM(i,j);
            DEM = tmp;
            int aux = DIMY; DIMY = DIMX; DIMX = aux;
        }
    
        if (flip_cols) {
            NumericMatrix tmp(DIMY, DIMX);
            for (int i = 0; i < DIMY; i++)
                for (int j = 0; j < DIMX; j++)
                    tmp(i, DIMX-1-j) = DEM(i,j);
            DEM = tmp;
        }
    
        // ---- CÁLCULO DE SOMBRAS ----
        NumericMatrix shadows(DIMY, DIMX);
        double slope = dx * dlightZ / dlightY;
    
        for (int y = 0; y < DIMX - 1; y++) {
    
            int LEN = DIMX - y;
    
            for (int x = 0; x < DIMY; x++) {
    
                double start = DEM(x, y);
                bool shadow = false;
    
                for (int k = 1; k < LEN; k++) {
                    double lightH = start + slope * k;
                    double terrain = DEM(x, y + k);
    
                    if (terrain > lightH) {
                        shadow = true;
                        break;
                    }
                }
    
                if (shadow)
                    shadows(x, y) = 1.0;
            }
        }
    
        // ---- DESHACER TRANSFORMACIONES ----
        if (flip_cols) {
            NumericMatrix tmp(DIMY, DIMX);
            for (int i = 0; i < DIMY; i++)
                for (int j = 0; j < DIMX; j++)
                    tmp(i, DIMX-1-j) = shadows(i,j);
            shadows = tmp;
        }
    
        if (transpose) {
            NumericMatrix tmp(DIMX, DIMY);
            for (int i = 0; i < DIMY; i++)
                for (int j = 0; j < DIMX; j++)
                    tmp(j,i) = shadows(i,j);
            shadows = tmp;
        }
    
        // ---- NORMALIZACIÓN FINAL ----
        int SY = shadows.nrow();
        int SX = shadows.ncol();
        double m = 0;
    
        for (int i = 0; i < SY; i++)
            for (int j = 0; j < SX; j++)
                if (shadows(i,j) > m) m = shadows(i,j);
    
        NumericMatrix out(SY, SX);
    
        if (m == 0) {
            for (int i = 0; i < SY; i++)
                for (int j = 0; j < SX; j++)
                    out(i,j) = 1.0;
            return out;
        }
    
        for (int i = 0; i < SY; i++)
            for (int j = 0; j < SX; j++)
                out(i,j) = 1.0 - shadows(i,j) / m;
    
        return out;
    }
')


#################################################

# 1. READ TENERIFE RASTER DATA

RESOLUTION=25
MAXIMO=3710.062  # Teide
DEM=readTIFF("tenerifecomposite.tif")*MAXIMO
hist(DEM[DEM>0], breaks=800)


#################################################

# 2. RESCALE DEM TO FULL HD (1920 X 1080)

f=1
fscale=0.347758887171561*f  # downsampling scaling (chosen to fit well in Full HD)
RESOLUTION=RESOLUTION/fscale  # downsampling also reduces RESOLUTION
DIMY=round(nrow(DEM)*fscale)
DIMX=round(ncol(DEM)*fscale)
DEMrs=arrayresample(DEM, DIMX, DIMY)
# DEMrs[DEMrs<0]=0  # clip to 0 Lanczos artifactsn (bilinear doesn't have)

DIMYfullHD=1080*f
DIMXfullHD=1920*f
DEM=matrix(0, nrow=DIMYfullHD, ncol=DIMXfullHD)
DEM[(DIMYfullHD/2-DIMY/2):(DIMYfullHD/2+DIMY/2-1),
    (DIMXfullHD/2-DIMX/2):(DIMXfullHD/2+DIMX/2-1)]=DEMrs
rm(DEMrs)

writeTIFF(DEM/max(DEM), "tenerifecomposite_fullHD.tif",
          bits.per.sample=16, compression="LZW")


#################################################

# 3. GENERATE STANDARD HILLSHADE

hillshade=hillshademap(DEM, dx=RESOLUTION, dlight=c(0, 1, 1))

# Save hillshade
writeTIFF(hillshade, "hillshade.tif", bits.per.sample=16, compression="LZW")

# Display hillshade
image(t(hillshade[nrow(hillshade):1,]), useRaster=TRUE,
      col=c(gray.colors(256, start=0, end=1, gamma=1)),
      asp=nrow(hillshade)/ncol(hillshade), axes=FALSE)


#################################################

# 4. CALCULATE SHADOWS AND APPLY THEM TO HILLSHADE

a=system.time(shadowmap(DEM, dx=RESOLUTION, dlight=c(0, 30, 5)))[3]  # 55.50s
b=system.time(shadowmap_cpp(DEM, dx=RESOLUTION, dlight=c(0, 30, 5)))[3]  # 2.21s
print(paste0("Improvement from ", round(a,2), "s to ", round(b,2), "s: ",
             round(log(a/b, 10), 2), " orders of magnitude"))

i=1
for (z in seq(from=1, to=10, length.out=40)) {
    name=paste0("hillshade_shadows_", ifelse(i<10,"0",""), i,".tif")
    shadows=shadowmap_cpp(DEM, dx=RESOLUTION, dlight=c(0, 30, z))  
    shadowsfinal=shadows
    shadowsfinal[shadowsfinal==0]=0.65  # shadow blend factor
    writeTIFF(hillshade*shadowsfinal, name, bits.per.sample=16, compression="LZW")
    i=i+1
}
