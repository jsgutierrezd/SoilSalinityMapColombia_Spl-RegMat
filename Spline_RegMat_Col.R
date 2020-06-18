setwd("~/SoilSalinityMapColombia_Spl-RegMat")
library(vip) #variable importance plot
library(pdp) #partial dependence plot
library(cowplot) #plots 
library(plyr) 
library(readxl)
library(raster)
library(rgdal)
library(sp)
library(caret)
library(dplyr)
library(doParallel)
library(ranger)
library(GSIF)
library(aqp)
library(lime)
library(caret)#neural network
library(reticulate)
library(magrittr) #pipe operator
library(randomForest)

### Profiles and sites datasets loading 

prof <- read_excel('G:\\My Drive\\IGAC_2020\\SALINIDAD\\INSUMOS\\BASES\\BASE_NAL_2020.xlsx',sheet = "PROFILES") %>% data.frame
names(prof)
prof <- prof[,c(1,2,3,4,5,10,11,14,15,27)]
sites <- read_excel('G:\\My Drive\\IGAC_2020\\SALINIDAD\\INSUMOS\\BASES\\BASE_NAL_2020.xlsx',sheet = "SITES") %>% data.frame

### Soil profile collection
dat_aqp <- prof
depths(dat_aqp) <- user.profile.id ~ top + bottom
site(dat_aqp) <- sites
coordinates(dat_aqp) <- ~ longitude + latitude
str(dat_aqp)
dat_aqp@sp@proj4string
dat_aqp@sp@proj4string <- CRS("+init=epsg:4326")
dat_aqp@sp@proj4string
class(dat_aqp)
coordinates(sites) <- ~ longitude + latitude
plot(sites)


### pH and ESP vertical distribution
x11()
agg <- slab(dat_aqp, fm= ~ phaq1+PSI)
xyplot(top ~ p.q50 | variable, data=agg, ylab='Profundidad (cm)',
       xlab='Mediana delimitada por los percentiles 25 y 75',
       lower=agg$p.q25, upper=agg$p.q75, ylim=c(190,-2),
       panel=panel.depth_function,
       alpha=0.25, sync.colors=TRUE,
       par.settings=list(superpose.line=list(col='RoyalBlue', lwd=2)),
       prepanel=prepanel.depth_function,
       cf=agg$contributing_fraction, cf.col='black', cf.interval=10, 
       layout=c(3,1), strip=strip.custom(bg=grey(0.8)),
       scales=list(x=list(tick.number=4, alternating=1, relation='free'))
)

dev.off()
dev.off()
dev.off()
```

### Spline pH and ESP

try(pH.0_30 <- mpspline(dat_aqp, 'phaq1', d = t(c(0,30))))
try(ESP.0_30 <- mpspline(dat_aqp, 'PSI', d = t(c(0,30))))
try(pH.30_100 <- mpspline(dat_aqp, 'phaq1', d = t(c(30,100))))
try(ESP.30_100 <- mpspline(dat_aqp, 'PSI', d = t(c(30,100))))

dat_subset <- data.frame(PERFIL = dat_aqp@site$user.profile.id,
                         latitude = dat_aqp@sp@coords[,2],
                         longitude = dat_aqp@sp@coords[,1],
                         pH.0_30 = pH.0_30$var.std[,1],
                         pH.30_100 = pH.30_100$var.std[,1],
                         ESP.0_30= ESP.0_30$var.std[,1],
                         ESP.30_100 = ESP.30_100$var.std[,1])


### Covariates loading and extraction by points

cov <-stack("G:\\My Drive\\IGAC_2020\\SALINIDAD\\INSUMOS\\COVARIABLES\\COV_SSMAP\\CovSSMAP_16062020.tif")
names(cov) <- readRDS("G:\\My Drive\\IGAC_2020\\SALINIDAD\\INSUMOS\\COVARIABLES\\COV_SSMAP\\NamesCovSSMAP_16062020.rds")
names(cov)
dat_subset_sp <- dat_subset
coordinates(dat_subset_sp) <- ~ longitude + latitude
str(dat_subset_sp)
dat_subset_sp
proy <- CRS("+proj=longlat +datum=WGS84")
dat_subset_sp@proj4string <- proy
dat_subset_sp <- spTransform (dat_subset_sp, CRS=projection(cov))



start <- Sys.time()
#dat <- extract(COV84, dat_subset_sp, sp = TRUE)
dat_subset <- cbind(dat_subset, extract(cov, dat_subset_sp))
write.csv(dat_subset, 'G:\\My Drive\\IGAC_2020\\SALINIDAD\\R_CODES\\PRUEBAPILOTOCVC\\RegMatrix_COL.csv')
(nzv <- caret::nearZeroVar(dat_subset[,8:50], saveMetrics = TRUE))
print(Sys.time() - start)


### Removing columns without data or zero(binary variables)
```{r}
names(dat_subset)
dat_subset1 <- dat_subset[,52:155]
dat_subset1 <- dat_subset1[, colSums(dat_subset1,na.rm = T) != 0]
dat_subset <- data.frame(dat_subset[,1:50],dat_subset1)
write.csv(dat_subset, 'G:\\My Drive\\IGAC_2020\\SALINIDAD\\R_CODES\\PRUEBAPILOTOCVC\\RegMatrix_VF_COL.csv',row.names = F)
```