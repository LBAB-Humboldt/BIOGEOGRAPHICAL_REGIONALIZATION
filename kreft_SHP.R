rm(list=ls(all=T))
gc() ##para borrar todo lo que quede oculto en memoria

memory.limit(size = 1000000) ## para asiganar 20 gigas del disco duro a ram
library("raster")
library("vegan")
library("maptools")
library("clustsig")
library(PBSmapping)
library(epitools)
library(RColorBrewer)
library(classInt)
library(sp)
library("cluster")
library("ape")
library("phangorn")
library(rgdal)
library(svGUI)
library("svDialogs")

# #### 0. DEFINIR RUTAS Y FUNCIONES ---------------------------------------------------


ruta_modelos<-(dlgDir(default = getwd(), title="ESPECIFIQUE LA RUTA DE LOS MODELOS")$res)
#"C:/Users/GIC 9/Documents/GBIF3/MODELOS"
ruta_elegidos<-(dlgDir(default = getwd(), title="ESPECIFIQUE LA RUTA DE MODELOS VERIFICADOS")$res)
ruta_salida<-(dlgDir(default = getwd(), title="ESPECIFIQUE LA RUTA DE LOS RESULTADOS")$res)
#"~/GBIF3/BETA"

### FUNCIONES

grafica_recambio=function(grilla,mascara,distancia2,xy1){
  # ver el area de estudio
  grilla2=mask(grilla,mascara)
  # selccionar celda
  celda=cellFromXY(grilla2,as.numeric(xy1))
  columna=distancia2[,which(colnames(distancia2)==as.character(celda))]
  colum=as.data.frame(cbind(as.numeric(row.names(distancia2)),columna))
  if(ncol(colum)==2){  
    names(colum)=c("id","dist")
    valores=grilla2[values(grilla2)%in%colum[,1]] 
    for (j in valores){
      grilla2[j]=colum[which(colum[,1]==j),2]
    }
    grilla2[which(values(grilla2)>=2)]=NA
    xy=as.data.frame(xy1)
    names(xy)=c("x","y")
    coordinates(xy)=~x+y
    resultados=c(grilla2,xy1) 
    return(resultados)}
} 

grafica_recambio_shp=function(mascara,distancia,xy1){
  # ver el area de estudio
  mascara2=mascara
  mascara2$recambio=as.numeric(NA)
  # selccionar celda
  xy1=SpatialPoints(xy1)
  celda=over(xy1,mascara2)$ID
  columna=distancia[,which(colnames(distancia)==as.character(celda))]
  colum=as.data.frame(cbind(as.numeric(row.names(distancia)),columna))
  if(ncol(colum)==2){  
    names(colum)=c("id","dist")
     valores=colum$id
    for (j in valores){
      mascara2$recambio[mascara2$ID==j]=colum$dist[colum$id==j]
      }
  mascara2$recambio=as.numeric(mascara2$recambio)
  }
  plot.new()
  spplot(mascara2,zcol="recambio")
  plot(xy1, add=T, col="red")
  resultados=c(mascara2,xy1) 
  return(resultados)} 

ESPECIES=function(xy){
  xy=t(as.matrix(xy))
  celda=cellFromXY(grilla2,as.numeric(xy))
  fila=DF2[which(row.names(DF2)==as.character(celda)),]
  sp=fila[,which(colSums(fila)!=0)]
  return(colnames(sp))
}

ESPECIES_shp=function(xy){
  xy=SpatialPoints(xy)
  celda=over(xy,mascara2)$ID
  fila=DF[which(row.names(DF)==as.character(celda)),]
  sp=fila[,which(colSums(fila)!=0)]
  return(colnames(sp))
}

col.br <- colorRampPalette(c("blue", "cyan", "yellow","gray"))


EVA1=function(g){
  clust=g
  alt=rev(sort(clust$height))
  ceros=which(alt==0)
  if (length(ceros)>1){ceros=ceros[-1]}
  alturas=alt[-ceros]
  resultado=NULL
  clust$height=sort(clust$height)
  for (i in 1:length(alturas)){
    corte=cutree(clust,h=alturas[i])
    result=cbind(alturas[i],max(corte))
    resultado=rbind(resultado,result)
  }
  resultado=as.data.frame(resultado)
  colnames(resultado)=c("metrica","grupos")
  return(resultado)
}



EVA2=function(cluster){
  cluster
  cluster$height=sort(cluster$height)
  ngrupos=max(cutree(cluster,h=0))
  DF_E=as.data.frame(as.numeric(celdas))
  row.names(DF_E)=celdas
  DF_E$grupo=NA
  proportion=NULL
  for (w in 2:ngrupos){
    cut=as.data.frame(cutree(cluster,k=w))
    names(cut)=c("group")
    DF_E$grupo[row.names(DF_E)==row.names(cut)]=cut$group[row.names(DF_E)==row.names(cut)]
    niveles=unique(DF_E$grupo)
    resumen=NULL
    for (z in niveles){
      sumfilas=colSums(DF_E[which(DF_E$grupo==z),2:(ncol(DF_E)-1)])
      sumfilas[sumfilas>0]=1
      res=c(z,sumfilas)
      resumen=rbind(resumen,res)
    }
    endemic=length(which(colSums(resumen)==1))
    propor=c(w,endemic/(ncol(resumen)-1))
    proportion=rbind(proportion,propor)
  }
  proportion=as.data.frame(proportion, row.names<-F)
  names(proportion)=c("grupos","metrica")
  return(proportion)
}



grup=function(c,eva){
  mod1=lm(eva$metrica[1:c]~eva$grupos[1:c])
  mod2=lm(eva$metrica[c+1:b]~eva$grupos[c+1:b])
  red1=sqrt(mean((mod1$residuals)^2))
  red2=sqrt(mean((mod2$residuals)^2))
  RMSE=(((c-1)/(b-1))*red1) + (((b-c)/(b-1))*red2)
  return(RMSE)
}

# ### 1. AREA DE ESTUDIO --------------------------------------------------

#MASCARAS # 
cat("seleccione .shp del area de estudio")
mascara=readShapePoly(file.choose()) ### define para que area se va ha realizar el analisis
mascara$ID=1:nrow(mascara)
extent=extent(mascara)

#area de estucio
cat("seleccione raster de area de estudio")
orobioma=raster(file.choose())#"C:/SIG/CLMATCAS COL/CLMATCAS COL/bio1_col.asc") 

factor=10 #derfine nivel de agregacion 
aoi=aggregate(orobioma,factor)
aoi=mask(aoi,mascara)

#crear grilla
grilla=aoi
names(grilla)="grilla"
grilla[1:ncell(grilla)]<-1:ncell(grilla)



#popo=matrix(0,ncell(aoi),ncell(aoi))

# ### 2. ESAMBLE DE MODELOS -----------------------------------------------

setwd(ruta_modelos)
MODELOS<-stack(list.files(pattern="*10p_cut.tif$")) ## para la prueba tome solo los 10 priemros modelos
#ojo esto es olo par ala prueba

nombres=list.files(pattern="*10p_cut.tif$")
nombres=as.data.frame(strsplit(nombres,"_"))
nombres=as.data.frame(t(nombres))
nombres2=paste(nombres$V1,nombres$V2,sep="_")


# cortar modelos para el area de estudio
MODELOS=mask(MODELOS,mascara)

#adjuntar modelos depurados 



mascara2=MODELOS[[1]]
##MAGNOLIAS y YUCA ARROZ
setwd(ruta_elegidos)
elegidos=list.files(pattern=".tif$")
# stack  sobreponer capas solo para raster

for (i in 2:length(elegidos)){
  especie=raster(elegidos[[i]])
  capa=resample(especie,mascara2,method="ngb")
  MODELOS=addLayer(MODELOS,capa)
}


nombres3=names(MODELOS2[[2409:2432]])

nombrestodos=c(nombres2,nombres3)

#### subir shp de mapas especies adicionales en formato raster
## INVEMAR

##UICN 

# #### 3. MATRIZ DE DISTANCIAS --------------------------------------------
setwd(ruta_salida)


#GENERAR TABLA 

### PARA SHP

## NOTA , DF se puede hacer aca con la siguiente linea , el rollo es que demore mucho ,
#si desisten  jorge lo corre 

MODELOS2=mask(MODELOS,mascara)
DF=extract(MODELOS2, mascara,df=T, fun=max)

DF=na.omit(DF)
head(DF)
 

#para reclasificar si es necesario 
# for (i in 2: ncol(DF)){
#   DF[which(DF[,i]!=0),i]=1
# }

#remover filas vacias
vacias=c(which(rowMeans(DF[,2:ncol(DF)])==0)) # solo estamos limpiando vacias , se puede menos de 5 (kreft)
if (length(vacias)==0){ DF=DF} else {DF=DF[-vacias,]}

celdas=DF[,1]
nombres2=c("ID", nombrestodos)
names(DF)=nombres2

#Matriz de distancias para shp
DISTAN=betadiver(DF[,-1], "sim")
distancia=as.matrix(DISTAN)
### PARA RASTER 
# es necesario sacar las dos matrices y las dos distancias

#matriz de distancias para raster
length((getValues(aoi)))
DF2=as.data.frame(matrix(NA,length((getValues(aoi))),1)) # cambiar el numero de filas 
DF2[,1]=getValues(grilla)

for(i in 2:nlayers(MODELOS2)){
  DF2[,i+1]=NA
  print(i)  
  agragado=aggregate(MODELOS2[[i]],factor)
  primera=resample(agragado,aoi)
  VALORES=getValues(primera)
  VALORES[which(VALORES>=0.25)]=1
  VALORES[which(VALORES<0.25)]=0
  VALORES=as.integer(VALORES)
  DF2[,i]=VALORES
}


DF2=DF2[,1:(ncol(DF2)-1)]



DF2=na.omit(DF2[,2:ncol(DF2)])
head(DF2)
celdas=row.names(DF2) 

#remover filas vacias
vacias2=c(which(rowMeans(DF2[,2:ncol(DF2)])==0)) # solo estamos limpiando vacias , se puede menos de 5 (kreft)
if (length(vacias2)==0){ DF2=DF2} else {DF2=DF2[-vacias2,]}

names(DF2)=nombres2


DISTAN2=betadiver(DF2[,2:ncol(DF2)], "sim")
distancia2=as.matrix(DISTAN2)

# ##########. 4 PATRONES DIVERSIDAD ALFA, GAMA, BETA ----------------------


##### ALFA 
setwd(ruta_salida)
alfa=sum(MODELOS2)
writeRaster(alfa,"alfa",overwrite=TRUE, format="GTiff")

#### GAMA
DF.RAS=as.data.frame(MODELOS2)
GAMA=grilla
for(d in 1:nlayers(MODELOS2)){
  print(d)  
  modelo=MODELOS2[[d]]
  agragar=aggregate(modelo,factor)
  CAPA=resample(agragar,aoi)
  GAMA=addLayer(GAMA,CAPA)
}

GAMA=aggregate(alfa,factor,fun=sum)
ALFA=aggregate(alfa,factor,fun=mean)


# BETA MULTIPLICATIVA
BETAm=(GAMA/ALFA)
#BETAm[BETAm>2]=2
plot(BETAm)

## BETA ADITIVA
BETAa=GAMA-ALFA
plot(BETAa)


## BETA EFECTIVE TURNOVER
BETAmt=(GAMA-ALFA)/ALFA
plot(BETAmt)

## BETA EFECTIVE TURNOVER
BETAmg=(GAMA-ALFA)/GAMA
plot(BETAmg)


writeRaster(ALFA,"ALFAmean",overwrite=TRUE)
writeRaster(BETAm,"BETAm.tif",overwrite=TRUE)
writeRaster(BETAa,"BETAa.tif",overwrite=TRUE)
writeRaster(BETAmg,"BETAmg.tif",overwrite=TRUE)


##PARA SACAR PATRONES POR POLIGONOS
# se remuestrea la media de la diversidad beta para los pixeles presentes en cada poligono,
# al intentar hacerlo por poligonos directamente tenemso el problema de como definir  alfa y gama
# de todas formas creo que la mejor  forma de mostrar eso creo que lo mejor es por pixel 

patrones.shp=mascara

BETAa.P=extract(BETAa, mascara,df=T,fun=mean,na.rm=TRUE)
BETAm.P=extract(BETAm, mascara,df=T,fun=mean,na.rm=TRUE)
BETAmg.P=extract(BETAmg, mascara,df=T,fun=mean, na.rm=TRUE)
BETAmt.P=extract(BETAmt, mascara,df=T,fun=mean, na.rm=TRUE)


### para graficar
patrones.shp$aditiva=NA
patrones.shp$multiplicativa=NA
patrones.shp$efective1=NA
patrones.shp$efective2=NA

for (i in patrones.shp$ID){
  #aditiva
  patrones.shp$aditiva[patrones.shp$ID==i]=BETAa.P$layer[BETAa.P$ID==i]
  #multiplicativa
  patrones.shp$multiplicativa[patrones.shp$ID==i]=BETAm.P$layer[BETAm.P$ID==i]
  #efective turnover1
  patrones.shp$efective1[patrones.shp$ID==i]=BETAmg.P$layer[BETAmg.P$ID==i]
  #efective trunover 2
  patrones.shp$efective2[patrones.shp$ID==i]=BETAmt.P$layer[BETAmt.P$ID==i]
  
}


spplot(patrones.shp, zcol="aditiva", main="aditiva")
spplot(patrones.shp, zcol="multiplicativa", main="multiplicativa")
spplot(patrones.shp, zcol="efective1", main="efective1")
spplot(patrones.shp, zcol="efective2", main="efective2")


## Guardar el resultado de los patrones por poligono

writePolyShape(patrones.shp, "patrones.shp", factor2char = TRUE, max_nchar=254)

# #### 5. VISUALIZACION GEOGRAFICA  RECAMBIO  -----------------------------

#selccionar punto de area de trabajo 
# lat = c(extent@ymin,extent@ymax)
# lon=c(extent@xmin,extent@xmax)
# center = c(lat=mean(lat), lon=mean(lon));
# zoom <- min(MaxZoom(range(lat), range(lon)));
#MyMap <- GetMap(center=center, zoom=zoom,maptype="terrain")
#tmp <- PlotOnStaticMap(MyMap, lat = c(40.702147,40.711614,40.718217), lon = c(-74.015794,-74.012318,-73.998284), destfile = "MyTile1.png", cex=1.5,pch=20,col=c('red', 'blue', 'green'), add=FALSE);
#tmp <- PlotPolysOnStaticMap(MyMap, polys=paramo);
#dev.off()


## PARA  SHP
plot.new()
plot(extent,main="SELECCIONE POLIGONO DE INTERES")
plot(mascara,add=T)
cat("selcciones el nuemro de areas para las cuales desea ver su distancia")

n=5
xy=locator(n=n)
#xy=as.data.frame((xy))
dev.off()

#graficas de recambio
for (l in 1:n) {
  print(l)
  xy1=cbind(xy[[1]][l],xy[[2]][l])
  RECAMBIO=grafica_recambio_shp(mascara,distancia,xy1)
  prueba=RECAMBIO[[1]]
  spplot(RECAMBIO[[1]],zcol="recambio") # grafica si lo haces uno por uno pero no se porque no grafica en el for  igual el shp queda guardado
  #plot(RECAMBIO[[2]],col="red")
  writePolyShape(prueba, paste0("dist",l,".shp"), factor2char = TRUE, max_nchar=254)
  
  ##especies presentes en poligono seleccionado 
  sp_en_celda=ESPECIES_shp(xy1)
  print(sp_en_celda)
}




### PARA RASTER

plot.new()
plot(extent,main="SELECCIONE AREA DE INTERES")
plot(grilla2,add=T)
cat("selcciones el nuemro de areas para las cuales desea ver su distancia")

n=5
xy=locator(n=n)
#xy=as.data.frame((xy))
dev.off()

#graficas de recambio
for (l in 1:n) {
  xy1=cbind(xy[[1]][l],xy[[2]][l])
  RECAMBIO=grafica_recambio(grilla,mascara,distancia2,xy1)
  plot(RECAMBIO[[1]],col=c("red",rev(col.br(10))))
  #plot(RECAMBIO[[2]],add=T,col="red")
  writeRaster(RECAMBIO[[1]],paste("dist",l),overwrite=TRUE, format="GTiff")
  
  grilla2=mask(grilla,mascara)
  ##especies presentes en celda seleccionada 
  sp_en_celda=ESPECIES(xy1)
  print(sp_en_celda)
}

# ### 6. ORDENACION -------------------------------------------------------
# CARO NO ESTAMOS HACIENDO ESTO  POR QUE REQUIERE MUCHO COMPUTO 
sol <- metaMDS(DF, distfun = betadiver, distance = "sim", k=2,  trymax=100)
sol <- metaMDS(distancia, k=2,  trymax=100)

 plot(sol,display = c("species"), choices = c(1, 2))
# ##### 6.1 HACER MATRIZ DE COLOR -----------------------------------------
# upper left: red - upper right: blue
# lower left: yellow - lower right: green

################---------------------------------------------------



my.grid<-expand.grid(x=seq(.01,1,.01),y=seq(.01,1,.01))
my.data<-seq(0,1,.01)

par(mfrow=c(1,1))
#- lower left - upper left
my.class<-classIntervals(my.data,n=100,style="equal")
my.col<-c("yellow","red") # choose colors
my.pal.1<-findColours(my.class,my.col)
plot(rep(0,101),my.data,pch=15,col=my.pal.1, cex=2, xlim=c(0,1),ylim=c(0,1))
#- lower right - upper right
my.class<-classIntervals(my.data,n=100,style="equal")
my.col<-c("green","blue") # choose colors
my.pal.2<-findColours(my.class,my.col)
points(rep(1,101),my.data,pch=15,col=my.pal.2, cex=2)
#- upper left - upper right
my.class<-classIntervals(my.data,n=100,style="equal")
my.col<-c("red","blue") # choose colors
my.pal.3<-findColours(my.class,my.col)
points(my.data,rep(1,101),pch=15,col=my.pal.3, cex=2)
#- lower left - lower right
my.class<-classIntervals(my.data,n=100,style="equal")
my.col<-c("yellow","green") # choose colors
my.pal.4<-findColours(my.class,my.col)
points(my.data,rep(0,101),pch=15,col=my.pal.4, cex=2)
#- horizontal
my.class<-classIntervals(my.data,n=100,style="equal")
my.col<-c(paste(my.pal.1[50]),paste(my.pal.2[50])) # choose colors
my.pal.h<-findColours(my.class,my.col)
points(my.data,rep(0.5,101),pch=15,col=my.pal.h, cex=2)
#- vertical
my.class<-classIntervals(my.data,n=100,style="equal")
my.col<-c(paste(my.pal.4[50]),paste(my.pal.3[50])) # choose colors
my.pal.v<-findColours(my.class,my.col)
points(rep(0.5,101),my.data,pch=15,col=my.pal.v, cex=2)

#------------------------------------
# loop: use left and right vertical ramp
# and interpolate horizontal lines
#------------------------------------

col.matrix<-matrix(nrow = 101, ncol = 101, NA)

for(i in 1:101){
  my.col<-c(paste(my.pal.1[i]),paste(my.pal.2[i])) # choose colors
  col.matrix[102-i,]<-findColours(my.class,my.col)
}

#------------------------------------
# plot full grid
#------------------------------------
par(mfrow=c(1,1))
par(mai=c(0,0,0,0))
plot(rep(0,101),my.data,pch=15,col=my.pal.1, cex=0.5, xlim=c(0,1),ylim=c(0,1),axes=F,ylab="",xlab="")

for(i in 1:101){
  col.temp<-col.matrix[i-1,]
  points(my.data,rep((i-1)/100,101),pch=15,col=col.temp, cex=2)
}
box()

str(col.matrix)


#Then you have to assign the right color value to your NMDS point by doing:

#-------------------------------------------|1
# transform NMDS results into 0-1 space
#-------------------------------------------

NMDS.points.t<-as.data.frame(sol$points)

NMDS.points.t$X<-NMDS.points.t$MDS1+abs(min(NMDS.points.t$MDS1))
NMDS.points.t$Y<-NMDS.points.t$MDS2+abs(min(NMDS.points.t$MDS2))

plot(NMDS.points.t$X,NMDS.points.t$Y)

m=max(NMDS.points.t$X,NMDS.points.t$Y)
NMDS.points.t$X<-NMDS.points.t$X/m
NMDS.points.t$Y<-NMDS.points.t$Y/m

#windows();
plot(NMDS.points.t$X,NMDS.points.t$Y)

#-------------------------------------------
# plot NMDS results using the 2D colors
# write color in HBWID link table
#-------------------------------------------

#windows();
plot(2,2,xlim=c(0,1),ylim=c(0,1),xlab="NMDS 1",ylab="NMDS 2",las=1,axes=F)
text(0.8,0,paste("stress: ",round((sol$stress),2)))
box()
col.HBWID<-matrix(nrow = length(NMDS.points.t$X), ncol = 2, NA)

for(i in 1:length(NMDS.points.t$X)){
  (posX=round(NMDS.points.t$X[i],2)*100)
  (posY=round(NMDS.points.t$Y[i],2)*100)
  col.HBWID[i,1]<-row.names(NMDS.points.t)[i] #$Quad360ID[i] #NMDS.points.t$Quad360ID[i]
  ifelse(noquote(posX)==0 | noquote(posY)==0 , 0, col.HBWID[i,2]<-col.matrix[noquote(posY),noquote(posX)])
  points(NMDS.points.t$X[i],NMDS.points.t$Y[i],col=col.matrix[noquote(posY),noquote(posX)],pch="+")
}

names(col.HBWID)=c("ID","col")
final=cbind(NMDS.points.t,col.HBWID[,2])


#.P O------------------------------------------

grilla2=mask(grilla,mascara)
breaks=as.numeric(row.names(final))
plot(colombia)
plot(grilla2,breaks=breaks,col=as.character(final[,5]),add=T,horizontal=T)

rgb=col2rgb(as.character(final[,5]))
rgb2=as.data.frame(cbind(as.numeric(row.names(final)),t(rgb)))
names(rgb2)=c("pixel_value", names(rgb2)[-1])

writeRaster(grilla2,"MDS",overwrite=TRUE, format="HFA",datatype="INT4S")
write.table(rgb2,"rgb_MDS.clr",sep=" ", row.names=F)

# #### 7. CLUSTER ANALYSIS ------------------------------------------------------


metodos=c("average","single" , "complete","ward","weighted")
######### TODO ESTA SOLO PARA SHP


## cluster con los diferentes metodos
#CLUSTER1=agnes(DISTAN, diss = T, method ="average", keep.diss = F, keep.data = F)
CLUSTER1=hclust(DISTAN,  method ="average")
CLUSTER2=hclust(DISTAN, method ="single" )
CLUSTER3=hclust(DISTAN, method ="complete")
CLUSTER4=hclust(DISTAN, method ="ward")
CLUSTER5=(agnes(DISTAN, diss = T, method ="weighted", keep.diss = F, keep.data = F))
CLUSTER5$height=sort(CLUSTER5$height)
CLUSTER6=hclust(DISTAN, method ="median")
CLUSTER6$height=sort(CLUSTER6$height)
CLUSTER7=hclust(DISTAN, method ="centroid")
CLUSTER8=diana(DISTAN, diss = T, keep.diss = F, keep.data = F)
#arbol Neighbour join tree
NJ=nj(DISTAN)
rotated=midpoint(NJ)
ultrametric=compute.brlen(rotated,power=1)
CLUSTER9=as.hclust(ultrametric)


##COPHENETIC CORRELATION COEFFICIENT

#The cophenetic distance between two observations that have been clustered is defined to be 
#the intergroup dissimilarity at which the two observations are first combined into a single cluster. 
#Note that this distance has many ties and restrictions.

#It can be argued that a dendrogram is an appropriate summary of some 
#data if the correlation between the original distances and the cophenetic distances is high. 
#Otherwise, it should simply be viewed as the description of the output of the clustering algorithm. 

#Cophenetic Distance
#The estimated distance between two points is the level at which they are fused in
#the dendrogram, or the height of the root. A good clustering method correctly
#reproduces the actual dissimilarities. The distance estimated from a dendrogram 
#is called cophenetic distance. The name echoes the origins of hierarchic
#clustering in old fashioned numeric taxonomy.1|+85/9*

#function cophenetic estimates the distances among all points from a dendrogram


#The cophenetic correlation coefficient
##structure from the dendrogram represents the actual distances. This measure is dened
#AS the correlation between the n(n-1)=2 pairwise dissimilarities between observations
#and their cophenetic dissimilarities from the dendrogram, i.e., the between cluster dis-
#similarities at which two observations are rst joined together in the same cluster.
#This parameter measures the correlation between distance values calculated during tree 
#building and the observed distance. The CCC is a measure of how faithfully a dendrogram 
#maintains the original pairwise distances. 


coph1 <- cor(cophenetic(CLUSTER1), DISTAN)
coph2 <- cor(cophenetic(CLUSTER2), DISTAN,use="pairwise.complete.obs") ##NA
coph3 <- cor(cophenetic(CLUSTER3), DISTAN)
coph4 <- cor(cophenetic(CLUSTER4), DISTAN)
coph5 <- cor(cophenetic(CLUSTER5), DISTAN)
coph6 <- cor(cophenetic(CLUSTER6), DISTAN)
coph7 <- cor(cophenetic(CLUSTER7), DISTAN) 
coph8 <- cor(cophenetic(CLUSTER8), DISTAN) 
coph9 <- cor(cophenetic(CLUSTER9), DISTAN) ### dimensiones incompatibles


# ##### 8. NUMERO DE GRUPOS --------------------------------------------------


### 1 particion


CLUSTERS=c("CLUSTER1","CLUSTER2","CLUSTER3","CLUSTER4","CLUSTER5","CLUSTER6","CLUSTER7","CLUSTER8","CLUSTER9") # "CLUSTER2","CLUSTER9"

EVALUACION1=NULL
for (g in 1:(length(CLUSTERS))){
  clust=get(CLUSTERS[g])
  evaluacion=EVA1(clust)
  evaluacion=unique(evaluacion)
  nombre=rep(CLUSTERS[g],nrow(evaluacion))
  eva=cbind(nombre,evaluacion)
  b=length(unique(eva$grupos))
  optimo=optimize(grup,c(2,b),eva,maximum = FALSE)
  c=floor(optimo[[1]])
  opt=as.character(c)
  plot(eva$metrica~eva$grupos, main=clust$method, sub=opt, xlab="N grupos",ylab="altura")
  legend("center",opt)
  abline(lm(eva$metrica[1:c]~eva$grupos[1:c]),col="red")
  abline(lm(eva$metrica[(c+1):b]~eva$grupos[(c+1):b]),col="red")
  eva=cbind(nombre,opt,evaluacion)
  EVALUACION1=rbind(EVALUACION1,eva)
  
}

#######CARO ESTA METRICA NO LA ESTAMOS USANDO POR ESO LA COMENTE

#   EVALUACION2=NULL
# 
# 
# for (g in 1:(length(CLUSTERS))){
#   evaluacion2=NULL
#   clust=get(CLUSTERS[g])
#   evaluacion2=EVA2(clust)
#   nombre2=rep(CLUSTERS[g],nrow(evaluacion2))
#   eva2=evaluacion2
#   b=length(unique(eva2$grupos))
#   optimo2=optimize(grup,c(2,b),eva2,maximum = FALSE)
#   c2=floor(optimo2[[1]])
#   opt2=as.character(c2)
#   plot(eva2$metrica~eva2$grupos, main=clust$method, sub=opt2, xlab="N grupos",ylab="% endemismo")
#   legend("center",opt2)
#   abline(lm(eva2$metrica[1:c]~eva2$grupos[1:c]),col="red")
#   abline(lm(eva2$metrica[(c+1):b]~eva2$grupos[(c+1):b]),col="red")
#   eva2=cbind(nombre2,opt2,evaluacion2)
#   EVALUACION2=rbind(EVALUACION2,eva2)
# }
# 
# 

# ############  9. GRAFICA DE RESULTADOS ----------------------------------


# ####para avarage y ward cortes sucesivos --------------------------------
#revisa los coeficientes y haz los dos mejores, normalmente son UPGMA (1) y Ward(4) 
LEVEL1=unique(EVALUACION1[,1:2])
DF_FINAL=as.data.frame(DF[,1])
row.names(DF_FINAL)=DF[,1]

ngrupos=as.numeric(as.character(LEVEL1[which(LEVEL1$nombre==CLUSTERS[4]),2])) ## se tomo el obtimo de ward (CLUSTER[4]), por ser el minimo de los optimos 
# tamebin puedes tomar el optimo de el mejor cluster y hacer cortes sucesivos hasta ese corte
tmp=mascara
for ( p in c(1,4)){ # p=1 para UPGMA 4 para ward, 
  cluster=get(CLUSTERS[p])
  cluster$height=sort(cluster$height)
  nombre=cluster$method
  
  
  #pdf(file=paste(nombre,"pdf",sep="."))
  for (t in 2:ngrupos){                #c(80,100,500)){ cambiar esto asi  si se quiere solo cierto corte  eg, 80,100,500    
    finalcuts=as.data.frame(cutree(cluster,k=t))
    names(finalcuts)=c("group")
    
    DF_FINAL[,ncol(DF_FINAL)+1]=NA
    DF_FINAL[row.names(DF_FINAL)==row.names(finalcuts),ncol(DF_FINAL)]=finalcuts$group[row.names(finalcuts)==row.names(DF_FINAL)]
    
    
    tmp@data[,ncol(tmp@data)+1]=NA
    filas2=(DF_FINAL[,1])
    for (i in filas2){
      tmp@data[tmp@data$ID==i, ncol(tmp@data)]=DF_FINAL[DF_FINAL[,1]==i,ncol(DF_FINAL)]
    }
    
    titulo=paste0(nombre,t)
    names(tmp)=c(names(tmp)[-ncol(tmp)],titulo)
  }
  #dev.off()
}

writePolyShape(tmp, "grupos.shp", factor2char = TRUE, max_nchar=254)


######################## HASTA ACA !!!!!!!!!
############## SI PLOTEAMOS POLIGOMOS



# si OROBIOMA es shp
tmp$clase=0

for (i in DF$ID){
  tmp$clase[tmp$Unit_ID==i]=DF$V1952[DF$ID==i]
}

spplot(tmp,zcol="clase")