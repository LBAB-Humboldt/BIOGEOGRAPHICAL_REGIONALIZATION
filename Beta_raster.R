Beta_raster=function(){
  
  rm(list=ls(all=T))
  gc() ##para borrar todo lo que quede oculto en memoria
  memory.limit(size = 1000000) ## para asiganar 20 gigas del disco duro a ram
  
  source("w:/Beta/beta-functions.R") #Load auxilliary functions
  
  ############# CARGAR  PAQUETES --------------------------------
  
  pkgs<-c("raster","vegan","maptools","clustsig","PBSmapping", "epitools","RColorBrewer",
          "classInt","sp","cluster","ape","phangorn","rgdal","svGUI","svDialogs","snowfall")
  
  lapply(pkgs,loadLibraries) #Load required libraries
  
  # #### 0. DEFINIR RUTAS, CORES Y FUNCIONES ---------------------------------------------------
  
  ## DEFINIR RUTAS 
  
  ruta_modelos<-(dlgDir(default = getwd(), title="ESPECIFIQUE LA RUTA DE LOS MODELOS")$res)
  ruta_iucn<-(dlgDir(default = getwd(), title="ESPECIFIQUE LA RUTA DE MAPAS DE EXPERTO")$res)
  ruta_salida<-(dlgDir(default = getwd(), title="ESPECIFIQUE LA RUTA DE LOS RESULTADOS")$res)
  
  
  ncores=as.numeric(dlgInput("defina el numero de cores a usar")$res)
  
  
  # ### 1. AREA DE ESTUDIO --------------------------------------------------
  
  #MASCARAS # 
  
  #Seleccionar .shp del area de estudio
  spFile<-dlgOpen(title = "Seleccione .shp del ?rea de estudio", filters = "*.shp")$res
  mascara=readShapePoly(spFile) #"F:/InfoGeo/COLOMBIA.shp"
  extent=extent(mascara)
  
  # Seleccionar ?rea de estudio en raster para crear grilla
  area.Ras<-dlgOpen(title = "Seleccione raster del ?rea de estudio", filters = "*.shp")$res
  area.raster<-raster(area.Ras)
  celdas<-Which(!is.na(area.raster),cells=T)
  
  #Para agrerar celdas para disminuir resoluci?n  
  factor=5 #derfine nivel de agregacion 
  aoi=aggregate(area.raster,factor)
  aoi=mask(aoi,mascara)
  
  #crear grilla
  grilla=aoi
  names(grilla)="grilla"
  grilla[1:ncell(grilla)]<-1:ncell(grilla)
  
  
  # ### 2. LISTA DE ESPECIES -----------------------------------------------
  
  #Modelos 
  
  sp.files<-list.files(ruta_modelos,pattern="*_10p_cut.tif$") # Aqu? se seleccionan todos los archivos en la carpeta. Todas las especies que est?n aqu? son incluidas en el an?lisis. 
  #Con los modelos disponibles a Julio de 2014 tienen problemas en la especie Apayana amydalina por lo que se saca del an?lsis  
  sp.files<-sp.files[-which(sp.files=="Ayapana_amygdalina_10p_cut.tif")]
  sp.names<-t(data.frame(strsplit(sp.files, "_")))
  sp.names<-paste0(sp.names[,1],"_",sp.names[,2])
  
  #Mapas IUCN 
  
  sp.files2<-list.files(ruta_iucn,pattern="*.tif")  
  sp.names2<-as.matrix(t(data.frame(strsplit(sp.files2,".tif"))))
  sp.names2<-as.vector(sp.names2[,1])
  
  
  #Todos 
  
  nombres.sp<-as.character(c(sp.names,sp.names2))
  
  
  ### 3. GENERAR TABLA  DE PRESENCIAS -------------------------------------
  
  ## 3.1 Tabla de presencia para especies con modelos. 
  
  setwd(ruta_modelos)
  
  clasificacion=matrix(c(0,0.25,0,0.25,1,1),nrow=2,ncol=3,byrow=T)
  
  breaks=seq(1,length(sp.files),by=round(length(sp.files)/(ncores)))
  
  Grupos= list(sp.files[1:breaks[2]])
  
  for (i in 2:ncores){
    if (i <ncores) {
      grupo<-sp.files[(breaks[i]+1):(breaks[i+1])]
    } 
    else {
      grupo<-sp.files[(breaks[i]+1):length(sp.files)]
    }
    Grupos[i] <- list(grupo)
  }
  
  Grupos_iucn= list(sp.files2[1:breaks[2]])
  
  for (i in 2:ncores){
    if (i <ncores) {
      grupo<-sp.files2[(breaks[i]+1):(breaks[i+1])]
    } 
    else {
      grupo<-sp.files2[(breaks[i]+1):length(sp.files2)]
    }
    Grupos_iucn[i] <- list(grupo)
  }
  
  ##3.1.1 PARELELIZADO: Crear Stacks para cada grupo de especies y pasar a data frame
  
  setwd(ruta_modelos)
  sfInit(parallel=T,cpus=ncores) # slaveOutfile=paste0(ruta_salida,"/LogFile_Corrida_Julio2014.txt")
  sfExport(list=c("staky_Agre","grilla", "clasificacion", "Grupos","factor","ruta_salida"))
  sfLibrary(rgdal)
  sfLibrary(raster)
  sfLibrary(snowfall)
  DF_pre<-sfClusterApplyLB(Grupos,staky_Agre)
  sfStop()
  
  setwd(ruta_iucn)
  sfInit(parallel=T,cpus=ncores) # slaveOutfile=paste0(ruta_salida,"/LogFile_Corrida_Julio2014.txt")
  sfExport(list=c("staky_Agre","grilla", "clasificacion", "Grupos_iucn","factor","ruta_salida"))
  sfLibrary(rgdal)
  sfLibrary(raster)
  sfLibrary(snowfall)
  DF_iucn<-sfClusterApplyLB(Grupos_iucn,staky_Agre)
  sfStop()
  
  
  #convertir a dataframe
  DF=as.data.frame(DF_pre)
  DF=DF[,-(agrep("grilla.", names(DF), max=1,ignore.case = F)[-1])] #limpiar  ID de mas generado al pegar los DF's preliminares
  names(DF)=c("ID",sp.names)
  
  save.image(paste0(ruta_salida,"/Corrida_", format(Sys.Date(), "%b_%d_%Y"),".RData"))
  
  
  write.csv(DF,paste0(ruta_salida,"/Tabla_presencias_sp_modelos.csv"))

  DF_prueba=apply(DF_prueba,c(1,2),function(x) {if (is.na(x)){x=0} else{x=x}})
  vacias2=c(which(rowMeans(DF_prueba[,2:ncol(DF_prueba)],na.rm=T)==0)) # solo estamos limpiando vacias , se puede menos de 5 (kreft)
  if (length(vacias2)==0){ DF_prueba=DF_prueba} else {DF_prueba=DF_prueba[-vacias2,]}
  DF_prueba<-as.data.frame(DF_prueba)
  
  DF=apply(DF,c(1,2),function(x) {if (is.na(x)){x=0} else{x=x}}) #remover filas vacias
  vacias2=c(which(rowMeans(DF[,2:ncol(DF)],na.rm=T)==0)) # solo estamos limpiando vacias , se puede menos de 5 (kreft)
  if (length(vacias2)==0){ DF=DF} else {DF=DF[-vacias2,]}
  
  ###A?adir paso para revisar que no haya n?meros diferentes a 1 y 0 ###
  
  # ### 4. GENERAR MATRIZ DE DISTANCIA --------------------------------------
  
  # Este paso aun no esta paralelizado este es el paso mas largo
  
  DISTAN2=betadiver(DF, "sim")
  distancia=as.matrix(DISTAN2)
  
  
  #5. PATRONES DE DIVERSIDAD  -----------------------------------------------
  
  setwd(ruta_modelos)
    
  ## ALFA  para col y amazonia
  sfInit(parallel=T,cpus=ncores)
  sfLibrary(raster)
  sfLibrary(rgdal)
  sfExport(list=c("groups2","alfas"))
  system.time(ALFAS<-sfClusterApplyLB(groups2,alfas))
  sfStop()
  #Revisar funci?n porque no est? funcionando. 
    
  setwd(ruta_salida)
  alfa1km=sum(stack(ALFAS))
  save(alfa1km,file="alfa1km.RData")
  ALFA=aggregate(alfa1km,factor,fun=mean, na.rm=T)
  plot(ALFA)
  
  
  ## GAMA para col y amazonia
  sfInit(parallel=T,cpus=ncores)
  sfLibrary(raster)
  sfLibrary(rgdal)
  sfExport(list=c("groups2","gamas","grilla", "clasificacion","factor"))
  system.time(GAMAS<-sfClusterApplyLB(groups2,gamas))
  sfStop()
  

  setwd(ruta_salida)
  GAMA=sum(stack(GAMAS))
  save(GAMA,file="GAMA.RData")
    
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
  
  writeRaster(alfa1km,"alfa1km",overwrite=TRUE)  
  writeRaster(ALFA,"ALFAmean",overwrite=TRUE)
  writeRaster(GAMA,"GAMA",overwrite=TRUE)
  writeRaster(BETAm,"BETAm.tif",overwrite=TRUE)
  writeRaster(BETAa,"BETAa.tif",overwrite=TRUE)
  writeRaster(BETAmg,"BETAmg.tif",overwrite=TRUE)
  
  
  #####  6. GRAFICAS DE RECAMBIO Y LISTAS SP -----------------------------------------
  n=as.numeric(dlgInput("defina el numero de areas para las cuales desea ver su recambio")$res)
  
  plot.new()
  plot(extent,main="HAGA CLICK SOBRE SUS POLIGONOS DE INTERES", sub="AVISO: esto puede demorarse,espere que aparezca el cursor y el mapa")
  plot(mascara,add=T)
  xy=locator(n=n)
  dev.off()
  
  xy1=SpatialPoints(xy)
  celda=over(xy1,mascara)$Unit_ID
  columna=distancia[,which(colnames(distancia)%in%as.character(celda))]
  colum=split(columna,col(columna))
  
  ##paralelo para  shp con distancias 
  
  sfInit(parallel=T,cpus=n)
  sfLibrary(raster)
  sfLibrary(rgdal)
  sfLibrary(sp)
  sfLibrary(maptools)
  sfLibrary(snow)
  sfExport(list=c("colum","mascara","grafica_recambio_shp"))
  system.time(RECAMBIOS<-sfClusterApplyLB(colum,grafica_recambio_shp))
  sfStop()
  
  
  ##paralelo para  listado de especies en los poligonos seleccionados
  
  sfInit(parallel=T,cpus=n)
  sfLibrary(raster)
  sfLibrary(rgdal)
  sfLibrary(snow)
  sfExport(list=c("celda","DF","ESPECIES_shp"))
  system.time(SP<-sfClusterApplyLB(celda,ESPECIES_shp))
  sfStop()
  
  ##GUARDAR RESULTADOS
  
  setwd(ruta_salida)
  for (i in 1:length(RECAMBIOS)){
    nombre=paste0("_",i)
    writePolyShape(RECAMBIOS[[i]], paste0("dist",nombre,".shp"), factor2char = TRUE, max_nchar=254)
    write.csv(SP[[i]],paste0("sp",nombre,".csv"))
  }
  
  
  
  # ###### 7. CLUSTERS ------------------------------------------------------
  
  tmp=mascara
  tmp$ID=1:nrow(tmp)
  lista.metodos=list("average","single" , "complete","ward","median","centroid")
  metodos<-dlgList(lista.metodos,multiple=TRUE,title="Seleccione el(los) metodos de agrupamiento que desea probar")$res
  
    
  ## cluster con los diferentes metodos
  
  for (t in 1:length(metodos)){
    nombre.clu=paste0("CLUSTER",t)
    metodo.usar=metodos[t]
    assign(nombre.clu,hclust(DISTAN2,  method =metodo.usar))
    
  }
  
  #NOTAS: Falta incorporar el metodo weighted por que ese solo esta disponible en agnes,
  #tocaria pasar todo a agnes y hacer sort de las height, tambien falta incluir el metodo
  #diana (arriba a  abajo)
  #CLUSTER1=agnes(DISTAN, diss = T, method ="average", keep.diss = F, keep.data = F)
  #CLUSTER5=(agnes(DISTAN2, diss = T, method ="weighted", keep.diss = F, keep.data = F))
  #CLUSTER5$height=sort(CLUSTER5$height)
  
  
  CLUSTERS=vector("list",length(metodos))
  for (j in 1:length(metodos)){
    CLUSTERS[[j]]=get(paste0("CLUSTER",j))
  }
  
  DF_FINAL=as.data.frame(row.names(DF))
  row.names(DF_FINAL)=row.names(DF)
  
  ngrupos=as.numeric(dlgInput("Especifique el numero maximo de grupos que quiere ver")$res)
  
  
  sfInit(parallel=T,cpus=length(CLUSTERS))
  sfLibrary(raster)
  sfLibrary(stats)
  sfLibrary(sp)
  sfLibrary(maptools)
  sfLibrary(snow)
  sfExport(list=c("CLUSTERS","DF_FINAL","CORTES_shp", "tmp", "ngrupos","ruta_salida"))
  system.time(CORTES_CLUST<-sfClusterApplyLB(CLUSTERS,CORTES_shp))
  sfStop()
  
  
  
  
  #Hallar correlaciones
  
  ## correlaciones para los cluster con los diferentes metodos
  
  for (k in 1:length(metodos)){
    nombre.clu=paste0("coph",k)
    metodo.usar=CLUSTERS[[k]]
    assign(nombre.clu,cor(cophenetic(CLUSTER1),DISTAN2))
    
  }
  
  coph1 <- cor(cophenetic(CLUSTER1), DISTAN)
  #coph2 <- cor(cophenetic(CLUSTER2), DISTAN,use="pairwise.complete.obs") ##NA
  coph3 <- cor(cophenetic(CLUSTER3), DISTAN)
  coph4 <- cor(cophenetic(CLUSTER4), DISTAN)
  coph5 <- cor(cophenetic(CLUSTER5), DISTAN)
  coph6 <- cor(cophenetic(CLUSTER6), DISTAN)
  coph7 <- cor(cophenetic(CLUSTER7), DISTAN) 
  
  
}