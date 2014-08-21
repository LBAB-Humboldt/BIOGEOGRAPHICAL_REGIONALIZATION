### FUNCIONES

# cargar paquetes
loadLibraries<-function(pkg){
  loadPkg<-require(pkg,,character.only=TRUE)
  if(loadPkg){
    print(paste0(pkg," is loaded correctly"))
  } else {
    print(paste0("trying to install ",pkg))
    install.packages(pkg)
    if(require(pkg,character.only=TRUE)){
      print(paste0(pkg," installed and loaded"))
    } else {
      stop(paste0("could not install ",pkg))
    }
  }
}

#grafica recambio.... grafica el recambio de una celda o poligono contra todas las demas

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
grafica_recambio_shp=function(colum){#,colum){
  # ver el area de estudio
  mascara2=mascara
  mascara2$recambio=as.numeric(0)
  colum1=as.data.frame(colum)
  #names(colum1)=c("dist")
  valores=as.numeric(row.names(colum1))
  mascara2@data[which(mascara2$Unit_ID%in%valores), ncol(mascara2@data)]=colum1[which(as.numeric(row.names(colum1))%in%valores),1]
  mascara2$recambio=as.numeric(mascara2$recambio)
  
  #resultados=c(mascara2) 
  return(mascara2)
} 

# ESPECIES... identifica que especies estan presentes en una celda o poligono dado
ESPECIES=function(xy){
  xy=t(as.matrix(xy))
  celda=cellFromXY(grilla2,as.numeric(xy))
  fila=DF2[which(row.names(DF2)==as.character(celda)),]
  sp=fila[,which(colSums(fila)!=0)]
  return(colnames(sp))
}
ESPECIES_shp=function(celda){
  celda2=celda
  fila=DF[which(row.names(DF)==as.character(celda2)),]
  sp=fila[,which(colSums(fila)!=0)]
  return(colnames(sp))
}

#colores de la paleta para  la grafica del ordenamiento multidimensional 
col.br <- colorRampPalette(c("blue", "cyan", "yellow","gray"))

# metricas para evaluacion 
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


### FUNCIONES PARALELO
countByUnit<-function(spFile,zones,ignoreCells,outFolder){
  spRaster<-raster(spFile)
  tmpRaster<-zones*spRaster
  unitCount<-as.data.frame(table(tmpRaster[ignoreCells]))
  if(nrow(unitCount)<=1){return}
  colnames(unitCount)<-c("Unit_ID",names(spRaster))
  write.csv(unitCount,paste0(outFolder,"/",names(spRaster),".csv"),row.names=F)
}

remExt<-function(input){
  return(unlist(strsplit(input, "\\_10p_cut.csv"))[1])
}

concatTables<-function(dbfs,field,rowIDs){
  concatTable<-as.data.frame(matrix(0,nrow=length(rowIDs),ncol=(length(dbfs)+1)))
  colnames(concatTable)<-c("ID",sapply(dbfs,remExt,USE.NAMES=FALSE))
  concatTable[,1]<-rowIDs
  for (i in 1:length(dbfs)){
    oneDbf<-read.csv(paste(path,dbfs[i],sep="\\"),as.is=TRUE)
    tmp_df1<-data.frame(ID=rowIDs)
    tmp_df2<-merge(tmp_df1,oneDbf,by.x="ID",by.y="Unit_ID",all.x=TRUE,incomparables=NA)
    concatTable[,i+1]<-tmp_df2[,field]
  }
  return(concatTable)
}

staky_original=function(sp.files){
  mod1=grilla
  for (i in 1:length(sp.files)){
    sfCat(paste0("inicia ",sp.files[i]," hora=",date()),sep="\n")
    mod2=aggregate(raster(sp.files[i]),factor)
    mod3=resample(mod2,grilla)
    mod4=reclassify(mod3,clasificacion)
    mod1=addLayer(mod1,mod4)
    sfCat(paste0("termina ",sp.files[i]," hora=",date()),sep="\n")
  }
  DF.pre=as.data.frame(mod1)
  return(DF.pre)
}

staky1=function(sp.files){
  trace<-tryCatch(
{
  mod1=grilla
  for (i in 1:length(sp.files)){
    sfCat(paste0("inicia ",sp.files[i]," hora=",date()),sep="\n",file=paste0(ruta_salida,"/Res_core_Sp_",sp.files[1]))
    mod2=aggregate(raster(sp.files[i]),factor)      
    mod3=resample(mod2,grilla)
    mod4=reclassify(mod3,clasificacion)
    mod1=addLayer(mod1,mod4)
    sfCat(paste0("termina ",sp.files[i]," hora=",date()),sep="\n")
  }
  DF.pre=as.data.frame(mod1)
  return(DF.pre)
},
error = function(e){cat(paste0(sp.files[i]," ",e, date()),sep="\n")})

}

staky_Agre=function(sp.files){
  sink(paste0(ruta_salida,"/Res_cores",sp.files[1],".txt"))
  on.exit(sink())
  trace<-tryCatch(
{
  mod1=grilla
  for (i in 1:length(sp.files)){
    sfCat(paste0("inicia ",sp.files[i]," hora=",date()),sep="\n")
    mod2=aggregate(raster(sp.files[i]),factor)      
    mod3=resample(mod2,grilla)
    mod4=reclassify(mod3,clasificacion)
    mod1=addLayer(mod1,mod4)
    sfCat(paste0("termina ",sp.files[i]," hora=",date()),sep="\n")
  }
  DF.pre=as.data.frame(mod1)
  return(DF.pre)
},
error = function(e){cat(paste0(sp.files[i]," ",e, date()),sep="\n")})
}

staky.h=function(Sp.humedales){
  mod1=grilla
  for (i in 1:length(Sp.humedales)){
    print(i)
    mod2=aggregate(raster(Sp.humedales[i]),factor)
    mod3=resample(mod2,grilla)
    mod4=reclassify(mod3,clasificacion)
    mod1=addLayer(mod1,mod4)
  }
  DF.pre=as.data.frame(mod1)
  return(DF.pre)
}

alfas=function(sp.files){
  mod1=raster(sp.files[1])
  for (i in 2:length(sp.files)){
    mod2=raster(sp.files[i])
    mod3=resample(mod2,mod1)
    extent(mod3)=extent(mod1)
    mod1=sum(mod1,mod3)
  }
  return(mod1)}

gamas=function(sp.files){
  mod1=aggregate(raster(sp.files[1]),factor)
  mod1=resample(mod1,grilla)
  mod1=reclassify(mod1,clasificacion)
  for (i in 2:length(sp.files)){
    mod2=aggregate(raster(sp.files[i]),factor)
    mod3=resample(mod2,grilla)
    mod4=reclassify(mod3,clasificacion)
    mod1=sum(mod1,mod4)
  }
  return(mod1)}


CORTES_shp=function(CLUSTERS){ 
  cluster=(CLUSTERS)
  cluster$height=sort(cluster$height)
  nombre=cluster$method
  
  for (t in 2:ngrupos){
    #print(t)#c(80,100,500)){ cambiar esto asi  si se quiere solo cierto corte  eg, 80,100,500    
    finalcuts=as.data.frame(cutree(cluster,k=t))
    names(finalcuts)=c("group")
    
    DF_FINAL[,ncol(DF_FINAL)+1]=NA
    DF_FINAL[which(row.names(DF_FINAL)%in%row.names(finalcuts)),ncol(DF_FINAL)]=finalcuts$group[which(row.names(finalcuts)%in%row.names(DF_FINAL))]
    
    tmp@data[,ncol(tmp@data)+1]=NA
    filas2=(DF_FINAL[,1])
    
    tmp@data[which(tmp@data$ID%in%filas2), ncol(tmp@data)]=DF_FINAL[which(DF_FINAL[,1]%in% filas2),ncol(DF_FINAL)]
    titulo=paste0(nombre,t)
    names(tmp)=c(names(tmp)[-ncol(tmp)],titulo)
  }
  writePolyShape(tmp, paste0(ruta_salida,"/",nombre,".shp"), factor2char = TRUE, max_nchar=254)
  return(tmp)
}

CORTES_raster=function(CLUSTERS){ 
  cluster=(CLUSTERS)
  cluster$height=sort(cluster$height)
  nombre=cluster$method
  
  for (t in 2:ngrupos){
    #print(t)#c(80,100,500)){ cambiar esto asi  si se quiere solo cierto corte  eg, 80,100,500    
    finalcuts=as.data.frame(cutree(cluster,k=t))
    names(finalcuts)=c("group")
    
    DF_FINAL[,ncol(DF_FINAL)+1]=NA
    DF_FINAL[which(row.names(DF_FINAL)%in%row.names(finalcuts)),ncol(DF_FINAL)]=finalcuts$group[which(row.names(finalcuts)%in%row.names(DF_FINAL))]
    
    tmp@data[,ncol(tmp@data)+1]=NA
    filas2=(DF_FINAL[,1])
    
    tmp@data[which(tmp@data$ID%in%filas2), ncol(tmp@data)]=DF_FINAL[which(DF_FINAL[,1]%in% filas2),ncol(DF_FINAL)]
    titulo=paste0(nombre,t)
    names(tmp)=c(names(tmp)[-ncol(tmp)],titulo)
  }
  writePolyShape(tmp, paste0(ruta_salida,"/",nombre,".shp"), factor2char = TRUE, max_nchar=254)
  return(tmp)
}
