  ## Curso de Genética del paisaje   INECOL AC Septiembre 2016
  ## Yessica Rico M 
  ## Ejercicio: Relación entre la conectividad del parche y la diversidad genética poblacional en Dianthus carthusianorum
  ## Datos & análisis Rico et al (2014) Con Biol 2,
--------------------------------------------------
   install.packages("MuMIn")
  
   require(MuMIn)
 # Importar datos diversidad genética por población (patch) 
   Diversity<-read.csv("Dc_diversity.csv", header=TRUE)   
   a <- c(1:nrow(Diversity))[order(Diversity$patch)]
   Diversity <-Diversity[a,]
   dim(Diversity)
   names(Diversity)
  
  ################################################################
  # Parte 1: Estimación de índices de conectividad por parche: Hanski
  # Para el cálculo se necesitan matrices de distancias entre pares de parches(dij), párametros de parches fuente (area Aj,tamaño poblacional Nj, presencia de la especie pj)  
  
  # Importar datos de coordenadas XY de todos los parches 
   patches <- read.csv("PATCH_XY_Dc.csv", header=TRUE )
   a <-c(1:nrow(patches))[order(patches$patch)] # ordena los parches ascendentemente
   a <-c(1:nrow(patches))[order(patches$patch)]
   patches <- patches[a,]
   dim(patches) # número de parches

  # Importar datos de distancias entre parches modelos predictivos (ver doc word para explicacion de los modelos dij) 
   #  pastore intensivo
   diShecte <- read.table("dist_sheep_cte_Dc_Sep011.txt", header=TRUE)  
   b <-c(1:nrow(diShecte))[order(names(diShecte))]
   diShecte<- diShecte[b,b]
   
   # pastore temporal
   diSheint <- read.table("dist_sheep_int_Dc_Sep011.txt", header=TRUE)  # Sheep intermitte
   a <- c(1:nrow(diSheint))[order(names(diSheint))]
   diSheint <- diSheint[a,a]

   # pastore presencia o ausencia 
   diShenu <-read.table("dist_Sheep_null_Dc_Sep011.txt", header=T) # Sheep null
   a <- c(1:nrow(diShenu))[order(names(diShenu))]
   diShenu<- diShenu[a,a]

   # efecto matriz (bosque)
   diveg <- read.table("dist_forest_Nov2011v2.txt", header=T)
   a<-c(1:ncol(diveg))[order(colnames(diveg))]
   b<-c(1:nrow(diveg))[order(rownames(diveg))]
   diveg<-diveg[b,a]
   diveg<-as.matrix(diveg, diag=T)

      dim(diShecte)
      dim(diSheint)
      dim(diShenu)
      dim(diveg)

    # distancias Euclidianas entre parches
   xy<-cbind(patches$x, patches$y)
   distance<-dist(xy, diag=T)
   distance2<-(distance/1000)
   distance2<-as.matrix(distance2)
   d<-as.matrix(distance2)
   
   # Matrices de distancias 
   diShecte <- as.matrix(diShecte, diag=T)
   diSheint <- as.matrix(diSheint, diag=T)
   diShenu <- as.matrix(diShenu, diag=T)
  
   data.frame(colnames(diShecte), colnames(diveg), colnames(diSheint),
   colnames(diShenu), patches$patch)

   # data frame que incluye todos los modelos de distancias dij
   dModels <- list(Eu = distance2, Shecte = diShecte, Sheint = diSheint, Shenu = diShenu, Forest=diveg)

   # Para calcular el indice de Hanski, primero es necesario optimizar el valor de alfa para cada modelo de distancias usando datos de presencia-ausencia en dos periodos 1989 y 2009
   # La función "get.alphafit" optimiza alfa para cada modelo y los valores se traspasan al cuado table.alpha
     
     table.alpha <- matrix(NA, nrow=5)
     dimnames(table.alpha) <- list(models=(c("Eu", "Shecte", "Sheint", "Shenu", "Forest")))

      nseq <-100
      alpha <- seq(0.1,2.5, length = nseq)

       get.alphafit <- function(alpha, d, pj, Op)
       {
         expo<-exp(-alpha* d )
         diag(expo)<-0
         matr<-sweep(expo,2, pj, "*")
         Si <-rowSums(sweep(matr, 2, Op/2, "*"), na.rm=TRUE)
         mod<- glm(cbind(Op,2 -Op) ~ Si, family=binomial)
         deviance(mod)
        }

     # Parameters
      pj <- patches$Dc.09
      Op <- (patches$Dc.89 + patches$Dc.09)

      for(m in 1:length(dModels))
      {
      table.alpha[m] <- (optimize(get.alphafit, interval=alpha, d=dModels[[m]],pj, Op)$minimum)
       }

    # Establecer parametros  Aj & Nj para los parches fuente
    # Parameters Si
     Aj <- patches$Ha
     Nj <- patches$pop09

   # Generación de cuadro con resultados de los indices de conectividad Si para cada parche y cada modelo dij
     Si <- data.frame(matrix(NA,nrow(patches),ncol=15))
     dimnames(Si) <- list(row.names(patches),
     paste(rep(names(dModels), rep(3,5)), rep(c("pj", "Aj", "Nj"),5), sep="_"))

     Source<- data.matrix(data.frame(pj=pj, Aj=pj*Aj, Nj=Nj))
     mod<-rep(1:5, rep(3,5))
     sb <- rep(1:3,5)

     ##### PARTE 2####################
     # Función que calcula el índice Si para cada parche y cada modelo dij
     get.Si <- function(alpha, d, Ap)
    {
      expo<-exp(-alpha*d)
      diag(expo)<-0
      matr<-sweep(expo,2, Ap, "*")
      S <- rowSums(sweep(matr, 2, Op/2, "*"), na.rm=TRUE)
    }

    for (n in 1:ncol(Si))
    {
    Si[,n] <- get.Si(alpha=table.alpha[mod[n]], d=dModels[[mod[n]]], Ap=Source[,sb[n]])
    }

    # Cuadro de resultados Si
    head(Si) # para abreviaturas ver doc word        
    
    #################################################
    #Parte 3 Selección sólo de los datos con presencia de la especie

    Si2 <-data.frame(Si, patch= patches$patch, Sampled = patches$Sampled)
    Si2 <- Si2[which(Si2$Sampled=="1"),]
    a <-order(Si2$patch)
    Si2 <- Si2[a,]
     
     # Cuadro son solo presencia de la especie
     dim(Si2)
     dim(Diversity)
     
    # Debido a que estamos usando índices de diversidad genética a nivel poblacional, es necesario eliminar parches con abundacia menor a 6 individuos, para evitar efectos de la deriva genética y estocasticidad demografica 
    Si2 <- Si2[(Si2$patch)!="A16",] 
    Si2 <- Si2[(Si2$patch)!="A39",] 
    Si2 <- Si2[(Si2$patch)!="A40",] 
    Si2 <- Si2[(Si2$patch)!="A46",] 
    Si2 <- Si2[(Si2$patch)!="A11",] 
    Si2 <- Si2[(Si2$patch)!="G21a",] 
    Si2 <- Si2[(Si2$patch)!="A20",] 
    Si2 <- Si2[(Si2$patch)!="A44",]  
    Si2 <- Si2[(Si2$patch)!="A02",] 

    # Eliminación de poblaciones pequeñas para diversidad genética

     Diversity <- Diversity[(Diversity$patch)!="A16",]
     Diversity <- Diversity[(Diversity$patch)!="A39",]
     Diversity <- Diversity[(Diversity$patch)!="A40",]
     Diversity <- Diversity[(Diversity$patch)!="A46",]
     Diversity <- Diversity[(Diversity$patch)!="A11",]
     Diversity <- Diversity[(Diversity$patch)!="G21a",]
     Diversity <- Diversity[(Diversity$patch)!="A44",]
     Diversity <- Diversity[(Diversity$patch)!="A20",]
     Diversity <- Diversity[(Diversity$patch)!="A02",]
     
     # Revisando el orden y dimension de los datos
     data.frame(patches= Si2$patch, patches= Diversity$patch)
     length(Si2$patch)
     length(Diversity$patch)

      Si <-Si2[,1:15]
      dim(Si)
     length(Diversity$patch)

    ##############################
    # Parte 4 Correlaciones de indices de conectividad para 15 modelos y diversidad genética
    ## Correlaciones con riqueza de alelos A, para analizar heterogocidad cambiar A por He 
    
       for (i in 1:15)
      {
        cat("\n", names(Si)[i],":", "\n")
        print(summary(lm(Diversity$A~ Si[,i])))
      }
    ##############################
    # Parte 5. Selección de modelos & valor acumulado por cada parametro 
      
      require(MuMIn)
     
     # Ajuste de modelos posibles a lso datos genéticos. Califica por orden de importancia todas las combinaciones posibles de predictores(5 distancias dij & tres parametros para parches fuente Aj, Pj, Nj) 
     ## El criterio utiliazado para la jerarquizacion de modelos es AICc y R2 ajustada
     ## Primero se generan modelos lineares (lm) y despues con la funcion dredge se hace la seleccion de modelos
     ## la funcion get models lista los mejores modelos & mientras que "importance" obtiene los modelos con mejores pesos acumulados wij. El valor de cada modelo va de 0 a 1, entre mas alto mayor soporte tiene el modelo
     # Nota: el argumento m.lim = c(NA, 1) especifica que solo se incluya un predictor, es decir un modelo ej. A~ Shenu_Aj , de lo contrario se incluirian combinaciones distintas entre los 15 modelos como predictores, es decir: A ~ Shenu_Aj + Shenu_Pj + Sheint_Aj....etc
      fm2 <- lm((Diversity$A) ~ ., data = data.frame(Si))
      options(na.action = "na.fail")
      (dd2 <- dredge(fm2,  m.lim = c(NA, 1), rank='AIC', extra="R^2"))
      top.models <- get.models(dd2, cumsum(c(0,weight)) <= .95)
      importance(top.models)

      # En este ejemplo, el interes es conocer cual parametro de conectividad dij (5 modelos) y cual parametro de parche fuente (tres parametros Aj, pj NJ) son los mejores para explicar la diversidad genetica, para ello se tiene que obtener para cada parametro de manera individual, su valor acumulado wij de entre todos los modelos.
      # Las siguientes lineas permiten obtener el valor para cada parametro

      test <- strsplit(names(dd2), "_")
      test <- test[sapply(test, length)==2]
      Parameters <- matrix(unlist(test), length(test), 2, byrow=TRUE)
      dd.01 <- !is.na(data.matrix(dd2[,c(1:nrow(Parameters)+1)]))

      Models<- list()
      Models$Eu <- apply(dd.01[,Parameters[,1]=="Eu"],1,max)
      Models$Matrix<- apply(dd.01[,Parameters[,1]=="Forest"],1,max)
      Models$Shecte <- apply(dd.01[,Parameters[,1]=="Shecte"],1,max)
      Models$Sheint<- apply(dd.01[,Parameters[,1]=="Sheint"],1,max)
      Models$Shenu<- apply(dd.01[,Parameters[,1]=="Shenu"],1,max)

      Models$pj  <- apply(dd.01[,Parameters[,2]=="pj"],1,max)
      Models$Aj <- apply(dd.01[,Parameters[,2]=="Aj"],1,max)
      Models$Nj <- apply(dd.01[,Parameters[,2]=="Nj"],1,max)

      # Importancia relativa por parametro

      Importance <- rep(NA,length(Models))
      names(Importance) <- names(Models)

      for(i in 1:length(Models))
      {
        Importance[i] <- sum(dd2$weight * unlist(Models[i]))
      }

      # Cuandro con los valores relativos de importancia wij para cada parametro en el modelo
      # Entre más alto el valor mejor el soporte del parametro
      round(Importance,2)

      # Grafica con los valores de importancia
      barplot((Importance), ylim=c(0,1.0))
       lines(rep(6,2), c(0,1), lwd=2, lty=8)
      
      # Parte 6: Análisis del mejor modelo seleccionado
      # ----------
      # only diversity and Si connectivity
   
      sol <- lm((Diversity$A)~(Si$Sheint_pj))
      summary(sol)
      
