
if (!require("shiny"))
  install.packages("shiny")
library(shiny)

if (!require("shinythemes"))
  install.packages("shinythemes")
library(shinythemes)


if (!require("vcfR"))
  install.packages("vcfR")
library(vcfR)


if (!require("dplyr"))
  install.packages("dplyr")
library(dplyr)


if (!require("tidyverse"))
  install.packages("tidyverse")
library(tidyverse)




shinyServer(function(input,output){
  
   observeEvent(input$aplicar, {
  
     
     
     #################### IMPORTACION Y LECTURA DEL ARCHIVO VCF ########################
     
    # Importacion de los datos del VCF de entrada
    if (is.null (input$archivo)) return (NULL)
    vcfdatos <- read.vcfR(input$archivo$datapath)

 
    # Extraer la infomacion del archivo
    FIJO<- as.data.frame(vcfdatos@fix[,-8])
    INFOvcf<- extract_info_tidy (vcfdatos, info_types= FALSE)
    datos<- cbind(FIJO, INFOvcf)
  
  
    # Selecciono solo las entradas en las que el gen no es OMG
    
    datos <- datos [!(datos$Gene.refGene == "OMG"),]
    
    
    # Selecciono algunas columnas y creo datos2
    
    datos2<- select(datos, avsnp147, CHROM, REF, ALT, Gene.refGene, Func.refGene, ExonicFunc.refGene, 
                    ExAC_ALL, gnomAD_exome_ALL, CLNALLELEID, CLNSIG, SIFT_score, SIFT_pred, Polyphen2_HVAR_score, 
                    Polyphen2_HVAR_pred, MutationTaster_score, MutationTaster_pred, 
                    phyloP20way_mammalian, phastCons20way_mammalian, CLNDN)
    
    


    # Cambio el nombre de las columnas de datos2
    
    colnames(datos2)<- c("dbSNP", "CHR", "REF", "ALT", "Gen", "Funcion", "Tipo", "ExAC", "gnomAD_exome", "alleleID", "CLNSIG",
                          "SIFT_score", "SIFT_pred", "PolyPhen2_score", 
                         "PolyPhen2_pred", "MT_score", "MT_pred", "phyloP", "phastCons", "CLNDN")
    
    
    # Divido AAChange.refGene en columnas. Contiene el efecto en el DNA y proteina en cada transcrito separados por ","
    
     aachange<- strsplit (as.character(datos$AAChange.refGene), split = ",", fixed = FALSE)
     transcritoA<- sapply(aachange,  "[", 1)
     transcritoB<- sapply(aachange, "[", 2)
     
     
     

     # Incluyo las nuevas columnas al archivo de datos2
    
     datos2<- cbind(datos2, transcritoA, transcritoB)
    
      
     
     ######################## APLICACION DE FILTROS ################### 
     
     
     #### Seleccion del GEN ####

    # Creo el subset Resultado. Muestra solo las variantes exonicas de cambio de sentido, para cada gen 
   
     resultado <- datos2 [(datos2$Gen == input$gen & datos2$Funcion == "exonic" & datos2$Tipo == "nonsynonymous_SNV"), ]
      
     
      
    # Si se seleccionan todos los genes en input$gen:
      
    if (input$gen == "todos") {
      resultado<- datos2 [(datos2$Funcion == "exonic" & datos2$Tipo == "nonsynonymous_SNV"),]
    }
     
      

     #### seleccion de la FRECUENCIA POBLACIONAL en ExAC ####

    # Filtro Resultado, para filtrar segun la frecuencia poblacional elegida:  
     
  
     
     if (input$frec !="") {
       resultado <- resultado[(resultado$ExAC>1 | resultado$ExAC<=input$frec),]
     }
     
     
     
     
     
     
     #### Seleccion de la clasificacion en CLINVAR ####
     
     if (input$clinvar!="") {
       
       if (input$clinvar == "todos" )  { 
         resultado }
       
       if (input$clinvar == "no") {
         resultado <- resultado [(resultado$CLNSIG == "."),]}
       
       if (input$clinvar == "b" )  { 
         resultado <- resultado [(resultado$CLNSIG == "Benign"),]}
       
       if (input$clinvar == "lb" )  { 
         resultado <- resultado [(resultado$CLNSIG == "Likely_benign"),]}
       
       if (input$clinvar == "blb" )  { 
         resultado <- resultado [(resultado$CLNSIG == "Benign/Likely_benign"),]}
       
       if (input$clinvar == "us" )  { 
         resultado <- resultado [(resultado$CLNSIG == "Uncertain_significance"),]}
       
       if (input$clinvar == "lp" )  { 
         resultado <- resultado [(resultado$CLNSIG == "Likely_pahogenic"),]}
       
       if (input$clinvar == "plp" )  { 
         resultado <- resultado [(resultado$CLNSIG == "Pathogenic/Likely_pathogenic"),]}
       
       if (input$clinvar == "p" )  { 
         resultado <- resultado [(resultado$CLNSIG == "Pathogenic"),]}
       
         }
     
     
     
      
    # Res. Selecciono las columnas de Resultado que quiero que aparezcan:
    
     res<- select(resultado, dbSNP, CHR, REF, ALT, Gen, ExAC, gnomAD_exome, alleleID, CLNSIG, SIFT_score,
                  SIFT_pred, PolyPhen2_score, PolyPhen2_pred, MT_score, MT_pred, 
                  phyloP, phastCons, CLNDN, transcritoA, transcritoB)
     
    

     
     explicacion<- paste (
       ("Informacion mostrada en cada columna:
       dbSNP: codigo de la variante en dbSNP.
       CHR: cromosoma en el que se encuentra la variante.
       REF: alelo de referencia.
       ALT: alelo alternativo.
       Gen: gen en el que se encuentra.
       ExAC: frecuencia de la variante en ExAC en todas las poblaciones.
       gnomAD_exome: frecuencia de la variante en gnomAD exome en todas las poblaciones.
       alleleID: codigo de la variante en ClinVar.
       CLNSIG: significado clinico de la variante segun ClinVar.
       SIFT_score: puntuacion en IFT, de 0 a 0.05 deletereas, >0.05 toleradas.
       SIFT_pred: prediccion en SIFT, D (deletereas) o T (toleradas).
       PolyPhen2_score: puntuacion en PolyPhen2, de 0 a 0.446 benignas, >0.446 y <0.908 posiblemente patogenicas, >0.908 probablemente patogenicas.
       PolyPhen2_pred: prediccion en PolyPhen 2, B (benigna) o D (deleterea).
       MT_score: puntuacion en MutationTaster, probabilidad de que la clasificacion sea cierta, <0.5 es dudoso.
       MT_pred: prediccion en MutationTaster, D (probablemente perjudicial), A (descrita como patogenica), N (probablemente benigna) o P (polimorfismo, descrito como benigna).
       phyloP: puntuacion de regiones conservadas,de -14 a +6, valores positivos indican region conservada y valores negativos region poco conservada.
       phastCons: puntuacion de regiones conservadas, de 0 a 1, cuanto mas cercano a 1 mas probable es que la region sea conservada.
       CLNDN: Interpretacion de ClinVar, con que patologias se encuentra asociada. 
       transcritoA: informacion de uno de los transcritos posibles sobre el gen, nombre del transcrito, exon en el que aparece la variante, cambio en la secuencia de DNA y cambio en la secuencia de la proteina. 
       transcritoB: informacion de otro de los transcritos posibles sobre el gen, nombre del transcrito, exon en el que aparece la variante, cambio en la secuencia de DNA y cambio en la secuencia de la proteina. "

        
        
        ))
     
  

     # Output: tabla con los resultados 
     
  output$variantes<-renderDataTable(res)
  output$texto<- renderText ("Solo se muestran las variantes exonicas de cambio de sentido
                                                          ")
  output$explicacion<- renderText(explicacion)
 
  
  
})
  
})





