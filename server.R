
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
                    phyloP20way_mammalian, phastCons20way_mammalian)
    
    


    # Cambio el nombre de las columnas de datos2
    
    colnames(datos2)<- c("rsID", "CHR", "REF", "ALT", "Gen", "Funcion", "Tipo", "ExAC_ALL", "gnomAD_exome", "alleleID", "CLNSIG",
                          "SIFT_csore", "SIFT_pred", "PolyPhen2_score", 
                         "PolyPhen2_pred", "MT_score", "MT_pred", "phyloP", "phastCons")
    
    
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
    
    if (input$frec !=""){
      resultado<- resultado [(resultado$ExAC_ALL <= input$frec),]}
           

     
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
    
     res<- select(resultado, rsID, CHR, REF, ALT, Gen, ExAC_ALL, gnomAD_exome, alleleID, CLNSIG, SIFT_csore,
                  SIFT_pred, PolyPhen2_score, PolyPhen2_pred, MT_score, MT_pred, 
                  phyloP, phastCons, transcritoA, transcritoB)
     
    


     # Output: tabla con los resultados 
     
  output$variantes<-renderDataTable(res)
  output$texto<- renderText ("Solo se muestran las variantes exónicas de cambio de sentido")

  
  
})
  
})





