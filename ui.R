if (!require("shinythemes"))
install.packages("shinythemes")
library(shinythemes)


shinyUI (fluidPage(
  theme = shinytheme ("sandstone"),
  
  titlePanel ("FiltVar"),
  
  sidebarLayout(
    sidebarPanel( 
      title = "Cargar el archivo",
      
      # Cargar el archivo VCF de entrada
      fileInput("archivo", "Seleccionar archivo VCF", multiple = FALSE, 
                accept = c("text", ".vcf")),
      
      # Seleccion del gen o genes de interes
      selectizeInput("gen", label= "Seleccion del gen", 
                     choices = c("Todos" = "todos", "A2ML1" = "A2ML1", "BRAF" = "BRAF", "CBL" = "CBL", "HRAS" = "HRAS", 
                                 "KRAS" = "KRAS", "LZTR1" = "LZTR1", "MAP2K1" = "MAP2K1", "MAP2K2" = "MAP2K2", 
                                 "MRAS" = "MRAS", "NF1" = "NF1", "NRAS" = "NRAS", "PPP1CB" = "PPP1CB", 
                                 "PTPN11" = "PTPN11", "RAF1" = "RAF1", "RASA1" = "RASA1", "RASA2"= "RASA2", 
                                 "RIT1" = "RIT1", "RRAS" = "RRAS", "SHOC2"= "SHOC2", "SOS1" = "SOS1", "SOS2" = "SOS2")), 
      
      
      # Seleccionar la frecuencia en ExAC
      sliderInput("frec", label= "Frecuencia en ExAC menor o igual a:", value = 0.05, min=0, max=1, step=0.05),
      
      
      # Seleccion del significado clinico
      selectizeInput("clinvar", label= "Significado clinico", 
                     choices = c("Todos" = "todos", "No reportado" = "no", "Patogenica" = "p", "Probablemente patogenica" = "lp", 
                                 "Patogenica/ probablemente patogenica" ="plp", "Significado incierto" = "us", 
                                 "Probablemente benigna" = "lb", "Benigna/ probablemente benigna"= "blb", "Benigna" = "b"
                                 )),
    

      # Boton para subir el archivo
      
     actionButton ("aplicar", "Aplicar los filtros seleccionados" )),
    
    
    
    mainPanel(
      
      textOutput("texto"), 
      dataTableOutput("variantes")
      
      
    )
    
  )
) )
  
  
  
