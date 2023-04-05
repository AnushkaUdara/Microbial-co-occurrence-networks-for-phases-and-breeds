# source("C:\\Users\\Anushka\\Desktop\\Anushka\\Research\\Methodology\\my_app\\TestShiny.R")
# 
# #Print User Input
# #install.packages("shiny")
# library(shiny)
# #install.packages("shinythemes")
# library(shinythemes)
# #library(shinyjs)
# 
# #jscode <- "shinyjs.closeWindow = function() { window.close(); }"
# 
# # Loading required packages
# library(tidyverse)
# library(phyloseq)
# library(qiime2R)
# library(BiocManager)
# library(microbiome)
# library(NetCoMi)
# library(WGCNA)
# 
# setwd("C:\\Users\\Anushka\\Desktop\\Anushka\\Research\\Methodology\\my_app")
# 
# #Define UI
# ui<-fluidPage(
#   theme = shinythemes::shinytheme("cerulean"),
#   titlePanel(h1(strong("Microbial Co-Occurrence Networks"))),
# 
#   sidebarLayout(
#     position = "left",
#     sidebarPanel(h2(strong("Construct Networks")),
#                  br(),
#                  textInput(inputId="Me",label="Network Construction Method"),
#                  textInput(inputId="Pa",label="Select Variable"),
#                  textInput(inputId="Sp",label="Select Sub-Variable"),
# 
#                  #titlePanel("Construct Network"),
#                  #useShinyjs(),
#                  #extendShinyjs(text = jscode, functions = c("closeWindow")),
#                  actionButton("script", "Build Network"),
#                  br(),br(),
# 
# 
#                  h2(strong("Visualize Networks")),
# 
#                  radioButtons(inputId="Taxa",label="Select Taxa",
#                               choices=c("Genus"="Genus","Family"="Family","Order"="Order","Class"="Class")),
# 
#                  radioButtons(inputId="Folder2",label="Select Image Type for Cattle Breed",
#                                     choices=c("R(All)"="I_R_All","R(Top100)"="I_R_Top100","Cytoscape(All)"="I_Cytoscape_All","Cytoscape(Top100)"="I_Cytoscape_Top100")),
# 
#                  selectInput(inputId = "Method2",label="Select Network Construction Method for Cattle Breed",
#                              choices=c("SparCC"="Sparcc","SPRING"="Spring")),
# 
#                  selectInput(inputId = "CattleBreed",label="Select Cattle Breed",
#                              choices=c("FJC","FC","SC","JSC","BC","FSC","AFS","JC","S","A","J"),selected = NULL),
# 
#                  radioButtons(inputId="Folder1",label="Select Image Type for Lactation Phase",
#                               choices=c("R(All)"="I_R_All","R(Top100)"="I_R_Top100","Cytoscape(All)"="I_Cytoscape_All","Cytoscape(Top100)"="I_Cytoscape_Top100")),
# 
#                  selectInput(inputId = "Method1",label="Select Network Construction Method for Lactation Phase",
#                              choices=c("SparCC"="Sparcc","SPRING"="Spring","SPIEC-EASI"="SPIEC-EASI","CCLasso"="CCLasso")),
# 
#                  selectInput(inputId = "LactationPhase",label="Select Lactation Phase",
#                              choices=c("Early Phase"="Early","Mid Phase"="Mid","LatePhase"="Late","Dry Phase"="Dry","All"="All"),selected = NULL),
# 
#                  textOutput("result")
#                  ),
# 
#     mainPanel(
#               imageOutput(outputId = "network2"),
#               br(),br(),br(),br(),br(),br(),br(),br(),
#               imageOutput(outputId = "network1"),
#               br(),br(),br(),br(),br(),br(),br(),br(),
#               imageOutput(outputId="network3")
#               )
#   )
# 
# 
# )
# 
# #Define server
# server<-function(input,output){
# 
#   output$network1 <- renderImage({
#     filename <- normalizePath(file.path('./',paste(input$Folder1,"/",input$Taxa,"_",input$Method1,"_Phase_",input$LactationPhase,'.png', sep='')))
#     outputArgs=list(src = filename)
#   },deleteFile=FALSE)
# 
#   output$network2 <- renderImage({
#     filename <- normalizePath(file.path('./',paste(input$Folder2,"/",input$Taxa,"_",input$Method2,"_Breed_",input$CattleBreed,'.png', sep='')))
#     outputArgs=list(src = filename)
#   },deleteFile=FALSE)
# 
#   output$network3 <- renderImage({
#     filename <- normalizePath(file.path('./Images3/res.png'))
#     outputArgs=list(src = filename)
#   },deleteFile=TRUE)
# 
# 
#   mylist <- reactiveVal() # we will store the inputs in a reactive list
# 
#   observe({ # create the list
#     mylist(list(
#       Me=input$Me,
#       Pa=input$Pa,
#       Sp=input$Sp))
#   })
# 
#   observeEvent(input$script, { # "runScript" is an action button
#     source("TestShiny.R", local = list2env(mylist()))
#   })
# 
# 
# }
# 
# #Run the app
# shinyApp(ui=ui,server=server)


# runExample("01_hello")      # a histogram
# runExample("02_text")       # tables and data frames
# runExample("03_reactivity") # a reactive expression
# runExample("04_mpg")        # global variables
# runExample("05_sliders")    # slider bars
# runExample("06_tabsets")    # tabbed panels
# runExample("07_widgets")    # help text and submit buttons
# runExample("08_html")       # Shiny app built from HTML
# runExample("09_upload")     # file upload wizard
# runExample("10_download")   # file download wizard
# runExample("11_timer")      # an automated timer
# 


#source("C:\\Users\\Anushka\\Desktop\\Anushka\\Research\\Methodology\\my_app\\TestShiny.R")

#Print User Input
#install.packages("shiny")
library(shiny)
#install.packages("shinythemes")
library(shinythemes)
#library(shinyjs)

#jscode <- "shinyjs.closeWindow = function() { window.close(); }"

# Loading required packages
library(tidyverse)
library(phyloseq)
library(qiime2R)
library(BiocManager)
library(microbiome)
library(NetCoMi)
library(WGCNA)

#shiny::runGitHub("shiny-phyloseq","joey711")

setwd("C:\\Users\\Anushka\\Desktop\\Anushka\\Research\\Methodology\\my_app")
#install.packages("shinydashboard")
library(shinydashboard)




ui <- dashboardPage(
  dashboardHeader(title = "MiCoNeV"),
  
  dashboardSidebar(sidebarMenu(
    menuItem("Visualize", tabName = "Visualize", icon = icon("circle-nodes")),
    menuItem("Associations", tabName = "Associations", icon = icon("th")),
    menuItem("Construct", tabName = "Construct", icon = icon("wrench"))
    
  )),
  
  dashboardBody(tabItems(
    
    
    
    # First tab content
    tabItem(tabName = "Construct",
            fluidRow(
              box(
                h2(strong("Construct Networks")),
                br(),
                textInput(inputId="Me",label="Network Construction Method"),
                textInput(inputId="Pa",label="Select Variable"),
                textInput(inputId="Sp",label="Select Sub-Variable"),
                
                #titlePanel("Construct Network"),
                #useShinyjs(),
                #extendShinyjs(text = jscode, functions = c("closeWindow")),
                actionButton("script", "Build Network"),
                br(),br()
              )
              
              
            )
          ),
    
    
    # Third tab content
    tabItem(tabName="Associations",
            h2(strong("Edge Data")),
            
            fluidRow(
              box(
                title = "Taxonomic Level", background = "navy", solidHeader = TRUE,
                radioButtons(inputId="Taxa12",label="Select Taxa",
                             choices=c("Genus"="Genus","Family"="Family","Order"="Order","Class"="Class"))
              )
            ),
            
            fluidRow(
              box(
                title = "Cattle Breed", background = "orange", solidHeader = TRUE,
                radioButtons(inputId="TaxaNum12",label="Select Taxa Number for Cattle Breed",
                             choices=c("All"="All","Top100"="Top100")),
                
                selectInput(inputId = "Method12",label="Select Network Construction Method for Cattle Breed",
                            choices=c("SparCC"="Sparcc","SPRING"="Spring")),
                
                selectInput(inputId = "CattleBreed12",label="Select Cattle Breed",
                            choices=c("FJC","FC","SC","JSC","BC","FSC","AFS","JC","S","A","J"),selected = NULL),
                
                actionButton(inputId = "go12", 
                             label = "Update")
                
              )
              
            ),
            
            fluidRow(
              box(
                tableOutput("contents12"),
                width=1000,
                style = "overflow-x: scroll"
              )

            )
            
            
            
            ),
    
   

    # Second tab content
    tabItem(tabName = "Visualize",
            h2(strong("Visualize Networks")),
            
            fluidRow(
              box(
                title = "Taxonomic Level", solidHeader = TRUE,
                radioButtons(inputId="Taxa",label="Select Taxa",
                             choices=c("Genus"="Genus","Family"="Family","Order"="Order","Class"="Class"))
              )
            ),
            
            
            fluidRow(
              
              # box(
              #   title = "Cattle Breed", background = "maroon", solidHeader = TRUE,
              #   radioButtons(inputId="Folder2",label="Select Image Type for Cattle Breed",
              #                choices=c("R(All)"="I_R_All","R(Top100)"="I_R_Top100","Cytoscape(All)"="I_Cytoscape_All","Cytoscape(Top100)"="I_Cytoscape_Top100")),
              #   
              #   selectInput(inputId = "Method2",label="Select Network Construction Method for Cattle Breed",
              #               choices=c("SparCC"="Sparcc","SPRING"="Spring")),
              #   
              #   selectInput(inputId = "CattleBreed",label="Select Cattle Breed",
              #               choices=c("FJC","FC","SC","JSC","BC","FSC","AFS","JC","S","A","J"),selected = NULL),
              #   
              #   actionButton(inputId = "go2", 
              #                label = "Update") 
              #   
              #   
              # ),
              # 
              # box(
              #   title = "Lactation Phase", background = "blue", solidHeader = TRUE,
              #   radioButtons(inputId="Folder1",label="Select Image Type for Lactation Phase",
              #                choices=c("R(All)"="I_R_All","R(Top100)"="I_R_Top100","Cytoscape(All)"="I_Cytoscape_All","Cytoscape(Top100)"="I_Cytoscape_Top100")),
              #   
              #   selectInput(inputId = "Method1",label="Select Network Construction Method for Lactation Phase",
              #               choices=c("SparCC"="Sparcc","SPRING"="Spring","SPIEC-EASI"="SPIEC-EASI","CCLasso"="CCLasso")),
              #   
              #   selectInput(inputId = "LactationPhase",label="Select Lactation Phase",
              #               choices=c("Early Phase"="Early","Mid Phase"="Mid","LatePhase"="Late","Dry Phase"="Dry","All"="All"),selected = NULL),
              #   
              #   #textOutput("result")
              #   
              #   actionButton(inputId = "go1", 
              #                label = "Update") 
              
              
              tabBox(
                title="Select Variable",
                id="tabset1",
                width="100px",
                side="right",
                tabPanel(
                  "Breed",
                  #title = "Cattle Breed", background = "maroon", solidHeader = TRUE,
                  radioButtons(inputId="Folder2",label="Select Image Type for Cattle Breed",
                               choices=c("R(All)"="I_R_All","R(Top100)"="I_R_Top100","Cytoscape(All)"="I_Cytoscape_All","Cytoscape(Top100)"="I_Cytoscape_Top100")),
                  
                  selectInput(inputId = "Method2",label="Select Network Construction Method for Cattle Breed",
                              choices=c("SparCC"="Sparcc","SPRING"="Spring")),
                  
                  selectInput(inputId = "CattleBreed",label="Select Cattle Breed",
                              choices=c("FJC","FC","SC","JSC","BC","FSC","AFS","JC","S","A","J"),selected = NULL),
                  
                  # actionButton(inputId = "go2",
                  #              label = "Update")
                      
                    
                ),
                
                tabPanel(
                  "Phase",
                  radioButtons(inputId="Folder1",label="Select Image Type for Lactation Phase",
                               choices=c("R(All)"="I_R_All","R(Top100)"="I_R_Top100","Cytoscape(All)"="I_Cytoscape_All","Cytoscape(Top100)"="I_Cytoscape_Top100")),

                  selectInput(inputId = "Method1",label="Select Network Construction Method for Lactation Phase",
                              choices=c("SparCC"="Sparcc","SPRING"="Spring","SPIEC-EASI"="SPIEC-EASI","CCLasso"="CCLasso")),

                  selectInput(inputId = "LactationPhase",label="Select Lactation Phase",
                              choices=c("Early Phase"="Early","Mid Phase"="Mid","LatePhase"="Late","Dry Phase"="Dry","All"="All"),selected = NULL),

                                 #textOutput("result")

                  # actionButton(inputId = "go1",
                  #              label = "Update")


                  
                  
                )
                
                
              )
              
              
                
              
              
            ),
            
            
            fluidRow(
              
              box(
                width=12,
                imageOutput(outputId = "network1"),
                br(),br(),br(),br(),
                style = "overflow-x: scroll"
              )
                
                
              
            ),
            
            # fluidRow(
            # 
            # 
            #   imageOutput(outputId = "network2"),
            #   br(),br(),br(),br(),br(),br()
            # 
            # ),
            
            
            

                             

                             

                             
    )
  )),
  
  tags$head(tags$style(HTML('
                        

                                /* body */
                                .content-wrapper, .right-side {
                                background-color: #ffffff;
                                }
                                
                                ')))
  
  
  
)

server <- function(input, output) {
  
  
      # data1 <- eventReactive(input$go1, { 
      #   filename <- normalizePath(file.path('./',paste(input$Folder1,"/",input$Taxa,"_",input$Method1,"_Phase_",input$LactationPhase,'.png', sep='')))
      #   outputArgs=list(src = filename) 
      # })
      # 
      # data2 <- eventReactive(input$go2, { 
      #   filename <- normalizePath(file.path('./',paste(input$Folder2,"/",input$Taxa,"_",input$Method2,"_Breed_",input$CattleBreed,'.png', sep='')))
      #   outputArgs=list(src = filename) 
      # })
      
      data12 <- eventReactive(input$go12, { 
        filename <- normalizePath(file.path('./Associations',paste(input$Taxa12,"_",input$Method12,"_Breed_",input$CattleBreed12,"_",input$TaxaNum12,'.csv', sep='')))
        #outputArgs=list(src = filename)
        df <- read.csv(filename,sep=",")
      })
      
      
      output$contents12<-renderTable({
        
        data12()
        
      })
      
  
      output$network1 <- renderImage({
        
        if(input$tabset1 == "Phase") {
          filename <- normalizePath(file.path('./',paste(input$Folder1,"/",input$Taxa,"_",input$Method1,"_",input$tabset1,"_",input$LactationPhase,'.png', sep='')))
        }
        else {
          filename <- normalizePath(file.path('./',paste(input$Folder2,"/",input$Taxa,"_",input$Method2,"_",input$tabset1,"_",input$CattleBreed,'.png', sep='')))
        }

        outputArgs=list(src = filename) 
      },deleteFile=FALSE)

      # output$network2 <- renderImage({
      #   filename <- normalizePath(file.path('./',paste(input$Folder2,"/",input$Taxa,"_",input$Method2,"_",input$tabset1,"_",input$CattleBreed,'.png', sep='')))
      #   outputArgs=list(src = filename)
      # },deleteFile=FALSE)

      output$network3 <- renderImage({
        filename <- normalizePath(file.path('./Images3/res.png'))
        outputArgs=list(src = filename)
      },deleteFile=TRUE)


      mylist <- reactiveVal() # we will store the inputs in a reactive list

      observe({ # create the list
        mylist(list(
          Me=input$Me,
          Pa=input$Pa,
          Sp=input$Sp))
      })

      observeEvent(input$script, { # "runScript" is an action button
        source("TestShiny.R", local = list2env(mylist()))
      })


  
  
}









#Run the app
shinyApp(ui=ui,server=server)