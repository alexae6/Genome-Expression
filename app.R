

## Setup
## load required packages and data import here
library(tidyverse)
library(magrittr)
library(dplyr)
library(shiny)
library(DT)
library(pheatmap)
library(data.table)

##### Commented out code Creates dataframes loaded into the app later
# Data for tab 1
#data_file <- "PRMT5.txt"
#if(file.exists(data_file)){
#  PRMT5 <- read_delim(file=data_file, 
#                      delim="\t", 
#                      col_types="cnnnnfc")
#}

#wdr77 <- "WDR77.txt"
#if(file.exists(wdr77)){
#  WDR77 <- read_delim(file=wdr77, 
#                      delim="\t", 
#                      col_types="cnnnnfc")
#}

# Data for tab 2
# Data Table of All Enrichment Results
#data_file <- "~/BCHM421/assignments/Final/go_enrichment.txt"
#if(file.exists(data_file)){
#  go_enrichment <- read_delim(file=data_file, 
#                              delim="\t")
#}

#data_file <- "~/BCHM421/assignments/Final/genes_wdown.txt"
#if(file.exists(data_file)){
#  genes_wdown <- read_delim(file=data_file, 
#                            delim="\t")
#}

#data_file <- "~/BCHM421/assignments/Final/genes_wup.txt"
#if(file.exists(data_file)){
#  genes_wup <- read_delim(file=data_file, 
#                          delim="\t")
#}
#data_file <- "~/BCHM421/assignments/Final/genes_pdown.txt"
#if(file.exists(data_file)){
#  genes_pdown <- read_delim(file=data_file, 
#                            delim="\t")
#}
#data_file <- "~/BCHM421/assignments/Final/genes_pup.txt"
#if(file.exists(data_file)){
#  genes_pup <- read_delim(file=data_file, 
#                          delim="\t")
#}

#save(PRMT5,WDR77,go_enrichment, genes_wdown, genes_pdown,genes_wup,genes_pup, file = "mydata.rda")


# Tab 1 Stuff
#tibble <- data.frame(PRMT5)
#prmt5 <- PRMT5 %>% mutate(.,Silenced = "PRMT")
#wdr77 <- WDR77 %>% mutate(.,Silenced = "WDR77")

# Tab 2 Stuff
# Create tibble with all gene information
#genes <- rbind(genes_pup,genes_pdown,genes_wup,genes_wdown)
#rna <- dplyr::union(wdr77,prmt5)
#symbol_FC <- rna %>% select(SYMBOL,Silenced,logFC)
#colnames(symbol_FC) <- c('Symbol','Silenced','logFC')

# Create tibble with all information
#symbol_FCP <- PRMT5 %>% select(SYMBOL,logFC)
#colnames(symbol_FCP) <- c('Symbol','logFC_PRMT')
#symbol_FCW <- WDR77 %>% select(SYMBOL,logFC)
#colnames(symbol_FCW) <- c('Symbol','logFC_WDR')
#all_info <- full_join(symbol_FCW,symbol_FCP, by='Symbol')

# Determine the biggest FC value between both silenced genes
#biggest <- all_info %>% group_by(Symbol) %>% mutate(maxFCvalue=max((abs(logFC_PRMT)),(abs(logFC_WDR))))

# Create tibble with all information in descending order
#hold1 <- full_join(biggest,symbol_FC)
#final_df <- full_join(hold1,genes, keep=FALSE) %>% unique()

# Eliminate repeats in ordering so make new tibble uisng modulus
#row_odd <- seq_len(nrow(final_df))%%2

#save(tibble,prmt5,wdr77,genes,rna,symbol_FC,symbol_FCP,symbol_FCW,all_info,biggest,hold1,final_df,row_odd,file ="tab12.rda")

# Tab 3 Stuff
#data_file <- "~/BCHM421/assignments/Final/cpm.txt"
#if(file.exists(data_file)){
#  cpm <- read_delim(file=data_file, 
#                    delim="\t")
#}

#data_file <- "~/BCHM421/assignments/Final/genes_pd.txt"
#if(file.exists(data_file)){
#  genes1 <- read_delim(file=data_file, 
#                    delim="\t")
#}

#data_file <- "~/BCHM421/assignments/Final/genes_pu.txt"
#if(file.exists(data_file)){
#  genes2 <- read_delim(file=data_file, 
#                       delim="\t")
#}

#data_file <- "~/BCHM421/assignments/Final/genes_wd.txt"
#if(file.exists(data_file)){
#  genes3 <- read_delim(file=data_file, 
#                       delim="\t")
#}

#data_file <- "~/BCHM421/assignments/Final/genes_wu.txt"
#if(file.exists(data_file)){
#  genes4 <- read_delim(file=data_file, 
#                       delim="\t")
#}

# Create tibble with all gene information
#allgenes <- rbind(genes_pup,genes_pdown,genes_wup,genes_wdown)


#colnames(cpm) <- c("Symbol","gfp_1","gfp_2","gfp_3","mep50_1","mep50_2","mep50_3","prmt5_1","prmt5_2","prmt5_3")
#cpm
#cpm2 <- full_join(cpm,allgenes)

#save(cpm,genes1,genes2,genes3,genes4,allgenes,cpm2, file="tab3stuff.rda")

###############################################################################
## Application

load('mydata.rda')
load('tab12.rda')
load('tab3stuff.rda')


## UI
ui <- navbarPage("Gene Explorer", 
                 tabPanel("Instructions", 
                          mainPanel(
                            h3("RNA seq Data"),
                            h5("This section displays the results of an RNAseq analysis. Data was collected for cells whose PRMT5 or WDR77 genes were silenced. You can explore the effects of these silenced genes through manipulating the RNAseq Results data table. If you wish for a visual representation, you can select gene symbols to be displayed in a column plot. Type or select your desired gene name(s) and you can see how their expression was affected in cells silenced for PRMT5 and WDR77."),
                            
                            h3("Gene Ontology "),
                            h5("A gene ontology enrichment analysis was performed using Gorilla. Data from genes that were up- and down-regulated in PRMT5 and WDR77 were combined and presented in the data table found on the Gene Ontology page. You can search for specific GO Terms (or GO_ID as presented in the table) and investigate the gene enrichment for different groups such as downregulated PRMT5 cells. Additionally, you can find the top number genes for a specific GO ID and have them displayed in a column plot. You can choose the maximum number of genes displayed and it will show you the genes and their expression in cells silenced for PRMT5 and WDR77."),
                          
                            h3("Heat Map"),
                            h5("This page displays a heat map with hierarchical clustering for genes of a selected GO Term. You can choose a GO Term to explore, and a heat map will be produced. You can also choose the distance method, linkage method, and whether to cluster by gene or sample. This heat map displays the log of CPM values for each sample.")
                            
                            )),
                 tabPanel("RNAseq Data", 
                          sidebarLayout(
                            sidebarPanel(
                              p(),
                              hr(),
                              selectInput('symbol',label=h3("Select gene symbols (start typing to find gene)"),choices="tibble$SYMBOL",selectize=TRUE, multiple=TRUE)
                            ),
                             mainPanel(
                              h3("Column Plot of Select Genes"),
                              plotOutput("table2"),
                              h3("RNAseq Results"),
                              DT::dataTableOutput(outputId = "table")
                            ))), 
                 tabPanel("Gene Ontology Data", 
                          sidebarLayout(
                            sidebarPanel(
                              p("Select a GO ID (start typing to find a GO ID)"),
                              hr(),
                              selectInput('ID',label=NULL,choices="go_enrichment",selectize=TRUE, multiple=FALSE),
                              sliderInput("integer", "Integer:",min = 1, max = 50,value=1)
                             
                            ),
                            mainPanel(
                              h3("Column Plot for GO Category"),
                              plotOutput("table3"),
                              h3("Data Table of GO Enrichment Results"),
                              DT::dataTableOutput(outputId = "GO")
        
                            ))),
                 tabPanel("Heat Map", 
                          sidebarLayout(
                            sidebarPanel(
                              p("Select a GO ID (start typing to find a GO ID)"),
                              hr(),
                              selectInput("tab3",label=NULL,choices="go_enrichment",selectize=TRUE, multiple=TRUE),
                              selectInput("TF",label=h3("Cluster by Sample?"),choices = list("TRUE", "FALSE"),selected = NULL,multiple = FALSE),
                              selectInput("TF2",label=h3("Cluster by Gene?"),choices = list("TRUE", "FALSE"),selected = NULL,multiple = FALSE),
                              selectInput("distance1",label=h3("Enter distance method for samples"),choices = list("euclidean", "manhattan","correlation"),selected = NULL,multiple = FALSE),
                              selectInput("distance2",label=h3("Enter distance method for genes"),choices = list("euclidean", "manhattan","correlation"),selected = NULL),
                              selectInput("linkage",label=h3("Enter linkage method for all clustering"),choices = list("average", "complete","single"),selected = NULL),
                              selectInput("scale",label=h3("Enter scale method for all clustering"),choices = list("column", "row","none"),selected = NULL)
                            ), 
                            mainPanel(
                              plotOutput("table4")
                              #DT::dataTableOutput("testtable")
                              
                          
                            )))
          )

## Server

server <- function(input, output, session) {
  
  # RNAseq Data Tab- allows user to select a gene by typing or using a drop down menu
  # Displays a bar plot showing gene expression amount and providing a table of detailed gene expression
  
  output$table <- DT::renderDataTable({
    datatable({rna},filter = 'top',extensions = c('Buttons'),
              options = list(
                dom = 'Blfrtip',
                buttons = list(list(extend = 'csv', filename= 'RNAseq_Explorer'),
                               list(extend = 'excel', filename= 'RNAseq_Explorer'))
                )
              
           ) %>% formatRound(columns=c('logFC', 'logCPM',"PValue","FDR"), digits=4)
  })
  output$symbol <- renderPrint(input$symbol)
  a <- reactive({
    req(input$symbol)
    p <- prmt5 %>% filter(SYMBOL %in% input$symbol)
    w <- wdr77 %>% filter(SYMBOL %in% input$symbol)
    a <- rbind(p,w)
    })
  output$table2 <- renderPlot({
   ggplot(a(),aes(x=SYMBOL,y=logFC,fill=Silenced))+geom_bar(stat="identity", position="dodge")+labs(title=paste("Results for Select Genes"))
  })
 observe({
  updateSelectInput(session,"symbol",choices=tibble$SYMBOL)
 }) 

 
 # Gene Ontology Data Tab- allows user to select a gene by typing or using a drop down menu
 # Uses a slider to say how many genes to display in bar graph displaying gene expression amount-
 # Sorted in ascending order
 
 output$ID <- renderPrint(input$ID)
 observe({
   updateSelectInput(session,"ID",choices=go_enrichment$GO_ID)
 }) 
 output$GO <- DT::renderDataTable({
   datatable({go_enrichment},filter = 'top',extensions = c('Buttons'),
             options = list(
               dom = 'Blfrtip',
               lengthMenu = c(5,10,25,50,100),
               buttons = list(list(extend = 'csv', filename= 'Geo_Explorer'),
                              list(extend = 'excel', filename= 'Geo_Explorer'))
             )
   )%>% formatRound(columns=c("FDR","Enrichment"), digits=4)  
 })
 
 output$ID <- renderPrint(input$ID)
 b <- reactive({
   req(input$ID)
   holder2 <- arrange(final_df,desc(maxFCvalue)) %>% filter(GO_ID %in% input$ID) %>% filter(row_number() %% 2 == 1) %>% head(n=input$integer)
   holder_p <- holder2 %>% select(Symbol,logFC_PRMT) %>% mutate(Silenced = "PRMT")
   colnames(holder_p) <- c("Symbol","logFC","Silenced")
   holder_w <- holder2 %>% select(Symbol,logFC_WDR) %>% mutate(Silenced = "WDR")
   colnames(holder_w) <- c("Symbol","logFC","Silenced")
   full_join(holder_p,holder_w)
 })
 
 #Plot
 output$table3 <- renderPlot({
   ggplot(b(),aes(x=Symbol, y= logFC, fill = Silenced))+ 
     geom_bar(stat="identity",position = "dodge2")+
     coord_flip()+labs(title=paste("Results for Top",input$integer,"Genes for",input$ID))
   
 })

 
 # Heatmap Tab- Allows user to select the gene and different clustering parameters to create
 # and displays a heatmap linking gene similarity
 
 
output$tab3 <- renderPrint(input$tab3)
observe({
   updateSelectInput(session,"tab3",choices=go_enrichment$GO_ID)
 }) 
output$TF <- renderPrint({})
output$TF2 <- renderPrint({})
output$distance1 <- renderPrint({})
output$distance2 <- renderPrint({})
output$linkage <- renderPrint({})
output$scale <- renderPrint({})

d <- reactive({
   cpm_values <- cpm2 %>% filter(GO_ID %in% input$tab3) %>% select(Symbol:prmt5_3)
   cpm_mat <- cpm_values %>% select(gfp_1:prmt5_3) %>% as.matrix() %>% na.omit
   #rownames(cpm_mat) <- cpm_values %>% .[,1]

 # dimnames(cpm_mat) <- NULL
 # colnames(cpm_mat) <- c("gfp_1","gfp_2","gfp_3","mep50_1","mep50_2","mep50_3","prmt5_1","prmt5_2","prmt5_3")
   # labels <- cpm2 %>% filter(GO_ID  %in% input$tab3)  %>% .[,1] %>% list()

   # labels <- cpm2 %>% filter(GO_ID  %in% input$tab3)  %>% .[,1] %>% list()
    #rownames(d()) <- labels
  })
 # rownames(cpm_mat) <- labels
  
 # cpm3 <- left_join(select_GO, cpm)
 # cpm_mat <- cpm3 %>% select(2:10) %>% as.matrix()  
  #rownames(cpm_mat) <- cpm3  %>% pull(gene)
  
  #Note: if interactive clustering method & clustering distance is up to input. If clustering by gene cluster_row=TRUE, if clustering by sample cluster_col = TRUE
 # pheatmap(cpm_mat, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", cluster_rows=TRUE, cluster_cols =TRUE, show_rownames = TRUE, clustering_method="average")
  
  
 # cpm_mat <- cpm2  %>% filter(GO_ID %in% input$ID2) %>% .[2:10,] %>% as.numeric() %>%
    #as.matrix() 
 # rownames(cpm_mat) <- cpm2 %>% filter(GO_ID %in% input$ID2) %>% dplyr::pull(Symbol)
  
 


output$testtable <- DT::renderDataTable({
  datatable({d()})
})


output$table4 <- renderPlot({
  pheatmap(d(),
           scale=input$scale,
           ColV = input$TF,
           RowV =input$TF2,
           cluster_distance_cols = input$distance1,
           clustering_distance_rows = input$distance2,  
           clustering_method=input$linkage
  
)}) 
#res = 150, 
#height = function(){
#  max(1200, length(d()) * 12)
#}, 
#width = 1200
#)

 
  }


# Run the application 

shinyApp(ui = ui, server = server)