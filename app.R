setwd('/Users/mrh/Desktop/ncsu_project/data')
library(dplyr)
library(matrixStats)
library(ggplot2)
library(stringr)
library(shiny)
library(plotly)
library(DT)
library(tidyr)
library(tidyverse)

df<- read.csv('processed_data_complete.csv')

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      checkboxGroupInput("substrate1", "Select Substrate 1", choices = names(df[4:22]), selected = names(df[4:6])),
      checkboxGroupInput("substrate2", "Select Substrate 2", choices = names(df[4:22]), selected = names(df[7:9]))
    ),
    mainPanel(
      tabsetPanel(
        tabPanel( 'Proteome',
          plotly::plotlyOutput("volcano_plot_1"),
          plotly::plotlyOutput("volcano_clust_1"),
          numericInput("percent_protein", "Protein NSAF Must Be Greater Than", value = 0.00, min = 0, max = 100),
          sliderInput("gene_range", "Adjust Target Gene Range", min = 1, max = max(df$gene, na.rm = T), value = c(1, max(df$gene, na.rm = T))),
          numericInput('clusters', 'Gene Clusters', 3, min = 1, max = max(df$gene, na.rm = T))
        ),
        tabPanel('Proteome table',
          dataTableOutput('working_dataframe')
        ),
        
        tabPanel('Specific Protein',
          numericInput("geneID", "Enter Gene ID", value = 5533, min = 1, max = max(df$gene, na.rm = T)),
          checkboxGroupInput('gene_comp', 'Select Substrates', choices = names(df[23:29]), selected = names(df[23:29])),
          plotly::plotlyOutput("gene_comparison"),
          plotly::plotlyOutput("gene_heat_map")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  t_test_fun <- function(x) {
    stop1 <- 1+length(input$substrate1)
    start2 <- stop1+1
    stop2 <- start2 + length(input$substrate2)
    t <- t.test(
     as.numeric(as.character(x[2:stop1])),
     as.numeric(as.character(x[start2:stop2])))
     return(t$p.value)
  }
  # data 1 contains only data on the two substrates you are comparing
  data1 <- reactive({
    data <- df%>%
      select(gene, input$substrate1, input$substrate2)%>% 
      mutate(sub1_mean = rowMeans(cbind(df[input$substrate1]), na.rm = T))%>%
      mutate(sub2_mean = rowMeans(cbind(df[input$substrate2]), na.rm = T))
    
    return(data)
  })
  
  cut_sub1 <- reactive({sum(data1()[,'sub1_mean'], na.rm = T) * (input$percent_protein / 100)})
  cut_sub2 <- reactive({sum(data1()[,'sub2_mean'], na.rm = T) * (input$percent_protein / 100)})
  # data 2 is a trimmed dataframe that only shows poteins that exceed a certain threshold of total protein
  data2 <- reactive({
    data_2 <- data1()%>%
      filter((sub1_mean > cut_sub1())&(sub2_mean > cut_sub2()))%>%
      filter(gene >= input$gene_range[1] & gene <= input$gene_range[2])%>%
      mutate(change = sub2_mean - sub1_mean)%>%
      mutate(move = ifelse(change > 0, 'increase', 'decrease'))
    return(data_2)
  })
  
 
  p <- reactive({apply(data2(), 1, t_test_fun)})
 
  
  data3 <- reactive({
    data_3 <- data2()%>%
    mutate(pvalues = p())%>%
    mutate(neglog10pvalues = -log(p()))
    return(data_3)
    })
  
  volcano_1 <- function() {
    ggplot(data3(), aes(x = change, y = neglog10pvalues, color = gene))+
      geom_point()+
      labs(title = sprintf("Proteom Change from %s to %s", input$substrate1[1], input$substrate2[1]), x = "Change (NSAF)", y = "-log10(p-value)")+
      theme(plot.title = element_text(hjust = 0.5))
  }
  
  cluster_setup <- reactive({
    data_clust <- data3()%>%
    select(gene)
    return(data_clust)
  })
  
  clusters <- reactive({
    kmeans(cluster_setup(), input$clusters)
  })
  
  
  volcano_clust <- function() {
    ggplot(data3(), aes(x = change, y = neglog10pvalues, color = as.factor(clusters()$cluster)))+
      geom_point(aes(text = data3()$gene))+
      labs( x = "Change (NSAF)", y = "-log10(p-value)", color = 'Cluster')
  }
  
  comp_gene <- reactive({
    data_gene <- df%>%
      filter(gene == input$geneID)%>%
      select(input$gene_comp)
    datalong <- pivot_longer(data_gene, cols = colnames(data_gene), names_to = "name")
    return(datalong)
  })
  
  heatmap <- reactive({
    data_gene2 <- df%>%
      filter(gene == input$geneID)%>%
      select(c(1,4:22))
    
    datalong2 <- pivot_longer(data_gene2, cols = colnames(data_gene2), names_to = 'name')
    datalong2['substrate'] <- 'unknown'
    datalong2 <- datalong2[2:20, 2:3]
    datalong2[1:3, 2] <- "Acetoin"
    datalong2[4:6, 2] <- 'C2B'
    datalong2[7:9, 2] <- 'EtOH'
    datalong2[10:12, 2] <- 'HIBA'
    datalong2[13:15, 2] <- 'T2B'
    datalong2[16:17, 2] <- 'Tetradecane'
    datalong2[18:19, 2] <- 'n.Heptane'
    attach(datalong2)
    comparison <- pairwise.t.test(value, substrate, p.adjust.method = "bonf")
    comp <- as.data.frame(comparison$p.value)
    comp <- comp%>%
      rownames_to_column()%>%
      gather(colname, value, -rowname)
    return(comp)
  })
  
  gene_heat_plot <- function() {
    ggplot(heatmap(), aes(x = rowname, y = colname, fill = value))+
      geom_tile()+
      labs(title = "Comparative P-values Between NSAF Values", x = "", y = "")+
      theme(plot.title = element_text(hjust = 0.5))
  }

  
  gene_comp_plot <- function() {
    ggplot(comp_gene())+
      geom_col(aes(x= name, y = value, color = name, fill = name))+
      labs(title = sprintf('Average NSAF for Gene %s', input$geneID), x = 'Substrate', y = 'Average NSAF')+
      theme(axis.text.x = element_text(angle = 45))
    }
  
  output$working_dataframe <- renderDataTable(data3())
  
  output$volcano_plot_1 <- plotly::renderPlotly({volcano_1()})
  
  output$volcano_clust_1 <- plotly::renderPlotly({volcano_clust()})
  
  output$gene_comparison <- plotly::renderPlotly({gene_comp_plot()})
  
  output$gene_heat_map <- plotly::renderPlotly({gene_heat_plot()})
  
}

shinyApp(ui = ui, server = server)
