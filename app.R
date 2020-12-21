setwd('/Users/mrh/Desktop/ncsu_project/data')
library(dplyr)
library(matrixStats)
library(ggplot2)
library(stringr)
library(shiny)
library(plotly)
library(DT)

df<- read.csv('processed_data_complete.csv')

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      checkboxGroupInput("substrate1", "Select Substrate 1", choices = names(df[4:22]), selected = names(df[4:6])),
      checkboxGroupInput("substrate2", "Select Substrate 2", choices = names(df[4:22]), selected = names(df[7:9])),
      numericInput("percent_protein", "Input Percent Total Protein to Display", value = 0.01, min = 0, max = 100),
    ),
    mainPanel(
      dataTableOutput('working_dataframe'),
      plotly::plotlyOutput("volcano_plot_1")
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
      select(gene, input$substrate1, input$substrate2)%>% # added gene name here
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
      xlab("Change (NSAF)")+
      ylab("-log10(p-value)")
  }
  
  
  output$working_dataframe <- renderDataTable(data3())
  
  output$volcano_plot_1 <- plotly::renderPlotly({volcano_1()}) # add function to create plot
  
}

shinyApp(ui = ui, server = server)



# next steps
# make a gene range selector
# make a table that compares the difference in specific protein across all substrates 