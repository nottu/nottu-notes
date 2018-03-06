library(shiny)
library(moments)
library(entropy)

dists <- c('Uniforme', 'Normal', 'Exponential')

ui <- fluidPage(
  titlePanel('Kurtosis y Negentropía'),
  sidebarPanel(
    selectInput('xcol', 'Distribucion X', dists),
    numericInput("nvars", "Número de Variables", 1, min=1, max=1000, step=1),

    tags$hr(),
    textOutput('kurt'),
    textOutput('negentropia')
    ),
  mainPanel(
    plotOutput('plot1')
    )
  )

server <- function(input, output) {
  nv = 100
  xdist <- reactive({
    y <- seq(0, 0, length=nv)
    switch( input$xcol,
      "Normal" = {
        for (i in 1:input$nvars) 
          y = y + rnorm(nv)
      },
      "Exponential" = {
        for (i in 1:input$nvars) 
          y = y + rexp(nv)
      },
      "Uniforme" = {
        for (i in 1:input$nvars) 
          y = y + runif(nv)
      }
      )
    return(y/input$nvars)
  })

  output$plot1 <- renderPlot({
    xvals <-xdist()
    hist(xvals, freq = TRUE, breaks=nv/10)
  })

  output$kurt <- renderText({
    xvals <-xdist()
    sprintf("Kurtosis %0.5g", kurtosis(xvals))
  })
  output$negentropia <- renderText({
    xvals <-xdist()
    v = var(xvals)
    m = mean(xvals)
    n = rnorm(nv, mean=m, sd=sqrt(v))
    e = entropy(n) - entropy(xvals)
    sprintf("Negentropia %0.5g", e)
  })
}

shinyApp(ui = ui, server = server)
