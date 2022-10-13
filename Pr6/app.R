library(shiny)
library(DT)

# Define UI for application working with power
ui <- fluidPage(

    # Application title
    titlePanel("Computation of required survey effort, followed by power to detect change"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("R",
                        "Cumulative change over study:",
                        min = 0.1,
                        max = 0.9,
                        value = 0.5),
            sliderInput("nyears",
                        "Number of annual surveys:",
                        min = 2,
                        max = 15,
                        value = 5),
            sliderInput("power",
                        "Power to detect change:",
                        min = 0.5,
                        max = 0.9,
                        value = .8),
            sliderInput("encrate",
                        "Encounter rate (detects per km):",
                        min = 0.1,
                        max = 2,
                        value = .25),
            sliderInput("b",
                        "Dispersion factor for design:",
                        min = 1,
                        max = 20,
                        value = 3),
            sliderInput("alpha",
                        "Type I error rate:",
                        min = 0.01,
                        max = .5,
                        value = .05),
        ),

        # Show results in five tabbed panels
        mainPanel(
          tabsetPanel(id="tabsets", 
            tabPanel("Description",
                     tags$h4("Description"),
                     tags$p("At the most basic, given an encounter rate from a pilot study, this application can help design a line transect survey.  Using Eqn. 2.4.2.2. of Buckland et al. (2015), the application can estimate the amount of line transect effort needed to estimate precision with a desired level of precision (measured by coefficient of variation)."),
                     tags$p("Subsequent to this calculation, computation of CV needed to achieve desired power based on Eqn. 10 of Gerrodette (1987). That formula assumes"),
                     tags$ul(
                       tags$li("CV is a constant over abundance"),
                       tags$li("Rate of change is linear over time"),
                     ),
                     tags$hr(),
                     tags$p("If interest lies only in designing a survey with a desired level of precision, the only slider that needs to be used is the Encounter Rate slider.  Consequently, the only result tab of interest is that labelled 'Eqn. 2.4'.  To conduct the prospective power calculation, the other sliders must be used and the other output tabs consulted."),
                     tags$h4("References"),
                     tags$ul(
                       tags$li("Buckland, S. T., Rexstad, E. A., Marques, T. A., & Oedekoven, C. S. (2015). Distance sampling: Methods and applications. Springer.
"),
                       tags$li("Gerrodette, T. (1987). A power analysis for detecting trends. Ecology, 68(5), 1364–1372. https://doi.org/10.2307/1939220
"),
                     )
                     ),
            tabPanel("Eqn. 2.4", value="eqn24Plot",
                     plotOutput("eqn24Plot")),
            tabPanel("CV graph", value="cvPlot",
                     plotOutput("cvPlot")),
            tabPanel("Effort graph", value="effortPlot",
                     plotOutput("effortPlot")),
            tabPanel("Numerical result for power calculation",
                     verbatimTextOutput("answer"),
                     textOutput("cv"),
                     textOutput("effort"))
          )
        )
    )
)

# Define server logic required to calculate power
server <- function(input, output) {
  
  lin.constn.cv <- function(R=0.1, n=10, alpha=0.05, power=0.8) {
    r <- R/(n-1)
    t1 <- r^2*n*(n-1)*(n+1)
    t2 <- 1+r*(n-1)*(1+r/6*(2*n-1))
    beta <- 1-power
    t3 <- 12 * (qnorm(alpha/2) + qnorm(beta))^2
    cv <- sqrt(t1/(t3*t2))
    return(cv)
  }
  eqn2.4 <- function(b, encrate, cv) {
    effort <- b/cv^2 * 1/encrate
    return(effort)
  }

  returnData <- 
    output$answer <- renderText({
      paste("Given:\n  Cumulative change=", input$R, 
            "\n  Years of surveys=", input$nyears,
            "\n  Power=", input$power,
            "\n  Encounter rate=", input$encrate,
            "\n  b=", input$b,
            "\n  alpha=", input$alpha)
                               })
      output$cv <- renderText({ paste("Required CV=", 
                                      round(lin.constn.cv(R=input$R, n=input$nyears,
                                                 alpha=input$alpha, power=input$power), 3))
                                })
      output$effort <- renderText({paste("Required effort=",
                                         round(eqn2.4(b=input$b, encrate=input$encrate,
                                                cv=lin.constn.cv(R=input$R, n=input$nyears, 
                                                                 alpha=input$alpha, power=input$power))))
                                })
      output$eqn24Plot <- renderPlot({
        cvseq <- seq(0.1, 0.5, by=0.01)
        resulteffort <- eqn2.4(b = input$b, encrate = input$encrate, cv=cvseq)
        plot(resulteffort, cvseq, main=paste("Effort to achieve desired CV\ngiven encounter rate=", input$encrate),
             ylab="Target CV", xlab="Effort (km)", type="l", lwd=2, 
             sub="Independent of popn change, number of years and power", font.sub=3)
      output$cvPlot <- renderPlot({
        rseq <- seq(0.1, 0.9, by=0.01)
        thecv <- lin.constn.cv(R=rseq, n=input$nyears, alpha=input$alpha, power=input$power)        
        plot(rseq, thecv, type="l", lwd=4,
             main=paste("CV with power=", input$power, "Number surveys=", input$nyears),
             xlab="Overall change in abundance",
             ylab="Needed coefficient of variation",
             ylim=c(0, 0.3),
             sub="Encounter rate and dispersion factor play no role in this calculation"
        )
        points(input$R, lin.constn.cv(R=input$R, n=input$nyears, alpha=input$alpha, power=input$power), 
               pch=19, col="red", cex=2)
        segments(0.0, lin.constn.cv(R=input$R, n=input$nyears, alpha=input$alpha, power=input$power),
                 input$R, lin.constn.cv(R=input$R, n=input$nyears, alpha=input$alpha, power=input$power),
                 lty=3)

      output$effortPlot <- renderPlot({
        rseq <- seq(0.1, 0.9, by=0.01)
        effort <- eqn2.4(b=input$b, encrate=input$encrate,
                         cv=lin.constn.cv(R=rseq, n=input$nyears, alpha=input$alpha, power=input$power))
        plot(rseq, effort, type="l", lwd=4,
             main=paste("Effort required with power=", input$power, 
                        "Number surveys=", input$nyears,
                        "Enc rate=", input$encrate,
                        "\n b=", input$b,
                        "alpha=", input$alpha),
             xlab="Overall change in abundance",
             ylab="Target effort per survey (km)",
             ylim=c(0, 20000)
        )
        points(input$R, 
               eqn2.4(b=input$b, encrate=input$encrate,
                      cv=lin.constn.cv(R=input$R, n=input$nyears, power=input$power)), 
               pch=19, col="red", cex=2)
        segments(0.0, eqn2.4(b=input$b, encrate=input$encrate,
                             cv=lin.constn.cv(R=input$R, n=input$nyears, power=input$power)), 
                 input$R, eqn2.4(b=input$b, encrate=input$encrate,
                                 cv=lin.constn.cv(R=input$R, n=input$nyears, power=input$power)), 
                 lty=3)
      })
    })
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
