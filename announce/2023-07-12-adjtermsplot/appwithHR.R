library(shiny)

# Define UI for application that draws detection functions
ui <- fluidPage(

    # Application title
    titlePanel("Key functions with cosine adjustments"),

    sidebarLayout(
        sidebarPanel(
          radioButtons("key",
                       "Key function:",
                       choices = c("HN", "HR"),
                       selected = "HN",
                       inline = TRUE),
          sliderInput("truncation",
                      "truncation distance:",
                      min = 100,
                      max = 200, step=10,
                      value = 200),
          sliderInput("sigma",
                        "sigma:",
                        min = 1,
                        max = 200, step=5,
                        value = 30),
          conditionalPanel(
            condition = "input.key == 'HR'",
              sliderInput("beta",
                          "beta:",
                          min = .1,
                          max = 10, step=.5,
                          value = 3)
          ),
          numericInput("nadj",
                       "Number of adjustment terms",
                       value = 0,
                       min = 0,
                       max = 2,
                       step = 1,
                       width = '25%'),
          # radioButtons("nadj",
          #              "Number of adjustment terms",
          #              choices = c("0", "1", "2"),
          #              selected = "0",
          #              inline = TRUE),
          conditionalPanel(
            condition = "input.nadj > 0",
              sliderInput("adj1",
                          "Cosine adjustment coeff 1:",
                          min = -0.25,
                          max = 0.25, step=0.05,
                          value = 0)
          ),
          conditionalPanel(
            condition = "input.nadj > 1",
              sliderInput("adj2",
                        "Cosine adjustment coeff 2:",
                        min = -0.25,
                        max = 0.25, step=0.05,
                        value = 0)
        )
        ),
        mainPanel(
           plotOutput("distPlot"),
           htmlOutput("description")
        )
    )
)

# Define server logic required to draw detection functions
server <- function(input, output) {
  gz<-function(z,
               sigma, adj1, adj2, beta,
               key="HN",adj_num=0,adjt="cos",w=max(z)){
    #this is a generic detection function that returns the probability of detecting an animal given the arguments:
    #    z               generic distance (radial or perpendicular) - can be scalar or vector
    #    key             the detection function key, default is Half Normal (HR and UNI might be alternatives)
    #    pars            the detection function parameter values (HN:s,adj1,adj2,...)(HR:s,b,adj1,adj2,...)(UNI:adj1,adj2,...)
    #    adjn            if 0, no adjustment terms are used, if > 0, number of adjustments to use
    #    adjt            type of adjustments to use, defaults to cosine series
    #    w               truncation distance, only fundamental for adjustments, by default the max of the distances
    #
    #RETURNS: a probability
    #
    #the standard HN and HR
    if(adj_num == 1) {
      if(key=="HR") {
        pars <- c(sigma, beta, adj1)
      } else {  # HN
        pars <- c(sigma, adj1)
      }
    }
    if(adj_num == 2) {
      if(key=="HR") {
        pars <- c(sigma, beta, adj1, adj2)
      } else {  # HN
        pars <- c(sigma, adj1, adj2)
      }
    }
    
    if(adj_num==0 & key=="HN") {
      pars <- sigma
      p<-exp(-z^2/(2*pars[1]^2))
    }
    if(adj_num==0 & key=="HR") {
      pars <- c(sigma, beta)
      scale.dist <- z/pars[1]
      inside <- -(scale.dist)^(-pars[2])
      p <- 1 - exp(inside)
    }
    #including cosines in the HN key
    if(adj_num>0 & key=="HN") {
      hn.bit<-exp(-z^2/(2*pars[1]^2))
      cos.bits<-matrix(nrow=length(z),ncol=adj_num)#storage for the cosine terms
      for(j in 1:adj_num) {
        cos.bits[,j]<-pars[j+1]*cos((j+1)*pi*z/w)
      }
      p<-(hn.bit*(1+apply(cos.bits,1,sum)))/(1+sum(pars[(length(pars)-adj_num+1):length(pars)]))
    }
    if(adj_num>0 & key=="HR") {
      scale.dist <- z/pars[1]
      inside <- -(scale.dist)^(-pars[2])
      hr.bit <- 1 - exp(inside)
      cos.bits<-matrix(nrow=length(z),ncol=adj_num)#storage for the cosine terms
      for(j in 1:adj_num) {
        cos.bits[,j]<-pars[j+2]*cos((j+1)*pi*z/w)
      }
      p<-(hr.bit*(1+apply(cos.bits,1,sum)))/(1+sum(pars[(length(pars)-adj_num+1):length(pars)]))
    }  
    #including cosines in the Uniform key
    if(adj_num>0 & key=="UNI") {
      cos.bits<-matrix(nrow=length(z),ncol=adj_num)#storage for the cosine terms
      for(j in 1:adj_num) {
        cos.bits[,j]<-pars[j+1]*cos(j*pi*z/w)
      }
      p<-(1*(1+apply(cos.bits,1,sum)))/(1+sum(pars[(length(pars)-adj_num+1):length(pars)]))
    }
    return(p)
  }
    output$distPlot <- renderPlot({
      distances <- seq(0, input$truncation, by=0.1)
      plot(distances, 
           gz(distances, key=input$key, 
              sigma = input$sigma, 
              adj1 = input$adj1,
              adj2 = input$adj2, 
              beta = input$beta, adjt="cos", adj_num=input$nadj),
           type="l", lwd=4, ylim=c(0,1),
           xlab="Distance", ylab="Detection probability")
    })
   output$description <- renderUI({
     str1 <- "Visual illustration of the effect of adjustments upon half normal and hazard rate detection function shape."
     str2 <- "This code has no error checking; for example, it makes no sense to have a second adjustment term without a first term"
     str4 <- "It is for illustrative purposes."
     HTML(paste(str1, str2, str4, sep='<br/>'))
   })
}

# Run the application 
shinyApp(ui = ui, server = server)
