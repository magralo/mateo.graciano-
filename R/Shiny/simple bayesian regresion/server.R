library(tidyverse)
library(MCMCpack)
library(rhandsontable)
library(DT)
source("Models/Normal.R")
source("Models/Probit.R")
modelos=c('Normal Model','Probit Model')
sel_model=selectInput("sel_model", "Select Model type",
                      modelos)

#### UI For each model
it1<- sliderInput("it", 
                  "MCMC Iterations:", 
                  value = 10000,
                  min = 10000, 
                  max = 100000,
                  step = 10000)
it2<- sliderInput("burnin", 
                  "Burn-in Sample:", 
                  value = 1000,
                  min = 1000, 
                  max = 10000,
                  step = 1000)

it3<- selectInput("keep", "Thinning parameter:", 
                  choices = c("1", "10", "20", "50", "100"), selected = "1")
MCMC_general=fluidRow(column(4,it1),column(4,it2),column(4,it3))

linear_normal=fluidPage(
  fluidRow(
    
    column(3,h6('Prior mean for your model'),rHandsontableOutput("priorMean1")),
    column(4,h6('Prior covariance for your model'),rHandsontableOutput("priorCov1"))      
           ),
  br(),
   fluidRow(
     column(3,numericInput('normal_scale','Prior scale parameter',value = 0.001,step = 0.001)),
     column(3,numericInput('normal_shape','Prior shape parameter',value = 0.001,step = 0.001))
   ),
  br(),
  h3('Specify MCMC parameters'),
  MCMC_general,
  actionButton("doModel", "Go!"),
  plotOutput("plotBeta")
)

uni_probit=fluidPage(
  fluidRow(
    
    column(3,h6('Prior mean for your model'),rHandsontableOutput("priorMean1")),
    column(4,h6('Prior covariance for your model'),rHandsontableOutput("priorCov1"))      
  ),
  br(),
  h3('Specify MCMC parameters'),
  MCMC_general,
  actionButton("doModel", "Go!"),
  plotOutput("plotBeta")
  
  
  
)


server <- function(input, output){
  ####Lectura de datos
  dataInput <- reactive({
    inFile1 <- input$file_name
    if (is.null(inFile1))
      return(NULL)
    read.csv(inFile1$datapath, header=input$header_file, sep=input$sep_file)
  })
  
  output$conditional=renderUI({
    if (!is.null(dataInput())){
      if(input$'prev_data'){
        fluidPage(h4('Preview of loaded data'),
                  dataTableOutput('preview'),
                  h3('Build your model'),
                  sel_model,
                  fluidRow(
                    column(3, selectInput("sel_y", "Select Model objective (y)",
                                          colnames(dataInput()))),
                    column(4,selectInput("sel_X", "Select Model objective (X)",
                                         colnames(dataInput()),multiple = TRUE))
                  ),
                  uiOutput('model_UI')
                  
        )
      }else{
        fluidPage(
              h3('Build your model'),
              sel_model,
              fluidRow(
              column(3, selectInput("sel_y", "Select Model objective (y)",
                          colnames(dataInput()))),
              column(4,selectInput("sel_X", "Select Model regressors (X)",
                          colnames(dataInput()),multiple = TRUE))
              ),
              uiOutput('model_UI')
              
              
        )
      }
      
    } else{
      fluidPage(h3('No loaded data'))
    }
  })
  
  output$model_UI=renderUI({
    switch(input$sel_model,
           'Normal Model' =linear_normal ,
           'Probit Model' = uni_probit
   )
  })
  
  output$preview=renderDataTable({
    datatable(head(dataInput()),rownames= FALSE)
  })
  
  output$priorMean1=renderRHandsontable({
    if(is.null(input$priorMean1) ){
      if(length(input$sel_X)>0){
        names=c('cte',input$sel_X)
        DF=data.frame(mean=rep(0,(length(names))))
       
        rownames(DF)=names
      }else{
        DF=data.frame(mean=0)
        rownames(DF)='cte'
      }
      
    }else{
      DF=hot_to_r(input$priorMean1)
      rn=rownames(DF)[-1]
      if(!identical(rn,input$sel_X)){
        names=c('cte',input$sel_X)
        DF=data.frame(mean=rep(0,(length(names))))
        
        rownames(DF)=names
      }
    }
    
    rhandsontable(DF)%>% 
      hot_col("mean",format="0.01")
    
  })
  
  
  output$priorCov1=renderRHandsontable({
    if(is.null(input$priorCov1) ){
      if(length(input$sel_X)>0){
        names=c('cte',input$sel_X)
        DF=as.data.frame(diag(length(names)))*0.001
        rownames(DF)=names
        colnames(DF)=names
      }else{
        DF=data.frame(mean=0.001)
        rownames(DF)='cte'
      }
      
    }else{
      DF=hot_to_r(input$priorCov1)
      
      rn=rownames(DF)[-1]
      if(!identical(rn,input$sel_X)){
        names=c('cte',input$sel_X)
        DF=as.data.frame(diag(length(names)))*0.001
        rownames(DF)=names
        colnames(DF)=names
      }
    }
    
    rhandsontable(DF)%>% 
      hot_col(colnames(DF),format="0.001")
    
  })
  rv <- reactiveValues(
    Posteriors=NULL
  )
  
  observeEvent(input$doModel, {
    print('hago algo')
    Bmean<- hot_to_r(input$priorMean1)[,1]
    
    Bvar<- solve(as.matrix(hot_to_r(input$priorCov1)))
    
    if (input$sel_model=='Normal Model'){
      a<- input$normal_scale
      b<- input$normal_shape
    }
    MCMC<- list(R=input$it,keep=as.numeric(input$keep),burnin=input$burnin)
    Data=list()
    data=dataInput()
    X=cbind(data.frame(cte=rep(0,nrow(data))),dplyr::select(data,input$sel_X))
    Data$X=as.matrix(X)
    Data$y=dplyr::select(data,input$sel_y)[,1]
    args <- switch(input$sel_model,
                   'Normal Model' = list(Data,list(betabar=Bmean,A=Bvar,a=a,b=b),MCMC),
                   'Probit Model' = list(Data,list(betabar=Bmean,A=Bvar),MCMC))
     outputModel=switch(input$sel_model,
                        'Normal Model' = do.call(Normal, args),
                        'Probit Model' = do.call(Probit, args)
     )
    
    rv$Posteriors=outputModel
    print(head(outputModel$betadraw))
  })
  
  output$plotBeta <- renderPlot({
    if(!is.null(rv$Posteriors)){
    data=as.data.frame(rv$Posteriors$betadraw)
    colnames(data)=c('cte',input$sel_X)
    data%>%
      gather(name,value)%>%
      ggplot(aes(value,fill=name))+geom_histogram()+ facet_grid(~ name)
    }
  })
  
  
}
