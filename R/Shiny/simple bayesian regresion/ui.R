### Introduction


file_name<- fileInput('file_name', 'Choose File',
                      accept=c('text/csv', 
                               'text/comma-separated-values,text/plain', 
                               '.csv'))

header_file<- checkboxInput('header_file', 'Header', TRUE)


sep_file<- radioButtons('sep_file', 'Separator',
                        c(Comma=',',
                          Semicolon=';',
                          Tab='\t'),
                        selected=',')
prev_data=checkboxInput('prev_data','preview of loaded data?',value = FALSE)
upload_file=fluidRow(column(4,file_name),column(2,header_file),column(2,sep_file),column(2,prev_data))



ui=fluidPage(
  includeCSS('mbie-styles.css'),
  includeCSS('tdstyles.css'),
  div(class='front-banner',
    div(class='imgcon'),
    div(class='hcon',h1('Mateo Graciano'),h6('Simple app fo bayesian modeling'))
  ),
  tabsetPanel(
    tabPanel('Introduction',
             h2('Here goes your own intro')
    ),
  tabPanel('Models',
    h4("Select file with your data"),
    upload_file,
    uiOutput('conditional')
  )
  
  )
)
  
