source("./mapRender.R")
source("./renderPFEtable.R")
source("./renderIFEtable.R")

get.state <- function(input) {
  if(input$countrySelector == "Japan") {
    return(input$prefectureSelector)
  }
  return(input$departmentSelector)
}

get.pca <- function(input) {
  if(input$pcaSelector == "EVR") {
    return(input$pcaSelector)
  } 
  return("K")
}


SERVER <- function(input, output) {
  output$map <- renderLeaflet({
    render.map(input$countrySelector, get.state(input))
  })
  output$PFETable <- DT::renderDataTable({
    
    sketch.upper <- htmltools::withTags(table(
      class = 'display',
      thead(
        tr(
          th(colspan = 1, '', style = 'text-align:center;'),
          th(colspan = 2, 'UFTS', style = 'text-align:center;'),
          th(colspan = 2, 'MFTS', style = 'text-align:center;'),
          th(colspan = 2, 'MLFTS', style = 'text-align:center;'),
          th(colspan = 2, 'FANOVA', style = 'text-align:center;'),
          th(colspan = 2, 'HDFPCA', style = 'text-align:center;')
        ),
        tr(
          lapply(c('h', rep(c('KLD', 'JSD'), 5)), function(x) th(x, style = 'text-align:center;'))
        )
      )
    ))
    DT::datatable(gen_pfe_table(input$countrySelector, get.state(input), 
                                input$genderSelector, get.pca(input)),
                  options = list(pageLength = 8,
                                 ordering = F), 
                  rownames = F, container = sketch.upper) %>% 
      DT::formatStyle('h',
                      fontWeight = DT::styleEqual(c("Mean","Median"), "bold"))
  })
  output$IFETable <- DT::renderDataTable({
    sketch.upper <- htmltools::withTags(table(
      class = 'display',
      thead(
        tr(
          th(colspan = 1, '', style = 'text-align:center;'),
          th(colspan = 3, 'UFTS', style = 'text-align:center;'),
          th(colspan = 3, 'MFTS', style = 'text-align:center;'),
          th(colspan = 3, 'MLFTS', style = 'text-align:center;'),
          th(colspan = 3, 'FANOVA', style = 'text-align:center;'),
          th(colspan = 3, 'HDFPCA', style = 'text-align:center;')
        ),
        tr(
          lapply(c('h', rep(c('ECP', 'CPD', 'IS'), 5)), function(x) th(x, style = 'text-align:center;'))
        )
      )
    ))
    DT::datatable(gen_ife_table(input$countrySelector, get.state(input), 
                                input$genderSelector, input$coverageSelector, get.pca(input)),
                  options = list(pageLength =8,
                                 ordering = F),
                  rownames = F, container = sketch.upper) %>% 
      DT::formatStyle('h',
                      fontWeight = DT::styleEqual(c("Mean","Median"), "bold"))
    
  })
  
}





