library(shiny)
library(tidyverse)
library(leaflet)

Ui <- fluidPage(
  fluidRow(
    column(4,
           wellPanel(
             div(id = "hidden",
                 selectInput(
                   "countrySelector",
                   "Country",
                   c("Japan")
                 ),
             ),
             conditionalPanel(
               "input.countrySelector == 'Japan'",
               selectInput("prefectureSelector",
                           "Prefecture",
                           sapply(readRDS("./names/names_prefectures.rds"), 
                                  function(x) x)
               )
             ),
             conditionalPanel(
               "input.countrySelector == 'France'",
               selectInput("departmentSelector",
                           "Department",
                           sapply(readRDS("./names/names_departments.rds"), 
                                  function(x) x)
               )
             ),
             selectInput("genderSelector",
                         "Gender",
                         c("Male", "Female")
             ),
             selectInput(
               "pcaSelector",
               "Selection method for the FPCA",
               c("EVR", "K = 6")
             ),
             selectInput(
               "coverageSelector",
               "Nominal Coverage",
               c("80%", "95%")
             )
           ),
           div(
             style = "display: flex; justify-content: center; align-items: center",
             h3(
               "Modeling and forecasting subnational age distribution of death counts",
               style = "font-weight: bold; text-align: center;"
             )),
           div(
             style = "display: flex; justify-content: center; align-items: center; padding-top:5%",
             h4(
               "Authors",
               br(),
               br(),
               "Han Lin Shang",
               br(),
               "Department of Actuarial Studies and Business Analytics",
               "Macquarie University",
               br(),
               br(),
               "Cristian F. Jiménez-Varón",
               br(),
               "University of York",
               br(),
               "Department of Mathematics",
             )
           )
    ),
    column(4,
           leafletOutput("map")
    ),
    column(4,
           div(
             style = "display: flex; justify-content: center; align-items: center",
             h4("Point forecast evaluation", style = "font-weight: bold;")
           ),
           div(DT::dataTableOutput("PFETable")),
           DT::dataTableOutput("summaryMetricsTable"),
           conditionalPanel(
             "input.coverageSelector == '80%'",
             div(
               style = "display: flex; justify-content: center; align-items: center",
               h4("Interval forecast evaluation 80% nominal coverage (sd approach)",
                  style = "font-weight: bold;")
             )
           ),
           conditionalPanel(
             "input.coverageSelector == '95%'",
             div(
               style = "display: flex; justify-content: center; align-items: center",
               h4("Interval forecast evaluation 95% nominal coverage (sd approach)",
                  style = "font-weight: bold;")
             )
           ),
           div(DT::dataTableOutput("IFETable")),
           div(
             style = "display: flex; justify-content: center; align-items: center",
             h6(style = "font-weight: bold;")
           )
    )
  ),
  tags$style(type = "text/css",
             ".dataTables_length, .dataTables_filter {display:none;}
             h4 {text-align:center;}
             #map {height: calc(100vh - 30px) !important;}
             td {text-align:center !important;}
             #hidden{display: none")
)