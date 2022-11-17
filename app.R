#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
# Find out more about building applications with Shiny here:
#    http://shiny.rstudio.com/
#

library(shiny)
library(rpart.plot)
library(ggplot2)
library(GGally)

### set home directory, source files, set link to the data files
setwd("~/Desktop/Study_KCL/PhD Projects/Ensemblemethods/GitHub_current")
source("Simulating_data.R")
source("EnsembleMethods_SeparateCodesByMethod.R")
elsa_file= "~/Desktop/Study_KCL/PhD Projects/Ensemblemethods/diabetes_data_for_method.csv"
foot_file = "~/Desktop/Study_KCL/PhD Projects/Ensemblemethods/5yearfoot_ensemble.csv"
###

######## user interface #########
# Define UI for application that draws a histogram
ui <- fluidPage(
  # Application title
  titlePanel("Simulated examples for the survival ensemble methods"),
  
  ### side panel ############################################
  sidebarLayout(
    sidebarPanel(
      #sliderInput("N", "Sample size:",
      #            min = 100, max = 20000,value = 150),
      numericInput(inputId = "N",
                   label = "Sample size:",
                   value = 150),
      hr(),
      inputPanel("Simulating data",
                 selectInput(
                   "data_type",
                   label = "Select data type",
                   choices = c("Diabetes_depression",
                               "Cross-terms", 
                               "Non-linear", 
                               "Linear",
                               "ELSA_Diabetes"),
                 )
                ),
      hr(),
      numericInput(inputId = "randomseed",
                   label = "Random seed for simulations:",
                   value = 42),
      hr(),
      numericInput(inputId = "observation_time",
                   label = "Observation time",
                   value = 5.0),
      hr(),
      numericInput(inputId = "fixed_time",
                   label = "Predicted time point",
                   value = 5.0),
      hr(),
      numericInput(inputId = "percent_censored",
                   label = "Simulated drop out rate",
                   value = 0.2),
      ),
    
    ### main panel ############################################
    
    mainPanel(
     tabsetPanel(
       id = 'methods_results',
       
       tabPanel("Sample statistics", 
                tags$b("Predictors stats:"),
                DT::dataTableOutput("data_stats_table"),
                hr(),
                plotOutput("distPlot1"),
                hr(),
                plotOutput("distPlot2"),
                hr(),
                tags$b("All stats:"),
                verbatimTextOutput("data_summary"),
       ),
       
       tabPanel("CoxPH", 
                DT::dataTableOutput("cox_traintest_table"),
                DT::dataTableOutput("cox_coef"),
                hr(),
                verbatimTextOutput("cox_table"),
                ),
       
       tabPanel("SRF", 
                DT::dataTableOutput("srf_traintest_table"),
                verbatimTextOutput("srf_table")),
       
       tabPanel("Ens1", 
                DT::dataTableOutput("ens1_traintest_table"),
                verbatimTextOutput("ens1_table"),
                ),
       
       tabPanel("Ens2", 
                DT::dataTableOutput("ens2_traintest_table"),
                plotOutput("ens2_plot_tree"),
                verbatimTextOutput("ens2_display_clusters"),
                hr(),
                verbatimTextOutput("ens2_table")),
       
       tabPanel("Ens3",                 
                DT::dataTableOutput("ens3_traintest_table"),
                plotOutput("ens3_plot_tree"),
                plotOutput("ens3_plot_tree_cv"),
                DT::dataTableOutput("ens3_single_cox_coef"),
                hr(),
                verbatimTextOutput("ens3_table")),
       
       tabPanel("Performance", 
                DT::dataTableOutput("performance_cox"),
                DT::dataTableOutput("performance_table_std"),
                DT::dataTableOutput("performance_table_diff_to_cox")),
       
       tabPanel("Conclusions", 
                verbatimTextOutput("conclusions_text")),
      )
    )
  )
)

### Server ############################################

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  ### load / simulate  data set 
  datasetInput <- reactive({
    switch(input$data_type,
           "Linear" = simulatedata_linear(input$N, 
                                          observe_time = input$observation_time,
                                          randomseed = input$randomseed,
                                          percentcensored = input$percent_censored),
           "Non-linear" = simulatedata_nonlinear(input$N,
                                                 observe_time = input$observation_time,
                                                 randomseed = input$randomseed,
                                                 percentcensored = input$percent_censored),
           "Cross-terms" = simulatedata_crossterms(input$N,
                                                   observe_time = input$observation_time,
                                                   randomseed = input$randomseed,
                                                   percentcensored = input$percent_censored),
           "ELSA_Diabetes" = read.csv(elsa_file), 
           "Diabetes_depression" = read.csv(foot_file)
    )
  })
  
  ### define predictors 
  predictfactors <- reactive({
    switch(input$data_type,
           "Linear" = c("age", "bmi", "hyp", "gender"),
           "Non-linear" = c("age", "bmi", "hyp", "gender"),
           "Cross-terms" = c("age", "bmi", "hyp", "gender"),
           "ELSA_Diabetes" = c("baseline_age_", "genderdum", "bmi_0_", "cvd_0", "hyp_0",
                               "baseline_exercise",  "B_wealth","baseline_B_smokstatus",
                               "t2dm_", "pc1_", "pc2_", "pc3_"),
           "Diabetes_depression"  = c("age", "sex", "texsev", "mean_hbamva",
                                      "anydep",  "sdscas", "socclas3")
           )
  })
  
  ### Describe the data ##################
   output$distPlot1 <- renderPlot({
        x <- datasetInput()
        bins = ifelse(input$N<=500, 
                      round(input$N/20,0), 
                      round(100,0))
        times_event = x[x$event==1, "time"]
        bins <- seq(min(times_event), 
                    max(times_event), 
                    length.out = bins)
        # draw the histogram with the specified number of bins
        hist(times_event, breaks = bins, 
             col = 'darkgray', border = 'white',
             xlab = 'Time',
             main = 'Event times')
    })
  
   output$data_summary <- renderPrint({
     x <- datasetInput()
     populationstats(x, input$fixed_time, input$data_type)
   })
   
   output$data_stats_table <- DT::renderDataTable({
     x <- datasetInput()
     DT::datatable(round(describe(x),4))
   })
   
   output$distPlot2 <- renderPlot({
     x <- datasetInput()
     predict.factors <- predictfactors()
     ggpairs(x[,predict.factors])
   })
     
   output$distPlot0 <- renderPlot({
      x <- datasetInput()
      bins = ifelse(input$N<=500, 
                    round(input$N/20,0), 
                    round(100,0))
      times_no_event = x[x$event==0, "time"]
      bins <- seq(min(times_no_event), 
                  max(times_no_event), 
                  length.out = bins)
      # draw the histogram with the specified number of bins
      hist(times_no_event, breaks = bins, 
           col = 'darkgray', border = 'white',
           xlab = 'Time',
           main = 'Censored times',
           ylim = c(0, 100.0))
    })
   
    #' End (describe the data) ##################
    
    ### Apply methods  ############################################
    
    CoxPH_cv <- reactive({
      x <- datasetInput()
      predict.factors <- predictfactors()
      # Cox model
      method_cox_cv(x, 
                    predict.factors, 
                    fixed_time = input$fixed_time, 
                    cv_number = 3, 
                    seed_to_fix = input$randomseed)
      })
    
    Ens2_train_on_all <- reactive({
      x <- datasetInput()
      predict.factors <- predictfactors()
      
      method_2A_train(x, predict.factors, fixed_time = input$fixed_time, 
                          internal_cv_k =3, 
                      seed_to_fix = input$randomseed)
    })
    
    Ens3_train_on_all <- reactive({
      x <- datasetInput()
      predict.factors <- predictfactors()
      
      method_3_train(x, predict.factors, fixed_time = input$fixed_time,  
                          internal_cv_k = 3, seed_to_fix = input$randomseed)
    })

    CoxPH_train_on_all <- reactive({
      x <- datasetInput()
      predict.factors <- predictfactors()
      method_cox_train(x, predict.factors)
    })
    
    output$cox_table <- renderPrint({
      CoxPH_cv()
     })
    
    output$performance_cox <- DT::renderDataTable({
      results_on_test_sets = performance_table_reactive()[c(1,2,5,6,7), ]
      DT::datatable(results_on_test_sets)
    })
    
    results_table_diff_to_cox <- reactive({
    # computes table which has differences to Cox performance
      results_diff_to_cox = round(cbind(  
        "Cox" = CoxPH_cv()$testaverage,   
        "SRF_diff" = SRF_cv()$testaverage - CoxPH_cv()$testaverage, 
        "1_diff" = Ens1_cv()$testaverage- CoxPH_cv()$testaverage, 
        "2_diff" = Ens2_cv()$testaverage-CoxPH_cv()$testaverage,   
        "3_diff" =  Ens3_cv()$testaverage-CoxPH_cv()$testaverage),
        6)
    })
    
    SRF_cv <- reactive({
      x <- datasetInput()
      predict.factors <- predictfactors()
      
      #data descriptive tables
      output$srf_traintest_table <- DT::renderDataTable({
        x <- datasetInput()
        rrr <- round(rbind("test" = SRF_cv()$testaverage, 
                           "train" = SRF_cv()$trainaverage), 4)
        DT::datatable(rrr[, c(2,5,6,7,1)])
      })
      
    method_srf_cv(x, 
                    predict.factors, 
                    fixed_time = input$fixed_time, 
                    cv_number = 3, 
                    internal_cv_k = 3,
                    seed_to_fix = input$randomseed)
    })
    
    output$srf_table <- renderPrint({
      x <- datasetInput()
      SRF_cv()
    })
    
    output$srf_traintest_table <- DT::renderDataTable({
      x <- datasetInput()
      rrr <- round(rbind("test" = SRF_cv()$testaverage, 
                  "train" = SRF_cv()$trainaverage), 4)
      DT::datatable(rrr[, c(2,5,6,7,1)])
    })
    
    output$cox_traintest_table <- DT::renderDataTable({
      x <- datasetInput()
      rrr<- round(rbind("test" = CoxPH_cv()$testaverage, 
                          "train" = CoxPH_cv()$trainaverage),4)
      DT::datatable(rrr[, c(2,5,6,7,1)])
    })
    
    output$ens1_traintest_table <- DT::renderDataTable({
      x <- datasetInput()
      rrr<- round(rbind("test" = Ens1_cv()$testaverage, 
                                "train" = Ens1_cv()$trainaverage),4)
      DT::datatable(rrr[, c(2,5,6,7,1)])
    })
    output$ens2_traintest_table <- DT::renderDataTable({
      x <- datasetInput()
      rrr<- round(rbind("test" = Ens2_cv()$testaverage, 
                                "train" = Ens2_cv()$trainaverage),4)
      DT::datatable(rrr[, c(2,5,6,7,1)])
    })
    output$ens3_traintest_table <- DT::renderDataTable({
      x <- datasetInput()
      rrr<- round(rbind("test" = Ens3_cv()$testaverage, 
                        "train" = Ens3_cv()$trainaverage),4)
      DT::datatable(rrr[, c(2,5,6,7,1)])
    })
    
    Ens1_cv <- reactive({
      x <- datasetInput()
      predict.factors <- predictfactors()
      # Cox model
      method_1A_cv(x, 
                    predict.factors, 
                    fixed_time = input$fixed_time, 
                    cv_number = 3, 
                    internal_cv_k = 3,
                    seed_to_fix = input$randomseed)
    })
    
    output$ens1_table <- renderPrint({
      x <- datasetInput()
      Ens1_cv()
    })
    
    Ens2_cv <- reactive({
      x <- datasetInput()
      predict.factors <- predictfactors()
      # Cox model
      method_2A_cv(x, 
                   predict.factors, 
                   fixed_time = input$fixed_time, 
                   cv_number = 3, 
                   internal_cv_k = 3,
                   seed_to_fix = input$randomseed)
    })
    
    output$ens2_table <- renderPrint({
      Ens2_cv()
    })
  
    Ens3_cv <- reactive({
      x <- datasetInput()
      predict.factors <- predictfactors()
      # Cox model
      method_3_cv(x, 
                   predict.factors, 
                   fixed_time = input$fixed_time, 
                   cv_number = 3, 
                   internal_cv_k = 3,
                   seed_to_fix = input$randomseed)
    })
    
    output$ens3_table <- renderPrint({
      x <- datasetInput()
      Ens3_cv()
    })
    
    output$ens2_display_clusters <- renderPrint({
    # method 2: display Cox models in each cluster 
      m2 <- Ens2_train_on_all()
      m2$clusters #clusters correspond to the top values in cluster plot 
      m2$coxmodels
    })
    
    output$ens2_plot_tree <- renderPlot({
      # method 2: display single tree
      m2 <- Ens2_train_on_all()
      rpart.plot(m2$treemodel, nn=TRUE, roundint = FALSE, digits = 4)
    })
    
    output$ens3_plot_tree <- renderPlot({
      # method 2: display single tree
      m3 <- Ens3_train_on_all()
      rpart.plot(m3$treemodel, nn=TRUE, roundint = FALSE, digits = 4)
    })
    
    output$ens3_plot_tree_cv <- renderPlot({
      # method 2: display single tree
      m3cv <- Ens3_cv()
      cv_number = length(m3cv$tuned_cv_models)
      par(mfrow=c(cv_number+1,1))
      for (i in 1:cv_number){
        rpart.plot(m3cv$tuned_cv_models[[i]]$treemodel, 
                   nn=TRUE, roundint = FALSE, digits = 4)
      }
    })
    
    output$ens3_single_cox_coef <- DT::renderDataTable({
      #display Cox coefficients for the MODIFIED cox model 
      m3 <- Ens3_train_on_all()
      ens3_coef_df = as.data.frame(round(summary(m3$modcoxmodel)$coefficients,4))
      DT::datatable(ens3_coef_df)
    })
    
    output$cox_coef <- DT::renderDataTable({
      #display Cox coefficients for the BASELINE cox model 
      m0 <- CoxPH_train_on_all()
      coxph_coef_df = round(summary(m0)$coefficients,4)
      DT::datatable(coxph_coef_df)
    })
    
    ### final performance table (MEAN) - for all methods
    performance_table_reactive <- reactive({ 
      round(cbind(  
        "Cox" = CoxPH_cv()$testaverage,   
        "SRF" = SRF_cv()$testaverage,   
        "1A" = Ens1_cv()$testaverage,  
        "2A" = Ens2_cv()$testaverage,   
        "3" =  Ens3_cv()$testaverage),6)
    })
   
    ### final performance table (STD) - for all methods
    #-> reactive
    performance_table_std_reactive <- reactive({ 
      round(cbind(  
        "Cox_std" = apply(CoxPH_cv()$test,  2, sd),
        "SRF_std" = apply(SRF_cv()$test,  2, sd),  
        "Ens1_std" = apply(Ens1_cv()$test,  2, sd),  
        "Ens2_std" =apply(Ens2_cv()$test,  2, sd),   
        "Ens3_std" = apply(Ens3_cv()$test,  2, sd)),6)
    })
    #-> output
    output$performance_table_std <- DT::renderDataTable({
      # outcomes and conclusions on linear/non-linear effects 
      DT::datatable(
        performance_table_std_reactive()[c(1,2,5,6,7),])
    })
    
    output$performance_table_diff_to_cox <- DT::renderDataTable({
      # outcomes and conclusions on linear/non-linear effects 
      DT::datatable(
        results_table_diff_to_cox()[c(1,2,5,6,7),])
    })
    
    output$conclusions_text<- renderPrint({
        # outcomes and conclusions on linear/non-linear effects 
        #  results_table_diff_to_cox()
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
