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
### set home directory, source files, set link to the data files

#laptop:
#setwd("~/Desktop/Study_KCL/PhD Projects/Ensemblemethods/GitHub_current")
#elsa_file= "~/Desktop/Study_KCL/PhD Projects/Ensemblemethods/diabetes_data_for_method.csv"
#foot_file = "~/Desktop/Study_KCL/PhD Projects/Ensemblemethods/5yearfoot_ensemble.csv"

#window pc
setwd("C:/Users/dinab/Desktop/PhD Projects/Ensemble methods/R Code")
elsa_file= "C:/Users/dinab/Desktop/PhD Projects/Ensemble methods/R Code/diabetes_data_for_method.csv"
foot_file = "C:/Users/dinab/Desktop/PhD Projects/Ensemble methods/R Code/5yearfoot_ensemble.csv"
hnscc_file = "C:/Users/dinab/Desktop/PhD Projects/Ensemble methods/R Code/hnscc_merged.csv"


source("Simulating_data.R")
source("EnsembleMethods_SeparateCodesByMethod.R")

###

st = ' "baseline_age_", "genderdum", "bmi_0_", "cvd_0", "hyp_0", "baseline_exercise",  "B_wealth", 
                             "baseline_B_smokstatus",  "t2dm_", "pc1_", "pc2_", "pc3_" '
Clean_String <- function(string){
  temp = string
  # Remove everything that is not a number or letter (may want to keep more 
  # stuff in your actual analyses). 
  temp <- stringr::str_replace_all(temp,",", " ")
  temp <- stringr::str_replace_all(temp,"'", " ")
  temp <- stringr::str_replace_all(temp,'["]', '')
  # Shrink down to just one white space
  temp <- stringr::str_replace_all(temp,"[\\s]+", " ")
  # Split it
  temp <- stringr::str_split(temp, " ")[[1]]
  # Get rid of trailing "" if necessary
  indexes <- which(temp == "")
  if(length(indexes) > 0){
    temp <- temp[-indexes]
  } 
  return(temp)
}

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
      selectInput(
        "data_type",
        label = "Select data type",
        choices = c("Simulated: Cross-terms", 
                    "Simulated: Non-linear", 
                    "Simulated: Linear",
                    "Diabetes_depression",
                    "ELSA_Diabetes",
                    "Hnscc",
                    "Custom"),
      ),
      numericInput(inputId = "fixed_time",
                   label = "Time point for event prediction:",
                   value = 5.0),
      
      numericInput(inputId = "randomseed_validation",
                   label = "Random seed for calibration and validation",
                   value = 42),
      
      numericInput(inputId = "k_outer",
                   label = "K_Outer loop CV (for validation)",
                   value = 3),
      
      numericInput(inputId = "k_inner",
                   label = "K_ inner CV folds (model tuning)",
                   value = 3),
      
      inputPanel("",
                 numericInput(inputId = "randomseed",
                              label = "Simulated data: random seed (generation):",
                              value = 42),
                 
                 numericInput(inputId = "N",
                              label = "Sample size:",
                              value = 150),
                 
                 numericInput(inputId = "observation_time",
                              label = "Observation time",
                              value = 5.0),
                 
                 numericInput(inputId = "percent_censored",
                              label = "Simulated drop out rate",
                              value = 0.2),
      ),
      
      inputPanel("",
                 textInput(
                   "custom_file",
                   label = "Custom data: path to data file",
                   value = "~/Desktop/Study_KCL/PhD Projects/Ensemblemethods/diabetes_data_for_method.csv",
                   placeholder = "~/Desktop/Study_KCL/PhD Projects/Ensemblemethods/diabetes_data_for_method.csv"
                 ),
                 
                 textInput(inputId = "custom_predictors",
                           label = "Predictors to use in the model",
                           value = ' "baseline_age_", "genderdum", "bmi_0_", "cvd_0", "hyp_0", "baseline_exercise",  "B_wealth", 
                             "baseline_B_smokstatus",  "t2dm_", "pc1_", "pc2_", "pc3_" ',
                           placeholder = "age, bmi, gender, hyp"),
                 
                 textInput(inputId = "custom_time",
                           label = "Time variable name",
                           value = "time"),
                 
                 textInput(inputId = "custom_event",
                           label = "Event indicator variable name",
                           value = "event"),
      ),
      
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
                 tags$b("Population stats:"),
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
           "Simulated: Linear" = simulatedata_linear(input$N, 
                                                     observe_time = input$observation_time,
                                                     randomseed = input$randomseed,
                                                     percentcensored = input$percent_censored),
           "Simulated: Non-linear" = simulatedata_nonlinear(input$N,
                                                            observe_time = input$observation_time,
                                                            randomseed = input$randomseed,
                                                            percentcensored = input$percent_censored),
           "Simulated: Cross-terms" = simulatedata_crossterms(input$N,
                                                              observe_time = input$observation_time,
                                                              randomseed = input$randomseed,
                                                              percentcensored = input$percent_censored),
           "ELSA_Diabetes" = read.csv(elsa_file), 
           "Diabetes_depression" = read.csv(foot_file),
           "Custom" = read.csv(input$custom_file),
           "Hnscc" = read.csv(hnscc_file)
    )
  })
  
  ### define predictors 
  predictfactors <- reactive({
    switch(input$data_type,
           "Simulated: Linear" = c("age", "bmi", "hyp", "gender"),
           "Simulated: Non-linear" = c("age", "bmi", "hyp", "gender"),
           "Simulated: Cross-terms" = c("age", "bmi", "hyp", "gender"),
           "ELSA_Diabetes" = c("baseline_age_", "genderdum", "bmi_0_", "cvd_0", "hyp_0",
                               "baseline_exercise",  "B_wealth","baseline_B_smokstatus",
                               "t2dm_", "pc1_", "pc2_", "pc3_"),
           "Diabetes_depression"  = c("age", "sex", "texsev", "mean_hbamva",
                                      "anydep",  "sdscas", "socclas3"),
           "Hnscc" = c("original_firstorder_10Percentile" ,"original_firstorder_90Percentile" ,"original_firstorder_Energy" ,"original_firstorder_Entropy" ,"original_firstorder_InterquartileRange" ,"original_firstorder_Kurtosis" ,"original_firstorder_Maximum" ,"original_firstorder_Mean" ,"original_firstorder_MeanAbsoluteDeviation" ,"original_firstorder_Median" ,"original_firstorder_Minimum" ,"original_firstorder_Range" ,"original_firstorder_RobustMeanAbsoluteDeviation" ,"original_firstorder_RootMeanSquared" ,"original_firstorder_Skewness" ,"original_firstorder_TotalEnergy" ,"original_firstorder_Uniformity" ,"original_firstorder_Variance" ,"original_glcm_Autocorrelation" ,"original_glcm_ClusterProminence" ,"original_glcm_ClusterShade" ,"original_glcm_ClusterTendency" ,"original_glcm_Contrast" ,"original_glcm_Correlation" ,"original_glcm_DifferenceAverage" ,"original_glcm_DifferenceEntropy" ,"original_glcm_DifferenceVariance" ,"original_glcm_Id" ,"original_glcm_Idm" ,"original_glcm_Idmn" ,"original_glcm_Idn" ,"original_glcm_Imc1" ,"original_glcm_Imc2" ,"original_glcm_InverseVariance" ,"original_glcm_JointAverage" ,"original_glcm_JointEnergy" ,"original_glcm_JointEntropy" ,"original_glcm_MCC" ,"original_glcm_MaximumProbability" ,"original_glcm_SumAverage" ,"original_glcm_SumEntropy" ,"original_glcm_SumSquares" ,"original_gldm_DependenceEntropy" ,"original_gldm_DependenceNonUniformity" ,"original_gldm_DependenceNonUniformityNormalized" ,"original_gldm_DependenceVariance" ,"original_gldm_GrayLevelNonUniformity" ,"original_gldm_GrayLevelVariance" ,"original_gldm_HighGrayLevelEmphasis" ,"original_gldm_LargeDependenceEmphasis" ,"original_gldm_LargeDependenceHighGrayLevelEmphasis" ,"original_gldm_LargeDependenceLowGrayLevelEmphasis" ,"original_gldm_LowGrayLevelEmphasis" ,"original_gldm_SmallDependenceEmphasis" ,"original_gldm_SmallDependenceHighGrayLevelEmphasis" ,"original_gldm_SmallDependenceLowGrayLevelEmphasis" ,"original_glrlm_GrayLevelNonUniformity" ,"original_glrlm_GrayLevelNonUniformityNormalized" ,"original_glrlm_GrayLevelVariance" ,"original_glrlm_HighGrayLevelRunEmphasis" ,"original_glrlm_LongRunEmphasis" ,"original_glrlm_LongRunHighGrayLevelEmphasis" ,"original_glrlm_LongRunLowGrayLevelEmphasis" ,"original_glrlm_LowGrayLevelRunEmphasis" ,"original_glrlm_RunEntropy" ,"original_glrlm_RunLengthNonUniformity" ,"original_glrlm_RunLengthNonUniformityNormalized" ,"original_glrlm_RunPercentage" ,"original_glrlm_RunVariance" ,"original_glrlm_ShortRunEmphasis" ,"original_glrlm_ShortRunHighGrayLevelEmphasis" ,"original_glrlm_ShortRunLowGrayLevelEmphasis" ,"original_glszm_GrayLevelNonUniformity" ,"original_glszm_GrayLevelNonUniformityNormalized" ,"original_glszm_GrayLevelVariance" ,"original_glszm_HighGrayLevelZoneEmphasis" ,"original_glszm_LargeAreaEmphasis" ,"original_glszm_LargeAreaHighGrayLevelEmphasis" ,"original_glszm_LargeAreaLowGrayLevelEmphasis" ,"original_glszm_LowGrayLevelZoneEmphasis" ,"original_glszm_SizeZoneNonUniformity" ,"original_glszm_SizeZoneNonUniformityNormalized" ,"original_glszm_SmallAreaEmphasis" ,"original_glszm_SmallAreaHighGrayLevelEmphasis" ,"original_glszm_SmallAreaLowGrayLevelEmphasis" ,"original_glszm_ZoneEntropy" ,"original_glszm_ZonePercentage" ,"original_glszm_ZoneVariance" ,"original_ngtdm_Busyness" ,"original_ngtdm_Coarseness" ,"original_ngtdm_Complexity" ,"original_ngtdm_Contrast" ,"original_ngtdm_Strength" ,"original_shape_Elongation" ,"original_shape_Flatness" ,"original_shape_LeastAxisLength" ,"original_shape_MajorAxisLength" ,"original_shape_Maximum2DDiameterColumn" ,"original_shape_Maximum2DDiameterRow" ,"original_shape_Maximum2DDiameterSlice" ,"original_shape_Maximum3DDiameter" ,"original_shape_MeshVolume" ,"original_shape_MinorAxisLength" ,"original_shape_Sphericity" 
                       ,"original_shape_SurfaceArea" ,"original_shape_SurfaceVolumeRatio" ,"original_shape_VoxelVolume"),
           "Custom" = Clean_String(input$custom_predictors),
    )
  })
  
  ### Describe the data ##################
  output$distPlot1 <- renderPlot({
    x <- datasetInput()
    nbins = ifelse(input$N<=500, 
                   round(input$N/20,0), 
                   round(100,0))
    times_event = x[x$event==1, "time"]
    bins <- seq(min(times_event), 
                max(times_event), 
                length.out = nbins)
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
    plot_len = ifelse(length(predict.factors)>10, 10, length(predict.factors))
    predict.factors.plot <- predict.factors[1:plot_len]
    ggpairs(x[,predict.factors.plot])
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
                  cv_number = input$k_outer,  
                  seed_to_fix = input$randomseed_validation)
  })
  
  Ens2_train_on_all <- reactive({
    x <- datasetInput()
    predict.factors <- predictfactors()
    
    method_2A_train(x, predict.factors, 
                    fixed_time = input$fixed_time, 
                    internal_cv_k = input$k_inner, 
                    seed_to_fix = input$randomseed_validation)
  })
  
  Ens3_train_on_all <- reactive({
    x <- datasetInput()
    predict.factors <- predictfactors()
    
    method_3_train(x, predict.factors, 
                   fixed_time = input$fixed_time,  
                   internal_cv_k = input$k_inner, 
                   seed_to_fix = input$randomseed_validation)
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
                  cv_number = input$k_outer, 
                  internal_cv_k = input$k_inner,
                  seed_to_fix = input$randomseed_validation)
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
                 cv_number = input$k_outer, 
                 internal_cv_k = input$k_inner,
                 seed_to_fix = input$randomseed_validation)
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
                 cv_number = input$k_outer, 
                 internal_cv_k = input$k_inner,
                 seed_to_fix = input$randomseed_validation)
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
                cv_number = input$k_outer, 
                internal_cv_k = input$k_inner,
                seed_to_fix = input$randomseed_validation)
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
