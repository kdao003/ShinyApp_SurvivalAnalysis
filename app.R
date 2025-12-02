#install.packages(c("shiny", "survival", "survminer", "dplyr", "ggplot2", "readr"))

# ==========================================================
# Shiny App for Survival Analysis (Stable Core Version)
# ==========================================================

library(shiny)
library(survival)
library(survminer)
library(readr)
library(dplyr)
library(rsconnect)
library(googlesheets4)

# Avoid requiring authentication for public sheets
gs4_deauth()

ui <- fluidPage(
  titlePanel("Survival Analysis App"),
  
  sidebarLayout(
    sidebarPanel(
      radioButtons(
        "data_source", "Select Data Source:",
        choices = c("Upload CSV" = "csv", "Import Google Sheet" = "sheet"),
        selected = "csv"
      ),
      conditionalPanel(
        condition = "input.data_source == 'csv'",
        fileInput("file", "Upload CSV File", accept = ".csv")
      ),
      conditionalPanel(
        condition = "input.data_source == 'sheet'",
        textInput("sheet_url", "Enter Google Sheet URL:", placeholder = "https://docs.google.com/spreadsheets/d/.../edit"),
        actionButton("load_sheet", "Load Google Sheet")
      ),
      uiOutput("column_selector"),
      uiOutput("cox_selector"),
      actionButton("run_analysis", "Run Survival Analysis"),
      actionButton("run_cox", "Run Cox PH Analysis"),
      width = 3
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Cleaned Data", tableOutput("data_preview")),
        tabPanel("Survival Plot", plotOutput("surv_plot"), verbatimTextOutput("fit_summary")),
        tabPanel("Cox Proportional Hazards",
                 verbatimTextOutput("cox_summary"),
                 plotOutput("cox_forest"))
      )
    )
  )
)

server <- function(input, output, session) {
  
  # ---- Load Data (CSV or Google Sheet) ----
  data_input <- reactive({
    if (input$data_source == "csv") {
      req(input$file)
      df <- read_csv(input$file$datapath, show_col_types = FALSE)
    } else {
      req(input$sheet_url)
      tryCatch({
        df <- read_sheet(input$sheet_url)
      }, error = function(e) {
        showNotification("Failed to load Google Sheet. Check URL and permissions.", type = "error")
        return(NULL)
      })
    }
    
    # --- Convert list-columns to atomic vectors ---
    df <- df %>%
      mutate(across(where(is.list), ~ sapply(., function(x) if (length(x) == 0) NA else as.character(x))))
    
    df[df == ""] <- NA
    df
  })
  
  # ---- Dynamic selectors ----
  output$column_selector <- renderUI({
    req(data_input())
    df <- data_input()
    
    tagList(
      selectInput("time_col", "Select Time Column", names(df)),
      selectInput("status_col", "Select Status Column (0=censored,1=event)", names(df)),
      selectInput("group_col", "Select Group Column (optional)", 
                  choices = c("None", names(df)), selected = "None")
    )
  })
  
  # ---- Dynamic Cox variable selector ----
  output$cox_selector <- renderUI({
    req(data_input())
    df <- data_input()
    
    tagList(
      helpText("Select variables for Cox model (besides time/status)"),
      checkboxGroupInput("cox_vars", "Select Predictor Variables", 
                         choices = names(df))
    )
  })
  
  # ---- Data Preview ----
  output$data_preview <- renderTable({
    req(data_input())
    head(data_input())
  })
  
  # ---- Survival analysis ----
  surv_results <- eventReactive(input$run_analysis, {
    req(input$time_col, input$status_col, input$group_col)
    df <- data_input()
    
    #make character
    time_col <- as.character(input$time_col)
    status_col <- as.character(input$status_col)
    group_col  <- as.character(input$group_col)
    
    # Replace blanks with NA and remove rows with missing values in key columns
    key_cols <- if (group_col != "None") {
      c(time_col, status_col, group_col)
    } else {
      c(time_col, status_col)
    }
    
    df <- df %>%
      mutate(across(where(is.character), ~ na_if(., " "))) %>%
      filter(if_all(all_of(key_cols), ~ !is.na(.)))
    
    #validate(need(nrow(df) > 0, "No valid rows remaining after removing NAs."))
    
    #Remove the rows that have an NA value
    #key_cols <- c(time_col, status_col, group_col)
    
    #df <- df %>%
      #filter(if_all(all_of(key_cols), ~ !is.na(.)))
    
    #validate(need(nrow(df) > 0, "No valid rows remaining after removing NAs."))
    
    #Create Survival Object
    if(group_col != "None")
    {
      df[[time_col]]   <- as.numeric(df[[time_col]])
      df[[status_col]] <- as.numeric(df[[status_col]])
      df[[group_col]] <- as.factor(df[[group_col]])
      
      surv_obj <- as.formula(paste0("Surv(", time_col, ",", status_col, ") ~", group_col))
      fit <- eval(substitute(survfit(surv_obj, data = df, type = "kaplan-meier")))
    }
    else
    {
      df[[time_col]]   <- as.numeric(df[[time_col]])
      df[[status_col]] <- as.numeric(df[[status_col]])
      
      surv_obj <- as.formula(paste0("Surv(", time_col, ",", status_col, ") ~", 1))
      fit <- eval(substitute(survfit(surv_obj, data = df, type = "kaplan-meier")))
    }
    
    # Create survival object
    #surv_obj <- as.formula(paste0("Surv(", time_col, ",", status_col, ") ~", group_col))
    
    #fit <- eval(substitute(survfit(surv_obj, data = df, type = "kaplan-meier")))
    
    # Fit model safely (explicit argument passing)
    #if (group_col != "None") {
    #   fit <- survfit(surv_obj ~ group_var)
    # } else {
    #  fit <- survfit(surv_obj ~ 1)
    #}
    
    list(fit = fit, df = df, group_col = group_col)
  })
  
  # ---- Plot ----
  output$surv_plot <- renderPlot({
    req(surv_results())
    res <- surv_results()
    
    # Ensure grouping column is treated as factor
    if (res$group_col != "None") {
      res$df[[res$group_col]] <- as.factor(res$df[[res$group_col]])
    }
    
    ggsurvplot(
      res$fit,
      data = res$df,
      risk.table = TRUE,
      conf.int = FALSE,
      pval = TRUE,
      ggtheme = theme_minimal(),
      palette = "Dark2",
      linetype = 1,           # fixed, not mapped
      linewidth = 0.8,        # fixed width
      conf.int.fill = "gray80",  # use fill instead of mapping CI by group
      conf.int.alpha = 0.25,     # transparency for CI
      conf.int.style = "ribbon"  # keep CI bands simple
    )
    
  })
  
  # ---- Survival Summary ----
  output$fit_summary <- renderPrint({
    req(surv_results())
    summary(surv_results()$fit)
  })
  
  # ---- Cox Proportional Hazards Model ----
  cox_results <- eventReactive(input$run_cox, {
    req(input$time_col, input$status_col, input$cox_vars)
    df <- data_input()
    
    time_col <- input$time_col
    status_col <- input$status_col
    predictors <- input$cox_vars
    
    # Clean and prepare data
    df <- df %>%
      mutate(across(where(is.character), ~ na_if(., " "))) %>%
      filter(if_all(all_of(c(time_col, status_col, predictors)), ~ !is.na(.)))
    
    df[[time_col]] <- as.numeric(df[[time_col]])
    df[[status_col]] <- as.numeric(df[[status_col]])
    
    # Convert character predictors to factors
    df <- df %>%
      mutate(across(all_of(predictors), ~ if (is.character(.)) as.factor(.) else .))
    
    # Fit Cox model
    cox_formula <- as.formula(
      paste0("Surv(", time_col, ",", status_col, ") ~ ", paste(predictors, collapse = " + "))
    )
    
    cox_fit <- coxph(cox_formula, data = df)
    list(fit = cox_fit, df = df)
  })
  
  # ---- Cox Model Summary ----
  output$cox_summary <- renderPrint({
    req(cox_results())
    summary(cox_results()$fit)
  })
  
  # ---- Cox Forest Plot ----
  output$cox_forest <- renderPlot({
    req(cox_results())
    res <- cox_results()
    ggforest(res$fit, data = res$df)
  })
  
}

shinyApp(ui, server)


