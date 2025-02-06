library(shiny)
library(lubridate)
library(tidyverse)
library(JMbayes2)
library(shinythemes)  # For themes
library(bslib)
load("data/jointFit.RData")

# load git data ####
# github_token <- Sys.getenv("GITHUB_PAT")
# git_url <- Sys.getenv("GITHUB_URL") - https://raw.githubusercontent.com/doratiass/private_data/main/6_min/jointFit.RData
# tmp <- tempfile()
# response <- GET(github_url, authenticate("doratiass", github_token))
# writeBin(content(response, "raw"), tmp)
# load(tmp)

ui <- navbarPage(
  theme = shinytheme("flatly"),  # Choose a theme
  collapsible = TRUE,
  fillable_mobile = TRUE,
  "Mortality Prediction from 6MWT",
  selected = "Mortality Risk",
  # Tab 1: about
  tabPanel(
    "About",
    about_card
  ),
  # Tab 2: Mortality Risk
  tabPanel(
    "Mortality Risk",
    sidebarLayout(
      sidebarPanel(
        # Basic patient info
        h3("Patient Information"),
        selectInput("gender", "Gender", choices = c("Male", "Female")),
        dateInput("dob", "Date of Birth", value = "1960-01-01"),
        # checkboxInput("ischemic", "Ischemic Etiology", value = TRUE),
        selectInput("nyha", "NYHA", choices = c("I", "II", "III", "IV", "Undetermined")),
        
        # UI for adding/removing multiple 6MWT measurements
        hr(),
        h4("6MWT Measurements"),
        helpText("Add measurements for 6-Minute Walk Test (6MWT)."),
        
        actionButton("addMeas", "Add Measurement", class = "btn btn-primary"),
        br(), br(),
        
        # Placeholder for the multiple measurements inputs
        uiOutput("measInputs"),
        
        hr(),
        actionButton("goButton", "Generate Prediction", class = "btn btn-success btn-lg btn-block"),
        width = 4
      ),
      
      mainPanel(
        h3("Prediction Results"),
        plotOutput("dynPredPlot", height = "500px"),
        hr(),
        h4("Graph Explanation"),
        htmlOutput("graphExplanation")  
      )
    )
  )
)

server <- function(input, output, session) {
  # ------------------------------------------------------------------
  # 1) Reactive data storage for multiple measurements, with row_id
  # ------------------------------------------------------------------
  # We'll store 6MWT info in a reactive data frame with columns:
  # row_id, x6mw_dist_meter, meas_date
  
  measurements <- reactiveVal(
    data.frame(
      row_id           = integer(0),
      x6mw_dist_meter  = numeric(0),
      meas_date        = as.Date(character(0)),
      stringsAsFactors = FALSE
    )
  )
  
  # Utility function to generate a new row_id that is unique
  get_next_row_id <- function() {
    df <- measurements()
    if (nrow(df) == 0) {
      return(1L)
    } else {
      return(max(df$row_id) + 1L)
    }
  }
  
  # When user clicks "Add Measurement", append a new row with defaults
  observeEvent(input$addMeas, {
    df <- measurements()
    new_row <- data.frame(
      row_id          = get_next_row_id(),
      x6mw_dist_meter = 400,
      meas_date       = Sys.Date(),
      stringsAsFactors = FALSE
    )
    measurements(bind_rows(df, new_row))
  })
  
  # ------------------------------------------------------------------
  # 2) Dynamic UI for each measurement row
  # ------------------------------------------------------------------
  # For each row in our measurements data frame, we create:
  # - A numericInput for x6mw_dist_meter
  # - A dateInput for meas_date
  # - An actionButton for removing that row
  output$measInputs <- renderUI({
    df <- measurements()
    
    if (nrow(df) == 0) {
      return(helpText("No measurements added yet. Click 'Add Measurement' above."))
    }
    
    # For each row, build a fluidRow with the three inputs
    lapply(seq_len(nrow(df)), function(i) {
      this_row_id <- df$row_id[i]
      dist_val    <- df$x6mw_dist_meter[i]
      date_val    <- df$meas_date[i]
      
      fluidRow(
        column(5,
               numericInput(
                 inputId = paste0("dist_", this_row_id),
                 label   = paste("6MW (m) ID =", this_row_id),
                 value   = dist_val,
                 min     = 0,
                 max     = 800
               )
        ),
        column(5,
               dateInput(
                 inputId = paste0("date_", this_row_id),
                 label   = "Date",
                 value   = date_val
               )
        ),
        column(2,
               actionButton(
                 inputId = paste0("removeMeas_", this_row_id),
                 label   = "X",
                 class   = "btn btn-danger btn-sm"
               )
        )
      )
    })
  })
  
  # ------------------------------------------------------------------
  # 3) Observe changes to each measurement's inputs
  # ------------------------------------------------------------------
  # We'll watch for:
  # - numericInput changes
  # - dateInput changes
  # - remove button clicks
  
  # A generic observer that updates the data frame for ANY row’s
  # distance or date changes:
  observe({
    df <- measurements()
    if (nrow(df) == 0) return(NULL)
    
    # For each row, retrieve updated inputs
    for (i in seq_len(nrow(df))) {
      rid <- df$row_id[i]
      dist_val <- input[[paste0("dist_", rid)]]
      date_val <- input[[paste0("date_", rid)]]
      
      if (!is.null(dist_val) && !is.null(date_val)) {
        df$x6mw_dist_meter[i] <- dist_val
        df$meas_date[i]       <- date_val
      }
    }
    measurements(df)
  })
  
  # For removing a row, we set up an observer for each row’s remove button
  # using an lapply. Another approach is a single observeEvent with
  # pattern-matching input IDs, but the below is more explicit.
  observe({
    df <- measurements()
    if (nrow(df) == 0) return(NULL)
    
    # For each row, watch removeMeas_rowid
    for (i in seq_len(nrow(df))) {
      rid <- df$row_id[i]
      
      # Use isolate or an if-check to see if the button was pressed
      remove_btn_id <- paste0("removeMeas_", rid)
      # If remove button doesn't exist yet, skip
      if (is.null(input[[remove_btn_id]])) next
      
      # We want to trigger whenever input[[remove_btn_id]] changes
      observeEvent(input[[remove_btn_id]], {
        # Remove the row where row_id == rid
        current_df <- measurements()
        new_df <- filter(current_df, row_id != rid)
        measurements(new_df)
      }, ignoreInit = TRUE)
    }
  })
  
  # ------------------------------------------------------------------
  # 4) Reaction to "Generate Prediction" button
  # ------------------------------------------------------------------
  # We build the final data frame for dynamic prediction
  # and create the plot.
  
  inputData <- eventReactive(input$goButton, {
    req(input$dob)
    df_meas <- measurements()
    
    validate(
      need(nrow(df_meas) > 0,
           "Please add at least one 6MWT measurement before generating prediction.")
    )
    
    first_meas_date <- min(df_meas$meas_date)
    
    # Build a data frame with required columns:
    # id, gender, age, ischemic_etiology, time, x6mw_dist_meter
    # "time" is months from the earliest measurement date
    df_out <- df_meas %>%
      mutate(
        id    = 1,
        gender = factor(input$gender, levels = c("Female", "Male")),
        age   = interval(ymd(input$dob), first_meas_date) / years(1),
        # ischemic_etiology = as.logical(input$ischemic),
        nyha = factor(input$nyha, c("II", "I", "III", "IV", "Undetermined")),
        time = round(interval(first_meas_date, meas_date) / months(1))
      ) %>%
      select(id, gender, age, nyha, time, x6mw_dist_meter) %>% #ischemic_etiology
      arrange(time)  
    
    df_out
  })
  
  # ------------------------------------------------------------------
  # 5) Create the dynamic prediction plot
  # ------------------------------------------------------------------
  output$dynPredPlot <- renderPlot({
    df <- inputData()
    req(df)
    
    t0 <- max(df$time) + 0.1
    
    dp <- plot_dyn_pred(jointFit, df, id = 1, t0 = t0)
    plot(dp)
  })
  
  output$graphExplanation <- renderUI({
    HTML(
      "<p>This graph illustrates the relationship between 6-minute walk test (6MWT) measurements and the predicted hazard of an event (e.g., mortality) over time, based on a joint modeling approach combining a linear mixed model (LMM) and a survival model.</p>
     <ul>
       <li><strong>6MWT Trajectory:</strong> The points on the left represent actual 6MWT distances (in meters) recorded during the follow-up period for the patient. The blue line shows the predicted trajectory based on the Linear Mixed Model (LMM), and the shaded blue area represents the confidence interval.</li>
       <li><strong>End of Follow-Up:</strong> The dotted vertical line marks the end of the follow-up period (e.g., the last 6MWT measurment). Measurements and predictions to the left of this line represent observed and modeled values during follow-up.</li>
       <li><strong>Predicted Hazard:</strong> The red line represents the predicted hazard (risk of an event, such as mortality) over time, as derived from the joint model. The shaded red area indicates the confidence interval for the hazard prediction. The increasing trend after the follow-up suggests rising risk.</li>
     </ul>"
    )
  })
}

shinyApp(ui = ui, server = server)