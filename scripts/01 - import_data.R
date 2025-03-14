# ============================================================================ #
# Data Preparation and Cleaning Script for 6-Minute Walk Test & Mortality
#
# This script performs the following tasks:
#   1. Loads required functions and libraries.
#   2. Reads in mortality and clinical data from Excel files.
#   3. Filters out excluded subjects based on criteria.
#   4. Cleans and recodes variables (both numeric and categorical).
#   5. Merges mortality data with clinical data and computes follow-up times.
#   6. Generates a baseline summary table (Table 1) and saves it as an HTML file.
#   7. Saves the final cleaned dataset and scaling factors for future analyses.
#
# Note: This project evaluates the impact of repeated 6-minute walk test 
#       measurements on mortality prediction using joint modeling methods.
# ============================================================================ #

# ============================================================================ #
# Load Functions and Libraries
# ============================================================================ #
source("scripts/00 - funcs.R")  # Load custom functions (e.g., merge_rows, label_get, etc.)

# Load required packages
library(tidyverse)    # Data manipulation and visualization
library(readxl)       # Reading Excel files
library(zoo)          # Working with time series (e.g., na.approx)
library(survival)     # Survival analysis

# ============================================================================ #
# Definitions -----------------------------------------------------------------
# ============================================================================ #
# -------------------------------------------------------------------------- #
## Define Exclusion IDs and Scaling Parameters ####
# -------------------------------------------------------------------------- #
# Read mortality Excel file, clean column names, and filter out rows that have an "exclude" flag 
# or a specific mortality_status (9999). Then, extract the 'id_enc' values to exclude.
ex_id <- readxl::read_xlsx("data/CHF_mortality.xlsx") %>%
  janitor::clean_names() %>%
  filter(!is.na(exclude) | mortality_status == 9999) %>%
  pull(id_enc)

# -------------------------------------------------------------------------- #
## Define Keywords and Variable Sets ####
# -------------------------------------------------------------------------- #
# Cardiac-related keywords (used later to flag cardiac diagnoses)
cardiac_keywords <- c("cardiac", "myocardial", "heart", "coronary", "angina", 
                      "arrhythmia", "cardiomyopathy", "ischemic", "infarction",
                      "valve", "aortic", "mitral", "pulmonary", "tricuspid",
                      "congestive", "chf", "heart failure", "hf", "systolic",
                      "pci", "cabg", "stent", "bypass", "defibrillator", "icd",
                      "Atrial", "ihd", "PTCA", "pacemaker", "LV D", "af", "pfo",
                      "Pericarditis", "pericardial", "pericardium", "pericardial",
                      "chest pain", "angina", "dyspnea", "palpitations", "edema",
                      "syncope", "cardiomegaly", "cardiomegaly", "cardiomegaly",
                      "Ventricular", "ventricular", "ventricle", "ventricles",
                      "myocarditis", "myocardial", "myocardium", "myocardial",
                      "NYHA", "לב")

# Variables to be included in the Seattle model (or other specific analyses)
seattle_vars <- c("age", "gender", "ischemic_etiology", "sbp", "nyha", 
                  "fusid", "k_spare", "allopurinol",
                  "bmi", "statin", "ace_inhibitor", "arbs",
                  "beta_blocker", "na", "hgb")

# ============================================================================ #
# Import and Clean Data -------------------------------------------------------
# ============================================================================ #

# -------------------------------------------------------------------------- #
## Mortality Data ####
# -------------------------------------------------------------------------- #
mortality <- readxl::read_xlsx("data/CHF_mortality.xlsx") %>%
  janitor::clean_names() %>%
  # Exclude IDs flagged earlier
  filter(!(id_enc %in% ex_id)) %>%
  # Create a follow-up date:
  #   - For deceased (mortality_status == 1): use the 'last_fup' date.
  #   - For survivors: assign a fixed censoring date ("2024-08-01 00:00:00").
  mutate(fup_date = as_datetime(ifelse(mortality_status == 1, 
                                       as_datetime(last_fup), 
                                       as_datetime("2024-08-01 00:00:00")))) %>%
  # Select and rename relevant columns
  select(id = id_enc, fup_date, mortality_status)


# -------------------------------------------------------------------------- #
## Clinical Data ####
# -------------------------------------------------------------------------- #
raw_df <- read_excel("data/CHF_final_wide.xlsx") %>%
  janitor::clean_names() %>%
  # Remove excluded subjects
  filter(!(id %in% ex_id)) %>%
  # Convert columns containing "_res" to numeric type
  mutate_at(vars(contains("_res")), as.numeric) %>%
  # Recode certain categorical variables (e.g., "כן" -> TRUE, "לא" -> FALSE)
  mutate_at(c("edema", "jvp", "night_dyspnea",
              "alcohol", "pnd", "exercise_dyspnea",
              "palpitations", "osa", 
              "dm", "htn", "lipid", "drugs"),
            function(x) case_when(
              x == "כן" ~ TRUE,
              x == "לא" ~ FALSE)) %>%
  # Clean and validate numerical variables by setting implausible values to NA
  mutate(
    x6mw_dist_meter = ifelse(x6mw_dist_meter > 800 | x6mw_dist_meter < 10, NA, x6mw_dist_meter),
    x6mw_hr_pre     = ifelse(x6mw_hr_pre > 300 | x6mw_hr_pre < 40, NA, x6mw_hr_pre),
    x6mw_hr_post    = ifelse(x6mw_hr_post > 300 | x6mw_hr_post < 40, NA, x6mw_hr_post),
    x6mw_sat_pre    = ifelse(x6mw_sat_pre > 100 | x6mw_sat_pre < 60, NA, x6mw_sat_pre),
    x6mw_sat_post   = ifelse(x6mw_sat_post > 100 | x6mw_sat_post < 60, NA, x6mw_sat_post),
    ht              = ifelse(ht < 100, NA, ht),
    wt              = ifelse(wt < 30 | wt > 200, NA, wt),
    hr              = ifelse(hr < 40 | hr > 300, NA, hr),
    sbp             = ifelse(sys < 50 | sys > 250, NA, sys),
    dbp             = ifelse(dias < 40 | dias > 180, NA, dias),
    jvp_cm          = ifelse(jvp_cm < 0 | jvp_cm > 20, NA, jvp_cm),
    alb             = ifelse(alb_res > 10, NA, alb_res),
    cr              = ifelse(cr_res > 10, NA, cr_res),
    hgb             = ifelse(hgb_res > 20 | hgb_res < 3, NA, hgb_res),
    na              = ifelse(na_res > 150 | na_res < 120, NA, na_res),
    alkp            = ifelse(alkp_res > 2000, NA, alkp_res),
    ggt             = ifelse(ggt_res > 1000, NA, ggt_res),
    bil             = ifelse(bil_res > 15, NA, bil_res),
    ua              = ifelse(ua_res > 100, NA, ua_res),
    k               = ifelse(k_res > 15, NA, k_res),
    glu             = ifelse(glu_res > 500 | glu_res < 20, NA, glu_res),
    
    # Recode categorical variables:
    # Convert gender from Hebrew character "ז" to "Male" (others become "Female")
    gender          = factor(ifelse(gender == "ז", "Male", "Female"),
                             levels = c("Female", "Male")),
    # Recode family status using provided Hebrew codes
    fam_status      = factor(case_when( 
      fam_status == "ר" ~ "Single",
      fam_status == "נ" ~ "Married",
      fam_status == "ג" ~ "Divorced",
      fam_status == "א" ~ "Widowed"
    )),
    # Recode healthcare provider (kupa) into standardized categories
    kupa            = factor(case_when(
      kupa %in% c("כללית", "שרותי בריאות כללית",
                  "כללית מושלם - שירותי בריאות כל",
                  "לית                         03")  ~ "Clalit",
      kupa %in% c("מכבי", "שב\"ן מכבי זהב",
                  "קופת חולים מכבי") ~ "Maccabi",
      kupa %in% c("מאוחדת", "קופת חולים מאוחדת") ~ "Meuhedet",
      kupa %in% c("לאומית", "קופת חולים לאומית", "לאומית עתיד") ~ "Leumit",
      TRUE ~ "Other"
    )),
    # Recode visit type: classify several descriptions as "Virtual" or "Regular"
    visit_type      = factor(case_when(
      visit_type == "ביקור וירטואלי (וידאו / טלפון)" ~ "Virtual",
      visit_type == "נרשם ללא נוכחות המטופל"           ~ "Virtual",
      visit_type == "רגיל"                                ~ "Regular",
      visit_type == "שיחה טלפונית"                        ~ "Virtual",
      visit_type == "שיחת וידאו"                          ~ "Virtual",
      TRUE ~ as.character(visit_type)
    )),
    # Set NYHA class as a factor with ordered levels
    nyha            = factor(nyha, levels = c("II", "I", "III", "IV", "Undetermined")),
    # Recode weight category into standardized groups
    weight_category = factor(case_when(
      weight_category %in% c("תקין", "טוב", "סביר", "תקין (מעט רזה)") ~ "Normal",
      weight_category %in% c("רזה", "רזה .", "רזון", "רזה מאד", "מעט רזה", "ירוד",
                             "רזה (ירד במשקל לאחר התקף הלב להערכתו מעל 10 ק\"ג)") ~ "Thin",
      weight_category %in% c("כחקטי", "קכקטי", "קכקסיה", "כחקטית") ~ "Cachexia",
      weight_category %in% c("עודף משקל", "עודף משקל קל", "השמנה יתירה", "עודף משקל ניכר",
                             "עודף משקל קשה", "השמנה ממארת", "השמנה", "השמנה מרכזית",
                             "עודף משקל בינוני", "השמנה ניתירה", "השמנה קלה", "עדיין השמנה יתירה", "עודף משקל .") ~ "Overweight",
      TRUE ~ as.character(weight_category)
    )),
    # Recode heart failure class as a factor with levels A-D
    hf_class        = factor(case_when(
      hf_class == "A" ~ "A",
      hf_class == "B" ~ "B",
      hf_class == "C" ~ "C",
      hf_class == "D" ~ "D"
    ), levels = c("A", "B", "C", "D")), 
    # Convert edema_detail to a factor (further decisions needed on categorization)
    edema_detail    = factor(edema_detail),
    # Flag for pillows-related dyspnea based on keywords in the free-text 'pillows_reason'
    pillows_reason_dyspnea = grepl("BIPAP|CIPAP|CPAP|נשימה|חמצן|משתעל", pillows_reason),
    # Impute jvp_cm: if missing and jvp is false then assume a default value of 3
    jvp_cm          = case_when(
      is.na(jvp_cm) & !jvp ~ 3,
      TRUE ~ as.numeric(jvp_cm)
    ),
    # Adjust smoke_stop_year if recorded in a two-digit format (adding 1900 if needed)
    smoke_stop_year = ifelse(smoke_stop_year < 1900, smoke_stop_year + 1900, smoke_stop_year),
    # Flag cardiac diagnosis if any of the cardiac keywords appear in the 'diag' field
    diag_cardiac    = grepl(paste0(cardiac_keywords, collapse = "|"), diag, ignore.case = TRUE),
    # Count the number of diagnoses by splitting the 'diag' field on '^'
    diag_num        = sapply(strsplit(diag, "\\^"), length),
    # Define ischemic etiology based on the presence of the Hebrew word for "ischemic"
    ischemic_etiology = grepl("\\bאיסכמי\\b", etiology_1, ignore.case = TRUE) & 
      !grepl("\\bלא איסכמי\\b", etiology_1, ignore.case = TRUE),
    # Create flags for various medications using pattern matching
    sglt2           = grepl("Jardiance|Empagliflozin|Forxiga|Dapagliflozin|Xigduo|Glyxambi|Empagliflozin", medications, ignore.case = TRUE),
    glp1            = grepl("Victoza|Trulicity|Ozempic|Saxenda|Xultophy|Liraglutide|Dulaglutide|Semaglutide", medications, ignore.case = TRUE),
    fusid           = grepl("\\b[Ff]usid\\b|\\b[Ff]urosemide\\b", medications, ignore.case = TRUE),
    torsemide       = grepl("\\bTorsemide\\b|\\bTorsemid\\b", medications, ignore.case = TRUE),
    bumetanide      = grepl("\\bBumetanide\\b", medications, ignore.case = TRUE),
    metolazone      = grepl("\\bMetolazone\\b|\\bZaroxylin\\b|\\bZaroxillin\\b|\\bZaroxyllyn\\b", medications, ignore.case = TRUE),
    k_spare         = grepl("\\bSpironolactone\\b|\\bAldactone\\b|\\bAldospirone\\b|\\bEplerenone\\b|\\bInspra\\b", medications, ignore.case = TRUE),
    thiazide        = grepl("\\bHydrochlorothiazide\\b", medications, ignore.case = TRUE),
    allopurinol     = grepl("\\bAllopurinol\\b|\\bAlloril\\b|\\bZylol\\b|\\bZylori\\b", medications, ignore.case = TRUE),
    beta_blocker    = grepl("Bisoprolol|Concor|Cardiloc|Metoprolol|Lopresor|Neobloc|Metopress|Carvedilol|Dimitone|Atenolol|Normiten|Normalol|Propranolol|Deralin|Prolol", medications, ignore.case = TRUE),
    anticoagulant   = grepl("Xarelto|Rivaroxaban|Coumadin|Warfarin|Eliquis|Apixaban|Pradaxa|Dabigatran|Sintrom|Acenocoumarol", medications, ignore.case = TRUE),
    antiplatlet     = grepl("Plavix|Clopidogrel|Brilinta|Ticagrelor|Effient|Prasugrel|Aspirin|Acetylsalicylic|Cardioaspirin|Cardiopirin|Aspocid|Aspicor|Aspimax|Aspirin", medications, ignore.case = TRUE),
    ace_inhibitor   = grepl("Tritace|Enalapril|Ramipril|Lisinopril|Captopril|Perindopril|Quinapril|Trandolapril|Moexipril|Benazepril", medications, ignore.case = TRUE),
    arbs            = grepl("Losartan|Ocsaar|Valsartan|Diovan|Irbesartan|Candesartan|Atacand|Olmesartan|Telmisartan|Eprosartan|Azilsartan", medications, ignore.case = TRUE),
    statin          = grepl("Atorvastatin|Simvastatin|Rosuvastatin|Pravastatin|Fluvastatin|Lovastatin|Pitavastatin", medications, ignore.case = TRUE),
    insulin         = grepl("Lantus|Levemir|Toujeo|Novolog|Humalog|Apidra|Fiasp|Novorapid|Humulin|Novolin|Tresiba|Basaglar", medications, ignore.case = TRUE),
    metformin       = grepl("Metformin|Glucophage|Glucovance|Janumet|Jentadueto|Kombiglyze|Synjardy|Xigduo|Invokamet|Qtern|Steglatro|Segluromet|Glyxambi|Steglujan", medications, ignore.case = TRUE),
    dpp4            = grepl("Januvia|Onglyza|Tradjenta|Nesina|Galvus|Zomelis|Qtern|Steglatro|Segluromet|Glyxambi|Steglujan", medications, ignore.case = TRUE),
    glinide         = grepl("Starlix|Prandin|Novonorm", medications, ignore.case = TRUE),
    sulfonylurea    = grepl("Glibenclamide|Glyburide|Glimepiride|Gliclazide|Tolbutamide|Diamicron|Amaryl|Daonil", medications, ignore.case = TRUE),
    # Count the number of medications by splitting the 'medications' string on '^'
    med_num         = sapply(strsplit(medications, "\\^"), length),
    # Flag if the 6MWT date matches the clinical date
    date_ok         = x6mw_date == clin_date
  ) %>%
  # Process free-text fields related to family history, hospitalizations, and emergency visits:
  mutate(
    # For family history, set detail and recode based on known negative responses
    family_history_detail = family_history, 
    family_history = case_when( 
      family_history %in% c("לא ידוע", "לא", "שולל", "אין", "ללא", "לא ידוע על איס\"ל במשפחה.",
                            "שולל מחלת לב במשפחה", "שולל מחלות לב או מוות פתאומי",
                            "שולל מחלת לב במשפחה.", "ללא רקע משפחתי", "לא ידועיש יל\"ד משפחתי.",
                            "לאלאב היתה שחפת ככל הנראה", "לא ידוע (ניצול שואה)", "שוללת",
                            "ללא סיפור משפחתי", "ללא חולי לב במשפחה", "הורים נפטרו משחפת.",
                            "שולל מחלות לב במשפחה", "ללא רקע") ~ FALSE,
      TRUE ~ !is.na(family_history)
    ),
    family_history_detail = ifelse(family_history, family_history_detail, NA_character_),
    .after = family_history
  ) %>%
  mutate(
    # For hospitalizations, create a flag and detail variable similarly
    hosp_detail = hosp, 
    hosp = case_when(
      hosp_detail %in% c("לא", "לאלא", "לא;", "לא. עבר ניתוח של טריגר פינגר באשפוז יום",
                         "לא.עקירה כירורגית של שן", "לא רלוונטי", "לא (אלקטיבי)",
                         "לא - למעט מרפאה המטולוגית", "לא (אלקטיבי לICD)", 
                         "לא (אשפוז יום).", "לא (היה מאושפז לפני כ8 חודשים בקרדיולוגיה ועבר השתלת ICD)",
                         "לא (למעט אלקטיבי)", "לא (למעט לביופסיה)", "לא (מטופל באשפוז יום המטולוגי)",
                         "לא (שדרוג לCRTD אלקטיבי.", "לא בארבע השנים האחרונות", "לא היה אשפוז",
                         "לא למעט אשפוז יום", "לא למעט ניתוח", "לא מעבר לאשפוז יום", 
                         "לא. היתה מאושפזת לפני שנתיים בקרדיולוגיה,  אז אבחנו לה CHF") ~ FALSE,
      TRUE ~ !is.na(hosp_detail)
    ),
    hosp_detail = ifelse(hosp, hosp_detail, NA_character_),
    .after = hosp
  ) %>%
  mutate(
    # Process emergency room (ER) visit information similarly
    er_detail = er, 
    er = case_when( 
      er_detail %in% c("לא", "לא;", "לא רלוונטי") ~ FALSE,
      TRUE ~ !is.na(er_detail)
    ),
    er_detail = ifelse(er, er_detail, NA_character_),
    .after = er
  ) %>%
  # Recode smoking status based on number of cigarettes per day, smoking age, and stop year.
  mutate(smoke = factor(case_when(
    smoke_cig_per_day > 10              ~ ">10",
    smoke_cig_per_day > 0               ~ "1-10",
    smoke_stop_year <= year(clin_date)  ~ "past smoker",
    !is.na(smoke_age)                   ~ "past smoker",
    TRUE                                ~ "never smoked"
  ), levels = c("never smoked", "past smoker", "1-10", ">10")),
  .before = smoke_cig_per_day) %>%
  # Fix Duplicated Rows and Impute Missing Values
  # Group by subject and clinical date to merge duplicate rows using the custom 'merge_rows' function.
  group_by(id, clin_date) %>%
  summarise(across(everything(), merge_rows), .groups = 'drop') %>%
  # For each subject, arrange by clinical date and impute missing height and weight using linear approximation
  group_by(id) %>%
  arrange(clin_date) %>%
  mutate(across(c("ht", "wt"), ~ if (n() > 1 & any(!is.na(ht))) {
    na.approx(.x, rule = 2)
  } else {
    .x
  })) %>%
  # Fill missing height/weight and laboratory values (columns containing '_res') in a downward direction
  fill(all_of(c("ht", "wt")), .direction = "downup") %>%
  fill(contains("_res"), .direction = "down") %>%
  ungroup() %>%
  # Calculate body mass index (BMI)
  mutate(bmi = wt / ((ht/100)^2)) %>%
  # Keep only records with matching dates, non-missing walk test date/distance, and from 2014 onwards
  filter(date_ok,
         !is.na(x6mw_date),
         !is.na(x6mw_dist_meter),
         year(x6mw_date) >= 2014)

# ============================================================================ #
# Create the Final Cleaned Dataset --------------------------------------------
# ============================================================================ #
clean_df_identified <- raw_df %>%
  # Merge the clinical data with mortality data by subject ID
  left_join(mortality, by = "id") %>%
  arrange(id, x6mw_date) %>%
  group_by(id) %>%
  mutate(
    # Scale the 6-minute walk test distance and age if necessary
    x6mw_dist_meter = x6mw_dist_meter / six_min_scale,
    age             = age / age_scale,
    # Define the first visit date and compute time (in months) from the first visit
    first_visit     = first(x6mw_date),
    # time            = round(time_length(x6mw_date - first_visit, "months")),
    time_raw            = time_length(x6mw_date - first_visit, "months"),
    time = case_when(
      time_raw <= 0 ~ 0,
      time_raw <= 1  ~ 1,
      TRUE         ~ round(time_raw)
    ),
    # Compute the time difference (in months) between the follow-up date and the first visit
    time_diff_months = time_length(fup_date - first_visit, "months"),
    # Define the last follow-up date based on mortality status and maximum follow-up time
    last_fup_date = case_when(
      mortality_status == 1 & time_diff_months <= max_fup_time ~ fup_date,
      mortality_status == 1 & time_diff_months > max_fup_time  ~ first_visit + months(max_fup_time),
      mortality_status == 0                                     ~ pmin(fup_date, first_visit + months(max_fup_time), na.rm = TRUE)
    ),
    # Adjust mortality status to 0 (censored) if the event occurred after the max follow-up time
    mortality_status = case_when(
      mortality_status == 1 & time_diff_months > max_fup_time ~ 0,
      TRUE ~ mortality_status
    ),
    # Calculate the follow-up time in months; if zero, set to 1
    fup_time = round(time_length(last_fup_date - first_visit, "months")),
    fup_time = ifelse(fup_time == 0, 1, fup_time)
  ) %>%
  ungroup() %>%
  # Keep only records within the maximum follow-up period
  filter(time <= max_fup_time) %>%
  # Select the final set of variables for analysis
  select(id, mortality_status, fup_time, time, x6mw_dist_meter, all_of(seattle_vars)) %>%
  group_by(id, mortality_status, fup_time, time) %>%
  summarise(
    across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),  # Take mean for numeric variables
    across(where(~ !is.numeric(.)), \(x) first(x))          # Take the first value for non-numeric variables
  ) %>%
  ungroup()

# -------------------------------------------------------------------------- #
## Anonymize the Data and Create a De-Identified Dataset ####
# -------------------------------------------------------------------------- #
# Create a mapping table
unique_ids <- unique(clean_df_identified$id)
anon_ids <- seq_along(unique_ids)
id_map <- data.frame(Original_ID = unique_ids, Anonymized_ID = anon_ids)

# Replace original IDs with anonymized IDs
clean_df <- id_map %>%
  right_join(clean_df_identified, by = c("Original_ID" = "id")) %>%
  select(-Original_ID) %>%
  rename(id = Anonymized_ID)

# Save the mapping table
write_csv(id_map, "data/id_mapping.csv")

# -------------------------------------------------------------------------- #
## Generate Baseline Summary Table (Table 1) ####
# -------------------------------------------------------------------------- #
tbl_1 <- clean_df %>%
  group_by(id) %>%
  # Compute the number of 6MWT tests per subject and the overall change from first to last measurement
  mutate(
    n_6mw      = length(unique(x6mw_dist_meter)),
    x6mw_change = x6mw_dist_meter[which.max(time)] - x6mw_dist_meter[which.min(time)]
  ) %>%
  ungroup() %>%
  # Filter to baseline records (time == 0)
  filter(time == 0) %>%
  # Create a categorical variable for 6MWT distance based on quantiles
  mutate(x6mw_dist_meter_q = cut(x6mw_dist_meter,
                                 c(0, quantile(x6mw_dist_meter, c(0.333, 0.667)), max(x6mw_dist_meter)),
                                 labels = c("Low percentile", "Middle percentile", "High percentile"))) %>%
  # Remove identifiers and time variable from the summary table
  select(-c(id, time, x6mw_change)) %>%
  # Rename variables using the custom label_get function for presentation
  rename_all(function(x) sapply(x, label_get, USE.NAMES = FALSE)) %>%
  # Create a summary table stratified by 6MWT distance quantiles using gtsummary
  tbl_summary(by = x6mw_dist_meter_q, 
              type = list(label_get("n_6mw") ~ "continuous"),
              statistic = list(all_continuous() ~ "{mean} ({sd})"),
              value = list(label_get("gender") ~ "Male"),
              missing = "no") %>%
  add_overall() %>%
  add_n(statistic = "{N_miss} ({p_miss})") %>%
  modify_header(n = "**Missing**") %>%
  add_p()

# Save the summary table as an HTML file
gt::gtsave(as_gt(tbl_1), file = "export/tbl_1.html")

# -------------------------------------------------------------------------- #
## Save the Final Cleaned Data and Clean Up Workspace ####
# -------------------------------------------------------------------------- #
# Remove temporary variables to clean the workspace
rm(cardiac_keywords, ex_id, mortality, raw_df)

# Save the final cleaned data and scaling factors for future analyses
save(clean_df, age_scale, six_min_scale, file = "data/CHF_data.RData")
