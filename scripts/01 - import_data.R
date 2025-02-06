source("scripts/00 - funcs.R")
library(tidyverse)
library(readxl)
library(zoo)
library(survival)

ex_id <- readxl::read_xlsx("data/CHF_mortality.xlsx") %>%
  janitor::clean_names() %>%
  filter(!is.na(exclude) |
           mortality_status == 9999) %>%
  pull(id_enc)

#readxl::excel_sheets("CHF_final.xlsx")
age_scale <- 1
six_min_scale <- 1
min_fup_time <- 6
max_fup_time <- 120

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

seattle_vars <- c("age", "gender", "ischemic_etiology", "sbp", "nyha", 
                  "fusid", "k_spare", "allopurinol",
                  "bmi", "statin", "ace_inhibitor", "arbs",
                  "beta_blocker", "na", "hgb")

# Import and clean data --------------------------------------------------------
mortality <- readxl::read_xlsx("data/CHF_mortality.xlsx") %>%
  janitor::clean_names() %>%
  filter(!(id_enc %in% ex_id)) %>%
  mutate(fup_date = as_datetime(ifelse(mortality_status == 1, 
                                       as_datetime(last_fup), 
                                       as_datetime("2024-08-01 00:00:00")))) %>%
  select(id = id_enc, fup_date, mortality_status)

raw_df <- read_excel("data/CHF_final_wide.xlsx") %>%
  janitor::clean_names() %>%
  filter(!(id %in% ex_id)) %>%
  mutate_at(vars(contains("_res")), as.numeric) %>%
  mutate_at(c("edema", "jvp", "night_dyspnea",
              "alcohol", "pnd", "exercise_dyspnea",
              "palpitations", "osa", 
              "dm", "htn", "lipid", "drugs"), function(x) case_when(
                x == "כן" ~ TRUE,
                x == "לא" ~ FALSE)) %>%
  mutate(
    # arrange ranges 
    x6mw_dist_meter = ifelse(x6mw_dist_meter > 800 | x6mw_dist_meter < 10, NA, x6mw_dist_meter),
    x6mw_hr_pre = ifelse(x6mw_hr_pre > 300 | x6mw_hr_pre < 40, NA, x6mw_hr_pre),
    x6mw_hr_post = ifelse(x6mw_hr_post > 300 | x6mw_hr_post < 40, NA, x6mw_hr_post),
    x6mw_sat_pre = ifelse(x6mw_sat_pre > 100 | x6mw_sat_pre < 60, NA, x6mw_sat_pre),
    x6mw_sat_post = ifelse(x6mw_sat_post > 100 | x6mw_sat_post < 60, NA, x6mw_sat_post),
    ht = ifelse(ht < 100, NA, ht),
    wt = ifelse(wt < 30 | wt > 200, NA, wt),
    hr = ifelse(hr < 40 | hr > 300, NA, hr),
    sbp = ifelse(sys < 50 | sys > 250, NA, sys),
    dbp = ifelse(dias < 40 | dias > 180, NA, dias),
    jvp_cm = ifelse(jvp_cm < 0 | jvp_cm > 20, NA, jvp_cm),
    alb = ifelse(alb_res > 10, NA, alb_res),
    cr = ifelse(cr_res > 10, NA, cr_res),
    hgb = ifelse(hgb_res > 20 | hgb_res < 3, NA, hgb_res),
    na = ifelse(na_res > 150 | na_res < 120, NA, na_res),
    alkp = ifelse(alkp_res > 2000, NA, alkp_res),
    ggt = ifelse(ggt_res > 1000, NA, ggt_res),
    bil = ifelse(bil_res > 15, NA, bil_res),
    ua = ifelse(ua_res > 100, NA, ua_res),
    k = ifelse(k_res > 15, NA, k_res),
    glu = ifelse(glu_res > 500 | glu_res < 20, NA, glu_res),
    # recode categorical variables
    gender = factor(ifelse(gender == "ז", "Male", "Female"),
                    levels = c("Female", "Male")),
    fam_status = factor(case_when( 
      fam_status == "ר" ~ "Single",
      fam_status == "נ" ~ "Married",
      fam_status == "ג" ~ "Divorced",
      fam_status == "א" ~ "Widowed"
    )),
    kupa = factor(case_when(
      kupa %in% c("כללית", "שרותי בריאות כללית",
                  "כללית מושלם - שירותי בריאות כל",
                  "לית                         03")  ~ "Clalit",
      kupa %in% c("מכבי","שב\"ן מכבי זהב",
                  "קופת חולים מכבי")~ "Maccabi",
      kupa %in% c("מאוחדת", "קופת חולים מאוחדת") ~ "Meuhedet",
      kupa %in% c("לאומית", "קופת חולים לאומית" ,
                  "לאומית עתיד") ~ "Leumit",
      TRUE ~ "Other"
    )),
    visit_type = factor(case_when(
      visit_type == "ביקור וירטואלי (וידאו / טלפון)" ~ "Virtual",
      visit_type == "נרשם ללא נוכחות המטופל" ~ "Virtual",
      visit_type == "רגיל" ~ "Regular",
      visit_type == "שיחה טלפונית" ~ "Virtual",
      visit_type == "שיחת וידאו" ~ "Virtual",
      TRUE ~ as.character(visit_type)
    )),
    nyha = factor(nyha, levels = c("II", "I", "III", "IV", "Undetermined")),
    weight_category = factor(case_when( #TODO - how to categorize this? should use many NAs
      weight_category %in% c("תקין", "טוב", "סביר", "תקין (מעט רזה)") ~ "Normal",
      weight_category %in% c("רזה", "רזה .", "רזון", "רזה מאד", "מעט רזה", "ירוד",
                             "רזה (ירד במשקל לאחר התקף הלב להערכתו מעל 10 ק\"ג)") ~ "Thin",
      weight_category %in% c("כחקטי", "קכקטי", "קכקסיה", "כחקטית") ~ "Cachexia",
      weight_category %in% c("עודף משקל", "עודף משקל קל", "השמנה יתירה", "עודף משקל ניכר",
                             "עודף משקל קשה", "השמנה ממארת", "השמנה", "השמנה מרכזית", "עודף משקל בינוני",
                             "השמנה ניתירה", "השמנה קלה", "עדיין השמנה יתירה", "עודף משקל .") ~ "Overweight",
      TRUE ~ as.character(weight_category)
    )),
    hf_class = factor(case_when(
      hf_class == "A" ~ "A",
      hf_class == "B" ~ "B",
      hf_class == "C" ~ "C",
      hf_class == "D" ~ "D",
    ), levels = c("A", "B", "C", "D")), 
    edema_detail = factor(edema_detail), #TODO - what to do with this?
    pillows_reason_dyspnea = grepl("BIPAP|CIPAP|CPAP|נשימה|חמצן|משתעל",pillows_reason),
    jvp_cm = case_when(
      is.na(jvp_cm) & !jvp ~ 3,
      TRUE ~ as.numeric(jvp_cm)), 
    smoke_stop_year = ifelse(smoke_stop_year < 1900, smoke_stop_year + 1900, smoke_stop_year),
    diag_cardiac = grepl(paste0(cardiac_keywords, collapse = "|"), diag, ignore.case = TRUE),
    diag_num = sapply(strsplit(diag, "\\^"), length),
    ischemic_etiology = grepl("\\bאיסכמי\\b",etiology_1, ignore.case = TRUE) & !grepl("\\bלא איסכמי\\b",etiology_1, ignore.case = TRUE),
    sglt2 = grepl("Jardiance|Empagliflozin|Forxiga|Dapagliflozin|Xigduo|Glyxambi|Empagliflozin",medications, ignore.case = TRUE),
    glp1 = grepl("Victoza|Trulicity|Ozempic|Saxenda|Xultophy|Liraglutide|Dulaglutide|Semaglutide",medications, ignore.case = TRUE),
    fusid = grepl("\\b[Ff]usid\\b|\\b[Ff]urosemide\\b",medications, ignore.case = TRUE),
    torsemide = grepl("\\bTorsemide\\b|\\bTorsemid\\b",medications, ignore.case = TRUE),
    bumetanide = grepl("\\bBumetanide\\b",medications, ignore.case = TRUE),
    metolazone = grepl("\\bMetolazone\\b|\\bZaroxylin\\b|\\bZaroxillin\\b|\\bZaroxyllyn\\b",medications, ignore.case = TRUE),
    k_spare = grepl("\\bSpironolactone\\b|\\bAldactone\\b|\\bAldospirone\\b|\\bEplerenone\\b|\\bInspra\\b",medications, ignore.case = TRUE),
    thiazide = grepl("\\bHydrochlorothiazide\\b",medications, ignore.case = TRUE),
    allopurinol = grepl("\\bAllopurinol\\b|\\bAlloril\\b|\\bZylol\\b|\\bZylori\\b",medications, ignore.case = TRUE),
    beta_blocker = grepl("Bisoprolol|Concor|Cardiloc|Metoprolol|Lopresor|Neobloc|Metopress|Carvedilol|Dimitone|Atenolol|Normiten|Normalol|Propranolol|Deralin|Prolol",medications, ignore.case = TRUE),
    anticoagulant = grepl("Xarelto|Rivaroxaban|Coumadin|Warfarin|Eliquis|Apixaban|Pradaxa|Dabigatran|Sintrom|Acenocoumarol",medications, ignore.case = TRUE),
    antiplatlet = grepl("Plavix|Clopidogrel|Brilinta|Ticagrelor|Effient|Prasugrel|Aspirin|Acetylsalicylic|Cardioaspirin|Cardiopirin|Aspocid|Aspicor|Aspimax|Aspirin",medications, ignore.case = TRUE),
    ace_inhibitor = grepl("Tritace|Enalapril|Ramipril|Lisinopril|Captopril|Perindopril|Quinapril|Trandolapril|Moexipril|Benazepril",medications, ignore.case = TRUE),
    arbs = grepl("Losartan|Ocsaar|Valsartan|Diovan|Irbesartan|Candesartan|Atacand|Olmesartan|Telmisartan|Eprosartan|Azilsartan",medications, ignore.case = TRUE),
    statin = grepl("Atorvastatin|Simvastatin|Rosuvastatin|Pravastatin|Fluvastatin|Lovastatin|Pitavastatin",medications, ignore.case = TRUE),
    insulin = grepl("Lantus|Levemir|Toujeo|Novolog|Humalog|Apidra|Fiasp|Novorapid|Humulin|Novolin|Tresiba|Basaglar|Lantus|Levemir|Toujeo|Novolog|Humalog|Apidra|Fiasp|Novorapid|Humulin|Novolin|Tresiba|Basaglar",medications, ignore.case = TRUE),
    metformin = grepl("Metformin|Glucophage|Glucovance|Janumet|Jentadueto|Kombiglyze|Synjardy|Xigduo|Invokamet|Qtern|Steglatro|Segluromet|Glyxambi|Steglujan|Glyxambi|Steglujan",medications, ignore.case = TRUE),
    dpp4 = grepl("Januvia|Onglyza|Tradjenta|Nesina|Galvus|Zomelis|Qtern|Steglatro|Segluromet|Glyxambi|Steglujan|Glyxambi|Steglujan",medications, ignore.case = TRUE),
    glinide = grepl("Starlix|Prandin|Novonorm|Starlix|Prandin|Novonorm",medications, ignore.case = TRUE),
    sulfonylurea = grepl("Glibenclamide|Glyburide|Glimepiride|Gliclazide|Tolbutamide|Diamicron|Amaryl|Daonil|Glynase|Glycomet|Glucotrol|Glyburide|Glimepiride|Gliclazide|Tolbutamide|Diamicron|Amaryl|Daonil|Glynase|Glycomet|Glucotrol",medications, ignore.case = TRUE),
    med_num = sapply(strsplit(medications, "\\^"), length),
    date_ok = x6mw_date == clin_date) %>%
  # arrange free text mwdical history and medications
  mutate(family_history_detail = family_history, 
         family_history = case_when( 
           family_history %in% c("לא ידוע", "לא", "שולל", "אין", "ללא", "לא ידוע על איס\"ל במשפחה.",
                                 "שולל מחלת לב במשפחה", "שולל מחלות לב או מוות פתאומי",
                                 "שולל מחלת לב במשפחה.", "ללא רקע משפחתי", "לא ידועיש יל\"ד משפחתי.",
                                 "לאלאב היתה שחפת ככל הנראה", "לא ידוע (ניצול שואה)", "שוללת",
                                 "ללא סיפור משפחתי", "ללא חולי לב במשפחה", "הורים נפטרו משחפת.",
                                 "שולל מחלות לב במשפחה", "ללא רקע") ~ FALSE,
           TRUE ~ !is.na(family_history)),
         family_history_detail = ifelse(family_history, family_history_detail, NA_character_),
         .after = family_history) %>%
  mutate(hosp_detail = hosp, 
         hosp = case_when(
           hosp_detail %in% c("לא", "לאלא", "לא;", "לא. עבר ניתוח של טריגר פינגר באשפוז יום",
                              "לא.עקירה כירורגית של שן", "לא רלוונטי", "לא (אלקטיבי)",
                              "לא - למעט מרפאה המטולוגית", "לא (אלקטיבי לICD)" , 
                              "לא (אשפוז יום).", "לא (היה מאושפז לפני כ8 חודשים בקרדיולוגיה ועבר השתלת ICD)",
                              "לא (למעט אלקטיבי)", "לא (למעט לביופסיה)", "לא (מטופל באשפוז יום המטולוגי)",
                              "לא (שדרוג לCRTD אלקטיבי.", "לא בארבע השנים האחרונות", "לא היה אשפוז",
                              "לא למעט אשפוז יום", "לא למעט ניתוח", "לא מעבר לאשפוז יום", 
                              "לא. היתה מאושפזת לפני שנתיים בקרדיולוגיה,  אז אבחנו לה CHF") ~ FALSE,
           TRUE ~ !is.na(hosp_detail)),
         hosp_detail = ifelse(hosp, hosp_detail, NA_character_),
         .after = hosp) %>%
  mutate(er_detail = er, 
         er = case_when( 
           er_detail %in% c("לא", "לא;", "לא רלוונטי") ~ FALSE,
           TRUE ~ !is.na(er_detail)
         ),
         er_detail = ifelse(er, er_detail, NA_character_),
         .after = er) %>%
  mutate(smoke = factor(case_when(
    smoke_cig_per_day > 10 ~ ">10",
    smoke_cig_per_day > 0 ~ "1-10",
    smoke_stop_year <= year(clin_date) ~ "past smoker",
    !is.na(smoke_age) ~ "past smoker",
    TRUE ~ "never smoked"), levels = c("never smoked", "past smoker", "1-10", ">10")),
    .before = smoke_cig_per_day) %>%
  # fix duplicated rows
  group_by(id, clin_date) %>%
  summarise(across(everything(), merge_rows), .groups = 'drop') %>%
  # fix missing height & labs
  group_by(id) %>%
  arrange(clin_date) %>%
  mutate(across(c("ht", "wt"), ~ if (n() > 1 & any(!is.na(ht))) {
    na.approx(.x, rule = 2)
  } else {
    .x
  })) %>%
  fill(all_of(c("ht", "wt")), .direction = "downup") %>%
  fill(contains("_res"), .direction = "down") %>%
  ungroup() %>%
  mutate(bmi = wt/((ht/100)^2)) %>%
  filter(date_ok,
         !is.na(x6mw_date),
         !is.na(x6mw_dist_meter),
         year(x6mw_date) >= 2014)

# create final DF ---------------------------------------------------------------
clean_df <- raw_df %>%
  left_join(mortality, by = "id") %>%
  arrange(id, x6mw_date) %>%
  group_by(id) %>%
  mutate(
    x6mw_dist_meter = x6mw_dist_meter/six_min_scale,
    age = age/age_scale,
    first_visit = first(x6mw_date),
    time = round(time_length(x6mw_date - first_visit, "months")),
    time_diff_months = time_length(fup_date - first_visit, "months"),
    last_fup_date = case_when(
      mortality_status == 1 & time_diff_months <= max_fup_time ~ fup_date,
      mortality_status == 1 & time_diff_months > max_fup_time ~ first_visit + months(max_fup_time),
      mortality_status == 0 ~ pmin(fup_date, first_visit + months(max_fup_time), na.rm = TRUE)
    ),
    mortality_status = case_when(
      mortality_status == 1 & time_diff_months > max_fup_time ~ 0,
      TRUE ~ mortality_status
    ),
    fup_time = round(time_length(last_fup_date - first_visit, "months")),
    fup_time = ifelse(fup_time == 0, 1, fup_time)) %>%
  ungroup() %>%
  filter(time <= max_fup_time) %>%
  select(id, mortality_status, fup_time, 
         time, x6mw_dist_meter, all_of(seattle_vars))

# Table 1 ---------------------------------------------------------------------
clean_df %>%
  group_by(id) %>%
  mutate(n_6mw = length(unique(x6mw_dist_meter)),
         x6mw_change = x6mw_dist_meter[which.max(time)] - x6mw_dist_meter[which.min(time)]) %>%
  ungroup() %>%
  filter(time == 0) %>%
  mutate(x6mw_dist_meter_q = cut(x6mw_dist_meter,c(
    0,quantile(x6mw_dist_meter, c(0.333, 0.667)),max(x6mw_dist_meter)),
    labels = c("Low percentile", "Middle percentile", "High percentile"))) %>%
  rename_all(function(x) sapply(x, label_get,USE.NAMES = FALSE)) %>%
  select(-c(id, time, )) %>%
  tbl_summary(by = x6mw_dist_meter_q, 
              type = list(label_get("n_6mw") ~ "continuous"),
              statistic = list(all_continuous() ~ "{mean} ({sd})"),
              missing = "no") %>%
  add_n(statistic = "{N_miss} ({p_miss})") %>%
  modify_header(n = "**Missing**") %>%
  add_p() -> tbl_1

gt::gtsave(as_gt(tbl_1), file = "export/tbl_1.html")

# save data --------------------------------------------------------------------
rm(cardiac_keywords, ex_id, mortality, raw_df)
save(clean_df, age_scale, six_min_scale, file = "data/CHF_data.RData")
