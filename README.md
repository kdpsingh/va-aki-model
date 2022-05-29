A repository to accompany our draft paper entitled “Reproducibility and
generalizability of a state-of-the-art acute kidney injury prediction
model”

# System Information

Our data preparation and analysis code were run on both VA and
University of Michigan Virtual Windows Machines, with different versions
of R (see details in paper).

The code we provide here as demo should work for:

Operating systems: Windows, Mac, Linux

Software: R 3.6+

The models we provide in this repository are MOJO objects which are
compatible with all versions of h2o.

Runtime: The runtime for preparing the data using the `gpmodels` depends on local computing resources and sample size, varying from hours to days. The runtime for generating predictions (Step 6) depends on local computing resources and test set size, varying from seconds to minutes.

# Instruction

## Step 1 - Prepare fixed and time-varying/temporal data

After gathering patients’ data, prepare one dataset for the fixed
variables and another for the time-varying/temporal variables.

The fixed dataset has one row for each hospital encounter, storing
baseline information of the encounter. The temporal dataset has multiple
rows for each hospital encounter, storing time-varying information of
the encounter with time stamps.

Please refer to our package
[gpmodels](https://github.com/ML4LHS/gpmodels) for sample fixed data and
sample temporal data.  
*Note: variables in sample\_fixed\_data and sample\_temporal\_data do
not reflect the full list of variables used in this study.*

Please refer to Supplementary Table 1 for the list of fixed and temporal
variables used in this study.

## Step 2 - Transform data

    # Install gpmodels from our github repository
    remotes::install_github('ML4LHS/gpmodels')
    library(gpmodels)
    library(tidyverse)
    library(lubridate)
    library(data.table)

    input = "[path_to_fixed_and_temporal_datasets]"
    output = "[data_output_path]"
    cohort = "2016to2020"

    ##################################################################
    ##                          Input Data                          ##
    ##################################################################

    fixed_data = 
      fread(file.path(input, 
                      paste("Inpatient", cohort, "fixed_data.csv", sep = "_"))) %>% 
      janitor::clean_names() %>% 
      mutate(admit_date = ymd_hms(admit_date), 
             discharge_date = ymd_hms(discharge_date)) %>% 
      rename(admit_time = admit_date, 
             discharge_time = discharge_date)

    temporal_data = fread(file.path(input, 
                                    paste("Inpatient", cohort, "temporal_data.csv", sep = "_"))) %>% 
      janitor::clean_names() %>% 
      mutate(time = ymd_hms(time))

    #################################################################
    ##                     Define time frame                      ##
    #################################################################

    tf_aki = 
      time_frame(
        fixed_data = fixed_data %>% select(-patient_id),
        temporal_data = temporal_data %>% filter(encounter_id %in% fixed_data$encounter_id) %>% select(-patient_id),
        fixed_id = "encounter_id", 
        fixed_start = "admit_time", 
        fixed_end = "discharge_time", 
        temporal_id = "encounter_id", 
        temporal_time = "time", 
        temporal_variable = "variable", 
        temporal_category = "category", 
        temporal_value = "value", 
        step = hours(6), 
        max_length = days(7), 
        chunk_size = 100, 
        output_folder = output, 
        create_folder = TRUE)

    tf_aki = 
      tf_aki %>% 
      pre_dummy_code()


    ##################################################################
    ##                        Add Predictors                        ##
    ##################################################################

    # Enable paralell computing
    future::plan('multisession')

    tf_aki %>% 
      add_rolling_predictors(category = c('lab', 'vitals'), 
                             lookback = hours(48), 
                             window = hours(6), 
                             stats = c(length = length,
                                       min = min, 
                                       mean = mean, 
                                       median = median, 
                                       max = max)) %>% 
      add_rolling_predictors(category = "medications", 
                             lookback = days(7), 
                             window = hours(24), 
                             stats = c(length = length))

    ##################################################################
    ##                         Add Outcomes                         ##
    ##################################################################

    tf_aki %>% 
      add_rolling_outcomes(variables = "sCr", 
                           lookahead = hours(48), 
                           stats = c(max = max))

    #################################################################
    ##                    Combine to Model Data                    ##
    #################################################################

    tf_combined = 
      tf_aki %>% 
      combine_output()

    tf_combined %>% write_csv(file.path(output, "Combined_dataset.csv"))

## Step 3 - Define outcome AKI stages and add calculated predictors

    library(tidyverse)

    aki_data = read_csv(file.path(output, "Combined_dataset.csv"))

    # Clean up dialysis outcome
    aki_data = 
      aki_data %>% 
      mutate(dialysis_time = lubridate::time_length(dialy_1stdatetime - admit_time, unit = "hours")) %>% 
      mutate(outcome_dialysis = if_else(!is.na(dialysis_time), 1, 0)) %>% 
      select(-dialy_1stdatetime)

    # Split outcome AKI into stages
    aki_data = 
      aki_data %>% 
      mutate(outcome_aki_stage = 
               case_when(outcome_dialysis == 1 & time >= (dialysis_time - 48) ~ "aki_3d", 
                         outcome_sCr_max_48 >= baseline_scr * 3 ~ "aki_3",
                         outcome_sCr_max_48 > baseline_scr & outcome_sCr_max_48 >= 4 ~ "aki_3",
                         outcome_sCr_max_48 >= baseline_scr * 2 & outcome_sCr_max_48 < baseline_scr * 3 ~ "aki_2",
                         outcome_sCr_max_48 >= baseline_scr * 1.5 & outcome_sCr_max_48 < baseline_scr * 2 ~ "aki_1",
                         outcome_sCr_max_48 >= baseline_scr + 0.3 ~ "aki_1",
                         is.na(outcome_sCr_max_48) ~ NA_character_,
                         TRUE ~ "no_aki")) %>% 
      mutate(outcome_aki_stage = outcome_aki_stage %>% as.factor())
    aki_data %>% count(outcome_aki_stage)

    # Current AKI stage
    aki_data = 
      aki_data %>% 
      mutate(current_aki_stage = 
               case_when(time >= dialysis_time ~ 'aki_3d', 
                         sCr_max_06 >= baseline_scr * 3 ~ 'aki_3',
                         sCr_max_06 > baseline_scr & sCr_max_06 >= 4 ~ 'aki_3',
                         sCr_max_06 >= baseline_scr * 2 & sCr_max_06 < baseline_scr * 3 ~ 'aki_2',
                         sCr_max_06 >= baseline_scr * 1.5 & sCr_max_06 < baseline_scr * 2 ~ 'aki_1',
                         sCr_max_06 >= baseline_scr + 0.3 ~ 'aki_1',
                         is.na(sCr_max_06) ~ NA_character_,
                         TRUE ~ 'no_aki')) %>% 
      mutate(current_aki_stage = current_aki_stage %>% as.factor())
    aki_data %>% count(current_aki_stage)

    # Add creatinine ratio and difference predictors
    aki_data = 
      aki_data %>% 
      mutate(cr_ratio_to_baseline_06 = sCr_max_06 / baseline_scr, 
             cr_diff_to_baseline_06 = sCr_max_06 - baseline_scr, 
             bun_to_cr_ratio_06 = BUN_max_06 / sCr_max_06)

    # Add CVP predictors - UM only step since UM does not collect CVP information
    cvp_predictors = 
      grep('^T_\\w+_\\d+$', names(aki_data), value = TRUE) %>% 
      str_replace("^T", "CVP")

    cvp_frame = matrix(ncol = length(cvp_predictors))
    colnames(cvp_frame) = cvp_predictors
    cvp_frame = cvp_frame %>% as_tibble()

    aki_data = 
      aki_data %>% 
      bind_cols(cvp_frame)


    # Inclusion/Exclusion
    cases_to_include = 
      aki_data %>% 
      select(encounter_id, time, baseline_scr, current_aki_stage, outcome_aki_stage)

    encounter_ids_to_include =
      cases_to_include %>% 
      filter(baseline_scr < 4) %>% 
      # filter(encounter_id %in% encounter_ids_without_aki_on_admission) %>% 
      # filter(!is.na(outcome_aki_stage)) %>% 
      pull(encounter_id) %>% 
      unique() %>% 
      as.character()

    aki_data = 
      aki_data %>% 
      filter(encounter_id %in% encounter_ids_to_include)

    # Propagate both current AKI stage and outcome AKI stage forward
    aki_stage_to_propagate = 
      aki_data %>% 
      select(encounter_id, time, current_aki_stage, outcome_aki_stage)

    aki_stage_to_propagate = 
      aki_stage_to_propagate %>% 
      group_by(encounter_id) %>% 
      fill(current_aki_stage, .direction = 'down') %>% 
      ungroup() %>% 
      mutate_at(vars(current_aki_stage), 
                . %>%
                  {
                    if_else(is.na(.),
                            factor('no_aki', levels=levels(aki_stage_to_propagate$current_aki_stage)),
                            .)
                  }) %>% 
      mutate(outcome_aki_stage = coalesce(outcome_aki_stage, current_aki_stage))

    aki_data$current_aki_stage = aki_stage_to_propagate$current_aki_stage
    aki_data$outcome_aki_stage = aki_stage_to_propagate$outcome_aki_stage

    # Remove remaining rows with missing outcomes - NONE
    aki_data = 
      aki_data %>% 
      filter(!is.na(outcome_aki_stage))

    aki_data = write_csv(file.path(output, "AKI_all.csv"))

A [sample
data](https://github.com/ML4LHS/va-aki-model/blob/main/sample_data/sample_data.csv)
resulting from Step 3 is provided in this repository.

## Step 4 - Split data into training/validation/test

VA - Training:Validation:Test = 64%:16%:20%

UM - Training:Validation:Test = 20%:20%:60%

Here we show the code to split UM data.

    library(h2o)
    h2o.init()

    encounters = read_csv(file.path("[path_to_raw_files]" ,"Inpatient_2016to2020_Encounters.csv"))
    h2o_aki_all = h2o.importFile(file.path(output, "AKI_all.csv"),  
                                 destination_frame = "h2o_aki_all")


    enc_ids = h2o_aki_all$encounter_id %>% as_tibble()

    # Add patient IDs
    enc_ids = 
      enc_ids %>% 
      left_join(encounters %>% select(EncounterID, PatientID),
                by = c("encounter_id" = "EncounterID"))

    pt_ids = enc_ids %>% distinct(PatientID) %>% pull()

    # Split into training/validatio/test - sample at patient level
    set.seed(1)
    pt_training = sample(pt_ids, 0.2*length(pt_ids))

    set.seed(1)
    pt_validation = sample(setdiff(pt_ids, pt_training), 0.2*length(pt_ids))

    pt_test = setdiff(pt_ids, c(pt_training, pt_validation))

    enc_ids = 
      enc_ids %>% 
      mutate(training = if_else(PatientID %in% pt_training, 1, 0), 
             validation = if_else(PatientID %in% pt_validation, 1, 0), 
             test = if_else(PatientID %in% pt_test, 1, 0))

    h2o_aki_training = h2o_aki_all[which(enc_ids$training == 1), ]
    h2o_aki_validation = h2o_aki_all[which(enc_ids$validation == 1), ]
    h2o_aki_test = h2o_aki_all[which(enc_ids$test == 1), ]

    h2o.exportFile(h2o_aki_training, file.path(output, "AKI_training_20_percent.csv"))
    h2o.exportFile(h2o_aki_validation, file.path(output, "AKI_validation_20_percent.csv"))
    h2o.exportFile(h2o_aki_test, file.path(output, "AKI_test_60_percent.csv"))

    h2o.shutdown()

## Step 5 - Original VA Model and Extended VA Model Training

Here is our code to train the original VA model using VA data.

    library(tidyverse)
    library(h2o)

    memory.limit(1e5)
    h2o.init(port = 1161, nthreads = 12, max_mem_size = '80G')

    h2o_aki_data = h2o.importFile("[va_data_input_path]", destination_frame = 'h2o_aki_data')
    va_encounters = read_rds('va_encounters.rds')
    va_training = read_rds('va_training_encounters.rds')
    va_validation = read_rds('va_validation_encounters.rds')
    va_test = read_rds('va_test_encounters.rds')

    # Split VA data
    h2o_training_rows = h2o_aki_data[, 'encounter_id'] %in% as.character(va_training)
    h2o_training_data = h2o_aki_data[h2o_training_rows, ]

    h2o_validation_rows = h2o_aki_data[, 'encounter_id'] %in% as.character(va_validation)
    h2o_validation_data = h2o_aki_data[h2o_validation_rows, ]

    h2o_test_rows = h2o_aki_data[, 'encounter_id'] %in% as.character(va_test)
    h2o_test_data = h2o_aki_data[h2o_test_rows, ]

    # Define predictors and outcome
    predictors = c([list of predictors])
    outcome = 'outcome_aki_stage'

    h2o_gbm = h2o.gbm(x = predictors,
                      y = outcome,
                      model_id = 'h2o_gbm_original_model',
                      training_frame = h2o_training_data,
                      validation_frame = h2o_validation_data,
                      seed = 1,
                      ntrees = 1000,
                      score_tree_interval = 10,
                      stopping_rounds = 5,
                      stopping_metric = 'logloss',
                      stopping_tolerance = 0.0005,
                      categorical_encoding = 'SortByResponse',
                      verbose = TRUE,
                      export_checkpoints_dir = "[model_output_path]")
    h2o.save_mojo(h2o_gbm, "[model_output_path]")

    h2o.shutdown()

Here is our code to train the extended VA model at UM. The training
strategy is the same as the original model except that we update the
early stopping criteria based on a lack of log loss improvement of
*0.0001* after 5 consecutive rounds.

    library(tidyverse)
    library(h2o)

    h2o.init(max_mem_size = "40G", nthreads = 6)

    h2o_aki_training = h2o.importFile(file.path("[um_data_input_path]", "AKI_training_20_percent.csv"))
    h2o_aki_validation = h2o.importFile(file.path("[um_data_input_path]", "AKI_validation_20_percent.csv"))
    va_aki_model = h2o.import_mojo(file.path("[model_input_path]", "h2o_gbm_original_model.zip"))

    predictors = va_aki_model@parameters$x
    outcome = va_aki_model@parameters$y

    h2o_gbm_mm = h2o.gbm(model_id = "h2o_gbm_extended_model",
                         x = predictors,
                         y = outcome,
                         training_frame = h2o_aki_training,
                         validation_frame = h2o_aki_validation,
                         checkpoint = va_aki_model@model_id,
                         seed = 1,
                         ntrees = 1000,
                         score_tree_interval = 10,
                         stopping_rounds = 5,
                         stopping_metric = "logloss",
                         stopping_tolerance = 0.0001,
                         categorical_encoding = "SortByResponse",
                         verbose = TRUE,
                         export_checkpoints_dir = "[model_output_path]")
    h2o.save_mojo(h2o_gbm_mm, "[model_output_path]")

    h2o.shutdown()

## Step 6 - Test data prediction

Here we show how to use the
[original](https://github.com/ML4LHS/va-aki-model/blob/main/models/h2o_gbm_original_model.zip)
and
[extended](https://github.com/ML4LHS/va-aki-model/blob/main/models/h2o_gbm_extended_model.zip)
model to make prediction on the UM test dataset.

    library(tidyverse)
    library(h2o)

    h2o.init(max_mem_size = '40G', nthreads = 6)

    # Load test data
    h2o_aki_test = h2o.importFile(file.path("[data_input_path]", "AKI_test_60_percent.csv"))

    # For demonstration purpose, use provided sample data
    #h2o_aki_test = h2o.importFile(file.path("[data_input_path]", "sample_data.csv"))

    # Load model
    va_model = h2o.import_mojo(file.path("[model_input_path]", "h2o_gbm_original_model.zip"))
    mm_model = h2o.import_mojo(file.path("[model_input_path]", "h2o_gbm_extended_model.zip"))

    # Use original VA model to make prediction on the test set
    va_model_test_preds = h2o.predict(va_model, h2o_aki_test) %>% as_tibble()

    # Use extended model to make prediction on the test set
    mm_Model_test_preds = h2o.predict(mm_model, h2o_aki_test) %>% as_tibble()

    h2o.shutdown()

We provide [sample predictions by original
model](https://github.com/ML4LHS/va-aki-model/blob/main/sample_predictions/sample_predictions_original_model.csv)
and [sample predictions by extended
model](https://github.com/ML4LHS/va-aki-model/blob/main/sample_predictions/sample_predictions_extended_model.csv)
in this repository.

*Note:*

*1. In our study, an unseen test dataset is used to generate the
predictions for evaluation. Here, for demonstration purpose, we use the
sample\_data.csv provided in this repository to generate the
predictions.*

*2. The “predict” column in the sample predictions are auto-generated by
h2o.predict(). In our study, we do not use this in our model evaluation.
For more details, see Methods - Model Evaluation in our paper.*

*3. “aki\_1”, “aki\_2”, “aki\_3”, “no\_aki” columns represent the
predicted risk values for progressing to each AKI stage in the next 48
hours by the model. For details, see Methods - Outcome Definition. *

Predictions by the original model:

<table style="width:100%;">
<colgroup>
<col style="width: 10%" />
<col style="width: 4%" />
<col style="width: 14%" />
<col style="width: 17%" />
<col style="width: 17%" />
<col style="width: 17%" />
<col style="width: 17%" />
</colgroup>
<thead>
<tr class="header">
<th>encounter_id</th>
<th>time</th>
<th>outcome_aki_stage</th>
<th>aki_1</th>
<th>aki_2</th>
<th>aki_3</th>
<th>no_aki</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>1</td>
<td>144</td>
<td>no_aki</td>
<td>0.0064551884955890235</td>
<td>5.222858922556572e-4</td>
<td>8.459502707443686e-4</td>
<td>1.412122574153119e-4</td>
</tr>
<tr class="even">
<td>2</td>
<td>48</td>
<td>no_aki</td>
<td>0.03600354069303235</td>
<td>0.0021072801356229573</td>
<td>0.0021641676522983254</td>
<td>1.1449415922994847e-4</td>
</tr>
<tr class="odd">
<td>3</td>
<td>60</td>
<td>no_aki</td>
<td>0.03384526925687363</td>
<td>5.596368040078928e-4</td>
<td>4.2588067841007406e-4</td>
<td>1.972534921280478e-4</td>
</tr>
<tr class="even">
<td>4</td>
<td>0</td>
<td>no_aki</td>
<td>0.06748869291831279</td>
<td>0.015829967056608125</td>
<td>0.004964146269355001</td>
<td>4.938794831467228e-4</td>
</tr>
<tr class="odd">
<td>4</td>
<td>6</td>
<td>aki_1</td>
<td>0.07718671685477586</td>
<td>0.011760862880380737</td>
<td>0.004462500716397853</td>
<td>3.671117370507371e-4</td>
</tr>
<tr class="even">
<td>5</td>
<td>54</td>
<td>aki_1</td>
<td>0.9972816792865564</td>
<td>3.0218589306194726e-4</td>
<td>8.564660049311132e-5</td>
<td>6.693205634094189e-4</td>
</tr>
<tr class="odd">
<td>5</td>
<td>60</td>
<td>aki_2</td>
<td>0.9957865319680796</td>
<td>4.88483694641826e-4</td>
<td>1.1565439792487673e-4</td>
<td>9.359069493640218e-4</td>
</tr>
<tr class="even">
<td>6</td>
<td>18</td>
<td>aki_2</td>
<td>0.9977068222930573</td>
<td>4.115113434234649e-4</td>
<td>3.839312808716371e-4</td>
<td>2.9839747478945265e-4</td>
</tr>
<tr class="odd">
<td>6</td>
<td>24</td>
<td>aki_3</td>
<td>0.9976837215491442</td>
<td>4.1573650617112747e-4</td>
<td>5.694113018531597e-4</td>
<td>4.1966539190284086e-4</td>
</tr>
<tr class="even">
<td>6</td>
<td>30</td>
<td>aki_3</td>
<td>0.9975796740956613</td>
<td>3.7793918425647877e-4</td>
<td>4.886861230560208e-4</td>
<td>3.6819038928775515e-4</td>
</tr>
</tbody>
</table>

Predictions by the extended model:

<table>
<colgroup>
<col style="width: 10%" />
<col style="width: 4%" />
<col style="width: 14%" />
<col style="width: 16%" />
<col style="width: 17%" />
<col style="width: 17%" />
<col style="width: 17%" />
</colgroup>
<thead>
<tr class="header">
<th>encounter_id</th>
<th>time</th>
<th>outcome_aki_stage</th>
<th>aki_1</th>
<th>aki_2</th>
<th>aki_3</th>
<th>no_aki</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>1</td>
<td>144</td>
<td>no_aki</td>
<td>0.006234659354318668</td>
<td>3.9396334357326586e-4</td>
<td>6.082752765154707e-4</td>
<td>8.65481823656329e-5</td>
</tr>
<tr class="even">
<td>2</td>
<td>48</td>
<td>no_aki</td>
<td>0.03181989816474396</td>
<td>0.0015259105059579803</td>
<td>0.0014938448792595441</td>
<td>6.736404516307144e-5</td>
</tr>
<tr class="odd">
<td>3</td>
<td>60</td>
<td>no_aki</td>
<td>0.03863906250008758</td>
<td>4.370226038962381e-4</td>
<td>3.1702490513051335e-4</td>
<td>1.25158505707926e-4</td>
</tr>
<tr class="even">
<td>4</td>
<td>0</td>
<td>no_aki</td>
<td>0.08524148096402387</td>
<td>0.014462842054231324</td>
<td>0.004323405924458679</td>
<td>3.6663392947475703e-4</td>
</tr>
<tr class="odd">
<td>4</td>
<td>6</td>
<td>aki_1</td>
<td>0.09565850641124808</td>
<td>0.010259884130698331</td>
<td>0.003710986546062025</td>
<td>2.6021934264484306e-4</td>
</tr>
<tr class="even">
<td>5</td>
<td>54</td>
<td>aki_1</td>
<td>0.9507510793806109</td>
<td>0.031460609108969334</td>
<td>7.493518809333979e-4</td>
<td>0.00151335194027806</td>
</tr>
<tr class="odd">
<td>5</td>
<td>60</td>
<td>aki_2</td>
<td>0.9486004134775309</td>
<td>0.020612851020189685</td>
<td>0.0010111274471720974</td>
<td>0.00211449429225202</td>
</tr>
<tr class="even">
<td>6</td>
<td>18</td>
<td>aki_2</td>
<td>1</td>
<td>0</td>
<td>0</td>
<td>0</td>
</tr>
<tr class="odd">
<td>6</td>
<td>24</td>
<td>aki_3</td>
<td>0.9106039396228075</td>
<td>0.0522658708963932</td>
<td>0.00591804537200534</td>
<td>0.0020818909115962406</td>
</tr>
<tr class="even">
<td>6</td>
<td>30</td>
<td>aki_3</td>
<td>0.9074636191598131</td>
<td>0.04735512019593059</td>
<td>0.005598886085908977</td>
<td>0.0018204227476364945</td>
</tr>
</tbody>
</table>
