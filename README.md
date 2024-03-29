# The transmission pattern of hand, foot and mouth disease 
Please cite as: Chen, Y., Nguyet, L. A., Qui, P. T., Hong, N. T. T., Ny, N. T. H., Anh, N. T., ... & Clapham, H. (2024). Age-time-specific transmission of hand-foot-and-mouth disease enterovirus serotypes in Vietnam: A catalytic model with maternal immunity. Epidemics, 100754

## Code
casedata_preprocess.R: the code preprocessing the case data.

models.R: the code for the implementation of the models.

result_analysis.R: the code generating estimations of Force of Infection, cases and seroprevalence based on the model output.

## Data
### Original data
The original case data.

The original serological survey data.

### Other related data
The population size data of Ho Chi Minh City, Vietnam.

The preprocessed case data (code available in casedata_preprocess.R).

The stan data as the input of the stan code to implement the models (code available in models.R).

The index of the age-time-specific Force of Infection.

### Model output data
The model fitting results (output of models.R).

### Estimations
The estimated Force of Infection, cases and seroprevalence from the models, which could be used to generate the figures (output of result_analysis.R).
