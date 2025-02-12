from functions import *
import os

def main():
    # get currrent working directory
    cwd = os.getcwd()

    # define file paht
    clinical_file = os.path.join(cwd, 'PDAC_Patient_clinical_data.xlsx')
    bca_file = os.path.join(cwd, 'PDAC_Patient_BCA_data.xlsx')


    df_clinical, sheets_dict = read_files(clinical_file, bca_file, header=1)

    # Process and merge data.
    merged_df = process_and_merge(df_clinical, sheets_dict)

    # save merged data.
    output_file = os.path.join(cwd, 'PDAC_Patient_clinical_BCA.xlsx')
    save_file(merged_df, output_file)

    df_imputed = impute_data(merged_df)

    imputation_check(merged_df, df_imputed, 'abdominal_cavity_pericardial_adipose_tiss_(mean)' )

    variables_to_test = [
            "whole_scan_bone_density_(mean_hu)",
            "abdominal_cavity_bone_density_(mean_hu)",
            "whole_scan_muscle_volume_(mean)",
            "abdominal_cavity_muscle_volume_(mean)",
            "whole_scan_total_adipose_tissue_vol_(mean)",
            "abdominal_cavity_total_adipose_tissue_vol_(mean)",
            "whole_scan_intramusc._adipose_tissu_(mean)",
            "abdominal_cavity_intramusc._adipose_tissu_(mean)",
            "whole_scan_subcutaneous_adip._tissu_(mean)",
            "abdominal_cavity_subcutaneous_adip._tissu_(mean)",
            "whole_scan_visceral_adipose_tissue__(mean)",
            "abdominal_cavity_visceral_adipose_tissue__(mean)",
            "whole_scan_pericardial_adipose_tiss_(mean)",
            "abdominal_cavity_pericardial_adipose_tiss_(mean)",
            'tumor_size_(mm)_'
            
    ]

    df_normality_results = run_normality_tests(df_imputed, variables_to_test)
    output_file = os.path.join(cwd, 'normality_test.xlsx')
    save_file(df_normality_results, output_file)

    for var in variables_to_test:
            if var in merged_df.columns:
                plot_qq(merged_df[var], var)
            else:
                logger.warning(f"Variable '{var}' not found in the DataFrame.")


    time_col = 'overall_survival_in_months_from_diagnosis'
    event_col = 'death_during_follow-up'
    therapy_col = 'neoadjuvant_therapy_(yes/no)'
    therapy_label = 'Neoadjuvant Therapy'


    # Preprocess the data to ensure proper types.
    df_survival = preprocess_survival_data(df_imputed, time_col, event_col, therapy_col)
    results = plot_km_curve(df_survival, time_col, event_col, therapy_col, therapy_label, save_path=Path("KM_neoadjuvant.png"))
    print("NeoAdjuvant Log-rank test p-value:", results.p_value)

    therapy_col = 'adjuvant_therapy_(yes/no)__'
    therapy_label = 'adjuvant Therapy'

    df_survival = preprocess_survival_data(df_imputed, time_col, event_col, therapy_col)
    results = plot_km_curve(df_survival, time_col, event_col, therapy_col, therapy_label, save_path=Path("KM_adjuvant.png"))
    print("Adjuvant Log-rank test p-value:", results.p_value)

    df_imputed = df_imputed.rename(columns={'sex_(male/female)': 'sex'})

    time_col = 'overall_survival_in_months_from_diagnosis'
    event_col = 'death_during_follow-up'
    therapy_cols = ['neoadjuvant_therapy_(yes/no)', 'adjuvant_therapy_(yes/no)__']

    additional_covariates = [
        'sex',
        'complications',
        'venous_invasion_yes/no_',
        't_(t-stage_pt_in_anderer_tabelle)',
        'grading',
        
    ]

    # Preprcess data
    df_cox = preprocess_cox_data(df_imputed, time_col, event_col, therapy_cols, additional_covariates)

    # fit cox model.
    cph_model = fit_cox_model(df_cox, time_col, event_col)

    body_composition_vars = [
        'abdominal_cavity_muscle_volume_(mean)',
        'abdominal_cavity_total_adipose_tissue_vol_(mean)',
        'abdominal_cavity_intramusc._adipose_tissu_(mean)',
        'abdominal_cavity_subcutaneous_adip._tissu_(mean)',
        'abdominal_cavity_visceral_adipose_tissue__(mean)',
        'whole_scan_muscle_volume_(mean)',
        'whole_scan_total_adipose_tissue_vol_(mean)',
        'whole_scan_intramusc._adipose_tissu_(mean)',
        'whole_scan_subcutaneous_adip._tissu_(mean)',
        'whole_scan_visceral_adipose_tissue__(mean)'
    ]

    # Therapy variables to test.
    therapy_vars = [ 'adjuvant_therapy_(yes/no)__', 'neoadjuvant_therapy_(yes/no)',]

    for therapy in therapy_vars:
        logger.info(f"Comparing body composition by therapy variable: {therapy}")
        compare_body_composition_by_therapy(df_imputed, therapy, body_composition_vars)


    ####


    intramuscular_col = 'whole_scan_intramusc._adipose_tissu_(mean)'
    muscle_col = 'whole_scan_muscle_volume_(mean)'
    complication_col = 'complications'  # assuming this is binary or coded as yes/no
    gender_col = 'sex'

    # convert complication column to binary 
    if not pd.api.types.is_numeric_dtype(df_imputed[complication_col]):
        df_imputed[complication_col] = df_imputed[complication_col].apply(lambda x: 1 if str(x).strip().lower() == "yes" else 0)

    # Cclculate the ratio
    df_imputed = calculate_ratio(df_imputed, intramuscular_col, muscle_col)

    # unadjusted Roc analysis to find  cutoff
    roc_analysis_for_ratio(df_imputed, 'im_to_muscle_ratio', complication_col)

    #  logistic regression adjusting for gender
    logistic_regression_adjusted(df_imputed, complication_col, 'im_to_muscle_ratio', gender_col)




    ###

    time_col = "overall_survival_in_months_from_diagnosis"
    event_col = "death_during_follow-up"

    # Define candidate body composition predictors.
    predictors = [
        'abdominal_cavity_muscle_volume_(mean)',
        'abdominal_cavity_total_adipose_tissue_vol_(mean)',
        'abdominal_cavity_intramusc._adipose_tissu_(mean)',
        'abdominal_cavity_subcutaneous_adip._tissu_(mean)',
        'abdominal_cavity_visceral_adipose_tissue__(mean)',
        'whole_scan_muscle_volume_(mean)',
        'whole_scan_total_adipose_tissue_vol_(mean)',
        'whole_scan_intramusc._adipose_tissu_(mean)',
        'whole_scan_subcutaneous_adip._tissu_(mean)',
        'whole_scan_visceral_adipose_tissue__(mean)'
    ]

    # Run univariable Cox regression analyses.
    results_df = analyze_body_composition(df_imputed, time_col, event_col, predictors)
    logger.info("Univariable Cox regression analysis complete.")
    print(results_df)

    # Save the summary table to Excel.
    output_file = Path("body_composition_survival_analysis.xlsx")
    results_df.to_excel(output_file, index=False)


if __name__ == '__main__':
    main()