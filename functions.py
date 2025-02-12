import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import logging
from sklearn.experimental import enable_iterative_imputer  # noqa
from sklearn.impute import IterativeImputer, SimpleImputer
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from pathlib import Path
from scipy.stats import shapiro
import statsmodels.api as sm
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines import CoxPHFitter
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_curve, auc
import statsmodels.api as sm
import statsmodels.formula.api as smf
import numpy as np

# Set up logging for tracking and debugging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


cwd = os.getcwd()



def read_files(clinical_file: str, bca_file: str, header: int = 1):
    """
    Reads the clinical data file and the BCA Excel file (all sheets) into DataFrames.
    
    Parameters:
        clinical_file (str): File path to the clinical data Excel file.
        bca_file (str): File path to the BCA Excel file.
        header (int): The header row index for the BCA sheets (default is 1).
    
    Returns:
        tuple: (df_clinical, sheets_dict)
            - df_clinical: DataFrame with clinical data.
            - sheets_dict: Dictionary of DataFrames from each sheet in the BCA file.
    """
    df_clinical = pd.read_excel(clinical_file)
    sheets_dict = pd.read_excel(bca_file, sheet_name=None, header=header)
    return df_clinical, sheets_dict

def process_and_merge(df_clinical: pd.DataFrame, sheets_dict: dict) -> pd.DataFrame:
    """
    Processes the BCA sheets by finding common IDs, filtering, renaming columns, merging,
    standardizing column names, and finally merging with the clinical data.
    
    parameters:
        df_clinical (pd.DataFrame): Clinical data DataFrame.
        sheets_dict (dict): Dictionary of DataFrames from the BCA file.
    
    returns:
        pd.DataFrame: The final merged DataFrame.
    """
    # identify common IDs across all BCA sheets.
    id_sets = {sheet: set(df['ID']) for sheet, df in sheets_dict.items()}
    common_ids = set.intersection(*id_sets.values())


    # filter each sheet to include only rows with common IDs.
    sheets_dict = {sheet: df[df['ID'].isin(common_ids)].copy() for sheet, df in sheets_dict.items()}
    
    # add a suffix to each column 
    for sheet, df in sheets_dict.items():
        df.rename(columns=lambda x: f"{x}_{sheet}" if x != 'ID' else x, inplace=True)
    
    # merge all BCA DataFrames on 'ID' using inner joins.
    merged_df = None
    for df in sheets_dict.values():
        if merged_df is None:
            merged_df = df
        else:
            merged_df = pd.merge(merged_df, df, on='ID', how='inner')
    
    # replace spaces with underscores and convert to lowercase.
    merged_df.columns = merged_df.columns.str.replace(' ', '_').str.lower()
    df_clinical.columns = df_clinical.columns.str.replace(' ', '_').str.lower()
    
    # Merge the BCA data with the clinical data on 'id'.
    final_df = pd.merge(merged_df, df_clinical, on='id', how='inner')
    return final_df

def save_file(df: pd.DataFrame, output_file_path: str):
    """
  
    Parameters:
        df (pd.DataFrame): The DataFrame to save.
        output_file_path (str): The path where the file will be saved.
    """
    try:
        df.to_excel(output_file_path, index=False)
        print(f"File saved successfully: {output_file_path}")
    except Exception as e:
        print(f"Error saving file: {e}")


def impute_data(df: pd.DataFrame) -> pd.DataFrame:
    """
    Impute missing values using IterativeImputer (MICE-like) for continuous variables.
    """
    def drop_columns(df: pd.DataFrame) -> pd.DataFrame:
        """
        drop columns with 100% missing values.
        """
        COLUMNS_DROP = [
        "date_of_birth", "pat_deathday", "end_of_follow-up_",
        "ct-scan_(intern/extern)", "layer_thickness_(mm)", "non-contrast_phase_(yes/no)",
        "arterial_resection_(yes/no)", "arterial_reconstruction_(yes/no)",
        "progression-free_survival_in_months_from_diagnosis",
        "progression_during_follow-up_(yes,no)", "site_of_progression"
        ]
        drop_cols = [col for col in COLUMNS_DROP if col in df.columns]
        if drop_cols:
            df = df.drop(columns=drop_cols)
        return df
    #drop columns with no data
    df = drop_columns(df)
    
  
    continuous_vars = [
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
        "tumor_size_(mm)_",
        "number_of_locoregional_lymph_nodes_",
        "ca19-9_before_surgery",
        "cea_before_surgery"
    ]
    

    # make sure cols exist in the DataFrame
    continuous_vars = [col for col in continuous_vars if col in df.columns]
    
    # build pipelines for different variable types
    continuous_pipeline = Pipeline(steps=[
        ('iter_imputer', IterativeImputer(random_state=42))
    ])
    
  
    
    # Combine pipelines with ColumnTransformer
    preprocessor = ColumnTransformer(
        transformers=[
            ('cont', continuous_pipeline, continuous_vars),
        ],
        remainder='passthrough'  # keep  remaining columns
    )
    
    logger.info("Starting imputation")
    imputed_array = preprocessor.fit_transform(df)
    
   
    remainder_cols = [col for col in df.columns if col not in continuous_vars]
    all_cols = continuous_vars + remainder_cols
    df_imputed = pd.DataFrame(imputed_array, columns=all_cols)
    
    logger.info("imputation completed")
    return df_imputed


def imputation_check(df_original: pd.DataFrame, df_imputed: pd.DataFrame, variable: str) :
    """
    sensitivy check comparing distribution of variables and ditribution of imputed values
    """
    
    # Extract observed values from the original data (non-missing) 
    observed = df_original[variable].dropna()

    # Extract the full imputed data for the same variable
    imputed = df_imputed[variable]

    # Plotting the distributions
    plt.figure(figsize=(10, 6))
    sns.kdeplot(observed, label='Observed', shade=True)
    sns.kdeplot(imputed, label='Imputed', shade=True)
    plt.title(f'Distribution Comparison for {variable}')
    plt.xlabel(variable)
    plt.legend()
    

def assess_normality(data: pd.Series, variable_name: str) -> dict:
    """
    asess the normal distribution of a variable using the Shapiroo-Wilk test.
    Parameters:
    -----------
    data : pd.Series
        seriees containing the data for the variable.
    variable_name : str
        (for logging and reporting
    
    Returns:
    --------
    dict
        dictionarry with the variable name, test statistic, and p-value.
        
    """
    # test requires complete cases.
    data_clean = data.dropna()
    
    try:
        statistic, p_value = shapiro(data_clean)
        return {"variable": variable_name, "statistic": statistic, "p_value": p_value}
    except Exception as e:
        logger.error(f"Error performing Shapiro-Wilk test for {variable_name}: {e}")
        return {"variable": variable_name, "statistic": None, "p_value": None, "error": str(e)}

def run_normality_tests(df: pd.DataFrame, variables: list) -> pd.DataFrame:
    """
    helper func. to run normality tests
    
    Parameters:
    -----------
    df : pd.DataFrame
        The DataFrame containing the variables.
    variables : list
        a lst of column names to test 
    
    Returns:
    --------
    pd.DataFrame
        A DataFrame summarizing test resultss
    """
    results = []
    for var in variables:
        if var in df.columns:
            result = assess_normality(df[var], var)
            results.append(result)
        else:
            logger.warning(f"Variable '{var}' not found in the DataFrame.")
    return pd.DataFrame(results)

def plot_qq(data: pd.Series, variable_name: str) -> None:
    """
    genrate a Q-Q plot f

    Parameters:
    -----------
    data : pd.Series
        data for which to generate theplot.
    variable_name : str
        (used for title and file naming).
  
    """
    # test requires complete dataa
    data_clean = data.dropna()
    

    fig = plt.figure()
    ax = fig.add_subplot(111)
    sm.qqplot(data_clean, line='s', ax=ax)
    ax.set_title(f"Q-Q Plot for {variable_name}")
    ax.set_xlabel("Theoretical Quantiles")
    ax.set_ylabel("Sample Quantiles")
    plt.tight_layout()
    
    save_path = os.path.join(cwd, 'Q_Q_plot.png')
    plt.savefig(save_path)
def create_overview_table(df: pd.DataFrame) -> pd.DataFrame:
    """
    creat an overview table sumarising the distribution of all variables in the DataFrame.
    for numeric variables, summary statistics include count, mean, std, min, quartiles, and max.
    For categorical variables, summary statistics include count, unique, top (most frequent), and frequency.
    
    
    Parameters:
    -----------
    df : pd.DataFrame
        imput DataFrame 
    
    Returns:
    --------
    pd.DataFrame
        summary DataFrame
    """
    summary_rows = []
    total_rows = len(df)
    
    for col in df.columns:
        col_data = df[col]
        missing_count = col_data.isnull().sum()
        missing_pct = (missing_count / total_rows) * 100
        
        # For numeric variables, use describe() to get summary statistics.
        if pd.api.types.is_numeric_dtype(col_data):
            desc = col_data.describe()
            row = {
                "Variable": col,
                "Type": "Numeric",
                "Count": desc.get("count", np.nan),
                "Missing": missing_count,
                "Missing (%)": round(missing_pct, 2),
                "Mean": round(desc.get("mean", np.nan), 2),
                "Std": round(desc.get("std", np.nan), 2),
                "Min": desc.get("min", np.nan),
                "25%": desc.get("25%", np.nan),
                "50%": desc.get("50%", np.nan),
                "75%": desc.get("75%", np.nan),
                "Max": desc.get("max", np.nan),
                "Unique": np.nan,
                "Top": np.nan,
                "Freq": np.nan
            }
        # For categorical variables, capture count, unique values, most frequent value, and its frequency.
        else:
            desc = col_data.describe()
            row = {
                "Variable": col,
                "Type": "Categorical",
                "Count": desc.get("count", np.nan),
                "Missing": missing_count,
                "Missing (%)": round(missing_pct, 2),
                "Mean": np.nan,
                "Std": np.nan,
                "Min": np.nan,
                "25%": np.nan,
                "50%": np.nan,
                "75%": np.nan,
                "Max": np.nan,
                "Unique": desc.get("unique", np.nan),
                "Top": desc.get("top", np.nan),
                "Freq": desc.get("freq", np.nan)
            }
        
        summary_rows.append(row)
    
    overview_df = pd.DataFrame(summary_rows)
    return overview_df

def preprocess_survival_data(df: pd.DataFrame, time_col: str, event_col: str, therapy_col: str) -> pd.DataFrame:
    """
    Preprocess the survival dataset:
    Convert the time col to numeric. convert the event col to binary convert the therapy indicator to binary.
    
    Parameters:
    -----------
    df : pd.DataFrame
         input DataFrame.
    time_col : str
        col name for time-toevent
    event_col : str
        col name for event indicator
    therapy_col : str
        col name for therapy indicator
    
    Returns:
    --------
    pd.DataFrame
        Processed datframe
    """
    # make sure time column is numeric.
    df[time_col] = pd.to_numeric(df[time_col], errors='coerce')
    
    if not pd.api.types.is_numeric_dtype(df[event_col]):
        df[event_col] = df[event_col].apply(lambda x: 1 if str(x).strip().lower() == "yes" else 0)
    
    # convert the therapy column if it is not already numeric (fixing bug)
    if not pd.api.types.is_numeric_dtype(df[therapy_col]):
        df[therapy_col] = df[therapy_col].apply(lambda x: 1 if str(x).strip().lower() == "yes" else 0)
    
    
    return df

def plot_km_curve(df: pd.DataFrame, time_col: str, event_col: str, therapy_col: str, therapy_label: str, save_path: Path = None):
    """
    generates Kaplan-Meier survival curves for two groups defined by the therapy variable,

    
    

    Returns:
    --------
    results : lifelines.statistics.StatisticalResult
        result of the log-rank test 
    """
    kmf = KaplanMeierFitter()
    
    # Split the data into two groups.
    group_yes = df[df[therapy_col] == 1]
    group_no = df[df[therapy_col] == 0]
    
    plt.figure(figsize=(10, 6))
    
    
    kmf.fit(durations=group_yes[time_col], event_observed=group_yes[event_col], label=f"{therapy_label}: Yes")
    ax = kmf.plot(ci_show=True)
    
   
    kmf.fit(durations=group_no[time_col], event_observed=group_no[event_col], label=f"{therapy_label}: No")
    kmf.plot(ax=ax, ci_show=True)
    
    plt.title(f"Kaplan-Meier Survival Curve by {therapy_label}")
    plt.xlabel("Time (months)")
    plt.ylabel("Survival Probability")
    plt.legend()
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path)
        logger.info(f"Kaplan-Meier plot saved to {save_path}")
    

    
    # Conduct a log-rank test to compare the survival curves.
    results = logrank_test(
        group_yes[time_col],
        group_no[time_col],
        event_observed_A=group_yes[event_col],
        event_observed_B=group_no[event_col]
    )
    logger.info(f"Log-rank test p-value for {therapy_label}: {results.p_value:.4f}")
    
    return results


def preprocess_cox_data(df: pd.DataFrame, time_col: str, event_col: str, therapy_cols: list, additional_covariates: list = None) -> pd.DataFrame:

    # drop rows with missing values
    df = df.dropna(subset=[time_col, event_col] + therapy_cols + (additional_covariates if additional_covariates else []))
    
    # time column to numeric.
    df[time_col] = pd.to_numeric(df[time_col], errors='coerce')
    
    # event indicator to numeric if not already.
    if not pd.api.types.is_numeric_dtype(df[event_col]):
        df[event_col] = df[event_col].apply(lambda x: 1 if str(x).strip().lower() == "yes" else 0)
    
    #therapy columns to numeric 
    for col in therapy_cols:
        if not pd.api.types.is_numeric_dtype(df[col]):
            df[col] = df[col].apply(lambda x: 1 if str(x).strip().lower() == "yes" else 0)
    
    
    cols = [time_col, event_col] + therapy_cols
    if additional_covariates:
        cols.extend(additional_covariates)
    df_cox = df[cols].copy()
    
    # Ccategorical covariates into dummy 
    if additional_covariates:
        df_cox = pd.get_dummies(df_cox, columns=additional_covariates, drop_first=True)
    
    return df_cox

def fit_cox_model(df: pd.DataFrame, time_col: str, event_col: str) -> CoxPHFitter:

    cph = CoxPHFitter()
    cph.fit(df, duration_col=time_col, event_col=event_col)
    cph.print_summary()  # Print hazard ratios, confidence intervals, and p-values.
    return cph

def compare_body_composition_by_therapy(df: pd.DataFrame, therapy_col: str, bc_vars: list) -> None:
    """
    cmpare body composition vars betw groups defined by  therapy using  Mann-Whitney U test.
    
    Parameters:
    -----------
    df : pd.DataFrame
        DataFramee containing body composition and therapy variables
    therapy_col : str
        Name of  therapy variable
    bc_vars : list
        List of body composition variables
    
    """
    # unique values and counts for therapy variable (bug)
    unique_vals = df[therapy_col].unique()
  
    
    for var in bc_vars:
        try:
            col_numeric = pd.to_numeric(df[var], errors='coerce')
            
            group_yes = col_numeric[df[therapy_col] == 1].dropna()
            group_no  = col_numeric[df[therapy_col] == 0].dropna()
            
            n_yes = group_yes.shape[0]
            n_no  = group_no.shape[0]
            
            # Perform the Mann-Whitney U test.
            stat, p_value = mannwhitneyu(group_yes, group_no, alternative='two-sided')
            logger.info(f"{var} - Mann-Whitney U test statistic: {stat:.4f}, p-value: {p_value:.4f}")
        except Exception as e:
            logger.error(f"Error testing {var} for {therapy_col}: {e}")

def calculate_ratio(df: pd.DataFrame, intramuscular_col: str, muscle_col: str) -> pd.DataFrame:

    
    df[intramuscular_col] = pd.to_numeric(df[intramuscular_col], errors='coerce')
    df[muscle_col] = pd.to_numeric(df[muscle_col], errors='coerce')
    
    df['im_to_muscle_ratio'] = df[intramuscular_col] / df[muscle_col]
    logger.info("Calculated intramuscular adipose to muscle ratio.")
    return df

def roc_analysis_for_ratio(df: pd.DataFrame, ratio_col: str, outcome_col: str):
    
    # Drop missing values
    data = df[[ratio_col, outcome_col]].dropna()
    y_true = data[outcome_col].values
    y_scores = data[ratio_col].values
    
    fpr, tpr, thresholds = roc_curve(y_true, y_scores)
    roc_auc = auc(fpr, tpr)
    
    # Calculate Youden's index and select optimal threshold
    youdens_index = tpr - fpr
    optimal_idx = np.argmax(youdens_index)
    optimal_threshold = thresholds[optimal_idx]
    logger.info(f"Optimal cutoff for {ratio_col} is {optimal_threshold:.4f}")
    
    # Plot ROC curve
    plt.figure(figsize=(8, 6))
    plt.plot(fpr, tpr, label=f'ROC curve (area = {roc_auc:.2f})')
    plt.plot([0, 1], [0, 1], 'k--')
    plt.scatter(fpr[optimal_idx], tpr[optimal_idx], color='red', label=f'Optimal cutoff = {optimal_threshold:.2f}')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve for Ratio Predicting Complications')
    plt.legend(loc='lower right')
    plt.tight_layout()
    save_path = os.path.join(cwd, 'ROC.png')
    plt.savefig(save_path)
    
    return fpr, tpr, thresholds, optimal_threshold

def logistic_regression_adjusted(df: pd.DataFrame, outcome_col: str, ratio_col: str, gender_col: str):

    # Ensure the outcome is numeric
    df[outcome_col] = pd.to_numeric(df[outcome_col], errors='coerce')
    # For gender, assume it's categorical; if needed, convert it.
    # We'll use a formula with C(gender_col) to indicate categorical.
    formula = f"{outcome_col} ~ {ratio_col} + C({gender_col})"
    model = smf.logit(formula=formula, data=df).fit(disp=False)
    print(model.summary())
    return model

def univariable_cox_analysis(df: pd.DataFrame, time_col: str, event_col: str, predictor: str) -> dict:
    
    cph = CoxPHFitter()
   
    df_temp = df[[time_col, event_col, predictor]].dropna().copy()
    
    
    df_temp[time_col] = pd.to_numeric(df_temp[time_col], errors='coerce')
    df_temp[event_col] = pd.to_numeric(df_temp[event_col], errors='coerce')
    df_temp[predictor] = pd.to_numeric(df_temp[predictor], errors='coerce')
    
    try:
        cph.fit(df_temp, duration_col=time_col, event_col=event_col)
        
        summary = cph.summary.loc[predictor]
        hr = np.exp(summary['coef'])
        lower_ci = np.exp(summary['coef lower 95%'])
        upper_ci = np.exp(summary['coef upper 95%'])
        p_value = summary['p']
        concordance = cph.concordance_index_
    except Exception as e:
        logger.error(f"Error analyzing predictor {predictor}: {e}")
        hr, lower_ci, upper_ci, p_value, concordance = np.nan, np.nan, np.nan, np.nan, np.nan

    return {
        "Predictor": predictor,
        "HR": hr,
        "95% CI Lower": lower_ci,
        "95% CI Upper": upper_ci,
        "p-value": p_value,
        "Concordance": concordance
    }

def analyze_body_composition(df: pd.DataFrame, time_col: str, event_col: str, predictors: list) -> pd.DataFrame:

    results = []
    for predictor in predictors:
        result = univariable_cox_analysis(df, time_col, event_col, predictor)
        results.append(result)
    return pd.DataFrame(results)

def plot_adipose(df, adipose_vars):
    # Melt the DataFrame to long format.
    df_melted = pd.melt(df, id_vars=['sex'], value_vars=adipose_vars,
                        var_name='Adipose_Type', value_name='Volume')

    # Convert the Volume column to numeric (in case it isn't already).
    df_melted['Volume'] = pd.to_numeric(df_melted['Volume'], errors='coerce')

    # Create a violin plot with box plot overlay, stratified by gender.
    plt.figure(figsize=(10, 6))
    sns.violinplot(x='Adipose_Type', y='Volume', hue='sex', data=df_melted, split=True,
                inner='box', palette='Set2')
    plt.title("Distribution of Adipose Tissue Compartments by Gender")
    plt.xlabel("Adipose Tissue Compartment")
    plt.ylabel("Volume")
    plt.xticks(rotation=45)
    plt.tight_layout()
    save_path = os.path.join(cwd, 'adipose_visualization.png')
    plt.savefig(save_path)