import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import os
from pathlib import Path
import torch
print(torch.__version__)


def wilcoxon_analysis(data_df, output_folder, filename_base):
    """
    Perform Wilcoxon rank-sum test between groups for each cell type
    
    Args:
        data_df (DataFrame): DataFrame containing group, cell_type_number, and proportion
        output_folder (str): Path to output folder
        filename_base (str): Base filename for output
    """
    # Get unique cell type numbers
    cell_types = sorted(data_df['cell_type_number'].unique())
    
    # Prepare results container
    test_results = []
    
    # Perform tests between each pair of groups
    group_pairs = [(1,2), (1,3), (2,3)]
    
    for cell_type in cell_types:
        cell_data = data_df[data_df['cell_type_number'] == cell_type]
        
        for group1, group2 in group_pairs:
            group1_data = cell_data[cell_data['group'] == group1]['proportion']
            group2_data = cell_data[cell_data['group'] == group2]['proportion']
            
            if len(group1_data) > 0 and len(group2_data) > 0:
                # Perform Wilcoxon rank-sum test
                statistic, pvalue = stats.ranksums(group1_data, group2_data)
                
                test_results.append({
                    'cell_type': cell_type,
                    'group_comparison': f'Group{group1}_vs_Group{group2}',
                    'statistic': statistic,
                    'pvalue': pvalue
                })
    
    # Convert to DataFrame
    results_df = pd.DataFrame(test_results)
    
    if not results_df.empty:
        # Apply Benjamini-Hochberg correction within each cell type
        results_df['padj'] = float('nan')
        for cell_type in results_df['cell_type'].unique():
            mask = results_df['cell_type'] == cell_type
            results_df.loc[mask, 'padj'] = multipletests(results_df.loc[mask, 'pvalue'], method='fdr_bh')[1]
        
        # Sort by cell type and adjusted p-value
        results_df = results_df.sort_values(['cell_type', 'padj'])
        
        # Save results
        output_filename = os.path.join(output_folder, f'{filename_base}_wilcoxon_results.csv')
        results_df.to_csv(output_filename, index=False)
        
        return results_df
    
    return None

def process_single_csv(input_file, output_folder):
    """
    Process a single CSV file and generate its results
    
    Args:
        input_file (str): Path to input CSV file
        output_folder (str): Path to output folder
    """
    # Read CSV file
    df = pd.read_csv(input_file)
    
    # Get cell types from first column and create mapping
    cell_types = df.iloc[:, 0].unique()
    cell_type_mapping = {cell_type: idx + 1 for idx, cell_type in enumerate(cell_types)}
    
    # Get base filename without extension
    base_filename = os.path.splitext(os.path.basename(input_file))[0]
    
    # Save cell type mapping for this file
    mapping_filename = f"{base_filename}_cell_type_mapping.txt"
    with open(os.path.join(output_folder, mapping_filename), 'w', encoding='utf-8') as f:
        f.write("Cell Type\tNumber\n")
        for cell_type, number in cell_type_mapping.items():
            f.write(f"{cell_type}\t{number}\n")
    
    # Initialize results list
    results = []
    
    # Process each column independently
    for column in df.columns[1:]:  # Skip first column (cell types)
        group = determine_group(column)
        if group is None:
            continue
            
        # Get column total for proportion calculation
        column_total = df[column].sum()
        
        if column_total > 0:
            # Process each cell type for this column
            for cell_type, value in zip(df.iloc[:, 0], df[column]):
                cell_type_number = cell_type_mapping[cell_type]
                proportion = value / column_total if column_total > 0 else 0
                
                results.append([
                    group,
                    cell_type_number,
                    proportion
                ])
    
    # Convert results to DataFrame
    results_df = pd.DataFrame(results, columns=['group', 'cell_type_number', 'proportion'])
    
    # Sort by group and cell type number
    results_df = results_df.sort_values(['group', 'cell_type_number'])
    
    # Save proportion results
    output_filename = f"processed_{base_filename}.csv"
    results_df.to_csv(os.path.join(output_folder, output_filename), index=False, float_format='%.6f')
    
    # Perform statistical analysis
    wilcoxon_analysis(results_df, output_folder, base_filename)
    
    return output_filename, mapping_filename

def determine_group(column_name):
    """
    Determine which group a column belongs to
    
    Args:
        column_name (str): Name of the column
    
    Returns:
        int: Group number (1, 2, or 3)
    """
    # Group 1 includes THD* columns and L_* columns
    if (column_name.startswith('THD') or 
        column_name in ['L_156_UP_CTRL', 'L_155_LOW_CTRL', 'L_45_CTRL', 'L_59_CTRL']):
        return 1
    elif column_name.startswith('GSM'):
        return 2
    elif column_name.startswith('Control_'):
        return 3
    return None

def process_csv_folder(input_folder, output_folder):
    """
    Process all CSV files in the input folder and save results to output folder
    
    Args:
        input_folder (str): Path to folder containing CSV files
        output_folder (str): Path to folder where results will be saved
    """
    # Create output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)
    
    # Process each CSV file
    processed_files = []
    for filename in os.listdir(input_folder):
        if filename.endswith('.csv'):
            print(f"Processing {filename}...")
            try:
                input_file = os.path.join(input_folder, filename)
                result_file, mapping_file = process_single_csv(input_file, output_folder)
                processed_files.append({
                    'input_file': filename,
                    'result_file': result_file,
                    'mapping_file': mapping_file
                })
            except Exception as e:
                print(f"Error processing {filename}: {str(e)}")
                continue
    
    # Create a summary file
    with open(os.path.join(output_folder, 'processing_summary.txt'), 'w', encoding='utf-8') as f:
        f.write("Processing Summary\n")
        f.write("=================\n\n")
        for file_info in processed_files:
            f.write(f"Input file: {file_info['input_file']}\n")
            f.write(f"Result file: {file_info['result_file']}\n")
            f.write(f"Mapping file: {file_info['mapping_file']}\n")
            f.write("-" * 50 + "\n")
    
    print(f"Processing complete. Results saved to {output_folder}")

def main():
    """
    Main function to run the CSV processing pipeline
    """
    # Get input and output folders from user
    input_folder = "E:/MasterDegree/Research/scRNAseq-LungNasal/1229"  # 替换为您的输入文件夹路径
    output_folder = "E:/MasterDegree/Research/scRNAseq-LungNasal/0101" # 替换为您的输出文件夹路径
    
    # Validate input folder exists
    if not os.path.exists(input_folder):
        print(f"Error: Input folder '{input_folder}' does not exist")
        return
    
    try:
        process_csv_folder(input_folder, output_folder)
    except Exception as e:
        print(f"Error processing files: {str(e)}")

if __name__ == "__main__":
    main()