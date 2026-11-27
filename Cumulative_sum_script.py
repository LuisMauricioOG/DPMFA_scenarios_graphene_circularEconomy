import pandas as pd
import numpy as np
import os

def read_csv_files(folder_path):
    """
    Reads Excel files containing 'Sink' in their filename from the specified folder
    and returns a dictionary of DataFrames.
    """
    # List all files in the directory
    files = os.listdir(folder_path)

    # Dictionary to store DataFrames
    dataframes = {}

    # Loop through the files and read those containing 'Sink' in the filename
    for file in files:
        if 'Sink' in file and file.endswith('.csv'):
            # Create a DataFrame from the Excel file
            df = pd.read_csv(os.path.join(folder_path, file))

            # Create a variable name from the file name (remove the extension)
            var_name = os.path.splitext(file)[0]

            # Store the DataFrame in the dictionary
            dataframes[var_name] = df

    return dataframes


def cumulative_sum(df_dict):
    """
    Computes the cumulative sum for each DataFrame in the provided dictionary,
    across the rows starting from the second column.
    """
    cumulative_sums = {}

    for name, df in df_dict.items():
        # Initialize the cumulative sum DataFrame with the same index and columns as the original
        cumulative_sum_df = pd.DataFrame(index=df.index, columns=df.columns)

        # Copy the first column as is
        cumulative_sum_df.iloc[:, 0] = df.iloc[:, 0]

        # Calculate cumulative sum row-wise starting from the second column
        for i in range(len(df)):
            row_cumsum = df.iloc[i, 1:].cumsum()
            cumulative_sum_df.iloc[i, 1:] = row_cumsum

        cumulative_sums[name] = cumulative_sum_df

    return cumulative_sums


def save_cumulative_sums(cumulative_sums, output_folder):
    """
    Saves the cumulative sum DataFrames into Excel files in the specified folder.
    """
    # Create the output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Save each cumulative sum DataFrame to an Excel file
    for name, df in cumulative_sums.items():
        output_file = os.path.join(output_folder, f"{name}_cumulative_sum.csv")
        df.to_csv(output_file, index=False)



