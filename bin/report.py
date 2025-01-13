import pandas as pd
import glob
import sys

# MERGE ALL DATA
# Get a list of all CSV files in the current directory
csv_files = glob.glob('*.csv')
samplesheet = sys.argv[1]

# LOAD SAMPLE SHEET
samplesheet_df = pd.read_csv(samplesheet, sep='\t')

# Rename the column
samplesheet_df.rename(columns={'SequenceID': 'Sample'}, inplace=True)
# Remove Barcode columns
samplesheet_df = samplesheet_df.drop(columns=['Barcode'])

# Initialize an empty DataFrame to store the merged data
merged_data = None

# Iterate over each CSV file
for i, file in enumerate(csv_files):
    # Read the CSV file into a DataFrame
    df = pd.read_csv(file)

    # If this is the first iteration, initialize merged_data with the first DataFrame
    if merged_data is None:
        merged_data = df
    else:
        # Merge the DataFrame with the merged_data DataFrame based on the "Sample" column
        merged_data = pd.concat([merged_data, df], axis=0, ignore_index=True)

# Group by 'Sample' and combine the rows
merged_data = merged_data.groupby('Sample', as_index=False).first()

# NIPH specific adjustment
# Replace ! with - in the Sample column
merged_data['Sample'] = merged_data['Sample'].str.replace('!', '-')

# Sort alphabetically with 'Sample' first
merged_data = merged_data[['Sample'] + sorted([col for col in merged_data.columns if col != 'Sample'])]

# Write the merged data to a new CSV file
merged_data.to_csv('merged_report.csv', index=False)
