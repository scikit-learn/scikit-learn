import re
import pandas as pd

import jpype
from jpype import JClass

from scipy.stats import kendalltau

def start_jvm():
    # Start the JVM and set up the classpath
    jvm_args = ["-Xmx1g"]  # Set maximum heap size for JVM
    cp = ["./", "./lib/*", "./weka"]  # Set the classpath for Weka and other Java dependencies
    jpype.startJVM(*jvm_args, classpath=cp, convertStrings=True)

def stop_jvm():
    jpype.shutdownJVM()

def load_xarff(file_path):
    """Custom loader to read XARFF files and convert them into a pandas DataFrame."""
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Separate metadata and data
    data_start_idx = lines.index('@data\n') + 1
    metadata = lines[:data_start_idx]
    data_lines = lines[data_start_idx:]

    # Extract column names and types from metadata
    column_names = []
    attribute_info = {}
    for line in metadata:
        if line.startswith('@attribute'):
            # Extract the column name and its type information
            match = re.findall(r'@attribute\s+([\w\d_-]+)\s+(.*)', line)
            if match:
                column_name, attribute_type = match[0]
                column_names.append(column_name)
                attribute_info[column_name] = attribute_type.strip()

    # Parse the data lines into a list of lists
    data = []
    for line in data_lines:
        # Split the line by commas and strip any extra spaces
        split_line = re.split(r',\s*(?![^{}]*\})', line.strip())
        # Handle the RANKING column by converting it to a list
        if split_line[-1].startswith('{'):
            ranking = split_line[-1][1:-1].split('>')
            split_line[-1] = ranking
        data.append(split_line)

    # Create a DataFrame from the parsed data
    df = pd.DataFrame(data, columns=column_names)

    return df, attribute_info

def save_to_xarff(df, file_path, relation_name="dataset", attribute_info=None):
    """Save a Pandas DataFrame as a XARFF file."""
    with open(file_path, 'w') as f:
        # Write the relation name
        f.write(f"@relation {relation_name}\n\n")

        # Write the attribute names and types using the provided attribute_info
        for column in df.columns:
            attribute_type = attribute_info.get(column, 'STRING')
            f.write(f"@attribute {column} {attribute_type}\n")

        # Write the data
        f.write("\n@data\n")
        for _, row in df.iterrows():
            row_data = []
            for value in row:
                if isinstance(value, list):
                    # Convert list back to ranking format
                    value = '>'.join(value)
                row_data.append(value)
            f.write(','.join(map(str, row_data)) + '\n')

def load_dataset_as_Instances(file):
    DataSource = JClass("weka.core.converters.ConverterUtils$DataSource")
    data = DataSource.read(file)
    data.setClassIndex(data.numAttributes() - 1)
    return data

def kendalls_tau(prefs, preds):
    # Using scipy's built-in Kendall's Tau implementation
    tau, _ = kendalltau(prefs, preds)
    return tau