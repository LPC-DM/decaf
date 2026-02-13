import json

# Input and output file names
input_name = 'Mz2000_mhs130_Mdm500_r1_impacts.json'
output_name = f'filtered_{input_name}'

# Load the JSON data from the input file
with open(input_name, 'r') as infile:
    data = json.load(infile)

# List of names to exclude
exclude_names = [
    'tt_norm2016pass',
    'tt_norm2017pass',
    'tt_norm2018pass',
    'wjets_norm2016pass',
    'wjets_norm2017pass',
    'wjets_norm2018pass'
]

# Filter the data, excluding elements with the names in the exclude_names list
filtered_params = [i for i in data['params'] if i['name'] not in exclude_names]

# Update the data with the filtered list
data['params'] = filtered_params

# Save the filtered data back to a new JSON file
with open(output_name, 'w') as outfile:
    json.dump(data, outfile, indent=4)

print(f"Filtered data saved to {output_name}")
