import os
import glob
import xml.etree.ElementTree as ET

def count_ecotypes_in_file(file_path):
    try:
        tree = ET.parse(file_path)
        root = tree.getroot()
        # Adjust the tag below if 'ecotype' is nested within other elements
        ecotypes = root.findall('.//ecotype')
        return len(ecotypes)
    except ET.ParseError as e:
        print(f"Error parsing {file_path}: {e}")
        return None

def summarize_ecotypes_in_folder(folder_path):
    xml_files = glob.glob(os.path.join(folder_path, '*.xml'))
    summary = {}

    for xml_file in xml_files:
        count = count_ecotypes_in_file(xml_file)
        if count is not None:
            summary[os.path.basename(xml_file)] = count

    return summary

if __name__ == "__main__":
    # Replace 'path_to_your_folder' with the actual path to your XML files
    folder_path = 'ecosim_results'
    summary = summarize_ecotypes_in_folder(folder_path)

    print("Ecotype Counts per File:")
    for filename, count in summary.items():
        print(f"{filename}: {count} ecotypes")
