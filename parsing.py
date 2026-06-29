import xml.etree.ElementTree as ET
from pathlib import Path

def count_ecotypes_in_file(file_path):
    try:
        return len(ET.parse(file_path).getroot().findall('.//ecotype'))
    except ET.ParseError as e:
        print(f"Error parsing {file_path}: {e}")
        return None

def summarize_ecotypes_in_folder(folder_path):
    summary = {}
    for xml_file in Path(folder_path).glob('*.xml'):
        count = count_ecotypes_in_file(xml_file)
        if count is not None:
            summary[xml_file.name] = count
    return summary

if __name__ == "__main__":
    # Replace 'path_to_your_folder' with the actual path to your XML files
    folder_path = 'ecosim_results'
    summary = summarize_ecotypes_in_folder(folder_path)

    print("Ecotype Counts per File:")
    for filename, count in summary.items():
        print(f"{filename}: {count} ecotypes")
