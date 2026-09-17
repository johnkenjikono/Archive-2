import xml.etree.ElementTree as ET
import pandas as pd
import os
import sys

def ecosim_xml_to_excel(xml_path, output_excel_path):
    tree = ET.parse(xml_path)
    root = tree.getroot()

    ecotypes = root.find(".//demarcation/ecotypes")
    if ecotypes is None:
        print(f"⚠️ No ecotypes found in {xml_path}")
        return

    data = []

    for eco in ecotypes.findall("ecotype"):
        ecotype_num = eco.attrib.get("number")
        size = eco.attrib.get("size")
        members = [m.attrib['name'] for m in eco.findall("member")]
        row = [ecotype_num, size] + members
        data.append(row)

    if not data:
        print(f"⚠️ No ecotype data to export in {xml_path}")
        return

    max_len = max(len(row) for row in data)
    for row in data:
        row += [""] * (max_len - len(row))

    columns = ["Ecotype", "Size"] + [f"Member_{i+1}" for i in range(max_len - 2)]
    df = pd.DataFrame(data, columns=columns)
    df.to_excel(output_excel_path, index=False)

def batch_parse_ecosim_xmls(input_folder):
    output_folder = os.path.join(input_folder, "parsed_results")
    if not os.path.exists(input_folder):
        print(f"❌ Input folder '{input_folder}' does not exist.")
        return

    os.makedirs(output_folder, exist_ok=True)

    processed_count = 0
    for file in os.listdir(input_folder):
        if file.endswith(".xml"):
            xml_path = os.path.join(input_folder, file)
            out_path = os.path.join(output_folder, f"{os.path.splitext(file)[0]}.xlsx")
            try:
                ecosim_xml_to_excel(xml_path, out_path)
                print(f"✅ Processed: {file}")
                processed_count += 1
            except Exception as e:
                print(f"❌ Failed to process {file}: {e}")

    if processed_count == 0:
        print("⚠️ No XML files found to process.")
    else:
        print(f"🎉 Successfully parsed {processed_count} XML files into {output_folder}")

if __name__ == "__main__":
    # If a folder path is given on the command line, use it!
    # Otherwise, default to "ecosim_output_rarefaction" based on the Rarefaction pipeline.
    input_folder = "ecosim_output_rarefaction"
    
    if len(sys.argv) > 1:
        input_folder = sys.argv[1]
    
    print(f"Starting post-processing on folder: '{input_folder}'")
    batch_parse_ecosim_xmls(input_folder)
