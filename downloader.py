import subprocess
import os
import zipfile
import shutil
import glob

def process_species_genome(species_name):
    # 1. Setup Naming
    # Example: "Edwardsiella tarda" -> "E.tarda"
    parts = species_name.split()
    if len(parts) >= 2:
        folder_name = f"{parts[0][0]}.{parts[1]}"
    else:
        folder_name = species_name.replace(" ", "_")
        
    zip_filename = f"{species_name.replace(' ', '_')}.zip"
    
    print(f"--- Processing {species_name} ---")

    # 2. Download via NCBI Datasets
    # Filters for Chromosome/Complete and Genbank (gbff)
    cmd = [
        "datasets", "download", "genome", "taxon", species_name,
        "--assembly-level", "chromosome,complete",
        "--include", "gbff",
        "--filename", zip_filename
    ]
    
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError:
        print(f"Failed to download {species_name}")
        return

    # 3. Unzip
    extract_path = "temp_ext"
    with zipfile.ZipFile(zip_filename, 'r') as zip_ref:
        zip_ref.extractall(extract_path)

    # 4. Define paths based on NCBI's folder structure
    # Standard structure is: ncbi_dataset/data/GCA_000.../genomic.gbff
    data_path = os.path.join(extract_path, "ncbi_dataset", "data")
    
    if os.path.exists(data_path):
        # 5. Create the final species directory (e.g., E.tarda)
        os.makedirs(folder_name, exist_ok=True)

        # 6. Move only the files (the .gbff files) into the species folder
        # This replaces the 'find -exec mv' logic
        for root, dirs, files in os.walk(data_path):
            for file in files:
                if file.endswith(".gbff"): # Targeting Genbank files
                    old_file_path = os.path.join(root, file)
                    new_file_path = os.path.join(folder_name, file)
                    shutil.move(old_file_path, new_file_path)
                    print(f"Moved: {file} to {folder_name}/")

    # 7. Cleanup
    # Remove the zip, the temp extraction folder, and any remaining 'data' artifacts
    if os.path.exists(zip_filename):
        os.remove(zip_filename)
    if os.path.exists(extract_path):
        shutil.rmtree(extract_path)
        
    print(f"Done. Files are in ./{folder_name}/")

# Example usage with your specific names
species_list = ["Bordetella bronchiseptica"]

for s in species_list:
    process_species_genome(s)