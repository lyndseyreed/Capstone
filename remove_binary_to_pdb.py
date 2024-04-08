import os
import shutil

def remove_binary(input_folder, output_folder):
    for root, dirs, files in os.walk(input_folder):
        # Get the relative path from the input folder to the current folder
        relative_path = os.path.relpath(root, input_folder)
        # Construct the corresponding output folder path
        output_root = os.path.join(output_folder, relative_path)
        # Create the corresponding output folder
        os.makedirs(output_root, exist_ok=True)
        for filename in files:
            input_file_path = os.path.join(root, filename)
            output_file_path = os.path.join(output_root, filename)
            with open(input_file_path, "rb") as input_file:
    