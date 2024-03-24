
def main(folder_path):
    # Load the ESM model
    model, alphabet = load_model_and_alphabet("esm1_t34_670M_UR50S")

    # Iterate through files in the folder
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        if os.path.isfile(file_path) and filename.endswith(".fasta"):
            print(f"Processing file: {filename}")

            # Get sequences from input file
            sequences = get_sequences(file_path)

            # Predict structures and compare
            for i, seq1 in enumerate(sequences):
                for j, seq2 in enumerate(sequences):
                    if i < j:  # Avoid redundant comparisons
                        print(f"Comparing structures for sequences {i+1} and {j+1}:")
                        pred_structure1 = predict_structure(seq1, model)
                        pred_structure2 = predict_structure(seq2, model)
                        rmsd = compare_structures(pred_structure1, pred_structure2)
                        print(f"RMSD between sequences {i+1} and {j+1}: {rmsd:.2f}")

if __name__ == "__main__":
    folder_path = "sequences"  # Provide the path to the folder containing input FASTA files
    main(folder_path)