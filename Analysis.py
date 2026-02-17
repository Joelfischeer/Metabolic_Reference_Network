def main():
    reference_file = "../reference_network.csv"
    given_file = "../metabolic_network.csv"

    metabolic_data_folder = "../metabolic_data"
    threshold = 0.3

    from Matrix_Comparison.Comparison import run_network_comparison
    run_network_comparison(reference_path=reference_file, 
                           given_path=given_file,
                           metabolic_data_folder=metabolic_data_folder,
                           threshold=threshold)



if __name__ == "__main__":
    main()
