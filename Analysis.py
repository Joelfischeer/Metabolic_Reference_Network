def main():
    reference_file = "../metabolic_data/reference_network.csv"
    organ_data = "../metabolic_data/organ_data.csv"
    connection_data = "../metabolic_data/connection_data.csv"


    given_file = "../metabolic_network.csv"

    threshold = 0.3

    from Matrix_Comparison.Comparison import run_network_comparison
    run_network_comparison(reference_path=reference_file, 
                           given_path=given_file,
                           organ_data=organ_data,
                           connection_data=connection_data,
                           threshold=threshold)



if __name__ == "__main__":
    main()
