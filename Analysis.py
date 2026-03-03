def main():
    organ_data = "metabolic_data/organ_data.csv"
    connection_data = "metabolic_data/connection_data.csv"


    given_file = "../metabolic_network.csv"

    threshold = 0.3

    from Matrix_Comparison.Comparison import run_network_comparison
    run_network_comparison(given_path=given_file,
                           organ_data=organ_data,
                           connection_data=connection_data,
                           threshold=threshold)



if __name__ == "__main__":
    main()
