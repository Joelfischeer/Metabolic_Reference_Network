def main():
    reference_file = "../reference_network.csv"
    given_file = "../metabolic_network.csv"

    threshold = 0.3

    from Matrix_Comparison.Comparison import run_network_comparison
    run_network_comparison(reference_file, given_file)



if __name__ == "__main__":
    main()
