import sys
import pandas as pd


def main():
    files = sys.argv[1:]
    nfiles = len(files)
    dfs = [pd.read_csv(countfile, sep="\t").drop(columns=["Chr", "Start", "End", "Strand", "Length"]) for countfile in files ] 
    print(dfs[0].columns)

    merged = pd.concat([df.set_index("Geneid") for df in dfs], axis=1)
    merged.to_csv("test.txt", sep="\t")


if __name__ == "__main__":
    main()