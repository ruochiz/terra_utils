"""
Combine Terra Tables

This script allows user to combine two or more data tables from Terra platform. It assumes the same samples exist in the tables supplied as inputs. It only acts on the columns containing the library names and fastqs. For the library name columns (i.e.: "ATAC_Lib" or "RNA_Lib"), the script copies the library information if it does not exist in all tables. For the fastq columns, the script combines the Google Bucket paths to the fastq files as an array.

This tool accepts tab separated value files (.tsv).

This script requires <requirements-here>

Usage:
    python combine_table.py </tables/to/combine/>

Author:
    Alia Meliki (aliameliki@gmail.com) / Ruochi Zhang (zhangruo@broadinstitute.org)
"""
import argparse
import pandas as pd
import re

def combine_table(tables):
    """
    Combines library and fastq columns from two or more TSV tables

    Args:
        tables (list): A list of paths to the TSV files

    Returns:
        combined_table: DataFrame of the first user supplied table with the fastq columns combined and the library columns filled in.
    """
    open_tables = [pd.read_csv(table, sep="\t") for table in tables]
    # Using the first column as the index, make sure they have the same rows and same order:
    # Set the first column as the index for each table
    for table in open_tables:
        table.index = list(table.iloc[:, 0].values)

    # Find the union of the indices across all tables
    union_index = open_tables[0].index
    for table in open_tables[1:]:
        union_index = union_index.union(table.index)
    print (union_index)
    # Reindex all tables to this union index and fill missing values with NaN
    open_tables = [table.reindex(union_index) for table in open_tables]
    for tab in open_tables:
        tab.iloc[:, 0] = list(tab.index)
    combined_table = open_tables[0] #TODO handle missing fastqs/library

    cols = ["Genome", "PKR", "R1_Subset", "whitelist"] # for those columns fill the first none nan results
    # For each column in cols_to_fill, fill the first non-NaN value across tables
    for col in cols:
        combined_table[col] = pd.concat([tab[col] for tab in open_tables], axis=1).bfill(axis=1).iloc[:, 0]

    # Compatible with outputs from release version and main version of the pipeline
    cols = {"lib": ["ATAC_Lib", "RNA_Lib"], "fastq": ["ATAC_fastq_R1", "ATAC_fastq_R2", "ATAC_raw_fastq_R1", "ATAC_raw_fastq_R2", "RNA_fastq_R1", "RNA_fastq_R2", "RNA_raw_fastq_R1", "RNA_raw_fastq_R2"]}

    for col in combined_table.columns:
        # Replace the lib if needed
        if col in cols["lib"]:
            for idx in combined_table[col].loc[combined_table[col].apply(isLib) == False].index:
                combined_table.loc[idx,col] = find_lib(idx, col, open_tables)
                
        # Combine fastq
        elif col in cols["fastq"]:
            for df in open_tables:
                df[col] =  df[col].apply(lambda x: "[\"\"]" if pd.isna(x) else x)


            for table in open_tables[1:]:
                combined_table[col] = combined_table[col] + table[col]
                combined_table[col] = combined_table[col].str.replace("[\"\"]", "", regex=False).str.replace("][",",",regex=False)
                
    return combined_table

def isLib(val): #val is a string
    val = str(val)
    if "i5" in val or "Ad1" in val:
        val = str(str)
        return True
    else:
        return False

def find_lib(idx, col, open_tables):
    """
    Finds the library name of a sample from a list of all supplied tables. It assumes samples are located on the same row for all tables

    Args:
        idx (int): Index or row num of the sample
        col (str): Column name
        open_tables (list): A list of DataFrames

    Returns:
        lib (str): Library name or None is the library does not exist
    """
    lib = None

    for table in open_tables:
        if isLib(table.loc[idx][col]) == True:
            lib = table.loc[idx][col]
        else:
            continue

    if lib:
        return lib
    else:
        return None # No lib for this modality, ie; an ATAC only or RNA only sample

def main(tables_to_combine, name):
    res = combine_table(tables_to_combine)

    if name is None:
        # Rename entity ID
        entity_name = res.columns[0].strip("entity:").strip("_id")
        entity_name = re.sub("_L00[1-9]", "", entity_name)
    else:
        entity_name = name
    res = res.rename(columns={res.columns[0]:("entity:" + entity_name + "_id")})
    print(entity_name, " " , res.columns[0])
    # exit()
    
    # Save table to current dir
    res.to_csv(f"./{entity_name}.tsv", sep = "\t", index=False)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine ATAC and RNA fastqs from two or more Terra entity tables")
    parser.add_argument("entity_tables", nargs = "*", help = "Two or more paths to Terra entity tables")
    parser.add_argument("--name", help="An optional name, can be None", default=None)
    args = parser.parse_args()
    
    main(args.entity_tables, args.name)
