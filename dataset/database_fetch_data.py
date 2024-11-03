import time

import requests  # 'requests' for making HTTP requests.
from bs4 import BeautifulSoup  # BeautifulSoup for web scraping HTML and XML files.
import psycopg2  # 'psycopg2' for connecting to PostgreSQL databases for fetching data from Chembl.
import mysql.connector as sql  # 'mysql.connector' for interfacing with MySQL databases for fetching data from Pharos.
import pandas as pd


def get_uniprot_species(uniprot_id):
    response = requests.get(f'https://www.uniprot.org/uniprot/{uniprot_id}.txt')
    for line in response.text.split('\n'):
        if line.startswith('OS'):
            return line[5:]


def get_protein_location(uniprot_id):
    """
    Given a UniProt ID, returns the subcellular location of the corresponding protein.

    This function scrapes data from the UniProt website, specifically from the XML data
    of the protein corresponding to the UniProt ID. It finds all subcellular locations
    listed and returns them as a joined string, with each location separated by a comma.

    Note that there is a delay of one second at the beginning of the function to prevent
    overwhelming the server with requests.

    Parameters:
    uniprot_id (str): The UniProt ID of the protein.

    Returns:
    str: A string of the subcellular locations of the protein, separated by commas.

    Example:
    # >>> get_protein_location("Q9Y2G9")
    'Mitochondrion matrix, Mitochondrion inner membrane'
    """
    time.sleep(1)
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.xml"
    response = requests.get(url)
    soup = BeautifulSoup(response.content, 'lxml-xml')

    subcellular_locations = soup.findAll('subcellularLocation')
    locations = [location.get_text(separator=' ', strip=True) for location in subcellular_locations]

    return ', '.join(locations)


def fetch_and_merge_data_pharos(table_names, merge_on):
    """
    Fetches data from the specified tables in a database and merges the resulting dataframes.

    Parameters:
    db_connection (sqlalchemy.engine.base.Connection): The database connection.
    table_names (list of str): The names of the tables to fetch data from.
    merge_on (str): The column name to merge the dataframes on.

    Returns:
    pandas.DataFrame: The merged dataframe.
    """
    db_connection = sql.connect(host='tcrd.kmc.io', db='tcrd540', user='tcrd')
    dataframes = []
    for table in table_names:
        query = f"SELECT * FROM {table}"
        df = pd.read_sql(query, con=db_connection)
        dataframes.append(df)

    merged_df = dataframes[0]
    for df in dataframes[1:]:
        merged_df = pd.merge(merged_df, df, left_on=merge_on, right_on='id', how='inner')

    return merged_df


# Usage example:
table_names = ["drug_activity", "target", "protein"]
merge_on = 'target_id'
df_pharos = fetch_and_merge_data_pharos(table_names, merge_on)
df_pharos.to_csv('pharos.csv')


def fetch_and_merge_data_chembl(database, user, password, host, port, query):
    """
    This function establishes a connection with a PostgreSQL database, executes a provided SQL query,
    and returns a DataFrame of the fetched results.

    :param database: String representing the name of the database in PostgreSQL to be connected to.
    :param user: String representing the username required to authenticate the PostgreSQL connection.
    :param password: String representing the password required to authenticate the PostgreSQL connection.
    :param host: String representing the IP address for establishing connection to the PostgreSQL database.
                 Use 'localhost' if the database is hosted locally.
    :param port: String representing the port number for establishing connection to the PostgreSQL database.
    :param query: String representing the SQL query to be executed on the database.

    :return: A pandas DataFrame containing data fetched from the executed query. Null entries in 'UniProt_ID'
             and 'Drug_Name' are filtered out and duplicate rows are dropped. This filtered DataFrame is also saved
             as a CSV file 'chembl.csv' in the current directory.

    Usage Example:
    --------------
    database='chembl_32'
    user='postgres'
    password='your_password'
    host='localhost'
    port='5432'
    query="SELECT * FROM your_table"

    df = fetch_and_merge_data_chembl(database, user, password, host, port, query)
    """
    conn = psycopg2.connect(
        database=database,
        user=user,
        password=password,
        host=host,
        port=port)

    cur = conn.cursor()

    cur.execute(query)

    rows = cur.fetchall()

    cur.close()
    conn.close()

    df = pd.DataFrame(rows, columns=["drug", "smiles", "uniprot"])
    df = df[(df['uniprot'].notnull()) & (df['drug'].notnull())]
    df.drop_duplicates(inplace=True)

    return df


# using the function (these variables are an example need to be modified according to the own datacent)

database = 'chembl_32'
user = 'postgres'
password = '1qaz'
host = 'localhost'
port = '5432'

query = """
    SELECT md.pref_name, cs.canonical_smiles, csq.accession
    FROM molecule_dictionary md 
    JOIN compound_structures cs ON md.molregno = cs.molregno
    JOIN activities act ON md.molregno = act.molregno
    JOIN assays ass ON act.assay_id = ass.assay_id
    JOIN target_dictionary td ON ass.tid = td.tid
    LEFT JOIN target_components tc ON td.tid = tc.tid
    LEFT JOIN component_sequences csq ON tc.component_id = csq.component_id
    WHERE td.organism = 'Homo sapiens'
"""

# fetch and save the from the chembl_32 dataset 
df_chembl = fetch_and_merge_data_chembl(database, user, password, host, port, query)
df_chembl.to_csv('chembl.csv', index=False)


def fetch_and_merge_data_drug_central(TSV_file_path):
    df = pd.read_csv(TSV_file_path)
    return df[['smiles', 'uniprot']]


df_drug_centeral = fetch_and_merge_data_drug_central('tchem_drugs_05122020.tsv') # it can be downloaded from the drugcenteral.org
df_drug_centeral.to_csv('drug_centeral.csv')
