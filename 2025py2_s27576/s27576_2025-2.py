#!/usr/bin/env python3
"""
Enhanced NCBI GenBank Data Retriever

Funkcjonalności:
- Filtrowanie długości sekwencji (min/max)
- Generowanie raportu CSV
- Tworzenie wykresu długości sekwencji
"""

from Bio import Entrez, SeqIO
import time
import pandas as pd
import matplotlib.pyplot as plt
import ssl
ssl._create_default_https_context = ssl._create_unverified_context


def run_search(tax_id, min_len, max_len, email, api_key):
    Entrez.email = email
    Entrez.api_key = api_key

    # Utworzenie zapytania z filtrem długości sekwencji
    query = f"txid{tax_id}[Organism] AND {min_len}:{max_len}[Sequence Length]"
    handle = Entrez.esearch(db="nucleotide", term=query, usehistory="y", retmax=0)
    result = Entrez.read(handle)

    # Zwracamy liczbę rekordów oraz informacje potrzebne do pobrania wyników
    return int(result["Count"]), result["WebEnv"], result["QueryKey"]


def fetch(webenv, query_key, start, batch_size=100):
    handle = Entrez.efetch(
        db="nucleotide", rettype="gb", retmode="text",
        retstart=start, retmax=batch_size,
        webenv=webenv, query_key=query_key
    )

    return list(SeqIO.parse(handle, "gb"))


def save_outputs(records, tax_id):
    # Tworzy DataFrame z danych i sortuje wg długości
    data = [{"accession": rec.id, "length": len(rec.seq), "description": rec.description} for rec in records]
    df = pd.DataFrame(data).sort_values("length", ascending=False)

    # Zapis danych do pliku CSV
    csv_path = f"taxid_{tax_id}.csv"
    df.to_csv(csv_path, index=False)

    # Generowanie wykresu długości sekwencji
    plt.figure(figsize=(12, 6))
    plt.plot(df["length"].values, marker='o')
    accessions = df["accession"].values
    spacing = max(1, len(accessions) // 20)
    plt.xticks(range(0, len(accessions), spacing), accessions[::spacing], rotation=90, fontsize=6)
    plt.xlabel("Accession")
    plt.ylabel("Sequence Length")
    plt.title(f"Sequence lengths - TaxID {tax_id}")
    plt.tight_layout()
    plt.savefig(f"taxid_{tax_id}.png")


def main():
    # Pobieranie danych od użytkownika
    email = input("Enter NCBI email: ")
    api_key = input("Enter NCBI API key: ")
    tax_id = input("Enter TaxID of the organism: ")
    min_len = input("Enter min sequence length: ")
    max_len = input("Enter max sequence length: ")

    # Wyszukiwanie rekordów
    count, webenv, qkey = run_search(tax_id, min_len, max_len, email, api_key)

    if count == 0:
        print("No records found:(")
        return

    print(f"Found {count} records, downloading...")

    all_records = []
    for i in range(0, min(count, 1000), 100):
        batch = fetch(webenv, qkey, i)
        all_records.extend(batch)
        time.sleep(0.4)

    save_outputs(all_records, tax_id)
    print("Done")

if __name__ == "__main__":
    main()


