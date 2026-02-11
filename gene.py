import json
import csv


def extract_uniprot_to_csv(json_filepath, csv_output_path):
    # Carica il file JSON fornito
    with open(json_filepath, 'r') as f:
        data = json.load(f)

    # Usa un set per raccogliere gli accession ID in modo da evitare duplicati
    uniprot_accessions = set()

    # Itera su ogni voce del dizionario principale
    for entry in data.values():
        # Itera sui valori del sotto-dizionario (es. {"A": "Q9Y653", "B": "P09471", ...})
        for accession in entry.values():
            if accession:
                uniprot_accessions.add(accession)

    # Ordina gli accession ID per coerenza
    sorted_accessions = sorted(list(uniprot_accessions))

    # Scrive i dati nel file CSV
    with open(csv_output_path, 'w', newline='') as csvfile:
        fieldnames = ['PROTEIN_AN']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for acc in sorted_accessions:
            writer.writerow({'PROTEIN_AN': acc})

    print(f"File CSV creato con successo: {csv_output_path}")
    print(f"Totale UniProt accession estratti: {len(sorted_accessions)}")


# Esecuzione dello script
extract_uniprot_to_csv('complex_mapping.json', 'proteins.csv')