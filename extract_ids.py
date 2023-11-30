# Description: Step1 - Extracts PubMed IDs and titles from PubMed search results and saves them to a CSV file

import os
import csv
import requests
from Bio import Entrez

email = "tscheidt@gmail.com"  # Replace with your actual email

def convert_pmids_to_pmcids(pmids):
    url = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"
    payload = {'ids': ','.join(pmids), 'format': 'json'}
    response = requests.post(url, data=payload)
    if response.status_code == 200:
        return response.json()
    else:
        print(f"Error with batch conversion: {response.status_code}")
        return None

def process_in_batches(id_list, batch_size=200):
    pmcid_mapping = {}
    for i in range(0, len(id_list), batch_size):
        batch = id_list[i:i + batch_size]
        result = convert_pmids_to_pmcids(batch)
        if result and 'records' in result:
            for record in result['records']:
                pmcid_mapping[record['pmid']] = record.get('pmcid', 'Not Available')
    return pmcid_mapping

def search_pubmed(keyword, email, batch_size=100):
    Entrez.email = email
    handle = Entrez.esearch(db="pubmed", term=keyword, retmax=batch_size, usehistory="y")
    record = Entrez.read(handle)
    handle.close()

    # Fetch all results in batches
    count = int(record["Count"])
    print(f"Total articles found: {count}")
    id_list = []
    for start in range(0, count, batch_size):
        handle = Entrez.esearch(db="pubmed", term=keyword, retstart=start, retmax=batch_size, usehistory="y", webenv=record["WebEnv"], query_key=record["QueryKey"])
        record = Entrez.read(handle)
        id_list.extend(record["IdList"])
        handle.close()

    return id_list

def fetch_details(id_list, email):
    ids = ','.join(id_list)
    Entrez.email = email
    handle = Entrez.efetch(db="pubmed", id=ids, retmode="xml", rettype="full")
    results = Entrez.read(handle)
    handle.close()
    return results

def extract_article_info(article, pmcid_mapping):
    try:
        pmid = str(article['MedlineCitation']['PMID'])
        pmcid = pmcid_mapping.get(pmid, 'Not Available')
        pdf_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/pdf/" if pmcid != 'Not Available' else 'Not Available'
        article_info = {
            'pmcid': pmcid,
            'pmid': pmid,
            'title': article['MedlineCitation']['Article']['ArticleTitle'],
            'pdf_url': pdf_url,
            'download_status': 'Not Downloaded'
        }
        return article_info
    except Exception as e:
        print(f"Error processing article: {e}")
        return None

def save_to_csv(articles, directory, filename, pmcid_mapping):
    filepath = os.path.join(directory, filename)
    if os.path.exists(filepath):
        os.remove(filepath)
    
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    with open(filepath, 'w', newline='', encoding='utf-8') as file:
        writer = csv.writer(file)
        writer.writerow(['PMCID', 'PMID', 'Title', 'PDF URL', 'Download Status'])
        
        for article in articles:
            article_info = extract_article_info(article, pmcid_mapping)
            if article_info:
                writer.writerow([article_info['pmcid'], article_info['pmid'], article_info['title'], article_info['pdf_url'], article_info['download_status']])
            else:
                print("Error in extracting article information")

# Main execution
keyword = "('myositis'[MeSH Terms] OR 'dermatomyositis'[MeSH Terms] OR 'dermatomyositis'[All Fields] OR ('juvenile'[All Fields] AND 'dermatomyositis'[All Fields]) OR 'juvenile dermatomyositis'[All Fields]) AND (1990:2024[pdat]) pubmed pmc open access[filter]'"
pubmed_ids = search_pubmed(keyword, email)
pmcid_mapping = process_in_batches(pubmed_ids)
articles = fetch_details(pubmed_ids, email)
save_to_csv(articles['PubmedArticle'], 'data', 'articles.csv', pmcid_mapping)
print("Article details have been saved to data/articles.csv")

