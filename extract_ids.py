# Description: Step1 - Extracts PubMed IDs and titles from PubMed search results and saves them to a CSV file

import os
import csv
from Bio import Entrez

email = "tscheidt@gmail.com"  # Replace with your actual email

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

def extract_article_info(article):
    try:
        article_info = {
            'pmid': str(article['MedlineCitation']['PMID']),
            'title': article['MedlineCitation']['Article']['ArticleTitle'],
            # Add extraction of additional fields here if needed
        }
        return article_info
    except Exception as e:
        print(f"Error processing article: {e}")
        return None

def save_to_csv(articles, directory, filename):
    filepath = os.path.join(directory, filename)
    if os.path.exists(filepath):
        os.remove(filepath)
    
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    with open(filepath, 'w', newline='', encoding='utf-8') as file:
        writer = csv.writer(file)
        writer.writerow(['PMID', 'Title'])  # Add more header fields here if needed
        
        for article in articles:
            article_info = extract_article_info(article)
            if article_info:
                writer.writerow([article_info['pmid'], article_info['title']])  # Add more fields to write here
            else:
                print("Error in extracting article information")

# Search for articles (including only open access articles)
keyword = "('myositis'[MeSH Terms] OR 'dermatomyositis'[MeSH Terms] OR 'dermatomyositis'[All Fields] OR ('juvenile'[All Fields] AND 'dermatomyositis'[All Fields]) OR 'juvenile dermatomyositis'[All Fields]) AND (1990:2024[pdat]) pubmed pmc open access[filter]'"
pubmed_ids = search_pubmed(keyword, email)

# Fetch details
articles = fetch_details(pubmed_ids, email)

# Save article details to CSV
save_to_csv(articles['PubmedArticle'], 'data', 'articles.csv')

print("Article details have been saved to data/articles.csv")
