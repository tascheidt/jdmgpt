from selenium import webdriver
from selenium.webdriver.chrome.options import Options
import csv
import os

def setup_selenium(download_directory):
    chrome_options = Options()
    chrome_options.add_argument("--headless")  # Enable headless mode
    prefs = {"download.default_directory": os.path.abspath(download_directory)}
    chrome_options.add_experimental_option("prefs", prefs)
    driver = webdriver.Chrome(options=chrome_options)
    return driver

def download_pdf_with_selenium(driver, url):
    print(f"Downloading {url}")
    driver.get(url)
    # Add logic here to wait for the download to complete if necessary

def read_and_update_csv(csv_filepath, download_directory, driver, max_downloads=None):
    updated_rows = []
    count = 0

    # Read all rows and decide which to download
    with open(csv_filepath, 'r') as file:
        reader = csv.reader(file)
        header = next(reader)
        for row in reader:
            pmcid, pdf_url, download_status = row[0], row[3], row[4]
            if download_status != 'Downloaded' and (max_downloads is None or count < max_downloads):
                file_path = os.path.join(download_directory, f"{pmcid}.pdf")
                if not os.path.exists(file_path):
                    download_pdf_with_selenium(driver, pdf_url)
                    row[4] = 'Downloaded'
                    count += 1
                else:
                    row[4] = 'Already Downloaded'
            updated_rows.append(row)

    # Write all rows back to the file
    with open(csv_filepath, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(header)
        writer.writerows(updated_rows)

# Usage
csv_filepath = 'data/articles.csv'
download_directory = '/Users/tscheidt/Documents/code/jdmgpt/downloaded_articles'
driver = setup_selenium(download_directory)
max_downloads = 5000
read_and_update_csv(csv_filepath, download_directory, driver, max_downloads)
driver.quit()
