import argparse
import os
import requests
from bs4 import BeautifulSoup as bs
from urllib.parse import urlparse, urlunparse
import json


def normalize_url(url: str) -> str:
    """
    Normalize URL: remove trailing # or anything, ensure it ends with a slash.
    """
    # Parse the URL
    parsed = urlparse(url)

    # Remove fragment (#) and query parameters
    normalized = urlunparse((parsed.scheme, parsed.netloc, parsed.path, '', '', ''))  # params  # query  # fragment

    # Ensure it ends with a slash
    if not normalized.endswith('/'):
        normalized += '/'

    return normalized


def scrape_addgene_kit(url: str) -> str:
    """
    Fetch Addgene kit page content using a simple HTTP request.
    """
    response = requests.get(url)
    response.raise_for_status()  # Raise an exception for bad status codes
    return response.text


def sanitize_string(text: str) -> str:
    """Sanitize string by trimming whitespace."""
    return text.strip()


def extract_plasmids_from_page(page_content: str) -> list[tuple[str, str, str, str]]:
    """
    Extract plasmid information from Addgene kit page HTML.

    Replicates the JavaScript logic:
    - Finds the last table with class 'kit-inventory-table' in '#kit-contents'
    - Validates headers are 'Well', 'Plasmid', 'Resistance'
    - Extracts plasmid data from tbody rows
    """
    soup = bs(page_content, "html.parser")

    kit_title = soup.find("h1", id="kit-title")
    if not kit_title:
        raise Exception("The html document does not have the expected structure: #kit-title not found")
    kit_title = kit_title.get_text(strip=True)

    # Find all tables with class 'kit-inventory-table' inside '#kit-contents'
    kit_contents = soup.find(id="kit-contents")
    if not kit_contents:
        raise Exception("The html document does not have the expected structure: #kit-contents not found")

    elements = kit_contents.find_all("table", class_="kit-inventory-table")
    if not elements:
        raise Exception("The html document does not have the expected structure: kit-inventory-table not found")

    # Get the last table (sometimes there is a table with only headers on top)
    kit_table = elements[-1]

    # Check that the headers of the table are Well, Plasmid and Resistance
    headers = kit_table.select("thead th")
    if len(headers) != 3:
        raise Exception("The html document does not have the expected structure: expected 3 headers")

    header_text = [header.get_text(strip=True) for header in headers]
    if header_text[0] != "Well" or header_text[1] != "Plasmid" or header_text[2] != "Resistance":
        raise Exception("The table does not have the expected headers")

    # Get the plasmid info
    plasmids = []
    tbody = kit_table.find("tbody")
    if not tbody:
        raise Exception("The table does not have a tbody")

    rows = tbody.find_all("tr")
    for row in rows:
        cells = row.find_all("td")
        if len(cells) < 3:
            continue  # Skip incomplete rows

        # Extract addgene_id from the link in the second cell
        second_cell = cells[1]
        link = second_cell.find("a")
        if not link or not link.get("href"):
            continue  # Skip rows without links

        href = link.get("href")
        # Extract ID from URL (handle trailing slash)
        # For now, it seems that ids are /id/ (they have trailing slash)
        # But it is better to be safe and pop twice, since first pop will
        # be false if there is no trailing slash
        id_href = href.split("/")
        # Pop twice to handle trailing slash (first pop may be empty string)
        # Match JavaScript behavior: pop() returns undefined for empty array (falsy)
        first_pop = id_href.pop() if id_href else None
        addgene_id = first_pop or (id_href.pop() if id_href else None)

        if not addgene_id:
            continue  # Skip if we can't extract ID

        plasmids.append(
            (
                sanitize_string(cells[0].get_text()),
                sanitize_string(cells[1].get_text()),
                addgene_id,
                sanitize_string(cells[2].get_text()),
            )
        )
    return kit_title, plasmids


def get_kit_dirname(url: str) -> str:
    "Turn the url into a directory name"
    parsed = urlparse(url)
    path = parsed.path
    if path.startswith("/kits/"):
        path = path[6:]
    # Remove the leading and trailing slashes
    path = path.strip("/")
    # Replace slashes with underscores
    path = path.replace("/", "_")
    return path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Scrape Addgene kit page and extract plasmid information")
    parser.add_argument("url", help="URL of the Addgene kit page")
    args = parser.parse_args()

    if not args.url.startswith("https://www.addgene.org/"):
        raise ValueError("The URL must start with https://www.addgene.org/")

    normalized_url = normalize_url(args.url)
    kit_name = get_kit_dirname(normalized_url)

    page_content = scrape_addgene_kit(normalized_url)
    kit_title, plasmids = extract_plasmids_from_page(page_content)

    if not os.path.exists(f"kits/{kit_name}"):
        os.makedirs(f"kits/{kit_name}")
    with open(f"kits/{kit_name}/plasmids.tsv", "w") as f:
        f.write("\t".join(["well", "name", "addgene_id", "resistance"]) + "\n")
        for plasmid in plasmids:
            f.write("\t".join(plasmid) + "\n")

    with open(f"kits/{kit_name}/info.json", "w") as f:
        json.dump(
            {
                "title": kit_title,
                "url": normalized_url,
            },
            f,
            indent=4,
        )
