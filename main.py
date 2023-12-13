from PyPDF2 import PdfReader
import requests
from Bio import Entrez
import pandas as pd
import numpy as np
from requests_html import HTMLSession
from requests.exceptions import ConnectionError
import xml.etree.ElementTree as ET
from pathlib import Path
import re
from metapub import PubMedFetcher
import boto3
import random
import json
from datetime import datetime, timedelta

# get config varables
config = json.load(open('config.json'))
keys = json.load(open('keys.json'))
session = boto3.Session(
    aws_access_key_id=keys['aws_access_key_id'],
    aws_secret_access_key=keys['aws_secret_access_key']
)
s3 = session.resource('s3')
fetch = PubMedFetcher()
random.seed(10)
# UDFS
# search_with_exclusion_criteria
def search(query):
   Entrez.email = 'example@email.com'
   exclusion_criteria = config['exclusion_criteria']
   query = exclusion_criteria + query
   handle = Entrez.esearch(db='pubmed',
                           sort='relevance',
                           retmax='250000',
                           retmode='xml',
                           term=query)
   results = Entrez.read(handle)
   return results


def fetch_details(id_list):
        ids = ','.join(id_list)
        Entrez.email = 'example@email.com'
        handle = Entrez.efetch(db='pubmed',
                            retmode='xml',
                            id=ids)
        results = Entrez.read(handle)
        return results

def util_fun(dataframe):
        #Standardizing months
        dataframe['Month'].replace('Jan', '01', inplace=True)
        dataframe['Month'].replace('Feb', '02', inplace=True)
        dataframe['Month'].replace('Mar', '03', inplace=True)
        dataframe['Month'].replace('Apr', '04', inplace=True)
        dataframe['Month'].replace('May', '05', inplace=True)
        dataframe['Month'].replace('Jun', '06', inplace=True)
        dataframe['Month'].replace('Jul', '07', inplace=True)
        dataframe['Month'].replace('Aug', '08', inplace=True)
        dataframe['Month'].replace('Sep', '09', inplace=True)
        dataframe['Month'].replace('Oct', '10', inplace=True)
        dataframe['Month'].replace('Nov', '11', inplace=True)
        dataframe['Month'].replace('Dec', '12', inplace=True)
        dataframe['Month'].replace('No Data', np.nan, inplace=True)
        return dataframe

def pmid2pmcid(email, pmid):
    Entrez.email = email

    handle = Entrez.elink(dbfrom="pubmed", db="pmc", linkname="pubmed_pmc", id=pmid, retmode="text")

    handle_read = handle.read()
    handle.close()

    root = ET.fromstring(handle_read)

    pmcid = "No Article found"

    for link in root.iter('Link'):
        for id in link.iter('Id'):
            pmcid = id.text
    return pmcid

def extract_doi(elements):
    """
    Extracts the DOI value from a list of StringElement objects.

    Args:
        elements (list): A list of StringElement objects.

    Returns:
        str or None: The DOI value, or None if not found.
    """
    test=[]
    #print(elements)
    for element in elements:
        #print(element)
        if element.attributes['IdType'] == 'doi':
            #print(element)
            DOI = element
            test.append(str(DOI))
    #print(test)
    return test

def fetch_citations(pmid):
    citation=fetch.article_by_pmid(pmid).citation
    print(citation)
    return citation

def extract_info(text):
    """Extracts name, keyword, and date from a given text."""
    # Use a regex pattern to capture the name
    name_pattern = r"(?P<name>[\w\s]+),?"
    # Use a non-capturing group to ignore the comma
    name_regex = re.compile(name_pattern)
    name = name_regex.findall(text)[0]


    # Use a regex pattern to capture the keyword
    keyword_pattern = r"\b(?P<keyword>et\s+al\.?)\b"
    keyword_regex = re.compile(keyword_pattern)
    if keyword_regex.findall(text)!=[]:
        keyword = keyword_regex.findall(text)[0]
    else:
        keyword='et al'

    # Use a regex pattern to capture the date
    date_pattern = r";\s*(\d{4})$"
    pattern = r"(\d{4});\s*"
    match = re.search(pattern, text)
    if match:
        # Print the extracted year
        date=match.group(1)
    else:
        date=''
        print("Year not found!")
    short_reference=''
    short_reference= name+' '+keyword+';'+date

    return short_reference

#Making a DF with article information
def fetch_pubmed_metadata(keyword_to_be_searched):
    pubmed_ids=[]
    title_list= []
    abstract_list=[]
    journal_list = []
    language_list =[]
    pubdate_year_list = []
    pubdate_month_list = []
    author_list=[]
    pmc_ids=[]
    publication_types=[]
    dois=[]
    dops=[]
    long_ref=[]
    short_ref = []
    studies_with_keyword = search(keyword_to_be_searched)
    studiesIdList = studies_with_keyword['IdList']
    studies = fetch_details(studiesIdList)

    chunk_size = 10000
    for chunk_i in range(0, len(studiesIdList), chunk_size):
        chunk = studiesIdList[chunk_i:chunk_i + chunk_size]
        papers = fetch_details(chunk)
        #print(papers['PubmedArticle'][0]['PubmedData']['ReferenceList'][0]['Reference'][])
        for i, paper in enumerate (papers['PubmedArticle']):

            doi = extract_doi(paper['PubmedData']['ArticleIdList'])
            dois.append(doi)
            try:
                pubmed_ids.append(paper['MedlineCitation']['PMID'])
            except:
                pubmed_ids.append('No IDS')
            try:
                publication_types.append(paper['MedlineCitation']['Article']['PublicationTypeList'])
            except:
                publication_types.append('No Data')
            try:
                pmc_id=pmid2pmcid('example@email.com',paper['MedlineCitation']['PMID'])
                pmc_ids.append(pmc_id)
            except:
                pmc_ids.append("No Article Found.")
            title_list.append(paper['MedlineCitation']['Article']['ArticleTitle'])
            try:
                abstract_list.append(paper['MedlineCitation']['Article']['Abstract']['AbstractText'][0])
            except:
                abstract_list.append('No Abstract')
            journal_list.append(paper['MedlineCitation']['Article']['Journal']['Title'])
            language_list.append(paper['MedlineCitation']['Article']['Language'][0])
            dop = paper['MedlineCitation']['Article']['ArticleDate']
            if dop!=[]:
                date = datetime.strptime(f"{dop[0]['Year']}-{dop[0]['Month']}-{dop[0]['Day']}", "%Y-%m-%d")
                dops.append(date)
            try:
                pubdate_year_list.append(paper['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Year'])
            except:
                pubdate_year_list.append('No Data')
            try:
                pubdate_month_list.append(paper['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Month'])
            except:
                pubdate_month_list.append('No Data')
            try :
                author_names = paper['MedlineCitation']['Article']['AuthorList']
                test=[]
                for j in range(len(author_names)):
                    #print(author_names[j]['ForeName']+' '+author_names[j]['LastName'])
                    test.append(author_names[j]['ForeName']+' '+author_names[j]['LastName'])

                author_list.append(test)
            except:
                author_list.append('No Author')
            try:
                citation=fetch_citations(paper['MedlineCitation']['PMID'])
                reference = extract_info(citation)
                short_ref.append(reference)
                long_ref.append(citation)
            except:
                pass

    #print(long_ref)
    df = pd.DataFrame(list(zip(
        pubmed_ids,pmc_ids,dops,dois,title_list,author_list,abstract_list,long_ref,short_ref,journal_list,language_list, publication_types,pubdate_year_list, pubdate_month_list
        )),
        columns=['PubMed_ID','PMCID','DOP','DOI','Title','Authors', 'Abstract','Citation','Reference','Journal','Language','publication_types', 'Year','Month'])
    #print(df)
    df = util_fun(df)

    return df

def download_article_pdf(df):
  s=HTMLSession()
  headers={'user-agent' : 'Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/114.0.0.0 Safari/537.36'}
  try:
      for index, row in df.iterrows() :
          pubID = row['PubMed_ID']
          pmcid = row['PMCID']
          pmc_id=pmid2pmcid('example@email.com',pubID)
          #print(pdf_path+'/'+pmc_id+'.pdf')
          if pmc_id!=None and pmc_id!='':
              try:
                  base_url='https://www.ncbi.nlm.nih.gov/pmc/articles/'
                  r=s.get(base_url+pmc_id+'/',headers=headers,timeout=5)
                  #print(r.status_code)
                  if r.status_code==200:

                      pdf_url='https://www.ncbi.nlm.nih.gov/'+r.html.find('a.int-view',first=True).attrs['href']
                      r=s.get(pdf_url,stream=True)
                      pdf_file = r.content
                      # s3.Object('gilead-datalake', f'pubmed/pdf/{pubID}.pdf').put(Body=pdf_file)
                      s3.Object(config['bucket_name'], f'pub-med/{pubID}/{pmcid}_{current_date}.pdf').put(Body=pdf_file)

              except ConnectionError as e:
                  pass
  except Exception as e:
      print("An Exception occured : {}".format(e))

def download_article_text(df):
  success = []
  failed = []
  try:
      for index, row in df.iterrows():
          pubID = row['PubMed_ID']
          pmcid = row['PMCID']
          url = "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json/"+pubID+"/unicode"
          response = requests.get(url)
          # print(response.status_code)
          if response.status_code==200:
              s3.Object(config['bucket_name'], f'pub-med/{pubID}/{pmcid}_{current_date}.json').put(Body=response.text)
              success.append(pubID)
          else:
              failed.append(pubID)
              print(f"HTTPS Error: {response.status_code} For pub-med id: {pubID}")
      reponse = {'status': 200,
                 'message': f"Successully downloaded {success} pubids\n Failed to Download {failed}"}
      return reponse
  except Exception as e:
      print("An Exception occured : {}".format(e))

def fetch_additional_metadata_helper_2(PMCID):
  if(PMCID=='No Article found'):
    return("No Article found'")
  else:

    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pmc&id={PMCID}&retmode=json"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        pmc_live_date = data['result'][PMCID]['pmclivedate']
        return pmc_live_date
    else:
        return 'Request failed'

def fetch_additional_metadata_helper(pubmed_id):

  url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={pubmed_id}&retmode=json"
  response = requests.get(url)
  if response.status_code == 200:
      data = response.json()

      # uid = data['result']['uids'][0]
      record = data['result'][pubmed_id]

      sortfirstauthor = record.get('sortfirstauthor', '')
      pages = record.get('pages', '')
      volume = record.get('volume', '')
      issue = record.get('issue', '')
      source = record.get('source', '')
      pubdate = record.get('sortpubdate', '').split(' ')[0]
      availablefromurl = record.get('availablefromurl', '')
      attributes = record.get('attributes', [])


      return sortfirstauthor,pages,volume,issue,source,pubdate,availablefromurl,attributes
  else:
      return 'Request failed', 'Request failed', 'Request failed', 'Request failed', 'Request failed', 'Request failed', 'Request failed', 'Request failed'

def fetch_additional_metadata(df):
  df['pmclivedate'] = df['PMCID'].apply(lambda x: pd.Series(fetch_additional_metadata_helper_2(x)))
  df[['sortfirstauthor', 'pages', 'volume', 'issue', 'source', 'pubdate', 'availablefromurl', 'attributes']] = df['PubMed_ID'].apply(lambda x: pd.Series(fetch_additional_metadata_helper(x)))
  return(df)

# Get the metadata
# Get the current date and time
current_datetime = datetime.now()

# Extract the current date
current_date = current_datetime.date()

# Calculate the date 15 days ago
date_15_days_ago = current_date - timedelta(days=config["past_days"])

# Format the dates in "YYYY/MM/DD"
formatted_current_date = current_date.strftime("%Y/%m/%d")
formatted_date_15_days_ago = date_15_days_ago.strftime("%Y/%m/%d")

for i in config["therapeutic area"]:
  search_keyword=f'''("{formatted_date_15_days_ago}"[Date - Publication] : "{formatted_current_date}"[Date - Publication]) AND ({i})'''
  pubmedMetaData=fetch_pubmed_metadata(search_keyword)
  pubmedMetaData=fetch_additional_metadata(pubmedMetaData)
  # pubmedMetaData.to_csv(f'{i}.csv', index=False)


  csv_buffer = pubmedMetaData.to_csv(index=False)
  s3.Object(config['bucket_name'], f'pub-med/metadata/metadata_{current_date}.csv').put(Body=csv_buffer)
  download_article_pdf(pubmedMetaData[['PubMed_ID','PMCID']])
  download_article_text(pubmedMetaData[['PubMed_ID','PMCID']])










print("End Check")