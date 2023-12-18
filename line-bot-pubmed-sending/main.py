import base64
import functions_framework
import os
import json
from datetime import datetime, timedelta
from flask import abort, jsonify
from Bio import Entrez
from openai import OpenAI
from linebot import LineBotApi
from linebot.models import TextSendMessage, FlexSendMessage

Entrez.email = "masataka.miki@gmail.com"

def search_articles(keyword, max_results=10):

    yesterday = (datetime.now() - timedelta(1)).strftime('%Y-%m-%d')
    print(yesterday)
    handle = Entrez.esearch(db='pubmed',
                            sort='relevance',
                            retmax=str(max_results),
                            mindate=yesterday,
                            maxdate=yesterday,
                            term=keyword)
    results = Entrez.read(handle)
    return results


def get_summary(pmid):
    handle = Entrez.esummary(db='pubmed', id=pmid, retmode='xml')
    summary = Entrez.read(handle)
    return summary


def fetch_article(pmid):
    handle = Entrez.efetch(db='pubmed', id=pmid, retmode='xml')
    article = Entrez.read(handle)['PubmedArticle'][0]
    # abstract = article['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
    return article


LINE_CHANNEL_ACCESS_TOKEN = os.environ.get("LINE_CHANNEL_ACCESS_TOKEN", "Specified environment variable is not set.")
line_bot_api = LineBotApi(LINE_CHANNEL_ACCESS_TOKEN)

# Triggered from a message on a Cloud Pub/Sub topic.
@functions_framework.cloud_event
def sending_pubmed(cloud_event):
    # Print out the data from Pub/Sub, to prove that it worked
    # print(base64.b64decode(cloud_event.data["message"]["data"]))
    
    # PubMed
    pubmed_pubtypes = ["Journal Article-aaa", "Books and Documents", "Clinical Trial", "Meta-Analysis", "Randomized Controlled Trial", "Review", "Systematic Review"]
    pubmed_results = search_articles(keyword="periodontal")
    print(pubmed_results)
    index = 0
    paper_title = ""
    paper_abstract = ""
    paper_uri = ""
    for pmid in pubmed_results["IdList"]:
        index += 1
        print(index)
        # print(pmid)
        summary = get_summary(pmid)
        print(summary)
        article = fetch_article(pmid)
        print(article)
        # check pubmed-pubtype
        print(summary[0]['PubTypeList'])
        if not set(summary[0]['PubTypeList']).isdisjoint(pubmed_pubtypes):
            paper_title = article['MedlineCitation']['Article']['ArticleTitle']
            paper_abstract = article['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
            paper_uri = "https://pubmed.ncbi.nlm.nih.gov/{0}/".format(pmid)
            break
    # print(index)

    if index == 0:
        line_bot_api.broadcast(TextSendMessage(text='本日の新着論文はありません'))
        return

    f = open('message.json', 'r')
    flex_message_json_dict = json.load(f)
    # print(flex_message_json_dict["body"]["contents"][0]["text"])
    # print(flex_message_json_dict["body"]["contents"][1]["contents"][0]["contents"][0]["text"])
    # print(flex_message_json_dict["footer"]["contents"][0]["action"]["uri"])

    flex_message_json_dict["body"]["contents"][0]["text"] = paper_title
    flex_message_json_dict["body"]["contents"][1]["contents"][0]["contents"][0]["text"] = paper_abstract
    flex_message_json_dict["footer"]["contents"][0]["action"]["uri"] = paper_uri

    # line_bot_api.broadcast(TextSendMessage(text='Hello World!'))
    line_bot_api.broadcast(FlexSendMessage(contents=flex_message_json_dict, alt_text="論文検索結果"))
