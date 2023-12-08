import base64
import functions_framework
import os
from openai import OpenAI
import json
from flask import abort, jsonify
from linebot import LineBotApi
from linebot.models import TextSendMessage, FlexSendMessage

LINE_CHANNEL_ACCESS_TOKEN = os.environ.get("LINE_CHANNEL_ACCESS_TOKEN", "Specified environment variable is not set.")

line_bot_api = LineBotApi(LINE_CHANNEL_ACCESS_TOKEN)

# Triggered from a message on a Cloud Pub/Sub topic.
@functions_framework.cloud_event
def sending_pubmed(cloud_event):
    # Print out the data from Pub/Sub, to prove that it worked
    # print(base64.b64decode(cloud_event.data["message"]["data"]))
    
    f = open('message.json', 'r')
    flex_message_json_dict = json.load(f)

    print(flex_message_json_dict["body"]["contents"][0]["text"])

    # line_bot_api.broadcast(TextSendMessage(text='Hello World!'))
    line_bot_api.broadcast(FlexSendMessage(contents=flex_message_json_dict, alt_text="論文検索結果"))
