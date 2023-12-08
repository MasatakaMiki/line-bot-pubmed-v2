import functions_framework
import os
from openai import OpenAI
import json
from flask import abort, jsonify
from linebot import (
    LineBotApi, WebhookHandler
)
from linebot.exceptions import (
    InvalidSignatureError
)
from linebot.models import (
    MessageEvent, TextMessage, TextSendMessage
)

LINE_CHANNEL_SECRET = os.environ.get("LINE_CHANNEL_SECRET", "Specified environment variable is not set.")
LINE_CHANNEL_ACCESS_TOKEN = os.environ.get("LINE_CHANNEL_ACCESS_TOKEN", "Specified environment variable is not set.")

line_bot_api = LineBotApi(LINE_CHANNEL_ACCESS_TOKEN)
handler = WebhookHandler(LINE_CHANNEL_SECRET)

@functions_framework.http
def main(request):
   
    signature = request.headers['X-Line-Signature']

    body = request.get_data(as_text=True)

    try:
        handler.handle(body, signature)
    except InvalidSignatureError:
        abort(400)
 
    return 'OK'

@handler.add(MessageEvent, message=TextMessage)
def handle_message(event):

    client = OpenAI(api_key=os.environ.get("OPENAI_API_KEY", "Specified environment variable is not set."))

    prompt = event.message.text

    response = client.chat.completions.create(
        messages=[
            {
                "role": "system",
                "content": "受け取った言語で返答をしてください。",
            },
            {
                "role": "user",
                "content": prompt,
            }
        ],
        # model="gpt-3.5-turbo",
        model="gpt-4-1106-preview",
    )
    print(f'RESPONSE:{response}')

    line_bot_api.reply_message(
        event.reply_token,
        TextSendMessage(text=response.choices[0].message.content))

@functions_framework.http
def hello_http(request):

    client = OpenAI(api_key=os.environ.get("OPENAI_API_KEY", "Specified environment variable is not set."))
    
    request_json = request.get_json(silent=True)
    request_args = request.args

    prompt = request_json['message']

    response = client.chat.completions.create(
        messages=[
            {
                "role": "system",
                "content": "受け取った言語で返答をしてください。",
            },
            {
                "role": "user",
                "content": prompt,
            }
        ],
        model="gpt-3.5-turbo",
    )
    print(f'RESPONSE:{response}')
    return f'ver.0.0.1:{response.choices[0].message.content}'

    # if request.args and 'message' in request.args:
    #     return request.args.get('message')
    # elif request_json and 'message' in request_json:
    #     return json.dumps(response['choices'][0]['text'], ensure_ascii=False, indent=2)
    # else:
    #     return f'Enter something'

    # if request_json and 'name' in request_json:
    #     name = request_json['name']
    # elif request_args and 'name' in request_args:
    #     name = request_args['name']
    # else:
    #     name = 'World'
    # return 'Hello {}!'.format(name)
