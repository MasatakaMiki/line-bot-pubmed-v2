import functions_framework
import os
from openai import OpenAI
import json

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
    return response.choices[0].message.content

    # response = openai.ChatCompletion.create(
    #     model="gpt-4",
    #     messages=[
    #         {"role": "system", "content": "日本語で文章が提供されるので英語に翻訳してください。"},
    #         {"role": "user", "content": "私は今日とても幸せです。"}
    #     ]
    # )

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
