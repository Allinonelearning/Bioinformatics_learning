#导入各种库
import requests
from lxml import etree
from lxml import html
from html.parser import HTMLParser
from lxml.html import fromstring, tostring
from certifi import where
import urllib3
http = urllib3.PoolManager(num_pools=50)
requests.packages.urllib3.disable_warnings()

# 怎么发送
# 发送给谁
url = "https://book.douban.com/top250"

# 伪装自己
# while True:
headers = {
    'User-Agent':'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/128.0.0.0 Safari/537.36',
    'Accept-Language': 'zh-CN,zh;q=0.9',
    'Accept':'text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.7'
}

# 发送请求
Resp = requests.get(url,headers= headers)
Resp.encoding = "utf-8"

# 响应应信息
e = etree.HTML(Resp.text)
e1 = html.tostring(e[0],encoding='utf-8').decode('utf-8')
print(Resp.text)
# print(Resp.text)
Book_name = e.xpath("//div[@class = \"indent\"]/table[1]/tbody[1]/tr[1]/td[2]/div[1]/a[1]/@title")
print(Book_name)
    # info = '\n'.join(e.xpath("//div[@class =\"txtnav\"]/p/text()"))
    # title = e.xpath("//h1/text()")[0]
    # url = "https://www.douluo1.com{}".format(e.xpath("//div[@class='page1']/a[4]/@href")[0])


# 保存
#     with open("斗罗大陆.txt", "a", encoding='utf-8') as f:
#         f.write(title + '\n\n' + info + '\n\n')
#     if url == "https://book.douban.com/top250?start=225":
#         break





