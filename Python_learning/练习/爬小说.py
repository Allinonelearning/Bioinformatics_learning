# 编写爬虫脚本，首先确定爬取目标，减少重复性工作
# 怎么发送
import requests
from lxml import etree
from certifi import where
import urllib3
http = urllib3.PoolManager(num_pools=50)
requests.packages.urllib3.disable_warnings()

# # urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
url = "https://www.douluo1.com/books/20539/10147194.html"# 发送给谁
while True:# 伪装自己
    headers = {
        'User-Agent':'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/128.0.0.0 Safari/537.36',
        'Accept-Language': 'zh-CN,zh;q=0.9',
        # 'Referer': 'https://www.douluo1.com/',
        'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.7'
    }
# 发送请求
    Resp = requests.get(url,headers= headers,verify=False)
# 设置编码
    Resp.encoding = "utf-8"
# 响应信息
# print(Resp.text)
    e = etree.HTML(Resp.text)
    info = '\n'.join(e.xpath("//div[@class =\"txtnav\"]/p/text()"))
    print(info)
    title = e.xpath("//h1/text()")[0]
    url = "https://www.douluo1.com{}".format(e.xpath("//div[@class='page1']/a[4]/@href")[0])
# 保存
    with open("斗罗大陆.txt","a",encoding='utf-8')as f:
        f.write(title+'\n\n'+info+'\n\n')
    if i * 25 >= 225:
        break
