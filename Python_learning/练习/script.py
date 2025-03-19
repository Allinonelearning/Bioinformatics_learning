from selenium import webdriver
import requests
from lxml import etree
import pandas as pd

# 设置 Selenium WebDriver
driver = webdriver.Chrome()

# 访问网页
url ="https://book.douban.com/top250"
i = 1
books_info = []
while True:
# 等待页面加载完成
    driver.get(url)
    driver.execute_script("window.scrollTo(0, document.body.scrollHeight);")
    driver.implicitly_wait(20)
# 获取页面源代码
    page_source = driver.page_source
# 解析页面
    e = etree.HTML(page_source)
# 执行 XPath 查询
    Book_list = e.xpath("//div[@class=\'indent\']/table/tbody/tr/td[2]/div/a")
    Book_info = e.xpath("//div[@class = \"indent\"]/table/tbody/tr/td[2]/p[1]/text()")
    Book_score = e.xpath("//div[@class = \"indent\"]/table/tbody/tr/td[2]/div[2]/span[2]/text()")
    # Book_famous= e.xpath("//div[@class = \"indent\"]/table/tbody/tr/td[2]/p[2]/span/text()")
# 打印结果
# 打印书名
#     for book, info,score,famous in zip(Book_list, Book_info,Book_score,Book_famous):
#         print(f"书名: {book.text.strip()}, 作品信息: {info.strip()},得分：{score.strip()},评论：{famous.strip()}")
    url = "https://book.douban.com/top250?start="+str(i*25)
    # print(url)
    i = i+1
    for book, info, score in zip(Book_list, Book_info, Book_score):
            books_info.append({"书名": book.text.strip(),
            "作品信息": info.strip(),
            "得分": score.strip()})
    if i * 25 >= 275:
        break
driver.quit()# 关闭浏览器
# 将书籍信息列表转换为 DataFrame
df_books = pd.DataFrame(books_info)
# 保存 DataFrame 为 CSV 文件
# df_books.to_csv("books_info.csv", index=False, encoding="utf-8-sig")
# # df_books.to_excel("books_info.xlsx", index=False)
df_books.describe()
