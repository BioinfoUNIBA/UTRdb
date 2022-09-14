import os                                                   
import sys                                    
import pandas as pd                          
from bs4 import BeautifulSoup

path = 'Ires_table.html'
data = []
list_header = []                                            
soup = BeautifulSoup(open(path),'html.parser')
header = soup.find_all("table")[0].find("tr")

for items in header:                                        
    try:                                      
        list_header.append(items.get_text()) 
    except:                  
        continue     

HTML_data = soup.find_all("table")[0].find_all("tr")[1:] 
for element in HTML_data:                                   
    sub_data = []                             
    for sub_element in element:              
        try:                 
            sub_data.append(sub_element.get_text())
        except:
            continue
    data.append(sub_data)

dataFrame = pd.DataFrame(data = data, columns = list_header)

dataFrame.to_csv('Ires_table.csv', index=False) 
