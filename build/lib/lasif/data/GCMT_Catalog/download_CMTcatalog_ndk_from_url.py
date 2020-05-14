import requests
import os

parent_url = 'https://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/NEW_MONTHLY/'
year = 2016

month_list = ['jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec']

for month in month_list:
    year_str = str(year)
    
    ndk_file = '%s%s.ndk'%(month,year_str[-2:])
    child_url = os.path.join(parent_url, year_str, ndk_file)
    print(("Requesting CMT catalog at %s"%child_url))
    
    child_folder = year_str
    
    try:
        r = requests.get(child_url)
        
        if not os.path.exists(child_folder):
            os.makedirs(child_folder)
               
        if r.status_code ==200:
            with open(os.path.join(child_folder, ndk_file), 'w') as f:
                f.write(r.content)
    except:
        continue
        
    
    