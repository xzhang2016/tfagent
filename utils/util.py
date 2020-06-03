import urllib.request
import os

def merge_dict_sum(dict1, dict2):
    """
    Merge two dictionaries and add values of common keys.
    Values of the input dicts can be any addable objects, like numeric, str, list.
    """
    dict3 = {**dict1, **dict2}
    for key, value in dict3.items():
        if key in dict1 and key in dict2:
            dict3[key] = value + dict1[key]
    return dict3
    
def merge_dict_list(dict1, dict2):
    """
    Merge two dictionaries and merge values of common keys to a list.
    """
    dict3 = {**dict1, **dict2}
    for key, value in dict3.items():
        if key in dict1 and key in dict2:
            dict3[key] = [value, dict1[key]]
    return dict3

def download_file_dropbox(url, fout_name):
    """
    Download file from dropbox
    """
    u = urllib.request.urlopen(url)
    data = u.read()
    u.close()
    
    #save
    with open(fout_name, 'wb') as fw:
        fw.write(data)
        
def make_folder(data_folder):
    if not os.path.isfile(data_folder):
        # Emulate mkdir -p (no error if folder exists)
        try:
            os.mkdir(data_folder)
            return True
        except Exception:
            return False
