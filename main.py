import pandas as pd
import requests as re
import argparse
from tqdm import tqdm

parser = argparse.ArgumentParser()

parser.add_argument('-input', help='Path to table')
parser.add_argument('-col_name', help='Column name with MEROPS families in table')
parser.add_argument('--fm', '-families_names', help="Enter MEROPS families names like: S8 A2 S56", nargs='+')
parser.add_argument('output')
parser.add_argument('type', type=int, choices=[1, 2],
                    help='1 - use MEROPS names from command line, 2 - take names from data frame')

args = parser.parse_args()

nested_list = []
pre_tabs = pd.DataFrame()
df1 = pd.DataFrame()

fam = args.fm

def pad_dict_list(dict_list, padel):
    """
    https://stackoverflow.com/questions/40442014/python-pandas-valueerror-arrays-must-be-all-same-length/56757839#56757839
    by robertspierre

    :param dict_list:
    :param padel:
    :return:
    """

    lmax = 0
    for lname in dict_list.keys():
        lmax = max(lmax, len(dict_list[lname]))
    for lname in dict_list.keys():
        ll = len(dict_list[lname])
        if  ll < lmax:
            dict_list[lname] += [padel] * (lmax - ll)
    return dict_list


def downloader_lst(familis):
    nested_fam_list = []
    for fm in tqdm(familis, colour='green'):
        try:
            url = f'https://www.ebi.ac.uk/merops/cgi-bin/famsum?family={fm}'
            holder = re.get(url)
            if holder:
                pre_data = pd.read_html(holder.text)[0]
                func_data = pd.read_html(holder.text)[2] # Так как строка про функции в другой таблице
                pre_lst = []
                list_data = pre_data.values.tolist()
                func_list = func_data.values.tolist()

                for i,j in zip(list_data, func_list):
                    if 'Name' in i:
                        list_data.remove(i)
                    if 'Active site residues' in i:
                        d = i
                        list_data.remove(i)
                        list_data.append(d)
                    if 'Biological functions' in j:
                        list_data.append(j)

                for j in list_data:
                    pre_lst.append(j[1])
                nested_fam_list.append([str(f'{fm} -- {url}'), pre_lst])

            else:
                nested_fam_list.append('')
        except ImportError:
            print(f'{fm} -- ERROR')

    pre_dict = dict(nested_fam_list)
    pre_dict['Name'] = ['Family type peptidase',
                         'Content of family',
                         'History',
                         'Catalytic type',
                         'Active site',
                         'Activities and specificities',
                         'Inhibitors',
                         'Molecular structure',
                         'Clan',
                         'Basis of clan assignment',
                         'Biological functions',
                         'Active site residues',]

    diction = pad_dict_list(pre_dict, padel='')
    df1 = pd.DataFrame.from_dict(diction)
    df1 = df1.set_index('Name')
    df1_transpose = df1.transpose()
    return df1_transpose


def downloader_tab(table):
    familis = []
    ids = []
    for idx, row in table.iterrows():
        if 'subfamily' in row[f'{args.col_name}'] or 'family' in row[f'{args.col_name}']:
            familis.append(str(row[f'{args.col_name}']).split(' ')[1])
            ids.append(row['Gene id'])

    df2 = downloader_lst(familis)
    return df2


if __name__ == '__main__':

    if args.type == 1:
        downloader_lst(fam).to_csv(fr'{args.output}')

    elif args.type == 2:
        data_frame = pd.read_csv(rf'{args.input}')
        downloader_tab(data_frame).to_csv(fr'{args.output}')
