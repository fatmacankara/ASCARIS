import ssl
import requests as r

def get_pdb_ids(protein_id):
    # Fetch PDB IDs associated with given UniProtID
    ssl._create_default_https_context = ssl._create_unverified_context
    url = 'https://www.uniprot.org/uploadlists/'
    params_ = {
        'from': 'ACC+ID',
        'to': 'PDB_ID',
        'format': 'tab',
        'query': protein_id.strip()
    }

    response = r.get(url, params=params_)
    response = response.text.split('\n')
    response = list(filter(None, response))

    pdbs = {}
    pdbs_per_protein = []
    for i in range(len(response)):
        try:
            if response[i].split('\t')[0] in pdbs.keys():
                pdbs_per_protein.append(response[i].split('\t')[1])
            elif response[i].split('\t')[0] not in pdbs.keys():
                pdbs_per_protein = []
                pdbs_per_protein.append(response[i].split('\t')[1])
        except:
            IndexError

        pdbs[response[i].split('\t')[0]] = pdbs_per_protein
    try:
        pdbs.pop('From')
    except:
        KeyError
        print(pdbs)
    return pdbs