from collections import Counter
import pandas as pd

def add_domains(data, path_to_domains):
    domains = pd.read_csv(path_to_domains, delimiter=' ')
    data = data.merge(domains, right_on='proteinID', left_on='uniprotID', how='left')
    data = data.drop(['proteinID'], axis=1)
    # Label each data point as range or notRange based on the relative distance of mutation and domain boundaries.
    data = data.astype('str')
    data.domStart = data.domStart.astype('float')
    data.domEnd = data.domEnd.astype('float')

    for i in data.index:
        if data.at[i, 'domain'] != 'nan':
            if int(data.at[i, 'domStart']) <= int(data.at[i, 'pos']) <= int(data.at[i, 'domEnd']):
                data.at[i, 'distance'] = 0
            else:
                distance = min(abs(int(data.at[i, 'domStart']) - int(data.at[i, 'pos'])),
                               abs(int(data.at[i, 'domEnd']) - int(data.at[i, 'pos'])))
                data.at[i, 'distance'] = int(distance)
        else:
            data.at[i, 'distance'] = 'nan'

    data = data.sort_values(by=['datapoint', 'distance']).reset_index(drop=True)  # Distances will be sorted.

    # Keep the one with the least distance. But we may have more than one range domains for a datapoint if distance = 0.
    # For this reason first we need to separate range ones so that when we take the first occurance to get the closest one
    # for non range ones, other distance=0 ones wont disappear.

    data_range = data[data.distance == 0]
    data_out_range = data[data.distance != 0]

    # For the range ones, find the most occurance

    dom = []
    for i in data_range.index:
        dom.append(data_range.at[i, 'domain'])

    domainCount = Counter(dom)  # Occurance of domains.

    # For out of range ones, take the closest distance.
    data_out_range = data_out_range.drop_duplicates(['datapoint'], keep='first')  # Already sorted above.
    domain_counts = pd.DataFrame(domainCount.items(), columns=['domain', 'count'])
    data_range_counts = data_range.merge(domain_counts, on='domain')
    data_range_counts = data_range_counts.sort_values(['datapoint', 'count'])
    data_range_counts = data_range_counts.drop_duplicates(['datapoint'], keep='last')  # Take with the higher count.
    data_range_counts = data_range_counts.drop(['count'], axis=1)

    # Merge them back together

    frames = [data_range_counts, data_out_range]
    data = pd.concat(frames, sort=False)  # Here when you concat two data frames, we might have range and not range with
    # min distance for the same data point. Delete the one coming from notRange one.
    data = data.sort_values(['datapoint', 'distance']).reset_index(drop=True)
    data = data.drop_duplicates(['datapoint'], keep='first')
    data = data.astype(str)
    return data
