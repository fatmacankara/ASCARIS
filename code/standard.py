def standardize(df, get_columns):
    cols_to_change = ['sasa', 'domaindistance3D', 'disulfide', 'intMet', 'intramembrane',
                      'naturalVariant', 'dnaBinding', 'activeSite', 'nucleotideBinding',
                      'lipidation', 'site', 'transmembrane', 'crosslink', 'mutagenesis',
                      'strand', 'helix', 'turn', 'metalBinding', 'repeat', 'caBinding',
                      'topologicalDomain', 'bindingSite', 'region', 'signalPeptide',
                      'modifiedResidue', 'zincFinger', 'motif', 'coiledCoil', 'peptide',
                      'transitPeptide', 'glycosylation', 'propeptide']
    for col in cols_to_change:  # because in the other ones, they are 3D distance. Here, no distance calculated.
        df[col] = 'nan'
    df = df[get_columns.columns]

    return df
