import pandas as pd


import prointvar.library


def variant_report_table(variant_table, path=None):
    """Format a variant table for a report"""
    report = variant_table.copy()
    report.columns = [e[-1] if isinstance(e, tuple) else e for e in report.columns]
    _whitelist = ['ID', 'HGVSp', 'AC', 'AC_Male', 'AC_POPMAX', 'POPMAX']
    report = report[_whitelist].dropna(subset=['HGVSp'])
    report['ID'] = report['ID'].fillna('-')
    report.loc[:, _whitelist[2:-1]] = (report[_whitelist[2:-1]].fillna(0).astype(int))
    report['HGVSp'] = report['HGVSp'].str.replace('.*:', '')
    report = report.rename(columns={'HGVSp': 'Mutation',
                                    'AC': 'Allele Count',
                                    'AC_Male': 'AC Male',
                                    'AC_POPMAX': 'AC Popn. Max',
                                    'POPMAX': 'Popn. Max'}
                          )
    if path:
        report.to_csv(path, index=False)
    return report


def categorise_mutagenesis_assays(effects):
    """Parse UniProt mutagenesis annotation effect to categorical"""
    map_ = {'Abolishes interaction with SARS-CoV spike glycoprotein': 'Abolishes interaction',
            'Strongly inhibits interaction with SARS-CoV spike glycoprotein': 'Strongly inhibits',
            'Inhibits interaction with SARS-CoV spike glycoprotein': 'Inhibits',
            'Slightly inhibits interaction with SARS-CoV spike glycoprotein': 'Slightly inhibits',
            'No effect on interaction with SARS-CoV spike glycoprotein': 'No effect'
           }
    cats = effects.map(map_).astype('category')
    cats = cats.cat.add_categories([c for c in map_.values() if c not in cats.cat.categories]).cat.reorder_categories(map_.values(), ordered=True)
    return cats


def read_MutaBind2_table(path):
    return pd.read_csv(path, skiprows=1)


def format_mCSM_mutation(mCSM_table):
    return (mCSM_table['wild-type'].map(prointvar.library.aa_codes_3to1_common)
            + mCSM_table['res-number'].astype('str')
            + mCSM_table['mutant'].map(prointvar.library.aa_codes_3to1_common))


def read_mCSM_table(path):
    """Read mCSM-PPI2 prediction file and add mutation column"""
    t = pd.read_csv(path)
    t['mutation'] = format_mCSM_mutation(t)
    return t

