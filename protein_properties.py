'''
Determines the properties of proteins based on aaindex tables
'''
from aaindex import aaindex1 as aa

class aa_index:
    def __init__(self):
        self.hydrophobicity = aa['ARGP820101'].values
        self.polarisability = aa['CHAM820101'].values
        self.vdw_volume = aa['FAUJ880103'].values
        self.polarity = aa['GRAR740102'].values
        
        self.codes = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C',
                    'GLU':'E', 'GLN':'Q', 'GLY':'G', 'HIS':'H', 'HYP':'O',
                    'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F',
                    'PRO':'P', 'GLP':'U', 'SER':'S', 'THR':'T', 'TRP':'W',
                    'TYR':'Y', 'VAL':'V'}
        
    def get_properties(self, sequence):#type should be list
        total_hydrophobicity, total_polarisability, total_vdw_volume, total_polarity = 0, 0, 0, 0
        for a in sequence:
            a = self.codes.get(a, '-')
            total_hydrophobicity += self.hydrophobicity[a]
            total_polarisability += self.polarisability[a]
            total_vdw_volume += self.vdw_volume[a]
            total_polarity += self.polarity[a]
        seq_len = len(sequence)
        norm_hydrophobicity = total_hydrophobicity/seq_len 
        norm_polarisability = total_polarisability/seq_len
        norm_vdw_volume = total_vdw_volume/seq_len
        norm_polarity = total_polarity/seq_len
        return norm_hydrophobicity, norm_polarisability, norm_vdw_volume, norm_polarity
