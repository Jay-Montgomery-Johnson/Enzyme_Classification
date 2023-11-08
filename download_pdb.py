'''
Search for PDB files with the desired characteristics and download the pdb files
'''
import requests, json, os
from Bio.PDB import PDBList

#pdb query built using hte query GUI
rows = 100
pdb_ids = {}
for EC in range(1, 8):
    pdb_ids[EC] =  []
    
    query = '''{
  "query": {
    "type": "group",
    "logical_operator": "and",
    "nodes": [
      {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "attribute": "exptl.method",
          "operator": "exact_match",
          "negation": false,
          "value": "X-RAY DIFFRACTION"
        }
      },
      {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "attribute": "rcsb_entry_info.diffrn_resolution_high.value",
          "operator": "less_or_equal",
          "negation": false,
          "value": 2.4
        }
      },
      {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "attribute": "rcsb_polymer_entity.rcsb_ec_lineage.id",
          "operator": "in",
          "negation": false,
          "value": [
            "''' + str(EC) + '''"
          ]
        }
      },
      {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "attribute": "rcsb_entry_info.selected_polymer_entity_types",
          "operator": "exact_match",
          "negation": false,
          "value": "Protein (only)"
        }
      },
      {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "attribute": "rcsb_entry_info.nonpolymer_entity_count",
          "operator": "equals",
          "negation": false,
          "value": 0
        }
      },
      {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "attribute": "rcsb_entry_info.deposited_polymer_monomer_count",
          "operator": "less_or_equal",
          "negation": false,
          "value": 1000
        }
      }
    ],
    "label": "text"
  },
  "return_type": "entry",
  "request_options": {
    "paginate": {
      "start": 0,
      "rows": ''' + str(rows) + '''
    },
    "results_content_type": [
      "experimental"
    ],
    "sort": [
      {
        "sort_by": "score",
        "direction": "desc"
      }
    ],
    "scoring_strategy": "combined"
  }
}'''
    
    response = requests.get(f'https://search.rcsb.org/rcsbsearch/v2/query?json={query}')
    data = json.loads(response.text)
    
    for p in data['result_set']:
        pdb_ids[EC].append(p['identifier'])

#record the pdb ids I use in a text file
with open('dataset/pdb_files_in_database.txt', 'w') as file:
     file.write(json.dumps(pdb_ids))

#download the pdb files
pdbl = PDBList()
for k in pdb_ids.keys():
    for i in pdb_ids[k]:
        pdbl.retrieve_pdb_file(i, pdir='dataset', file_format='pdb')
        os.rename(f'dataset/pdb{str(i).lower()}.ent', f'dataset/ec_{k}/pdb/{i}.pdb')
