from opencloning.syntax import Syntax
import json
from pydna.parsers import parse
import os

index_file = os.path.join("syntaxes", "index.json")
index = json.load(open(index_file))

for syntax_entry in index:
    # Load the syntax
    syntax_path = syntax_entry["path"]
    if "uses_syntax" in syntax_entry:
        syntax_file = os.path.join("syntaxes", syntax_entry["uses_syntax"], "syntax.json")
    else:
        syntax_file = os.path.join("syntaxes", syntax_path, "syntax.json")

    syntax = Syntax.model_validate_json(open(syntax_file).read())

    plasmids = list()
    for associated_kit in syntax_entry["kits"]:
        associated_kit: dict
        kit_path = os.path.join("kits", associated_kit["kit"])
        plasmid_names = associated_kit.get("names", None)
        associated_plasmids_tsv = os.path.join(kit_path, "plasmids.tsv")
        associated_plasmids_tsv_lines = open(associated_plasmids_tsv).readlines()
        for line in associated_plasmids_tsv_lines[1:]:
            well, name, addgene_id, resistance, content = line.split("\t")
            if plasmid_names is not None and name not in plasmid_names:
                continue
            content = content.strip()
            seq = parse(f'addgene_plasmids/{addgene_id}.gb')[0]
            resp = syntax.assign_plasmid_to_syntax_part(seq)
            if len(resp) == 0:
                print(f'Skipped plasmid {name} ({addgene_id}) because it could not be assigned to a syntax part')
                continue
            elif len(resp) > 1:
                print(f'Skipped plasmid {name} ({addgene_id}) because it could be assigned to multiple syntax parts')
                continue
            else:
                resp = resp[0]

            left_overhang, right_overhang = resp['key'].split('-')
            longest_feature_name = ''
            if resp['longest_feature'] is not None and 'label' in resp['longest_feature'].qualifiers:
                longest_feature_name = resp['longest_feature'].qualifiers['label'][0]
            if content != '':
                name += f" ({content})"

            plasmid = {
                'id': len(plasmids) + 1,
                'plasmid_name': name,
                'left_overhang': left_overhang,
                'right_overhang': right_overhang,
                'key': resp['key'],
                'type': 'AddgeneIdSource',
                'source': {
                    'id': 1,
                    'type': 'AddgeneIdSource',
                    'input': [],
                    'repository_id': addgene_id,
                },
            }
            plasmids.append(plasmid)

    if not os.path.exists(f'syntaxes/{syntax_path}'):
        os.makedirs(f'syntaxes/{syntax_path}')
    with open(f'syntaxes/{syntax_path}/plasmids.json', 'w') as f:
        json.dump(plasmids, f, indent=4)
