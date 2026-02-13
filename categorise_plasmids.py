from opencloning.syntax import Syntax
import json
from pydna.parsers import parse
import os
import shutil


def merge_syntaxes(syntax_file1: str, syntax_file2: str) -> Syntax:
    syntax1 = Syntax.model_validate_json(open(syntax_file1).read())
    syntax2 = Syntax.model_validate_json(open(syntax_file2).read())

    syntax1.parts = syntax2.parts + syntax1.parts
    syntax1.overhangNames.update(syntax2.overhangNames)
    # Reassing ids of parts
    for i, part in enumerate(syntax1.parts):
        part.id = i + 1
    return syntax1


index_file = os.path.join("syntaxes", "index.json")
index = json.load(open(index_file))

for syntax_entry in index:
    # Load the syntax
    syntax_path = syntax_entry["path"]
    syntax_file = os.path.join("syntaxes", syntax_path, "syntax.json")
    if "uses_syntax" in syntax_entry:
        parent_syntax = os.path.join("syntaxes", syntax_entry["uses_syntax"], "syntax.json")
        # Copy the syntax.json file to the syntaxes directory
        shutil.copy(parent_syntax, syntax_file)
    elif "extends_syntax" in syntax_entry:
        parent_syntax = os.path.join("syntaxes", syntax_entry["extends_syntax"], "syntax.json")
        sub_syntax = os.path.join("syntaxes", syntax_path, "sub_syntax.json")
        # Merge the syntaxes
        merged_syntax = merge_syntaxes(sub_syntax, parent_syntax)
        with open(syntax_file, 'w') as f:
            f.write(merged_syntax.model_dump_json(indent=4))

    syntax = Syntax.model_validate_json(open(syntax_file).read())

    plasmids = list()
    for associated_kit in syntax_entry["kits"]:
        print('>', associated_kit["kit"])
        associated_kit: dict
        kit_path = os.path.join("kits", associated_kit["kit"])
        plasmid_names = associated_kit.get("names", None)
        associated_plasmids_tsv = os.path.join(kit_path, "plasmids.tsv")
        associated_plasmids_tsv_lines = open(associated_plasmids_tsv).readlines()
        for line in associated_plasmids_tsv_lines[1:]:
            ls = line.split("\t")
            well, name, addgene_id, resistance = ls[:4]
            content = '' if len(ls) < 5 else ls[4]
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

# Write a minified index.json file
mini_index = dict()
for syntax_entry in index:
    mini_index[syntax_entry["name"]] = {
        "path": syntax_entry["path"],
        "description": syntax_entry["description"]
    }

with open('syntaxes/index.min.json', 'w') as f:
    json.dump(mini_index, f, indent=4)
