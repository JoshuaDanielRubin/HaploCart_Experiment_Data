cat ../src/pangenome_mapping | grep + | cut -d, -f1,4 | sed -E 's/("([^"]*)")?,/\2\t/g' >> ../src/parsed_pangenome_mapping


