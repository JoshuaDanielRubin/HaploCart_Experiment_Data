
for NODE in $(seq 1 1 11821)
do
odgi position -i ../src/graph_dir/graph.odgi -g $NODE -r "H2a2a1" >> ../data/pangenome_mapping
echo "Done with"
echo $BASE
done
