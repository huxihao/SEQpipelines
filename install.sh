echo "address http://159.226.118.212" > server.txt
echo "home home/"$USER >> server.txt
mkdir "/var/www/html/home/"$USER
ln -s /data/pipelines/transidx single-cell-rna-seq/
ln -s /data/pipelines/softwares/bin single-cell-rna-seq/
ln -s $PWD"/single-cell-rna-seq/work" "/var/www/html/home/"$USER"/single-cell-rna-seq"
ln -s /data/pipelines/transidx rna-seq/
ln -s /data/pipelines/softwares/bin rna-seq/
ln -s $PWD"/rna-seq/work" "/var/www/html/home/"$USER"/rna-seq"
#ln -s /data/pipelines/transidx rna-crosslink/
ln -s /data/pipelines/softwares/bin rna-crosslink/
ln -s $PWD"/rna-crosslink/work" "/var/www/html/home/"$USER"/rna-crosslink"
