## Please download the MiClip package the put this script to the same folder
library("MiClip")
mc = MiClip("data.sam")
rd = MiClip.read(mc)
en = MiClip.enriched(rd)
bd = MiClip.binding(en)
write.table(bd$enriched, "bd_enriched.txt", sep="\t")
write.table(bd$sites, "bd_sites.txt", sep="\t")
write.table(bd$clusters, "bd_clusters.txt", sep="\t")
