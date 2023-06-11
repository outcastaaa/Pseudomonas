# 使用hmmscan重新筛选一遍蛋白,在本地跑会更快，后期可以在本地跑
## 1. 使用Pfam数据库
```bash
#下载pfam数据库
mkdir -p ~/data/HMM/PFAM
cd  ~/data/HMM/PFAM
for basename in Pfam-A.hmm Pfam-A.hmm.dat active_site.dat; do
    wget -N --content-disposition ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/${basename}.gz
    wget -N --content-disposition ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/${basename}.gz
    wget -N --content-disposition ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/${basename}.gz
done

for basename in Pfam-A.hmm Pfam-A.hmm.dat active_site.dat; do
    echo "==> ${basename}"
    gzip -dcf ${basename}.gz > ${basename}
done

# 直接链接过来
mkdir -p ~/Shenwei/data/Pseudomonas_xrz/HMM/PFAM
ln -s ~/data/HMM/PFAM ~/Shenwei/data/Pseudomonas_xrz/HMM/
#格式化pfam数据库:Pfam 数据库中每个编号代表一个蛋白质家族。Pfam 分 A 和 B 两个数据库，其中 A 数据库是经过手工校正的高质量数据库
hmmpress ~/Shenwei/data/Pseudomonas_xrz/HMM/PFAM/Pfam-A.hmm


# 将提取的dNTP蛋白序列与pfam数据库比对
cd ~/Shenwei/data/Pseudomonas_xrz/PROTEINS
gzip -dc all.replace.fa.gz > all.replace.fa
faops some PROTEINS/all.replace.fa <(tsv-select -f 1 dNTP/dNTP.filter.replace.tsv)  dNTP/dNTP.pfam.fa

E_value=1e-10
NAME=dNTP
bsub -q mpi -n 24 -J "hmmscan" -o dNTP/ " hmmscan --cpu 12 -E 1e-10 --domE 1e-10 --noali --tblout ${NAME}/${NAME}.tbl  HMM/PFAM/Pfam-A.hmm  dNTP/dNTP.pfam.fa"

cp ~/shiyuqi/pseudomonas/abstract.pl ~/Shenwei/data/Pseudomonas_xrz/script/
perl script/abstract.pl dNTP/dNTP.tbl > dNTP/dNTP.abstract.tsv
# 提取其中关键信息：target name，accession，query name，full sequence E-value，best 1 domain E-value，description of target

#查看dNTP同家族蛋白的数据库登录号以及结构域描述
cat  dNTP/dNTP.abstract.tsv | tsv-summarize -g 2,6  --count
# PF13286.6       Phosphohydrolase-associated_domain      1756
# PF01966.22      HD_domain       1564
# PF10237.9       Probable_N6-adenine_methyltransferase   1

#查看domain以及domain的描述
cat  dNTP/dNTP.abstract.tsv | tsv-summarize -g 1,6 --count
# HD_assoc        Phosphohydrolase-associated_domain      1756
# HD      HD_domain       1564
# N6-adenineMlase Probable_N6-adenine_methyltransferase   1

#同一蛋白序列可以匹配到多条model序列,只保留e值最小的且description符合该famliy的菌株蛋白序列名
cp ~/shiyuqi/pseudomonas/compare.pl ~/Shenwei/data/Pseudomonas_xrz/script/
perl script/compare.pl dNTP/dNTP.abstract.tsv >dNTP/dNTP_minevalue.tsv
tsv-summarize -g 3 --count dNTP/dNTP_minevalue.tsv
# Phosphohydrolase-associated_domain      1756
# HD_domain       18 后续删除即可


#拼接属名等信息并统计拷贝数
cat dNTP/dNTP_minevalue.tsv | tsv-select -f 1,3 |
tsv-filter --str-in-fld 2:"Phosphohydrolase-associated_domain" | \
tsv-summarize -g 2 --count 
# Phosphohydrolase-associated_domain      1756
cat dNTP/dNTP_minevalue.tsv | tsv-select -f 1,3 |
tsv-filter --str-in-fld 2:"Phosphohydrolase-associated_domain" >dNTP/dNTP.hmmscan_filter.replace.tsv

cat dNTP/dNTP.hmmscan_filter.replace.tsv | tsv-select -f 1 | grep -Eo '([^_]+_[^_]+)$' | sed 's/$/.1/' > dNTP/dNTP.hmmscan_filter_WP.replace.tsv

paste dNTP/dNTP.hmmscan_filter_WP.replace.tsv dNTP/dNTP.hmmscan_filter.replace.tsv > dNTP/dNTP.hmmscan_filter.summary


perl script/make_table.pl -t summary/strains.taxon.tsv -i dNTP/dNTP.hmmscan_filter.summary -a summary/total.lst > dNTP/hmmscan_filter.statistics.tsv

```
