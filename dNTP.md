# dNTP triphosphohydrolase
1. 使用hmmsearch最大范围内统计不同菌株中蛋白的copy数  
2. 使用hmmscan保留最显著的domain为给定famliy的蛋白，并统计不同菌株中蛋白的copy数据,e值为1e-50  
3. 使用blastp将种子序列在所有蛋白里面检索直至没有新的结果出现  


+ 检索：使用hmmsearch抓取相同domain的基因或者基因家族  

```bash
cd ~/data/Pseudomonas
cat STRAINS/Pseudom_aeruginosa_PAO1/*.tsv | grep "IPR006261"
# 数据库登录号：发现只在tiger数据库检索到了
#  TIGRFAM TIGR01353 
# NP_249815.1     94ee528dd44d70777e0ed12875c07c17        498     TIGRFAM TIGR01353       dGTP_triPase: putative dGTPase     36      496     1.1E-107        T       15-02-2023      IPR006261       dNTP triphosphohydrolase
# NP_251733.1     26e9713ec4adde923c798dd755da30f3        443     TIGRFAM TIGR01353       dGTP_triPase: putative dGTPase     24      436     4.7E-98         T       15-02-2023      IPR006261       dNTP triphosphohydrolase
# IPR006261       dNTP triphosphohydrolase        2
# NP_249815.1     IPR006261       dNTP triphosphohydrolase
# NP_251733.1     IPR006261       dNTP triphosphohydrolase


cd ~/Shenwei/data
mkdir -p  Pseudomonas_xrz/dNTP/HMM
mkdir -p  ~/Shenwei/data/Pseudomonas_xrz/summary/
cp ~/Shenwei/data/Pseudomonas/dNTP/HMM/* ~/Shenwei/data/Pseudomonas_xrz/dNTP/HMM/
cp ~/Shenwei/data/Pseudomonas/summary/* ~/Shenwei/data/Pseudomonas_xrz/summary/
ln -s ~/Shenwei/data/Pseudomonas/ASSEMBLY  ~/Shenwei/data/Pseudomonas_xrz/ASSEMBLY
# 对比师兄的total.lst与师姐的taxon/${GENUS}中的STRAIN
# total.lst -> ~/Shenwei/data/Pseudomonas_xrz/summary/total.lst
# 
# ${GENUS}属 
# Acinetobacter
# Azotobacter
# Bordetella
# Burkholderia
# Pseudomonas
# Serratia
# Stenotrophomonas
# Stutzerimonas


#在蛋白库里搜索该基因的结构域
cd ~/Shenwei/data/Pseudomonas_xrz
E_VALUE=1e-20
for domain in TIGR01353;do
    >&2 echo "==> domain [${domain}]"

    cat summary/total.lst | 
        parallel --no-run-if-empty --linebuffer -k -j 24 "
            gzip -dcf ASSEMBLY/{}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw dNTP/HMM/${domain}.hmm - |
                grep '>>' |
                    perl -nl -e '
                        m{>>\s+(\S+)} or next;
                        \$n = \$1;
                        \$s = \$n;
                        \$s =~ s/\.\d+//;
                        printf qq{%s\t%s_%s\n}, \$n, {}, \$s;
                        '
        " > dNTP/${domain}.replace.tsv

        # -E ${E_VALUE}在结果中报告E值小于这个阈值的序列
        # --domE ${E_VALUE} 在结果中报告E值小于这个阈值的domain

done
#  4736 TIGR01353.replace.tsv
mv dNTP/TIGR01353.replace.tsv dNTP/dNTP.replace.tsv 


# 看一下输出文件
E_VALUE=1e-20
for domain in TIGR01353;do
    >&2 echo "==> domain [${domain}]"

    cat summary/total.lst | 
        parallel --no-run-if-empty --linebuffer -k -j 24 "
            gzip -dcf ASSEMBLY/{}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE}  --noali --notextw dNTP/HMM/${domain}.hmm - >> output.txt" 
        # -E ${E_VALUE}在结果中报告E值小于这个阈值的序列
        # --domE ${E_VALUE} 在结果中报告E值小于这个阈值的domain
done


# tsv-join brkb/PF03631.replace.tsv \
#     -f brkb/PTHR30213.replace.tsv |
#     tsv-join -f brkb/TIGR00765.replace.tsv
#     > brkb/brkb.replace.tsv
# # 6637

# 统计结果
cp -r ~/Shenwei/data/Pseudomonas/script ~/Shenwei/data/Pseudomonas_xrz/script
perl script/make_table.pl -t summary/strains.taxon.tsv -i dNTP/TIGR01353.replace.tsv -a summary/total.lst > dNTP/statistics.tsv


# 根据ncbi的annotation查看有哪些结构域
ln -s ~/data/Pseudomonas/PROTEINS ~/Shenwei/data/Pseudomonas_xrz/PROTEINS
cat dNTP/dNTP.replace.tsv | tsv-select -f 2,1 | tsv-join -d 1 \
-f PROTEINS/all.annotation.tsv -k 1 --append-fields 2 | tsv-summarize -g 3 --count
# deoxyguanosinetriphosphate triphosphohydrolase  3735 铜绿
# deoxyguanosinetriphosphate triphosphohydrolase family protein   17
# anti-phage deoxyguanosine triphosphatase        4
# dNTP triphosphohydrolase, partial       1
# dNTP triphosphohydrolase        78 部分铜绿
# hypothetical protein    2
# dGTPase 891 铜绿
# HD domain-containing protein    1
# 标注为抗噬菌体酶的物种：Acin_bau，Acin_pit，Pseudom_fragi，Ser_ply拷贝数均为1

# 整理结构域结果，不符合的被剔除
cat dNTP/dNTP.replace.tsv | tsv-select -f 2,1 | tsv-join -d 1 \
-f PROTEINS/all.annotation.tsv -k 1 --append-fields 2 | tsv-filter --str-not-in-fld 3:"hypothetical protein" | tsv-filter --str-not-in-fld 3:"HD domain-containing protein" > dNTP/dNTP.filter.replace.tsv
cat dNTP/dNTP.filter.replace.tsv | tsv-select -f 2,1  > dNTP/dNTP.filter.summary

cd ~/Shenwei/data/Pseudomonas_xrz
perl script/make_table.pl -t summary/strains.taxon.tsv -i dNTP/dNTP.filter.summary -a summary/total.lst > dNTP/filter.statistics.tsv
```

## 初步建树
+ Protein tree全部菌株建立蛋白树

```bash
mkdir -p protein_tree

pigz -dcf PROTEINS/all.replace.fa.gz |
    faops some stdin <(tsv-select -f 2 dNTP/dNTP.filter.summary | grep "Pseudom_aeru") protein_tree/dNTP.pseudom.fa

mkdir -p ~/Shenwei/data/Pseudomonas_xrz/dNTP/reference/
cp  ~/Shenwei/data/Pseudomonas/reference/dNTP/* ~/Shenwei/data/Pseudomonas_xrz/dNTP/reference/
cat dNTP/reference/TIGR01353.replace.tsv | grep "Acin_pit" > protein_tree/dNTP.outgroup.tsv


# 设置外类群
for i in $(sed 's/\t/-/g' protein_tree/dNTP.outgroup.tsv);do
    PROTEIN=$(echo $i | cut -d "-" -f 1)
    REPLACE=$(echo $i | cut -d "-" -f 2)
    SPECIES=$(echo $i | cut -d "-" -f 2 | rev | tr "_" "\t" | tsv-select --exclude 1,2 | rev | tr "\t" "_")
# rev是reverse，tr将_替换为空格

    gzip -dcf ASSEMBLY/${SPECIES}/*_protein.faa.gz |
        faops some stdin <(echo $PROTEIN) stdout |
        faops replace stdin <(echo $i | tr "-" "\t") stdout
done \
    > protein_tree/dNTP.outgroup.fa

# 找到铜绿假单胞菌的所有菌株建树
cd ~/Shenwei/data/Pseudomonas_xrz/
(cat protein_tree/dNTP.pseudom.fa && cat protein_tree/dNTP.outgroup.fa) > protein_tree/dNTP.tree.fa

muscle -in protein_tree/dNTP.tree.fa -out protein_tree/dNTP.tree.aln.fa

bsub -q serial -n 15 -J "iqtree" -o protein_tree/  "iqtree -s protein_tree/dNTP.tree.aln.fa -m MFP -bb 1000 -alrt 1000 -nt AUTO"

# 重新定根
nw_reroot protein_tree/dNTP.tree.aln.fa.treefile $(nw_labels protein_tree/dNTP.tree.aln.fa.treefile | grep -E "Acin_pit") |
    nw_order -c n - \
    > protein_tree/dNTP.tree.reoot.treefile
```

+ ref_tree
```bash
E_VALUE=1e-20

# 对模式物种菌株建树，10/15中找到了
# make tree
for domain in TIGR01353;do
    >&2 echo "==> domain [${domain}]"

    cat summary/reference.lst | 
        parallel --no-run-if-empty --linebuffer -k -j 24 "
            gzip -dcf ASSEMBLY/{}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw dNTP/HMM/${domain}.hmm - |
                grep '>>' |
                    perl -nl -e '
                        m{>>\s+(\S+)} or next;
                        \$n = \$1;
                        \$s = \$n;
                        \$s =~ s/\.\d+//;
                        printf qq{%s\t%s_%s\n}, \$n, {}, \$s;
                        '
        " > dNTP/reference/${domain}.replace.tsv
done

# tsv-join reference/PF03631.replace.tsv \
#     -f reference/PTHR30213.replace.tsv | 
# tsv-join -f reference/TIGR00765.replace.tsv > reference/brkb.replace.tsv

for i in $(sed 's/\t/-/g' dNTP/reference/TIGR01353.replace.tsv);do
    PROTEIN=$(echo $i | cut -d "-" -f 1)
    REPLACE=$(echo $i | cut -d "-" -f 2)
    SPECIES=$(echo $i | cut -d "-" -f 2 | rev | tr "_" "\t" | tsv-select --exclude 1,2 | rev | tr "\t" "_")
    
    gzip -dcf ASSEMBLY/${SPECIES}/*_protein.faa.gz |
        faops some stdin <(echo $PROTEIN) stdout |
        faops replace stdin <(echo $i | tr "-" "\t") stdout
done \
    > dNTP/reference/dNTP.fa

muscle -in dNTP/reference/dNTP.fa -out dNTP/reference/dNTP.aln.fa


bsub -q serial -n 15 -J "iqtree" -o dNTP/reference/ "iqtree -s dNTP/reference/dNTP.aln.fa -m MFP -bb 1000 -alrt 1000 -nt AUTO"
```


+ 分析模式物种共线性，还没做
```bash
PAO1(NC_002516.2):PA2751(3112878-3113777,899),PA0951(1037650-1038885,1235)
PA7(NC_009656.1):PSPA7_RS11930(2505736-2506635,899),PSPA7_RS21670(4680813-4682048,1235)
LESB58(NC_011770.1):PLES_RS12070(2489873-2490772,899),PLES_RS22480(4799500-4800735,1235)
UCBPP_PA14(NC_008463.1):PA14_RS11605(2460195-2461094,899),PA14_RS21155(4612649-4613884,1235)

P.anta(NZ_LT629704.1):BLQ27_RS25350(5529747-5530616,869),BLQ27_RS10310(2253379-2254335,956)
P.kna(NZ_CP047698.1):GVI47_RS08125(1710205-1711077,872),GVI47_RS06420(1370800-1372029,1229)
P.citro(NZ_CP014158.1):PcP3B5_RS20975(4836045-4836953,908),PcP3B5_RS08355(1837059-1838288,1229)
P.oti(NZ_AP022642.1):PotMrB4_RS04125(872834-873730,896),PotMrB4_RS21425(4640586-4641806,1220)

P.alcalig(NZ_CP014784.1):,A0T30_RS16330(3552696-3553913,1217)
P.mend(NZ_CP013124.1):,DW68_RS08615(1825499-1826716,1217)
P.lal(NZ_CP043311.1):,FXN65_RS21535(4663689-4664909,1220)

Azot.chroo1(NZ_CP011835.1):ACG10_RS12850(2804478-2805347,869),ACG10_RS14595(3152001-3153227,1226)
Azot.chroo(CP011835.1):ACG10_12555(2804478-2805347,869),ACG10_14255(3152001-3153227,1226)
Azot.vine(NC_012560.1):AVIN_RS14965(3370019-3370888,869),AVIN_RS16780(3738396-3739622,1226)
Stu.stut(NZ_AP024722.1):PszF2a_RS06855(1409687-1410562,875),PszF2a_RS14595(3090034-3091266,1232)
Cellv.japonicus(NZ_CP043306.1):,FY115_RS11055(2617501-2618784,1283)
```

```bash
# reference
# Clade 1
clinker gff/P.aeru.PAO1.gff gff/P.aeru.PA7.gff gff/P.aeru.LESB58.gff gff/P.aeru.UCBPP_PA14.gff -r NC_002516.2:3102878-3123777 NC_009656.1:2495736-2516635 NC_011770.1:2479873-2500772 NC_008463.1:2450195-2471094 -p test.html
# Clade 2
clinker gff/P.aeru.PAO1.gff gff/P.aeru.PA7.gff gff/P.aeru.LESB58.gff gff/P.aeru.UCBPP_PA14.gff -r NC_002516.2:1027650-1048885 NC_009656.1:4670813-4692048 NC_011770.1:4789500-4810735 NC_008463.1:4602649-4623884 -p test.html

# Clade 1 and Clade 2
python script/fetch_gbk.py -i genebank/PAO1.gbff -l PA2751 -e 10000 -o brkb
python script/fetch_gbk.py -i genebank/PAO1.gbff -l PA0951 -e 10000 -o brkb
clinker brkb/*.gbk -p test.html

# Pseudom 
# Clade 1
clinker gff/P.aeru.PAO1.gff gff/P.kna.gff gff/P.citro.gff gff/P.otitidis.gff gff/Azot.vine.gff gff/Stu.stut.gff -r NC_002516.2:3102878-3123777 NZ_CP047698.1:1700205-1721077 NZ_CP014158.1:4826045-4846953 NZ_AP022642.1:862834-883730 NC_012560.1:3360019-3380888 NZ_AP024722.1:1399687-1420562 -p test.html
# Clade 2
python script/fetch_gbk.py -i genebank/mend.gbff -l DW68_RS08615 -e 10000 -o brkb
clinker gff/P.aeru.PAO1.gff gff/P.kna.gff gff/P.citro.gff gff/P.alcalig.gff gff/P.lal.gff gff/P.otitidis.gff gff/Azot.vine.gff gff/Stu.stut.gff gff/Cellv.japonicus.gff brkb/P.mend.gbk -r NC_002516.2:1027650-1048885 NZ_CP047698.1:1360800-1382029 NZ_CP014158.1:1827059-1848288 NZ_CP014784.1:3542696-3563913 NZ_CP043311.1:4653689-4674909 NZ_AP022642.1:4630586-4651806 NC_012560.1:3728396-3749622 NZ_AP024722.1:3080034-3101266 NZ_CP043306.1:2607501-2628784 -p test.html
```

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
# HD_domain       18 要保留！要不然铜绿没有了


#拼接属名等信息并统计拷贝数
cat dNTP/dNTP_minevalue.tsv | tsv-select -f 1,3 |
tsv-summarize -g 2 --count 
# Phosphohydrolase-associated_domain      1756
# HD_domain       18

cat dNTP/dNTP_minevalue.tsv | tsv-select -f 1,3  >dNTP/dNTP.hmmscan_filter.replace.tsv
cat dNTP/dNTP.hmmscan_filter.replace.tsv | tsv-select -f 1 | grep -Eo '([^_]+_[^_]+)$' | sed 's/$/.1/' > dNTP/dNTP.hmmscan_filter_WP.replace.tsv
paste dNTP/dNTP.hmmscan_filter_WP.replace.tsv dNTP/dNTP.hmmscan_filter.replace.tsv > dNTP/dNTP.hmmscan_filter.summary
# 形成统计表格
perl script/make_table.pl -t summary/strains.taxon.tsv -i dNTP/dNTP.hmmscan_filter.summary -a summary/total.lst > dNTP/hmmscan_filter.statistics.tsv




# 因为筛选出来的太少，使用dNTP/dNTP.abstract.tsv代替dNTP/dNTP_minevalue.tsv重新处理一遍
#拼接属名等信息并统计拷贝数
cat dNTP/dNTP.abstract.tsv | tsv-select -f 3,6 | tsv-summarize -g 2 --count 
# Phosphohydrolase-associated_domain      1756
# HD_domain       1564
# Probable_N6-adenine_methyltransferase   1

cat dNTP/dNTP.abstract.tsv | tsv-select -f 3,6 | tsv-filter --str-not-in-fld 2:"Probable_N6-adenine_methyltransferase" >dNTP/dNTP.hmmscan_filter.replace2.tsv
cat dNTP/dNTP.hmmscan_filter.replace2.tsv | tsv-select -f 1 | grep -Eo '([^_]+_[^_]+)$' | sed 's/$/.1/' > dNTP/dNTP.hmmscan_filter_WP.replace2.tsv
paste dNTP/dNTP.hmmscan_filter_WP.replace2.tsv dNTP/dNTP.hmmscan_filter.replace2.tsv > dNTP/dNTP.hmmscan_filter2.summary


perl script/make_table.pl -t summary/strains.taxon.tsv -i dNTP/dNTP.hmmscan_filter2.summary -a summary/total.lst > dNTP/hmmscan_filter.statistics2.tsv

cat dNTP/hmmscan_filter.statistics2.tsv | awk '{sum += $4} END {print sum}'
# 3318
```

## 2. 使用tigrfams数据库，该蛋白只有TIGER中找到了蛋白数据，有可能只能在TIGER下查找
```bash
#下载tigerfam数据库
mkdir -p ~/data/HMM/TIGERFAM
cd ~/data/HMM/TIGERFAM
wget -N --content-disposition https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM.tgz
tar -xzvf hmm_PGAP.HMM.tgz
cat *.HMM >tigrfams.hmm
ln -s ~/shiyuqi/pseudomonas/TIGERFAM ~/Shenwei/data/Pseudomonas_xrz/HMM/

#格式化tigerfam数据库 
bsub -q mpi -n 24 -J "hmm" "hmmpress ~/Shenwei/data/Pseudomonas_xrz/HMM/TIGERFAM/tigerfam.hmm"


# 将提取的dNTP蛋白序列与数据库比对
cd ~/Shenwei/data/Pseudomonas_xrz/PROTEINS
gzip -dc all.replace.fa.gz > all.replace.fa
cd ~/Shenwei/data/Pseudomonas_xrz/
faops some PROTEINS/all.replace.fa <(tsv-select -f 1 dNTP/dNTP.filter.replace.tsv)  dNTP/dNTP.tiger.fa

E_value=1e-40
NAME=dNTP
bsub -q mpi -n 24 -J "hmmscan" -o dNTP/ " hmmscan --cpu 12 -E 1e-40 --domE 1e-40 --noali --tblout ${NAME}/${NAME}.tiger.tbl  HMM/TIGERFAM/tigerfam.hmm  ${NAME}/${NAME}.tiger.fa"

perl script/abstract.pl dNTP/dNTP.tiger.tbl > dNTP/dNTP_tiger.abstract.tsv
# 提取其中关键信息：target name，accession，query name，full sequence E-value，best 1 domain E-value，description of target

#查看dNTP同家族蛋白的数据库登录号以及结构域描述
cat  dNTP/dNTP_tiger.abstract.tsv | tsv-summarize -g 2,6  --count
# NF002205.0      NCBI_Protein_Cluster_(PRK):_deoxyguanosinetriphosphate_triphosphohydrolase      3222
# TIGR01353.1     JCVI:_dNTP_triphosphohydrolase  4726
# NF003429.0      NCBI_Protein_Cluster_(PRK):_dGTPase     3218
# NF003701.0      NCBI_Protein_Cluster_(PRK):_dGTPase     4704
# NF002829.0      NCBI_Protein_Cluster_(PRK):_deoxyguanosinetriphosphate_triphosphohydrolase      3826
# NF002326.0      NCBI_Protein_Cluster_(PRK):_deoxyguanosinetriphosphate_triphosphohydrolase      1567
# NF002328.0      NCBI_Protein_Cluster_(PRK):_deoxyguanosinetriphosphate_triphosphohydrolase      1567
# NF002327.0      NCBI_Protein_Cluster_(PRK):_deoxyguanosinetriphosphate_triphosphohydrolase      1567
# NF002330.0      NCBI_Protein_Cluster_(PRK):_deoxyguanosinetriphosphate_triphosphohydrolase      1567
# NF002329.0      NCBI_Protein_Cluster_(PRK):_deoxyguanosinetriphosphate_triphosphohydrolase      1508

#查看domain以及domain的描述
cat  dNTP/dNTP_tiger.abstract.tsv | tsv-summarize -g 1,6 --count


#同一蛋白序列可以匹配到多条model序列,只保留e值最小的且description符合该famliy的菌株蛋白序列名
perl script/compare.pl dNTP/dNTP_tiger.abstract.tsv >dNTP/dNTP_tiger_minevalue.tsv
tsv-summarize -g 3 --count dNTP/dNTP_tiger_minevalue.tsv
# NCBI_Protein_Cluster_(PRK):_deoxyguanosinetriphosphate_triphosphohydrolase      3763
# NCBI_Protein_Cluster_(PRK):_dGTPase     920 
# JCVI:_dNTP_triphosphohydrolase  43

#不筛选domain拼接属名等信息并统计拷贝数
cat dNTP/dNTP_tiger_minevalue.tsv | tsv-select -f 1,3 |
tsv-summarize -g 2 --count 
cat dNTP/dNTP_tiger_minevalue.tsv | tsv-select -f 1,3  >dNTP/dNTP.hmmscan_tiger_filter.replace.tsv
cat dNTP/dNTP.hmmscan_tiger_filter.replace.tsv | tsv-select -f 1 | grep -Eo '([^_]+_[^_]+)$' | sed 's/$/.1/' > dNTP/dNTP.hmmscan_tiger_filter_WP.replace.tsv
paste dNTP/dNTP.hmmscan_tiger_filter_WP.replace.tsv dNTP/dNTP.hmmscan_tiger_filter.replace.tsv > dNTP/dNTP.hmmscan_tiger_filter.summary
# 形成统计表格
perl script/make_table.pl -t summary/strains.taxon.tsv -i dNTP/dNTP.hmmscan_tiger_filter.summary -a summary/total.lst > dNTP/hmmscan_tiger_filter.statistics.tsv


# 去掉最后一个统计一下，因为PAO1含有前两个domain命名
cat dNTP/dNTP_tiger_minevalue.tsv | tsv-select -f 1,3 | tsv-filter --str-not-in-fld 2:"JCVI:_dNTP_triphosphohydrolase" >dNTP/dNTP.hmmscan_tiger_filter_new.replace.tsv
cat dNTP/dNTP.hmmscan_tiger_filter_new.replace.tsv | tsv-select -f 1 | grep -Eo '([^_]+_[^_]+)$' | sed 's/$/.1/' > dNTP/dNTP.hmmscan_tiger_filter_WP_new.replace.tsv
paste dNTP/dNTP.hmmscan_tiger_filter_WP_new.replace.tsv dNTP/dNTP.hmmscan_tiger_filter_new.replace.tsv > dNTP/dNTP.hmmscan_tiger_filter_new.summary
# 形成统计表格
perl script/make_table.pl -t summary/strains.taxon.tsv -i dNTP/dNTP.hmmscan_tiger_filter_new.summary -a summary/total.lst > dNTP/hmmscan_tiger_filter_new.statistics.tsv
```
!!! 结果和前面有些不同，有些物种两拷贝呗排除在外

# 使用blastp筛选一遍hmmscan结果蛋白
* 使用上面排除了最后一个domain的结果
```bash
cd ~/Shenwei/data/Pseudomonas_xrz
#提取hmmsearch和hmmscan结果
cat dNTP/dNTP.hmmscan_tiger_filter_new.replace.tsv | cut -f 1 > dNTP/dNTP_diamond1.tsv

#第一轮diamond
faops some PROTEINS/all.replace.fa  dNTP/dNTP_diamond1.tsv dNTP/dNTP_diamond1.fa
diamond makedb --in dNTP/dNTP_diamond1.fa --db dNTP/dNTP
bsub -q mpi -n 24 -J "blastp" -o dNTP/dNTP "
diamond blastp --db dNTP/dNTP.dmnd --query PROTEINS/all.replace.fa -e 1e-10 --outfmt 6 --threads 16 --out dNTP/dNTP_blastp_result1.tsv"

#第二轮diamond
faops some PROTEINS/all.replace.fa <(cat dNTP/dNTP_blastp_result1.tsv | cut -f 2 | sort -n | uniq) dNTP/dNTP_diamond2.fa
diamond makedb --in dNTP/dNTP_diamond2.fa --db dNTP/dNTP
bsub -q mpi -n 24 -J "blastp" -o dNTP/dNTP "
diamond blastp --db dNTP/dNTP.dmnd --query PROTEINS/all.replace.fa -e 1e-10 --outfmt 6 --threads 16 --out dNTP/dNTP_blastp_result2.tsv"


#hmmer结果
cat dNTP/dNTP_diamond1.tsv | wc -l  #

#第一轮diamond的query
cut -f 1 dNTP/dNTP_blastp_result1.tsv | sort -n | uniq | wc -l #
#第一轮diamond的target
cut -f 2 dNTP/dNTP_blastp_result1.tsv | sort -n | uniq | wc -l #

#第二轮diamond的query
cut -f 1 dNTP/dNTP_blastp_result2.tsv | sort -n | uniq | wc -l #
#第二轮diamond的target
cut -f 2 dNTP/dNTP_blastp_result2.tsv | sort -n | uniq | wc -l #
```


























+ Pseudom tree
```bash
E_VALUE=1e-20

# make tree
for domain in TIGR00765 PF03631 PTHR30213;do
    >&2 echo "==> domain [${domain}]"

    cat Pseudom.rep.lst | 
        parallel --no-run-if-empty --linebuffer -k -j 24 "
            gzip -dcf ASSEMBLY/{}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw brkb/hmm/${domain}.hmm - |
                grep '>>' |
                    perl -nl -e '
                        m{>>\s+(\S+)} or next;
                        \$n = \$1;
                        \$s = \$n;
                        \$s =~ s/\.\d+//;
                        printf qq{%s\t%s_%s\n}, \$n, {}, \$s;
                        '
        " > pseudom_tree/${domain}.replace.tsv
done

tsv-join pseudom_tree/PF03631.replace.tsv \
    -f pseudom_tree/PTHR30213.replace.tsv | 
tsv-join -f pseudom_tree/TIGR00765.replace.tsv > pseudom_tree/brkb.replace.tsv

pigz -dcf PROTEINS/all.replace.fa.gz |
    faops some stdin <(tsv-select -f 2 pseudom_tree/brkb.replace.tsv) pseudom_tree/brkb.fa

muscle -in pseudom_tree/brkb.fa -out pseudom_tree/brkb.aln.fa

iqtree -s pseudom_tree/brkb.aln.fa -m MFP -bb 1000 -alrt 1000 -nt AUTO

#nw_reroot pseudom_tree/brkb.aln.fa.treefile $(nw_labels pseudom_tree/brkb.aln.fa.treefile | grep -E "Acin_bau") |
#    nw_order -c n - \
#    > pseudom_tree/brkb.reoot.treefile
```

+ gamma tree
```bash
E_VALUE=1e-20

# make tree
for domain in TIGR00765 PF03631 PTHR30213;do
    >&2 echo "==> domain [${domain}]"

    cat gamma.rep.lst | 
        parallel --no-run-if-empty --linebuffer -k -j 24 "
            gzip -dcf ASSEMBLY/{}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw brkb/hmm/${domain}.hmm - |
                grep '>>' |
                    perl -nl -e '
                        m{>>\s+(\S+)} or next;
                        \$n = \$1;
                        \$s = \$n;
                        \$s =~ s/\.\d+//;
                        printf qq{%s\t%s_%s\n}, \$n, {}, \$s;
                        '
        " > gamma_tree/${domain}.replace.tsv
done

tsv-join gamma_tree/PF03631.replace.tsv \
    -f gamma_tree/PTHR30213.replace.tsv | 
tsv-join -f gamma_tree/TIGR00765.replace.tsv > gamma_tree/brkb.replace.tsv

for i in $(sed 's/\t/-/g' gamma_tree/brkb.replace.tsv);do
    PROTEIN=$(echo $i | cut -d "-" -f 1)
    REPLACE=$(echo $i | cut -d "-" -f 2)
    SPECIES=$(echo $i | cut -d "-" -f 2 | rev | tr "_" "\t" | tsv-select --exclude 1,2 | rev | tr "\t" "_")
    
    gzip -dcf ASSEMBLY/${SPECIES}/*_protein.faa.gz |
        faops some stdin <(echo $PROTEIN) stdout |
        faops replace stdin <(echo $i | tr "-" "\t") stdout
done \
    > gamma_tree/brkb.fa

muscle -in gamma_tree/brkb.fa -out gamma_tree/brkb.aln.fa

iqtree -s gamma_tree/brkb.aln.fa -m MFP -bb 1000 -alrt 1000 -nt AUTO

nw_reroot gamma_tree/brkb.aln.fa.treefile $(nw_labels gamma_tree/brkb.aln.fa.treefile | grep -E "Sta_aure") |
nw_order -c n - \
    > gamma_tree/brkb.reoot.treefile

# statistics info
for domain in TIGR00765 PF03631 PTHR30213;do
    >&2 echo "==> domain [${domain}]"

    cat gamma.rep.lst | 
        parallel --no-run-if-empty --linebuffer -k -j 24 "
            echo -e '{}'
            gzip -dcf ASSEMBLY/{}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw brkb/hmm/${domain}.hmm - |
                grep '>>'
        " > gamma_tree/${domain}_info.tsv
    
    perl script/statistics_info.pl gamma_tree/${domain}_info.tsv > tem&&
    mv tem gamma_tree/${domain}_info.tsv
done

tsv-join gamma_tree/PF03631_info.tsv \
    -f gamma_tree/PTHR30213_info.tsv | 
tsv-join -f gamma_tree/TIGR00765_info.tsv > gamma_tree/brkb_info.tsv
```
