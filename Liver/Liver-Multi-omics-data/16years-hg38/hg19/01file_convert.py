###bigBed文件转为bed文件###
# 1、下载小工具：http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/

# 2、Linux安装小程序：（1）安装教程：https://www.cnblogs.com/foreverycc/archive/2013/05/17/3083032.html
                                   （2）cd+程序所在位置；chmod +x bigWigToBedGraph；添加路径，vim ~/.bashrc；按“i”插入路径；按“Esc”退出，并按“：wq”保存生效；按“:q”退出
                                   （3）source ~/.bashrc激活
# 3、运行命令：bigWigToBedGraph H3K27ac.bigWig H3K27ac.BedGraph


###组学数据由hg38转为hg19###
# 1、下载坐标转换对应文件，hg38到hg19
         hg19到hg38网址：http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
                                      wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
         gunzip hg38ToHg19.over.chain.gz 
        #下载转换工具
           wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
           chmod 755 liftOver
        #查看帮助
           ./liftOver
# 2、#格式：liftOver oldFile map.chain newFile unMapped
        #  直接执行，若不在一个文件夹，用全路径，如：/home/user/tools/liftOver
        ./liftOver GSM3715311_CTCF_CT_peaks.bed hg38ToHg19.over.chain GSM3715311_CTCF_CT_peaks_hg19.bed unMapped


##hg38转换为hg19后进行排序