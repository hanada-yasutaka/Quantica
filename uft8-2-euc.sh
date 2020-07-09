for name in *.m; do
    if [ `nkf -g ${name}` = "EUC-JP" ]; then
        tmp=${name/.m/_utf8.m}
        nkf -w ${name} > ${tmp}
        mv ${tmp} ${name}
        echo "convert ${name} to euc" 
    fi
done
