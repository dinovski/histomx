HISTOMX_PATH=~/Desktop/histomx
HISTOMX=${HISTOMX_PATH}/bin/render-histomx_report

RMD_FILE=${HISTOMX_PATH}/scripts/histomx_kidney.Rmd
RMD_FILE=${HISTOMX_PATH}/scripts/histomx_kidney_full.Rmd

mkdir -p $HISTOMX_PATH/refRCCs
#symlink to ref RCCs
#ln -s ~/Dropbox/PTG/transcriptomics/nanostring/kidney/refset/current/RefSet-KTD1-FFPE $HISTOMX_PATH/refRCCs

$HISTOMX -m ${RMD_FILE} -f ${HISTOMX_PATH}/test_files/test.RCC -i 'test' -p ${HISTOMX_PATH}/test_files/patient-test.txt -r ${HISTOMX_PATH}/test_files/rna-test.txt

$HISTOMX -m ${RMD_FILE} -f ${HISTOMX_PATH}/test_files/test.RCC -i 'test' -p ${HISTOMX_PATH}/test_files/empty.json -r ${HISTOMX_PATH}/test_files/rna-test.json
$HISTOMX -m ${RMD_FILE} -f ${HISTOMX_PATH}/test_files/test.RCC -i 'test' -p ${HISTOMX_PATH}/test_files/empty.json -r ${HISTOMX_PATH}/test_files/rna-test.txt
$HISTOMX -m ${RMD_FILE} -f ${HISTOMX_PATH}/test_files/test.RCC -i 'test' -p ${HISTOMX_PATH}/test_files/patient-test.json -r ${HISTOMX_PATH}/test_files/rna-test.json
$HISTOMX -m ${RMD_FILE} -f ${HISTOMX_PATH}/test_files/amr.RCC -i 'amr-test'
$HISTOMX -m ${RMD_FILE} -f ${HISTOMX_PATH}/test_files/tcmr.RCC -i 'tcmr-test'
$HISTOMX -m ${RMD_FILE} -f ~/Dropbox/PTG/transcriptomics/nanostring/kidney/beta-testing/KTD1-RNAlater-run-3/RCC/20201215_KTD1-RNAlater-run-3_NG-1738_05.RCC -i "NG-1738" -p ~/Desktop/patient-test.txt

$HISTOMX -m ${RMD_FILE} -f ${HISTOMX_PATH}/test_files/20210401_RefSet-run-30_NS-413_11.RCC -i "NS-413" -p ${HISTOMX_PATH}/test_files/ns-413-patient.txt 
$HISTOMX -m ${RMD_FILE} -f ${HISTOMX_PATH}/test_files/20210804_RefSet-run-60_NS-748_12.RCC -i "NS-748" -p ${HISTOMX_PATH}/test_files/ns-748-patient.txt 

##--------------
## NS-344: Active AMR: g3, ptc3, v3, i0, t0, C4d3
$HISTOMX -m ${RMD_FILE} -f ${HISTOMX_PATH}/input_files/20210318_RefSet-run-25_NS-344_04.RCC -i 'amr-test' \
-p ${HISTOMX_PATH}/input_files/patient_info_NS344.txt -r ${HISTOMX_PATH}/input_files/rna_info_NS344.txt

DIR=~/Desktop/histomx/new_run/
for file in $(ls ${DIR}/*RCC);do
    ID=${file%.RCC}
    $HISTOMX -m ${RMD_FILE} -f ${file} -i '${ID}'
    #wkhtmltopdf --no-stop-slow-scripts ${DIR}/$(basename $file ".RCC")_histomx.html ${DIR}/$(basename $file "html")pdf
done

RCC=${DIR}/20220125_KTD2-RNAlater-run-8_NG-2304_12.RCC
$HISTOMX -m ${RMD_FILE} -f ${RCC} -i 'NG-2304'

file=${DIR}/20220125_KTD2-RNAlater-run-8_NG-2294_02_histomx.html
wkhtmltopdf --load-error-handling ignore --no-stop-slow-scripts --enable-local-file-access --encoding utf8 --javascript-delay 10 $(basename $file) ${DIR}/$(basename $file "html")pdf

# v1.7
HISTOMX_PATH=~/Desktop/histomx-2b834cc75d33680cc0697bc1e0d768d1236344b4
HISTOMX=${HISTOMX_PATH}/bin/render-histomx_report
RMD_FILE=${HISTOMX_PATH}/scripts/histomx_kidney.Rmd
DIR=~/Dropbox/PTG/transcriptomics/nanostring/kidney/cedars

$HISTOMX -m ${RMD_FILE} -f $DIR/RCC/20220713_RefSet-runCedars01-Kidneys_NS-384C_06.RCC -p $DIR/sample_info_NS-384.txt
$HISTOMX -m ${RMD_FILE} -f $DIR/RCC/20220713_RefSet-runCedars01-Kidneys_NS-391C_09.RCC -p $DIR/sample_info_NS-391.txt
$HISTOMX -m ${RMD_FILE} -f $DIR/RCC/20220713_RefSet-runCedars01-Kidneys_NS-394C_10.RCC -p $DIR/sample_info_NS-394.txt
$HISTOMX -m ${RMD_FILE} -f $DIR/RCC/20220713_RefSet-runCedars01-Kidneys_NS-403C_12.RCC -p $DIR/sample_info_NS-403.txt

for file in $(ls ${DIR}/RCC/*RCC);do
    ID=${file%.RCC}
    #$HISTOMX -m ${RMD_FILE} -f ${file}
    wkhtmltopdf --no-stop-slow-scripts ${DIR}/RCC/$(basename $file ".RCC")_histomx.html ${DIR}/$(basename $file ".RCC").pdf
done

## clone specific commit
git checkout <commit_sha>

## xeno
IDIR=~/Dropbox/PTG/transcriptomics/nanostring/kidney/xeno/RCCs/

$HISTOMX -m ${RMD_FILE} -f ${IDIR}/20211224_RefSet-Run-70_NP-01_01.RCC -p ${IDIR}/patient_info_NP-01.txt -r ${IDIR}/rna_info_NP-01.txt
$HISTOMX -m ${RMD_FILE} -f ${IDIR}/20211224_RefSet-Run-70_NP-02_02.RCC -p ${IDIR}/patient_info_NP-02.txt -r ${IDIR}/rna_info_NP-02.txt

$HISTOMX -m ${RMD_FILE} -f ${IDIR}/20220301_Nano-Pig-run1_NP-03_07.RCC -p ${IDIR}/patient_info_NP-03.txt -r ${IDIR}/rna_info_NP-03.txt
$HISTOMX -m ${RMD_FILE} -f ${IDIR}/20220301_Nano-Pig-run1_NP-07_11.RCC -p ${IDIR}/patient_info_NP-07.txt
$HISTOMX -m ${RMD_FILE} -f ${IDIR}/20220301_Nano-Pig-run1_NP-08_12.RCC -p ${IDIR}/patient_info_NP-08.txt -r ${IDIR}/rna_info_NP-08.txt
$HISTOMX -m ${RMD_FILE} -f ${IDIR}/20220422_KTD2-RNAlater-Run11_NP-09_10.RCC -p ${IDIR}/patient_info_NP-09.txt -r ${IDIR}/rna_info_NP-09.txt
