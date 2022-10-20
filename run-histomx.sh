HISTOMX_PATH=~/Desktop/histomx
HISTOMX=${HISTOMX_PATH}/bin/render-histomx_report

RMD_FILE=${HISTOMX_PATH}/scripts/histomx_kidney.Rmd
RMD_FILE=${HISTOMX_PATH}/scripts/histomx_kidney_bkv.Rmd
RMD_FILE=${HISTOMX_PATH}/scripts/histomx_kidney_simple.Rmd

mkdir -p $HISTOMX_PATH/refRCCs
#symlink to ref RCCs
ln -s ~/Dropbox/PTG/transcriptomics/nanostring/kidney/refset/current/RefSet-KTD1-FFPE $HISTOMX_PATH/refRCCs

$HISTOMX -m ${RMD_FILE} -f ${HISTOMX_PATH}/test_files/test.RCC -i 'test' -p ${HISTOMX_PATH}/test_files/patient-test.json -r ${HISTOMX_PATH}/test_files/rna-test.json
$HISTOMX -m ${RMD_FILE} -f ${HISTOMX_PATH}/test_files/test.RCC -i 'test' -p ${HISTOMX_PATH}/test_files/patient-test.txt -r ${HISTOMX_PATH}/test_files/rna-test.txt
$HISTOMX -m ${RMD_FILE} -f ${HISTOMX_PATH}/test_files/amr.RCC -i 'amr-test'
$HISTOMX -m ${RMD_FILE} -f ${HISTOMX_PATH}/test_files/tcmr.RCC -i 'tcmr-test'


# LUNG
HISTOMX_PATH=~/Desktop/histomx-lung
HISTOMX=${HISTOMX_PATH}/bin/render-histomx_report

$HISTOMX -m ${HISTOMX_PATH}/scripts/histomx_lung_simple.Rmd -f ${HISTOMX_PATH}/test_files/tcmr.RCC -i 'test' -p ${HISTOMX_PATH}/test_files/patient-test.txt -r ${HISTOMX_PATH}/test_files/rna-test.txt

##--------------
## NS-344: Active AMR: g3, ptc3, v3, i0, t0, C4d3
$HISTOMX -m ${RMD_FILE} -f ${HISTOMX_PATH}/input_files/20210318_RefSet-run-25_NS-344_04.RCC -i 'amr-test' \
-p ${HISTOMX_PATH}/input_files/patient_info_NS344.txt -r ${HISTOMX_PATH}/input_files/rna_info_NS344.txt

## NS-141: Active AMR: g3, cg0, ptc3, v3, i0, t0, C4d3
$HISTOMX -m ${RMD_FILE} -f ${HISTOMX_PATH}/input_files/20200505_RefSet-run-8_NS-141_10.RCC -i 'NS-141-amr' \
-p ${HISTOMX_PATH}/input_files/patient_info_NS141.txt

## NS-85: Active AMR: g3, cg0, ptc3, i1, t0, C4d0
$HISTOMX -m ${RMD_FILE} -f ${HISTOMX_PATH}/input_files/20191231_RefSet-run-4_NS-85_05.RCC -i 'NS-85-amr' \
-p ${HISTOMX_PATH}/input_files/patient_info_NS85.txt -r ${HISTOMX_PATH}/input_files/rna_info_NS85.txt

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

## xeno
RMD_FILE=${HISTOMX_PATH}/scripts/histomx_xeno.Rmd
DIR=~/Dropbox/PTG/transcriptomics/nanostring/kidney/xeno/RCCs/

$HISTOMX -m ${RMD_FILE} -f ${DIR}/20211224_RefSet-Run-70_NP-01_01.RCC -i "NP-01" -r ${DIR}/rna_info_NP-01.txt -p ${DIR}/patient_info_NP-01.txt
$HISTOMX -m ${RMD_FILE} -f ${DIR}/20211224_RefSet-Run-70_NP-02_02.RCC -i "NP-02" -r ${DIR}/rna_info_NP-02.txt -p ${DIR}/patient_info_NP-02.txt
$HISTOMX -m ${RMD_FILE} -f ${DIR}/20220301_Nano-Pig-run1_NP-03_07.RCC -i "NP-03" -r ${DIR}/rna_info_NP-03.txt -p ${DIR}/patient_info_NP-03.txt
$HISTOMX -m ${RMD_FILE} -f ${DIR}/20220301_Nano-Pig-run1_NP-07_11.RCC -i "NP-07" -p ${DIR}/patient_info_NP-07.txt
$HISTOMX -m ${RMD_FILE} -f ${DIR}/20220301_Nano-Pig-run1_NP-08_12.RCC -i "NP-08" -r ${DIR}/rna_info_NP-08.txt -p ${DIR}/patient_info_NP-08.txt
$HISTOMX -m ${RMD_FILE} -f ${DIR}/20220422_KTD2-RNAlater-Run11_NP-09_10.RCC -i "NP-09" -r ${DIR}/rna_info_NP-09.txt -p ${DIR}/patient_info_NP-09.txt

for file in $(ls ${DIR}/*html);do
    wkhtmltopdf --no-stop-slow-scripts $file ${DIR}/$(basename $file ".html").pdf
done

open ${DIR}/RCCs/*pdf

## pediatric
HISTOMX_PATH=~/Desktop/histomx
HISTOMX=${HISTOMX_PATH}/bin/render-histomx_report

RMD_FILE=${HISTOMX_PATH}/scripts/histomx_kidney_bkv.Rmd

DIR=~/Dropbox/PTG/transcriptomics/nanostring/kidney/pediatric/RCC/
$HISTOMX -m ${RMD_FILE} -f ${DIR}/*/*/20200618_208901300821_HJU12290-0010-S15-138_09.RCC -i "S15-138"
$HISTOMX -m ${RMD_FILE} -f ${DIR}/*/*/20200620_208901370821_HJU12290-0029-S17-48_03.RCC -i "S17-48"
$HISTOMX -m ${RMD_FILE} -f ${DIR}/*/20210602_RefSet-run-49_NS-634_10.RCC -i "NS-634"

wkhtmltopdf /Users/dina.zielinski/Dropbox/PTG/transcriptomics/nanostring/kidney/pediatric/RCC/round1-emory/20200618_208901300821_RCC/S15-138_histomx.html ~/Desktop/ped/S15-138_histomx.pdf
wkhtmltopdf /Users/dina.zielinski/Dropbox/PTG/transcriptomics/nanostring/kidney/pediatric/RCC/round1-emory/20200620_208901370821_RCC/S17-48_histomx.html ~/Desktop/ped/S17-48_histomx.pdf
wkhtmltopdf /Users/dina.zielinski/Dropbox/PTG/transcriptomics/nanostring/kidney/pediatric/RCC/round2-borderline/NS-634_histomx.html ~/Desktop/ped/NS-634_histomx.pdf

# CEDARS
# v1.7
HISTOMX_PATH=~/Desktop/histomx-2b834cc75d33680cc0697bc1e0d768d1236344b4
HISTOMX=${HISTOMX_PATH}/bin/render-histomx_report
RMD_FILE=${HISTOMX_PATH}/scripts/histomx_kidney.Rmd
DIR=~/Dropbox/PTG/transcriptomics/nanostring/kidney/cedars

$HISTOMX -m ${RMD_FILE} -f $DIR/RCC/20220713_RefSet-runCedars01-Kidneys_NS-384C_06.RCC -p $DIR/sample_info_NS-384.txt
$HISTOMX -m ${RMD_FILE} -f $DIR/RCC/20220713_RefSet-runCedars01-Kidneys_NS-391C_09.RCC -p $DIR/sample_info_NS-391.txt
$HISTOMX -m ${RMD_FILE} -f $DIR/RCC/20220713_RefSet-runCedars01-Kidneys_NS-394C_10.RCC -p $DIR/sample_info_NS-394.txt
$HISTOMX -m ${RMD_FILE} -f $DIR/RCC/20220713_RefSet-runCedars01-Kidneys_NS-403C_12.RCC -p $DIR/sample_info_NS-403.txt

$HISTOMX -m ${RMD_FILE} -f $DIR/RCC/20220713_RefSet-runCedars01-Kidneys_NS-156C_01.RCC
$HISTOMX -m ${RMD_FILE} -f $DIR/RCC/20220713_RefSet-runCedars01-Kidneys_NS-336C_02.RCC
$HISTOMX -m ${RMD_FILE} -f $DIR/RCC/20220713_RefSet-runCedars01-Kidneys_NS-344C_03.RCC
$HISTOMX -m ${RMD_FILE} -f $DIR/RCC/20220713_RefSet-runCedars01-Kidneys_NS-364C_05.RCC
$HISTOMX -m ${RMD_FILE} -f $DIR//RCC/20220713_RefSet-runCedars01-Kidneys_NS-387C_07.RCC
$HISTOMX -m ${RMD_FILE} -f $DIR/RCC/20220713_RefSet-runCedars01-Kidneys_NS-389C_08.RCC

for file in $(ls ${DIR}/RCC/*RCC);do
    ID=${file%.RCC}
    #$HISTOMX -m ${RMD_FILE} -f ${file}
    wkhtmltopdf --no-stop-slow-scripts ${DIR}/RCC/$(basename $file ".RCC")_histomx.html ${DIR}/$(basename $file ".RCC").pdf
done

## clone specific commit
git checkout <commit_sha>


