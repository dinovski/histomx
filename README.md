**Generate HistoMx report for a single RCC file**

The following R packages must be installed:  
"archetypes", "DESeq2", "dplyr", "ggplot2", "ggrepel", "gridExtra", "htmltools", "kableExtra", "knitr", "MASS", "plyr", "prettydoc", "pROC", "rmarkdown", "RUVSeq", "RCRnorm", "stringi", "stringr", "xml2", "yaml"

1. Open the 'terminal' application on your computer.

2. Change to wherever you want to use as your working directory:
```
cd ~/Desktop
```

> **OPTION 1 (without git):**
> * Click the 'Code' button at the top of this repository and then 'Download ZIP'
> * Unzip and copy this directory ('histomx-main') to your working directory (eg. ~/Desktop) and rename it to 'histomx'
>
> **OPTION 2 (with git version control):**
> * Install git by entering the following in the terminal window:
> ```
> brew install git
> ```
> Other installation options: https://www.atlassian.com/fr/git/tutorials/install-git
>
> * Type the following command in order to clone this repository (ie. copy everything to your working directory):
> ```
> git clone https://github.com/dinovski/histomx.git
> ```

3. Create a directory called 'refRCCs'.
```
mkdir -p $HISTOMX_PATH/refRCCs
```

4. Download all reference RCC files:  
https://drive.google.com/drive/folders/1Wzi9LCof7QMcYyx7kLOiKWuBFEY2o8Zk?usp=sharing

5. Move all RCC files to this directory:
```
mv ~/Downloads/RefSet-KTD1-FFPE/* refRCCs/
```

6. Set variables to define paths to scripts/files
```
HISTOMX_PATH=~/Desktop/histomx 

HISTOMX=${HISTOMX_PATH}/bin/render-histomx_report  

RMD_FILE=${HISTOMX_PATH}/scripts/histomx_report.Rmd
```

7. Run histomx to generate the report
* This will output an .html file (you can specify the name of the file with the '-i' argument).
* The -m (or --rmd) and -f (or --rcc) arguments are required

Running the HISTOMX command will show all possible input arguments:
```
$HISTOMX
```

The following are 2 examples (with and without patient or RNA QC files)
```
$HISTOMX -m ${RMD_FILE} -f ${HISTOMX_PATH}/test_files/test.RCC -i 'test-run'
```
With optional patient and RNA sample files:
```
$HISTOMX -m ${RMD_FILE} -f ${HISTOMX_PATH}/test_files/test.RCC -i 'test-run' -p ${HISTOMX_PATH}/test_files/patient-test.txt -r ${HISTOMX_PATH}/test_files/rna-test.txt
```
All output files are written to the same directory as the input RCC file (in this case 'test_files')

*Troubleshooting:*.  
If you see an error with 'pandoc' you may need to update your version of RStudio. See here:  
https://bookdown.org/yihui/rmarkdown-cookbook/install-pandoc.html

