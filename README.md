**Generate a HistoMx report for a single RCC file**

R must be installed on your machine. You can download it here:  
https://cran.r-project.org/bin/macosx/

NOTE: These instructions are for Mac only (it should theoretically work in windows bash. Another option is to run the R markdown command directly in R which can be found in the bin file) 

1. Open the 'Terminal' application on your computer.

2. Change to wherever you want to use as your working directory:
```
cd ~/Desktop
```

3. Dowload the histomx code:  
* Click the 'Code' button at the top of this repository and then 'Download ZIP'
* Unzip and copy this directory ('histomx-main') to your working directory (eg. ~/Desktop)

4. Install all R packages and dependencies.  

Run the following commands in the terminal or in RStudio. This will only install missing packages and dependencies.  

```
R

source("~/Desktop/histomx-main/static/histomx_pckgs.Rdmpd")
ip <- as.data.frame(installed.packages())
to_install <- setdiff(histomx_pckgs$package, ip$Package)
BiocManager::install(to_install) 

quit(save="no")
```

5. Set variables to define paths to scripts/files
```
HISTOMX_PATH=~/Desktop/histomx-main  

HISTOMX=${HISTOMX_PATH}/bin/render-histomx_report  
```

The type of report is determined by the Markdown file used (choose one):
```
RMD_FILE=histomx_kidney.Rmd (standard BHOT report)
RMD_FILE=histomx_kidney_bkv (standard BHOT report with BKV expression)
RMD_FILE=histomx_kidney_simple (simplified report w/o pathways)
RMD_FILE=histomx_xeno.Rmd (BPOT report)
```

6. Create a directory called 'refRCCs'.
```
mkdir -p $HISTOMX_PATH/refRCCs
```

7. Download all reference RCC files from the server.  

8. Move all RCC files to this directory:
```
mv ~/Downloads/RefSet-KTD1-FFPE/* refRCCs/
```

9. Run histomx on a single RCC file to generate the report
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


**TROUBLESHOOTING:**  
If you see an error with 'pandoc' you may need to update your version of RStudio. See here:  
https://bookdown.org/yihui/rmarkdown-cookbook/install-pandoc.html

> **To use git version control:**
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

