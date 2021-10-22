**Generate a HistoMx report for a single RCC file**

R must be installed on your machine. You can dowload it here:  
https://cran.r-project.org/bin/macosx/

1. Open the 'Terminal' application on your computer.

2. Change to wherever you want to use as your working directory:
```
cd ~/Desktop
```

3. Install all R packages and dependencies.  

First, download the following file:  
https://drive.google.com/file/d/1V89z4OYOJxdMUtL54jKeg5CrP-SIR32q/view?usp=sharing

Then run the following commands from the commandline in or in RStudio. This will only install missing packages and dependencies.  
NOTE: some packages may need to be installed from Bioconductor: eg.```BiocManager::install("DESeq2")```
```
R

source("~/Downloads/histomx_pckgs.Rdmpd")
ip <- as.data.frame(installed.packages())
to_install <- setdiff(histomx_pckgs$package, ip$Package)
install.packages(to_install)

quit(save="no")
```

4. Dowload the histomx code:  
* Click the 'Code' button at the top of this repository and then 'Download ZIP'
* Unzip and copy this directory ('histomx-main') to your working directory (eg. ~/Desktop) and rename it to 'histomx'

5. In the terminal window, navigate to the histomx directory:
```
cd ~/Desktop/histomx
```

6. Set variables to define paths to scripts/files
```
HISTOMX_PATH=~/Desktop/histomx 

HISTOMX=${HISTOMX_PATH}/bin/render-histomx_report  

RMD_FILE=${HISTOMX_PATH}/scripts/histomx_report.Rmd
```

7. Create a directory called 'refRCCs'.
```
mkdir -p $HISTOMX_PATH/refRCCs
```

8. Download all reference RCC files:  
https://drive.google.com/drive/folders/1Wzi9LCof7QMcYyx7kLOiKWuBFEY2o8Zk?usp=sharing

9. Move all RCC files to this directory:
```
mv ~/Downloads/RefSet-KTD1-FFPE/* refRCCs/
```

10. Run histomx on a single RCC file to generate the report
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

*Troubleshooting:*  
If you see an error with 'pandoc' you may need to update your version of RStudio. See here:  
https://bookdown.org/yihui/rmarkdown-cookbook/install-pandoc.html

> **To us git version control:**
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

