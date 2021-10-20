**Generate HistoMx report for a single RCC file**

1. Open the 'terminal' application on your computer.

2. Change to wherever you want to use as your working directory:
```
cd ~/Desktop
```

**OPTION 1 (without git):**
* Click the 'Code' button at the top of this repository and then 'Download ZIP'
* Unzip and copy this directory (histomx-main) to your working directory (eg. ~/Desktop)

**OPTION 2 (with git version control):**
* Install git by entering the following in the terminal window:
```
brew install git
```
Other installation options:https://www.atlassian.com/fr/git/tutorials/install-git

* Type the following command in order to clone this repository (ie. copy everything to your working directory):
```
git clone https://github.com/dinovski/histomx.git
```

3. Set variables to define paths to scripts/files
```
HISTOMX_PATH='~/Desktop/histomx/'
HISTOMX='~/Desktop/histomx/bin/render-histomx_report'
RMD_FILE='~/Desktop/scripts/histomx_report.Rmd'
```

4. Run histomx to generate the report
* This will output an .html file (you can specify the name of the file with the '-i' argument).
* The -m (or --rmd) and -f (or --rcc) arguments are required

Running the HISTOMX command will show all possible input arguments:
```
$HISTOMX
```

The following are 2 examples (with and without patient or RNA QC files)
```
$HISTOMX -m ${RMD_FILE} -f ${HISTOMX_PATH}/sample_files/test.RCC -i 'test-run'
```
With optional patient and RNA sample files:
```
$HISTOMX -m ${RMD_FILE} -f ${HISTOMX_PATH}/sample_files/test.RCC -i 'test-run' -p ${HISTOMX_PATH}/sample_files/test-patient.txt -r ${HISTOMX_PATH}/sample_files/test-rna.txt
```
