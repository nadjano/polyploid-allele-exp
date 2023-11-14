## nextflow-draft:

In this repository 
I have created template for Nextflow -code . which can be use as base for Omics Analysis using Nextflow



**Step1: Set up AWS cloud9 :**

First we have to setup aws account as explained in setup.md file present in setup folder



**Step2: Create nextflow scripts :**

According to our need we have to write nectflow scripts . 
I have splitted nextflow scripts according to my wish as below :

**step3: nextflow.config :**

In this i have added My infomration manifest and how nextflow script should run either locally or cloud or through clusters.

- Added Manifest block
- Added Profiles - like local or AWS 
- Added process block 
- 

**step4 : main nextflow script :**

- I have added main.nf script which includes or call different modules from modules folder
- Also this main.nf has input output directories path

**step4 : modules nextflow scripts :**

- For each tasks i have created seperate nextflow scripts so i can use them in main nextflow script
- This is helpful to edit or add more scripts according to our needs 
