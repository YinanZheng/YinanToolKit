### Use Bioc-devel
### This is R Code
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "devel")
BiocManager::valid()              # checks for out of date packages


## Ctrl+Alt+Enter submit code
## passphrase for REMP : zyn1029

## Origin: Github
## upstream: Bioconductor
## Master: devel repository

## For a fresh start:
## First clone repository from Github (using Github Desktop)

## In the local repository directory, Sync an existing repositories
cd D:/Github/REMP
cd /Users/yzk256/Documents/GitHub

git remote add upstream git@git.bioconductor.org:packages/REMP.git

git fetch --all

git checkout master
git merge origin/master
git merge upstream/master
git push upstream master
git push origin master

git checkout RELEASE_3_10
git merge upstream/RELEASE_3_10
git merge origin/RELEASE_3_10
git push upstream RELEASE_3_10
git push origin RELEASE_3_10

###############################################
###############################################

## To fix bugs in devel and release:

git fetch --all

git checkout master
git merge origin/master
git merge upstream/master

git checkout RELEASE_3_10
git merge upstream/RELEASE_3_10
git merge origin/RELEASE_3_10


git checkout master

### Start Making changes to devel version
### ....
### ....
### Remember to update README.md and NEWS

### Add modified files (-A: all files)
git add -A

### after version bump  
git add DESCRIPTION


git commit -m "REMP 1.11.1 updates"



### If necessary, update release verion simutanuously
git checkout RELEASE_3_10
git cherry-pick master

### For release version (y minus 1)
### Fix version bump by editing 'Version:' field of DESCRIPTION, then
git add DESCRIPTION
git commit -m "REMP 1.10.1 updates"
### 

### Push and update devel
git checkout master
git push upstream master
git push origin master

### Push and update release (only current release can be updated)
git checkout RELEASE_3_10
git push upstream RELEASE_3_10
git push origin RELEASE_3_10








