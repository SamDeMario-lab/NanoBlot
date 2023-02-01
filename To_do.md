# To Do

## Code

### Check YML file
MAC
LINUX
WINDOWS

### Add step checking probe specificity
Just before writing the final bam we should check for a multimapped flag.
  Actually I was thinking about this more and we should put it in the R Package and do the check after we read in the samples.

### DESEQ 
Review 3 suggested we use it for replicates. The only reasonable way to do this is to just check if the overall level of a transcript changes. 
so the psuedocode for this is to make the normalizeNanoblot take an argument that is basically the DeSeq info argument --> and basically helps them call DeSeq but they are literally calling DeSeq --> and then have the DeSeq just output the Deseq adjusted or whatever table it normally outputs 
And then make the plotting function take different plotting infos 

### Fix Rsamtools BamFileList()
We have a function called Scanbamfiles() which we should replace with BamFileList()
I think we only have to change 2 functions

### Check integrity of user specified regions.  

### Clean GitHub Repo

### Test the checkForMultimapping() function
I want to make sure it behaves as expected

### README
## Remove dependencies
## Add Function Descriptions

## Genes for multimapping test
TDH1
UBI UB4
