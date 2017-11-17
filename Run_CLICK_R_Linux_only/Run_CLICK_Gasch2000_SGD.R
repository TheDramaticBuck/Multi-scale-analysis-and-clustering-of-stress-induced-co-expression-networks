# Code for running CLICK from R as an executable
# See also http://acgt.cs.tau.ac.il/expander/

# Thomas Leja and Nuno Rocha Nene


#Load data 

setwd('.') # Change working directory to folder where current file is held

stress='./Gasch_SGD/diamideTreatment/' ### Other stresses: HOtimeCourse; HS25-37; HS29-33; HS30-37; HS37-25; HSmild1; HSmild2; Hyper-osmotic; Hypo-osmotic; menadione.

stress=paste(stress,"Data_Gasch2000_SGD.R",sep="")

source(stress) ### expression data (SGD)

genenames='./Gasch_SGD/diamideTreatment'

genenames=paste(genenames,"/Data_Gasch2000_SGD_genenames.R",sep="")

source(genenames)

data=datamatrix_select

rownames(data) <- genenames_select

data.table <- as.table(data)

# Write data to clickInput.orig file

# CLICK requires the first line to contain the dimensions (space separated) of the data matrix (rows x columns). 
# The following lines should contain data.matrix (tab separated). 
# If row names contain special characters like ".", "/" or space the click will freze, I usally replace any special characters with underscore before running click.

write(dim(data), file="clickInput.orig", sep=" ")

write.table(data, file="clickInput.orig", sep="\t", col.names=F,row.names=T,quote=F,append=T)


# CLICK LINUX server


path = "./clickLinux" # do not use "/" at the end of this line, it will be added in the loop function below.
setwd(path) # change path to where exe is located

dir.create("./clickLinux/output")

exp.hom = format(seq(0.5,1,by=0.01),digits=3) # specify the homogeneity range since CLICK requires this input

for (i in 1:length(exp.hom)) {
	p = read.delim("clickParams.txt",header=F,stringsAsFactors=F)[,1]
	p[6] = paste(path,"/clickLinux/output/click_",exp.hom[i]," ",sep="")
	p[10] = paste(exp.hom[i]," ",sep="")
	write.table(p,file="./clickParams.txt",col.names=F,row.names=F,quote=F)
	
	system("./clickLinux/click.exe clickParams.txt")
}


### Import solutions to a matrix format

setwd("./clickLinux/output")

res = list.files(pattern="[.]res[.]sol")

m = read.delim(res[1],header=F,row.names=1,stringsAsFactors=F)
for (i in 2:length(res)) {
	tmp = read.delim(res[i],header=F,row.names=1,stringsAsFactors=F)
	m = cbind(m, tmp)
}
colnames(m) = gsub("[.]res[.]sol","",res)

setwd("../")

write.csv(m,paste("./Gasch_SGD/",Gasch_SGD_exp,file="clickSolutions.csv"," ",sep=""), quote=F)

### CLICK solutions can be evaluated under the Hav, Sav, GOTOav performance functions too.


