library(xlsx)
library(lsa)
library(RTextTools)

#Genrating FoodEx2 hierarchy needed for training Pincare embeddings
foodEx2<-read.csv("../Foodex2-revision2.csv",sep=";")
 
 data<-foodEx2[foodEx2$statef %in% c("r","d","s","c"),]
 
 code<-as.character(data$code)
 parent<-as.character(data$reportParentCode)
 
 pairs<-cbind(parent,code)
 
 indices<-c()
 counter<-1
 
 for(i in 1:nrow(pairs)){
 	if(any(pairs[i,]=="")){
 		indices[counter]<-i
 		counter<-counter+1
 	}
 }
 
 pairs<-pairs[-indices,]
 write.csv(pairs,"../foodEx2_hierarchy.csv",sep=",")
 
#clustering with tSNE

library(xlsx)
library(VGAM)
library(caret)
library(dplyr)
library(cluster)
library(Rtsne) # for t-SNE plot
library(ggplot2) 

#Please change the path where your embeddings are stored 
data1<-read.table("../poincare_foodex2_100D.txt",sep=" ",header=FALSE)
data1<-data1[,2:ncol(data1)]
data1<-data.matrix(data1)
data1<-data.frame(data1)
 gower_dist <- daisy(data1,
                    metric = "gower",
                    type = list(logratio = 3))
                    gower_mat <- as.matrix(gower_dist)
                    
#Estiamte the number of clusters
sil_width <- c(NA)
for(i in seq(150,200,10)){
  pam_fit <- pam(gower_dist,
                 diss = TRUE,
                 k = i)   
 sil_width[i] <- pam_fit$silinfo$avg.width
 
}



plot(1:300, sil_width,
     xlab = "Number of clusters",
     ylab = "Silhouette Width")
lines(1:300, sil_width)

#Use the number of optimal clusters to perform the clustering
 pam_fit <- pam(gower_dist, diss = TRUE, k = 230)
 data1[pam_fit$medoids, ]
 tsne_obj <- Rtsne(gower_dist,check_duplicates = TRUE, is.distance=TRUE)
 tsne_data <- tsne_obj$Y %>%
  data.frame() %>%
     setNames(c("X", "Y")) %>%
  mutate(cluster = factor(pam_fit$clustering)
 
)
ggplot(aes(x = X, y = Y), data = tsne_data) +
geom_point(aes(color = cluster))
pam_results <- data1 %>%
  # dplyr::select(GO.) %>%
   mutate(cluster = pam_fit$clustering) %>%
   group_by(cluster) %>%
   do(the_summary = summary(.))
pam_results$the_summary



 
 
 #classificaiton
 
 
 labels<-as.character(data$statef)
 code<-as.character(data$code)
 #Please change the path where your embeddings are stored 
 data1<-read.table("../poincare_foodex2_100D.txt",sep=" ",header=FALSE)
 
list_term<-as.character(data1$V1)

data1<-data1[,2:ncol(data1)]
data1<-data.matrix(data1)

target<-c()
for(i in 1:length(list_term)){
	
	index<-which(code==list_term[i])
	if(length(index)!=0){
		target[i]<-labels[index]
		}
	if(length(index)==0){
		target[i]<-"XX"
	}
}
 

data_class<-cbind(data1,target)

indices<-which(target=="XX")

data_class<-data_class[-indices,]




container <- create_container(data1[-indices,], data_class[,ncol(data_class)], trainSize=1:5003, virgin=FALSE)

 SVM <- cross_validate(container, 10, "SVM")
 SLDA <- cross_validate(container, 10, "SLDA") 
 BAGGING <- cross_validate(container, 10, "BAGGING")
 BOOSTING <- cross_validate(container, 10, "BOOSTING")
 RF <- cross_validate(container, 10, "RF") 
 NNET <- cross_validate(container, 10, "NNET") 
 TREE <- cross_validate(container, 10, "TREE")

  
#food matching
library(coreNLP)
initCoreNLP()

library(NLP)
library(openNLP)

names<-as.character(data$name)
names<-gsub("[[:punct:]]", "", names)
names<-iconv(names, "ISO_8859-2", "UTF-8")
names<-tolower(names)



chunk_annotator<-Maxent_Chunk_Annotator()
 
sent_token_annotator <- Maxent_Sent_Token_Annotator() 
word_token_annotator <- Maxent_Word_Token_Annotator() 
pos_tag_annotator <- Maxent_POS_Tag_Annotator() 


POStagging<-function(s){
  annotObj <- annotateString(s)
   tokens<-getToken(annotObj)
   
   words<-tokens$token
   tags<-tokens$POS
   lemmas<-tokens$lemma
  
   df<-data.frame(words,tags,lemmas)
   return(df)
   
   
 } 


 Default_chunking<-function(s){
	a2 <- annotate(s, list(sent_token_annotator, word_token_annotator))
 	a3 <- annotate(s, pos_tag_annotator, a2)
 	a3w <- subset(a3, type == "word") 
	tags <- sapply(a3w$features, `[[`, "POS")
 	
 	a4<-annotate(s, chunk_annotator, a3)
 	a4w <- subset(a4, type == "word") 
 	chunks<-sapply(a4w$features, `[[`, "chunk_tag")
 	
 	beginnings1<- grep("^B",chunks)
 	beginnings2<-grep("^O",chunks)
 	beginnings<-sort(c(beginnings1,beginnings2))
 	
 	
 	return(list(chunks=chunks,beginnings=beginnings))
 	
 }

Nouns<-c()
Adjectives<-c()
Verbs<-c()
Length<-c()
NP_number<-c()
for(i in 1:length(names)){
  nouns<-c()
  adjectives<-c()
  verbs<-c()
  if(names[i]!=""){
    annoObj<-annotateString(names[i])
    annotations<-getToken(annoObj)
    chunks<-Default_chunking(names[i])$chunks
    words<-annotations$token
    tags <- annotations$POS
    lemmas<-annotations$lemma
    nouns<-lemmas[grep("NN[SP(PS)]?",tags)]
    adjectives<-lemmas[grep("JJ[RS]?",tags)]
    verbs<-lemmas[grep("VB[DGNPZ]?",tags)]
    Nouns[i]<-paste(nouns,collapse=".",sep=".")
    Adjectives[i]<-paste(adjectives,collapse=".",sep=".")
    Verbs[i]<-paste(verbs,collapse=".",sep=".")
    Length[i]<-nrow(annotations)
	NP_number[i]<-length(grep("B-NP",chunks))

  }
  if(names[i]==""){
    Nouns[i]<-paste(nouns,collapse=".",sep=".")
    Adjectives[i]<-paste(adjectives,collapse=".",sep=".")
    Verbs[i]<-paste(verbs,collapse=".",sep=".")
    Length[i]<-0
	NP_number[i]<-0
  }
}
Nouns_foodex2_lemma<-Nouns
Adjectives_foodex2_lemma<-Adjectives
Verbs_foodex2_lemma<-Verbs
Length_foodex2<-Length
NP_number_foodex2<-NP_number




indices<-c(990,1846,53,2557,1864,215,2807,231)
ENGFDNAM<-names[indices]

Nouns<-c()
Adjectives<-c()
Verbs<-c()
Length<-c()
NP_number<-c()
for(i in 1:length(ENGFDNAM)){
	nouns<-c()
	adjectives<-c()
	verbs<-c()
	if(ENGFDNAM[i]!="" && !is.na(ENGFDNAM[i])){
	annoObj<-annotateString(ENGFDNAM[i])
	annotations<-getToken(annoObj)
	chunks<-Default_chunking(ENGFDNAM[i])$chunks
	
	words<-annotations$token
	tags <- annotations$POS
	lemmas<-annotations$lemma
	nouns<-lemmas[grep("NN[SP(PS)]?",tags)]
	adjectives<-lemmas[grep("JJ[RS]?",tags)]
	verbs<-lemmas[grep("VB[DGNPZ]?",tags)]
	Nouns[i]<-paste(nouns,collapse=".",sep=".")
	Adjectives[i]<-paste(adjectives,collapse=".",sep=".")
	Verbs[i]<-paste(verbs,collapse=".",sep=".")
	Length[i]<-nrow(annotations)
	NP_number[i]<-length(grep("B-NP",chunks))
	}
	if(ENGFDNAM[i]=="" || is.na(ENGFDNAM[i])){
		Nouns[i]<-paste(nouns,collapse=".",sep=".")
	Adjectives[i]<-paste(adjectives,collapse=".",sep=".")
	Verbs[i]<-paste(verbs,collapse=".",sep=".")
	Length[i]<-0
	NP_number[i]<-0
	}
}

Nouns_fir_food_lemma<-Nouns
Adjectives_fir_food_lemma<-Adjectives
Verbs_fir_food_lemma<-Verbs
Length_fir_food<-Length
NP_number_fir_food<-NP_number


codes<-c()
p_max<-c()

prob<-list()

	names_matched<-c()
for(i in 1:length(Nouns_fir_food_lemma)){
 	if(Nouns_fir_food_lemma[i]!="" && !is.na(Nouns_fir_food_lemma[i])){
 	
 		nouns_temp<-strsplit(Nouns_fir_food_lemma[i],"\\.")[[1]]
 		adjectives_temp<-strsplit(Adjectives_fir_food_lemma[i],"\\.")[[1]]
 		verbs_temp<-strsplit(Verbs_fir_food_lemma[i],"\\.")[[1]]
 		
 		index_DS<-c()
 		for(j in 1:length(nouns_temp)){
 			index_DS<-c(index_DS,grep(nouns_temp[j],Nouns_foodex2_lemma,fixed=TRUE))
 		}
 		
 		index_DS<-setdiff(index_DS,c(i))
 		
 		index_DS<-unique(index_DS)
 		
 		sort(index_DS,decreasing=FALSE)
 		
 		p<-c()
 		
 		if(length(index_DS)!=0){
 		for(k in 1:length(index_DS)){
 			nouns_candidate<-strsplit(Nouns_foodex2_lemma[index_DS[k]],"\\.")[[1]]
 			adjectives_candidate<-strsplit(Adjectives_foodex2_lemma[index_DS[k]],"\\.")[[1]]
 			verbs_candidate<-strsplit(Verbs_foodex2_lemma[index_DS[k]],"\\.")[[1]]
 			
 			
 			p[k]<-(length(intersect(nouns_temp,nouns_candidate)))/(length(union(nouns_temp,nouns_candidate)))*(length(intersect(union(adjectives_temp,verbs_temp),union(adjectives_candidate,verbs_candidate)))+1)/(length(union(union(adjectives_temp,verbs_temp),union(adjectives_candidate,verbs_candidate)))+2)  
 					}
 		
 		
 		
 		prob[[i]]<-p
 			indices<-sort(p,decreasing=TRUE,index.return=TRUE)$ix
 			p_max[i]<-max(p)
 			if(length(indices)>0){
 		#	codes[i]<-paste(codes_foodex2[index_DS[indices]],collapse=".",sep=".")
 			names_matched[i]<-paste(names[index_DS[indices[1:6]]],collapse=".",sep=".")
 			 			}
 			
 		}
 		
 		
 		}
 		
 	
 	if(Nouns_fir_food_lemma[i]=="" || is.na(Nouns_fir_food_lemma[i])){
 		#codes[i]<-0
 	}
 	
 		
	
}

 A<-rbind(strsplit(names_matched[1],"\\.")[[1]],strsplit(names_matched[2],"\\.")[[1]],strsplit(names_matched[3],"\\.")[[1]],strsplit(names_matched[4],"\\.")[[1]],strsplit(names_matched[5],"\\.")[[1]],strsplit(names_matched[6],"\\.")[[1]],strsplit(names_matched[7],"\\.")[[1]],strsplit(names_matched[8],"\\.")[[1]])

rownames(A)<-ENGFDNAM


#embedding matching
test_code<-code[indices]

#Please change the path where your embeddings are stored
data1<-read.table("../poincare_foodex2_200D.txt",sep=" ",header=FALSE)

list_terms<-as.character(data1$V1)
data1<-data1[,2:ncol(data1)]
data1<-data.matrix(data1)


match<-list()
for (i in 1:length(test_code)){
	index1<-which(list_terms==test_code[i])
	sim<-c()
	for(j in 1:length(list_terms)){
		
		sim[j]<-cosine(data1[index1,],data1[j,])
	}
	
	ind<-sort(sim,decreasing=TRUE,index.return=TRUE)$ix
	match[[i]]<-ind[1:6]
	
}

B<-rbind(list_terms[match[[1]]],list_terms[match[[2]]],list_terms[match[[3]]],list_terms[match[[4]]],list_terms[match[[5]]],list_terms[match[[6]]],list_terms[match[[7]]],list_terms[match[[8]]])

for(i in 1:nrow(B)){
	for(j in 1:ncol(B)){
		index<-which(code==B[i,j])
		B[i,j]<-names[index]
	}
}


r<-sample(which(data$statef=="r"),100)
d<-sample(which(data$statef=="d"),100)
c<-sample(which(data$statef=="c"),100)
s<-sample(which(data$statef=="s"),100)

indices<-c(r,d,c,s)
indices<-sort(indices)

x_temp<-indices

indices<-x_temp


A<-matrix(0,2400,2)
counter1<-1
for(i in 1:400){
	x<-strsplit(names_matched[i],"\\.")[[1]]
	for(j in 1:length(x)){
		A[counter1,1]<-ENGFDNAM[i]
		A[counter1,2]<-x[j]
		counter1<-counter1+1
	}
}

#matching

data_temp<-cbind(list_term,target)

r<-sample(which(target=="r"),100)
d<-sample(which(target=="d"),100)
c<-sample(which(target=="c"),100)
s<-sample(which(target=="s"),100)
indices<-c(r,d,c,s)
indices<-sort(indices)

x_temp<-indices

indices<-x_temp

test_code<-list_term[indices]



A<-matrix(0,2400,2)
counter1<-1

for(i in 1:length(match)){
	index1<-which(code==test_code[i])
	first_temp<-names[index1]
	for(j in 1:length(match[[i]])){
		A[counter1,1]<-first_temp
		index2<-which(code==list_terms[match[[i]][j]])
		A[counter1,2]<-names[index2]
		counter1<-counter1+1
	}
	
	
	
}

