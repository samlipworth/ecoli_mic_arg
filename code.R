
##############libraries##############

library(data.table)
library(tidyverse)
library(survival)
library(patchwork)

ec_match<-read_csv('./main_data.csv') 
tem<-read_csv('./tem.csv')


###########gent###########################
#read in data

ec_match2<-select(ec_match,guuid,Gentamicin_lower,Gentamicin_upper) %>% filter(!is.na(Gentamicin_upper))
amrfinder<-read_tsv('amrfinder.tsv')
ec_match2<-filter(ec_match2,guuid %in% amrfinder$Name)
#select relevant genes from amrfinder
amrfinder<-amrfinder%>% filter(grepl('AMINOGLYCOSIDE',Class )) %>%  select(guuid=Name,gene=`Gene symbol`)
length(unique(amrfinder$gene))


#if gene occurs >=10 times it is called as it is in amrfinder, 5-10 as "gene_other" and <5 as "other gene"
gent<-read_tsv('amrfinder.tsv')  %>% filter(Name %in% ec_match2$guuid) %>% filter(grepl('GENTAMICIN',Subclass )) %>%  select(guuid=Name,gene=`Gene symbol`)
length(unique(gent$gene))
count<-amrfinder %>% filter(str_length(guuid) ==36) %>%  group_by(gene) %>% count() %>% filter(n>=10)

common_genes<-filter(amrfinder,gene %in% count$gene)
length(unique(amrfinder$gene)) 
length(unique(common_genes$gene)) 

deal_with_rare_genes <- function(gene_name) {
  if (startsWith(gene_name, "ampC")) {
    gsub("^ampC_.*", "ampC_other", gene_name)
  } else if (!grepl("_", gene_name)) {
    # Append '_other' at the end
    paste0(gene_name, "_other")
  } else {
    # Replace after the last '_'
    gsub("_[^_]+$", "_other", gene_name)
  }
}

amrfinder$gene <- ifelse(amrfinder$gene %in% count$gene,
                         amrfinder$gene,
                         sapply(amrfinder$gene, deal_with_rare_genes))

singletons<-amrfinder %>%  group_by(gene) %>% count() %>% filter(n<5)

amrfinder$gene<-ifelse(amrfinder$gene %in% singletons$gene,"other gene",amrfinder$gene)
amrfinder_save1<-amrfinder

#amrfinder$gene<-ifelse(amrfinder$gene %in% count$gene,amrfinder$gene,'other gene')
x<-data.frame(ec_match2$guuid)
x<-left_join(x,amrfinder,by=c("ec_match2.guuid"="guuid"))
x<-select(x,guuid=ec_match2.guuid,gene)
amrfinder<-data.frame(table(x)) %>% spread(gene,Freq)
names(amrfinder)<-str_replace_all(names(amrfinder),'-','_')
amrfinder[,2:ncol(amrfinder)][amrfinder[,2:ncol(amrfinder)]>1]<-1


genes<-names(amrfinder)
genes<-genes[!genes =='guuid']

amrfinder<-left_join(amrfinder,ec_match2,by=c("guuid"="guuid"))

mlst<-read_tsv('./mlst.tsv',col_names = F) %>% select(guuid=X1,spec=X2,mlst=X3)



big<-c('131','95','69','73')
mlst$mlst<-ifelse(mlst$mlst %in% big, mlst$mlst, 'other')

phylogroups<-read_tsv('./phylogroups.tsv')
mlst<-left_join(mlst,phylogroups,by=c("guuid"))
mlst$mlst<-ifelse(mlst$mlst =='other',mlst$phylogroup,mlst$mlst)
other<-c('U','U/cryptic','cryptic','E','G')
mlst$mlst<-ifelse(mlst$mlst %in% other,'other',mlst$mlst)


mlst<-select(mlst,guuid,mlst)

amrfinder<-left_join(amrfinder,mlst,by=c("guuid"))

out=NULL

names(amrfinder)<-janitor::make_clean_names(names(amrfinder))
genes<-janitor::make_clean_names(genes)


genes<-append(genes,"mlst")
amrfinder$mlst<-factor(amrfinder$mlst,levels=c('B2','131','95','73','69','A','B1','C','D','F','other'))
(Y <- with(amrfinder, Surv(gentamicin_lower, gentamicin_upper, type = "interval2")))
x<-amrfinder$guuid
amrfinder<-select(amrfinder,-guuid,-gentamicin_lower,-gentamicin_upper)

a<-amrfinder[,grep("^V",colnames(amrfinder))]

#univariable models

for(gene in genes){
  print(which(genes ==gene))
  
  print(gene)
  
  b<-select(amrfinder,gene)
  
  #c<-cbind(b,a)
  c<-b
  t<-as.formula(paste("Y", paste(names(c), collapse = " + "), sep=' ~ '))
  
  
  
  m<-broom::tidy(survreg(t ,data=c,dist="gaussian"),conf.int=T)
  
  m<-if(gene=='mlst'){m[2:11,]}else{m[2,]}    
  out=rbind(out,data.frame(m))
}
out_save<-out
include<-out

amrfinder_save<-amrfinder
big<-c('st131','st95','st73','st69','other')
include<-filter(include,!gene %in% big)

#multivariable models

amrfinder<-amrfinder[names(amrfinder) %in% include$term | names(amrfinder) %in% names(a) ]
t<-as.formula(paste("Y", paste(names(amrfinder), collapse = " + "), sep=' ~ '))
model<-survreg(t, data=amrfinder,dist="gaussian")
out2<-broom::tidy(survreg(t, data=amrfinder,dist="gaussian"),conf.int=T)
out2$which<-"multivariable no population"

amrfinder<-amrfinder_save
include<-out
amrfinder<-amrfinder[names(amrfinder) %in% include$term | names(amrfinder) %in% names(a) | names(amrfinder) =='mlst' ]
t<-as.formula(paste("Y", paste(names(amrfinder), collapse = " + "), sep=' ~ '))
model<-survreg(t, data=amrfinder,dist="gaussian")
out3<-broom::tidy(survreg(t, data=amrfinder,dist="gaussian"),conf.int=T)
out3$which<-"multivariable population"

#here an "r" gene is associated with resistance whereas a "line" gene is independently resistance conferring
#we worry about significance later
gent_r_genes<-filter(out3,2^(estimate) >1) %>% filter(!grepl('mlst',term)) 
gent_line<- 2^(log2(2) - model$coefficients[1])
gent_line_genes<-filter(out3,2^estimate >gent_line) %>% filter(!grepl('mlst',term))
gent_all_genes<-filter(out3,!grepl('mlst',term))

out$which<-"univariable"
out<-rbind(out,out2,out3)

count$gene<-janitor::make_clean_names(count$gene)

out<-filter(out,term %in% count$gene| grepl('mlst',term) | term=='other_gene')



ceph_genes<-read_tsv('amrfinder.tsv') %>% filter(grepl('GENTAMICIN',Subclass))
ceph_genes<-unique(ceph_genes$`Gene symbol`)
ceph_genes<-janitor::make_clean_names(ceph_genes)

a<-ifelse(out$term %in% ceph_genes,'red','black')

out<-filter(out,!grepl('blaEC',gene))
out$estimate<-as.numeric(out$estimate)
out$conf.high<-as.numeric(out$conf.high)
out$conf.low<-as.numeric(out$conf.low)

X<-c("mlstB1"  ,     "mlstA"    ,   "mlstC"  ,      "mlstD" ,       "mlstE",       
     "mlstF",        "mlstG" ,       "mlstother") 
out$term<-ifelse(out$term %in% X, str_replace(out$term,'mlst','phylogroup '),ifelse(grepl('mlst',out$term),str_replace(out$term,'mlst','ST '),out$term))

out<- out %>% 
  add_row(term = "GAP") %>% 
  add_row(term = "GAP2")

out$estimate<-2^out$estimate
out$conf.high<-2^out$conf.high
out$conf.low<-2^out$conf.low


out$which<-ifelse(grepl("GAP",out$term),'univariable',out$which)
out$which<-factor(out$which,levels=c("univariable","multivariable no population","multivariable population"))

amrfinder<-amrfinder_save1
m<-select(mlst,guuid,gene=mlst) %>% filter(!is.na(gene))
amrfinder<-rbind(amrfinder,m)
count<-amrfinder %>% group_by(gene) %>% distinct(guuid) %>%  count()
count$gene<-ifelse(str_length(count$gene) >3,janitor::make_clean_names(count$gene),count$gene)
count$gene <- case_when(count$gene == '131' ~ 'ST 131',
                        count$gene == '69' ~ 'ST 69',
                        count$gene == '95' ~ 'ST 95',
                        count$gene == '73' ~ 'ST 73',
                        count$gene == 'A' ~ 'phylogroup A',
                        count$gene == 'B1' ~ 'phylogroup B1',
                        count$gene == 'B2' ~ 'phylogroup B2',
                        count$gene == 'C' ~ 'phylogroup C',
                        count$gene == 'D' ~ 'phylogroup D',
                        count$gene == 'F' ~ 'phylogroup F',
                        count$gene == 'other' ~ 'phylogroup other',
                        TRUE ~ count$gene)


out<-left_join(out,count,by=c("term"="gene"))
out$term<-paste0(out$term, ' n=(',out$n,')')


out$term<-str_replace_all(out$term," n=[(]NA[)]",'')
out$term<-factor(out$term, levels=c("aac_3_i_id n=(93)","aac_3_i_ie n=(77)","aac_6_ib_cr5 n=(126)","aad_a1 n=(342)",
                                    "aad_a2 n=(48)","aad_a5 n=(311)","ant_2_ia n=(10)","aph_3_ib n=(749)",      
                                    "aph_3_ia n=(89)","aph_6_id n=(746)","other_gene n=(24)","GAP","GAP2","ST 131 n=(646)",
                                    "ST 95 n=(463)","ST 73 n=(683)","ST 69 n=(413)","phylogroup A n=(246)",
                                    "phylogroup B1 n=(267)","phylogroup C n=(101)","phylogroup D n=(285)","phylogroup F n=(136)",    
                                    "phylogroup other n=(118)"))

gentr<-c("aac_3_i_id n=(93)","aac_3_i_ie n=(77)","ant_2_ia n=(10)")

a<-ifelse(out$term %in% gentr,"red","black")
a[10:55]<-"black"
gent_plot<-ggplot(out) +
  aes(x=term,y=estimate,color=which,group=which) +
  geom_point(position = position_dodge(0.9)) + geom_errorbar(data=out,aes(ymin=conf.low,ymax=conf.high),position=position_dodge(width=0.9),width=0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90,color=a),panel.grid.major.x = element_blank()) +
  xlab('Gene') + ylab('Fold-change in MIC') + 
  labs(color="Model") +
  geom_hline(yintercept =2^(log2(2) - model$coefficients[1]),linetype="dashed",color="dark red") +
  geom_hline(yintercept=1)+
  scale_x_discrete(labels = function(x) ifelse(x == "GAP" | x=="GAP2", "", x)) +
  scale_color_viridis_d() +
  scale_y_log10()


gent<-out
gent$drug<-"gentamicin"

marginaleffects::hypotheses(model, "aac_3_i_id = aac_3_i_ie")

amrfinder<-amrfinder_save
predictors <- out_save$term
predictors<-ifelse(grepl("mlst",predictors),"mlst",predictors)
predictors<-unique(predictors)
predictors<-predictors[!is.na(predictors)  & ! predictors =="(Intercept)" & ! predictors == "Log(scale)"]

out=NULL

#here we test all inveractions in a pairwise manner

for (excluded_gene in predictors) {
  print(excluded_gene)
  
  other_predictors <- setdiff(predictors, excluded_gene)
  
  for (current_gene in other_predictors) {
    
    interaction_term <- paste(excluded_gene, current_gene, sep="*")
   
    formula_str <- as.formula(paste( "Y ~", paste(setdiff(other_predictors,current_gene), collapse = "+"), "+", interaction_term))
    summary_model <- broom::tidy(survreg(formula_str, data = amrfinder, dist = "gaussian"))

    
    out <- rbind(out,data.frame(summary_model))
    
  }
}




replace_mlst_patterns <- function(texts) {
  # Replace mlst* before the colon
  texts <- gsub("mlst[A-Za-z0-9_]*:", "mlst:", texts)
  # Replace mlst* after the colon
  texts <- gsub(":[A-Za-z0-9_]*mlst[A-Za-z0-9_]*", ":mlst", texts)
  return(texts)
}


#only keep interaction terms where p<0.01
interaction_terms<-filter(out,p.value < 0.01) %>% filter(grepl(':',term))
interaction_terms$term<-replace_mlst_patterns(interaction_terms$term)
interaction_terms<-unique(interaction_terms$term)
interaction_terms

model<-survreg(Y ~ aac_3_i_id + aac_3_i_ie + aac_3_i_va_other + aac_6_ib_cr5 + 
                 aad_a1 + aad_a2 + aad_a5 + ant_2_ia + aph_3_ib + aph_3_ia + 
                 aph_6_id + other_gene   
                 ,data=amrfinder,dist="gaussian")

elimination<-  step(model, direction="backward")

selected_vars <- names(coef(elimination))
selected_vars<-selected_vars[!selected_vars=="(Intercept)"]
selected_vars<-c(selected_vars,"mlst")

interaction_string <- paste(interaction_terms, collapse = " + ")
new_formula <- as.formula(paste("Y ~ ", paste(selected_vars, collapse = " + "), " + ", interaction_string))
model<-survreg(new_formula,data=amrfinder,dist="gaussian")
elimination<-step(model,direction="backward",k=3.8)

###leave one out cross validation


ceph_gene<-read_tsv('amrfinder.tsv') %>% filter(Subclass=='GENTAMICIN') %>% filter(Scope=='core')


out<-NULL
amrfinder$guuid<-x
for (i in 1:nrow(amrfinder)) {
  print(i)
  model<-survreg( Y[-i,] ~ aac_3_i_id + aac_3_i_ie + aac_3_i_va_other + 
                    aac_6_ib_cr5 + aad_a1 + aad_a5 + ant_2_ia + other_gene + 
                    mlst + aad_a1:other_gene + aad_a5:ant_2_ia, data = amrfinder[-i,], 
                  dist = "gaussian")
  
  prediction <- predict(model, newdata = amrfinder[i, ])
  correct<-ifelse(prediction >= as.numeric(Y[i,1]) & prediction <= as.numeric(Y[i,2]),1,0)
  
  out<-rbind(out,data.frame(correct,prediction,as.numeric(Y[i,1]), as.numeric(Y[i,2]),amrfinder$guuid[i]))
}


out$amrfinder_prediction<-ifelse(out$amrfinder.guuid.i. %in% ceph_gene$Name,'R','S')
out$amrfinder_prediction<-as.factor(out$amrfinder_prediction)
#breakpoint >2 = log2(2) = 1
out$binary<-ifelse(out$as.numeric.Y.i..2.. >1,'R','S')
out$binary_prediction<-ifelse(out$prediction>1,'R','S')
out$binary<-as.factor(out$binary)
out$binary_prediction<-as.factor(out$binary_prediction)

extend_margin <- function(out, n){
  out$margin <- NA # Initialize margin column
  for(i in 0:n){
    condition <- out$prediction >= (out$as.numeric.Y.i..1.. - i) & 
      out$prediction <= (out$as.numeric.Y.i..2.. + i)
    out$margin[is.na(out$margin) & condition] <- i
  }
  return(out)
}


out2<-extend_margin(out,10)
out2$margin<-ifelse(out2$prediction < out2$as.numeric.Y.i..1..,out2$margin *-1,out2$margin)

out2$which<- ifelse(out2$binary_prediction==out2$binary & out2$binary=='S', 'Correct S',
                  ifelse(out2$binary_prediction==out2$binary & out2$binary =='R','Correct R',
                  ifelse(out2$binary_prediction=='S' & out2$binary =='R','VM',
                         ifelse(out2$binary_prediction =='R' & out2$binary =='S','M',NA))))

out$within_1<-ifelse(out$prediction >= (out$as.numeric.Y.i..1.. -1) & out$prediction <= (out$as.numeric.Y.i..2.. +1),1,0)
out$within_2<-ifelse(out$prediction >= (out$as.numeric.Y.i..1.. -2) & out$prediction <= (out$as.numeric.Y.i..2.. +2),1,0)


caret::confusionMatrix(out$binary_prediction,out$binary,positive='R')
table(out$correct)
prop.test(2412,(2412+457))
2412+457


table(out$within_1)
prop.test(2785,(2785+84))
2785+84

cm<-caret::confusionMatrix(as.factor(out$binary_prediction),reference=as.factor(out$binary),positive="R")
tn<-cm$table[2,2]
fn<-cm$table[2,1]
tp<-cm$table[1,1]
fp<-cm$table[1,2]
#concordance
prop.test(x=(tp+tn),n=(tp+tn+fp+fn))
tp+tn
tp+tn+fp+fn
#sens
prop.test(x=(tp),n=(tp+fn))
tp
tp+fn
#spec
prop.test(x=(tn),n=(tn+fp))
tn
tn+fp

#Maj
prop.test(x=(fp),n=(tn+fp))
fp
tn+fp
#Vmj
prop.test(x=(fn),n=(tp+fn))
fn
tp + fn

#PPV
prop.test(x=(tp),n=(tp+fp))
tp
tp+fp

#NPV
prop.test(x=(tn),n=(tn+fn))
tn
tn+fn

colors<-c("#09820B" ,"#03FF25" ,"#FFA703", "#ED3A37")
names(colors)<-c("Correct R","Correct S","M","VM")
gent_residual<-ggplot(out2) +
  aes(x=margin, fill=which) +
  geom_bar(stat="count") + theme_minimal() +
  #scale_fill_manual(values=colors) + 
  scale_x_continuous(breaks = seq(min(out2$margin),max(out2$margin),1)) + 
  xlab("True - Predicted MIC") +
  labs(fill="") +
  ggtitle("Gentamicin") +
  scale_fill_viridis_d()


amrfinder<-read_tsv('amrfinder.tsv')
gent<-select(ec_match,guuid,Gentamicin) %>% filter(guuid %in% amrfinder$Name)

amrfinder<-amrfinder%>% filter(grepl('AMINOGLYCOSIDE',Class )) %>%  select(guuid=Name,gene=`Gene symbol`)
length(unique(amrfinder$gene))

length(unique(gent$gene))
count<-amrfinder %>% filter(str_length(guuid) ==36) %>%  group_by(gene) %>% count() %>% filter(n>=10)


deal_with_rare_genes <- function(gene_name) {
  if (startsWith(gene_name, "ampC")) {
    gsub("^ampC_.*", "ampC_other", gene_name)
  } else if (!grepl("_", gene_name)) {
    # Append '_other' at the end
    paste0(gene_name, "_other")
  } else {
    # Replace after the last '_'
    gsub("_[^_]+$", "_other", gene_name)
  }
}

amrfinder$gene <- ifelse(amrfinder$gene %in% count$gene,
                         amrfinder$gene,
                         sapply(amrfinder$gene, deal_with_rare_genes))

singletons<-amrfinder %>%  group_by(gene) %>% count() %>% filter(n<5)

amrfinder$gene<-ifelse(amrfinder$gene %in% singletons$gene,"other gene",amrfinder$gene)


gent_r_genes<-filter(gent_r_genes,p.value<0.05)
amrfinder$gene<-janitor::make_clean_names(amrfinder$gene,allow_dupes = T)
amrfinder_r<-filter(amrfinder,gene %in% gent_r_genes$term)
amrfinder_r<- amrfinder_r %>% group_by(guuid) %>% count()

gent$Gentamicin<-factor(gent$Gentamicin,levels=c('<=1','2','4','>4'))
gent<-filter(gent,!is.na(Gentamicin))
amrfinder_r<-left_join(gent,amrfinder_r,by=c("guuid"="guuid"))
amrfinder_r$n<-ifelse(is.na(amrfinder_r$n),0,amrfinder_r$n)
amrfinder_r<-filter(amrfinder_r,!is.na(Gentamicin))

r<-c('4','>4')

amrfinder_r$R<-ifelse(amrfinder_r$Gentamicin %in% r,'R','S')
amrfinder_r$R<-factor(amrfinder_r$R,levels=c('S','R'))
gent_gene_dist<-ggplot(amrfinder_r) +
  aes(x=n,fill=R) +
  geom_histogram(bins=length(unique(amrfinder_r$n))) +
  facet_wrap(~Gentamicin,nrow = 1) + theme_minimal() +
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),strip.text = element_blank()) +
  xlab("Number of ARGs/mutations")+
  scale_fill_viridis_d() +
  labs(fill="")

ng<-amrfinder_r %>% 
  group_by(Gentamicin) %>% 
  summarise(prop_no_gene = sum(n==0)/n(),
            lower=prop.test(sum(n==0),sum(n()))$conf.int[1],
            upper=prop.test(sum(n==0),sum(n()))$conf.int[2])


gent_ng_plot<-ggplot(ng) +
  aes(x=Gentamicin,y=prop_no_gene) +
  geom_point() +
  geom_errorbar(data=ng,aes(ymin=lower,ymax=upper),width=0) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank()) +
  ylab("Proportion with no ARG/mutation") +
  xlab("MIC")


gent_eucast_plot<-ng_plot / gent_gene_dist  +
  plot_annotation(title = "Gentamicin")





###################cotrim ecoli#####################


ec_match2<-select(ec_match,guuid,Cotrim_lower,Cotrim_upper) %>% filter(!is.na(Cotrim_upper))
amrfinder<-read_tsv('amrfinder.tsv')
ec_match2<-filter(ec_match2,guuid %in% amrfinder$Name)
amrfinder<-amrfinder%>% filter(grepl('TRIMETHOPRIM',Subclass )|grepl('SULFONAMIDE',Subclass)) %>%  select(guuid=Name,gene=`Gene symbol`)
length(unique(amrfinder$gene))

gent<-read_tsv('amrfinder.tsv')  %>% filter(Name %in% ec_match2$guuid) %>% filter(grepl('TRIMETHOPRIM',Subclass )| grepl('SULFONAMIDE',Subclass)) %>%  select(guuid=Name,gene=`Gene symbol`)
length(unique(gent$gene))
count<-amrfinder %>% filter(str_length(guuid) ==36) %>%  group_by(gene) %>% count() %>% filter(n>=10)

deal_with_rare_genes <- function(gene_name) {
  if (startsWith(gene_name, "ampC")) {
    gsub("^ampC_.*", "ampC_other", gene_name)
  } else if (!grepl("_", gene_name)) {
    # Append '_other' at the end
    paste0(gene_name, "_other")
  } else {
    # Replace after the last '_'
    gsub("_[^_]+$", "_other", gene_name)
  }
}

amrfinder$gene <- ifelse(amrfinder$gene %in% count$gene,
                         amrfinder$gene,
                         sapply(amrfinder$gene, deal_with_rare_genes))

singletons<-amrfinder %>%  group_by(gene) %>% count() %>% filter(n<5)

amrfinder$gene<-ifelse(amrfinder$gene %in% singletons$gene,"other gene",amrfinder$gene)
amrfinder_save1<-amrfinder

common_genes<-filter(amrfinder,gene %in% count$gene)
length(unique(amrfinder$gene)) # 28 aminoglycoside genes
length(unique(common_genes$gene)) # 10 occur >=10
#amrfinder$gene<-ifelse(amrfinder$gene %in% count$gene,amrfinder$gene,'other gene')
x<-data.frame(ec_match2$guuid)
x<-left_join(x,amrfinder,by=c("ec_match2.guuid"="guuid"))
x<-select(x,guuid=ec_match2.guuid,gene)
amrfinder<-data.frame(table(x)) %>% spread(gene,Freq)
names(amrfinder)<-str_replace_all(names(amrfinder),'-','_')
amrfinder[,2:ncol(amrfinder)][amrfinder[,2:ncol(amrfinder)]>1]<-1


genes<-names(amrfinder)
genes<-genes[!genes =='guuid']



amrfinder<-left_join(amrfinder,ec_match2,by=c("guuid"="guuid"))

mlst<-read_tsv('./mlst.tsv',col_names = F) %>% select(guuid=X1,spec=X2,mlst=X3)



big<-c('131','95','69','73')
mlst$mlst<-ifelse(mlst$mlst %in% big, mlst$mlst, 'other')

phylogroups<-read_tsv('./phylogroups.tsv')
mlst<-left_join(mlst,phylogroups,by=c("guuid"))
mlst$mlst<-ifelse(mlst$mlst =='other',mlst$phylogroup,mlst$mlst)
other<-c('U','U/cryptic','cryptic','E','G')
mlst$mlst<-ifelse(mlst$mlst %in% other,'other',mlst$mlst)

mlst<-select(mlst,guuid,mlst)

amrfinder<-left_join(amrfinder,mlst,by=c("guuid"))

out=NULL

names(amrfinder)<-janitor::make_clean_names(names(amrfinder))
genes<-janitor::make_clean_names(genes)


genes<-append(genes,"mlst")
amrfinder$mlst<-factor(amrfinder$mlst,levels=c('B2','131','95','73','69','A','B1','C','D','F','other'))
(Y <- with(amrfinder, Surv(cotrim_lower, cotrim_upper, type = "interval2")))
x<-amrfinder$guuid
amrfinder<-select(amrfinder,-guuid,-cotrim_lower,-cotrim_upper)

a<-amrfinder[,grep("^V",colnames(amrfinder))]
for(gene in genes){
  print(which(genes ==gene))
  
  print(gene)
  
  b<-select(amrfinder,gene)
  
  #c<-cbind(b,a)
  c<-b
  t<-as.formula(paste("Y", paste(names(c), collapse = " + "), sep=' ~ '))
  
  
  
  m<-broom::tidy(survreg(t ,data=c,dist="gaussian"),conf.int=T)
  
  m<-if(gene=='mlst'){m[2:11,]}else{m[2,]}    
  out=rbind(out,data.frame(m))
}
include<-out
out_save<-out

amrfinder_save<-amrfinder
big<-c('st131','st95','st73','st69','other')
include<-filter(include,!gene %in% big)

amrfinder<-amrfinder[names(amrfinder) %in% include$term | names(amrfinder) %in% names(a) ]
t<-as.formula(paste("Y", paste(names(amrfinder), collapse = " + "), sep=' ~ '))
model<-survreg(t, data=amrfinder,dist="gaussian")
out2<-broom::tidy(survreg(t, data=amrfinder,dist="gaussian"),conf.int=T)
out2$which<-"multivariable no population"

amrfinder<-amrfinder_save
include<-out
amrfinder<-amrfinder[names(amrfinder) %in% include$term | names(amrfinder) %in% names(a) | names(amrfinder) =='mlst' ]
t<-as.formula(paste("Y", paste(names(amrfinder), collapse = " + "), sep=' ~ '))
model<-survreg(t, data=amrfinder,dist="gaussian")
out3<-broom::tidy(survreg(t, data=amrfinder,dist="gaussian"),conf.int=T)
out3$which<-"multivariable population"

cotrim_r_genes<-filter(out3,2^(estimate) >1) %>% filter(!grepl('mlst',term))
cotrim_line<-2^(log2(4) - model$coefficients[1])
cotrim_line_genes<-filter(out3,2^estimate>cotrim_line) %>% filter(!grepl('mlst',term))
cotrim_all_genes<-filter(out3,!grepl('mlst',term))

out$which<-"univariable"
out<-rbind(out,out2,out3)

count$gene<-janitor::make_clean_names(count$gene)

out<-filter(out,term %in% count$gene| grepl('mlst',term) | term=='other_gene')


ceph_genes<-read_tsv('amrfinder.tsv') %>% filter(grepl('TRIMETHOPRIM',Subclass)| grepl('SULFONAMIDE',Subclass))
ceph_genes<-unique(ceph_genes$`Gene symbol`)
ceph_genes<-janitor::make_clean_names(ceph_genes)

a<-ifelse(out$term %in% ceph_genes,'red','black')

out<-filter(out,!grepl('blaEC',gene))
out$estimate<-as.numeric(out$estimate)
out$conf.high<-as.numeric(out$conf.high)
out$conf.low<-as.numeric(out$conf.low)

X<-c("mlstB1"  ,     "mlstA"    ,   "mlstC"  ,      "mlstD" ,       "mlstE",       
     "mlstF",        "mlstG" ,       "mlstother") 
out$term<-ifelse(out$term %in% X, str_replace(out$term,'mlst','phylogroup '),ifelse(grepl('mlst',out$term),str_replace(out$term,'mlst','ST '),out$term))


out<- out %>% 
  add_row(term = "GAP") %>% 
  add_row(term = "GAP2")

out$estimate<-2^out$estimate
out$conf.high<-2^out$conf.high
out$conf.low<-2^out$conf.low


out$which<-ifelse(grepl("GAP",out$term),'univariable',out$which)
out$which<-factor(out$which,levels=c("univariable","multivariable no population","multivariable population"))

amrfinder<-amrfinder_save1
m<-select(mlst,guuid,gene=mlst) %>% filter(!is.na(gene))
amrfinder<-rbind(amrfinder,m)
count<-amrfinder %>% group_by(gene) %>% distinct(guuid) %>%  count()
count$gene<-ifelse(str_length(count$gene) >3,janitor::make_clean_names(count$gene),count$gene)
count$gene <- case_when(count$gene == '131' ~ 'ST 131',
                        count$gene == '69' ~ 'ST 69',
                        count$gene == '95' ~ 'ST 95',
                        count$gene == '73' ~ 'ST 73',
                        count$gene == 'A' ~ 'phylogroup A',
                        count$gene == 'B1' ~ 'phylogroup B1',
                        count$gene == 'B2' ~ 'phylogroup B2',
                        count$gene == 'C' ~ 'phylogroup C',
                        count$gene == 'D' ~ 'phylogroup D',
                        count$gene == 'F' ~ 'phylogroup F',
                        count$gene == 'other' ~ 'phylogroup other',
                        TRUE ~ count$gene)


out<-left_join(out,count,by=c("term"="gene"))
out$term<-paste0(out$term, ' n=(',out$n,')')


out$term<-str_replace_all(out$term," n=[(]NA[)]",'')
out$term<-factor(out$term, levels=c("dfr_a1 n=(195)","dfr_a12 n=(44)","dfr_a14 n=(117)","dfr_a17 n=(331)",
                                    "dfr_a36 n=(16)","dfr_a5 n=(124)","dfr_a7 n=(61)","dfr_a8 n=(10)",
                                    "other_gene n=(19)","sul1 n=(613)","sul2 n=(785)","sul3 n=(12)","GAP","GAP2",
                                    "ST 131 n=(646)","ST 95 n=(463)","ST 73 n=(683)","ST 69 n=(413)",
                                    "phylogroup A n=(246)","phylogroup B1 n=(267)","phylogroup C n=(101)","phylogroup D n=(285)",
                                    "phylogroup F n=(136)","phylogroup other n=(118)"))



cotrim_plot<-ggplot(out) +
  aes(x=term,y=estimate,color=which,group=which) +
  geom_point(position = position_dodge(0.9)) + geom_errorbar(data=out,aes(ymin=conf.low,ymax=conf.high),position=position_dodge(width=0.9),width=0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90),panel.grid.major.x = element_blank()) +
  scale_y_log10() +
  xlab('Gene') + ylab('Fold-change in MIC') + 
  labs(color="Model") +
  ggtitle("Co-trimoxazole") + geom_hline(yintercept =2^(log2(4) - model$coefficients[1]),linetype="dashed",color="dark red") +
  geom_hline(yintercept=1)+
  scale_x_discrete(labels = function(x) ifelse(x == "GAP" | x=="GAP2", "", x)) +
  scale_color_viridis_d()

cotrim<-out
cotrim$drug<-"cotrim"


amrfinder<-amrfinder_save
predictors <- out_save$term
predictors<-ifelse(grepl("mlst",predictors),"mlst",predictors)
predictors<-unique(predictors)
predictors<-predictors[!is.na(predictors)  & ! predictors =="(Intercept)" & ! predictors == "Log(scale)"]

out=NULL



for (excluded_gene in predictors) {
  print(excluded_gene)
  
  other_predictors <- setdiff(predictors, excluded_gene)
  
  for (current_gene in other_predictors) {
    
    interaction_term <- paste(excluded_gene, current_gene, sep="*")
    
    
    formula_str <- as.formula(paste( "Y ~", paste(setdiff(other_predictors,current_gene), collapse = "+"), "+", interaction_term))
    summary_model <- broom::tidy(survreg(formula_str, data = amrfinder, dist = "gaussian"))
    
  
    
    out <- rbind(out,data.frame(summary_model))
    
  }
}




replace_mlst_patterns <- function(texts) {
  # Replace mlst* before the colon
  texts <- gsub("mlst[A-Za-z0-9_]*:", "mlst:", texts)
  # Replace mlst* after the colon
  texts <- gsub(":[A-Za-z0-9_]*mlst[A-Za-z0-9_]*", ":mlst", texts)
  return(texts)
}



interaction_terms<-filter(out,p.value < 0.01) %>% filter(grepl(':',term))
interaction_terms$term<-replace_mlst_patterns(interaction_terms$term)
interaction_terms<-unique(interaction_terms$term)
interaction_terms

model<-survreg(Y ~ dfr_a1 + dfr_a12 + dfr_a14 + dfr_a17 + dfr_a19_other + dfr_a36 + 
                 dfr_a5 + dfr_a7 + dfr_a8 + other_gene + sul1 + sul2 + sul3 + 
                 mlst 
               ,data=amrfinder,dist="gaussian")

elimination<-  step(model, direction="backward")

selected_vars <- names(coef(elimination))
selected_vars<-selected_vars[!selected_vars=="(Intercept)"]


interaction_string <- paste(interaction_terms, collapse = " + ")
new_formula <- as.formula(paste("Y ~ ", paste(selected_vars, collapse = " + "), " + ", interaction_string))
model<-survreg(new_formula,data=amrfinder,dist="gaussian")
elimination<-step(model,direction="backward",k=3.8)

###leave one out cross validation


ceph_gene<-read_tsv('amrfinder.tsv') %>% filter(grepl('TRIMETHOPRIM',Subclass) | grepl('SULFONAMIDE',Subclass)) %>% filter(Scope=='core')


out<-NULL
amrfinder$guuid<-x
for (i in 1:nrow(amrfinder)) {
  print(i)
  model<-survreg( Y[-i,] ~ dfr_a1 + dfr_a12 + dfr_a14 + dfr_a17 + 
                    dfr_a19_other + dfr_a36 + dfr_a5 + dfr_a7 + dfr_a8 + other_gene + 
                    sul1 + sul2 + sul3 + dfr_a1:sul1 + dfr_a1:sul2 + dfr_a12:dfr_a17 + 
                    dfr_a14:dfr_a17 + dfr_a14:sul2 + dfr_a17:dfr_a7 + dfr_a17:sul2 + 
                    dfr_a5:sul1 + dfr_a5:sul2 + sul1:sul2 + sul1:sul3 + sul1:mlst, data = amrfinder[-i,], 
                  dist = "gaussian")
  
  prediction <- predict(model, newdata = amrfinder[i, ])
  correct<-ifelse(prediction >= as.numeric(Y[i,1]) & prediction <= as.numeric(Y[i,2]),1,0)
  
  out<-rbind(out,data.frame(correct,prediction,as.numeric(Y[i,1]), as.numeric(Y[i,2]),amrfinder$guuid[i]))
}


out$amrfinder_prediction<-ifelse(out$amrfinder.guuid.i. %in% ceph_gene$Name,'R','S')
out$amrfinder_prediction<-as.factor(out$amrfinder_prediction)
#breakpoint >4 = log2(4) = 2
out$binary<-ifelse(out$as.numeric.Y.i..2.. >2,'R','S')
out$binary_prediction<-ifelse(out$prediction>2,'R','S')
out$binary<-as.factor(out$binary)
out$binary_prediction<-as.factor(out$binary_prediction)

extend_margin <- function(out, n){
  out$margin <- NA # Initialize margin column
  for(i in 0:n){
    condition <- out$prediction >= (out$as.numeric.Y.i..1.. - i) & 
      out$prediction <= (out$as.numeric.Y.i..2.. + i)
    out$margin[is.na(out$margin) & condition] <- i
  }
  return(out)
}


out2<-extend_margin(out,10)
out2$margin<-ifelse(out2$prediction < out2$as.numeric.Y.i..1..,out2$margin *-1,out2$margin)

out2$which<- ifelse(out2$binary_prediction==out2$binary & out2$binary=='S', 'Correct S',
                    ifelse(out2$binary_prediction==out2$binary & out2$binary =='R','Correct R',
                           ifelse(out2$binary_prediction=='S' & out2$binary =='R','VM',
                                  ifelse(out2$binary_prediction =='R' & out2$binary =='S','M',NA))))

out$within_1<-ifelse(out$prediction >= (out$as.numeric.Y.i..1.. -1) & out$prediction <= (out$as.numeric.Y.i..2.. +1),1,0)
out$within_2<-ifelse(out$prediction >= (out$as.numeric.Y.i..1.. -2) & out$prediction <= (out$as.numeric.Y.i..2.. +2),1,0)


caret::confusionMatrix(out$binary_prediction,out$binary)
table(out$correct)
prop.test(2635,(2635+231))
2635+231


table(out$within_1)
prop.test(2742,(2742+124))
2742+124

cm<-caret::confusionMatrix(as.factor(out$binary_prediction),reference=as.factor(out$binary),positive="R")
tn<-cm$table[2,2]
fn<-cm$table[2,1]
tp<-cm$table[1,1]
fp<-cm$table[1,2]
#concordance
prop.test(x=(tp+tn),n=(tp+tn+fp+fn))
tp+tn
tp+tn+fp+fn
#sens
prop.test(x=(tp),n=(tp+fn))
tp
tp+fn
#spec
prop.test(x=(tn),n=(tn+fp))
tn
tn+fp

#Maj
prop.test(x=(fp),n=(tn+fp))
fp
tn+fp
#Vmj
prop.test(x=(fn),n=(tp+fn))
fn
tp + fn

#PPV
prop.test(x=(tp),n=(tp+fp))
tp
tp+fp

#NPV
prop.test(x=(tn),n=(tn+fn))
tn
tn+fn


colors<-c("#09820B" ,"#03FF25" ,"#FFA703", "#ED3A37")
names(colors)<-c("Correct R","Correct S","M","VM")
cotrim_residual<-ggplot(out2) +
  aes(x=margin, fill=which) +
  geom_bar(stat="count") + theme_minimal() +
  #scale_fill_manual(values=colors) + 
  scale_x_continuous(breaks = seq(min(out2$margin),max(out2$margin),1)) + 
  xlab("True - Predicted MIC") +
  labs(fill="") +
  ggtitle("Co-trimoxazole") +
  scale_fill_viridis_d()

out2$cotrim_pred<-
  case_when(out2$prediction <=0 ~ '<=1/19',
            out2$prediction >0 & out2$prediction <=1 ~ '2/38',
            out2$prediction >1 & out2$prediction <=2 ~ '4/76',
            out2$prediction >2 ~ '>4/76')

em<-select(ec_match,guuid,Cotrim)
out2<-left_join(out2,em,by=c("amrfinder.guuid.i."="guuid"))

out2$cotrim_pred<-factor(out2$cotrim_pred,levels=c('<=1/19','2/38','4/76','>4/76'))
out2$Cotrim<-factor(out2$Cotrim,levels=c('<=1/19','2/38','4/76','>4/76'))

cotrim_sup_dist<-ggplot() +
  geom_histogram(data=out2,aes(x=cotrim_pred,fill=Cotrim),stat="count",) + scale_fill_viridis_d() +
  theme_minimal() + xlab("Predicted co-amoxiclav MIC") + labs(fill="Measured co-amoxiclav\nMIC") +
  ggplot() +
  geom_histogram(data=out2,aes(x=Cotrim,fill=cotrim_pred),stat="count") + scale_fill_viridis_d() +
  theme_minimal() + xlab("Measured co-amoxiclav MIC") + labs(fill="Predicted co-amoxiclav\nMIC") +
  plot_annotation(tag_levels = "A")





amrfinder<-read_tsv('amrfinder.tsv')
cotrim<-select(ec_match,guuid,Cotrim) %>% filter(guuid %in% amrfinder$Name)


amrfinder<-amrfinder%>% filter(grepl('TRIMETHOPRIM',Subclass )|grepl('SULFONAMIDE',Subclass)) %>%  select(guuid=Name,gene=`Gene symbol`)
length(unique(amrfinder$gene))


count<-amrfinder %>% filter(str_length(guuid) ==36) %>%  group_by(gene) %>% count() %>% filter(n>=10)



deal_with_rare_genes <- function(gene_name) {
  if (startsWith(gene_name, "ampC")) {
    gsub("^ampC_.*", "ampC_other", gene_name)
  } else if (!grepl("_", gene_name)) {
    # Append '_other' at the end
    paste0(gene_name, "_other")
  } else {
    # Replace after the last '_'
    gsub("_[^_]+$", "_other", gene_name)
  }
}

amrfinder$gene <- ifelse(amrfinder$gene %in% count$gene,
                         amrfinder$gene,
                         sapply(amrfinder$gene, deal_with_rare_genes))

singletons<-amrfinder %>%  group_by(gene) %>% count() %>% filter(n<5)

amrfinder$gene<-ifelse(amrfinder$gene %in% singletons$gene,"other gene",amrfinder$gene)


cotrim_r_genes<-filter(cotrim_r_genes,p.value<0.05)
amrfinder$gene<-janitor::make_clean_names(amrfinder$gene,allow_dupes = T)
amrfinder_r<-filter(amrfinder,gene %in% cotrim_r_genes$term)
amrfinder_r<- amrfinder_r %>% group_by(guuid) %>% count()

cotrim$Cotrim<-factor(cotrim$Cotrim,levels=c('<=1/19','2/38','4/76','>4/76'))
cotrim<-filter(cotrim,!is.na(Cotrim))
amrfinder_r<-left_join(cotrim,amrfinder_r,by=c("guuid"="guuid"))
amrfinder_r$n<-ifelse(is.na(amrfinder_r$n),0,amrfinder_r$n)
amrfinder_r<-filter(amrfinder_r,!is.na(Cotrim))

r<-c('>4/76')

amrfinder_r$R<-ifelse(amrfinder_r$Cotrim %in% r,'R','S')
amrfinder_r$R<-factor(amrfinder_r$R,levels=c('S','R'))
cotrim_gene_dist<-ggplot(amrfinder_r) +
  aes(x=n,fill=R) +
  geom_histogram(bins=length(unique(amrfinder_r$n))) +
  facet_wrap(~Cotrim,nrow = 1) + theme_minimal() +
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),strip.text = element_blank()) +
  xlab("Number of ARGs/mutations")+
  scale_fill_viridis_d() +
  labs(fill="")

ng<-amrfinder_r %>% 
  group_by(Cotrim) %>% 
  summarise(prop_no_gene = sum(n==0)/n(),
            lower=prop.test(sum(n==0),sum(n()))$conf.int[1],
            upper=prop.test(sum(n==0),sum(n()))$conf.int[2])


cotrim_ng_plot<-ggplot(ng) +
  aes(x=Cotrim,y=prop_no_gene) +
  geom_point() +
  geom_errorbar(data=ng,aes(ymin=lower,ymax=upper),width=0) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank()) +
  ylab("Proportion with no ARG/mutation") +
  xlab("MIC")


cotrim_eucast_plot<-ng_plot / cotrim_gene_dist  +
  plot_annotation(title = "Co-trimoxazole")



############ampicillin ecoli################
ec_match2<-select(ec_match,guuid,Ampicillin_lower,Ampicillin_upper) %>% filter(!is.na(Ampicillin_upper))
amrfinder<-read_tsv('amrfinder.tsv')
tem_promotor_mutations<-filter(amrfinder,grepl('TEMp',`Gene symbol`))
ec_match2<-filter(ec_match2,guuid %in% amrfinder$Name)
amrfinder<-amrfinder%>% filter(Class=='BETA-LACTAM' ) %>%  select(guuid=Name,gene=`Gene symbol`)
length(unique(amrfinder$gene))
amrfinder$gene<-str_replace_all(amrfinder$gene,'-','_')
amrfinder$gene<-ifelse(grepl('blaTEM_',amrfinder$gene) & amrfinder$guuid %in% tem_promotor_mutations$Name,
                       paste0(amrfinder$gene,'_',tem_promotor_mutations$`Gene symbol`[match(amrfinder$guuid,tem_promotor_mutations$Name)]),amrfinder$gene)
amrfinder <- amrfinder %>%
  filter(!grepl("(?<!\\w)blaTEMp_", gene, perl = TRUE))


ceph<-read_tsv('amrfinder.tsv')  %>% filter(Name %in% ec_match2$guuid) %>% filter(grepl('BETA-LACTAM',Class )) %>%  select(guuid=Name,gene=`Gene symbol`) 
ceph$gene<-str_replace_all(ceph$gene,'-','_')


count<-amrfinder %>% filter(str_length(guuid) ==36) %>%  group_by(gene) %>% count() %>% filter(n>=10)
singletons<-amrfinder %>% filter(str_length(guuid) ==36) %>%  group_by(gene) %>% count() %>% filter(n==1)


deal_with_rare_genes <- function(gene_name) {
  if (startsWith(gene_name, "ampC")) {
    gsub("^ampC_.*", "ampC_other", gene_name)
  } else if (!grepl("_", gene_name)) {
    # Append '_other' at the end
    paste0(gene_name, "_other")
  } else {
    # Replace after the last '_'
    gsub("_[^_]+$", "_other", gene_name)
  }
}

amrfinder$gene <- ifelse(amrfinder$gene %in% count$gene,
                         amrfinder$gene,
                                ifelse(!amrfinder$gene %in% singletons$gene,sapply(amrfinder$gene, deal_with_rare_genes),amrfinder$gene))

amrfinder$gene<-ifelse(amrfinder$gene %in% singletons$gene,"other gene",amrfinder$gene)
amrfinder_save1<-amrfinder


common_genes<-filter(ceph,gene %in% count$gene)
length(unique(ceph$gene)) # 49 amp genes
length(unique(common_genes$gene)) # 14 occur >=10

x<-data.frame(ec_match2$guuid)
x<-left_join(x,amrfinder,by=c("ec_match2.guuid"="guuid"))
x<-select(x,guuid=ec_match2.guuid,gene)
amrfinder<-data.frame(table(x)) %>% spread(gene,Freq)
names(amrfinder)<-str_replace_all(names(amrfinder),'-','_')
amrfinder[,2:ncol(amrfinder)][amrfinder[,2:ncol(amrfinder)]>1]<-1

genes<-names(amrfinder)
genes<-genes[!genes =='guuid']

amrfinder<-left_join(amrfinder,ec_match2,by=c("guuid"="guuid"))

mlst<-read_tsv('./mlst.tsv',col_names = F) %>% select(guuid=X1,spec=X2,mlst=X3)


big<-c('131','95','69','73')
mlst$mlst<-ifelse(mlst$mlst %in% big, mlst$mlst, 'other')

phylogroups<-read_tsv('./phylogroups.tsv')
mlst<-left_join(mlst,phylogroups,by=c("guuid"))
mlst$mlst<-ifelse(mlst$mlst =='other',mlst$phylogroup,mlst$mlst)
other<-c('U','U/cryptic','cryptic','E','G')
mlst$mlst<-ifelse(mlst$mlst %in% other,'other',mlst$mlst)

mlst<-select(mlst,guuid,mlst)

amrfinder<-left_join(amrfinder,mlst,by=c("guuid"))

out=NULL

names(amrfinder)<-janitor::make_clean_names(names(amrfinder))
genes<-janitor::make_clean_names(genes)


genes<-append(genes,"mlst")
amrfinder$mlst<-factor(amrfinder$mlst,levels=c('B2','131','95','73','69','A','B1','C','D','F','other'))
(Y <- with(amrfinder, Surv(ampicillin_lower, ampicillin_upper, type = "interval2")))
x<-amrfinder$guuid
amrfinder<-select(amrfinder,-guuid,-ampicillin_lower,-ampicillin_upper)

a<-amrfinder[,grep("^V",colnames(amrfinder))]
for(gene in genes){
  print(which(genes ==gene))
  
  print(gene)
  
  b<-select(amrfinder,gene)
  
  #c<-cbind(b,a)
  c<-b
  t<-as.formula(paste("Y", paste(names(c), collapse = " + "), sep=' ~ '))
  
  
  
  m<-broom::tidy(survreg(t ,data=c,dist="gaussian"),conf.int=T)
  
  m<-if(gene=='mlst'){m[2:11,]}else{m[2,]}    
  out=rbind(out,data.frame(m))
}
include<-out
out_save<-out
amrfinder_save<-amrfinder
big<-c('st131','st95','st73','st69','other')
include<-filter(include,!gene %in% big)

amrfinder<-amrfinder[names(amrfinder) %in% include$term | names(amrfinder) %in% names(a) ]
t<-as.formula(paste("Y", paste(names(amrfinder), collapse = " + "), sep=' ~ '))
model<-survreg(t, data=amrfinder,dist="gaussian")
out2<-broom::tidy(survreg(t, data=amrfinder,dist="gaussian"),conf.int=T)
out2$which<-"multivariable no population"

amrfinder<-amrfinder_save
include<-out
amrfinder<-amrfinder[names(amrfinder) %in% include$term | names(amrfinder) %in% names(a) | names(amrfinder) =='mlst' ]
t<-as.formula(paste("Y", paste(names(amrfinder), collapse = " + "), sep=' ~ '))
model<-survreg(t, data=amrfinder,dist="gaussian")
out3<-broom::tidy(survreg(t, data=amrfinder,dist="gaussian"),conf.int=T)
out3$which<-"multivariable population"

amp_r_genes<-filter(out3,2^(estimate) >1) %>% filter(!grepl('mlst',term))
amp_line<-2^(log2(8) - model$coefficients[1])
amp_line_genes<-filter(out3,2^estimate >amp_line) %>% filter(!grepl('mlst',term))
amp_all_genes<-filter(out3,!grepl('mlst',term))

out$which<-"univariable"
out<-rbind(out,out2,out3)

count$gene<-janitor::make_clean_names(count$gene)

out<-filter(out,term %in% count$gene| grepl('mlst',term) | term=='other_gene')



ceph_genes<-read_tsv('amrfinder.tsv') %>% filter(grepl('CEPHALOSPORIN',Subclass))
ceph_genes<-unique(ceph_genes$`Gene symbol`)
ceph_genes<-janitor::make_clean_names(ceph_genes)

a<-ifelse(out$term %in% ceph_genes,'red','black')

out<-filter(out,!grepl('blaEC',gene))
out$estimate<-as.numeric(out$estimate)
out$conf.high<-as.numeric(out$conf.high)
out$conf.low<-as.numeric(out$conf.low)

out$term <- case_when(out$term == "bla_tem_1_bla_te_mp_c32t" ~ "bla_tem_1 + blaTEMp_c32t",
                      out$term == "bla_tem_1_bla_te_mp_g162t" ~ "bla_tem_1 + blaTEMp_g162t",
                      TRUE ~ out$term)

X<-c("mlstB1"  ,     "mlstA"    ,   "mlstC"  ,      "mlstD" ,       "mlstE",       
     "mlstF",        "mlstG" ,       "mlstother") 
out$term<-ifelse(out$term %in% X, str_replace(out$term,'mlst','phylogroup '),ifelse(grepl('mlst',out$term),str_replace(out$term,'mlst','ST '),out$term))



out<- out %>% 
  add_row(term = "GAP") %>% 
  add_row(term = "GAP2")

out$estimate<-2^out$estimate
out$conf.high<-2^out$conf.high
out$conf.low<-2^out$conf.low

out$which<-ifelse(grepl("GAP",out$term),'univariable',out$which)
out$which<-factor(out$which,levels=c("univariable","multivariable no population","multivariable population"))

amrfinder_save1$gene <- case_when(amrfinder_save1$gene == "blaTEM_1_blaTEMp_C32T" ~ "bla_tem_1 + blaTEMp_c32t",
                                  amrfinder_save1$gene== "blaTEM_1_blaTEMp_G162T" ~ "bla_tem_1 + blaTEMp_g162t",
                      TRUE ~ amrfinder_save1$gene)

amrfinder<-amrfinder_save1
m<-select(mlst,guuid,gene=mlst) %>% filter(!is.na(gene))
amrfinder<-rbind(amrfinder,m)
count<-amrfinder %>% group_by(gene) %>% distinct(guuid) %>%  count()
count$gene<-ifelse(str_length(count$gene) >3,janitor::make_clean_names(count$gene),count$gene)
count$gene <- case_when(count$gene == '131' ~ 'ST 131',
                        count$gene == '69' ~ 'ST 69',
                        count$gene == '95' ~ 'ST 95',
                        count$gene == '73' ~ 'ST 73',
                        count$gene == 'A' ~ 'phylogroup A',
                        count$gene == 'B1' ~ 'phylogroup B1',
                        count$gene == 'B2' ~ 'phylogroup B2',
                        count$gene == 'C' ~ 'phylogroup C',
                        count$gene == 'D' ~ 'phylogroup D',
                        count$gene == 'F' ~ 'phylogroup F',
                        count$gene == 'other' ~ 'phylogroup other',
                        TRUE ~ count$gene)

count$gene <- case_when(count$gene== "bla_tem_1_bla_te_mp_c32t" ~ "bla_tem_1 + blaTEMp_c32t",
                        count$gene== "bla_tem_1_bla_te_mp_g162t" ~ "bla_tem_1 + blaTEMp_g162t",
                      TRUE ~ count$gene)



out<-left_join(out,count,by=c("term"="gene"))
out$term<-paste0(out$term, ' n=(',out$n,')')


out$term<-str_replace_all(out$term," n=[(]NA[)]",'')
out$term<-factor(out$term, levels=c("amp_c_c_11t n=(15)","amp_c_c_42t n=(10)","amp_c_t_32a n=(34)","bla_ctx_m_14 n=(16)",     
                                    "bla_ctx_m_15 n=(171)","bla_ctx_m_27 n=(32)","bla_ec n=(2216)","bla_ec_5 n=(667)",        
                                    "bla_oxa_1 n=(163)","bla_shv_1 n=(56)","bla_tem n=(34)","bla_tem_1 n=(1033)",  
                                    "bla_tem_1 + blaTEMp_c32t n=(90)" , "bla_tem_1 + blaTEMp_g162t n=(15)" ,"bla_tem_40 n=(13)","other_gene n=(22)","GAP","GAP2",       
                                    "ST 131 n=(646)","ST 95 n=(463)","ST 73 n=(683)","ST 69 n=(413)",
                                    "phylogroup A n=(246)","phylogroup B1 n=(267)","phylogroup C n=(101)","phylogroup D n=(285)",
                                    "phylogroup F n=(136)","phylogroup other n=(118)"))




amp_plot<-ggplot(out) +
  aes(x=term,y=estimate,color=which,group=which) +
  geom_point(position = position_dodge(0.9)) + geom_errorbar(data=out,aes(ymin=conf.low,ymax=conf.high),position=position_dodge(width=0.9),width=0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90),panel.grid.major.x = element_blank()) +
  scale_y_log10() +
  xlab('Gene') + ylab('Fold-change in MIC') + 
  labs(color="Model") +
  ggtitle("Ampicillin") + geom_hline(yintercept =2^(log2(8) - model$coefficients[1]),linetype="dashed",color="dark red") +
  geom_hline(yintercept=1)+
  scale_x_discrete(labels = function(x) ifelse(x == "GAP" | x=="GAP2", "", x)) +
  scale_color_viridis_d()


ampicillin<-out
ampicillin$drug<-"ampicillin"


amrfinder<-amrfinder_save
predictors <- out_save$term
predictors<-ifelse(grepl("mlst",predictors),"mlst",predictors)
predictors<-unique(predictors)
predictors<-predictors[!is.na(predictors)  & ! predictors =="(Intercept)" & ! predictors == "Log(scale)"]

out=NULL



for (excluded_gene in predictors) {
  print(excluded_gene)
  other_predictors <- setdiff(predictors, excluded_gene)
  
  for (current_gene in other_predictors) {
    interaction_term <- paste(excluded_gene, current_gene, sep="*")
    
    formula_str <- as.formula(paste("Y ~", paste(setdiff(other_predictors,current_gene), collapse = "+"), "+", interaction_term))
    summary_model <- broom::tidy(survreg(formula_str, data = amrfinder, dist = "gaussian"))
    
    
    

    out <- rbind(out,data.frame(summary_model))
    
  }
}




replace_mlst_patterns <- function(texts) {
  texts <- gsub("mlst[A-Za-z0-9_]*:", "mlst:", texts)
  texts <- gsub(":[A-Za-z0-9_]*mlst[A-Za-z0-9_]*", ":mlst", texts)
  return(texts)
}



interaction_terms<-filter(out,p.value < 0.01) %>% filter(grepl(':',term))
interaction_terms$term<-replace_mlst_patterns(interaction_terms$term)
interaction_terms<-unique(interaction_terms$term)
interaction_terms

model<-survreg(Y ~ amp_c_c_11t + amp_c_c_42t + amp_c_other + amp_c_t_32a + bla_carb_other + 
                 bla_cmy_other + bla_ctx_m_14 + bla_ctx_m_15 + bla_ctx_m_27 + 
                 bla_ctx_m_other + bla_ctx_other + bla_dha_other + bla_ec + 
                 bla_ec_5 + bla_mal_other + bla_oxa_1 + bla_oxa_other + bla_shv_1 + 
                 bla_shv_other + bla_tem + bla_tem_1 + bla_tem_1_bla_te_mp_c32t + 
                 bla_tem_1_bla_te_mp_g162t + bla_tem_1_bla_te_mp_other + bla_tem_30_bla_te_mp_other + 
                 bla_tem_33_bla_te_mp_other + bla_tem_40 + bla_tem_40_bla_te_mp_other + 
                 bla_tem_other + omp_c_other + omp_f_other + other_gene + 
                 mlst  ,data=amrfinder,dist="gaussian")

elimination<-  step(model, direction="backward")

selected_vars <- names(coef(elimination))
selected_vars<-selected_vars[!selected_vars=="(Intercept)"]
selected_vars<-selected_vars[!grepl('mlst',selected_vars)]
selected_vars<-c(selected_vars,"mlst")


interaction_string <- paste(interaction_terms, collapse = " + ")
new_formula <- as.formula(paste("Y ~ ", paste(selected_vars, collapse = " + "), " + ", interaction_string))
model<-survreg(new_formula,data=amrfinder,dist="gaussian")
elimination<-step(model,direction="backward",k=3.8)

###leave one out cross validation


ceph_gene<-read_tsv('amrfinder.tsv') %>% filter(Subclass=='BETA-LACTAM') %>% filter(Scope=='core')


out<-NULL
amrfinder$guuid<-x
for (i in 1:nrow(amrfinder)) {
  print(i)
  model<-survreg( Y[-i,] ~ amp_c_c_11t + amp_c_c_42t + amp_c_other + 
                    amp_c_t_32a + bla_carb_other + bla_cmy_other + bla_ctx_m_14 + 
                    bla_ctx_m_15 + bla_ctx_m_27 + bla_ctx_m_other + bla_ctx_other + 
                    bla_dha_other + bla_mal_other + bla_oxa_1 + bla_shv_1 + bla_shv_other + 
                    bla_tem + bla_tem_1 + bla_tem_1_bla_te_mp_c32t + bla_tem_1_bla_te_mp_g162t + 
                    bla_tem_1_bla_te_mp_other + bla_tem_30_bla_te_mp_other + 
                    bla_tem_33_bla_te_mp_other + bla_tem_40 + bla_tem_40_bla_te_mp_other + 
                    bla_tem_other + other_gene + amp_c_c_42t:bla_tem_1 + amp_c_t_32a:bla_tem_1 + 
                    amp_c_t_32a:bla_tem_30_bla_te_mp_other + bla_cmy_other:bla_oxa_1 + 
                    bla_ctx_m_14:bla_tem_1 + bla_ctx_m_15:bla_oxa_1 + bla_ctx_m_15:bla_tem_1 + 
                    bla_ctx_m_27:bla_tem_1 + bla_ctx_m_27:bla_tem_1_bla_te_mp_c32t + 
                    bla_tem_1:bla_ec + other_gene:bla_ec + bla_tem_1:bla_ec_5 + 
                    bla_oxa_1:bla_tem_1 + bla_shv_1:bla_tem_1 + bla_shv_other:bla_tem_1 + 
                    bla_tem_1:other_gene, data = amrfinder[-i,], 
                  dist = "gaussian")
  
  prediction <- predict(model, newdata = amrfinder[i, ])
  correct<-ifelse(prediction >= as.numeric(Y[i,1]) & prediction <= as.numeric(Y[i,2]),1,0)
  
  out<-rbind(out,data.frame(correct,prediction,as.numeric(Y[i,1]), as.numeric(Y[i,2]),amrfinder$guuid[i]))
}


out$amrfinder_prediction<-ifelse(out$amrfinder.guuid.i. %in% ceph_gene$Name,'R','S')
out$amrfinder_prediction<-as.factor(out$amrfinder_prediction)
out$binary<-ifelse(out$as.numeric.Y.i..2.. >3,'R','S')
out$binary_prediction<-ifelse(out$prediction>3,'R','S')
out$binary<-as.factor(out$binary)
out$binary_prediction<-as.factor(out$binary_prediction)

extend_margin <- function(out, n){
  out$margin <- NA # Initialize margin column
  for(i in 0:n){
    condition <- out$prediction >= (out$as.numeric.Y.i..1.. - i) & 
      out$prediction <= (out$as.numeric.Y.i..2.. + i)
    out$margin[is.na(out$margin) & condition] <- i
  }
  return(out)
}


out2<-extend_margin(out,10)
out2$margin<-ifelse(out2$prediction < out2$as.numeric.Y.i..1..,out2$margin *-1,out2$margin)

out2$which<- ifelse(out2$binary_prediction==out2$binary & out2$binary=='S', 'Correct S',
                    ifelse(out2$binary_prediction==out2$binary & out2$binary =='R','Correct R',
                           ifelse(out2$binary_prediction=='S' & out2$binary =='R','VM',
                                  ifelse(out2$binary_prediction =='R' & out2$binary =='S','M',NA))))

out$within_1<-ifelse(out$prediction >= (out$as.numeric.Y.i..1.. -1) & out$prediction <= (out$as.numeric.Y.i..2.. +1),1,0)
out$within_2<-ifelse(out$prediction >= (out$as.numeric.Y.i..1.. -2) & out$prediction <= (out$as.numeric.Y.i..2.. +2),1,0)


caret::confusionMatrix(out$binary_prediction,out$binary)
table(out$correct)
prop.test(2552,(2552+317))
2552+317


table(out$within_1)
prop.test(2734,(2734+135))
2734+135

cm<-caret::confusionMatrix(as.factor(out$binary_prediction),reference=as.factor(out$binary),positive="R")
tn<-cm$table[2,2]
fn<-cm$table[2,1]
tp<-cm$table[1,1]
fp<-cm$table[1,2]
#concordance
prop.test(x=(tp+tn),n=(tp+tn+fp+fn))
tp+tn
tp+tn+fp+fn
#sens
prop.test(x=(tp),n=(tp+fn))
tp
tp+fn
#spec
prop.test(x=(tn),n=(tn+fp))
tn
tn+fp

#Maj
prop.test(x=(fp),n=(tn+fp))
fp
tn+fp
#Vmj
prop.test(x=(fn),n=(tp+fn))
fn
tp + fn

#PPV
prop.test(x=(tp),n=(tp+fp))
tp
tp+fp

#NPV
prop.test(x=(tn),n=(tn+fn))
tn
tn+fn

colors<-c("#09820B" ,"#03FF25" ,"#FFA703", "#ED3A37")
names(colors)<-c("Correct R","Correct S","M","VM")
amp_residual<-ggplot(out2) +
  aes(x=margin, fill=which) +
  geom_bar(stat="count") + theme_minimal() +
  #scale_fill_manual(values=colors) + 
  scale_x_continuous(breaks = seq(min(out2$margin),max(out2$margin),1)) + 
  xlab("True - Predicted MIC") +
  labs(fill="") +
  ggtitle("Ampicillin") +
  scale_fill_viridis_d()



amrfinder<-read_tsv('amrfinder.tsv') %>% filter(grepl('BETA-LACTAM',Class))%>%  select(guuid=Name,gene=`Gene symbol`)

amp_r_genes<-filter(amp_r_genes,p.value<0.05)
amrfinder$gene<-janitor::make_clean_names(amrfinder$gene,allow_dupes = T)
amrfinder_r<-filter(amrfinder,gene %in% amp_r_genes$term)
amrfinder_r<- amrfinder_r %>% group_by(guuid) %>% count()
amp<-select(ec_match,guuid,Ampicillin) 
amp$Ampicillin<-factor(amp$Ampicillin,levels=c('<=2','4','8','>8'))
amrfinder_r<-left_join(amp,amrfinder_r,by=c("guuid"="guuid"))
amrfinder_r$n<-ifelse(is.na(amrfinder_r$n),0,amrfinder_r$n)
amrfinder_r<-filter(amrfinder_r,!is.na(Ampicillin))

r<-c('>8')

amrfinder_r$R<-ifelse(amrfinder_r$Ampicillin %in% r,'R','S')
amp_gene_dist<-ggplot(amrfinder_r) +
  aes(x=n,fill=R) +
  geom_histogram(bins=4) +
  facet_wrap(~Ampicillin,nrow = 1) + theme_minimal() +
  scale_fill_brewer(palette = "Dark2") +
  labs(fill="")


amrfinder<-read_tsv('amrfinder.tsv')
amp<-select(ec_match,guuid,Ampicillin) %>% filter(guuid %in% amrfinder$Name)
tem_promotor_mutations<-filter(amrfinder,grepl('TEMp',`Gene symbol`))
ec_match2<-filter(ec_match2,guuid %in% amrfinder$Name)
amrfinder<-amrfinder%>% filter(Class=='BETA-LACTAM' ) %>%  select(guuid=Name,gene=`Gene symbol`)
length(unique(amrfinder$gene))
amrfinder$gene<-str_replace_all(amrfinder$gene,'-','_')
amrfinder_keep<-amrfinder
amrfinder$gene<-ifelse(grepl('blaTEM_',amrfinder$gene) & amrfinder$guuid %in% tem_promotor_mutations$Name,
                       paste0(amrfinder$gene,'_',tem_promotor_mutations$`Gene symbol`[match(amrfinder$guuid,tem_promotor_mutations$Name)]),amrfinder$gene)
amrfinder <- amrfinder %>%
  filter(!grepl("(?<!\\w)blaTEMp_", gene, perl = TRUE))
coamox<-read_tsv('amrfinder.tsv')  %>% filter(Name %in% ec_match2$guuid) %>% filter(grepl('BETA-LACTAM',Class )) %>%  select(guuid=Name,gene=`Gene symbol`) 
coamox$gene<-str_replace_all(coamox$gene,'-','_')

count<-amrfinder %>% filter(str_length(guuid) ==36) %>%  group_by(gene) %>% count() %>% filter(n>=10)


deal_with_rare_genes <- function(gene_name) {
  if (startsWith(gene_name, "ampC")) {
    gsub("^ampC_.*", "ampC_other", gene_name)
  } else if (!grepl("_", gene_name)) {
    paste0(gene_name, "_other")
  } else {
    gsub("_[^_]+$", "_other", gene_name)
  }
}

amrfinder$gene <- ifelse(amrfinder$gene %in% count$gene,
                         amrfinder$gene,
                         sapply(amrfinder$gene, deal_with_rare_genes))

singletons<-amrfinder %>%  group_by(gene) %>% count() %>% filter(n<5)

amrfinder$gene<-ifelse(amrfinder$gene %in% singletons$gene,"other gene",amrfinder$gene)


amp_r_genes<-filter(amp_r_genes,p.value<0.05)
amrfinder$gene<-janitor::make_clean_names(amrfinder$gene,allow_dupes = T)
amrfinder_r<-filter(amrfinder,gene %in% amp_r_genes$term)
amrfinder_r<- amrfinder_r %>% group_by(guuid) %>% count()

amp$Ampicillin<-factor(amp$Ampicillin,levels=c('<=2','4','8','>8'))
amp<-filter(amp,!is.na(Ampicillin))
amrfinder_r<-left_join(amp,amrfinder_r,by=c("guuid"="guuid"))
amrfinder_r$n<-ifelse(is.na(amrfinder_r$n),0,amrfinder_r$n)
amrfinder_r<-filter(amrfinder_r,!is.na(Ampicillin))

r<-c('>8')

amrfinder_r$R<-ifelse(amrfinder_r$Ampicillin %in% r,'R','S')
amrfinder_r$R<-factor(amrfinder_r$R,levels=c('S','R'))
amp_gene_dist<-ggplot(amrfinder_r) +
  aes(x=n,fill=R) +
  geom_histogram(bins=length(unique(amrfinder_r$n))) +
  facet_wrap(~Ampicillin,nrow = 1) + theme_minimal() +
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),strip.text = element_blank()) +
  xlab("Number of ARGs/mutations")+
  scale_fill_viridis_d() +
  labs(fill="")

ng<-amrfinder_r %>% 
  group_by(Ampicillin) %>% 
  summarise(prop_no_gene = sum(n==0)/n(),
            lower=prop.test(sum(n==0),sum(n()))$conf.int[1],
            upper=prop.test(sum(n==0),sum(n()))$conf.int[2])


amp_ng_plot<-ggplot(ng) +
  aes(x=Ampicillin,y=prop_no_gene) +
  geom_point() +
  geom_errorbar(data=ng,aes(ymin=lower,ymax=upper),width=0) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank()) +
  ylab("Proportion with no ARG/mutation") +
  xlab("MIC")


amp_eucast_plot<-ng_plot / amp_gene_dist  +
  plot_annotation(title = "Ampicillin")




############co-amox ecoli################
ec_match2<-select(ec_match,guuid,Coamox_lower,Coamox_upper) %>% filter(!is.na(Coamox_upper))
amrfinder<-read_tsv('amrfinder.tsv')
tem_promotor_mutations<-filter(amrfinder,grepl('TEMp',`Gene symbol`))
ec_match2<-filter(ec_match2,guuid %in% amrfinder$Name)
amrfinder<-amrfinder%>% filter(Class=='BETA-LACTAM' ) %>%  select(guuid=Name,gene=`Gene symbol`)
length(unique(amrfinder$gene))
amrfinder$gene<-str_replace_all(amrfinder$gene,'-','_')
amrfinder_keep<-amrfinder
amrfinder$gene<-ifelse(grepl('blaTEM_',amrfinder$gene) & amrfinder$guuid %in% tem_promotor_mutations$Name,
                       paste0(amrfinder$gene,'_',tem_promotor_mutations$`Gene symbol`[match(amrfinder$guuid,tem_promotor_mutations$Name)]),amrfinder$gene)
amrfinder <- amrfinder %>%
  filter(!grepl("(?<!\\w)blaTEMp_", gene, perl = TRUE))
coamox<-read_tsv('amrfinder.tsv')  %>% filter(Name %in% ec_match2$guuid) %>% filter(grepl('BETA-LACTAM',Class )) %>%  select(guuid=Name,gene=`Gene symbol`) 
coamox$gene<-str_replace_all(coamox$gene,'-','_')

count<-amrfinder %>% filter(str_length(guuid) ==36) %>%  group_by(gene) %>% count() %>% filter(n>=10)


deal_with_rare_genes <- function(gene_name) {
  if (startsWith(gene_name, "ampC")) {
    gsub("^ampC_.*", "ampC_other", gene_name)
  } else if (!grepl("_", gene_name)) {
    # Append '_other' at the end
    paste0(gene_name, "_other")
  } else {
    # Replace after the last '_'
    gsub("_[^_]+$", "_other", gene_name)
  }
}

amrfinder$gene <- ifelse(amrfinder$gene %in% count$gene,
                         amrfinder$gene,
                         sapply(amrfinder$gene, deal_with_rare_genes))

singletons<-amrfinder %>%  group_by(gene) %>% count() %>% filter(n<5)

amrfinder$gene<-ifelse(amrfinder$gene %in% singletons$gene,"other gene",amrfinder$gene)
amrfinder_save1<-amrfinder

common_coamox<-filter(coamox,gene %in% count$gene)
length(unique(coamox$gene)) # 49 amp genes
length(unique(common_coamox$gene)) # 14 occur >=10
x<-data.frame(ec_match2$guuid)
x<-left_join(x,amrfinder,by=c("ec_match2.guuid"="guuid"))
x<-select(x,guuid=ec_match2.guuid,gene)
amrfinder<-data.frame(table(x)) %>% spread(gene,Freq)
names(amrfinder)<-str_replace_all(names(amrfinder),'-','_')
amrfinder[,2:ncol(amrfinder)][amrfinder[,2:ncol(amrfinder)]>1]<-1

genes<-names(amrfinder)
genes<-genes[!genes =='guuid']

amrfinder<-left_join(amrfinder,ec_match2,by=c("guuid"="guuid"))

mlst<-read_tsv('./mlst.tsv',col_names = F) %>% select(guuid=X1,spec=X2,mlst=X3)


big<-c('131','95','69','73')
mlst$mlst<-ifelse(mlst$mlst %in% big, mlst$mlst, 'other')

phylogroups<-read_tsv('./phylogroups.tsv')
mlst<-left_join(mlst,phylogroups,by=c("guuid"))
mlst$mlst<-ifelse(mlst$mlst =='other',mlst$phylogroup,mlst$mlst)
other<-c('U','U/cryptic','cryptic','E','G')
mlst$mlst<-ifelse(mlst$mlst %in% other,'other',mlst$mlst)
#mlst<-select(mlst,-spec)
#mlst<-pivot_wider(mlst,id_cols = 'guuid',names_from = 'mlst',values_fill = 0,values_from = 'present')
mlst<-select(mlst,guuid,mlst)

amrfinder<-left_join(amrfinder,mlst,by=c("guuid"))

out=NULL

names(amrfinder)<-janitor::make_clean_names(names(amrfinder))

genes<-janitor::make_clean_names(genes)



genes<-append(genes,"mlst")
amrfinder$mlst<-factor(amrfinder$mlst,levels=c('B2','131','95','73','69','A','B1','C','D','F','other'))
(Y <- with(amrfinder, Surv(coamox_lower, coamox_upper, type = "interval2")))
x<-amrfinder$guuid
amrfinder<-select(amrfinder,-guuid,-coamox_lower,-coamox_lower)

a<-amrfinder[,grep("^V",colnames(amrfinder))]
for(gene in genes){
  print(which(genes ==gene))
  
  print(gene)
  
  b<-select(amrfinder,gene)
  
  #c<-cbind(b,a)
  c<-b
  t<-as.formula(paste("Y", paste(names(c), collapse = " + "), sep=' ~ '))
  
  
  
  m<-broom::tidy(survreg(t ,data=c,dist="gaussian"),conf.int=T)
  
  m<-if(gene=='mlst'){m[2:13,]}else{m[2,]}    
  out=rbind(out,data.frame(m))
}
include<-out
save_out<-out
amrfinder_save<-amrfinder
include<-filter(include,!gene %in% mlst$mlst)

amrfinder<-amrfinder[names(amrfinder) %in% include$term | names(amrfinder) %in% names(a) ]
t<-as.formula(paste("Y", paste(names(amrfinder), collapse = " + "), sep=' ~ '))
model<-survreg(t, data=amrfinder,dist="gaussian")
out2<-broom::tidy(survreg(t, data=amrfinder,dist="gaussian"),conf.int=T)
out2$which<-"multivariable no population"

amrfinder<-amrfinder_save
include<-out
amrfinder<-amrfinder[names(amrfinder) %in% include$term | names(amrfinder) %in% names(a) | names(amrfinder) =='mlst' ]
t<-as.formula(paste("Y", paste(names(amrfinder), collapse = " + "), sep=' ~ '))
model<-survreg(t, data=amrfinder,dist="gaussian")
out3<-broom::tidy(survreg(t, data=amrfinder,dist="gaussian"),conf.int=T)
out3$which<-"multivariable population"
coamox_t<-t

coamox_r_genes<-filter(out3,2^(estimate) >1) %>% filter(!grepl('mlst',term))
coamox_line<-2^(log2(8) - model$coefficients[1])
coamox_line_genes<-filter(out3,2^estimate > coamox_line)%>% filter(!grepl('mlst',term))
coamox_all_genes<-filter(out3,!grepl('mlst',term))

out$which<-"univariable"
out<-rbind(out,out2,out3)
save_out<-out
count$gene<-janitor::make_clean_names(count$gene)

out<-filter(out,term %in% count$gene| grepl('mlst',term) | term=='other_gene')



coamox_genes<-read_tsv('amrfinder.tsv') %>% filter(grepl('BETA-LACTAM',Class))
coamox_genes<-unique(coamox_genes$`Gene symbol`)
coamox_genes<-janitor::make_clean_names(coamox_genes)

a<-ifelse(out$term %in% coamox_genes,'red','black')

out<-filter(out,!grepl('blaEC',gene))
out$estimate<-as.numeric(out$estimate)
out$conf.high<-as.numeric(out$conf.high)
out$conf.low<-as.numeric(out$conf.low)


out$term <- case_when(out$term == "bla_tem_1_bla_te_mp_c32t" ~ "bla_tem_1 + blaTEMp_c32t",
                      out$term == "bla_tem_1_bla_te_mp_g162t" ~ "bla_tem_1 + blaTEMp_g162t",
                      TRUE ~ out$term)

X<-c("mlstB1"  ,     "mlstA"    ,   "mlstC"  ,      "mlstD" ,       "mlstE",       
     "mlstF",        "mlstG" ,       "mlstother") 
out$term<-ifelse(out$term %in% X, str_replace(out$term,'mlst','phylogroup '),ifelse(grepl('mlst',out$term),str_replace(out$term,'mlst','ST '),out$term))

out<- out %>% 
  add_row(term = "GAP") %>% 
  add_row(term = "GAP2")

out$estimate<-2^out$estimate
out$conf.high<-2^out$conf.high
out$conf.low<-2^out$conf.low

out$which<-ifelse(grepl("GAP",out$term),'univariable',out$which)
out$which<-factor(out$which,levels=c("univariable","multivariable no population","multivariable population"))

amrfinder_save1$gene <- case_when(amrfinder_save1$gene == "blaTEM_1_blaTEMp_C32T" ~ "bla_tem_1 + blaTEMp_c32t",
                                  amrfinder_save1$gene== "blaTEM_1_blaTEMp_G162T" ~ "bla_tem_1 + blaTEMp_g162t",
                                  TRUE ~ amrfinder_save1$gene)

amrfinder<-amrfinder_save1
m<-select(mlst,guuid,gene=mlst) %>% filter(!is.na(gene))
amrfinder<-rbind(amrfinder,m)
count<-amrfinder %>% group_by(gene) %>% distinct(guuid) %>%  count()
count$gene<-ifelse(str_length(count$gene) >3,janitor::make_clean_names(count$gene),count$gene)
count$gene <- case_when(count$gene == '131' ~ 'ST 131',
                        count$gene == '69' ~ 'ST 69',
                        count$gene == '95' ~ 'ST 95',
                        count$gene == '73' ~ 'ST 73',
                        count$gene == 'A' ~ 'phylogroup A',
                        count$gene == 'B1' ~ 'phylogroup B1',
                        count$gene == 'B2' ~ 'phylogroup B2',
                        count$gene == 'C' ~ 'phylogroup C',
                        count$gene == 'D' ~ 'phylogroup D',
                        count$gene == 'F' ~ 'phylogroup F',
                        count$gene == 'other' ~ 'phylogroup other',
                        TRUE ~ count$gene)

count$gene <- case_when(count$gene== "bla_tem_1_bla_te_mp_c32t" ~ "bla_tem_1 + blaTEMp_c32t",
                        count$gene== "bla_tem_1_bla_te_mp_g162t" ~ "bla_tem_1 + blaTEMp_g162t",
                        TRUE ~ count$gene)



out<-left_join(out,count,by=c("term"="gene"))
out$term<-paste0(out$term, ' n=(',out$n,')')


out$term<-str_replace_all(out$term," n=[(]NA[)]",'')
out$term<-factor(out$term, levels=c("amp_c_c_11t n=(15)","amp_c_c_42t n=(10)","amp_c_t_32a n=(34)",
                                    "bla_ctx_m_14 n=(16)","bla_ctx_m_15 n=(171)","bla_ctx_m_27 n=(32)",
                                    "bla_ec n=(2216)","bla_ec_5 n=(667)","bla_oxa_1 n=(163)",
                                    "bla_shv_1 n=(56)","bla_tem n=(34)","bla_tem_1 n=(1033)",
                                    "bla_tem_1 + blaTEMp_c32t n=(90)","bla_tem_1 + blaTEMp_g162t n=(15)","bla_tem_40 n=(13)",              
                                    "other_gene n=(22)","GAP","GAP2","ST 131 n=(646)","ST 95 n=(463)",
                                    "ST 73 n=(683)","ST 69 n=(413)","phylogroup A n=(246)",
                                    "phylogroup B1 n=(267)","phylogroup C n=(101)","phylogroup D n=(285)",
                                    "phylogroup F n=(136)","phylogroup other n=(118)"))



coamox_plot<-ggplot(out) +
  aes(x=term,y=estimate,color=which,group=which) +
  geom_point(position = position_dodge(0.9)) + geom_errorbar(data=out,aes(ymin=conf.low,ymax=conf.high),position=position_dodge(width=0.9),width=0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90),panel.grid.major.x = element_blank()) +
  scale_y_log10() +
  xlab('Gene') + ylab('Fold-change in MIC') + 
  labs(color="Model") +
  ggtitle("Co-amoxiclav") + geom_hline(yintercept =2^(log2(8) - model$coefficients[1]),linetype="dashed",color="dark red") +
  geom_hline(yintercept=1)+
  scale_x_discrete(labels = function(x) ifelse(x == "GAP" | x=="GAP2", "", x)) +
  scale_color_viridis_d()


coamoxiclav<-out
coamoxiclav$drug<-"coamox"

amrfinder<-amrfinder_save
predictors <- save_out$term
predictors<-ifelse(grepl("mlst",predictors),"mlst",predictors)
predictors<-unique(predictors)
predictors<-predictors[!is.na(predictors)  & ! predictors =="(Intercept)" & ! predictors == "Log(scale)"]

out=NULL



for (excluded_gene in predictors) {
  print(excluded_gene)
  other_predictors <- setdiff(predictors, excluded_gene)
  
  for (current_gene in other_predictors) {
    interaction_term <- paste(excluded_gene, current_gene, sep="*")
    
    formula_str <- as.formula(paste("Y", "~", paste(setdiff(other_predictors,current_gene), collapse = "+"), "+", interaction_term))
    summary_model <- broom::tidy(survreg(formula_str, data = amrfinder, dist = "gaussian"))
    out <- rbind(out,data.frame(summary_model))
    
  }
}




replace_mlst_patterns <- function(texts) {
  texts <- gsub("mlst[A-Za-z0-9_]*:", "mlst:", texts)
  texts <- gsub(":[A-Za-z0-9_]*mlst[A-Za-z0-9_]*", ":mlst", texts)
  return(texts)
}



interaction_terms<-filter(out,p.value < 0.01) %>% filter(grepl(':',term))
interaction_terms$term<-replace_mlst_patterns(interaction_terms$term)
interaction_terms$term
interaction_terms<-unique(interaction_terms$term)

model<-survreg(Y ~ amp_c_c_11t + amp_c_c_42t + amp_c_other + amp_c_t_32a + bla_carb_other + 
                 bla_cmy_other + bla_ctx_m_14 + bla_ctx_m_15 + bla_ctx_m_27 + 
                 bla_ctx_m_other + bla_ec + bla_ec_5 + bla_oxa_1 + bla_oxa_other + 
                 bla_shv_1 + bla_shv_other + bla_tem + bla_tem_1 + bla_tem_1_bla_te_mp_c32t + 
                 bla_tem_1_bla_te_mp_g162t + bla_tem_30_bla_te_mp_other + 
                 bla_tem_33_bla_te_mp_other + bla_tem_40 + bla_tem_other + 
                 omp_c_other + other_gene + mlst,data=amrfinder,dist="gaussian")

elimination<-  step(model, direction="backward")

selected_vars <- names(coef(elimination))
selected_vars<-selected_vars[!selected_vars=="(Intercept)"]
selected_vars<-gsub('mlst.*','mlst',selected_vars)
selected_vars<-unique(selected_vars)

interaction_string <- paste(interaction_terms, collapse = " + ")
new_formula <- as.formula(paste("Y ~ ", paste(selected_vars, collapse = " + "), " + ", interaction_string))
model<-survreg(new_formula,data=amrfinder,dist="gaussian")
elimination<-step(model,direction="backward",k=3.8)



coamox_gene<-read_tsv('amrfinder.tsv') %>% filter(Class=='BETA-LACTAM') %>% filter(Scope=='core')


out<-NULL
amrfinder$guuid<-x
for (i in 1:nrow(amrfinder)) {
  print(i)
  model<-survreg( Y[-i,] ~   amp_c_c_11t + amp_c_c_42t + amp_c_other + 
                    amp_c_t_32a + bla_carb_other + bla_cmy_other + bla_ctx_m_14 + 
                    bla_ctx_m_27 + bla_ctx_m_other + bla_oxa_1 + bla_oxa_other + 
                    bla_shv_1 + bla_shv_other + bla_tem + bla_tem_1 + bla_tem_1_bla_te_mp_c32t + 
                    bla_tem_1_bla_te_mp_g162t + bla_tem_30_bla_te_mp_other + 
                    bla_tem_33_bla_te_mp_other + bla_tem_40 + bla_tem_other + 
                    omp_c_other + other_gene + mlst + amp_c_t_32a:bla_tem_30_bla_te_mp_other + 
                    bla_cmy_other:bla_oxa_1 + bla_oxa_1:bla_ctx_m_15 + bla_oxa_1:bla_tem_1 + 
                    bla_oxa_1:bla_tem_1_bla_te_mp_c32t + bla_oxa_other:bla_tem + 
                    bla_tem_1:other_gene + other_gene:mlst + amp_c_t_32a:bla_ec_5, data = amrfinder[-i,], 
                  dist = "gaussian")
  
  prediction <- predict(model, newdata = amrfinder[i, ])
  correct<-ifelse(prediction >= as.numeric(Y[i,1]) & prediction <= as.numeric(Y[i,2]),1,0)
  
  out<-rbind(out,data.frame(correct,prediction,as.numeric(Y[i,1]), as.numeric(Y[i,2]),amrfinder$guuid[i]))
}


out$amrfinder_prediction<-ifelse(out$amrfinder.guuid.i. %in% coamox_gene$Name,'R','S')
out$amrfinder_prediction<-as.factor(out$amrfinder_prediction)
out$binary<-ifelse(out$as.numeric.Y.i..2.. >3,'R','S')
out$binary_prediction<-ifelse(out$prediction>3,'R','S')
out$binary<-as.factor(out$binary)
out$binary_prediction<-as.factor(out$binary_prediction)

extend_margin <- function(out, n){
  out$margin <- NA # Initialize margin column
  for(i in 0:n){
    condition <- out$prediction >= (out$as.numeric.Y.i..1.. - i) & 
      out$prediction <= (out$as.numeric.Y.i..2.. + i)
    out$margin[is.na(out$margin) & condition] <- i
  }
  return(out)
}


out2<-extend_margin(out,10)
out2$margin<-ifelse(out2$prediction < out2$as.numeric.Y.i..1..,out2$margin *-1,out2$margin)

out2$which<- ifelse(out2$binary_prediction==out2$binary & out2$binary=='S', 'Correct S',
                    ifelse(out2$binary_prediction==out2$binary & out2$binary =='R','Correct R',
                           ifelse(out2$binary_prediction=='S' & out2$binary =='R','VM',
                                  ifelse(out2$binary_prediction =='R' & out2$binary =='S','M',NA))))

out$within_1<-ifelse(out$prediction >= (out$as.numeric.Y.i..1.. -1) & out$prediction <= (out$as.numeric.Y.i..2.. +1),1,0)
out$within_2<-ifelse(out$prediction >= (out$as.numeric.Y.i..1.. -2) & out$prediction <= (out$as.numeric.Y.i..2.. +2),1,0)


caret::confusionMatrix(out$binary_prediction,out$binary)
table(out$correct)
prop.test(1467,(1467+1403))
1403+1467


table(out$within_1)
prop.test(2512,(2512+358))
2512+358

cm<-caret::confusionMatrix(as.factor(out$binary_prediction),reference=as.factor(out$binary),positive="R")
tn<-cm$table[2,2]
fn<-cm$table[2,1]
tp<-cm$table[1,1]
fp<-cm$table[1,2]
#concordance
prop.test(x=(tp+tn),n=(tp+tn+fp+fn))
tp+tn
tp+tn+fp+fn
#sens
prop.test(x=(tp),n=(tp+fn))
tp
tp+fn
#spec
prop.test(x=(tn),n=(tn+fp))
tn
tn+fp

#Maj
prop.test(x=(fp),n=(tn+fp))
fp
tn+fp
#Vmj
prop.test(x=(fn),n=(tp+fn))
fn
tp + fn

#PPV
prop.test(x=(tp),n=(tp+fp))
tp
tp+fp

#NPV
prop.test(x=(tn),n=(tn+fn))
tn
tn+fn

colors<-c("#09820B" ,"#03FF25" ,"#FFA703", "#ED3A37")
names(colors)<-c("Correct R","Correct S","M","VM")
coamox_residual<-ggplot(out2) +
  aes(x=margin, fill=which) +
  geom_bar(stat="count") + theme_minimal() +
  #scale_fill_manual(values=colors) + 
  scale_x_continuous(breaks = seq(min(out2$margin),max(out2$margin),1)) + 
  xlab("True - Predicted MIC") +
  labs(fill="") +
  ggtitle("Co-amoxiclav") +
  scale_fill_viridis_d()

out2$Coamox_pred<-
  case_when(out2$prediction <=1 ~ '<=2/2',
            out2$prediction >1 & out2$prediction <=2 ~ '4/2',
            out2$prediction >2 & out2$prediction <=3 ~ '8/2',
            out2$prediction >3 & out2$prediction <=4 ~ '16/2',
            out2$prediction >4 & out2$prediction <=5 ~ '32/2',
            out2$prediction >5 ~ '>32/2')

em<-select(ec_match,guuid,Coamox)
out2<-left_join(out2,em,by=c("amrfinder.guuid.i."="guuid"))

out2$Coamox_pred<-factor(out2$Coamox_pred,levels=c('<=2/2','4/2','8/2','16/2','32/2','>32/2'))
out2$Coamox<-factor(out2$Coamox,levels=c('<=2/2','4/2','8/2','16/2','32/2','>32/2'))
out2<-filter(out2,!Coamox==">8/2")
coamox_sup_dist<-ggplot() +
  geom_histogram(data=out2,aes(x=Coamox_pred,fill=Coamox),stat="count",) + scale_fill_viridis_d() +
  theme_minimal() + xlab("Predicted co-amoxiclav MIC") + labs(fill="Measured co-amoxiclav\nMIC") +
  ggplot() +
  geom_histogram(data=out2,aes(x=Coamox,fill=Coamox_pred),stat="count") + scale_fill_viridis_d() +
  theme_minimal() + xlab("Measured co-amoxiclav MIC") + labs(fill="Predicted co-amoxiclav\nMIC") +
  plot_annotation(tag_levels = "A")

amrfinder<-read_tsv('amrfinder.tsv')
coamox<-select(ec_match,guuid,Coamox) %>% filter(guuid %in% amrfinder$Name)
tem_promotor_mutations<-filter(amrfinder,grepl('TEMp',`Gene symbol`))
ec_match2<-filter(ec_match2,guuid %in% amrfinder$Name)
amrfinder<-amrfinder%>% filter(Class=='BETA-LACTAM' ) %>%  select(guuid=Name,gene=`Gene symbol`)
length(unique(amrfinder$gene))
amrfinder$gene<-str_replace_all(amrfinder$gene,'-','_')
amrfinder_keep<-amrfinder
amrfinder$gene<-ifelse(grepl('blaTEM_',amrfinder$gene) & amrfinder$guuid %in% tem_promotor_mutations$Name,
                       paste0(amrfinder$gene,'_',tem_promotor_mutations$`Gene symbol`[match(amrfinder$guuid,tem_promotor_mutations$Name)]),amrfinder$gene)
amrfinder <- amrfinder %>%
  filter(!grepl("(?<!\\w)blaTEMp_", gene, perl = TRUE))


count<-amrfinder %>% filter(str_length(guuid) ==36) %>%  group_by(gene) %>% count() %>% filter(n>=10)


deal_with_rare_genes <- function(gene_name) {
  if (startsWith(gene_name, "ampC")) {
    gsub("^ampC_.*", "ampC_other", gene_name)
  } else if (!grepl("_", gene_name)) {
    paste0(gene_name, "_other")
  } else {
    gsub("_[^_]+$", "_other", gene_name)
  }
}

amrfinder$gene <- ifelse(amrfinder$gene %in% count$gene,
                         amrfinder$gene,
                         sapply(amrfinder$gene, deal_with_rare_genes))

singletons<-amrfinder %>%  group_by(gene) %>% count() %>% filter(n<5)

amrfinder$gene<-ifelse(amrfinder$gene %in% singletons$gene,"other gene",amrfinder$gene)


coamox_r_genes<-filter(coamox_r_genes,p.value<0.05)
amrfinder$gene<-janitor::make_clean_names(amrfinder$gene,allow_dupes = T)
amrfinder_r<-filter(amrfinder,gene %in% coamox_r_genes$term)
amrfinder_r<- amrfinder_r %>% group_by(guuid) %>% count()

coamox$Coamox<-factor(coamox$Coamox,levels=c('<=2/2','4/2','8/2','16/2','32/2','>32/2'))
coamox<-filter(coamox,!is.na(Coamox))
amrfinder_r<-left_join(coamox,amrfinder_r,by=c("guuid"="guuid"))
amrfinder_r$n<-ifelse(is.na(amrfinder_r$n),0,amrfinder_r$n)
amrfinder_r<-filter(amrfinder_r,!is.na(Coamox))

r<-c('16/2','32/2','>32/2')

amrfinder_r$R<-ifelse(amrfinder_r$Coamox %in% r,'R','S')
amrfinder_r$R<-factor(amrfinder_r$R,levels=c('S','R'))
coamox_gene_dist<-ggplot(amrfinder_r) +
  aes(x=n,fill=R) +
  geom_histogram(bins=length(unique(amrfinder_r$n))) +
  facet_wrap(~Coamox,nrow = 1) + theme_minimal() +
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),strip.text = element_blank()) +
  xlab("Number of ARGs/mutations")+
  scale_fill_viridis_d() +
  labs(fill="")

ng<-amrfinder_r %>% 
  group_by(Coamox) %>% 
  summarise(prop_no_gene = sum(n==0)/n(),
            lower=prop.test(sum(n==0),sum(n()))$conf.int[1],
            upper=prop.test(sum(n==0),sum(n()))$conf.int[2])

            
coamox_ng_plot<-ggplot(ng) +
  aes(x=Coamox,y=prop_no_gene) +
  geom_point() +
  geom_errorbar(data=ng,aes(ymin=lower,ymax=upper),width=0) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank()) +
  ylab("Proportion with no ARG/mutation") +
  xlab("MIC")


coamox_eucast_plot<-ng_plot / coamox_gene_dist  +
  plot_annotation(title = "Co-amoxiclav")






###############piptaz#############


ec_match2<-select(ec_match,guuid,piptaz_lower,piptaz_upper) %>% filter(!is.na(piptaz_upper))
amrfinder<-read_tsv('amrfinder.tsv')
tem_promotor_mutations<-filter(amrfinder,grepl('TEMp',`Gene symbol`))
ec_match2<-filter(ec_match2,guuid %in% amrfinder$Name)
amrfinder<-amrfinder%>% filter(Class=='BETA-LACTAM' ) %>% filter(Scope=='core') %>%  select(guuid=Name,gene=`Gene symbol`) 
length(unique(amrfinder$gene))
amrfinder$gene<-str_replace_all(amrfinder$gene,'-','_')
amrfinder$gene<-ifelse(grepl('blaTEM_',amrfinder$gene) & amrfinder$guuid %in% tem_promotor_mutations$Name,
                       paste0(amrfinder$gene,'_',tem_promotor_mutations$`Gene symbol`[match(amrfinder$guuid,tem_promotor_mutations$Name)]),amrfinder$gene)
amrfinder <- amrfinder %>%
  filter(!grepl("(?<!\\w)blaTEMp_", gene, perl = TRUE))


piptaz<-read_tsv('amrfinder.tsv')  %>% filter(Name %in% ec_match2$guuid) %>% filter(grepl('BETA-LACTAM',Class )) %>%  select(guuid=Name,gene=`Gene symbol`) 
piptaz$gene<-str_replace_all(piptaz$gene,'-','_')

count<-amrfinder %>% filter(str_length(guuid) ==36) %>%  group_by(gene) %>% count() %>% filter(n>=10)


deal_with_rare_genes <- function(gene_name) {
  if (startsWith(gene_name, "ampC")) {
    gsub("^ampC_.*", "ampC_other", gene_name)
  } else if (!grepl("_", gene_name)) {
    paste0(gene_name, "_other")
  } else {
    gsub("_[^_]+$", "_other", gene_name)
  }
}

amrfinder$gene <- ifelse(amrfinder$gene %in% count$gene,
                         amrfinder$gene,
                         sapply(amrfinder$gene, deal_with_rare_genes))

singletons<-amrfinder %>%  group_by(gene) %>% count() %>% filter(n<10)

amrfinder$gene<-ifelse(amrfinder$gene %in% singletons$gene,"other gene",amrfinder$gene)
amrfinder_save1<-amrfinder
table(amrfinder$gene)

common_piptaz<-filter(piptaz,gene %in% count$gene)
length(unique(piptaz$gene)) # 49 amp genes
length(unique(common_piptaz$gene)) # 12 occur >=10
x<-data.frame(ec_match2$guuid)
x<-left_join(x,amrfinder,by=c("ec_match2.guuid"="guuid"))
x<-select(x,guuid=ec_match2.guuid,gene)
amrfinder<-data.frame(table(x)) %>% spread(gene,Freq)
names(amrfinder)<-str_replace_all(names(amrfinder),'-','_')
amrfinder[,2:ncol(amrfinder)][amrfinder[,2:ncol(amrfinder)]>1]<-1

genes<-names(amrfinder)
genes<-genes[!genes =='guuid']

amrfinder<-left_join(amrfinder,ec_match2,by=c("guuid"="guuid"))

mlst<-read_tsv('./mlst.tsv',col_names = F) %>% select(guuid=X1,spec=X2,mlst=X3)

#mlst$present<-1

big<-c('131','95','69','73')
mlst$mlst<-ifelse(mlst$mlst %in% big, mlst$mlst, 'other')

phylogroups<-read_tsv('./phylogroups.tsv')
mlst<-left_join(mlst,phylogroups,by=c("guuid"))
mlst$mlst<-ifelse(mlst$mlst =='other',mlst$phylogroup,mlst$mlst)
other<-c('U','U/cryptic','cryptic','E','G')
mlst$mlst<-ifelse(mlst$mlst %in% other,'other',mlst$mlst)
mlst<-select(mlst,guuid,mlst)

amrfinder<-left_join(amrfinder,mlst,by=c("guuid"))

out=NULL

names(amrfinder)<-janitor::make_clean_names(names(amrfinder))
genes<-janitor::make_clean_names(genes)


genes<-append(genes,"mlst")
amrfinder$mlst<-factor(amrfinder$mlst,levels=c('B2','131','95','73','69','A','B1','C','D','F','other'))
(Y <- with(amrfinder, Surv(piptaz_lower, piptaz_upper, type = "interval2")))
x<-amrfinder$guuid
amrfinder<-select(amrfinder,-guuid,-piptaz_lower,-piptaz_upper)

a<-amrfinder[,grep("^V",colnames(amrfinder))]
for(gene in genes){
  print(which(genes ==gene))
  
  print(gene)
  
  b<-select(amrfinder,gene)
  
  #c<-cbind(b,a)
  c<-b
  t<-as.formula(paste("Y", paste(names(c), collapse = " + "), sep=' ~ '))
  
  
  
  m<-broom::tidy(survreg(t ,data=c,dist="gaussian"),conf.int=T)
  
  m<-if(gene=='mlst'){m[2:11,]}else{m[2,]}    
  out=rbind(out,data.frame(m))
}
include<-out
save_out<-out
amrfinder_save<-amrfinder
big<-c('st131','st95','st73','st69','other')
include<-filter(include,!gene %in% big)

amrfinder<-amrfinder[names(amrfinder) %in% include$term | names(amrfinder) %in% names(a) ]
t<-as.formula(paste("Y", paste(names(amrfinder), collapse = " + "), sep=' ~ '))
model<-survreg(t, data=amrfinder,dist="gaussian")
out2<-broom::tidy(survreg(t, data=amrfinder,dist="gaussian"),conf.int=T)
out2$which<-"multivariable no population"

amrfinder<-amrfinder_save
include<-out
amrfinder<-amrfinder[names(amrfinder) %in% include$term | names(amrfinder) %in% names(a) | names(amrfinder) =='mlst' ]
t<-as.formula(paste("Y", paste(names(amrfinder), collapse = " + "), sep=' ~ '))
model<-survreg(t, data=amrfinder,dist="gaussian")
out3<-broom::tidy(survreg(t, data=amrfinder,dist="gaussian"),conf.int=T)
out3$which<-"multivariable population"

piptaz_r_genes<-filter(out3,round(estimate,1) >0) %>% filter(!grepl('mlst',term))
piptaz_line<-2^(log2(8) - model$coefficients[1])
piptaz_line_genes<-filter(out3,2^estimate > piptaz_line)%>% filter(!grepl('mlst',term))
piptaz_all_genes<-filter(out3,!grepl('mlst',term))

out$which<-"univariable"
out<-rbind(out,out2,out3)
count$gene<-janitor::make_clean_names(count$gene)

out<-filter(out,term %in% count$gene| grepl('mlst',term) | term=='other_gene')



out<-filter(out,!grepl('blaEC',gene))
out$estimate<-as.numeric(out$estimate)
out$conf.high<-as.numeric(out$conf.high)
out$conf.low<-as.numeric(out$conf.low)



out$term <- case_when(out$term == "bla_tem_1_bla_te_mp_c32t" ~ "bla_tem_1 + blaTEMp_c32t",
                      out$term == "bla_tem_1_bla_te_mp_g162t" ~ "bla_tem_1 + blaTEMp_g162t",
                      TRUE ~ out$term)

X<-c("mlstB1"  ,     "mlstA"    ,   "mlstC"  ,      "mlstD" ,       "mlstE",       
     "mlstF",        "mlstG" ,       "mlstother") 
out$term<-ifelse(out$term %in% X, str_replace(out$term,'mlst','phylogroup '),ifelse(grepl('mlst',out$term),str_replace(out$term,'mlst','ST '),out$term))
 
out<- out %>% 
  add_row(term = "GAP") %>% 
  add_row(term = "GAP2")

out$estimate<-2^out$estimate
out$conf.high<-2^out$conf.high
out$conf.low<-2^out$conf.low

out$which<-ifelse(grepl("GAP",out$term),'univariable',out$which)
out$which<-factor(out$which,levels=c("univariable","multivariable no population","multivariable population"))

amrfinder_save1$gene <- case_when(amrfinder_save1$gene == "blaTEM_1_blaTEMp_C32T" ~ "bla_tem_1 + blaTEMp_c32t",
                                  amrfinder_save1$gene== "blaTEM_1_blaTEMp_G162T" ~ "bla_tem_1 + blaTEMp_g162t",
                                  TRUE ~ amrfinder_save1$gene)

amrfinder<-amrfinder_save1
m<-select(mlst,guuid,gene=mlst) %>% filter(!is.na(gene))
amrfinder<-rbind(amrfinder,m)
count<-amrfinder %>% group_by(gene) %>% distinct(guuid) %>%  count()
count$gene<-ifelse(str_length(count$gene) >3,janitor::make_clean_names(count$gene),count$gene)
count$gene <- case_when(count$gene == '131' ~ 'ST 131',
                        count$gene == '69' ~ 'ST 69',
                        count$gene == '95' ~ 'ST 95',
                        count$gene == '73' ~ 'ST 73',
                        count$gene == 'A' ~ 'phylogroup A',
                        count$gene == 'B1' ~ 'phylogroup B1',
                        count$gene == 'B2' ~ 'phylogroup B2',
                        count$gene == 'C' ~ 'phylogroup C',
                        count$gene == 'D' ~ 'phylogroup D',
                        count$gene == 'F' ~ 'phylogroup F',
                        count$gene == 'other' ~ 'phylogroup other',
                        TRUE ~ count$gene)

count$gene <- case_when(count$gene== "bla_tem_1_bla_te_mp_c32t" ~ "bla_tem_1 + blaTEMp_c32t",
                        count$gene== "bla_tem_1_bla_te_mp_g162t" ~ "bla_tem_1 + blaTEMp_g162t",
                        TRUE ~ count$gene)



out<-left_join(out,count,by=c("term"="gene"))
out$term<-paste0(out$term, ' n=(',out$n,')')


out$term<-str_replace_all(out$term," n=[(]NA[)]",'')
out$term<-factor(out$term, levels=c( "amp_c_c_11t n=(15)","amp_c_c_42t n=(10)","amp_c_t_32a n=(34)",
                                     "bla_ctx_m_14 n=(16)","bla_ctx_m_15 n=(171)","bla_ctx_m_27 n=(32)",
                                     "bla_oxa_1 n=(163)","bla_shv_1 n=(56)","bla_tem n=(34)",
                                     "bla_tem_1 n=(1033)","bla_tem_1 + blaTEMp_c32t n=(90)","bla_tem_1 + blaTEMp_g162t n=(15)",
                                     "bla_tem_40 n=(13)","other_gene n=(71)","GAP","GAP2","ST 131 n=(646)",
                                     "ST 95 n=(463)","ST 73 n=(683)","ST 69 n=(413)",
                                     "phylogroup A n=(246)","phylogroup B1 n=(267)","phylogroup C n=(101)",
                                     "phylogroup D n=(285)","phylogroup F n=(136)","phylogroup other n=(118)"))


piptaz_plot<-ggplot(out) +
  aes(x=term,y=estimate,color=which,group=which) +
  geom_point(position = position_dodge(0.9)) + geom_errorbar(data=out,aes(ymin=conf.low,ymax=conf.high),position=position_dodge(width=0.9),width=0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90),panel.grid.major.x = element_blank()) +
  scale_y_log10() +
  xlab('Gene') + ylab('Fold-change in MIC') + 
  labs(color="Model") +
  ggtitle("Piperacillin-Tazobactam") + geom_hline(yintercept =2^(log2(8) - model$coefficients[1]),linetype="dashed",color="dark red") +
  geom_hline(yintercept=1)+
  scale_x_discrete(labels = function(x) ifelse(x == "GAP" | x=="GAP2", "", x)) +
  scale_color_viridis_d()

piptaz<-out
piptaz$drug<-"piptaz"


amrfinder<-amrfinder_save
predictors <- save_out$term
predictors<-ifelse(grepl("mlst",predictors),"mlst",predictors)
predictors<-unique(predictors)
predictors<-predictors[!is.na(predictors)  & ! predictors =="(Intercept)" & ! predictors == "Log(scale)"]

out=NULL



for (excluded_gene in predictors) {
  print(excluded_gene)
  other_predictors <- setdiff(predictors, excluded_gene)
  
  for (current_gene in other_predictors) {
    interaction_term <- paste(excluded_gene, current_gene, sep="*")
    
    formula_str <- as.formula(paste("Y ~", paste(setdiff(other_predictors,current_gene), collapse = "+"), "+", interaction_term))
    summary_model <- broom::tidy(survreg(formula_str, data = amrfinder, dist = "gaussian"))
    out <- rbind(out,data.frame(summary_model))
    
  }
}




replace_mlst_patterns <- function(texts) {
  texts <- gsub("mlst[A-Za-z0-9_]*:", "mlst:", texts)
  texts <- gsub(":[A-Za-z0-9_]*mlst[A-Za-z0-9_]*", ":mlst", texts)
  return(texts)
}



interaction_terms<-filter(out,p.value < 0.01) %>% filter(grepl(':',term))
interaction_terms$term<-replace_mlst_patterns(interaction_terms$term)
interaction_terms<-unique(interaction_terms$term)
interaction_terms

model<-survreg(Y ~ amp_c_c_11t + amp_c_c_42t + amp_c_t_32a + bla_cmy_other + 
                 bla_ctx_m_14 + bla_ctx_m_15 + bla_ctx_m_27 + bla_ctx_m_other + 
                 bla_oxa_1 + bla_shv_1 + bla_tem + bla_tem_1 + bla_tem_1_bla_te_mp_c32t + 
                 bla_tem_1_bla_te_mp_g162t + bla_tem_40 + bla_tem_other + 
                 other_gene + mlst,data=amrfinder,dist="gaussian")

elimination<-  step(model, direction="backward")

selected_vars <- names(coef(elimination))
selected_vars<-selected_vars[!selected_vars=="(Intercept)"]
selected_vars <- gsub("^mlst.*", "mlst", selected_vars)
selected_vars<-unique(selected_vars)

interaction_string <- paste(interaction_terms, collapse = " + ")
new_formula <- as.formula(paste("Y ~ ", paste(selected_vars, collapse = " + "), " + ", interaction_string))
model<-survreg(new_formula,data=amrfinder,dist="gaussian")
elimination<-step(model,direction="backward",k=3.8)

###leave one out cross validation


piptaz_gene<-read_tsv('amrfinder.tsv') %>% filter(Class=='BETA-LACTAM') %>% filter(Scope=='core')


out<-NULL
amrfinder$guuid<-x
for (i in 1:nrow(amrfinder)) {
  print(i)
  model<-survreg( Y[-i,] ~ amp_c_c_42t + bla_ctx_m_15 + bla_ctx_m_other + 
                    bla_oxa_1 + bla_shv_1 + bla_tem_1 + bla_tem_1_bla_te_mp_c32t + 
                    bla_tem_1_bla_te_mp_g162t + bla_tem_other + other_gene + 
                    mlst + bla_ctx_m_15:bla_oxa_1 + bla_ctx_m_15:bla_tem_1 + 
                    bla_ctx_m_15:other_gene + bla_oxa_1:bla_tem_1_bla_te_mp_c32t + 
                    other_gene:mlst, data = amrfinder[-i,], 
                  dist = "gaussian")
  
  prediction <- predict(model, newdata = amrfinder[i, ])
  correct<-ifelse(prediction >= as.numeric(Y[i,1]) & prediction <= as.numeric(Y[i,2]),1,0)
  
  out<-rbind(out,data.frame(correct,prediction,as.numeric(Y[i,1]), as.numeric(Y[i,2]),amrfinder$guuid[i]))
}

out$amrfinder_prediction<-ifelse(out$amrfinder.guuid.i. %in% piptaz_gene$Name,'amrfinder R','amrfinder S')
out$binary<-ifelse(out$as.numeric.Y.i..2.. >3,'R','S')
out$binary_prediction<-ifelse(out$prediction>3,'R','S')
out$binary<-as.factor(out$binary)
out$binary_prediction<-as.factor(out$binary_prediction)


extend_margin <- function(out, n){
  out$margin <- NA # Initialize margin column
  for(i in 0:n){
    condition <- out$prediction >= (out$as.numeric.Y.i..1.. - i) & 
      out$prediction <= (out$as.numeric.Y.i..2.. + i)
    out$margin[is.na(out$margin) & condition] <- i
  }
  return(out)
}


out2<-extend_margin(out,10)
out2$margin<-ifelse(out2$prediction < out2$as.numeric.Y.i..1..,out2$margin *-1,out2$margin)

out2$which<- ifelse(out2$binary_prediction==out2$binary & out2$binary=='S', 'Correct S',
                    ifelse(out2$binary_prediction==out2$binary & out2$binary =='R','Correct R',
                           ifelse(out2$binary_prediction=='S' & out2$binary =='R','VM',
                                  ifelse(out2$binary_prediction =='R' & out2$binary =='S','M',NA))))

out$within_1<-ifelse(out$prediction >= (out$as.numeric.Y.i..1.. -1) & out$prediction <= (out$as.numeric.Y.i..2.. +1),1,0)
out$within_2<-ifelse(out$prediction >= (out$as.numeric.Y.i..1.. -2) & out$prediction <= (out$as.numeric.Y.i..2.. +2),1,0)


caret::confusionMatrix(out$binary_prediction,out$binary)
table(out$correct)
prop.test(2535,(2535+333))
2535+333


table(out$within_1)
prop.test(2698,(2698+170))
2698+170

cm<-caret::confusionMatrix(as.factor(out$binary_prediction),reference=as.factor(out$binary),positive="R")
tn<-cm$table[2,2]
fn<-cm$table[2,1]
tp<-cm$table[1,1]
fp<-cm$table[1,2]
#concordance
prop.test(x=(tp+tn),n=(tp+tn+fp+fn))
tp+tn
tp+tn+fp+fn
#sens
prop.test(x=(tp),n=(tp+fn))
tp
tp+fn
#spec
prop.test(x=(tn),n=(tn+fp))
tn
tn+fp

#Maj
prop.test(x=(fp),n=(tn+fp))
fp
tn+fp
#Vmj
prop.test(x=(fn),n=(tp+fn))
fn
tp + fn

#PPV
prop.test(x=(tp),n=(tp+fp))
tp
tp+fp

#NPV
prop.test(x=(tn),n=(tn+fn))
tn
tn+fn

colors<-c("#09820B" ,"#03FF25" ,"#FFA703", "#ED3A37")
names(colors)<-c("Correct R","Correct S","M","VM")
piptaz_residual<-ggplot(out2) +
  aes(x=margin, fill=which) +
  geom_bar(stat="count") + theme_minimal() +
  #scale_fill_manual(values=colors) + 
  scale_x_continuous(breaks = seq(min(out2$margin,na.rm = T),max(out2$margin,na.rm = T),1)) + 
  xlab("True - Predicted MIC") +
  labs(fill="") +
  ggtitle("Piperacillin-Tazobactam") +
  scale_fill_viridis_d()

out2$piptaz_pred<-case_when(
  out2$prediction <= 2 ~ '<=4/4',
  out2$prediction >2 & out2$prediction <= 3 ~ '8/4',
  out2$prediction >3 & out2$prediction <=4 ~ '16/4',
  out2$prediction >4 & out2$prediction <=5 ~  '32/4',
  out2$prediction >5 & out2$prediction <=6 ~ '64/4',
  out2$prediction >6 ~ '>64/4'
)

em<-select(ec_match,guuid,Piptaz)

out2<-left_join(out2,em,by=c("amrfinder.guuid.i."="guuid"))
out2$piptaz_pred<-factor(out2$piptaz_pred,levels=c('<=4/4','8/4','16/4','32/4','64/4','>64/4'))
out2$Piptaz<-factor(out2$Piptaz,levels=c('<=4/4','8/4','16/4','32/4','64/4','>64/4'))

piptaz_sup_dist<-ggplot() +
  geom_histogram(data=out2,aes(x=piptaz_pred,fill=Piptaz),stat="count",) + scale_fill_viridis_d() +
  theme_minimal() + xlab("Predicted piperacillin-tazobactam MIC") + labs(fill="Measured piperacillin-tazobactam\nMIC") +
  ggplot() +
  geom_histogram(data=out2,aes(x=Piptaz,fill=piptaz_pred),stat="count") + scale_fill_viridis_d() +
  theme_minimal() + xlab("Measured piperacillin-tazobactam MIC") + labs(fill="Predicted piperacillin-tazobactam\nMIC") +
  plot_annotation(tag_levels = "A")


amrfinder<-read_tsv('amrfinder.tsv')
piptaz<-select(ec_match,guuid,Piptaz)  %>% filter(guuid %in% amrfinder$Name)
tem_promotor_mutations<-filter(amrfinder,grepl('TEMp',`Gene symbol`))
ec_match2<-filter(ec_match2,guuid %in% amrfinder$Name)
amrfinder<-amrfinder%>% filter(Class=='BETA-LACTAM' ) %>%  select(guuid=Name,gene=`Gene symbol`)
length(unique(amrfinder$gene))
amrfinder$gene<-str_replace_all(amrfinder$gene,'-','_')
amrfinder_keep<-amrfinder
amrfinder$gene<-ifelse(grepl('blaTEM_',amrfinder$gene) & amrfinder$guuid %in% tem_promotor_mutations$Name,
                       paste0(amrfinder$gene,'_',tem_promotor_mutations$`Gene symbol`[match(amrfinder$guuid,tem_promotor_mutations$Name)]),amrfinder$gene)
amrfinder <- amrfinder %>%
  filter(!grepl("(?<!\\w)blaTEMp_", gene, perl = TRUE))

count<-amrfinder %>% filter(str_length(guuid) ==36) %>%  group_by(gene) %>% count() %>% filter(n>=10)


deal_with_rare_genes <- function(gene_name) {
  if (startsWith(gene_name, "ampC")) {
    gsub("^ampC_.*", "ampC_other", gene_name)
  } else if (!grepl("_", gene_name)) {
    # Append '_other' at the end
    paste0(gene_name, "_other")
  } else {
    # Replace after the last '_'
    gsub("_[^_]+$", "_other", gene_name)
  }
}

amrfinder$gene <- ifelse(amrfinder$gene %in% count$gene,
                         amrfinder$gene,
                         sapply(amrfinder$gene, deal_with_rare_genes))

singletons<-amrfinder %>%  group_by(gene) %>% count() %>% filter(n<5)

amrfinder$gene<-ifelse(amrfinder$gene %in% singletons$gene,"other gene",amrfinder$gene)


piptaz_r_genes<-filter(piptaz_r_genes,p.value<0.05)
amrfinder$gene<-janitor::make_clean_names(amrfinder$gene,allow_dupes = T)
amrfinder_r<-filter(amrfinder,gene %in% coamox_r_genes$term)
amrfinder_r<- amrfinder_r %>% group_by(guuid) %>% count()

piptaz$Piptaz<-factor(piptaz$Piptaz,levels=c('<=4/4','8/4','16/4','32/4','64/4','>64/4'))
pipaz<-filter(piptaz,!is.na(Piptaz))
amrfinder_r<-left_join(piptaz,amrfinder_r,by=c("guuid"="guuid"))
amrfinder_r$n<-ifelse(is.na(amrfinder_r$n),0,amrfinder_r$n)
amrfinder_r<-filter(amrfinder_r,!is.na(Piptaz))

r<-c('16/4','32/4','64/4','>64/4')

amrfinder_r$R<-ifelse(amrfinder_r$Piptaz %in% r,'R','S')
amrfinder_r$R<-factor(amrfinder_r$R,levels=c('S','R'))
piptaz_gene_dist<-ggplot(amrfinder_r) +
  aes(x=n,fill=R) +
  geom_histogram(bins=length(unique(amrfinder_r$n))) +
  facet_wrap(~Piptaz,nrow = 1) + theme_minimal() +
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),strip.text = element_blank()) +
  xlab("Number of ARGs/mutations")+
  scale_fill_viridis_d() +
  labs(fill="")

ng<-amrfinder_r %>% 
  group_by(Piptaz) %>% 
  summarise(prop_no_gene = sum(n==0)/n(),
            lower=prop.test(sum(n==0),sum(n()))$conf.int[1],
            upper=prop.test(sum(n==0),sum(n()))$conf.int[2])


piptaz_ng_plot<-ggplot(ng) +
  aes(x=Piptaz,y=prop_no_gene) +
  geom_point() +
  geom_errorbar(data=ng,aes(ymin=lower,ymax=upper),width=0) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank()) +
  ylab("Proportion with no ARG/mutation") +
  xlab("MIC")


piptaz_eucast_plot<-ng_plot / piptaz_gene_dist  +
  plot_annotation(title = "Piperacillin-Tazobactam")


############cephalexin################
ec_match2<-select(ec_match,guuid, Cephalexin_lower,Cephalexin_upper) %>% filter(!is.na(Cephalexin_upper))
amrfinder<-read_tsv('amrfinder.tsv')
tem_promotor_mutations<-filter(amrfinder,grepl('TEMp',`Gene symbol`))
ec_match2<-filter(ec_match2,guuid %in% amrfinder$Name)
amrfinder<-amrfinder%>% filter(Class=='BETA-LACTAM' ) %>%  select(guuid=Name,gene=`Gene symbol`)
length(unique(amrfinder$gene))
amrfinder$gene<-str_replace_all(amrfinder$gene,'-','_')

amrfinder$gene<-ifelse(grepl('blaTEM_',amrfinder$gene) & amrfinder$guuid %in% tem_promotor_mutations$Name,
                       paste0(amrfinder$gene,'_',tem_promotor_mutations$`Gene symbol`[match(amrfinder$guuid,tem_promotor_mutations$Name)]),amrfinder$gene)
amrfinder <- amrfinder %>%
  filter(!grepl("(?<!\\w)blaTEMp_", gene, perl = TRUE))

count<-amrfinder %>% filter(str_length(guuid) ==36) %>%  group_by(gene) %>% count() %>% filter(n>=10)



deal_with_rare_genes <- function(gene_name) {
  if (startsWith(gene_name, "ampC")) {
    gsub("^ampC_.*", "ampC_other", gene_name)
  } else if (!grepl("_", gene_name)) {
    # Append '_other' at the end
    paste0(gene_name, "_other")
  } else {
    # Replace after the last '_'
    gsub("_[^_]+$", "_other", gene_name)
  }
}

amrfinder$gene <- ifelse(amrfinder$gene %in% count$gene,
                         amrfinder$gene,
                         sapply(amrfinder$gene, deal_with_rare_genes))

singletons<-amrfinder %>%  group_by(gene) %>% count() %>% filter(n<5)

amrfinder$gene<-ifelse(amrfinder$gene %in% singletons$gene,"other gene",amrfinder$gene)



x<-data.frame(ec_match2$guuid)
x<-left_join(x,amrfinder,by=c("ec_match2.guuid"="guuid"))
x<-select(x,guuid=ec_match2.guuid,gene)
amrfinder<-data.frame(table(x)) %>% spread(gene,Freq)
names(amrfinder)<-str_replace_all(names(amrfinder),'-','_')
amrfinder[,2:ncol(amrfinder)][amrfinder[,2:ncol(amrfinder)]>1]<-1

genes<-names(amrfinder)
genes<-genes[!genes =='guuid']
amrfinder<-left_join(amrfinder,ec_match2,by=c("guuid"="guuid"))

mlst<-read_tsv('./mlst.tsv',col_names = F) %>% select(guuid=X1,spec=X2,mlst=X3)


big<-c('131','95','69','73')
mlst$mlst<-ifelse(mlst$mlst %in% big, mlst$mlst, 'other')

phylogroups<-read_tsv('./phylogroups.tsv')
mlst<-left_join(mlst,phylogroups,by=c("guuid"))
mlst$mlst<-ifelse(mlst$mlst =='other',mlst$phylogroup,mlst$mlst)
other<-c('U','U/cryptic','cryptic','E','G')
mlst$mlst<-ifelse(mlst$mlst %in% other,'other',mlst$mlst)


mlst<-select(mlst,guuid,mlst)

amrfinder<-left_join(amrfinder,mlst,by=c("guuid"))

out=NULL

names(amrfinder)<-janitor::make_clean_names(names(amrfinder))
genes<-janitor::make_clean_names(genes)


genes<-append(genes,"mlst")
amrfinder$mlst<-factor(amrfinder$mlst,levels=c('B2','131','95','73','69','A','B1','C','D','F','other'))
(Y <- with(amrfinder, Surv(cephalexin_lower, cephalexin_upper, type = "interval2")))
x<-amrfinder$guuid
amrfinder<-select(amrfinder,-guuid,-cephalexin_lower,-cephalexin_upper)

a<-amrfinder[,grep("^V",colnames(amrfinder))]
for(gene in genes){
  print(which(genes ==gene))
  
  print(gene)
  
  b<-select(amrfinder,gene)
  
  #c<-cbind(b,a)
  c<-b
  t<-as.formula(paste("Y", paste(names(c), collapse = " + "), sep=' ~ '))
  
  
  
  m<-broom::tidy(survreg(t ,data=c,dist="gaussian"),conf.int=T)
  
  m<-if(gene=='mlst'){m[2:13,]}else{m[2,]}    
  out=rbind(out,data.frame(m))
}
include<-out
amrfinder_save<-amrfinder
big<-c('st131','st95','st73','st69','other')
include<-filter(include,!gene %in% big)
amrfinder<-amrfinder[names(amrfinder) %in% include$term | names(amrfinder) %in% names(a) ]
t<-as.formula(paste("Y", paste(names(amrfinder), collapse = " + "), sep=' ~ '))
model<-survreg(t, data=amrfinder,dist="gaussian")
out2<-broom::tidy(survreg(t, data=amrfinder,dist="gaussian"),conf.int=T)
out2$which<-"multivariable no population"

amrfinder<-amrfinder_save
include<-out
amrfinder<-amrfinder[names(amrfinder) %in% include$term | names(amrfinder) %in% names(a) | names(amrfinder) =='mlst' ]
t<-as.formula(paste("Y", paste(names(amrfinder), collapse = " + "), sep=' ~ '))
model<-survreg(t, data=amrfinder,dist="gaussian")
out3<-broom::tidy(survreg(t, data=amrfinder,dist="gaussian"),conf.int=T)
out3$which<-"multivariable population"



out$which<-"univariable"
out<-rbind(out,out2,out3)

count$gene<-janitor::make_clean_names(count$gene)

out<-filter(out,term %in% count$gene| grepl('mlst',term) | term=='other_gene')


out<-filter(out,!grepl('blaEC',gene))
out$estimate<-as.numeric(out$estimate)
out$conf.high<-as.numeric(out$conf.high)
out$conf.low<-as.numeric(out$conf.low)


xx<-c("mlstB1"  ,     "mlstA"    ,   "mlstC"  ,      "mlstD" ,       "mlstE",       
      "mlstF",        "mlstG" ,       "mlstother") 
out$term<-ifelse(out$term %in% xx, str_replace(out$term,'mlst','phylogroup '),ifelse(grepl('mlst',out$term),str_replace(out$term,'mlst','ST '),out$term))
out$term<- case_when(
  out$term == 'bla_tem_1_bla_te_mp_c32t' ~ 'blaTEM-1 + blaTEMp-C32T',
  out$term == 'bla_tem_1_bla_te_mp_g162t'~ 'blaTEM-1 + blaTEMp-G162T',
  TRUE ~ out$term
)

out$estimate<-2^out$estimate
out$conf.high<-2^out$conf.high
out$conf.low<-2^out$conf.low

out$term<-factor(out$term, levels=c("amp_c_c_11t","amp_c_c_42t","amp_c_t_32a",             
                                    "bla_ctx_m_14","bla_ctx_m_15","bla_ctx_m_27",
                                    "bla_ec","bla_ec_5","bla_oxa_1",
                                    "bla_shv_1","bla_tem","bla_tem_1",
                                    "blaTEM-1 + blaTEMp-C32T", "blaTEM-1 + blaTEMp-G162T", "bla_tem_40",             
                                    "other_gene","ST 131","ST 95",
                                    "ST 73","ST 69","phylogroup A",            
                                    "phylogroup B1","phylogroup C","phylogroup D",             
                                    "phylogroup F","phylogroup other" ))

cephalexin_plot<-ggplot(out) +
  aes(x=term,y=estimate,color=which,group=which) +
  geom_point(position = position_dodge(0.9)) + geom_errorbar(data=out,aes(ymin=conf.low,ymax=conf.high),position="dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90)) +
  xlab('Gene') + ylab('Effect on log2 MIC') + 
  labs(color="Model") + geom_hline(yintercept =2^(log2(16) - model$coefficients[1]),linetype="dashed")+
  ggtitle("Cephalexin") + ylim(-1,6)


cephalexin<-filter(out,which=='multivariable population')
cephalexin$which<-"cephalexin"
cephalexin_line<-2^(log2(16) - model$coefficients[1])

############cefuroxime################
ec_match2<-select(ec_match,guuid, cefuroxime_lower,cefuroxime_upper) %>% filter(!is.na(cefuroxime_upper))
amrfinder<-read_tsv('amrfinder.tsv')
tem_promotor_mutations<-filter(amrfinder,grepl('TEMp',`Gene symbol`))
ec_match2<-filter(ec_match2,guuid %in% amrfinder$Name)
amrfinder<-amrfinder%>% filter(Class=='BETA-LACTAM' ) %>%  select(guuid=Name,gene=`Gene symbol`)
length(unique(amrfinder$gene))
amrfinder$gene<-str_replace_all(amrfinder$gene,'-','_')

amrfinder$gene<-ifelse(grepl('blaTEM_',amrfinder$gene) & amrfinder$guuid %in% tem_promotor_mutations$Name,
                       paste0(amrfinder$gene,'_',tem_promotor_mutations$`Gene symbol`[match(amrfinder$guuid,tem_promotor_mutations$Name)]),amrfinder$gene)
amrfinder <- amrfinder %>%
  filter(!grepl("(?<!\\w)blaTEMp_", gene, perl = TRUE))
coamox<-read_tsv('amrfinder.tsv')  %>% filter(Name %in% ec_match2$guuid) %>% filter(grepl('BETA-LACTAM',Class )) %>%  select(guuid=Name,gene=`Gene symbol`) 
coamox$gene<-str_replace_all(coamox$gene,'-','_')


count<-amrfinder %>% filter(str_length(guuid) ==36) %>%  group_by(gene) %>% count() %>% filter(n>=10)

deal_with_rare_genes <- function(gene_name) {
  if (startsWith(gene_name, "ampC")) {
    gsub("^ampC_.*", "ampC_other", gene_name)
  } else if (!grepl("_", gene_name)) {
    paste0(gene_name, "_other")
  } else {
    gsub("_[^_]+$", "_other", gene_name)
  }
}

amrfinder$gene <- ifelse(amrfinder$gene %in% count$gene,
                         amrfinder$gene,
                         sapply(amrfinder$gene, deal_with_rare_genes))

singletons<-amrfinder %>%  group_by(gene) %>% count() %>% filter(n<5)

amrfinder$gene<-ifelse(amrfinder$gene %in% singletons$gene,"other gene",amrfinder$gene)
amrfinder_save1<-amrfinder


x<-data.frame(ec_match2$guuid)
x<-left_join(x,amrfinder,by=c("ec_match2.guuid"="guuid"))
x<-select(x,guuid=ec_match2.guuid,gene)
amrfinder<-data.frame(table(x)) %>% spread(gene,Freq)
names(amrfinder)<-str_replace_all(names(amrfinder),'-','_')
amrfinder[,2:ncol(amrfinder)][amrfinder[,2:ncol(amrfinder)]>1]<-1

genes<-names(amrfinder)
genes<-genes[!genes =='guuid']
amrfinder<-left_join(amrfinder,ec_match2,by=c("guuid"="guuid"))

mlst<-read_tsv('./mlst.tsv',col_names = F) %>% select(guuid=X1,spec=X2,mlst=X3)


big<-c('131','95','69','73')
mlst$mlst<-ifelse(mlst$mlst %in% big, mlst$mlst, 'other')

phylogroups<-read_tsv('./phylogroups.tsv')
mlst<-left_join(mlst,phylogroups,by=c("guuid"))
mlst$mlst<-ifelse(mlst$mlst =='other',mlst$phylogroup,mlst$mlst)
other<-c('U','U/cryptic','cryptic','E','G')
mlst$mlst<-ifelse(mlst$mlst %in% other,'other',mlst$mlst)


mlst<-select(mlst,guuid,mlst)

amrfinder<-left_join(amrfinder,mlst,by=c("guuid"))

out=NULL

names(amrfinder)<-janitor::make_clean_names(names(amrfinder))
genes<-janitor::make_clean_names(genes)


genes<-append(genes,"mlst")
amrfinder$mlst<-factor(amrfinder$mlst,levels=c('B2','131','95','73','69','A','B1','C','D','F','other'))
(Y <- with(amrfinder, Surv(cefuroxime_lower, cefuroxime_upper, type = "interval2")))
x<-amrfinder$guuid
amrfinder<-select(amrfinder,-guuid,-cefuroxime_lower,-cefuroxime_upper)

a<-amrfinder[,grep("^V",colnames(amrfinder))]
for(gene in genes){
  print(which(genes ==gene))
  
  print(gene)
  
  b<-select(amrfinder,gene)
  
  #c<-cbind(b,a)
  c<-b
  t<-as.formula(paste("Y", paste(names(c), collapse = " + "), sep=' ~ '))
  
  
  
  m<-broom::tidy(survreg(t ,data=c,dist="gaussian"),conf.int=T)
  
  m<-if(gene=='mlst'){m[2:13,]}else{m[2,]}    
  out=rbind(out,data.frame(m))
}
include<-out
out_save<-out
amrfinder_save<-amrfinder
big<-c('st131','st95','st73','st69','other')
include<-filter(include,!gene %in% big)
amrfinder<-amrfinder[names(amrfinder) %in% include$term | names(amrfinder) %in% names(a) ]
t<-as.formula(paste("Y", paste(names(amrfinder), collapse = " + "), sep=' ~ '))
model<-survreg(t, data=amrfinder,dist="gaussian")
out2<-broom::tidy(survreg(t, data=amrfinder,dist="gaussian"),conf.int=T)
out2$which<-"multivariable no population"

amrfinder<-amrfinder_save
include<-out
amrfinder<-amrfinder[names(amrfinder) %in% include$term | names(amrfinder) %in% names(a) | names(amrfinder) =='mlst' ]
t<-as.formula(paste("Y", paste(names(amrfinder), collapse = " + "), sep=' ~ '))
model<-survreg(t, data=amrfinder,dist="gaussian")
out3<-broom::tidy(survreg(t, data=amrfinder,dist="gaussian"),conf.int=T)
out3$which<-"multivariable population"

cefuroxime_r_genes<-filter(out3,2^(estimate) >1) %>% filter(!grepl('mlst',term))
cefuroxime_line<-2^(log2(8) - model$coefficients[1])
cefuroxime_line_genes<-filter(out3,2^estimate > cefuroxime_line)
cefuroxime_all_genes<-filter(out3,!grepl('mlst',term))

out$which<-"univariable"
out<-rbind(out,out2,out3)
#out<-left_join(out,out2,by=c("term"="term")) 

#out<-filter(out,!grepl('V',term))
count$gene<-janitor::make_clean_names(count$gene)

out<-filter(out,term %in% count$gene| grepl('mlst',term) | term=='other_gene')



ceph_genes<-read_tsv('amrfinder.tsv') %>% filter(grepl('CEPHALOSPORIN',Subclass))
ceph_genes<-unique(ceph_genes$`Gene symbol`)
ceph_genes<-janitor::make_clean_names(ceph_genes)

a<-ifelse(out$term %in% ceph_genes,'red','black')

out<-filter(out,!grepl('blaEC',gene))
out$estimate<-as.numeric(out$estimate)
out$conf.high<-as.numeric(out$conf.high)
out$conf.low<-as.numeric(out$conf.low)


xx<-c("mlstB1"  ,     "mlstA"    ,   "mlstC"  ,      "mlstD" ,       "mlstE",       
     "mlstF",        "mlstG" ,       "mlstother") 
out$term<-ifelse(out$term %in% xx, str_replace(out$term,'mlst','phylogroup '),ifelse(grepl('mlst',out$term),str_replace(out$term,'mlst','ST '),out$term))

out<- out %>% 
  add_row(term = "GAP") %>% 
  add_row(term = "GAP2")

out$estimate<-2^out$estimate
out$conf.high<-2^out$conf.high
out$conf.low<-2^out$conf.low


out$which<-ifelse(grepl("GAP",out$term),'univariable',out$which)
out$which<-factor(out$which,levels=c("univariable","multivariable no population","multivariable population"))

amrfinder_save1$gene <- case_when(amrfinder_save1$gene == "blaTEM_1_blaTEMp_C32T" ~ "bla_tem_1 + blaTEMp_c32t",
                                  amrfinder_save1$gene== "blaTEM_1_blaTEMp_G162T" ~ "bla_tem_1 + blaTEMp_g162t",
                                  TRUE ~ amrfinder_save1$gene)

amrfinder<-amrfinder_save1
m<-select(mlst,guuid,gene=mlst) %>% filter(!is.na(gene))
amrfinder<-rbind(amrfinder,m)
count<-amrfinder %>% group_by(gene) %>% distinct(guuid) %>%  count()
count$gene<-ifelse(str_length(count$gene) >3,janitor::make_clean_names(count$gene),count$gene)
count$gene <- case_when(count$gene == '131' ~ 'ST 131',
                        count$gene == '69' ~ 'ST 69',
                        count$gene == '95' ~ 'ST 95',
                        count$gene == '73' ~ 'ST 73',
                        count$gene == 'A' ~ 'phylogroup A',
                        count$gene == 'B1' ~ 'phylogroup B1',
                        count$gene == 'B2' ~ 'phylogroup B2',
                        count$gene == 'C' ~ 'phylogroup C',
                        count$gene == 'D' ~ 'phylogroup D',
                        count$gene == 'F' ~ 'phylogroup F',
                        count$gene == 'other' ~ 'phylogroup other',
                        TRUE ~ count$gene)

count$gene <- case_when(count$gene== "bla_tem_1_bla_te_mp_c32t" ~ "bla_tem_1 + blaTEMp_c32t",
                        count$gene== "bla_tem_1_bla_te_mp_g162t" ~ "bla_tem_1 + blaTEMp_g162t",
                        TRUE ~ count$gene)



out<-left_join(out,count,by=c("term"="gene"))
out$term<-paste0(out$term, ' n=(',out$n,')')


out$term<-str_replace_all(out$term," n=[(]NA[)]",'')
out$term<-factor(out$term, levels=c(  "amp_c_c_11t n=(15)","amp_c_c_42t n=(10)","amp_c_t_32a n=(34)","bla_ctx_m_14 n=(16)" ,     
                                      "bla_ctx_m_15 n=(171)","bla_ctx_m_27 n=(32)","bla_ec n=(2216)","bla_ec_5 n=(667)",        
                                      "bla_oxa_1 n=(163)","bla_shv_1 n=(56)","bla_tem n=(34)","bla_tem_1 n=(1033)",
                                      "bla_tem_1_bla_te_mp_c32t", "bla_tem_1_bla_te_mp_g162t" ,"bla_tem_40 n=(13)","other_gene n=(24)","GAP","GAP2",
                                      "ST 131 n=(646)","ST 95 n=(463)","ST 73 n=(683)","ST 69 n=(413)",
                                      "phylogroup A n=(246)","phylogroup B1 n=(267)","phylogroup C n=(101)","phylogroup D n=(285)",  
                                      "phylogroup F n=(136)","phylogroup other n=(118)"))



cef_plot<-ggplot(out) +
  aes(x=term,y=estimate,color=which,group=which) +
  geom_point(position = position_dodge(0.9)) + geom_errorbar(data=out,aes(ymin=conf.low,ymax=conf.high),position=position_dodge(width=0.9),width=0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90),panel.grid.major.x = element_blank()) +
  scale_y_log10() +
  xlab('Gene') + ylab('Fold-change in MIC') + 
  labs(color="Model") +
  ggtitle("Cefuroxime") + geom_hline(yintercept =2^(log2(8) - model$coefficients[1]),linetype="dashed",color="dark red") +
  geom_hline(yintercept=1)+
  scale_x_discrete(labels = function(x) ifelse(x == "GAP" | x=="GAP2", "", x)) +
  scale_color_viridis_d()
cefuroxime<-filter(out,which=='multivariable population')
cefuroxime$drug<-"cefuroxime"
cefuroxime_line<-2^(log2(8) - model$coefficients[1])


amrfinder<-amrfinder_save
predictors <- out_save$term
predictors<-ifelse(grepl("mlst",predictors),"mlst",predictors)
predictors<-unique(predictors)
predictors<-predictors[!is.na(predictors)  & ! predictors =="(Intercept)" & ! predictors == "Log(scale)"]

out=NULL



for (excluded_gene in predictors) {
  print(excluded_gene)
  other_predictors <- setdiff(predictors, excluded_gene)
  
  for (current_gene in other_predictors) {
    interaction_term <- paste(excluded_gene, current_gene, sep="*")
    
    formula_str <- as.formula(paste("Y", "~", paste(setdiff(other_predictors,current_gene), collapse = "+"), "+", interaction_term))
    summary_model <- broom::tidy(survreg(formula_str, data = amrfinder, dist = "gaussian"))
    out <- rbind(out,data.frame(summary_model))
    
  }
}




replace_mlst_patterns <- function(texts) {
  texts <- gsub("mlst[A-Za-z0-9_]*:", "mlst:", texts)
  texts <- gsub(":[A-Za-z0-9_]*mlst[A-Za-z0-9_]*", ":mlst", texts)
  return(texts)
}



interaction_terms<-filter(out,p.value < 0.01) %>% filter(grepl(':',term))
interaction_terms$term<-replace_mlst_patterns(interaction_terms$term)
interaction_terms<-unique(interaction_terms$term)
interaction_terms

model<-survreg(Y ~ amp_c_c_11t + amp_c_c_42t + amp_c_other + amp_c_t_32a + bla_carb_other + 
                 bla_cmy_other + bla_ctx_m_14 + bla_ctx_m_15 + bla_ctx_m_27 + 
                 bla_ctx_m_other + bla_ec + bla_ec_5 + bla_oxa_1 + bla_oxa_other + 
                 bla_shv_1 + bla_shv_other + bla_tem + bla_tem_1 + bla_tem_1_bla_te_mp_c32t + 
                 bla_tem_1_bla_te_mp_g162t + bla_tem_30_bla_te_mp_other + 
                 bla_tem_33_bla_te_mp_other + bla_tem_40 + bla_tem_other + 
                 omp_c_other + other_gene + mlst
                 ,data=amrfinder,dist="gaussian")

elimination<-  step(model, direction="backward")

selected_vars <- names(coef(elimination))
selected_vars<-selected_vars[!selected_vars=="(Intercept)"]
selected_vars<-selected_vars[!grepl('mlst',selected_vars)]
selected_vars<-c(selected_vars,'mlst')

interaction_string <- paste(interaction_terms, collapse = " + ")
new_formula <- as.formula(paste("Y ~ ", paste(selected_vars, collapse = " + "), " + ", interaction_string))
model<-survreg(new_formula,data=amrfinder,dist="gaussian")
elimination<-step(model,direction="backward",k=3.8)

###leave one out cross validation


ceph_gene<-read_tsv('amrfinder.tsv') %>% filter(Subclass=='CEPHALOSPORIN') %>% filter(Scope=='core')

out<-NULL
amrfinder$guuid<-x
for (i in 1:nrow(amrfinder)) {
  print(i)
  model<-survreg( Y[-i,] ~  amp_c_c_11t + amp_c_c_42t + amp_c_other + 
                    amp_c_t_32a + bla_cmy_other + bla_ctx_m_14 + bla_ctx_m_15 + 
                    bla_ctx_m_27 + bla_ctx_m_other + bla_ec_5 + bla_oxa_1 + bla_oxa_other + 
                    bla_shv_other + bla_tem_1_bla_te_mp_c32t + bla_tem_1_bla_te_mp_g162t + 
                    bla_tem_other + other_gene + mlst + bla_cmy_other:bla_oxa_1 + 
                    bla_ctx_m_15:bla_oxa_1 + bla_ec_5:bla_oxa_1 + bla_oxa_1:bla_tem_1_bla_te_mp_c32t, data = amrfinder[-i,], 
                  dist = "gaussian")
  
  prediction <- predict(model, newdata = amrfinder[i, ])
  correct<-ifelse(prediction >= as.numeric(Y[i,1]) & prediction <= as.numeric(Y[i,2]),1,0)
  
  out<-rbind(out,data.frame(correct,prediction,as.numeric(Y[i,1]), as.numeric(Y[i,2]),amrfinder$guuid[i]))
}


out$amrfinder_prediction<-ifelse(out$amrfinder.guuid.i. %in% ceph_gene$Name,'R','S')
out$amrfinder_prediction<-as.factor(out$amrfinder_prediction)
out$binary<-ifelse(out$as.numeric.Y.i..2.. >3,'R','S')
out$binary_prediction<-ifelse(out$prediction>3,'R','S')
out$binary<-as.factor(out$binary)
out$binary_prediction<-as.factor(out$binary_prediction)

extend_margin <- function(out, n){
  out$margin <- NA # Initialize margin column
  for(i in 0:n){
    condition <- out$prediction >= (out$as.numeric.Y.i..1.. - i) & 
      out$prediction <= (out$as.numeric.Y.i..2.. + i)
    out$margin[is.na(out$margin) & condition] <- i
  }
  return(out)
}


out2<-extend_margin(out,10)
out2$margin<-ifelse(out2$prediction < out2$as.numeric.Y.i..1..,out2$margin *-1,out2$margin)

out2$which<- ifelse(out2$binary_prediction==out2$binary & out2$binary=='S', 'Correct S',
                    ifelse(out2$binary_prediction==out2$binary & out2$binary =='R','Correct R',
                           ifelse(out2$binary_prediction=='S' & out2$binary =='R','VM',
                                  ifelse(out2$binary_prediction =='R' & out2$binary =='S','M',NA))))

out$within_1<-ifelse(out$prediction >= (out$as.numeric.Y.i..1.. -1) & out$prediction <= (out$as.numeric.Y.i..2.. +1),1,0)
out$within_2<-ifelse(out$prediction >= (out$as.numeric.Y.i..1.. -2) & out$prediction <= (out$as.numeric.Y.i..2.. +2),1,0)


caret::confusionMatrix(out$binary_prediction,out$binary)
table(out$correct)
prop.test(1413,(992+1413))
1413+992


table(out$within_1)
prop.test(2314,(2314+91))
2314+91

cm<-caret::confusionMatrix(as.factor(out$binary_prediction),reference=as.factor(out$binary),positive="R")
tn<-cm$table[2,2]
fn<-cm$table[2,1]
tp<-cm$table[1,1]
fp<-cm$table[1,2]
#concordance
prop.test(x=(tp+tn),n=(tp+tn+fp+fn))
tp+tn
tp+tn+fp+fn
#sens
prop.test(x=(tp),n=(tp+fn))
tp
tp+fn
#spec
prop.test(x=(tn),n=(tn+fp))
tn
tn+fp

#Maj
prop.test(x=(fp),n=(tn+fp))
fp
tn+fp
#Vmj
prop.test(x=(fn),n=(tp+fn))
fn
tp + fn

#PPV
prop.test(x=(tp),n=(tp+fp))
tp
tp+fp

#NPV
prop.test(x=(tn),n=(tn+fn))
tn
tn+fn

colors<-c("#09820B" ,"#03FF25" ,"#FFA703", "#ED3A37")
names(colors)<-c("Correct R","Correct S","M","VM")
cefuroxime_residual<-ggplot(out2) +
  aes(x=margin, fill=which) +
  geom_bar(stat="count") + theme_minimal() +
  #scale_fill_manual(values=colors) + 
  scale_x_continuous(breaks = seq(min(out2$margin),max(out2$margin),1)) + 
  xlab("True - Predicted MIC") +
  labs(fill="") +
  ggtitle("Cefuroxime") +
  scale_fill_viridis_d()

em<-select(ec_match,guuid,Cefuroxime)
out2<-left_join(out2,em,by=c("amrfinder.guuid.i."="guuid"))

out2$Cefur_pred<- case_when(
  out2$prediction <=1 ~ '<=2',
  out2$prediction >1 & out2$prediction <= 2 ~ '4',
  out2$prediction >2 & out2$prediction <=3 ~ '8',
  out2$prediction >3 ~ '>8'
)

out2$Cefur_pred<-factor(out2$Cefur_pred,levels=c('<=2','4','8','>8'))
out2$Cefuroxime<-factor(out2$Cefuroxime,levels=c('<=2','4','8','>8'))


cefur_sup_dist<-ggplot() +
  geom_histogram(data=out2,aes(x=Cefur_pred,fill=Cefuroxime),stat="count",) + scale_fill_viridis_d() +
  theme_minimal() + xlab("Predicted cefuroxime MIC") + labs(fill="Measured cefuroxime\nMIC") +
  ggplot() +
  geom_histogram(data=out2,aes(x=Cefuroxime,fill=Cefur_pred),stat="count") + scale_fill_viridis_d() +
  theme_minimal() + xlab("Measured cefuroxime MIC") + labs(fill="Predicted cefuroxime\nMIC") +
  plot_annotation(tag_levels = "A")


amrfinder<-read_tsv('amrfinder.tsv')
cefur<-select(ec_match,guuid,Cefuroxime)  %>% filter(guuid %in% amrfinder$Name)
tem_promotor_mutations<-filter(amrfinder,grepl('TEMp',`Gene symbol`))
ec_match2<-filter(ec_match2,guuid %in% amrfinder$Name)
amrfinder<-amrfinder%>% filter(Class=='BETA-LACTAM' ) %>%  select(guuid=Name,gene=`Gene symbol`)
length(unique(amrfinder$gene))
amrfinder$gene<-str_replace_all(amrfinder$gene,'-','_')
amrfinder_keep<-amrfinder
amrfinder$gene<-ifelse(grepl('blaTEM_',amrfinder$gene) & amrfinder$guuid %in% tem_promotor_mutations$Name,
                       paste0(amrfinder$gene,'_',tem_promotor_mutations$`Gene symbol`[match(amrfinder$guuid,tem_promotor_mutations$Name)]),amrfinder$gene)
amrfinder <- amrfinder %>%
  filter(!grepl("(?<!\\w)blaTEMp_", gene, perl = TRUE))

count<-amrfinder %>% filter(str_length(guuid) ==36) %>%  group_by(gene) %>% count() %>% filter(n>=10)


deal_with_rare_genes <- function(gene_name) {
  if (startsWith(gene_name, "ampC")) {
    gsub("^ampC_.*", "ampC_other", gene_name)
  } else if (!grepl("_", gene_name)) {
    paste0(gene_name, "_other")
  } else {
    gsub("_[^_]+$", "_other", gene_name)
  }
}

amrfinder$gene <- ifelse(amrfinder$gene %in% count$gene,
                         amrfinder$gene,
                         sapply(amrfinder$gene, deal_with_rare_genes))

singletons<-amrfinder %>%  group_by(gene) %>% count() %>% filter(n<5)

amrfinder$gene<-ifelse(amrfinder$gene %in% singletons$gene,"other gene",amrfinder$gene)


cefuroxime_r_genes<-filter(cefuroxime_r_genes,p.value<0.05)
amrfinder$gene<-janitor::make_clean_names(amrfinder$gene,allow_dupes = T)
amrfinder_r<-filter(amrfinder,gene %in% cefuroxime_r_genes$term)
amrfinder_r<- amrfinder_r %>% group_by(guuid) %>% count()

cefur$Cefuroxime<-factor(cefur$Cefuroxime,levels=c('<=2','4','8','>8'))
cefur<-filter(cefur,!is.na(Cefuroxime))
amrfinder_r<-left_join(cefur,amrfinder_r,by=c("guuid"="guuid"))
amrfinder_r$n<-ifelse(is.na(amrfinder_r$n),0,amrfinder_r$n)
amrfinder_r<-filter(amrfinder_r,!is.na(Cefuroxime))

r<-c('>8')

amrfinder_r$R<-ifelse(amrfinder_r$Cefuroxime %in% r,'R','S')
amrfinder_r$R<-factor(amrfinder_r$R,levels=c('S','R'))
cefur_gene_dist<-ggplot(amrfinder_r) +
  aes(x=n,fill=R) +
  geom_histogram(bins=length(unique(amrfinder_r$n))) +
  facet_wrap(~Cefuroxime,nrow = 1) + theme_minimal() +
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),strip.text = element_blank()) +
  xlab("Number of ARGs/mutations")+
  scale_fill_viridis_d() +
  labs(fill="")

ng<-amrfinder_r %>% 
  group_by(Cefuroxime) %>% 
  summarise(prop_no_gene = sum(n==0)/n(),
            lower=prop.test(sum(n==0),sum(n()))$conf.int[1],
            upper=prop.test(sum(n==0),sum(n()))$conf.int[2])


cefur_ng_plot<-ggplot(ng) +
  aes(x=Cefuroxime,y=prop_no_gene) +
  geom_point() +
  geom_errorbar(data=ng,aes(ymin=lower,ymax=upper),width=0) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank()) +
  ylab("Proportion with no ARG/mutation") +
  xlab("MIC")


cefur_eucast_plot<-ng_plot / cefur_gene_dist  +
  plot_annotation(title = "Cefuroxime")




############ceftriaxone################
ec_match2<-select(ec_match,guuid,Ceftriaxone_lower,Ceftriaxone_upper) %>% filter(!is.na(Ceftriaxone_upper))
amrfinder<-read_tsv('amrfinder.tsv')
tem_promotor_mutations<-filter(amrfinder,grepl('TEMp',`Gene symbol`))
ec_match2<-filter(ec_match2,guuid %in% amrfinder$Name)
amrfinder<-amrfinder%>% filter(Class=='BETA-LACTAM' ) %>%  select(guuid=Name,gene=`Gene symbol`)
length(unique(amrfinder$gene))
amrfinder$gene<-str_replace_all(amrfinder$gene,'-','_')

amrfinder$gene<-ifelse(grepl('blaTEM_',amrfinder$gene) & amrfinder$guuid %in% tem_promotor_mutations$Name,
                       paste0(amrfinder$gene,'_',tem_promotor_mutations$`Gene symbol`[match(amrfinder$guuid,tem_promotor_mutations$Name)]),amrfinder$gene)
amrfinder <- amrfinder %>%
  filter(!grepl("(?<!\\w)blaTEMp_", gene, perl = TRUE))

count<-amrfinder %>% filter(str_length(guuid) ==36) %>%  group_by(gene) %>% count() %>% filter(n>=10)


ceph<-read_tsv('amrfinder.tsv')  %>% filter(Name %in% ec_match2$guuid) %>% filter(grepl('CEPHALOSPORIN',Subclass )) %>%  select(guuid=Name,gene=`Gene symbol`) 
ceph$gene<-str_replace_all(ceph$gene,'-','_')

deal_with_rare_genes <- function(gene_name) {
  if (startsWith(gene_name, "ampC")) {
    gsub("^ampC_.*", "ampC_other", gene_name)
  } else if (!grepl("_", gene_name)) {
    paste0(gene_name, "_other")
  } else {
    gsub("_[^_]+$", "_other", gene_name)
  }
}


amrfinder$gene <- ifelse(amrfinder$gene %in% count$gene,
                         amrfinder$gene,
                         sapply(amrfinder$gene, deal_with_rare_genes))

singletons<-amrfinder %>%  group_by(gene) %>% count() %>% filter(n<5)

amrfinder$gene<-ifelse(amrfinder$gene %in% singletons$gene,"other gene",amrfinder$gene)
amrfinder_save1<-amrfinder



common_genes<-filter(ceph,gene %in% count$gene)
length(unique(ceph$gene)) # 24 ceph genes
length(unique(common_genes$gene)) # 8 occur >=10
x<-data.frame(ec_match2$guuid)
x<-left_join(x,amrfinder,by=c("ec_match2.guuid"="guuid"))
x<-select(x,guuid=ec_match2.guuid,gene)
amrfinder<-data.frame(table(x)) %>% spread(gene,Freq)
names(amrfinder)<-str_replace_all(names(amrfinder),'-','_')
amrfinder[,2:ncol(amrfinder)][amrfinder[,2:ncol(amrfinder)]>1]<-1

genes<-names(amrfinder)
genes<-genes[!genes =='guuid']
amrfinder<-left_join(amrfinder,ec_match2,by=c("guuid"="guuid"))

mlst<-read_tsv('./mlst.tsv',col_names = F) %>% select(guuid=X1,spec=X2,mlst=X3)

big<-c('131','95','69','73')
mlst$mlst<-ifelse(mlst$mlst %in% big, mlst$mlst, 'other')
phylogroups<-read_tsv('./phylogroups.tsv')
mlst<-left_join(mlst,phylogroups,by=c("guuid"))
mlst$mlst<-ifelse(mlst$mlst =='other',mlst$phylogroup,mlst$mlst)
other<-c('U','U/cryptic','cryptic','E','G')
mlst$mlst<-ifelse(mlst$mlst %in% other,'other',mlst$mlst)

mlst<-select(mlst,guuid,mlst)

amrfinder<-left_join(amrfinder,mlst,by=c("guuid"))

out=NULL

names(amrfinder)<-janitor::make_clean_names(names(amrfinder))
genes<-janitor::make_clean_names(genes)


genes<-append(genes,"mlst")
amrfinder$mlst<-factor(amrfinder$mlst,levels=c('B2','131','95','73','69','A','B1','C','D','F','other'))
(Y <- with(amrfinder, Surv(ceftriaxone_lower, ceftriaxone_upper, type = "interval2")))
guuids<-amrfinder$guuid
amrfinder<-select(amrfinder,-guuid,-ceftriaxone_lower,-ceftriaxone_upper)

a<-amrfinder[,grep("^V",colnames(amrfinder))]
for(gene in genes){
  print(which(genes ==gene))
  
  print(gene)
  
  b<-select(amrfinder,gene)
  
  c<-b
  t<-as.formula(paste("Y", paste(names(c), collapse = " + "), sep=' ~ '))
  
  
  
  m<-broom::tidy(survreg(t ,data=c,dist="gaussian"),conf.int=T)
  
  m<-if(gene=='mlst'){m[2:13,]}else{m[2,]}    
  out=rbind(out,data.frame(m))
}
out_save<-out
include<-out
amrfinder_save<-amrfinder

include<-filter(include,!gene %in% mlst$mlst)

amrfinder<-amrfinder[names(amrfinder) %in% include$term | names(amrfinder) %in% names(a) ]
t<-as.formula(paste("Y", paste(names(amrfinder), collapse = " + "), sep=' ~ '))
model<-survreg(t, data=amrfinder,dist="gaussian")
out2<-broom::tidy(survreg(t, data=amrfinder,dist="gaussian"),conf.int=T)
out2$which<-"multivariable no population"

amrfinder<-amrfinder_save
include<-out
amrfinder<-amrfinder[names(amrfinder) %in% include$term | names(amrfinder) %in% names(a) | names(amrfinder) =='mlst' ]
t<-as.formula(paste("Y", paste(names(amrfinder), collapse = " + "), sep=' ~ '))
model<-survreg(t, data=amrfinder,dist="gaussian")
out3<-broom::tidy(survreg(t, data=amrfinder,dist="gaussian"),conf.int=T)
out3$which<-"multivariable population"

ceftriaxone_r_genes<-filter(out3,2^(estimate) >1) %>% filter(!grepl('mlst',term))
ceftriaxone_line<-2^(log2(2) - model$coefficients[1])
ceftriaxone_line_genes<-filter(out3,2^estimate > ceftriaxone_line)
ceftriaxone_all_genes<-filter(out3,!grepl('mlst',term))

out$which<-"univariable"
out<-rbind(out,out2,out3)
count$gene<-janitor::make_clean_names(count$gene)

out<-filter(out,term %in% count$gene| grepl('mlst',term) | term=='other_gene')


ceph_genes<-read_tsv('amrfinder.tsv') %>% filter(grepl('CEPHALOSPORIN',Subclass))
ceph_genes<-unique(ceph_genes$`Gene symbol`)
ceph_genes<-janitor::make_clean_names(ceph_genes)



out<-filter(out,!grepl('blaEC',gene))
out$estimate<-as.numeric(out$estimate)
out$conf.high<-as.numeric(out$conf.high)
out$conf.low<-as.numeric(out$conf.low)

x<-c("mlstB1"  ,     "mlstA"    ,   "mlstC"  ,      "mlstD" ,       "mlstE",       
     "mlstF",        "mlstG" ,       "mlstother") 
out$term<-ifelse(out$term %in% x, str_replace(out$term,'mlst','phylogroup '),ifelse(grepl('mlst',out$term),str_replace(out$term,'mlst','ST '),out$term))


out<- out %>% 
  add_row(term = "GAP") %>% 
  add_row(term = "GAP2")

out$estimate<-2^out$estimate
out$conf.high<-2^out$conf.high
out$conf.low<-2^out$conf.low


out$which<-ifelse(grepl("GAP",out$term),'univariable',out$which)
out$which<-factor(out$which,levels=c("univariable","multivariable no population","multivariable population"))

amrfinder_save1$gene <- case_when(amrfinder_save1$gene == "blaTEM_1_blaTEMp_C32T" ~ "bla_tem_1 + blaTEMp_c32t",
                                  amrfinder_save1$gene== "blaTEM_1_blaTEMp_G162T" ~ "bla_tem_1 + blaTEMp_g162t",
                                  TRUE ~ amrfinder_save1$gene)

amrfinder<-amrfinder_save1
m<-select(mlst,guuid,gene=mlst) %>% filter(!is.na(gene))
amrfinder<-rbind(amrfinder,m)
count<-amrfinder %>% group_by(gene) %>% distinct(guuid) %>%  count()
count$gene<-ifelse(str_length(count$gene) >3,janitor::make_clean_names(count$gene),count$gene)
count$gene <- case_when(count$gene == '131' ~ 'ST 131',
                        count$gene == '69' ~ 'ST 69',
                        count$gene == '95' ~ 'ST 95',
                        count$gene == '73' ~ 'ST 73',
                        count$gene == 'A' ~ 'phylogroup A',
                        count$gene == 'B1' ~ 'phylogroup B1',
                        count$gene == 'B2' ~ 'phylogroup B2',
                        count$gene == 'C' ~ 'phylogroup C',
                        count$gene == 'D' ~ 'phylogroup D',
                        count$gene == 'F' ~ 'phylogroup F',
                        count$gene == 'other' ~ 'phylogroup other',
                        TRUE ~ count$gene)

count$gene <- case_when(count$gene== "bla_tem_1_bla_te_mp_c32t" ~ "bla_tem_1 + blaTEMp_c32t",
                        count$gene== "bla_tem_1_bla_te_mp_g162t" ~ "bla_tem_1 + blaTEMp_g162t",
                        TRUE ~ count$gene)



out<-left_join(out,count,by=c("term"="gene"))
out$term<-paste0(out$term, ' n=(',out$n,')')


out$term<-str_replace_all(out$term," n=[(]NA[)]",'')
out$term<-factor(out$term, levels=c("amp_c_c_11t n=(15)","amp_c_c_42t n=(10)","amp_c_t_32a n=(34)","bla_ctx_m_14 n=(16)",     
                                    "bla_ctx_m_15 n=(171)","bla_ctx_m_27 n=(32)","bla_ec n=(2216)","bla_ec_5 n=(667)",
                                    "bla_oxa_1 n=(163)","bla_shv_1 n=(56)","bla_tem n=(34)","bla_tem_1 n=(1033)",
                                    "bla_tem_1_bla_te_mp_c32t","bla_tem_1_bla_te_mp_g162t","bla_tem_40 n=(13)","other_gene n=(24)", "GAP","GAP2",        
                                    "ST 131 n=(646)","ST 95 n=(463)","ST 73 n=(683)","ST 69 n=(413)",
                                    "phylogroup A n=(246)","phylogroup B1 n=(267)","phylogroup C n=(101)","phylogroup D n=(285)",    
                                    "phylogroup F n=(136)","phylogroup other n=(118)" ))

ceph_genes<- c("amp_c_c_11t n=(15)","amp_c_c_42t n=(10)","amp_c_t_32a n=(34)","bla_ctx_m_14 n=(16)","bla_ctx_m_15 n=(171)",  
               "bla_ctx_m_27 n=(32)","bla_ec_5 n=(667)","bla_oxa_1 n=(163)")

a<-ifelse(out$term %in% ceph_genes,"red","black")
a[10:70]<-'black' #gap messes things up
ceftriaxone_plot<-ggplot(out) +
  aes(x=term,y=estimate,color=which,group=which) +
  geom_point(position = position_dodge(0.9)) + geom_errorbar(data=out,aes(ymin=conf.low,ymax=conf.high),position=position_dodge(width=0.9),width=0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90,color=a),panel.grid.major.x = element_blank()) +
  xlab('Gene') + ylab('Fold-change in MIC') + 
  labs(color="Model") +
  geom_hline(yintercept =2^(log2(2) - model$coefficients[1]),linetype="dashed",color="dark red") +
  geom_hline(yintercept=1)+
  scale_x_discrete(labels = function(x) ifelse(x == "GAP" | x=="GAP2", "", x)) +
  scale_color_viridis_d() +
  scale_y_log10(breaks=c(1,3,10,30)) 


ceftriaxone<-filter(out,which=='multivariable population')
ceftriaxone$drug<-"ceftriaxone"
ceftriaxone_line<-2^(log2(2) - model$coefficients[1])

amrfinder<-amrfinder_save


amrfinder<-amrfinder_save
predictors <- out_save$term
predictors<-ifelse(grepl("mlst",predictors),"mlst",predictors)
predictors<-unique(predictors)
predictors<-predictors[!is.na(predictors)  & ! predictors =="(Intercept)" & ! predictors == "Log(scale)"]

out=NULL



for (excluded_gene in predictors) {
  print(excluded_gene)
  other_predictors <- setdiff(predictors, excluded_gene)
  
  for (current_gene in other_predictors) {
    interaction_term <- paste(excluded_gene, current_gene, sep="*")
    
    formula_str <- as.formula(paste("Y ~", paste(setdiff(other_predictors,current_gene), collapse = "+"), "+", interaction_term))
    summary_model <- broom::tidy(survreg(formula_str, data = amrfinder, dist = "gaussian"))
    out <- rbind(out,data.frame(summary_model))
    
  }
}




replace_mlst_patterns <- function(texts) {
  texts <- gsub("mlst[A-Za-z0-9_]*:", "mlst:", texts)
  texts <- gsub(":[A-Za-z0-9_]*mlst[A-Za-z0-9_]*", ":mlst", texts)
  return(texts)
}



interaction_terms<-filter(out,p.value < 0.01) %>% filter(grepl(':',term))
interaction_terms$term<-replace_mlst_patterns(interaction_terms$term)
interaction_terms<-unique(interaction_terms$term)
interaction_terms

model<-survreg(Y ~  amp_c_c_11t + amp_c_c_42t + amp_c_other + amp_c_t_32a + bla_carb_other + 
                 bla_cmy_other + bla_ctx_m_14 + bla_ctx_m_15 + bla_ctx_m_27 + 
                 bla_ctx_m_other + bla_ec + bla_ec_5 + bla_oxa_1 + bla_oxa_other + 
                 bla_shv_1 + bla_shv_other + bla_tem + bla_tem_1 + bla_tem_1_bla_te_mp_c32t + 
                 bla_tem_1_bla_te_mp_g162t + bla_tem_30_bla_te_mp_other + 
                 bla_tem_33_bla_te_mp_other + bla_tem_40 + bla_tem_other + 
                 omp_c_other + other_gene + mlst,data=amrfinder,dist="gaussian")

elimination<-  step(model, direction="backward")

selected_vars <- names(coef(elimination))
selected_vars<-selected_vars[!selected_vars=="(Intercept)"]

interaction_string <- paste(interaction_terms, collapse = " + ")
new_formula <- as.formula(paste("Y ~ ", paste(selected_vars, collapse = " + "), " + ", interaction_string))
model<-survreg(new_formula,data=amrfinder,dist="gaussian")
elimination<-step(model,direction="backward",k=3.8)

###leave one out cross validation


ceph_gene<-read_tsv('amrfinder.tsv') %>% filter(Subclass=='CEPHALOSPORIN') %>% filter(Scope=='core')
ceftriaxone_genes<-c('blaCTX-M-101',
                     'blaCTX-M-14',
                     'blaCTX-M-15',
                     'blaCTX-M-189',
                     'blaCTX-M-1',
                     'blaCTX-M-27',
                     'blaCTX-M-32',
                     'blaCTX-M-3',
                     'blaCTX-M-55',
                     'blaCTX-M-9',
                     'blaSHV-102',
                     'blaSHV-12',
                     'blaTEM-101',
                     'blaTEM-102',
                     'blaTEM-106',
                     'blaTEM-52',
                     'blaTEM-63')

ceftriaxone_genes<-read_tsv('amrfinder.tsv') %>% filter(`Gene symbol` %in% ceftriaxone_genes)
out<-NULL
amrfinder$guuid<-guuids
for (i in 1:nrow(amrfinder)) {
  print(i)
  model<-survreg( Y[-i,] ~ amp_c_c_11t + amp_c_c_42t + bla_cmy_other + 
                    bla_ctx_m_14 + bla_ctx_m_15 + bla_ctx_m_27 + bla_ctx_m_other + 
                    bla_ec + bla_oxa_1 + bla_shv_other + bla_tem_1 + bla_tem_1_bla_te_mp_c32t + 
                    bla_tem_other + other_gene + amp_c_c_11t:bla_ec + bla_ctx_m_15:bla_ec_5 + 
                    bla_ctx_m_other:bla_ec + bla_ec:bla_ec_5 + bla_ec:bla_shv_other + 
                    bla_oxa_1:bla_tem_1_bla_te_mp_c32t + bla_tem_1:other_gene, 
                  data = amrfinder[-i,], 
                  dist = "gaussian")
  
  prediction <- predict(model, newdata = amrfinder[i, ])
  correct<-ifelse(prediction >= as.numeric(Y[i,1]) & prediction <= as.numeric(Y[i,2]),1,0)
  
  out<-rbind(out,data.frame(correct,prediction,as.numeric(Y[i,1]), as.numeric(Y[i,2]),amrfinder$guuid[i]))
}

out$amrfinder_prediction<-ifelse(out$amrfinder.guuid.i. %in% ceph_gene$Name,'amrfinder R','amrfinder S')
out$resfinder_prediction<-ifelse(out$amrfinder.guuid.i. %in% ceftriaxone_genes$Name,'resfinder R','refinder S')
#breakpoint >2 = log2(2) = 1
out$binary<-ifelse(out$as.numeric.Y.i..2.. >1,'R','S')
out$binary_prediction<-ifelse(out$prediction>1,'R','S')
out$binary<-as.factor(out$binary)
out$binary_prediction<-as.factor(out$binary_prediction)

car::linearHypothesis(model, "bla_ctx_m_27 = bla_ctx_m_15")



extend_margin <- function(out, n){
  out$margin <- NA
  for(i in 0:n){
    condition <- out$prediction >= (out$as.numeric.Y.i..1.. - i) & 
      out$prediction <= (out$as.numeric.Y.i..2.. + i)
    out$margin[is.na(out$margin) & condition] <- i
  }
  return(out)
}


out2<-extend_margin(out,10)
out2$margin<-ifelse(out2$prediction < out2$as.numeric.Y.i..1..,out2$margin *-1,out2$margin)

out2$which<- ifelse(out2$binary_prediction==out2$binary & out2$binary=='S', 'Correct S',
                    ifelse(out2$binary_prediction==out2$binary & out2$binary =='R','Correct R',
                           ifelse(out2$binary_prediction=='S' & out2$binary =='R','VM',
                                  ifelse(out2$binary_prediction =='R' & out2$binary =='S','M',NA))))

out$within_1<-ifelse(out$prediction >= (out$as.numeric.Y.i..1.. -1) & out$prediction <= (out$as.numeric.Y.i..2.. +1),1,0)
out$within_2<-ifelse(out$prediction >= (out$as.numeric.Y.i..1.. -2) & out$prediction <= (out$as.numeric.Y.i..2.. +2),1,0)


caret::confusionMatrix(out$binary_prediction,out$binary)
table(out$correct)
prop.test(2788,(2788+82))
2788+82


table(out$within_1)
prop.test(2810,(2810+60))
2810+60

cm<-caret::confusionMatrix(as.factor(out$binary_prediction),reference=as.factor(out$binary),positive="R")
tn<-cm$table[2,2]
fn<-cm$table[2,1]
tp<-cm$table[1,1]
fp<-cm$table[1,2]
#concordance
prop.test(x=(tp+tn),n=(tp+tn+fp+fn))
tp+tn
tp+tn+fp+fn
#sens
prop.test(x=(tp),n=(tp+fn))
tp
tp+fn
#spec
prop.test(x=(tn),n=(tn+fp))
tn
tn+fp

#Maj
prop.test(x=(fp),n=(tn+fp))
fp
tn+fp
#Vmj
prop.test(x=(fn),n=(tp+fn))
fn
tp + fn

#PPV
prop.test(x=(tp),n=(tp+fp))
tp
tp+fp

#NPV
prop.test(x=(tn),n=(tn+fn))
tn
tn+fn


table(out$within_1)
prop.test(2807,(2807+63))
2807+63


colors<-c("#09820B" ,"#03FF25" ,"#FFA703", "#ED3A37")
names(colors)<-c("Correct R","Correct S","M","VM")
ceftriaxone_residual<-ggplot(out2) +
  aes(x=margin, fill=which) +
  geom_bar(stat="count") + theme_minimal() +
  #scale_fill_manual(values=colors) + 
  scale_x_continuous(breaks = seq(min(out2$margin),max(out2$margin),1)) + 
  xlab("True - Predicted MIC") +
  labs(fill="") +
  ggtitle("Ceftriaxone") +
  scale_fill_viridis_d()


amrfinder<-read_tsv('amrfinder.tsv')
ceft<-select(ec_match,guuid,Ceftriaxone)  %>% filter(guuid %in% amrfinder$Name)
tem_promotor_mutations<-filter(amrfinder,grepl('TEMp',`Gene symbol`))
ec_match2<-filter(ec_match2,guuid %in% amrfinder$Name)
amrfinder<-amrfinder%>% filter(Class=='BETA-LACTAM' ) %>%  select(guuid=Name,gene=`Gene symbol`)
length(unique(amrfinder$gene))
amrfinder$gene<-str_replace_all(amrfinder$gene,'-','_')
amrfinder_keep<-amrfinder
amrfinder$gene<-ifelse(grepl('blaTEM_',amrfinder$gene) & amrfinder$guuid %in% tem_promotor_mutations$Name,
                       paste0(amrfinder$gene,'_',tem_promotor_mutations$`Gene symbol`[match(amrfinder$guuid,tem_promotor_mutations$Name)]),amrfinder$gene)
amrfinder <- amrfinder %>%
  filter(!grepl("(?<!\\w)blaTEMp_", gene, perl = TRUE))

count<-amrfinder %>% filter(str_length(guuid) ==36) %>%  group_by(gene) %>% count() %>% filter(n>=10)


deal_with_rare_genes <- function(gene_name) {
  if (startsWith(gene_name, "ampC")) {
    gsub("^ampC_.*", "ampC_other", gene_name)
  } else if (!grepl("_", gene_name)) {
    paste0(gene_name, "_other")
  } else {
    gsub("_[^_]+$", "_other", gene_name)
  }
}

amrfinder$gene <- ifelse(amrfinder$gene %in% count$gene,
                         amrfinder$gene,
                         sapply(amrfinder$gene, deal_with_rare_genes))

singletons<-amrfinder %>%  group_by(gene) %>% count() %>% filter(n<5)

amrfinder$gene<-ifelse(amrfinder$gene %in% singletons$gene,"other gene",amrfinder$gene)


ceft_r_genes<-filter(ceftriaxone_r_genes,p.value<0.05)
amrfinder$gene<-janitor::make_clean_names(amrfinder$gene,allow_dupes = T)
amrfinder_r<-filter(amrfinder,gene %in% ceftriaxone_r_genes$term)
amrfinder_r<- amrfinder_r %>% group_by(guuid) %>% count()

ceft$Ceftriaxone<-case_when(
  ceft$Ceftriaxone == '<=0.5' ~ '<=1',
  ceft$Ceftriaxone == '1' ~ '<=1',
  ceft$Ceftriaxone == '<=1' ~ '<=1',
  TRUE ~ ceft$Ceftriaxone
)
ceft$Ceftriaxone<-factor(ceft$Ceftriaxone,levels=c('<=1','2','4','>4'))
ceft<-filter(ceft,!is.na(Ceftriaxone))
amrfinder_r<-left_join(ceft,amrfinder_r,by=c("guuid"="guuid"))
amrfinder_r$n<-ifelse(is.na(amrfinder_r$n),0,amrfinder_r$n)
amrfinder_r<-filter(amrfinder_r,!is.na(Ceftriaxone))

r<-c('4','>4')

amrfinder_r$R<-ifelse(amrfinder_r$Ceftriaxone %in% r,'R','S')
amrfinder_r$R<-factor(amrfinder_r$R,levels=c('S','R'))
ceft_gene_dist<-ggplot(amrfinder_r) +
  aes(x=n,fill=R) +
  geom_histogram(bins=length(unique(amrfinder_r$n))) +
  facet_wrap(~Ceftriaxone,nrow = 1) + theme_minimal() +
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),strip.text = element_blank()) +
  xlab("Number of ARGs/mutations")+
  scale_fill_viridis_d() +
  labs(fill="")

ng<-amrfinder_r %>% 
  group_by(Ceftriaxone) %>% 
  summarise(prop_no_gene = sum(n==0)/n(),
            lower=prop.test(sum(n==0),sum(n()))$conf.int[1],
            upper=prop.test(sum(n==0),sum(n()))$conf.int[2])


ceft_ng_plot<-ggplot(ng) +
  aes(x=Ceftriaxone,y=prop_no_gene) +
  geom_point() +
  geom_errorbar(data=ng,aes(ymin=lower,ymax=upper),width=0) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank()) +
  ylab("Proportion with no ARG/mutation") +
  xlab("MIC")


ceft_eucast_plot<-ng_plot / ceft_gene_dist  +
  plot_annotation(title = "Ceftriaxone")





############cefepime################
ec_match2<-select(ec_match,guuid, cefepime_lower,cefepime_upper) %>% filter(!is.na(cefepime_upper))
amrfinder<-read_tsv('amrfinder.tsv')
tem_promotor_mutations<-filter(amrfinder,grepl('TEMp',`Gene symbol`))
ec_match2<-filter(ec_match2,guuid %in% amrfinder$Name)
amrfinder<-amrfinder%>% filter(Class=='BETA-LACTAM' ) %>%  select(guuid=Name,gene=`Gene symbol`)
length(unique(amrfinder$gene))
amrfinder$gene<-str_replace_all(amrfinder$gene,'-','_')

amrfinder$gene<-ifelse(grepl('blaTEM_',amrfinder$gene) & amrfinder$guuid %in% tem_promotor_mutations$Name,
                       paste0(amrfinder$gene,'_',tem_promotor_mutations$`Gene symbol`[match(amrfinder$guuid,tem_promotor_mutations$Name)]),amrfinder$gene)
amrfinder <- amrfinder %>%
  filter(!grepl("(?<!\\w)blaTEMp_", gene, perl = TRUE))

count<-amrfinder %>% filter(str_length(guuid) ==36) %>%  group_by(gene) %>% count() %>% filter(n>=10)



deal_with_rare_genes <- function(gene_name) {
  if (startsWith(gene_name, "ampC")) {
    gsub("^ampC_.*", "ampC_other", gene_name)
  } else if (!grepl("_", gene_name)) {
    paste0(gene_name, "_other")
  } else {
    gsub("_[^_]+$", "_other", gene_name)
  }
}

amrfinder$gene <- ifelse(amrfinder$gene %in% count$gene,
                         amrfinder$gene,
                         sapply(amrfinder$gene, deal_with_rare_genes))

singletons<-amrfinder %>%  group_by(gene) %>% count() %>% filter(n<5)

amrfinder$gene<-ifelse(amrfinder$gene %in% singletons$gene,"other gene",amrfinder$gene)



x<-data.frame(ec_match2$guuid)
x<-left_join(x,amrfinder,by=c("ec_match2.guuid"="guuid"))
x<-select(x,guuid=ec_match2.guuid,gene)
amrfinder<-data.frame(table(x)) %>% spread(gene,Freq)
names(amrfinder)<-str_replace_all(names(amrfinder),'-','_')
amrfinder[,2:ncol(amrfinder)][amrfinder[,2:ncol(amrfinder)]>1]<-1

genes<-names(amrfinder)
genes<-genes[!genes =='guuid']
amrfinder<-left_join(amrfinder,ec_match2,by=c("guuid"="guuid"))

mlst<-read_tsv('./mlst.tsv',col_names = F) %>% select(guuid=X1,spec=X2,mlst=X3)


big<-c('131','95','69','73')
mlst$mlst<-ifelse(mlst$mlst %in% big, mlst$mlst, 'other')

phylogroups<-read_tsv('./phylogroups.tsv')
mlst<-left_join(mlst,phylogroups,by=c("guuid"))
mlst$mlst<-ifelse(mlst$mlst =='other',mlst$phylogroup,mlst$mlst)
other<-c('U','U/cryptic','cryptic','E','G')
mlst$mlst<-ifelse(mlst$mlst %in% other,'other',mlst$mlst)

mlst<-select(mlst,guuid,mlst)

amrfinder<-left_join(amrfinder,mlst,by=c("guuid"))

out=NULL

names(amrfinder)<-janitor::make_clean_names(names(amrfinder))
genes<-janitor::make_clean_names(genes)


genes<-append(genes,"mlst")
amrfinder$mlst<-factor(amrfinder$mlst,levels=c('B2','131','95','73','69','A','B1','C','D','F','other'))
(Y <- with(amrfinder, Surv(cefepime_lower, cefepime_upper, type = "interval2")))
guuids<-amrfinder$guuid
amrfinder<-select(amrfinder,-guuid,-cefepime_lower,-cefepime_upper)

a<-amrfinder[,grep("^V",colnames(amrfinder))]
for(gene in genes){
  print(which(genes ==gene))
  
  print(gene)
  
  b<-select(amrfinder,gene)
  
  c<-b
  t<-as.formula(paste("Y", paste(names(c), collapse = " + "), sep=' ~ '))
  
  
  
  m<-broom::tidy(survreg(t ,data=c,dist="gaussian"),conf.int=T)
  
  m<-if(gene=='mlst'){m[2:13,]}else{m[2,]}    
  out=rbind(out,data.frame(m))
}
include<-out
out_save<-out
amrfinder_save<-amrfinder
big<-c('st131','st95','st73','st69','other')
include<-filter(include,!gene %in% big)
amrfinder<-amrfinder[names(amrfinder) %in% include$term | names(amrfinder) %in% names(a) ]
t<-as.formula(paste("Y", paste(names(amrfinder), collapse = " + "), sep=' ~ '))
model<-survreg(t, data=amrfinder,dist="gaussian")
out2<-broom::tidy(survreg(t, data=amrfinder,dist="gaussian"),conf.int=T)
out2$which<-"multivariable no population"

amrfinder<-amrfinder_save
include<-out
amrfinder<-amrfinder[names(amrfinder) %in% include$term | names(amrfinder) %in% names(a) | names(amrfinder) =='mlst' ]
t<-as.formula(paste("Y", paste(names(amrfinder), collapse = " + "), sep=' ~ '))
model<-survreg(t, data=amrfinder,dist="gaussian")
out3<-broom::tidy(survreg(t, data=amrfinder,dist="gaussian"),conf.int=T)
out3$which<-"multivariable population"



out$which<-"univariable"
out<-rbind(out,out2,out3)
count$gene<-janitor::make_clean_names(count$gene)

out<-filter(out,term %in% count$gene| grepl('mlst',term) | term=='other_gene')




ceph_genes<-read_tsv('amrfinder.tsv') %>% filter(grepl('CEPHALOSPORIN',Subclass))
ceph_genes<-unique(ceph_genes$`Gene symbol`)
ceph_genes<-janitor::make_clean_names(ceph_genes)

a<-ifelse(out$term %in% ceph_genes,'red','black')

out<-filter(out,!grepl('blaEC',gene))
out$estimate<-as.numeric(out$estimate)
out$conf.high<-as.numeric(out$conf.high)
out$conf.low<-as.numeric(out$conf.low)
out$estimate<-2^out$estimate
out$conf.high<-2^out$conf.high
out$conf.low<-2^out$conf.low

xx<-c("mlstB1"  ,     "mlstA"    ,   "mlstC"  ,      "mlstD" ,       "mlstE",       
      "mlstF",        "mlstG" ,       "mlstother") 
out$term<-ifelse(out$term %in% xx, str_replace(out$term,'mlst','phylogroup '),ifelse(grepl('mlst',out$term),str_replace(out$term,'mlst','ST '),out$term))

out$term <- case_when(out$term== "bla_tem_1_bla_te_mp_c32t" ~ "blaTEM-1 + blaTEMp-C32T",
                        out$term== "bla_tem_1_bla_te_mp_g162t" ~ "blaTEM-1 + blaTEMp-G162T",
                        TRUE ~ out$term)



out$term<-factor(out$term, levels=c("amp_c_c_11t","amp_c_c_42t","amp_c_t_32a",
                                    "bla_ctx_m_14","bla_ctx_m_15","bla_ctx_m_27",
                                    "bla_ec","bla_ec_5","bla_oxa_1",
                                    "bla_shv_1","bla_tem","bla_tem_1",
                                    "blaTEM-1 + blaTEMp-C32T","blaTEM-1 + blaTEMp-G162T" ,"bla_tem_40",              
                                    "other_gene","ST 131","ST 95",
                                    "ST 73","ST 69","phylogroup A",
                                    "phylogroup B1","phylogroup C","phylogroup D",
                                    "phylogroup F","phylogroup other"))

cefepime_plot<-ggplot(out) +
  aes(x=term,y=estimate,color=which,group=which) +
  geom_point(position = position_dodge(0.9)) + geom_errorbar(data=out,aes(ymin=conf.low,ymax=conf.high),position="dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90,color=a)) +
  xlab('Gene') + ylab('Effect on log2 MIC') + 
  labs(color="Model") + geom_hline(yintercept =2^(log2(4) - model$coefficients[1]),linetype="dashed")+
  ggtitle("Cefepime") + ylim(-2,6)

cefepime_line=2^(log2(4)-model$coefficients[1])
cefepime<-filter(out,which=='multivariable population')
cefepime$which<-"cefepime"

amrfinder<-amrfinder_save
predictors <- out_save$term
predictors<-ifelse(grepl("mlst",predictors),"mlst",predictors)
predictors<-unique(predictors)
predictors<-predictors[!is.na(predictors)  & ! predictors =="(Intercept)" & ! predictors == "Log(scale)"]

out=NULL



for (excluded_gene in predictors) {
  print(excluded_gene)
  other_predictors <- setdiff(predictors, excluded_gene)
  
  for (current_gene in other_predictors) {
    interaction_term <- paste(excluded_gene, current_gene, sep="*")
    
    formula_str <- as.formula(paste("Y ~", paste(setdiff(other_predictors,current_gene), collapse = "+"), "+", interaction_term))
    summary_model <- broom::tidy(survreg(formula_str, data = amrfinder, dist = "gaussian"))
    out <- rbind(out,data.frame(summary_model))
    
  }
}




replace_mlst_patterns <- function(texts) {
  texts <- gsub("mlst[A-Za-z0-9_]*:", "mlst:", texts)
  texts <- gsub(":[A-Za-z0-9_]*mlst[A-Za-z0-9_]*", ":mlst", texts)
  return(texts)
}



interaction_terms<-filter(out,p.value < 0.01) %>% filter(grepl(':',term))
interaction_terms$term<-replace_mlst_patterns(interaction_terms$term)
interaction_terms<-unique(interaction_terms$term)
interaction_terms


model<-survreg(Y ~  amp_c_c_11t + amp_c_c_42t + amp_c_other + amp_c_t_32a + bla_carb_other + 
                 bla_cmy_other + bla_ctx_m_14 + bla_ctx_m_15 + bla_ctx_m_27 + 
                 bla_ctx_m_other + bla_ec + bla_ec_5 + bla_oxa_1 + bla_oxa_other + 
                 bla_shv_1 + bla_shv_other + bla_tem + bla_tem_1 + bla_tem_1_bla_te_mp_c32t + 
                 bla_tem_1_bla_te_mp_g162t + bla_tem_30_bla_te_mp_other + 
                 bla_tem_33_bla_te_mp_other + bla_tem_40 + bla_tem_other + 
                 omp_c_other + other_gene + mlst,data=amrfinder,dist="gaussian")

elimination<-  step(model, direction="backward")

selected_vars <- names(coef(elimination))
selected_vars<-selected_vars[!selected_vars=="(Intercept)"]

interaction_string <- paste(interaction_terms, collapse = " + ")
new_formula <- as.formula(paste("Y ~ ", paste(selected_vars, collapse = " + "), " + ", interaction_string))
model<-survreg(new_formula,data=amrfinder,dist="gaussian")
elimination<-step(model,direction="backward",k=3.8)

###leave one out cross validation


ceph_gene<-read_tsv('amrfinder.tsv') %>% filter(Subclass=='CEPHALOSPORIN') %>% filter(Scope=='core')


out<-NULL
amrfinder$guuid<-guuids
for (i in 1:nrow(amrfinder)) {
  print(i)
  model<-survreg( Y[-i,] ~ amp_c_c_11t + bla_ctx_m_14 + bla_ctx_m_15 + 
                    bla_ctx_m_27 + bla_ctx_m_other + bla_oxa_1 + bla_shv_other + 
                    bla_tem_33_bla_te_mp_other + bla_tem_other + other_gene + 
                    amp_c_c_11t:bla_ec + amp_c_c_11t:bla_ec_5 + bla_ec:amp_c_t_32a + 
                    bla_ec_5:amp_c_t_32a + amp_c_t_32a:mlst + bla_ctx_m_14:mlst + 
                    bla_ctx_m_15:bla_oxa_1 + bla_ctx_m_15:bla_tem_1 + bla_ctx_m_15:mlst + 
                    bla_ctx_m_27:bla_ec + bla_ctx_m_27:bla_ec_5 + bla_ctx_m_27:bla_tem_1 + 
                    bla_ctx_m_27:mlst + bla_ctx_m_other:bla_ec + bla_ctx_m_other:mlst + 
                    bla_ec:bla_ec_5 + bla_ec:mlst + bla_oxa_1:bla_ec_5 + bla_ec_5:mlst + 
                    bla_oxa_1:mlst + mlst:bla_tem_1_bla_te_mp_c32t + bla_tem_other:mlst, data = amrfinder[-i,], 
                  dist = "gaussian")
  
  prediction <- predict(model, newdata = amrfinder[i, ])
  correct<-ifelse(prediction >= as.numeric(Y[i,1]) & prediction <= as.numeric(Y[i,2]),1,0)
  
  out<-rbind(out,data.frame(correct,prediction,as.numeric(Y[i,1]), as.numeric(Y[i,2]),amrfinder$guuid[i]))
}

out$amrfinder_prediction<-ifelse(out$amrfinder.guuid.i. %in% ceph_gene$Name,'amrfinder R','amrfinder S')
out$resfinder_prediction<-ifelse(out$amrfinder.guuid.i. %in% ceftriaxone_genes$Name,'resfinder R','refinder S')
#breakpoint >4 = log2(2) = 2
out$binary<-ifelse(out$as.numeric.Y.i..2.. >2,'R','S')
out$binary_prediction<-ifelse(out$prediction>2,'R','S')
out$binary<-as.factor(out$binary)
out$binary_prediction<-as.factor(out$binary_prediction)

car::linearHypothesis(model, "bla_ctx_m_27 = bla_ctx_m_15")



extend_margin <- function(out, n){
  out$margin <- NA # Initialize margin column
  for(i in 0:n){
    condition <- out$prediction >= (out$as.numeric.Y.i..1.. - i) & 
      out$prediction <= (out$as.numeric.Y.i..2.. + i)
    out$margin[is.na(out$margin) & condition] <- i
  }
  return(out)
}


out2<-extend_margin(out,10)
out2$margin<-ifelse(out2$prediction < out2$as.numeric.Y.i..1..,out2$margin *-1,out2$margin)

out2$which<- ifelse(out2$binary_prediction==out2$binary & out2$binary=='S', 'Correct S',
                    ifelse(out2$binary_prediction==out2$binary & out2$binary =='R','Correct R',
                           ifelse(out2$binary_prediction=='S' & out2$binary =='R','VM',
                                  ifelse(out2$binary_prediction =='R' & out2$binary =='S','M',NA))))

out$within_1<-ifelse(out$prediction >= (out$as.numeric.Y.i..1.. -1) & out$prediction <= (out$as.numeric.Y.i..2.. +1),1,0)
out$within_2<-ifelse(out$prediction >= (out$as.numeric.Y.i..1.. -2) & out$prediction <= (out$as.numeric.Y.i..2.. +2),1,0)


caret::confusionMatrix(out$binary_prediction,out$binary)
table(out$correct)
prop.test(2188,(2188+216))
2188+216


table(out$within_1)
prop.test(2298,(2298+106))
2298+106

cm<-caret::confusionMatrix(as.factor(out$binary_prediction),reference=as.factor(out$binary),positive="R")
tn<-cm$table[2,2]
fn<-cm$table[2,1]
tp<-cm$table[1,1]
fp<-cm$table[1,2]
#concordance
prop.test(x=(tp+tn),n=(tp+tn+fp+fn))
tp+tn
tp+tn+fp+fn
#sens
prop.test(x=(tp),n=(tp+fn))
tp
tp+fn
#spec
prop.test(x=(tn),n=(tn+fp))
tn
tn+fp

#Maj
prop.test(x=(fp),n=(tn+fp))
fp
tn+fp
#Vmj
prop.test(x=(fn),n=(tp+fn))
fn
tp + fn

#PPV
prop.test(x=(tp),n=(tp+fp))
tp
tp+fp

#NPV
prop.test(x=(tn),n=(tn+fn))
tn
tn+fn





colors<-c("#09820B" ,"#03FF25" ,"#FFA703", "#ED3A37")
names(colors)<-c("Correct R","Correct S","M","VM")
cefepime_residual<-ggplot(out2) +
  aes(x=margin, fill=which) +
  geom_bar(stat="count") + theme_minimal() +
  #scale_fill_manual(values=colors) + 
  scale_x_continuous(breaks = seq(min(out2$margin),max(out2$margin),1)) + 
  xlab("True - Predicted MIC") +
  labs(fill="") +
  ggtitle("Cefepime") +
  scale_fill_viridis_d()



######ceph combined############

cefuroxime<-select(cefuroxime,-n,-drug)
ceftriaxone<-select(ceftriaxone,-n,-drug)
cefuroxime$which<-"cefuroxime"
ceftriaxone$which<-"ceftriaxone"
ceph_all<-rbind(cephalexin,cefuroxime,ceftriaxone,cefepime)
ceph_all$term<-str_replace_all(ceph_all$term,' N[=].*','')
ceph_all$term<-str_replace_all(ceph_all$term,' n[=][(].*','')
ceph_all<-filter(ceph_all,!grepl('ST',term) & !grepl('phylogroup',term))

ceph_all$term <- case_when(ceph_all$term== "bla_tem_1_bla_te_mp_c32t" ~ "blaTEM-1 + blaTEMp-C32T",
                      ceph_all$term== "bla_tem_1_bla_te_mp_g162t" ~ "blaTEM-1 + blaTEMp-G162T",
                      TRUE ~ ceph_all$term)



palette=c( "#CF2F2C" ,"#449419", "#1A36B0" ,"#4F0A09")
ceph_all$line<-ifelse(ceph_all$which=='cephalexin',cephalexin_line,
                      ifelse(ceph_all$which == 'cefuroxime', cefuroxime_line,
                             ifelse(ceph_all$which == 'ceftriaxone', ceftriaxone_line,
                                    cefepime_line)))

ceph_all$which<-factor(ceph_all$which,levels=c("cephalexin","cefuroxime","ceftriaxone","cefepime"))

format_names_nicely <- function(label) {
  if (startsWith(label, "bla")) {
    parts <- strsplit(label, "bla", fixed = TRUE)[[1]]
    
    new_label <- toupper(gsub("_", "-", sub("_", "", parts[2]), fixed = TRUE))
    return(paste("bla", new_label, sep=""))
  } 
  else if (startsWith(label, "amp")) {
    parts <- strsplit(label, "amp_", fixed = TRUE)[[1]]
    new_label <- paste0(toupper(substring(parts[2], 1, 1)), substring(parts[2], 2))
    new_label <- gsub("_", "-", new_label, fixed = TRUE)
    return(paste("amp", new_label, sep=""))
  } else {
    return(gsub("-", " ", label, fixed = TRUE))
  }
}


ceph_all$term<-as.character(ceph_all$term)
ceph_all$term <- sapply(ceph_all$term, format_names_nicely)

ggplot(ceph_all) +
  aes(x=term,y=estimate,color=which) +
  geom_point(position = position_dodge(0.9)) + 
  geom_errorbar(aes(ymin=conf.low,ymax=conf.high),position="dodge",width=0) +
  scale_color_viridis_d() +
  geom_hline(aes(yintercept = line),linetype="dashed",color="dark red") +
  theme_minimal() + theme(axis.text.x = element_text(angle=90),strip.text = element_blank(),panel.grid.major.x = element_blank()) +
  labs(color="Antibiotic") + ylab("Fold-change in MIC") + xlab("Gene/mutation") +
  facet_wrap(~which)  + 
  scale_y_log10()
  

library(patchwork)

cephalexin_plot + cefuroxime_plot + ceftriaxone_plot + cefepime_plot  +
  plot_layout(guides="collect")
###################cipro#####################

ec_match2<-select(ec_match,guuid,Ciprofloxacin_lower,Ciprofloxacin_upper) %>% filter(!is.na(Ciprofloxacin_upper))
amrfinder<-read_tsv('amrfinder.tsv') 
ec_match2<-filter(ec_match2,guuid %in% amrfinder$Name)
amrfinder<-amrfinder %>% filter(grepl('QUINOLONE',Class )) %>%  select(guuid=Name,gene=`Gene symbol`)


a_count<-amrfinder %>% group_by(guuid) %>% count()
count<-amrfinder %>% filter(str_length(guuid) ==36) %>%  group_by(gene) %>% count() %>% filter(n>=10)

deal_with_rare_genes <- function(gene_name) {
  if (startsWith(gene_name, "ampC")) {
    gsub("^ampC_.*", "ampC_other", gene_name)
  } else if (!grepl("_", gene_name)) {
    paste0(gene_name, "_other")
  } else {
    gsub("_[^_]+$", "_other", gene_name)
  }
}

amrfinder$gene <- ifelse(amrfinder$gene %in% count$gene,
                         amrfinder$gene,
                         sapply(amrfinder$gene, deal_with_rare_genes))

singletons<-amrfinder %>%  group_by(gene) %>% count() %>% filter(n<5)

amrfinder$gene<-ifelse(amrfinder$gene %in% singletons$gene,"other gene",amrfinder$gene)
amrfinder_save1<-amrfinder



x<-data.frame(ec_match2$guuid)
x<-left_join(x,amrfinder,by=c("ec_match2.guuid"="guuid"))
x<-select(x,guuid=ec_match2.guuid,gene)
amrfinder<-data.frame(table(x)) %>% spread(gene,Freq)
names(amrfinder)<-str_replace_all(names(amrfinder),'-','_')
amrfinder[,2:ncol(amrfinder)][amrfinder[,2:ncol(amrfinder)]>1]<-1
genes<-names(amrfinder)
genes<-genes[!genes =='guuid']


amrfinder<-left_join(amrfinder,ec_match2,by=c("guuid"="guuid"))

mlst<-read_tsv('./mlst.tsv',col_names = F) %>% select(guuid=X1,spec=X2,mlst=X3)


big<-c('131','95','69','73')
mlst$mlst<-ifelse(mlst$mlst %in% big, mlst$mlst, 'other')

phylogroups<-read_tsv('./phylogroups.tsv')
mlst<-left_join(mlst,phylogroups,by=c("guuid"))
mlst$mlst<-ifelse(mlst$mlst =='other',mlst$phylogroup,mlst$mlst)
other<-c('U','U/cryptic','cryptic','E','G')
mlst$mlst<-ifelse(mlst$mlst %in% other,'other',mlst$mlst)

mlst<-select(mlst,guuid,mlst)

amrfinder<-left_join(amrfinder,mlst,by=c("guuid"))


out=NULL

names(amrfinder)<-janitor::make_clean_names(names(amrfinder))
genes<-janitor::make_clean_names(genes)


genes<-append(genes,"mlst")
amrfinder$mlst<-factor(amrfinder$mlst,levels=c('B2','131','95','73','69','A','B1','C','D','F','other'))
(Y <- with(amrfinder, Surv(ciprofloxacin_lower, ciprofloxacin_upper, type = "interval2")))
x<-amrfinder$guuid

amrfinder<-select(amrfinder,-guuid,-ciprofloxacin_lower,-ciprofloxacin_upper)

a<-amrfinder[,grep("^V",colnames(amrfinder))]
for(gene in genes){
  print(which(genes ==gene))
  
  print(gene)
  
  b<-select(amrfinder,gene)
  
  #c<-cbind(b,a)
  c<-b
  t<-as.formula(paste("Y", paste(names(c), collapse = " + "), sep=' ~ '))
  
  
  
  m<-broom::tidy(survreg(t ,data=c,dist="gaussian"),conf.int=T)
  
  m<-if(gene=='mlst'){m[2:11,]}else{m[2,]}    
  out=rbind(out,data.frame(m))
}
include<-out
out_save<-out

amrfinder_save<-amrfinder
big<-c('st131','st95','st73','st69','other')
include<-filter(include,!gene %in% big)
amrfinder<-amrfinder[names(amrfinder) %in% include$term | names(amrfinder) %in% names(a) ]
t<-as.formula(paste("Y", paste(names(amrfinder), collapse = " + "), sep=' ~ '))
model<-survreg(t, data=amrfinder,dist="gaussian")
out2<-broom::tidy(survreg(t, data=amrfinder,dist="gaussian"),conf.int=T)
out2$which<-"multivariable no population"

amrfinder<-amrfinder_save
include<-out
amrfinder<-amrfinder[names(amrfinder) %in% include$term | names(amrfinder) %in% names(a) | names(amrfinder) =='mlst' ]
t<-as.formula(paste("Y", paste(names(amrfinder), collapse = " + "), sep=' ~ '))
model<-survreg(t, data=amrfinder,dist="gaussian")
out3<-broom::tidy(survreg(t, data=amrfinder,dist="gaussian"),conf.int=T)
out3$which<-"multivariable population"

cipro_r_genes<-filter(out3,2^(estimate) >1) %>% filter(!grepl('mlst',term))
cipro_line<-2^(log2(0.5) - model$coefficients[1])
cipro_line_genes<-filter(out3,2^(estimate) > cipro_line)
cipro_all_genes<-filter(out3,!grepl('mlst',term))


out$which<-"univariable"
out<-rbind(out,out2,out3)
count$gene<-janitor::make_clean_names(count$gene)

out<-filter(out,term %in% count$gene| grepl('mlst',term))



ceph_genes<-read_tsv('amrfinder.tsv') %>% filter(grepl('QUINOLONE',Class))
ceph_genes<-unique(ceph_genes$`Gene symbol`)
ceph_genes<-janitor::make_clean_names(ceph_genes)

a<-ifelse(out$term %in% ceph_genes,'red','black')


X<-c("mlstB1"  ,     "mlstA"    ,   "mlstC"  ,      "mlstD" ,       "mlstE",       
     "mlstF",        "mlstG" ,       "mlstother") 
out$term<-ifelse(out$term %in% X, str_replace(out$term,'mlst','phylogroup '),ifelse(grepl('mlst',out$term),str_replace(out$term,'mlst','ST '),out$term))


out<- out %>% 
  add_row(term = "GAP") %>% 
  add_row(term = "GAP2")

out$estimate<-2^out$estimate
out$conf.high<-2^out$conf.high
out$conf.low<-2^out$conf.low

out$which<-ifelse(grepl("GAP",out$term),'univariable',out$which)
out$which<-factor(out$which,levels=c("univariable","multivariable no population","multivariable population"))


amrfinder<-amrfinder_save1
m<-select(mlst,guuid,gene=mlst) %>% filter(!is.na(gene))
amrfinder<-rbind(amrfinder,m)
count<-amrfinder %>% group_by(gene) %>% distinct(guuid) %>%  count()
count$gene<-ifelse(str_length(count$gene) >3,janitor::make_clean_names(count$gene),count$gene)
count$gene <- case_when(count$gene == '131' ~ 'ST 131',
                        count$gene == '69' ~ 'ST 69',
                        count$gene == '95' ~ 'ST 95',
                        count$gene == '73' ~ 'ST 73',
                        count$gene == 'A' ~ 'phylogroup A',
                        count$gene == 'B1' ~ 'phylogroup B1',
                        count$gene == 'B2' ~ 'phylogroup B2',
                        count$gene == 'C' ~ 'phylogroup C',
                        count$gene == 'D' ~ 'phylogroup D',
                        count$gene == 'F' ~ 'phylogroup F',
                        count$gene == 'other' ~ 'phylogroup other',
                        TRUE ~ count$gene)


out<-left_join(out,count,by=c("term"="gene"))
out$term<-paste0(out$term, ' n=(',out$n,')')


out$term<-str_replace_all(out$term," n=[(]NA[)]",'')
out$term<-factor(out$term, levels=c("aac_6_ib_cr5 n=(126)","gyr_a_d87g n=(15)","gyr_a_d87n n=(361)","gyr_a_d87y n=(15)",      
                                    "gyr_a_s83a n=(11)","gyr_a_s83l n=(610)","par_c_e84g n=(11)","par_c_e84v n=(202)",
                                    "par_c_s57t n=(42)","par_c_s80i n=(380)","par_e_d475e n=(62)","par_e_i355t n=(73)",
                                    "par_e_i529l n=(374)","par_e_l416f n=(48)","par_e_s458a n=(49)","qnr_b19 n=(13)",
                                    "qnr_s1 n=(27)","GAP","GAP2","ST 131 n=(646)","ST 95 n=(463)","ST 73 n=(683)",
                                    "ST 69 n=(413)","phylogroup A n=(246)","phylogroup B1 n=(267)","phylogroup C n=(101)",   
                                    "phylogroup D n=(285)","phylogroup F n=(136)","phylogroup other n=(118)"))


cipro_plot<-ggplot(out) +
  aes(x=term,y=estimate,color=which,group=which) +
  geom_point(position = position_dodge(0.9)) + geom_errorbar(data=out,aes(ymin=conf.low,ymax=conf.high),position=position_dodge(width=0.9),width=0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90),panel.grid.major.x = element_blank()) +
  scale_y_log10() +
  xlab('Gene') + ylab('Fold-change in MIC') + 
  labs(color="Model") +
  ggtitle("Ciprofloxacin") + geom_hline(yintercept =2^(log2(0.5) - model$coefficients[1]),linetype="dashed",color="dark red") +
  geom_hline(yintercept=1)+
  scale_x_discrete(labels = function(x) ifelse(x == "GAP" | x=="GAP2", "", x)) +
  scale_color_viridis_d()

cipro<-out
cipro$drug<-"ciprofloxacin"


amrfinder<-amrfinder_save

amrfinder$guuid<-x
amrfinder<-select(amrfinder,-mlst)
amrfinder<-select(amrfinder,guuid,everything())
x<-names(amrfinder)
x<-x[!grepl('qnr',x)]
x<-x[!grepl('aac',x)]
amrfinder2<-select(amrfinder,x)
amrfinder$n<-rowSums(amrfinder2[,2:ncol(amrfinder2)])

amrfinder<-left_join(amrfinder,mlst,by=c("guuid"))
model<-survreg(Y ~ n + qnr_b19 + qnr_s1 + aac_6_ib_cr5  + mlst,data=amrfinder,dist="gaussian")


predictors <- c('n','qnr_b19','qnr_s1','aac_6_ib_cr5','mlst')
predictors<-ifelse(grepl("mlst",predictors),"mlst",predictors)
predictors<-unique(predictors)
predictors<-predictors[!is.na(predictors)  & ! predictors =="(Intercept)" & ! predictors == "Log(scale)"]

out=NULL



for (excluded_gene in predictors) {
  print(excluded_gene)
  other_predictors <- setdiff(predictors, excluded_gene)
  
  for (current_gene in other_predictors) {
    interaction_term <- paste(excluded_gene, current_gene, sep="*")
    formula_str <- as.formula(paste("Y ~", paste(setdiff(other_predictors,current_gene), collapse = "+"), "+", interaction_term))
    summary_model <- broom::tidy(survreg(formula_str, data = amrfinder, dist = "gaussian"))
    out <- rbind(out,data.frame(summary_model))
    
  }
}




replace_mlst_patterns <- function(texts) {
  texts <- gsub("mlst[A-Za-z0-9_]*:", "mlst:", texts)
  texts <- gsub(":[A-Za-z0-9_]*mlst[A-Za-z0-9_]*", ":mlst", texts)
  return(texts)
}



interaction_terms<-filter(out,p.value < 0.01) %>% filter(grepl(':',term))
interaction_terms$term<-replace_mlst_patterns(interaction_terms$term)
interaction_terms<-unique(interaction_terms$term)
interaction_terms


elimination<-  step(model, direction="backward")

selected_vars <- names(coef(elimination))
selected_vars<-selected_vars[!selected_vars=="(Intercept)"]
selected_vars<-selected_vars[!grepl('mlst',selected_vars)]
selected_vars<-c(selected_vars,"mlst")

interaction_string <- paste(interaction_terms, collapse = " + ")
new_formula <- as.formula(paste("Y ~ ", paste(selected_vars, collapse = " + ")))
model<-survreg(new_formula,data=amrfinder,dist="gaussian")
elimination<-step(model,direction="backward",k=3.8)

###leave one out cross validation


ceph_gene<-read_tsv('amrfinder.tsv') %>% filter(Subclass=='QUINOLONE') %>% filter(Scope=='core')

out<-NULL

for (i in 1:nrow(amrfinder)) {
  print(i)
  model<-survreg( Y[-i,] ~  n + qnr_b19 + qnr_s1 + mlst, data = amrfinder[-i,], 
                  dist = "gaussian")
  
  prediction <- predict(model, newdata = amrfinder[i, ])
  correct<-ifelse(prediction >= as.numeric(Y[i,1]) & prediction <= as.numeric(Y[i,2]),1,0)
  
  out<-rbind(out,data.frame(correct,prediction,as.numeric(Y[i,1]), as.numeric(Y[i,2]),amrfinder$guuid[i]))
}


out$amrfinder_prediction<-ifelse(out$amrfinder.guuid.i. %in% ceph_gene$Name,'R','S')
out$amrfinder_prediction<-as.factor(out$amrfinder_prediction)
#breakpoint >0.5 = log2(0.5) = -1
out$binary<-ifelse(out$as.numeric.Y.i..2.. >-1,'R','S')
out$binary_prediction<-ifelse(out$prediction>-1,'R','S')
out$binary<-as.factor(out$binary)
out$binary_prediction<-as.factor(out$binary_prediction)

extend_margin <- function(out, n){
  out$margin <- NA # Initialize margin column
  for(i in 0:n){
    condition <- out$prediction >= (out$as.numeric.Y.i..1.. - i) & 
      out$prediction <= (out$as.numeric.Y.i..2.. + i)
    out$margin[is.na(out$margin) & condition] <- i
  }
  return(out)
}


out2<-extend_margin(out,10)
out2$margin<-ifelse(out2$prediction < out2$as.numeric.Y.i..1..,out2$margin *-1,out2$margin)

out2$which<- ifelse(out2$binary_prediction==out2$binary & out2$binary=='S', 'Correct S',
                    ifelse(out2$binary_prediction==out2$binary & out2$binary =='R','Correct R',
                           ifelse(out2$binary_prediction=='S' & out2$binary =='R','VM',
                                  ifelse(out2$binary_prediction =='R' & out2$binary =='S','M',NA))))

out$within_1<-ifelse(out$prediction >= (out$as.numeric.Y.i..1.. -1) & out$prediction <= (out$as.numeric.Y.i..2.. +1),1,0)
out$within_2<-ifelse(out$prediction >= (out$as.numeric.Y.i..1.. -2) & out$prediction <= (out$as.numeric.Y.i..2.. +2),1,0)


caret::confusionMatrix(out$binary_prediction,out$binary)
table(out$correct)
prop.test(2600,(2600+237))
2600+237


table(out$within_1)
prop.test(2782,(2782+55))
2782+55

cm<-caret::confusionMatrix(as.factor(out$binary_prediction),reference=as.factor(out$binary),positive="R")
tn<-cm$table[2,2]
fn<-cm$table[2,1]
tp<-cm$table[1,1]
fp<-cm$table[1,2]
#concordance
prop.test(x=(tp+tn),n=(tp+tn+fp+fn))
tp+tn
tp+tn+fp+fn
#sens
prop.test(x=(tp),n=(tp+fn))
tp
tp+fn
#spec
prop.test(x=(tn),n=(tn+fp))
tn
tn+fp

#Maj
prop.test(x=(fp),n=(tn+fp))
fp
tn+fp
#Vmj
prop.test(x=(fn),n=(tp+fn))
fn
tp + fn

#PPV
prop.test(x=(tp),n=(tp+fp))
tp
tp+fp

#NPV
prop.test(x=(tn),n=(tn+fn))
tn
tn+fn


colors<-c("#09820B" ,"#03FF25" ,"#FFA703", "#ED3A37")
names(colors)<-c("Correct R","Correct S","M","VM")
cipro_residual<-ggplot(out2) +
  aes(x=margin, fill=which) +
  geom_bar(stat="count") + theme_minimal() +
  #scale_fill_manual(values=colors) + 
  scale_x_continuous(breaks = seq(min(out2$margin),max(out2$margin),1)) + 
  xlab("True - Predicted MIC") +
  labs(fill="") +
  ggtitle("Ciprofloxacin") +
  scale_fill_viridis_d()



amrfinder<-read_tsv('amrfinder.tsv')
cipro<-select(ec_match,guuid,Ciprofloxacin)  %>% filter(guuid %in% amrfinder$Name)

amrfinder<-amrfinder %>% filter(grepl('QUINOLONE',Class )) %>%  select(guuid=Name,gene=`Gene symbol`)


a_count<-amrfinder %>% group_by(guuid) %>% count()
count<-amrfinder %>% filter(str_length(guuid) ==36) %>%  group_by(gene) %>% count() %>% filter(n>=10)


deal_with_rare_genes <- function(gene_name) {
  if (startsWith(gene_name, "ampC")) {
    gsub("^ampC_.*", "ampC_other", gene_name)
  } else if (!grepl("_", gene_name)) {
    paste0(gene_name, "_other")
  } else {
    gsub("_[^_]+$", "_other", gene_name)
  }
}

amrfinder$gene <- ifelse(amrfinder$gene %in% count$gene,
                         amrfinder$gene,
                         sapply(amrfinder$gene, deal_with_rare_genes))

singletons<-amrfinder %>%  group_by(gene) %>% count() %>% filter(n<5)

amrfinder$gene<-ifelse(amrfinder$gene %in% singletons$gene,"other gene",amrfinder$gene)


cipro_r_genes<-filter(cipro_r_genes,p.value<0.05)
amrfinder$gene<-janitor::make_clean_names(amrfinder$gene,allow_dupes = T)
amrfinder_r<-filter(amrfinder,gene %in% cipro_r_genes$term)
amrfinder_r<- amrfinder_r %>% group_by(guuid) %>% count()

cipro$Ciprofloxacin <- case_when(
  cipro$Ciprofloxacin == '<=0.125' ~ '<=0.25',
  cipro$Ciprofloxacin == '0.25' ~ '<=0.25',
  TRUE ~ cipro$Ciprofloxacin
)

cipro$Ciprofloxacin<-factor(cipro$Ciprofloxacin,levels=c('<=0.25','0.5','1','>1'))
cipro<-filter(cipro,!is.na(Ciprofloxacin))
amrfinder_r<-left_join(cipro,amrfinder_r,by=c("guuid"="guuid"))
amrfinder_r$n<-ifelse(is.na(amrfinder_r$n),0,amrfinder_r$n)
amrfinder_r<-filter(amrfinder_r,!is.na(Ciprofloxacin))

r<-c('1','>1')

amrfinder_r$R<-ifelse(amrfinder_r$Ciprofloxacin %in% r,'R','S')
amrfinder_r$R<-factor(amrfinder_r$R,levels=c('S','R'))
cipro_gene_dist<-ggplot(amrfinder_r) +
  aes(x=n,fill=R) +
  geom_histogram(bins=length(unique(amrfinder_r$n))) +
  facet_wrap(~Ciprofloxacin,nrow = 1) + theme_minimal() +
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),strip.text = element_blank()) +
  xlab("Number of ARGs/mutations")+
  scale_fill_viridis_d() +
  labs(fill="")

ng<-amrfinder_r %>% 
  group_by(Ciprofloxacin) %>% 
  summarise(prop_no_gene = sum(n==0)/n(),
            lower=prop.test(sum(n==0),sum(n()))$conf.int[1],
            upper=prop.test(sum(n==0),sum(n()))$conf.int[2])


cipro_ng_plot<-ggplot(ng) +
  aes(x=Ciprofloxacin,y=prop_no_gene) +
  geom_point() +
  geom_errorbar(data=ng,aes(ymin=lower,ymax=upper),width=0) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank()) +
  ylab("Proportion with no ARG/mutation") +
  xlab("MIC")


cipro_eucast_plot<-ng_plot / cipro_gene_dist  +
  plot_annotation(title = "Ciprofloxacin")


#############eucast plot all################

(((amp_ng_plot + ggtitle("Ampicillin")) / amp_gene_dist) | ((coamox_ng_plot + ggtitle("Co-amoxiclav")) / coamox_gene_dist)) /
  (((piptaz_ng_plot + ggtitle("Piperacillin-Tazobactam")) / piptaz_gene_dist) | ((cefur_ng_plot + ggtitle("Cefuroxime")) / cefur_gene_dist)) /
  (((ceft_ng_plot + ggtitle("Ceftriaxone")) / ceft_gene_dist) | ((gent_ng_plot + ggtitle("Gentamicin")) / gent_gene_dist)) /
  (((cipro_ng_plot + ggtitle("Ciprofloxacin")) / cipro_gene_dist)  |((cotrim_ng_plot + ggtitle("Co-trimoxazole")) / cotrim_gene_dist)) +plot_layout(guides="collect",ncol=2,nrow=2) +
  plot_annotation(tag_levels = 'A')
  


#####################levo#############################


ec_match2<-select(ec_match,guuid,Levofloxacin,levofloxacin_lower,levofloxacin_upper) %>% filter(!is.na(levofloxacin_upper))
amrfinder<-read_tsv('amrfinder.tsv') 
ec_match2<-filter(ec_match2,guuid %in% amrfinder$Name)
amrfinder<-amrfinder %>% filter(grepl('QUINOLONE',Class )) %>%  select(guuid=Name,gene=`Gene symbol`)


a_count<-amrfinder %>% group_by(guuid) %>% count()
count<-amrfinder %>% filter(str_length(guuid) ==36) %>%  group_by(gene) %>% count() %>% filter(n>=10)

deal_with_rare_genes <- function(gene_name) {
  if (startsWith(gene_name, "ampC")) {
    gsub("^ampC_.*", "ampC_other", gene_name)
  } else if (!grepl("_", gene_name)) {
    paste0(gene_name, "_other")
  } else {
    gsub("_[^_]+$", "_other", gene_name)
  }
}

amrfinder$gene <- ifelse(amrfinder$gene %in% count$gene,
                         amrfinder$gene,
                         sapply(amrfinder$gene, deal_with_rare_genes))

singletons<-amrfinder %>%  group_by(gene) %>% count() %>% filter(n<5)

amrfinder$gene<-ifelse(amrfinder$gene %in% singletons$gene,"other gene",amrfinder$gene)
amrfinder_save1<-amrfinder



x<-data.frame(ec_match2$guuid)
x<-left_join(x,amrfinder,by=c("ec_match2.guuid"="guuid"))
x<-select(x,guuid=ec_match2.guuid,gene)
amrfinder<-data.frame(table(x)) %>% spread(gene,Freq)
names(amrfinder)<-str_replace_all(names(amrfinder),'-','_')
amrfinder[,2:ncol(amrfinder)][amrfinder[,2:ncol(amrfinder)]>1]<-1
genes<-names(amrfinder)
genes<-genes[!genes =='guuid']


amrfinder<-left_join(amrfinder,ec_match2,by=c("guuid"="guuid"))

mlst<-read_tsv('./mlst.tsv',col_names = F) %>% select(guuid=X1,spec=X2,mlst=X3)


big<-c('131','95','69','73')
mlst$mlst<-ifelse(mlst$mlst %in% big, mlst$mlst, 'other')

phylogroups<-read_tsv('./phylogroups.tsv')
mlst<-left_join(mlst,phylogroups,by=c("guuid"))
mlst$mlst<-ifelse(mlst$mlst =='other',mlst$phylogroup,mlst$mlst)
other<-c('U','U/cryptic','cryptic','E','G')
mlst$mlst<-ifelse(mlst$mlst %in% other,'other',mlst$mlst)

mlst<-select(mlst,guuid,mlst)

amrfinder<-left_join(amrfinder,mlst,by=c("guuid"))


out=NULL

names(amrfinder)<-janitor::make_clean_names(names(amrfinder))
genes<-janitor::make_clean_names(genes)


genes<-append(genes,"mlst")
amrfinder$mlst<-factor(amrfinder$mlst,levels=c('B2','131','95','73','69','A','B1','C','D','F','other'))
(Y <- with(amrfinder, Surv(levofloxacin_lower, levofloxacin_upper, type = "interval2")))
x<-amrfinder$guuid

amrfinder<-select(amrfinder,-guuid,-levofloxacin_lower,-levofloxacin_upper)

a<-amrfinder[,grep("^V",colnames(amrfinder))]
for(gene in genes){
  print(which(genes ==gene))
  
  print(gene)
  
  b<-select(amrfinder,gene)
  
  #c<-cbind(b,a)
  c<-b
  t<-as.formula(paste("Y", paste(names(c), collapse = " + "), sep=' ~ '))
  
  
  
  m<-broom::tidy(survreg(t ,data=c,dist="gaussian"),conf.int=T)
  
  m<-if(gene=='mlst'){m[2:11,]}else{m[2,]}    
  out=rbind(out,data.frame(m))
}
include<-out
amrfinder_save<-amrfinder
big<-c('st131','st95','st73','st69','other')
include<-filter(include,!gene %in% big)
amrfinder<-amrfinder[names(amrfinder) %in% include$term | names(amrfinder) %in% names(a) ]
t<-as.formula(paste("Y", paste(names(amrfinder), collapse = " + "), sep=' ~ '))
model<-survreg(t, data=amrfinder,dist="gaussian")
out2<-broom::tidy(survreg(t, data=amrfinder,dist="gaussian"),conf.int=T)
out2$which<-"multivariable no population"

amrfinder<-amrfinder_save
include<-out
amrfinder<-amrfinder[names(amrfinder) %in% include$term | names(amrfinder) %in% names(a) | names(amrfinder) =='mlst' ]
t<-as.formula(paste("Y", paste(names(amrfinder), collapse = " + "), sep=' ~ '))
model<-survreg(t, data=amrfinder,dist="gaussian")
out3<-broom::tidy(survreg(t, data=amrfinder,dist="gaussian"),conf.int=T)
out3$which<-"multivariable population"

cipro_r_genes<-filter(out3,round(estimate,1) >0) %>% filter(!grepl('mlst',term))
cipro_line<-log2(2) - model$coefficients[1]
cipro_line_genes<-filter(out3,estimate > cipro_line)
cipro_all_genes<-filter(out3,!grepl('mlst',term))


out$which<-"univariable"
out<-rbind(out,out2,out3)
#out<-left_join(out,out2,by=c("term"="term")) 

#out<-filter(out,!grepl('V',term))
count$gene<-janitor::make_clean_names(count$gene)

out<-filter(out,term %in% count$gene| grepl('mlst',term))

ceph_genes<-read_tsv('amrfinder.tsv') %>% filter(grepl('QUINOLONE',Class))
ceph_genes<-unique(ceph_genes$`Gene symbol`)
ceph_genes<-janitor::make_clean_names(ceph_genes)

a<-ifelse(out$term %in% ceph_genes,'red','black')


x<-c("mlstB1"  ,     "mlstA"    ,   "mlstC"  ,      "mlstD" ,       "mlstE",       
     "mlstF",        "mlstG" ,       "mlstother") 
out$term<-ifelse(out$term %in% x, str_replace(out$term,'mlst','phylogroup '),ifelse(grepl('mlst',out$term),str_replace(out$term,'mlst','ST '),out$term))


out<- out %>% 
  add_row(term = "GAP") %>% 
  add_row(term = "GAP2")

out$estimate<-2^out$estimate
out$conf.high<-2^out$conf.high
out$conf.low<-2^out$conf.low

out$which<-ifelse(grepl("GAP",out$term),'univariable',out$which)
out$which<-factor(out$which,levels=c("univariable","multivariable no population","multivariable population"))


amrfinder<-amrfinder_save1
m<-select(mlst,guuid,gene=mlst) %>% filter(!is.na(gene))
amrfinder<-rbind(amrfinder,m)
count<-amrfinder %>% group_by(gene) %>% distinct(guuid) %>%  count()
count$gene<-ifelse(str_length(count$gene) >3,janitor::make_clean_names(count$gene),count$gene)
count$gene <- case_when(count$gene == '131' ~ 'ST 131',
                        count$gene == '69' ~ 'ST 69',
                        count$gene == '95' ~ 'ST 95',
                        count$gene == '73' ~ 'ST 73',
                        count$gene == 'A' ~ 'phylogroup A',
                        count$gene == 'B1' ~ 'phylogroup B1',
                        count$gene == 'B2' ~ 'phylogroup B2',
                        count$gene == 'C' ~ 'phylogroup C',
                        count$gene == 'D' ~ 'phylogroup D',
                        count$gene == 'F' ~ 'phylogroup F',
                        count$gene == 'other' ~ 'phylogroup other',
                        TRUE ~ count$gene)


out<-left_join(out,count,by=c("term"="gene"))
out$term<-paste0(out$term, ' n=(',out$n,')')


out$term<-str_replace_all(out$term," n=[(]NA[)]",'')
out$term<-factor(out$term, levels=c("aac_6_ib_cr5 n=(126)","gyr_a_d87g n=(15)","gyr_a_d87n n=(361)","gyr_a_d87y n=(15)",      
                                    "gyr_a_s83a n=(11)","gyr_a_s83l n=(610)","par_c_e84g n=(11)","par_c_e84v n=(202)",
                                    "par_c_s57t n=(42)","par_c_s80i n=(380)","par_e_d475e n=(62)","par_e_i355t n=(73)",
                                    "par_e_i529l n=(374)","par_e_l416f n=(48)","par_e_s458a n=(49)","qnr_b19 n=(13)",
                                    "qnr_s1 n=(27)","GAP","GAP2","ST 131 n=(646)","ST 95 n=(463)","ST 73 n=(683)",
                                    "ST 69 n=(413)","phylogroup A n=(246)","phylogroup B1 n=(267)","phylogroup C n=(101)",   
                                    "phylogroup D n=(285)","phylogroup F n=(136)","phylogroup other n=(118)"))


levo_plot<-ggplot(out) +
  aes(x=term,y=estimate,color=which,group=which) +
  geom_point(position = position_dodge(0.9)) + geom_errorbar(data=out,aes(ymin=conf.low,ymax=conf.high),position=position_dodge(width=0.9),width=0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90),panel.grid.major.x = element_blank()) +
  xlab('Gene') + ylab('Fold-change in MIC') + 
  labs(color="Model") +
  ggtitle("Levofloxacin") + geom_hline(yintercept =2^(log2(1) - model$coefficients[1]),linetype="dashed",color="dark red") +
  geom_hline(yintercept=1)+
  scale_x_discrete(labels = function(x) ifelse(x == "GAP" | x=="GAP2", "", x)) +
  scale_color_viridis_d()

levo_out<-out

###################cipro 2##########################


ec_match2<-select(ec_match,guuid,Cipro2,Cipro2_lower,Cipro2_upper) %>% filter(!is.na(Cipro2_upper))
amrfinder<-read_tsv('amrfinder.tsv') 
ec_match2<-filter(ec_match2,guuid %in% amrfinder$Name)
amrfinder<-amrfinder %>% filter(grepl('QUINOLONE',Class )) %>%  select(guuid=Name,gene=`Gene symbol`)


a_count<-amrfinder %>% group_by(guuid) %>% count()
count<-amrfinder %>% filter(str_length(guuid) ==36) %>%  group_by(gene) %>% count() %>% filter(n>=10)

deal_with_rare_genes <- function(gene_name) {
  if (startsWith(gene_name, "ampC")) {
    gsub("^ampC_.*", "ampC_other", gene_name)
  } else if (!grepl("_", gene_name)) {
    paste0(gene_name, "_other")
  } else {
    gsub("_[^_]+$", "_other", gene_name)
  }
}

amrfinder$gene <- ifelse(amrfinder$gene %in% count$gene,
                         amrfinder$gene,
                         sapply(amrfinder$gene, deal_with_rare_genes))

singletons<-amrfinder %>%  group_by(gene) %>% count() %>% filter(n<5)

amrfinder$gene<-ifelse(amrfinder$gene %in% singletons$gene,"other gene",amrfinder$gene)
amrfinder_save1<-amrfinder



x<-data.frame(ec_match2$guuid)
x<-left_join(x,amrfinder,by=c("ec_match2.guuid"="guuid"))
x<-select(x,guuid=ec_match2.guuid,gene)
amrfinder<-data.frame(table(x)) %>% spread(gene,Freq)
names(amrfinder)<-str_replace_all(names(amrfinder),'-','_')
amrfinder[,2:ncol(amrfinder)][amrfinder[,2:ncol(amrfinder)]>1]<-1
genes<-names(amrfinder)
genes<-genes[!genes =='guuid']


amrfinder<-left_join(amrfinder,ec_match2,by=c("guuid"="guuid"))

mlst<-read_tsv('./mlst.tsv',col_names = F) %>% select(guuid=X1,spec=X2,mlst=X3)


big<-c('131','95','69','73')
mlst$mlst<-ifelse(mlst$mlst %in% big, mlst$mlst, 'other')

phylogroups<-read_tsv('./phylogroups.tsv')
mlst<-left_join(mlst,phylogroups,by=c("guuid"))
mlst$mlst<-ifelse(mlst$mlst =='other',mlst$phylogroup,mlst$mlst)
other<-c('U','U/cryptic','cryptic','E','G')
mlst$mlst<-ifelse(mlst$mlst %in% other,'other',mlst$mlst)

mlst<-select(mlst,guuid,mlst)

amrfinder<-left_join(amrfinder,mlst,by=c("guuid"))


out=NULL

names(amrfinder)<-janitor::make_clean_names(names(amrfinder))
genes<-janitor::make_clean_names(genes)


genes<-append(genes,"mlst")
amrfinder$mlst<-factor(amrfinder$mlst,levels=c('B2','131','95','73','69','A','B1','C','D','F','other'))
(Y <- with(amrfinder, Surv(cipro2_lower, cipro2_upper, type = "interval2")))
x<-amrfinder$guuid

amrfinder<-select(amrfinder,-guuid,-cipro2_lower,-cipro2_upper)

a<-amrfinder[,grep("^V",colnames(amrfinder))]
for(gene in genes){
  print(which(genes ==gene))
  
  print(gene)
  
  b<-select(amrfinder,gene)
  
  #c<-cbind(b,a)
  c<-b
  t<-as.formula(paste("Y", paste(names(c), collapse = " + "), sep=' ~ '))
  
  
  
  m<-broom::tidy(survreg(t ,data=c,dist="gaussian"),conf.int=T)
  
  m<-if(gene=='mlst'){m[2:11,]}else{m[2,]}    
  out=rbind(out,data.frame(m))
}
include<-out
amrfinder_save<-amrfinder
big<-c('st131','st95','st73','st69','other')
include<-filter(include,!gene %in% big)
amrfinder<-amrfinder[names(amrfinder) %in% include$term | names(amrfinder) %in% names(a) ]
t<-as.formula(paste("Y", paste(names(amrfinder), collapse = " + "), sep=' ~ '))
model<-survreg(t, data=amrfinder,dist="gaussian")
out2<-broom::tidy(survreg(t, data=amrfinder,dist="gaussian"),conf.int=T)
out2$which<-"multivariable no population"

amrfinder<-amrfinder_save
include<-out
amrfinder<-amrfinder[names(amrfinder) %in% include$term | names(amrfinder) %in% names(a) | names(amrfinder) =='mlst' ]
t<-as.formula(paste("Y", paste(names(amrfinder), collapse = " + "), sep=' ~ '))
model<-survreg(t, data=amrfinder,dist="gaussian")
out3<-broom::tidy(survreg(t, data=amrfinder,dist="gaussian"),conf.int=T)
out3$which<-"multivariable population"

cipro_r_genes<-filter(out3,round(estimate,1) >0) %>% filter(!grepl('mlst',term))
cipro_line<-log2(2) - model$coefficients[1]
cipro_line_genes<-filter(out3,estimate > cipro_line)
cipro_all_genes<-filter(out3,!grepl('mlst',term))


out$which<-"univariable"
out<-rbind(out,out2,out3)
#out<-left_join(out,out2,by=c("term"="term")) 

#out<-filter(out,!grepl('V',term))
count$gene<-janitor::make_clean_names(count$gene)

out<-filter(out,term %in% count$gene| grepl('mlst',term))

ceph_genes<-read_tsv('amrfinder.tsv') %>% filter(grepl('QUINOLONE',Class))
ceph_genes<-unique(ceph_genes$`Gene symbol`)
ceph_genes<-janitor::make_clean_names(ceph_genes)

a<-ifelse(out$term %in% ceph_genes,'red','black')


x<-c("mlstB1"  ,     "mlstA"    ,   "mlstC"  ,      "mlstD" ,       "mlstE",       
     "mlstF",        "mlstG" ,       "mlstother") 
out$term<-ifelse(out$term %in% x, str_replace(out$term,'mlst','phylogroup '),ifelse(grepl('mlst',out$term),str_replace(out$term,'mlst','ST '),out$term))


out<- out %>% 
  add_row(term = "GAP") %>% 
  add_row(term = "GAP2")

out$estimate<-2^out$estimate
out$conf.high<-2^out$conf.high
out$conf.low<-2^out$conf.low

out$which<-ifelse(grepl("GAP",out$term),'univariable',out$which)
out$which<-factor(out$which,levels=c("univariable","multivariable no population","multivariable population"))


amrfinder<-amrfinder_save1
m<-select(mlst,guuid,gene=mlst) %>% filter(!is.na(gene))
amrfinder<-rbind(amrfinder,m)
count<-amrfinder %>% group_by(gene) %>% distinct(guuid) %>%  count()
count$gene<-ifelse(str_length(count$gene) >3,janitor::make_clean_names(count$gene),count$gene)
count$gene <- case_when(count$gene == '131' ~ 'ST 131',
                        count$gene == '69' ~ 'ST 69',
                        count$gene == '95' ~ 'ST 95',
                        count$gene == '73' ~ 'ST 73',
                        count$gene == 'A' ~ 'phylogroup A',
                        count$gene == 'B1' ~ 'phylogroup B1',
                        count$gene == 'B2' ~ 'phylogroup B2',
                        count$gene == 'C' ~ 'phylogroup C',
                        count$gene == 'D' ~ 'phylogroup D',
                        count$gene == 'F' ~ 'phylogroup F',
                        count$gene == 'other' ~ 'phylogroup other',
                        TRUE ~ count$gene)


out<-left_join(out,count,by=c("term"="gene"))
out$term<-paste0(out$term, ' n=(',out$n,')')


out$term<-str_replace_all(out$term," n=[(]NA[)]",'')
out$term<-factor(out$term, levels=c("aac_6_ib_cr5 n=(126)","gyr_a_d87g n=(15)","gyr_a_d87n n=(361)","gyr_a_d87y n=(15)",      
                                    "gyr_a_s83a n=(11)","gyr_a_s83l n=(610)","par_c_e84g n=(11)","par_c_e84v n=(202)",
                                    "par_c_s57t n=(42)","par_c_s80i n=(380)","par_e_d475e n=(62)","par_e_i355t n=(73)",
                                    "par_e_i529l n=(374)","par_e_l416f n=(48)","par_e_s458a n=(49)","qnr_b19 n=(13)",
                                    "qnr_s1 n=(27)","GAP","GAP2","ST 131 n=(646)","ST 95 n=(463)","ST 73 n=(683)",
                                    "ST 69 n=(413)","phylogroup A n=(246)","phylogroup B1 n=(267)","phylogroup C n=(101)",   
                                    "phylogroup D n=(285)","phylogroup F n=(136)","phylogroup other n=(118)"))


cipro2_plot<-ggplot(out) +
  aes(x=term,y=estimate,color=which,group=which) +
  geom_point(position = position_dodge(0.9)) + geom_errorbar(data=out,aes(ymin=conf.low,ymax=conf.high),position=position_dodge(width=0.9),width=0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90),panel.grid.major.x = element_blank()) +
  xlab('Gene') + ylab('Fold-change in MIC') + 
  labs(color="Model") +
  ggtitle("Cipro2") + geom_hline(yintercept =2^(log2(0.5) - model$coefficients[1]),linetype="dashed",color="dark red") +
  geom_hline(yintercept=1)+
  scale_x_discrete(labels = function(x) ifelse(x == "GAP" | x=="GAP2", "", x)) +
  scale_color_viridis_d() +
  ylim(0,30)

cipro2<-out


levo_out<-filter(levo_out,which=="multivariable population")
cipro2<-filter(cipro2,which=="multivariable population")
levo_out$which<-"levo"
cipro2$which<-"cipro"

quin<-left_join(cipro2,levo_out,by=c("term"))
quin$term<-as.character(quin$term)

sig<-c('qnr_s1 n=(27)','gyr_a_s83l n=(610)')
quin$sig<-ifelse(quin$term %in% sig ,1,0)
quin$alpha<-ifelse(quin$sig==1,1,0.4)
quin$term2<-ifelse(quin$sig==1,quin$term,'other')


quin$term2<-str_replace_all(quin$term2,' n=.*','')
quin$term2<-ifelse(quin$term2=='gyr_a_s83l','gyrA-S83L',ifelse(quin$term2=='qnr_s1','qnrS1','other'))
quin$term2<-factor(quin$term2,levels=c('gyrA-S83L','qnrS1','other'))
cols<-c( "#2A003D", "#337D57", "#080606")
cip_effects<-ggplot() +
  geom_point(data=quin,aes(x=estimate.x,y=estimate.y,color=term2,alpha=alpha)) +
  geom_abline(intercept = 0, slope=1, linetype="dashed") +
  geom_errorbar(data=quin,aes(xmin=conf.low.x,x=estimate.x,y=estimate.y,xmax=conf.high.x,alpha=alpha,color=term2)) +
  geom_errorbar(data=quin,aes(x=estimate.x,y=estimate.y,ymin=conf.low.y,ymax=conf.high.y,alpha=alpha,color=term2)) +
  theme_minimal() +
  xlab("Fold change in Ciprofloxacin MIC") +
  ylab("Fold change in Levofloxacin MIC") +
  scale_color_manual(values=cols) +
  guides(alpha=FALSE) +
  labs(color="Gene")

cip_effects /cipro_plot + plot_annotation(tag_levels="A")
  
#########Fig S3############

amp_residual + coamox_residual + piptaz_residual +cefuroxime_residual + ceftriaxone_residual + cefepime_residual + gent_residual +
   cipro_residual + cotrim_residual +
   plot_layout(guides="collect")


#########table S1#########

ts1<-rbind(gent, ceftriaxone, cipro, ampicillin, coamoxiclav, cefuroxime, cotrim,piptaz)


ts1<-filter(ts1,which== "multivariable population")
ts1$"95% CI"<-paste0(round(ts1$conf.low,2)," - ",round(ts1$conf.high,2))
ts1<-select(ts1,Gene=term,"Estimated log2(change in MIC)"=estimate,"p value"=p.value,`95% CI`,Antibiotic=drug)

x<-c("mlstB1"  ,     "mlstA"    ,   "mlstC"  ,      "mlstD" ,       "mlstE",       
     "mlstF",        "mlstG" ,       "mlstother") 
ts1$Gene<-as.character(ts1$Gene)
ts1$Gene<-ifelse(ts1$Gene %in% x, str_replace(ts1$Gene,'mlst','phylogroup '),ifelse(grepl('mlst',ts1$Gene),str_replace(ts1$Gene,'mlst','ST '),ts1$Gene))
ts1$`Estimated log2(change in MIC)`<-round(ts1$`Estimated log2(change in MIC)`,2)

format_pvalues <- function(pvals) {
  sapply(pvals, function(p) {
    if (p < 0.001) {
      return("<0.001")
    } else {
      return(formatC(p, format = "f", digits = 2))
    }
  })
}
ts1$`p value`<-format_pvalues(ts1$`p value`)
ts1<-select(ts1,Gene,`Estimated log2(change in MIC)`,`95% CI`,`p value`,Antibiotic)
ts1 <- ts1 %>% kableExtra::kable() 
write_lines(ts1,'ts1.html')



#######MIC to rule out ARG#######

######coamox##############

coamox_r_genes<-filter(coamox_r_genes,p.value<0.05)
amrfinder<-read_tsv('amrfinder.tsv') 
tem_promotor_mutations<-filter(amrfinder,grepl('TEMp',`Gene symbol`))
ec_match2<-filter(ec_match2,guuid %in% amrfinder$Name)
amrfinder<-amrfinder%>% filter(Class=='BETA-LACTAM' ) %>%  select(guuid=Name,gene=`Gene symbol`)
length(unique(amrfinder$gene))
amrfinder$gene<-str_replace_all(amrfinder$gene,'-','_')

amrfinder$gene<-ifelse(grepl('blaTEM_',amrfinder$gene) & amrfinder$guuid %in% tem_promotor_mutations$Name,
                       paste0(amrfinder$gene,'_',tem_promotor_mutations$`Gene symbol`[match(amrfinder$guuid,tem_promotor_mutations$Name)]),amrfinder$gene)
amrfinder <- amrfinder %>%
  filter(!grepl("(?<!\\w)blaTEMp_", gene, perl = TRUE))
coamox<-read_tsv('amrfinder.tsv')  %>% filter(Name %in% ec_match2$guuid) %>% filter(grepl('BETA-LACTAM',Class )) %>%  select(guuid=Name,gene=`Gene symbol`) 
coamox$gene<-str_replace_all(coamox$gene,'-','_')

count<-amrfinder %>% filter(str_length(guuid) ==36) %>%  group_by(gene) %>% count() %>% filter(n>=10)



deal_with_rare_genes <- function(gene_name) {
  if (startsWith(gene_name, "ampC")) {
    gsub("^ampC_.*", "ampC_other", gene_name)
  } else if (!grepl("_", gene_name)) {
    paste0(gene_name, "_other")
  } else {
    gsub("_[^_]+$", "_other", gene_name)
  }
}

amrfinder$gene <- ifelse(amrfinder$gene %in% count$gene,
                         amrfinder$gene,
                         sapply(amrfinder$gene, deal_with_rare_genes))

singletons<-amrfinder %>%  group_by(gene) %>% count() %>% filter(n<5)

amrfinder$gene<-ifelse(amrfinder$gene %in% singletons$gene,"other gene",amrfinder$gene)
amrfinder$gene<-janitor::make_clean_names(amrfinder$gene)
amrfinder$gene<- str_remove(amrfinder$gene, "_\\d+$")


amrfinder<-filter(amrfinder,gene %in% coamox_r_genes$term)
coamox_r_genes<-select(coamox_r_genes,term,estimate)
coamox_r_genes$estimate<-2^coamox_r_genes$estimate
amrfinder<-left_join(amrfinder,coamox_r_genes,by=c("gene"="term"))

pheno<-select(ec_match,guuid,Coamox) %>% filter(!is.na(Coamox)) %>% filter(!Coamox =='>8/2')
coamox_amrfinder<-filter(amrfinder,guuid %in%pheno$guuid)

pheno$Coamox<-as.numeric(factor(pheno$Coamox,levels=c('<=2/2','4/2','8/2','16/2','32/2','>32/2')))
mic_bands<-unique(pheno$Coamox)

genes<-unique(coamox_amrfinder$gene)

out=NULL

for(level in 6:1){
  pheno2<-filter(pheno,Coamox==level)
  pheno2<-left_join(pheno2,coamox_amrfinder,by=c("guuid"="guuid"))
  pheno2$estimate<-ifelse(is.na(pheno2$estimate),0,pheno2$estimate)
  pheno2<-pheno2 %>% group_by(guuid) %>% summarise(estimate=sum(estimate))
  

  
  out<-rbind(out,data.frame(level,pheno2))
}


names(out)<-c('MIC','guuid','estimate')



out$MIC<-case_when(out$MIC ==1 ~ '<=2/2',
                    out$MIC ==2 ~ '4/2',
                    out$MIC ==3 ~ '8/2',
                    out$MIC ==4 ~ '16/2',
                    out$MIC ==5 ~ '32/2',
                    out$MIC ==6 ~ '>32/2')

out$MIC<-factor(out$MIC,levels=c('<=2/2','4/2','8/2','16/2','32/2','>32/2'))

breaks <- seq(from = 0, to = max(out$estimate) + 1, by = 1)
out$estimate_bin<-cut(out$estimate,breaks = breaks, right = FALSE, labels = paste(head(breaks, -1), head(breaks, -1) +1 , sep="-"))

out<-out %>% group_by(MIC) %>% count(estimate_bin)
first_colors <- colorRampPalette(c("green", "darkgreen"))(3)
rest_colors <- colorRampPalette(c( "#F7BDBC", "#FA0A0A"))(length(unique(out$estimate_bin)) - 3)
color_palette <- c(first_colors, rest_colors)




line_position <- which(levels(out$MIC) == "8/2") + 0.5

sens<-c('<=2/2','4/2','8/2')

out$class<-ifelse(out$estimate_bin=='0-1'| out$estimate_bin=='1-2' | out$estimate_bin =='2-3','S','R')

out <- out %>%
  mutate(mid_point = sapply(strsplit(as.character(estimate_bin), "-"), function(x) {
    mean(as.numeric(x))
  }))


mic_coamox<-ggplot(out) +
  aes(x=MIC,y=n,fill=mid_point) +
  geom_bar(stat='identity')  +
  geom_vline(aes(xintercept = line_position), linetype="dotted", color = "black") +
  theme_minimal() +
  labs(fill="Total estimated\nfold-change in MIC") +
  scale_fill_gradientn(colors=color_palette)
  



calculate_ratio <- function(out, sens, conf.level = 0.95) {
  ratios <- out %>%
    group_by(MIC) %>%
    summarise(
      R_count = sum(n[class == "R"]),
      S_count = sum(n[class == "S"]),
      total = R_count + S_count,
      ratio = ifelse(MIC %in% sens,
                     R_count / total,
                     S_count / total)
    )
  
  # Compute the confidence intervals
  conf_int <- mapply(function(x, n) {
    test <- prop.test(x, n, conf.level = conf.level)
    c(lower = test$conf.int[1], upper = test$conf.int[2])
  },
  x = ifelse(ratios$MIC %in% sens, ratios$R_count, ratios$S_count),
  n = ratios$total,
  SIMPLIFY = FALSE
  )
  
  # Add the confidence intervals to the ratios
  ratios$lower_bound <- sapply(conf_int, "[[", "lower")
  ratios$upper_bound <- sapply(conf_int, "[[", "upper")
  
  # Join with the original data
  ratios %>%
    select(-R_count, -S_count, -total) %>%
    left_join(out, by = "MIC") %>%
    select(MIC, ratio, lower_bound, upper_bound) %>%
    distinct()
}

out2<-calculate_ratio(out,sens) 
out2$error<-ifelse(out2$MIC %in% sens,'Erroneous identification\nas wild-type','Erroneous identification\nas non wild-type')



margin_plot <- ggplot(out2) +
  aes(x=MIC, y=ratio*100, ymin=lower_bound*100, ymax=upper_bound*100,color=error) +
  geom_point() +
  geom_errorbar(width=0) +
  theme_minimal() +
  ylab("% Error") +
  xlab("") + 
  ggtitle("Co-amoxiclav") +
  labs(color='Error type')

combined_plot<-margin_plot/mic_coamox 
co_amox_combined_plot <- combined_plot + plot_layout(heights = c(1, 3))
#######gent#######

gent_r_genes<-filter(gent_r_genes,p.value<0.05)
amrfinder<-read_tsv('amrfinder.tsv')
amrfinder<-amrfinder%>% filter(grepl('AMINOGLYCOSIDE',Class)) %>%  select(guuid=Name,gene=`Gene symbol`)
amrfinder$gene<-str_replace_all(amrfinder$gene,'-','_')

count<-amrfinder %>% filter(str_length(guuid) ==36) %>%  group_by(gene) %>% count() %>% filter(n>=10)


deal_with_rare_genes <- function(gene_name) {
  if (startsWith(gene_name, "ampC")) {
    gsub("^ampC_.*", "ampC_other", gene_name)
  } else if (!grepl("_", gene_name)) {
    paste0(gene_name, "_other")
  } else {
    gsub("_[^_]+$", "_other", gene_name)
  }
}

amrfinder$gene <- ifelse(amrfinder$gene %in% count$gene,
                         amrfinder$gene,
                         sapply(amrfinder$gene, deal_with_rare_genes))

singletons<-amrfinder %>%  group_by(gene) %>% count() %>% filter(n<5)

amrfinder$gene<-ifelse(amrfinder$gene %in% singletons$gene,"other gene",amrfinder$gene)
amrfinder$gene<-janitor::make_clean_names(amrfinder$gene)
amrfinder$gene<- str_remove(amrfinder$gene, "_\\d+$")


amrfinder<-filter(amrfinder,gene %in% gent_r_genes$term)
gent_r_genes<-select(gent_r_genes,term,estimate)
gent_r_genes$estimate<-2^gent_r_genes$estimate
amrfinder<-left_join(amrfinder,gent_r_genes,by=c("gene"="term"))
#for each MIC interval what is the performance as a diagnostic test to exclude the ARG (e.g. confirm wild type)

pheno<-select(ec_match,guuid,Gentamicin) %>% filter(!is.na(Gentamicin)) 
coamox_amrfinder<-filter(amrfinder,guuid %in%pheno$guuid)

pheno$Gentamicin<-as.numeric(factor(pheno$Gentamicin,levels=c('<=1','2','4','>4')))
mic_bands<-unique(pheno$Gentamicin)

genes<-unique(coamox_amrfinder$gene)
out=NULL

for(level in length(mic_bands):1){
  pheno2<-filter(pheno,Gentamicin==level)
  pheno2<-left_join(pheno2,coamox_amrfinder,by=c("guuid"="guuid"))
  pheno2$estimate<-ifelse(is.na(pheno2$estimate),0,pheno2$estimate)
  pheno2<-pheno2 %>% group_by(guuid) %>% summarise(estimate=sum(estimate))
 
  #pheno2$estimate<-ifelse(pheno2$estimate ==0,'0',paste0(' >',round(pheno2$estimate),' - ','<=',(round(pheno2$estimate) +1)))
  
  
  out<-rbind(out,data.frame(level,pheno2))
  }


names(out)<-c('MIC','guuid','estimate')

out$MIC<-case_when(out$MIC ==1 ~ '<=1',
                    out$MIC ==2 ~ '2',
                    out$MIC ==3 ~ '4',
                    out$MIC ==4 ~ '>4')

out$MIC<-factor(out$MIC,levels=c('<=1','2','4','>4'))

breaks <- seq(from = 0, to = max(out$estimate) + 1, by = 1)
out$estimate_bin<-cut(out$estimate,breaks = breaks, right = FALSE, labels = paste(head(breaks, -1), head(breaks, -1) +1 , sep="-"))

out<-out %>% group_by(MIC) %>% count(estimate_bin)
first_colors <- colorRampPalette(c("green", "darkgreen"))(4)
rest_colors <- colorRampPalette(c( "#F7BDBC", "#FA0A0A"))(length(unique(out$estimate_bin)) - 4)
color_palette <- c(first_colors, rest_colors)




line_position <- which(levels(out$MIC) == "2") + 0.5

sens<-c('<=1','2')

out$class<-ifelse(out$estimate_bin=='0-1'| out$estimate_bin=='1-2' | out$estimate_bin =='2-3' | out$estimate_bin=='3-4','S','R')

out <- out %>%
  mutate(mid_point = sapply(strsplit(as.character(estimate_bin), "-"), function(x) {
    mean(as.numeric(x))
  }))


mic_gent<-ggplot(out) +
  aes(x=MIC,y=n,fill=mid_point) +
  geom_bar(stat='identity')  +
  geom_vline(aes(xintercept = line_position), linetype="dotted", color = "black") +
  theme_minimal() +
  labs(fill="Total estimated\nfold-change in MIC") +
  scale_fill_gradientn(colors=color_palette)


calculate_ratio <- function(out, sens, conf.level = 0.95) {
  ratios <- out %>%
    group_by(MIC) %>%
    summarise(
      R_count = sum(n[class == "R"]),
      S_count = sum(n[class == "S"]),
      total = R_count + S_count,
      ratio = ifelse(MIC %in% sens,
                     R_count / total,
                     S_count / total)
    )
  
  # Compute the confidence intervals
  conf_int <- mapply(function(x, n) {
    test <- prop.test(x, n, conf.level = conf.level)
    c(lower = test$conf.int[1], upper = test$conf.int[2])
  },
  x = ifelse(ratios$MIC %in% sens, ratios$R_count, ratios$S_count),
  n = ratios$total,
  SIMPLIFY = FALSE
  )
  
  # Add the confidence intervals to the ratios
  ratios$lower_bound <- sapply(conf_int, "[[", "lower")
  ratios$upper_bound <- sapply(conf_int, "[[", "upper")
  
  # Join with the original data
  ratios %>%
    select(-R_count, -S_count, -total) %>%
    left_join(out, by = "MIC") %>%
    select(MIC, ratio, lower_bound, upper_bound) %>%
    distinct()
}

out2<-calculate_ratio(out,sens) 
out2$error<-ifelse(out2$MIC %in% sens,'Erroneous identification\nas wild-type','Erroneous identification\nas non wild-type')


margin_plot <- ggplot(out2) +
  aes(x=MIC, y=ratio*100, ymin=lower_bound*100, ymax=upper_bound*100,color=error) +
  geom_point() +
  geom_errorbar(width=0) +
  theme_minimal() +
  ylab("% Error") +
  xlab("") + 
  ggtitle("Gentamicin") +
  labs(color='Error type')

combined_plot<-margin_plot/mic_gent
gent_combined_plot <- combined_plot + plot_layout(heights = c(1, 3))
########ceftriaxone##########
ceftriaxone_r_genes<-filter(ceftriaxone_r_genes,p.value<0.05)

amrfinder<-read_tsv('amrfinder.tsv') 
tem_promotor_mutations<-filter(amrfinder,grepl('TEMp',`Gene symbol`))
ec_match2<-filter(ec_match2,guuid %in% amrfinder$Name)
amrfinder<-amrfinder%>% filter(Class=='BETA-LACTAM' ) %>%  select(guuid=Name,gene=`Gene symbol`)
length(unique(amrfinder$gene))
amrfinder$gene<-str_replace_all(amrfinder$gene,'-','_')

amrfinder$gene<-ifelse(grepl('blaTEM_',amrfinder$gene) & amrfinder$guuid %in% tem_promotor_mutations$Name,
                       paste0(amrfinder$gene,'_',tem_promotor_mutations$`Gene symbol`[match(amrfinder$guuid,tem_promotor_mutations$Name)]),amrfinder$gene)
amrfinder <- amrfinder %>%
  filter(!grepl("(?<!\\w)blaTEMp_", gene, perl = TRUE))
coamox<-read_tsv('amrfinder.tsv')  %>% filter(Name %in% ec_match2$guuid) %>% filter(grepl('BETA-LACTAM',Class )) %>%  select(guuid=Name,gene=`Gene symbol`) 
coamox$gene<-str_replace_all(coamox$gene,'-','_')

count<-amrfinder %>% filter(str_length(guuid) ==36) %>%  group_by(gene) %>% count() %>% filter(n>=10)


deal_with_rare_genes <- function(gene_name) {
  if (startsWith(gene_name, "ampC")) {
    gsub("^ampC_.*", "ampC_other", gene_name)
  } else if (!grepl("_", gene_name)) {
    paste0(gene_name, "_other")
  } else {
    gsub("_[^_]+$", "_other", gene_name)
  }
}

amrfinder$gene <- ifelse(amrfinder$gene %in% count$gene,
                         amrfinder$gene,
                         sapply(amrfinder$gene, deal_with_rare_genes))

singletons<-amrfinder %>%  group_by(gene) %>% count() %>% filter(n<5)

amrfinder$gene<-ifelse(amrfinder$gene %in% singletons$gene,"other gene",amrfinder$gene)
amrfinder$gene<-janitor::make_clean_names(amrfinder$gene)
amrfinder$gene<- str_remove(amrfinder$gene, "_\\d+$")


amrfinder<-filter(amrfinder,gene %in% ceftriaxone_r_genes$term)
gent_r_genes<-select(ceftriaxone_r_genes,term,estimate)
ceftriaxone_r_genes$estimate<-2^ceftriaxone_r_genes$estimate
amrfinder<-left_join(amrfinder,ceftriaxone_r_genes,by=c("gene"="term"))


pheno<-select(ec_match,guuid,Ceftriaxone) %>% filter(!is.na(Ceftriaxone)) 
pheno$Ceftriaxone<-ifelse(pheno$Ceftriaxone=='<=0.5','<=1',ifelse(pheno$Ceftriaxone =='1','<=1',pheno$Ceftriaxone))
coamox_amrfinder<-filter(amrfinder,guuid %in%pheno$guuid)

pheno$Ceftriaxone<-as.numeric(factor(pheno$Ceftriaxone,levels=c('<=1','2','4','>4')))
mic_bands<-unique(pheno$Ceftriaxone)

out=NULL
for(level in length(mic_bands):1){
  pheno2<-filter(pheno,Ceftriaxone==level)
  pheno2<-left_join(pheno2,coamox_amrfinder,by=c("guuid"="guuid"))
  pheno2$estimate<-ifelse(is.na(pheno2$estimate),0,pheno2$estimate)
  pheno2<-pheno2 %>% group_by(guuid) %>% summarise(estimate=sum(estimate))
  
  #pheno2$estimate<-ifelse(pheno2$estimate ==0,'0',paste0(' >',round(pheno2$estimate),' - ','<=',(round(pheno2$estimate) +1)))
  
  
  out<-rbind(out,data.frame(level,pheno2))
}


names(out)<-c('MIC','guuid','estimate')

out$MIC<-case_when(out$MIC ==1 ~ '<=1',
                   out$MIC ==2 ~ '2',
                   out$MIC ==3 ~ '4',
                   out$MIC ==4 ~ '>4')

out$MIC<-factor(out$MIC,levels=c('<=1','2','4','>4'))

breaks <- seq(from = 0, to = max(out$estimate) + 1, by = 1)
out$estimate_bin<-cut(out$estimate,breaks = breaks, right = FALSE, labels = paste(head(breaks, -1), head(breaks, -1) +1 , sep="-"))

out<-out %>% group_by(MIC) %>% count(estimate_bin)
first_colors <- colorRampPalette(c("green", "darkgreen"))(10)
rest_colors <- colorRampPalette(c( "#F7BDBC", "#FA0A0A"))(length(unique(levels(out$estimate_bin))) - 10)
color_palette <- c(first_colors, rest_colors)



line_position <- which(levels(out$MIC) == "2") + 0.5


sens<-c('<=1','2')

out$class<-ifelse(out$estimate_bin=='0-1'|out$estimate_bin=='1-2' | out$estimate_bin=='2-3'|
                  out$estimate_bin=='3-4'|out$estimate_bin=='4-5'|out$estimate_bin=='5-6'|
                    out$estimate_bin=='6-7'|out$estimate_bin=='7-8'|out$estimate_bin=='8-9'|out$estimate_bin=='9-10','S','R')

out <- out %>%
  mutate(mid_point = sapply(strsplit(as.character(estimate_bin), "-"), function(x) {
    mean(as.numeric(x))
  }))

mic_ceftriaxone<-ggplot(out) +
  aes(x=MIC,y=n,fill=mid_point) +
  geom_bar(stat='identity')  +
  geom_vline(aes(xintercept = line_position), linetype="dotted", color = "black") +
  theme_minimal() +
  labs(fill="Total estimated\nfold-change in MIC") +
  scale_fill_gradientn(colors=color_palette)



calculate_ratio <- function(out, sens, conf.level = 0.95) {
  ratios <- out %>%
    group_by(MIC) %>%
    summarise(
      R_count = sum(n[class == "R"]),
      S_count = sum(n[class == "S"]),
      total = R_count + S_count,
      ratio = ifelse(MIC %in% sens,
                     R_count / total,
                     S_count / total)
    )
  
  # Compute the confidence intervals
  conf_int <- mapply(function(x, n) {
    test <- prop.test(x, n, conf.level = conf.level)
    c(lower = test$conf.int[1], upper = test$conf.int[2])
  },
  x = ifelse(ratios$MIC %in% sens, ratios$R_count, ratios$S_count),
  n = ratios$total,
  SIMPLIFY = FALSE
  )
  
  # Add the confidence intervals to the ratios
  ratios$lower_bound <- sapply(conf_int, "[[", "lower")
  ratios$upper_bound <- sapply(conf_int, "[[", "upper")
  
  # Join with the original data
  ratios %>%
    select(-R_count, -S_count, -total) %>%
    left_join(out, by = "MIC") %>%
    select(MIC, ratio, lower_bound, upper_bound) %>%
    distinct()
}

out2<-calculate_ratio(out,sens) 
out2$error<-ifelse(out2$MIC %in% sens,'Erroneous identification\nas wild-type','Erroneous identification\nas non wild-type')


margin_plot <- ggplot(out2) +
  aes(x=MIC, y=ratio*100, ymin=lower_bound*100, ymax=upper_bound*100,color=error) +
  geom_point() +
  geom_errorbar(width=0) +
  theme_minimal() +
  ylab("% Error") +
  xlab("") + 
  ggtitle("Ceftriaxone") +
  labs(color='Error type')

combined_plot<-margin_plot/mic_ceftriaxone
ceftriaxone_combined_plot <- combined_plot + plot_layout(heights = c(1, 3))

########cefuroxime##########
cefuroxime_r_genes<-filter(cefuroxime_r_genes,p.value <0.05)

amrfinder<-read_tsv('amrfinder.tsv') 
tem_promotor_mutations<-filter(amrfinder,grepl('TEMp',`Gene symbol`))
ec_match2<-filter(ec_match2,guuid %in% amrfinder$Name)
amrfinder<-amrfinder%>% filter(Class=='BETA-LACTAM' ) %>%  select(guuid=Name,gene=`Gene symbol`)
length(unique(amrfinder$gene))
amrfinder$gene<-str_replace_all(amrfinder$gene,'-','_')

amrfinder$gene<-ifelse(grepl('blaTEM_',amrfinder$gene) & amrfinder$guuid %in% tem_promotor_mutations$Name,
                       paste0(amrfinder$gene,'_',tem_promotor_mutations$`Gene symbol`[match(amrfinder$guuid,tem_promotor_mutations$Name)]),amrfinder$gene)
amrfinder <- amrfinder %>%
  filter(!grepl("(?<!\\w)blaTEMp_", gene, perl = TRUE))
coamox<-read_tsv('amrfinder.tsv')  %>% filter(Name %in% ec_match2$guuid) %>% filter(grepl('BETA-LACTAM',Class )) %>%  select(guuid=Name,gene=`Gene symbol`) 
coamox$gene<-str_replace_all(coamox$gene,'-','_')

count<-amrfinder %>% filter(str_length(guuid) ==36) %>%  group_by(gene) %>% count() %>% filter(n>=10)

deal_with_rare_genes <- function(gene_name) {
  if (startsWith(gene_name, "ampC")) {
    gsub("^ampC_.*", "ampC_other", gene_name)
  } else if (!grepl("_", gene_name)) {
    paste0(gene_name, "_other")
  } else {
    gsub("_[^_]+$", "_other", gene_name)
  }
}

amrfinder$gene <- ifelse(amrfinder$gene %in% count$gene,
                         amrfinder$gene,
                         sapply(amrfinder$gene, deal_with_rare_genes))

singletons<-amrfinder %>%  group_by(gene) %>% count() %>% filter(n<5)

amrfinder$gene<-ifelse(amrfinder$gene %in% singletons$gene,"other gene",amrfinder$gene)
amrfinder$gene<-janitor::make_clean_names(amrfinder$gene)
amrfinder$gene<- str_remove(amrfinder$gene, "_\\d+$")


amrfinder<-filter(amrfinder,gene %in% cefuroxime_r_genes$term)
cefuroxime_r_genes<-select(cefuroxime_r_genes,term,estimate)
cefuroxime_r_genes$estimate<-2^gent_r_genes$estimate
amrfinder<-left_join(amrfinder,cefuroxime_r_genes,by=c("gene"="term"))


pheno<-select(ec_match,guuid,Cefuroxime) %>% filter(!is.na(Cefuroxime)) 
coamox_amrfinder<-filter(amrfinder,guuid %in%pheno$guuid)

pheno$Cefuroxime<-as.numeric(factor(pheno$Cefuroxime,levels=c('<=2','4','8','>8')))
mic_bands<-unique(pheno$Cefuroxime)

out=NULL
for(level in length(mic_bands):1){
  pheno2<-filter(pheno,Cefuroxime==level)
  pheno2<-left_join(pheno2,coamox_amrfinder,by=c("guuid"="guuid"))
  pheno2$estimate<-ifelse(is.na(pheno2$estimate),0,pheno2$estimate)
  pheno2<-pheno2 %>% group_by(guuid) %>% summarise(estimate=sum(estimate))

  out<-rbind(out,data.frame(level,pheno2))
}


names(out)<-c('MIC','guuid','estimate')

out$MIC<-case_when(out$MIC ==1 ~ '<=2',
                   out$MIC ==2 ~ '4',
                   out$MIC ==3 ~ '8',
                   out$MIC ==4 ~ '>8')

out$MIC<-factor(out$MIC,levels=c('<=2','4','8','>8'))

breaks <- seq(from = 0, to = max(out$estimate) + 1, by = 1)
out$estimate_bin<-cut(out$estimate,breaks = breaks, right = FALSE, labels = paste(head(breaks, -1), head(breaks, -1) +1 , sep="-"))

out<-out %>% group_by(MIC) %>% count(estimate_bin)
first_colors <- colorRampPalette(c("green", "darkgreen"))(3)
rest_colors <- colorRampPalette(c( "#F7BDBC", "#FA0A0A"))(length(unique(levels(out$estimate_bin))) - 3)
color_palette <- c(first_colors, rest_colors)

line_position <- which(levels(out$MIC) == "8") + 0.5


sens<-c('<=2','4','8')

out$class<-ifelse(out$estimate_bin=='0-1'| out$estimate_bin=='1-2'|out$estimate_bin=='2-3','S','R')

out <- out %>%
  mutate(mid_point = sapply(strsplit(as.character(estimate_bin), "-"), function(x) {
    mean(as.numeric(x))
  }))


mic_cefuroxime<-ggplot(out) +
  aes(x=MIC,y=n,fill=mid_point) +
  geom_bar(stat='identity')  +
  geom_vline(aes(xintercept = line_position), linetype="dotted", color = "black") +
  theme_minimal() +
  labs(fill="Total estimated\nfold-change in MIC") +
  scale_fill_gradientn(colors=color_palette)



calculate_ratio <- function(out, sens, conf.level = 0.95) {
  ratios <- out %>%
    group_by(MIC) %>%
    summarise(
      R_count = sum(n[class == "R"]),
      S_count = sum(n[class == "S"]),
      total = R_count + S_count,
      ratio = ifelse(MIC %in% sens,
                     R_count / total,
                     S_count / total)
    )
  
  # Compute the confidence intervals
  conf_int <- mapply(function(x, n) {
    test <- prop.test(x, n, conf.level = conf.level)
    c(lower = test$conf.int[1], upper = test$conf.int[2])
  },
  x = ifelse(ratios$MIC %in% sens, ratios$R_count, ratios$S_count),
  n = ratios$total,
  SIMPLIFY = FALSE
  )
  
  # Add the confidence intervals to the ratios
  ratios$lower_bound <- sapply(conf_int, "[[", "lower")
  ratios$upper_bound <- sapply(conf_int, "[[", "upper")
  
  # Join with the original data
  ratios %>%
    select(-R_count, -S_count, -total) %>%
    left_join(out, by = "MIC") %>%
    select(MIC, ratio, lower_bound, upper_bound) %>%
    distinct()
}

out2<-calculate_ratio(out,sens) 
out2$error<-ifelse(out2$MIC %in% sens,'Erroneous identification\nas wild-type','Erroneous identification\nas non wild-type')


margin_plot <- ggplot(out2) +
  aes(x=MIC, y=ratio*100, ymin=lower_bound*100, ymax=upper_bound*100,color=error) +
  geom_point() +
  geom_errorbar(width=0) +
  theme_minimal() +
  ylab("% Error") +
  xlab("") + 
  ggtitle("Cefuroxime") +
  labs(color='Error type')

combined_plot<-margin_plot/mic_cefuroxime
cefuroxime_combined_plot <- combined_plot + plot_layout(heights = c(1, 3))

##########piptaz###########

piptaz_r_genes<-filter(piptaz_r_genes,p.value<0.05)
amrfinder<-read_tsv('amrfinder.tsv') 

tem_promotor_mutations<-filter(amrfinder,grepl('TEMp',`Gene symbol`))
ec_match2<-filter(ec_match2,guuid %in% amrfinder$Name)
amrfinder<-amrfinder%>% filter(Class=='BETA-LACTAM' ) %>%  select(guuid=Name,gene=`Gene symbol`)
length(unique(amrfinder$gene))
amrfinder$gene<-str_replace_all(amrfinder$gene,'-','_')

amrfinder$gene<-ifelse(grepl('blaTEM_',amrfinder$gene) & amrfinder$guuid %in% tem_promotor_mutations$Name,
                       paste0(amrfinder$gene,'_',tem_promotor_mutations$`Gene symbol`[match(amrfinder$guuid,tem_promotor_mutations$Name)]),amrfinder$gene)
amrfinder <- amrfinder %>%
  filter(!grepl("(?<!\\w)blaTEMp_", gene, perl = TRUE))
coamox<-read_tsv('amrfinder.tsv')  %>% filter(Name %in% ec_match2$guuid) %>% filter(grepl('BETA-LACTAM',Class )) %>%  select(guuid=Name,gene=`Gene symbol`) 
coamox$gene<-str_replace_all(coamox$gene,'-','_')

count<-amrfinder %>% filter(str_length(guuid) ==36) %>%  group_by(gene) %>% count() %>% filter(n>=10)


deal_with_rare_genes <- function(gene_name) {
  if (startsWith(gene_name, "ampC")) {
    gsub("^ampC_.*", "ampC_other", gene_name)
  } else if (!grepl("_", gene_name)) {
    paste0(gene_name, "_other")
  } else {
    gsub("_[^_]+$", "_other", gene_name)
  }
}

amrfinder$gene <- ifelse(amrfinder$gene %in% count$gene,
                         amrfinder$gene,
                         sapply(amrfinder$gene, deal_with_rare_genes))

singletons<-amrfinder %>%  group_by(gene) %>% count() %>% filter(n<5)

amrfinder$gene<-ifelse(amrfinder$gene %in% singletons$gene,"other gene",amrfinder$gene)
amrfinder$gene<-janitor::make_clean_names(amrfinder$gene)
amrfinder$gene<- str_remove(amrfinder$gene, "_\\d+$")


amrfinder<-filter(amrfinder,gene %in% piptaz_r_genes$term)
piptaz_r_genes<-select(piptaz_r_genes,term,estimate)
piptaz_r_genes$estimate<-2^piptaz_r_genes$estimate
amrfinder<-left_join(amrfinder,piptaz_r_genes,by=c("gene"="term"))

pheno<-select(ec_match,guuid,Piptaz) %>% filter(!is.na(Piptaz)) 
piptaz_amrfinder<-filter(amrfinder,guuid %in%pheno$guuid)

pheno$Piptaz<-as.numeric(factor(pheno$Piptaz,levels=c('<=4/4','8/4','16/4','32/4','64/4','>64/4')))
mic_bands<-unique(pheno$Piptaz)

out=NULL
for(level in length(mic_bands):1){
  pheno2<-filter(pheno,Piptaz==level)
  pheno2<-left_join(pheno2,piptaz_amrfinder,by=c("guuid"="guuid"))
  pheno2$estimate<-ifelse(is.na(pheno2$estimate),0,pheno2$estimate)
  pheno2<-pheno2 %>% group_by(guuid) %>% summarise(estimate=sum(estimate))
  
  out<-rbind(out,data.frame(level,pheno2))
}


names(out)<-c('MIC','guuid','estimate')


out$MIC<-case_when(out$MIC ==1 ~ '<=4/4',
                    out$MIC ==2 ~ '8/4',
                    out$MIC ==3 ~ '16/4',
                    out$MIC ==4 ~ '32/4',
                    out$MIC ==5 ~ '64/4',
                    out$MIC ==6 ~ '>64/4')

out$MIC<-factor(out$MIC,levels=c('<=4/4','8/4','16/4','32/4','64/4','>64/4'))


breaks <- seq(from = 0, to = max(out$estimate) + 1, by = 1)
out$estimate_bin<-cut(out$estimate,breaks = breaks, right = FALSE, labels = paste(head(breaks, -1), head(breaks, -1) +1 , sep="-"))

out<-out %>% group_by(MIC) %>% count(estimate_bin)
first_colors <- colorRampPalette(c("green", "darkgreen"))(6)
rest_colors <- colorRampPalette(c( "#F7BDBC", "#FA0A0A"))(length(unique(levels(out$estimate_bin))) - 6)
color_palette <- c(first_colors, rest_colors)


line_position <- which(levels(out$MIC) == "8/4") + 0.5


sens<-c('<=4/4','8/4')

out$class<-ifelse(out$estimate_bin=='0-1'| out$estimate_bin=='1-2' | out$estimate_bin =='2-3' |
                    out$estimate_bin=='3-4' | out$estimate_bin=='4-5'|out$estimate_bin=='5-6','S','R')

out <- out %>%
  mutate(mid_point = sapply(strsplit(as.character(estimate_bin), "-"), function(x) {
    mean(as.numeric(x))
  }))


mic_piptaz<-ggplot(out) +
  aes(x=MIC,y=n,fill=mid_point) +
  geom_bar(stat='identity')  +
  geom_vline(aes(xintercept = line_position), linetype="dotted", color = "black") +
  theme_minimal() +
  labs(fill="Total estimated\nfold-change in MIC") +
  scale_fill_gradientn(colors=color_palette)



calculate_ratio <- function(out, sens, conf.level = 0.95) {
  ratios <- out %>%
    group_by(MIC) %>%
    summarise(
      R_count = sum(n[class == "R"]),
      S_count = sum(n[class == "S"]),
      total = R_count + S_count,
      ratio = ifelse(MIC %in% sens,
                     R_count / total,
                     S_count / total)
    )
  
  # Compute the confidence intervals
  conf_int <- mapply(function(x, n) {
    test <- prop.test(x, n, conf.level = conf.level)
    c(lower = test$conf.int[1], upper = test$conf.int[2])
  },
  x = ifelse(ratios$MIC %in% sens, ratios$R_count, ratios$S_count),
  n = ratios$total,
  SIMPLIFY = FALSE
  )
  
  # Add the confidence intervals to the ratios
  ratios$lower_bound <- sapply(conf_int, "[[", "lower")
  ratios$upper_bound <- sapply(conf_int, "[[", "upper")
  
  # Join with the original data
  ratios %>%
    select(-R_count, -S_count, -total) %>%
    left_join(out, by = "MIC") %>%
    select(MIC, ratio, lower_bound, upper_bound) %>%
    distinct()
}

out2<-calculate_ratio(out,sens) 
out2$error<-ifelse(out2$MIC %in% sens,'Erroneous identification\nas wild-type','Erroneous identification\nas non wild-type')


margin_plot <- ggplot(out2) +
  aes(x=MIC, y=ratio*100, ymin=lower_bound*100, ymax=upper_bound*100,color=error) +
  geom_point() +
  geom_errorbar(width=0) +
  theme_minimal() +
  ylab("% Error") +
  xlab("") + 
  ggtitle("Piperacillin-Tazobactam") +
  labs(color='Error type')

combined_plot<-margin_plot/mic_piptaz
piptaz_combined_plot <- combined_plot + plot_layout(heights = c(1, 3))

#####amp#######


amp_r_genes<-filter(amp_r_genes,p.value<0.05)
amrfinder<-read_tsv('amrfinder.tsv') 

tem_promotor_mutations<-filter(amrfinder,grepl('TEMp',`Gene symbol`))
ec_match2<-filter(ec_match2,guuid %in% amrfinder$Name)
amrfinder<-amrfinder%>% filter(Class=='BETA-LACTAM' ) %>%  select(guuid=Name,gene=`Gene symbol`)
length(unique(amrfinder$gene))
amrfinder$gene<-str_replace_all(amrfinder$gene,'-','_')

amrfinder$gene<-ifelse(grepl('blaTEM_',amrfinder$gene) & amrfinder$guuid %in% tem_promotor_mutations$Name,
                       paste0(amrfinder$gene,'_',tem_promotor_mutations$`Gene symbol`[match(amrfinder$guuid,tem_promotor_mutations$Name)]),amrfinder$gene)
amrfinder <- amrfinder %>%
  filter(!grepl("(?<!\\w)blaTEMp_", gene, perl = TRUE))
coamox<-read_tsv('amrfinder.tsv')  %>% filter(Name %in% ec_match2$guuid) %>% filter(grepl('BETA-LACTAM',Class )) %>%  select(guuid=Name,gene=`Gene symbol`) 
coamox$gene<-str_replace_all(coamox$gene,'-','_')

count<-amrfinder %>% filter(str_length(guuid) ==36) %>%  group_by(gene) %>% count() %>% filter(n>=10)


deal_with_rare_genes <- function(gene_name) {
  if (startsWith(gene_name, "ampC")) {
    gsub("^ampC_.*", "ampC_other", gene_name)
  } else if (!grepl("_", gene_name)) {
    paste0(gene_name, "_other")
  } else {
    gsub("_[^_]+$", "_other", gene_name)
  }
}

amrfinder$gene <- ifelse(amrfinder$gene %in% count$gene,
                         amrfinder$gene,
                         sapply(amrfinder$gene, deal_with_rare_genes))

singletons<-amrfinder %>%  group_by(gene) %>% count() %>% filter(n<5)


amrfinder$gene<-ifelse(amrfinder$gene %in% singletons$gene,"other gene",amrfinder$gene)
amrfinder$gene<-janitor::make_clean_names(amrfinder$gene)
amrfinder$gene<- str_remove(amrfinder$gene, "_\\d+$")


amrfinder<-filter(amrfinder,gene %in% amp_r_genes$term)
amp_r_genes<-select(amp_r_genes,term,estimate)
amp_r_genes$estimate<-2^amp_r_genes$estimate
amrfinder<-left_join(amrfinder,amp_r_genes,by=c("gene"="term"))
pheno<-select(ec_match,guuid,Ampicillin) %>% filter(!is.na(Ampicillin)) 
coamox_amrfinder<-filter(amrfinder,guuid %in%pheno$guuid)

pheno$Ampicillin<-as.numeric(factor(pheno$Ampicillin,levels=c('<=2','4','8','>8')))
mic_bands<-unique(pheno$Ampicillin)


out=NULL
for(level in length(mic_bands):1){
  pheno2<-filter(pheno,Ampicillin==level)
  pheno2<-left_join(pheno2,coamox_amrfinder,by=c("guuid"="guuid"))
  pheno2$estimate<-ifelse(is.na(pheno2$estimate),0,pheno2$estimate)
  pheno2<-pheno2 %>% group_by(guuid) %>% summarise(estimate=sum(estimate))

  out<-rbind(out,data.frame(level,pheno2))
}


names(out)<-c('MIC','guuid','estimate')

out$MIC<-case_when(out$MIC ==1 ~ '<=2',
                    out$MIC ==2 ~ '4',
                    out$MIC ==3 ~ '8',
                    out$MIC ==4 ~ '>8')

out$MIC<-factor(out$MIC,levels=c('<=2','4','8','>8'))

breaks <- seq(from = 0, to = max(out$estimate) + 1, by = 1)
out$estimate_bin<-cut(out$estimate,breaks = breaks, right = FALSE, labels = paste(head(breaks, -1), head(breaks, -1) +1 , sep="-"))

out<-out %>% group_by(MIC) %>% count(estimate_bin)
first_colors <- colorRampPalette(c("green", "darkgreen"))(5)
rest_colors <- colorRampPalette(c( "#F7BDBC", "#FA0A0A"))(length(unique(levels(out$estimate_bin))) - 5)
color_palette <- c(first_colors, rest_colors)

line_position <- which(levels(out$MIC) == "8") + 0.5


sens<-c('<=2','4','8')

out$class<-ifelse(out$estimate_bin=='0-1'| out$estimate_bin=='1-2' | out$estimate_bin =='2-3'|
                  out$estimate_bin=='3-4'|out$estimate_bin=='4-5','S','R')

out <- out %>%
  mutate(mid_point = sapply(strsplit(as.character(estimate_bin), "-"), function(x) {
    mean(as.numeric(x))
  }))

mic_amp<-ggplot(out) +
  aes(x=MIC,y=n,fill=mid_point) +
  geom_bar(stat='identity')  +
  geom_vline(aes(xintercept = line_position), linetype="dotted", color = "black") +
  theme_minimal() +
  labs(fill="Total estimated\nfold-change in MIC") +
  scale_fill_gradientn(colors=color_palette)




calculate_ratio <- function(out, sens, conf.level = 0.95) {
  ratios <- out %>%
    group_by(MIC) %>%
    summarise(
      R_count = sum(n[class == "R"]),
      S_count = sum(n[class == "S"]),
      total = R_count + S_count,
      ratio = ifelse(MIC %in% sens,
                     R_count / total,
                     S_count / total)
    )
  
  # Compute the confidence intervals
  conf_int <- mapply(function(x, n) {
    test <- prop.test(x, n, conf.level = conf.level)
    c(lower = test$conf.int[1], upper = test$conf.int[2])
  },
  x = ifelse(ratios$MIC %in% sens, ratios$R_count, ratios$S_count),
  n = ratios$total,
  SIMPLIFY = FALSE
  )
  
  # Add the confidence intervals to the ratios
  ratios$lower_bound <- sapply(conf_int, "[[", "lower")
  ratios$upper_bound <- sapply(conf_int, "[[", "upper")
  
  # Join with the original data
  ratios %>%
    select(-R_count, -S_count, -total) %>%
    left_join(out, by = "MIC") %>%
    select(MIC, ratio, lower_bound, upper_bound) %>%
    distinct()
}

out2<-calculate_ratio(out,sens) 
out2$error<-ifelse(out2$MIC %in% sens,'Erroneous identification\nas wild-type','Erroneous identification\nas non wild-type')


margin_plot <- ggplot(out2) +
  aes(x=MIC, y=ratio*100, ymin=lower_bound*100, ymax=upper_bound*100,color=error) +
  geom_point() +
  geom_errorbar(width=0) +
  theme_minimal() +
  ylab("% Error") +
  xlab("") + 
  ggtitle("Ampicillin") +
  labs(color='Error type')

combined_plot<-margin_plot/mic_amp 
amp_combined_plot <- combined_plot + plot_layout(heights = c(1, 3))
#######cotrim##########

cotrim_r_genes<-filter(cotrim_r_genes,p.value<0.05)
amrfinder<-read_tsv('amrfinder.tsv')
amrfinder<-amrfinder%>% filter(grepl('TRIMETHOPRIM',Class) | grepl('SULFONAMIDE',Class)) %>%  select(guuid=Name,gene=`Gene symbol`)
amrfinder$gene<-str_replace_all(amrfinder$gene,'-','_')

count<-amrfinder %>% filter(str_length(guuid) ==36) %>%  group_by(gene) %>% count() %>% filter(n>=10)


deal_with_rare_genes <- function(gene_name) {
  if (startsWith(gene_name, "ampC")) {
    gsub("^ampC_.*", "ampC_other", gene_name)
  } else if (!grepl("_", gene_name)) {
    paste0(gene_name, "_other")
  } else {
    gsub("_[^_]+$", "_other", gene_name)
  }
}

amrfinder$gene <- ifelse(amrfinder$gene %in% count$gene,
                         amrfinder$gene,
                         sapply(amrfinder$gene, deal_with_rare_genes))

singletons<-amrfinder %>%  group_by(gene) %>% count() %>% filter(n<5)

amrfinder$gene<-ifelse(amrfinder$gene %in% singletons$gene,"other gene",amrfinder$gene)
amrfinder$gene<-janitor::make_clean_names(amrfinder$gene)
amrfinder$gene<- str_remove(amrfinder$gene, "_\\d+$")


amrfinder<-filter(amrfinder,gene %in% cotrim_r_genes$term)
cotrim_r_genes<-select(cotrim_r_genes,term,estimate)
  cotrim_r_genes$estimate<-2^cotrim_r_genes$estimate
amrfinder<-left_join(amrfinder,cotrim_r_genes,by=c("gene"="term"))
pheno<-select(ec_match,guuid,Cotrim) %>% filter(!is.na(Cotrim)) 
coamox_amrfinder<-filter(amrfinder,guuid %in%pheno$guuid)

pheno<-select(ec_match,guuid,Cotrim) %>% filter(!is.na(Cotrim))
coamox_amrfinder<-filter(amrfinder,guuid %in%pheno$guuid)

pheno$Cotrim<-as.numeric(factor(pheno$Cotrim,levels=c('<=1/19','2/38','4/76','>4/76')))
mic_bands<-unique(pheno$Cotrim)

out=NULL
for(level in length(mic_bands):1){
  pheno2<-filter(pheno,Cotrim==level)
  pheno2<-left_join(pheno2,coamox_amrfinder,by=c("guuid"="guuid"))
  pheno2$estimate<-ifelse(is.na(pheno2$estimate),0,pheno2$estimate)
  pheno2<-pheno2 %>% group_by(guuid) %>% summarise(estimate=sum(estimate))

  out<-rbind(out,data.frame(level,pheno2))
}


names(out)<-c('MIC','guuid','estimate')

out$MIC<-case_when(out$MIC ==1 ~ '<=1/19',
                    out$MIC ==2 ~ '2/38',
                    out$MIC ==3 ~ '4/76',
                    out$MIC ==4 ~ '>4/76')

out$MIC<-factor(out$MIC,levels=c('<=1/19','2/38','4/76','>4/76'))

breaks <- seq(from = 0, to = max(out$estimate) + 1, by = 1)
out$estimate_bin<-cut(out$estimate,breaks = breaks, right = FALSE, labels = paste(head(breaks, -1), head(breaks, -1) +1 , sep="-"))

out<-out %>% group_by(MIC) %>% count(estimate_bin)
first_colors <- colorRampPalette(c("green", "darkgreen"))(9)
rest_colors <- colorRampPalette(c( "#F7BDBC", "#FA0A0A"))(length(unique(levels(out$estimate_bin))) - 9)
color_palette <- c(first_colors, rest_colors)

line_position <- which(levels(out$MIC) == "4/76") + 0.5


sens<-c('<=1/19','2/38','4/76')

out$class<-ifelse(out$estimate_bin =='0-1'| out$estimate_bin == '1-2' | out$estimate_bin =='2-3'|
out$estimate_bin == '3-4'| out$estimate_bin == '4-5'| out$estimate_bin == '5-6' | 
out$estimate_bin == '6-7'| out$estimate_bin == '7-8'| out$estimate_bin == '8-9'| out$estimate_bin == '9-10'| 
out$estimate_bin == '10-11','S','R')


out <- out %>%
  mutate(mid_point = sapply(strsplit(as.character(estimate_bin), "-"), function(x) {
    mean(as.numeric(x))
  }))

mic_cotrim<-ggplot(out) +
  aes(x=MIC,y=n,fill=mid_point) +
  geom_bar(stat='identity')  +
  geom_vline(aes(xintercept = line_position), linetype="dotted", color = "black") +
  theme_minimal() +
  labs(fill="Total estimated\nfold-change in MIC") +
  scale_fill_gradientn(colors=color_palette)




calculate_ratio <- function(out, sens, conf.level = 0.95) {
  ratios <- out %>%
    group_by(MIC) %>%
    summarise(
      R_count = sum(n[class == "R"]),
      S_count = sum(n[class == "S"]),
      total = R_count + S_count,
      ratio = ifelse(MIC %in% sens,
                     R_count / total,
                     S_count / total)
    )
  
  # Compute the confidence intervals
  conf_int <- mapply(function(x, n) {
    test <- prop.test(x, n, conf.level = conf.level)
    c(lower = test$conf.int[1], upper = test$conf.int[2])
  },
  x = ifelse(ratios$MIC %in% sens, ratios$R_count, ratios$S_count),
  n = ratios$total,
  SIMPLIFY = FALSE
  )
  
  # Add the confidence intervals to the ratios
  ratios$lower_bound <- sapply(conf_int, "[[", "lower")
  ratios$upper_bound <- sapply(conf_int, "[[", "upper")
  
  # Join with the original data
  ratios %>%
    select(-R_count, -S_count, -total) %>%
    left_join(out, by = "MIC") %>%
    select(MIC, ratio, lower_bound, upper_bound) %>%
    distinct()
}

out2<-calculate_ratio(out,sens) 
out2$error<-ifelse(out2$MIC %in% sens,'Erroneous identification\nas wild-type','Erroneous identification\nas non wild-type')


margin_plot <- ggplot(out2) +
  aes(x=MIC, y=ratio*100, ymin=lower_bound*100, ymax=upper_bound*100,color=error) +
  geom_point() +
  geom_errorbar(width=0) +
  theme_minimal() +
  ylab("% Error") +
  xlab("") + 
  ggtitle("Co-trimoxazole") +
  labs(color='Error type')

combined_plot<-margin_plot/mic_cotrim 
cotrim_combined_plot <- combined_plot + plot_layout(heights = c(1, 3))


###########mic assembly#############

amp_combined_plot + co_amox_combined_plot #+ piptaz_combined_plot

###########cipro phylogeny############

tree<-ape::read.tree('./mic.tree')
x<-tree$tip.label
amrfinder<-readr::read_tsv('./amrfinder.tsv')
x<-x[x %in% amrfinder$Name]

pruned.tree<-ape::drop.tip(tree,tree$tip.label[-match(x, tree$tip.label)])
tree<-phytools::midpoint.root(pruned.tree)

tree<-ggtree::ggtree(tree)
amrfinder<-dplyr::filter(amrfinder,Name %in% x & grepl('QUINOLONE',Class))
amrfinder<-select(amrfinder,guuid=Name,gene=`Gene symbol`)
x<-data.frame(x)
a<-left_join(x,amrfinder,by=c("x"="guuid"))
names(a)<-c("guuid","gene")
amrfinder<-data.frame(table(a)) %>% spread(gene,Freq)
x<-amrfinder$guuid
amrfinder<-select(amrfinder,-guuid)
rownames(amrfinder)<-x
amrfinder[amrfinder>1]<-1
g<-tree +ylim(-300,5000)
g
ggtree::gheatmap(g,amrfinder,low = "white",high = "black",colnames_angle = 90,colnames_offset_y = -200,font.size = 2) +
   theme(legend.position = "none")


#######co-amox sub-breakpoint#######
ec_match2<-ec_match  %>%  select(guuid,Coamox_lower,Coamox_upper) %>% filter(!is.na(Coamox_upper))  %>% filter(Coamox_upper <=3)
amrfinder<-read_tsv('amrfinder.tsv') %>% filter(Scope=="core")
tem_promotor_mutations<-filter(amrfinder,grepl('TEMp',`Gene symbol`))
ec_match2<-filter(ec_match2,guuid %in% amrfinder$Name)
amrfinder<-amrfinder%>% filter(Class=='BETA-LACTAM' ) %>%  select(guuid=Name,gene=`Gene symbol`)
length(unique(amrfinder$gene))
amrfinder$gene<-str_replace_all(amrfinder$gene,'-','_')

amrfinder$gene<-ifelse(grepl('blaTEM_',amrfinder$gene) & amrfinder$guuid %in% tem_promotor_mutations$Name,
                       paste0(amrfinder$gene,'_',tem_promotor_mutations$`Gene symbol`[match(amrfinder$guuid,tem_promotor_mutations$Name)]),amrfinder$gene)
amrfinder <- amrfinder %>%
  filter(!grepl("(?<!\\w)blaTEMp_", gene, perl = TRUE))

count<-amrfinder %>% filter(str_length(guuid) ==36) %>%  group_by(gene) %>% count() %>% filter(n>=10)


deal_with_rare_genes <- function(gene_name) {
  if (startsWith(gene_name, "ampC")) {
    gsub("^ampC_.*", "ampC_other", gene_name)
  } else if (!grepl("_", gene_name)) {
    paste0(gene_name, "_other")
  } else {
    gsub("_[^_]+$", "_other", gene_name)
  }
}

amrfinder$gene <- ifelse(amrfinder$gene %in% count$gene,
                         amrfinder$gene,
                         sapply(amrfinder$gene, deal_with_rare_genes))

singletons<-amrfinder %>%  group_by(gene) %>% count() %>% filter(n<5)

amrfinder$gene<-ifelse(amrfinder$gene %in% singletons$gene,"other gene",amrfinder$gene)



x<-data.frame(ec_match2$guuid)
x<-left_join(x,amrfinder,by=c("ec_match2.guuid"="guuid"))
x<-select(x,guuid=ec_match2.guuid,gene)
amrfinder<-data.frame(table(x)) %>% spread(gene,Freq)
names(amrfinder)<-str_replace_all(names(amrfinder),'-','_')
amrfinder[,2:ncol(amrfinder)][amrfinder[,2:ncol(amrfinder)]>1]<-1

genes<-names(amrfinder)
genes<-genes[!genes =='guuid']

amrfinder<-left_join(amrfinder,ec_match2,by=c("guuid"="guuid"))

mlst<-read_tsv('./mlst.tsv',col_names = F) %>% select(guuid=X1,spec=X2,mlst=X3)
#mlst$present<-1

big<-c('131','95','69','73')
mlst$mlst<-ifelse(mlst$mlst %in% big, mlst$mlst, 'other')

phylogroups<-read_tsv('./phylogroups.tsv')
mlst<-left_join(mlst,phylogroups,by=c("guuid"))
mlst$mlst<-ifelse(mlst$mlst =='other',mlst$phylogroup,mlst$mlst)
other<-c('U','U/cryptic','cryptic','E','G')
mlst$mlst<-ifelse(mlst$mlst %in% other,'other',mlst$mlst)
mlst<-select(mlst,guuid,mlst)

amrfinder<-left_join(amrfinder,mlst,by=c("guuid"))

out=NULL

names(amrfinder)<-janitor::make_clean_names(names(amrfinder))
genes<-janitor::make_clean_names(genes)


genes<-append(genes,"mlst")
amrfinder$mlst<-factor(amrfinder$mlst,levels=c('B2','131','95','73','69','A','B1','C','D','F','other'))
(Y <- with(amrfinder, Surv(coamox_lower, coamox_upper, type = "interval2")))
x<-amrfinder$guuid
amrfinder<-select(amrfinder,-guuid,-coamox_lower,-coamox_lower)

a<-amrfinder[,grep("^V",colnames(amrfinder))]
for(gene in genes){
  print(which(genes ==gene))
  
  print(gene)
  
  b<-select(amrfinder,gene)
  
  #c<-cbind(b,a)
  c<-b
  t<-as.formula(paste("Y", paste(names(c), collapse = " + "), sep=' ~ '))
  
  
  
  m<-broom::tidy(survreg(t ,data=c,dist="gaussian"),conf.int=T)
  
  m<-if(gene=='mlst'){m[2:13,]}else{m[2,]}    
  out=rbind(out,data.frame(m))
}

tidy_uni <- out %>% 
  dplyr::select(term, estimate, conf.low, conf.high, p.value) %>%
  dplyr::mutate(estimate=2^estimate,conf.low=2^conf.low,conf.high=2^conf.high) %>% 
  dplyr::mutate_at(vars(estimate, conf.low, conf.high), round, 2) %>%
  dplyr::mutate(conf_interval = paste0("(", conf.low, "-", conf.high, ")"),
                estimate_and_conf_interval = paste(estimate, conf_interval),
                p_value_formatted = ifelse(p.value < 0.001, "<0.001", signif(p.value, 2))) %>%
  dplyr::select(-estimate, -conf.low, -conf.high, -conf_interval, -p.value) # Optional, remove the original columns

count$gene<-janitor::make_clean_names(count$gene)
tidy_uni<-filter(tidy_uni,term %in% count$gene| grepl('mlst',term) | term=='other_gene')

include<-out
amrfinder_save<-amrfinder
include<-filter(include,!gene %in% mlst$mlst)

amrfinder<-amrfinder[names(amrfinder) %in% include$term | names(amrfinder) %in% names(a) ]
t<-as.formula(paste("Y", paste(names(amrfinder), collapse = " + "), sep=' ~ '))
model<-survreg(t, data=amrfinder,dist="gaussian")
out2<-broom::tidy(survreg(t, data=amrfinder,dist="gaussian"),conf.int=T)
out2$which<-"multivariable no population"

amrfinder<-amrfinder_save
include<-out
amrfinder<-amrfinder[names(amrfinder) %in% include$term | names(amrfinder) %in% names(a) | names(amrfinder) =='mlst' ]
t<-as.formula(paste("Y", paste(names(amrfinder), collapse = " + "), sep=' ~ '))
model<-survreg(t, data=amrfinder,dist="gaussian")
out3<-broom::tidy(survreg(t, data=amrfinder,dist="gaussian"),conf.int=T)
out3$which<-"multivariable population"



out$which<-"univariable"
out<-rbind(out,out2,out3)
count$gene<-janitor::make_clean_names(count$gene)

out<-filter(out,term %in% count$gene| grepl('mlst',term) | term=='other_gene')

out<-filter(out,!grepl('blaEC',gene))
out$estimate<-as.numeric(out$estimate)
out$conf.high<-as.numeric(out$conf.high)
out$conf.low<-as.numeric(out$conf.low)

out$term<-ifelse(out$term =='other_gene','other gene',out$term)
ggplot(out) +
  aes(x=term,y=estimate,color=which,group=which) +
  geom_point(position = position_dodge(0.9)) + geom_errorbar(data=out,aes(ymin=conf.low,ymax=conf.high),position="dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90)) +
  xlab('Gene') + ylab('Effect on log2 MIC') + 
  labs(color="Model")# +
  #ggtitle("Co-amoxiclav sub-breakpoint") #+ geom_hline(yintercept =log2(8) - model$coefficients[1],linetype="dashed")


tidy_multi <- broom::tidy(model, conf.int = TRUE) %>% 
  dplyr::select(term, estimate, conf.low, conf.high, p.value) %>%
  dplyr::mutate(estimate=2^estimate,conf.low=2^conf.low,conf.high=2^conf.high) %>% 
  dplyr::mutate_at(vars(estimate, conf.low, conf.high), round, 2) %>%
  dplyr::mutate(conf_interval = paste0("(", conf.low, "-", conf.high, ")"),
                estimate_and_conf_interval = paste(estimate, conf_interval),
                p_value_formatted = ifelse(p.value < 0.001, "<0.001", signif(p.value, 2))) %>%
  dplyr::select(-estimate, -conf.low, -conf.high, -conf_interval, -p.value) # Optional, remove the original columns

tidy_multi<-filter(tidy_multi,term %in% count$gene| grepl('mlst',term) | term=='other_gene')
all<-left_join(tidy_uni,tidy_multi,by=c("term"))
names(all)<-c("Gene","Univariable","p","Multivariable","p")
kableExtra::kable(all)


full<-logLik(model)
intercept<-logLik(update(model, . ~ 1))

R2_cox_snell <- 1 - exp((-2 / n) * (full - intercept))

#############cefur sub breakpoint##############

ec_match2<-select(ec_match,guuid,cefuroxime_lower,cefuroxime_upper) %>% filter(!is.na(cefuroxime_upper)) %>% filter(cefuroxime_upper <=3)
amrfinder<-read_tsv('amrfinder.tsv') %>% filter(Scope=="core")
tem_promotor_mutations<-filter(amrfinder,grepl('TEMp',`Gene symbol`))
ec_match2<-filter(ec_match2,guuid %in% amrfinder$Name)
amrfinder<-amrfinder%>% filter(Class=='BETA-LACTAM' ) %>%  select(guuid=Name,gene=`Gene symbol`)
length(unique(amrfinder$gene))
amrfinder$gene<-str_replace_all(amrfinder$gene,'-','_')

amrfinder$gene<-ifelse(grepl('blaTEM_',amrfinder$gene) & amrfinder$guuid %in% tem_promotor_mutations$Name,
                       paste0(amrfinder$gene,'_',tem_promotor_mutations$`Gene symbol`[match(amrfinder$guuid,tem_promotor_mutations$Name)]),amrfinder$gene)
amrfinder <- amrfinder %>%
  filter(!grepl("(?<!\\w)blaTEMp_", gene, perl = TRUE))


count<-amrfinder %>% filter(str_length(guuid) ==36) %>%  group_by(gene) %>% count() %>% filter(n>=10)

deal_with_rare_genes <- function(gene_name) {
  if (startsWith(gene_name, "ampC")) {
    gsub("^ampC_.*", "ampC_other", gene_name)
  } else if (!grepl("_", gene_name)) {
    paste0(gene_name, "_other")
  } else {
    gsub("_[^_]+$", "_other", gene_name)
  }
}

amrfinder$gene <- ifelse(amrfinder$gene %in% count$gene,
                         amrfinder$gene,
                         sapply(amrfinder$gene, deal_with_rare_genes))

singletons<-amrfinder %>%  group_by(gene) %>% count() %>% filter(n<5)

amrfinder$gene<-ifelse(amrfinder$gene %in% singletons$gene,"other gene",amrfinder$gene)




common_genes<-filter(ceph,gene %in% count$gene)
length(unique(ceph$gene)) # 24 ceph genes
length(unique(common_genes$gene)) # 8 occur >=10
x<-data.frame(ec_match2$guuid)
x<-left_join(x,amrfinder,by=c("ec_match2.guuid"="guuid"))
x<-select(x,guuid=ec_match2.guuid,gene)
amrfinder<-data.frame(table(x)) %>% spread(gene,Freq)
names(amrfinder)<-str_replace_all(names(amrfinder),'-','_')
amrfinder[,2:ncol(amrfinder)][amrfinder[,2:ncol(amrfinder)]>1]<-1

genes<-names(amrfinder)
genes<-genes[!genes =='guuid']

amrfinder<-left_join(amrfinder,ec_match2,by=c("guuid"="guuid"))

mlst<-read_tsv('./mlst.tsv',col_names = F) %>% select(guuid=X1,spec=X2,mlst=X3)


big<-c('131','95','69','73')
mlst$mlst<-ifelse(mlst$mlst %in% big, mlst$mlst, 'other')
phylogroups<-read_tsv('./phylogroups.tsv')
mlst<-left_join(mlst,phylogroups,by=c("guuid"))
mlst$mlst<-ifelse(mlst$mlst =='other',mlst$phylogroup,mlst$mlst)
other<-c('U','U/cryptic','cryptic','E','G')
mlst$mlst<-ifelse(mlst$mlst %in% other,'other',mlst$mlst)

mlst<-select(mlst,guuid,mlst)

amrfinder<-left_join(amrfinder,mlst,by=c("guuid"))

out=NULL

names(amrfinder)<-janitor::make_clean_names(names(amrfinder))
genes<-janitor::make_clean_names(genes)


genes<-append(genes,"mlst")
amrfinder$mlst<-factor(amrfinder$mlst,levels=c('B2','131','95','73','69','A','B1','C','D','F','other'))
(Y <- with(amrfinder, Surv(cefuroxime_lower, cefuroxime_upper, type = "interval2")))
guuids<-amrfinder$guuid
amrfinder<-select(amrfinder,-guuid,-cefuroxime_lower,-cefuroxime_upper)

a<-amrfinder[,grep("^V",colnames(amrfinder))]
for(gene in genes){
  print(which(genes ==gene))
  
  print(gene)
  
  b<-select(amrfinder,gene)
  
  c<-b
  t<-as.formula(paste("Y", paste(names(c), collapse = " + "), sep=' ~ '))
  
  
  
  m<-broom::tidy(survreg(t ,data=c,dist="gaussian"),conf.int=T)
  
  m<-if(gene=='mlst'){m[2:13,]}else{m[2,]}    
  out=rbind(out,data.frame(m))
}
include<-out

tidy_uni <- out %>% 
  dplyr::select(term, estimate, conf.low, conf.high, p.value) %>%
  dplyr::mutate(estimate=2^estimate,conf.low=2^conf.low,conf.high=2^conf.high) %>% 
  dplyr::mutate_at(vars(estimate, conf.low, conf.high), round, 2) %>%
  dplyr::mutate(conf_interval = paste0("(", conf.low, "-", conf.high, ")"),
                estimate_and_conf_interval = paste(estimate, conf_interval),
                p_value_formatted = ifelse(p.value < 0.001, "<0.001", signif(p.value, 2))) %>%
  dplyr::select(-estimate, -conf.low, -conf.high, -conf_interval, -p.value) # Optional, remove the original columns

count$gene<-janitor::make_clean_names(count$gene)
tidy_uni<-filter(tidy_uni,term %in% count$gene| grepl('mlst',term) | term=='other_gene')


amrfinder_save<-amrfinder

include<-filter(include,!gene %in% mlst$mlst)

amrfinder<-amrfinder[names(amrfinder) %in% include$term | names(amrfinder) %in% names(a) ]
t<-as.formula(paste("Y", paste(names(amrfinder), collapse = " + "), sep=' ~ '))
model<-survreg(t, data=amrfinder,dist="gaussian")
out2<-broom::tidy(survreg(t, data=amrfinder,dist="gaussian"),conf.int=T)
out2$which<-"multivariable no population"

amrfinder<-amrfinder_save
include<-out
amrfinder<-amrfinder[names(amrfinder) %in% include$term | names(amrfinder) %in% names(a) | names(amrfinder) =='mlst' ]
t<-as.formula(paste("Y", paste(names(amrfinder), collapse = " + "), sep=' ~ '))
model<-survreg(t, data=amrfinder,dist="gaussian")
out3<-broom::tidy(survreg(t, data=amrfinder,dist="gaussian"),conf.int=T)
out3$which<-"multivariable population"



out$which<-"univariable"
out<-rbind(out,out2,out3)
count$gene<-janitor::make_clean_names(count$gene)

out<-filter(out,term %in% count$gene| grepl('mlst',term) | term=='other_gene')




out<-filter(out,!grepl('blaEC',gene))
out$estimate<-as.numeric(out$estimate)
out$conf.high<-as.numeric(out$conf.high)
out$conf.low<-as.numeric(out$conf.low)


tidy_multi <- broom::tidy(model, conf.int = TRUE) %>% 
  dplyr::select(term, estimate, conf.low, conf.high, p.value) %>%
  dplyr::mutate(estimate=2^estimate,conf.low=2^conf.low,conf.high=2^conf.high) %>% 
  dplyr::mutate_at(vars(estimate, conf.low, conf.high), round, 2) %>%
  dplyr::mutate(conf_interval = paste0("(", conf.low, "-", conf.high, ")"),
                estimate_and_conf_interval = paste(estimate, conf_interval),
                p_value_formatted = ifelse(p.value < 0.001, "<0.001", signif(p.value, 2))) %>%
  dplyr::select(-estimate, -conf.low, -conf.high, -conf_interval, -p.value) # Optional, remove the original columns

tidy_multi<-filter(tidy_multi,term %in% count$gene| grepl('mlst',term) | term=='other_gene')

all2<-left_join(tidy_uni,tidy_multi,by=c("term"))
names(all2)<-c("Gene","Univariable","p","Multivariable","p")

all_both<-rbind(all,all2)

x<-c("mlstB1"  ,     "mlstA"    ,   "mlstC"  ,      "mlstD" ,       "mlstE",       
     "mlstF",        "mlstG" ,       "mlstother") 
all_both$Gene<-ifelse(all_both$Gene %in% x, str_replace(all_both$Gene,'mlst','phylogroup '),ifelse(grepl('mlst',all_both$Gene),str_replace(all_both$Gene,'mlst','ST '),all_both$Gene))

all_both %>% kableExtra::kable()

full<-logLik(model)
intercept<-logLik(update(model, . ~ 1))

R2_cox_snell <- 1 - exp((-2 / n) * (full - intercept))
#############cephal sub breakpoint##############

ec_match2<-select(ec_match,guuid,Cephalexin_lower,Cephalexin_upper) %>% filter(!is.na(Cephalexin_upper)) %>% filter(Cephalexin_upper <=4)
amrfinder<-read_tsv('amrfinder.tsv') 
amrfinder$`Gene symbol`<-str_replace_all(amrfinder$`Gene symbol`,'-','_')
ec_match2<-filter(ec_match2,guuid %in% amrfinder$Name)
amrfinder<-amrfinder%>% filter(Class=='BETA-LACTAM' ) %>%  select(guuid=Name,gene=`Gene symbol`) %>% filter(guuid %in% ec_match2$guuid)
length(unique(amrfinder$gene))

ceph<-read_tsv('amrfinder.tsv')  %>% filter(Name %in% ec_match2$guuid) %>% filter(grepl('CEPHALOSPORIN',Subclass )) %>%  select(guuid=Name,gene=`Gene symbol`)
ceph$gene<-str_replace_all(ceph$gene,'-','_')
count<-amrfinder %>% filter(str_length(guuid) ==36) %>%  group_by(gene) %>% count() %>% filter(n>=10)

deal_with_rare_genes <- function(gene_name) {
  if (startsWith(gene_name, "ampC")) {
    gsub("^ampC_.*", "ampC_other", gene_name)
  } else if (!grepl("_", gene_name)) {
    paste0(gene_name, "_other")
  } else {
    gsub("_[^_]+$", "_other", gene_name)
  }
}

amrfinder$gene <- ifelse(amrfinder$gene %in% count$gene,
                         amrfinder$gene,
                         sapply(amrfinder$gene, deal_with_rare_genes))

singletons<-amrfinder %>%  group_by(gene) %>% count() %>% filter(n<5)

amrfinder$gene<-ifelse(amrfinder$gene %in% singletons$gene,"other gene",amrfinder$gene)




common_genes<-filter(ceph,gene %in% count$gene)
length(unique(ceph$gene)) # 24 ceph genes
length(unique(common_genes$gene)) # 8 occur >=10
x<-data.frame(ec_match2$guuid)
x<-left_join(x,amrfinder,by=c("ec_match2.guuid"="guuid"))
x<-select(x,guuid=ec_match2.guuid,gene)
amrfinder<-data.frame(table(x)) %>% spread(gene,Freq)
names(amrfinder)<-str_replace_all(names(amrfinder),'-','_')
amrfinder[,2:ncol(amrfinder)][amrfinder[,2:ncol(amrfinder)]>1]<-1

genes<-names(amrfinder)
genes<-genes[!genes =='guuid']

amrfinder<-left_join(amrfinder,ec_match2,by=c("guuid"="guuid"))

mlst<-read_tsv('./mlst.tsv',col_names = F) %>% select(guuid=X1,spec=X2,mlst=X3)

big<-c('131','95','69','73')
mlst$mlst<-ifelse(mlst$mlst %in% big, mlst$mlst, 'other')
phylogroups<-read_tsv('./phylogroups.tsv')
mlst<-left_join(mlst,phylogroups,by=c("guuid"))
mlst$mlst<-ifelse(mlst$mlst =='other',mlst$phylogroup,mlst$mlst)
other<-c('U','U/cryptic','cryptic','E','G')
mlst$mlst<-ifelse(mlst$mlst %in% other,'other',mlst$mlst)

mlst<-select(mlst,guuid,mlst)

amrfinder<-left_join(amrfinder,mlst,by=c("guuid"))

out=NULL

names(amrfinder)<-janitor::make_clean_names(names(amrfinder))
genes<-janitor::make_clean_names(genes)


genes<-append(genes,"mlst")
amrfinder$mlst<-factor(amrfinder$mlst,levels=c('B2','131','95','73','69','A','B1','C','D','F','other'))
(Y <- with(amrfinder, Surv(cephalexin_lower, cephalexin_upper, type = "interval2")))
guuids<-amrfinder$guuid
amrfinder<-select(amrfinder,-guuid,-cephalexin_lower,-cephalexin_upper)

a<-amrfinder[,grep("^V",colnames(amrfinder))]
for(gene in genes){
  print(which(genes ==gene))
  
  print(gene)
  
  b<-select(amrfinder,gene)
  
  #c<-cbind(b,a)
  c<-b
  t<-as.formula(paste("Y", paste(names(c), collapse = " + "), sep=' ~ '))
  
  
  
  m<-broom::tidy(survreg(t ,data=c,dist="gaussian"),conf.int=T)
  
  m<-if(gene=='mlst'){m[2:13,]}else{m[2,]}    
  out=rbind(out,data.frame(m))
}
include<-out

amrfinder_save<-amrfinder

include<-filter(include,!gene %in% mlst$mlst)

amrfinder<-amrfinder[names(amrfinder) %in% include$term | names(amrfinder) %in% names(a) ]
t<-as.formula(paste("Y", paste(names(amrfinder), collapse = " + "), sep=' ~ '))
model<-survreg(t, data=amrfinder,dist="gaussian")
out2<-broom::tidy(survreg(t, data=amrfinder,dist="gaussian"),conf.int=T)
out2$which<-"multivariable no population"

amrfinder<-amrfinder_save
include<-out
amrfinder<-amrfinder[names(amrfinder) %in% include$term | names(amrfinder) %in% names(a) | names(amrfinder) =='mlst' ]
t<-as.formula(paste("Y", paste(names(amrfinder), collapse = " + "), sep=' ~ '))
model<-survreg(t, data=amrfinder,dist="gaussian")
out3<-broom::tidy(survreg(t, data=amrfinder,dist="gaussian"),conf.int=T)
out3$which<-"multivariable population"



out$which<-"univariable"
out<-rbind(out,out2,out3)

count$gene<-janitor::make_clean_names(count$gene)

out<-filter(out,term %in% count$gene| grepl('mlst',term) | term=='other_gene')

ceph_genes<-read_tsv('amrfinder.tsv') %>% filter(grepl('CEPHALOSPORIN',Subclass))
ceph_genes<-unique(ceph_genes$`Gene symbol`)
ceph_genes<-janitor::make_clean_names(ceph_genes)


out<-filter(out,!grepl('blaEC',gene))
out$estimate<-as.numeric(out$estimate)
out$conf.high<-as.numeric(out$conf.high)
out$conf.low<-as.numeric(out$conf.low)

x<-c("mlstB1"  ,     "mlstA"    ,   "mlstC"  ,      "mlstD" ,       "mlstE",       
     "mlstF",        "mlstG" ,       "mlstother") 
out$term<-ifelse(out$term %in% x, str_replace(out$term,'mlst','phylogroup '),ifelse(grepl('mlst',out$term),str_replace(out$term,'mlst','ST '),out$term))
out$term<-factor(out$term, levels=c("amp_c_c_11t",  "amp_c_c_42t" , "amp_c_t_32a" , "bla_ctx_m_14", "bla_ctx_m_15", "bla_ctx_m_27",
                                    "bla_ec" ,      "bla_ec_5" ,    "bla_oxa_1" ,   "bla_shv_1",    "bla_tem" ,     "bla_tem_1" ,  
                                    "bla_tem_30",   "bla_tem_40",   "other_gene",   "ST 131",      "ST 95"   ,    "ST 73",      
                                    "ST 69","phylogroup A","phylogroup B1","phylogroup B2","phylogroup C","phylogroup D","phylogroup E","phylogroup F","phylogroup G","phylogroup other"))
ggplot(out) +
  aes(x=term,y=estimate,color=which,group=which) +
  geom_point(position = position_dodge(0.9)) + geom_errorbar(data=out,aes(ymin=conf.low,ymax=conf.high),position="dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90)) +
  xlab('Gene') + ylab('Effect on log2 MIC') + 
  labs(color="Model") +
  ggtitle("Ceftriaxone") + geom_hline(yintercept =log2(8) - model$coefficients[1],linetype="dashed")



tidy_multi <- broom::tidy(model, conf.int = TRUE) %>% 
  dplyr::select(term, estimate, conf.low, conf.high, p.value) %>%
  dplyr::mutate_at(vars(estimate, conf.low, conf.high), round, 2) %>%
  dplyr::mutate(conf_interval = paste0("(", conf.low, "-", conf.high, ")"),
                estimate_and_conf_interval = paste(estimate, conf_interval),
                p_value_formatted = ifelse(p.value < 0.001, "<0.001", signif(p.value, 2))) %>%
  dplyr::select(-estimate, -conf.low, -conf.high, -conf_interval, -p.value) # Optional, remove the original columns


all2<-left_join(tidy_uni,tidy_multi,by=c("term"))
names(all2)<-c("Univariable","p","Multivariable","p")


all_both<-rbind(all,all2)
kableExtra::kable(all)




#######co-amox no known BL###########

amrfinder<-read_tsv('amrfinder.tsv') %>% filter(Class=='BETA-LACTAM') %>% filter(Scope=='core')
pheno<-select(ec_match,guuid,Coamox) %>% filter(!is.na(Coamox)) %>% filter(!guuid %in% amrfinder$Name) %>% filter(Coamox!='>8/2')

mlst<-read_tsv('./mlst.tsv',col_names = F) %>% select(guuid=X1,spec=X2,mlst=X3)

big<-c('131','95','69','73')
mlst$mlst<-ifelse(mlst$mlst %in% big, mlst$mlst, 'other')
phylogroups<-read_tsv('./phylogroups.tsv')
mlst<-left_join(mlst,phylogroups,by=c("guuid"))
mlst$mlst<-ifelse(mlst$mlst =='other',mlst$phylogroup,mlst$mlst)
other<-c('U','U/cryptic','cryptic','E','G')
mlst$mlst<-ifelse(mlst$mlst %in% other,'other',mlst$mlst)


mlst<-select(mlst,guuid,mlst)

pheno<-left_join(pheno,mlst,by=c("guuid"="guuid"))

pheno$Coamox<-factor(pheno$Coamox,levels=c('<=2/2','4/2','8/2','16/2','32/2','>32/2'))

my_cols<-c("#36FF39" ,"#28BD2D", "#176E1A", "#F5C8C7" ,"#FF817D" ,"#FC0307")

coamox_no_bl<-pheno %>% 
  filter(!is.na(mlst)) %>%
  group_by(mlst, Coamox) %>%
  summarise(n = n()) %>%
  mutate(prop = n / sum(n)) %>% 
ggplot(aes(x=mlst,y=prop,fill=Coamox)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  scale_fill_manual(values = my_cols) +
  labs(fill="Co-amoxiclav MIC") +
  xlab("ST/Phylogroup") +
  ylab("Proportion") +
  ggtitle("Co-amoxiclav")

########Cefuroxime no known BL###########


amrfinder<-read_tsv('./amrfinder.tsv') %>% filter(Class=='BETA-LACTAM') %>% filter(Scope=='core')
pheno<-select(ec_match,guuid,Cefuroxime) %>% filter(!is.na(Cefuroxime)) %>% filter(!guuid %in% amrfinder$Name) 

mlst<-read_tsv('./mlst.tsv',col_names = F) %>% select(guuid=X1,spec=X2,mlst=X3)



big<-c('131','95','69','73')
mlst$mlst<-ifelse(mlst$mlst %in% big, mlst$mlst, 'other')
phylogroups<-read_tsv('./phylogroups.tsv')
mlst<-left_join(mlst,phylogroups,by=c("guuid"))
mlst$mlst<-ifelse(mlst$mlst =='other',mlst$phylogroup,mlst$mlst)
other<-c('U','U/cryptic','cryptic','E','G')
mlst$mlst<-ifelse(mlst$mlst %in% other,'other',mlst$mlst)


mlst<-select(mlst,guuid,mlst)

pheno<-left_join(pheno,mlst,by=c("guuid"="guuid"))

pheno$Cefuroxime<-factor(pheno$Cefuroxime,levels=c('<=2','4','8','>8'))

my_cols<-c("#36FF39" ,"#28BD2D", "#176E1A", "#FC0307")

cefur_no_bl<-pheno %>% 
  filter(!is.na(mlst)) %>%
  group_by(mlst, Cefuroxime) %>%
  summarise(n = n()) %>%
  mutate(prop = n / sum(n)) %>% 
  ggplot(aes(x=mlst,y=prop,fill=Cefuroxime)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  scale_fill_manual(values = my_cols) +
  labs(fill="Cefuroxime MIC") +
  xlab("ST/Phylogroup") +
  ylab("Proportion") + 
  ggtitle("Cefuroxime")


coamox_no_bl + cefur_no_bl + plot_layout(guides="collect")
##########coamox ecoli tem cov######


tem<-rename(tem,'guuid'=X32)

ec_match2<-select(ec_match,guuid,Coamox_lower,Coamox_upper) %>% filter(!is.na(Coamox_upper))
amrfinder<-read_tsv('amrfinder.tsv')
tem_promotor_mutations<-filter(amrfinder,grepl('TEMp',`Gene symbol`))
ec_match2<-filter(ec_match2,guuid %in% amrfinder$Name)
amrfinder<-amrfinder%>% filter(Class=='BETA-LACTAM' ) %>%  select(guuid=Name,gene=`Gene symbol`)
length(unique(amrfinder$gene))
amrfinder$gene<-str_replace_all(amrfinder$gene,'-','_')

amrfinder$gene<-ifelse(grepl('blaTEM_',amrfinder$gene) & amrfinder$guuid %in% tem_promotor_mutations$Name,
                       paste0(amrfinder$gene,'_',tem_promotor_mutations$`Gene symbol`[match(amrfinder$guuid,tem_promotor_mutations$Name)]),amrfinder$gene)
amrfinder <- amrfinder %>%
  filter(!grepl("(?<!\\w)blaTEMp_", gene, perl = TRUE))
coamox<-read_tsv('amrfinder.tsv')  %>% filter(Name %in% ec_match2$guuid) %>% filter(grepl('BETA-LACTAM',Class )) %>%  select(guuid=Name,gene=`Gene symbol`) 
coamox$gene<-str_replace_all(coamox$gene,'-','_')


count<-amrfinder %>% filter(str_length(guuid) ==36) %>%  group_by(gene) %>% count() %>% filter(n>=10)


deal_with_rare_genes <- function(gene_name) {
  if (startsWith(gene_name, "ampC")) {
    gsub("^ampC_.*", "ampC_other", gene_name)
  } else if (!grepl("_", gene_name)) {
    paste0(gene_name, "_other")
  } else {
    gsub("_[^_]+$", "_other", gene_name)
  }
}

amrfinder$gene <- ifelse(amrfinder$gene %in% count$gene,
                         amrfinder$gene,
                         sapply(amrfinder$gene, deal_with_rare_genes))

singletons<-amrfinder %>%  group_by(gene) %>% count() %>% filter(n<5)

amrfinder$gene<-ifelse(amrfinder$gene %in% singletons$gene,"other gene",amrfinder$gene)

amrfinder2<-filter(amrfinder,!grepl('blaEC',gene))
c<-amrfinder2 %>% group_by(guuid) %>% count() %>% filter(n==1)
bt<-filter(amrfinder,gene=='blaTEM_1')

common_coamox<-filter(coamox,gene %in% count$gene)
length(unique(coamox$gene)) # 49 amp genes
length(unique(common_coamox$gene)) # 14 occur >=10
x<-data.frame(ec_match2$guuid)
x<-left_join(x,amrfinder,by=c("ec_match2.guuid"="guuid"))
x<-select(x,guuid=ec_match2.guuid,gene)
amrfinder<-data.frame(table(x)) %>% spread(gene,Freq)
names(amrfinder)<-str_replace_all(names(amrfinder),'-','_')
amrfinder[,2:ncol(amrfinder)][amrfinder[,2:ncol(amrfinder)]>1]<-1

genes<-names(amrfinder)
genes<-genes[!genes =='guuid']

amrfinder<-left_join(amrfinder,ec_match2,by=c("guuid"="guuid"))

mlst<-read_tsv('./mlst.tsv',col_names = F) %>% select(guuid=X1,spec=X2,mlst=X3)


big<-c('131','95','69','73')
mlst$mlst<-ifelse(mlst$mlst %in% big, mlst$mlst, 'other')

phylogroups<-read_tsv('./phylogroups.tsv')
mlst<-left_join(mlst,phylogroups,by=c("guuid"))
mlst$mlst<-ifelse(mlst$mlst =='other',mlst$phylogroup,mlst$mlst)
other<-c('U','U/cryptic','cryptic','E','G')
mlst$mlst<-ifelse(mlst$mlst %in% other,'other',mlst$mlst)

mlst<-select(mlst,guuid,mlst)

amrfinder<-left_join(amrfinder,mlst,by=c("guuid"))

out=NULL

names(amrfinder)<-janitor::make_clean_names(names(amrfinder))
genes<-janitor::make_clean_names(genes)


genes<-append(genes,"mlst")
amrfinder$mlst<-factor(amrfinder$mlst,levels=c('B2','131','95','73','69','A','B1','C','D','F','other'))
(Y <- with(amrfinder, Surv(coamox_lower, coamox_upper, type = "interval2")))
x<-amrfinder$guuid

amrfinder<-left_join(amrfinder,tem,by=c("guuid"))
amrfinder$tem_cov_corrected<-ifelse(is.na(amrfinder$tem_cov_corrected),0,amrfinder$tem_cov_corrected)
amrfinder$tem_cov_corrected<-log2(amrfinder$tem_cov_corrected +1)

amrfinder<-select(amrfinder,-guuid,-coamox_lower,-coamox_lower)

a<-amrfinder[,grep("^V",colnames(amrfinder))]

genes<-c(genes,'tem_cov_corrected')
keep<-amrfinder
for(gene in genes){
  print(which(genes ==gene))
  
  print(gene)
  
  b<-select(amrfinder,gene)
  
  #c<-cbind(b,a)
  c<-b
  t<-as.formula(paste("Y", paste(names(c), collapse = " + "), sep=' ~ '))
  
  
  
  m<-broom::tidy(survreg(t ,data=c,dist="gaussian"),conf.int=T)
  
  m<-if(gene=='mlst'){m[2:13,]}else{m[2,]}    
  out=rbind(out,data.frame(m))
}
include<-out
out_save<-out

amrfinder_save<-amrfinder
include<-filter(include,!gene %in% mlst$mlst)

amrfinder<-amrfinder[names(amrfinder) %in% include$term | names(amrfinder) %in% names(a) ]
t<-as.formula(paste("Y", paste(names(amrfinder), collapse = " + "), sep=' ~ '))
model<-survreg(t, data=amrfinder,dist="gaussian")
out2<-broom::tidy(survreg(t, data=amrfinder,dist="gaussian"),conf.int=T)
out2$which<-"multivariable no population"

amrfinder<-amrfinder_save
include<-out
amrfinder<-amrfinder[names(amrfinder) %in% include$term | names(amrfinder) %in% names(a) | names(amrfinder) =='mlst' ]
t<-as.formula(paste("Y", paste(names(amrfinder), collapse = " + "), sep=' ~ '))
model<-survreg(t, data=amrfinder,dist="gaussian")
out3<-broom::tidy(survreg(t, data=amrfinder,dist="gaussian"),conf.int=T)
out3$which<-"multivariable population"

out$which<-"univariable"
out<-rbind(out,out2,out3)

count$gene<-janitor::make_clean_names(count$gene)

out<-filter(out,term %in% count$gene| grepl('mlst',term) | term=='other_gene' | term=="tem_cov_corrected")

coamox_genes<-read_tsv('amrfinder.tsv') %>% filter(grepl('BETA-LACTAM',Subclass))
coamox_genes<-unique(coamox_genes$`Gene symbol`)
coamox_genes<-janitor::make_clean_names(coamox_genes)

a<-ifelse(out$term %in% coamox_genes,'red','black')

out<-filter(out,!grepl('blaEC',gene))
out$estimate<-as.numeric(out$estimate)
out$conf.high<-as.numeric(out$conf.high)
out$conf.low<-as.numeric(out$conf.low)


coamox<-ggplot(out) +
  aes(x=term,y=estimate,color=which,group=which) +
  geom_point(position = position_dodge(0.9)) + geom_errorbar(data=out,aes(ymin=conf.low,ymax=conf.high),position="dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90,color=a)) +
  xlab('Gene') + ylab('Effect on log2 MIC') + 
  labs(color="Model") +
  ggtitle("Co-amoxiclav") + geom_hline(yintercept =log2(8) - model$coefficients[1],linetype="dashed")

amrfinder<-amrfinder_save
predictors <- out_save$term
predictors<-ifelse(grepl("mlst",predictors),"mlst",predictors)
predictors<-unique(predictors)
predictors<-predictors[!is.na(predictors)  & ! predictors =="(Intercept)" & ! predictors == "Log(scale)"]

out=NULL



for (excluded_gene in predictors) {
  print(excluded_gene)
  other_predictors <- setdiff(predictors, excluded_gene)
  
  for (current_gene in other_predictors) {
    interaction_term <- paste(excluded_gene, current_gene, sep="*")
    formula_str <- as.formula(paste( "Y ~", paste(setdiff(other_predictors,current_gene), collapse = "+"), "+", interaction_term))
    summary_model <- broom::tidy(survreg(formula_str, data = amrfinder, dist = "gaussian"))
    out <- rbind(out,data.frame(summary_model))
    
  }
}




replace_mlst_patterns <- function(texts) {
  texts <- gsub("mlst[A-Za-z0-9_]*:", "mlst:", texts)
  texts <- gsub(":[A-Za-z0-9_]*mlst[A-Za-z0-9_]*", ":mlst", texts)
  return(texts)
}



interaction_terms<-filter(out,p.value < 0.01) %>% filter(grepl(':',term))
interaction_terms$term<-replace_mlst_patterns(interaction_terms$term)
interaction_terms<-unique(interaction_terms$term)
interaction_terms

model<-survreg(Y ~ amp_c_c_11t + amp_c_c_42t + amp_c_other + amp_c_t_32a + bla_carb_other + 
                 bla_cmy_other + bla_ctx_m_14 + bla_ctx_m_15 + bla_ctx_m_27 + 
                 bla_ctx_m_other + bla_ec + bla_ec_5 + bla_oxa_1 + bla_oxa_other + 
                 bla_shv_1 + bla_shv_other + bla_tem +   bla_tem_1_bla_te_mp_c32t + 
                 bla_tem_1_bla_te_mp_g162t + bla_tem_30_bla_te_mp_other + 
                 bla_tem_33_bla_te_mp_other + bla_tem_40 + bla_tem_other + 
                 omp_c_other + other_gene + mlst + bla_tem_1+tem_cov_corrected,data=amrfinder,dist="gaussian")

elimination<-  step(model, direction="backward")

selected_vars <- names(coef(elimination))
selected_vars<-selected_vars[!selected_vars=="(Intercept)"]
selected_vars<-gsub('mlst.*','mlst',selected_vars)
selected_vars<-unique(selected_vars)

interaction_string <- paste(interaction_terms, collapse = " + ")
new_formula <- as.formula(paste("Y ~ ", paste(selected_vars, collapse = " + "), " + ", interaction_string))
model<-survreg(new_formula,data=amrfinder,dist="gaussian")
elimination<-step(model,direction="backward",k=3.8)

###leave one out cross validation


coamox_gene<-read_tsv('amrfinder.tsv') %>% filter(Class=='BETA-LACTAM') %>% filter(Scope=='core')


out<-NULL
amrfinder$guuid<-x
for (i in 1:nrow(amrfinder)) {
  print(i)
  model<-survreg( Y[-i,] ~  amp_c_c_11t + amp_c_c_42t + amp_c_other + 
                    amp_c_t_32a + bla_carb_other + bla_cmy_other + bla_ctx_m_14 + 
                    bla_ctx_m_15 + bla_ctx_m_27 + bla_ctx_m_other + bla_ec_5 + 
                    bla_oxa_1 + bla_oxa_other + bla_shv_1 + bla_shv_other + bla_tem + 
                    bla_tem_1_bla_te_mp_c32t + bla_tem_1_bla_te_mp_g162t + bla_tem_30_bla_te_mp_other + 
                    bla_tem_33_bla_te_mp_other + bla_tem_40 + bla_tem_other + 
                    omp_c_other + other_gene + mlst + bla_tem_1 + tem_cov_corrected + 
                    amp_c_t_32a:bla_ec + amp_c_t_32a:bla_tem_30_bla_te_mp_other + 
                    bla_cmy_other:bla_oxa_1 + bla_ctx_m_14:bla_tem_1 + bla_ctx_m_15:bla_oxa_1 + 
                    bla_oxa_1:bla_tem_1 + bla_oxa_1:bla_tem_1_bla_te_mp_c32t + 
                    bla_oxa_other:bla_tem + bla_oxa_other:bla_tem_1 + bla_tem:tem_cov_corrected + 
                    omp_c_other:bla_tem_1 + other_gene:bla_tem_1 + bla_tem_1:tem_cov_corrected + 
                    other_gene:mlst, data = amrfinder[-i,], 
                  dist = "gaussian")
  
  prediction <- predict(model, newdata = amrfinder[i, ])
  correct<-ifelse(prediction >= as.numeric(Y[i,1]) & prediction <= as.numeric(Y[i,2]),1,0)
  
  out<-rbind(out,data.frame(correct,prediction,as.numeric(Y[i,1]), as.numeric(Y[i,2]),amrfinder$guuid[i]))
}


out$amrfinder_prediction<-ifelse(out$amrfinder.guuid.i. %in% coamox_gene$Name,'R','S')
out$amrfinder_prediction<-as.factor(out$amrfinder_prediction)
out$binary<-ifelse(out$as.numeric.Y.i..2.. >3,'R','S')
out$binary_prediction<-ifelse(out$prediction>3,'R','S')
out$binary<-as.factor(out$binary)
out$binary_prediction<-as.factor(out$binary_prediction)

extend_margin <- function(out, n){
  out$margin <- NA # Initialize margin column
  for(i in 0:n){
    condition <- out$prediction >= (out$as.numeric.Y.i..1.. - i) & 
      out$prediction <= (out$as.numeric.Y.i..2.. + i)
    out$margin[is.na(out$margin) & condition] <- i
  }
  return(out)
}


out2<-extend_margin(out,10)
out2$margin<-ifelse(out2$prediction < out2$as.numeric.Y.i..1..,out2$margin *-1,out2$margin)

out2$which<- ifelse(out2$binary_prediction==out2$binary & out2$binary=='S', 'Correct S',
                    ifelse(out2$binary_prediction==out2$binary & out2$binary =='R','Correct R',
                           ifelse(out2$binary_prediction=='S' & out2$binary =='R','VM',
                                  ifelse(out2$binary_prediction =='R' & out2$binary =='S','M',NA))))

out$within_1<-ifelse(out$prediction >= (out$as.numeric.Y.i..1.. -1) & out$prediction <= (out$as.numeric.Y.i..2.. +1),1,0)
out$within_2<-ifelse(out$prediction >= (out$as.numeric.Y.i..1.. -2) & out$prediction <= (out$as.numeric.Y.i..2.. +2),1,0)


caret::confusionMatrix(out$binary_prediction,out$binary)
table(out$correct)
prop.test(1620,(1620+1250))
1620+1250


table(out$within_1)
prop.test(2600,(2600+270))
2600+270

cm<-caret::confusionMatrix(as.factor(out$binary_prediction),reference=as.factor(out$binary),positive="R")
tn<-cm$table[2,2]
fn<-cm$table[2,1]
tp<-cm$table[1,1]
fp<-cm$table[1,2]
#concordance
prop.test(x=(tp+tn),n=(tp+tn+fp+fn))
tp+tn
tp+tn+fp+fn
#sens
prop.test(x=(tp),n=(tp+fn))
tp
tp+fn
#spec
prop.test(x=(tn),n=(tn+fp))
tn
tn+fp

#Maj
prop.test(x=(fp),n=(tn+fp))
fp
tn+fp
#Vmj
prop.test(x=(fn),n=(tp+fn))
fn
tp + fn

#PPV
prop.test(x=(tp),n=(tp+fp))
tp
tp+fp

#NPV
prop.test(x=(tn),n=(tn+fn))
tn
tn+fn



colors<-c("#09820B" ,"#03FF25" ,"#FFA703", "#ED3A37")
names(colors)<-c("Correct R","Correct S","M","VM")
ggplot(out2) +
  aes(x=margin, fill=which) +
  geom_bar(stat="count") + theme_minimal() +
  scale_fill_manual(values=colors) + 
  scale_x_continuous(breaks = seq(min(out2$margin),max(out2$margin),1)) + 
  xlab("Residual") +
  labs(fill="") +
  ggtitle("Gentamicin")

###############piptaz tem cov#############

tem<-rename(tem,'guuid'=X32)

ec_match2<-select(ec_match,guuid,piptaz_lower,piptaz_upper) %>% filter(!is.na(piptaz_upper))
amrfinder<-read_tsv('amrfinder.tsv')
tem_promotor_mutations<-filter(amrfinder,grepl('TEMp',`Gene symbol`))
ec_match2<-filter(ec_match2,guuid %in% amrfinder$Name)
amrfinder<-amrfinder%>% filter(Class=='BETA-LACTAM' ) %>%  select(guuid=Name,gene=`Gene symbol`)
length(unique(amrfinder$gene))
amrfinder$gene<-str_replace_all(amrfinder$gene,'-','_')

amrfinder$gene<-ifelse(grepl('blaTEM_',amrfinder$gene) & amrfinder$guuid %in% tem_promotor_mutations$Name,
                       paste0(amrfinder$gene,'_',tem_promotor_mutations$`Gene symbol`[match(amrfinder$guuid,tem_promotor_mutations$Name)]),amrfinder$gene)
amrfinder <- amrfinder %>%
  filter(!grepl("(?<!\\w)blaTEMp_", gene, perl = TRUE))
coamox<-read_tsv('amrfinder.tsv')  %>% filter(Name %in% ec_match2$guuid) %>% filter(grepl('BETA-LACTAM',Class )) %>%  select(guuid=Name,gene=`Gene symbol`) 
coamox$gene<-str_replace_all(coamox$gene,'-','_')


count<-amrfinder %>% filter(str_length(guuid) ==36) %>%  group_by(gene) %>% count() %>% filter(n>=10)


deal_with_rare_genes <- function(gene_name) {
  if (startsWith(gene_name, "ampC")) {
    gsub("^ampC_.*", "ampC_other", gene_name)
  } else if (!grepl("_", gene_name)) {
    paste0(gene_name, "_other")
  } else {
    gsub("_[^_]+$", "_other", gene_name)
  }
}

amrfinder$gene <- ifelse(amrfinder$gene %in% count$gene,
                         amrfinder$gene,
                         sapply(amrfinder$gene, deal_with_rare_genes))

singletons<-amrfinder %>%  group_by(gene) %>% count() %>% filter(n<5)

amrfinder$gene<-ifelse(amrfinder$gene %in% singletons$gene,"other gene",amrfinder$gene)

amrfinder2<-filter(amrfinder,!grepl('blaEC',gene))
c<-amrfinder2 %>% group_by(guuid) %>% count() %>% filter(n==1)
bt<-filter(amrfinder,gene=='blaTEM_1')

common_coamox<-filter(coamox,gene %in% count$gene)
length(unique(coamox$gene)) # 49 amp genes
length(unique(common_coamox$gene)) # 14 occur >=10
x<-data.frame(ec_match2$guuid)
x<-left_join(x,amrfinder,by=c("ec_match2.guuid"="guuid"))
x<-select(x,guuid=ec_match2.guuid,gene)
amrfinder<-data.frame(table(x)) %>% spread(gene,Freq)
names(amrfinder)<-str_replace_all(names(amrfinder),'-','_')
amrfinder[,2:ncol(amrfinder)][amrfinder[,2:ncol(amrfinder)]>1]<-1

genes<-names(amrfinder)
genes<-genes[!genes =='guuid']

amrfinder<-left_join(amrfinder,ec_match2,by=c("guuid"="guuid"))

mlst<-read_tsv('./mlst.tsv',col_names = F) %>% select(guuid=X1,spec=X2,mlst=X3)

big<-c('131','95','69','73')
mlst$mlst<-ifelse(mlst$mlst %in% big, mlst$mlst, 'other')

phylogroups<-read_tsv('./phylogroups.tsv')
mlst<-left_join(mlst,phylogroups,by=c("guuid"))
mlst$mlst<-ifelse(mlst$mlst =='other',mlst$phylogroup,mlst$mlst)
other<-c('U','U/cryptic','cryptic','E','G')
mlst$mlst<-ifelse(mlst$mlst %in% other,'other',mlst$mlst)

mlst<-select(mlst,guuid,mlst)

amrfinder<-left_join(amrfinder,mlst,by=c("guuid"))

out=NULL

names(amrfinder)<-janitor::make_clean_names(names(amrfinder))
genes<-janitor::make_clean_names(genes)


genes<-append(genes,"mlst")
amrfinder$mlst<-factor(amrfinder$mlst,levels=c('B2','131','95','73','69','A','B1','C','D','F','other'))
(Y <- with(amrfinder, Surv(piptaz_lower, piptaz_upper, type = "interval2")))
x<-amrfinder$guuid

amrfinder<-left_join(amrfinder,tem,by=c("guuid"))
amrfinder$tem_cov_corrected<-ifelse(is.na(amrfinder$tem_cov_corrected),0,amrfinder$tem_cov_corrected)
amrfinder$tem_cov_corrected<-log2(amrfinder$tem_cov_corrected +1)

amrfinder<-select(amrfinder,-guuid,-piptaz_lower,-piptaz_upper)

a<-amrfinder[,grep("^V",colnames(amrfinder))]

genes<-c(genes,'tem_cov_corrected')
keep<-amrfinder
for(gene in genes){
  print(which(genes ==gene))
  
  print(gene)
  
  b<-select(amrfinder,gene)
  
  #c<-cbind(b,a)
  c<-b
  t<-as.formula(paste("Y", paste(names(c), collapse = " + "), sep=' ~ '))
  
  
  
  m<-broom::tidy(survreg(t ,data=c,dist="gaussian"),conf.int=T)
  
  m<-if(gene=='mlst'){m[2:13,]}else{m[2,]}    
  out=rbind(out,data.frame(m))
}
include<-out
out_save<-out

amrfinder_save<-amrfinder
include<-filter(include,!gene %in% mlst$mlst)

amrfinder<-amrfinder[names(amrfinder) %in% include$term | names(amrfinder) %in% names(a) ]
t<-as.formula(paste("Y", paste(names(amrfinder), collapse = " + "), sep=' ~ '))
model<-survreg(t, data=amrfinder,dist="gaussian")
out2<-broom::tidy(survreg(t, data=amrfinder,dist="gaussian"),conf.int=T)
out2$which<-"multivariable no population"

amrfinder<-amrfinder_save
include<-out
amrfinder<-amrfinder[names(amrfinder) %in% include$term | names(amrfinder) %in% names(a) | names(amrfinder) =='mlst' ]
t<-as.formula(paste("Y", paste(names(amrfinder), collapse = " + "), sep=' ~ '))
model<-survreg(t, data=amrfinder,dist="gaussian")
out3<-broom::tidy(survreg(t, data=amrfinder,dist="gaussian"),conf.int=T)
out3$which<-"multivariable population"

out$which<-"univariable"
out<-rbind(out,out2,out3)

count$gene<-janitor::make_clean_names(count$gene)

out<-filter(out,term %in% count$gene| grepl('mlst',term) | term=='other_gene' | term=="tem_cov_corrected")


coamox_genes<-read_tsv('amrfinder.tsv') %>% filter(grepl('BETA-LACTAM',Subclass))
coamox_genes<-unique(coamox_genes$`Gene symbol`)
coamox_genes<-janitor::make_clean_names(coamox_genes)

a<-ifelse(out$term %in% coamox_genes,'red','black')

out<-filter(out,!grepl('blaEC',gene))
out$estimate<-as.numeric(out$estimate)
out$conf.high<-as.numeric(out$conf.high)
out$conf.low<-as.numeric(out$conf.low)


piptaz_tem_copy<-ggplot(out) +
  aes(x=term,y=estimate,color=which,group=which) +
  geom_point(position = position_dodge(0.9)) + geom_errorbar(data=out,aes(ymin=conf.low,ymax=conf.high),position="dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90,color=a)) +
  xlab('Gene') + ylab('Effect on log2 MIC') + 
  labs(color="Model") +
  ggtitle("Co-amoxiclav") + geom_hline(yintercept =log2(8) - model$coefficients[1],linetype="dashed")

amrfinder<-amrfinder_save
predictors <- out_save$term
predictors<-ifelse(grepl("mlst",predictors),"mlst",predictors)
predictors<-unique(predictors)
predictors<-predictors[!is.na(predictors)  & ! predictors =="(Intercept)" & ! predictors == "Log(scale)"]

out=NULL



for (excluded_gene in predictors) {
  print(excluded_gene)
  other_predictors <- setdiff(predictors, excluded_gene)
  
  for (current_gene in other_predictors) {
    interaction_term <- paste(excluded_gene, current_gene, sep="*")
    
    formula_str <- as.formula(paste("Y ~", paste(setdiff(other_predictors,current_gene), collapse = "+"), "+", interaction_term))
    summary_model <- broom::tidy(survreg(formula_str, data = amrfinder, dist = "gaussian"))
    
    out <- rbind(out,data.frame(summary_model))
    
  }
}




replace_mlst_patterns <- function(texts) {
  texts <- gsub("mlst[A-Za-z0-9_]*:", "mlst:", texts)
  texts <- gsub(":[A-Za-z0-9_]*mlst[A-Za-z0-9_]*", ":mlst", texts)
  return(texts)
}



interaction_terms<-filter(out,p.value < 0.01) %>% filter(grepl(':',term))
interaction_terms$term<-replace_mlst_patterns(interaction_terms$term)
interaction_terms<-unique(interaction_terms$term)
interaction_terms

model<-survreg(Y ~ amp_c_c_11t + amp_c_c_42t + amp_c_other + amp_c_t_32a + bla_carb_other + 
                 bla_cmy_other + bla_ctx_m_14 + bla_ctx_m_15 + bla_ctx_m_27 + 
                 bla_ctx_m_other + bla_ec + bla_ec_5 + bla_oxa_1 + bla_oxa_other + 
                 bla_shv_1 + bla_shv_other + bla_tem + bla_tem_1 + bla_tem_1_bla_te_mp_c32t + 
                 bla_tem_1_bla_te_mp_g162t + bla_tem_30_bla_te_mp_other + 
                 bla_tem_33_bla_te_mp_other + bla_tem_40 + bla_tem_other + 
                 omp_c_other + other_gene + mlst + tem_cov_corrected,data=amrfinder,dist="gaussian")

elimination<-  step(model, direction="backward")

selected_vars <- names(coef(elimination))
selected_vars<-selected_vars[!selected_vars=="(Intercept)"]
selected_vars<-gsub('mlst.*','mlst',selected_vars)
selected_vars<-unique(selected_vars)

interaction_string <- paste(interaction_terms, collapse = " + ")
new_formula <- as.formula(paste("Y ~ ", paste(selected_vars, collapse = " + "), " + ", interaction_string))
model<-survreg(new_formula,data=amrfinder,dist="gaussian")
elimination<-step(model,direction="backward",k=3.8)

###leave one out cross validation


coamox_gene<-read_tsv('amrfinder.tsv') %>% filter(Class=='BETA-LACTAM') %>% filter(Scope=='core')


out<-NULL
amrfinder$guuid<-x
for (i in 1:nrow(amrfinder)) {
  print(i)
  model<-survreg( Y[-i,] ~  amp_c_c_42t + bla_ctx_m_15 + bla_ctx_m_other + 
                    bla_ec + bla_oxa_1 + bla_oxa_other + bla_shv_1 + bla_tem_1 + 
                    bla_tem_1_bla_te_mp_c32t + bla_tem_1_bla_te_mp_g162t + bla_tem_33_bla_te_mp_other + 
                    bla_tem_other + omp_c_other + other_gene + tem_cov_corrected + 
                    bla_oxa_1:bla_shv_1 + bla_oxa_other:bla_tem_1 + bla_tem_1:omp_c_other + 
                    bla_tem_1:other_gene + bla_tem_1_bla_te_mp_c32t:tem_cov_corrected + 
                    omp_c_other:tem_cov_corrected, data = amrfinder[-i,], 
                  dist = "gaussian")
  
  prediction <- predict(model, newdata = amrfinder[i, ])
  correct<-ifelse(prediction >= as.numeric(Y[i,1]) & prediction <= as.numeric(Y[i,2]),1,0)
  
  out<-rbind(out,data.frame(correct,prediction,as.numeric(Y[i,1]), as.numeric(Y[i,2]),amrfinder$guuid[i]))
}


out$amrfinder_prediction<-ifelse(out$amrfinder.guuid.i. %in% coamox_gene$Name,'R','S')
out$amrfinder_prediction<-as.factor(out$amrfinder_prediction)
#breakpoint >8 = log2(8) = 3
out$binary<-ifelse(out$as.numeric.Y.i..2.. >3,'R','S')
out$binary_prediction<-ifelse(out$prediction>3,'R','S')
out$binary<-as.factor(out$binary)
out$binary_prediction<-as.factor(out$binary_prediction)

extend_margin <- function(out, n){
  out$margin <- NA # Initialize margin column
  for(i in 0:n){
    condition <- out$prediction >= (out$as.numeric.Y.i..1.. - i) & 
      out$prediction <= (out$as.numeric.Y.i..2.. + i)
    out$margin[is.na(out$margin) & condition] <- i
  }
  return(out)
}


out2<-extend_margin(out,10)
out2$margin<-ifelse(out2$prediction < out2$as.numeric.Y.i..1..,out2$margin *-1,out2$margin)

out2$which<- ifelse(out2$binary_prediction==out2$binary & out2$binary=='S', 'Correct S',
                    ifelse(out2$binary_prediction==out2$binary & out2$binary =='R','Correct R',
                           ifelse(out2$binary_prediction=='S' & out2$binary =='R','VM',
                                  ifelse(out2$binary_prediction =='R' & out2$binary =='S','M',NA))))

out$within_1<-ifelse(out$prediction >= (out$as.numeric.Y.i..1.. -1) & out$prediction <= (out$as.numeric.Y.i..2.. +1),1,0)
out$within_2<-ifelse(out$prediction >= (out$as.numeric.Y.i..1.. -2) & out$prediction <= (out$as.numeric.Y.i..2.. +2),1,0)


caret::confusionMatrix(out$binary_prediction,out$binary)
table(out$correct)
prop.test(2546,(2546+322))
2546+322


table(out$within_1)
prop.test(2713,(2713+155))
2713+155

cm<-caret::confusionMatrix(as.factor(out$binary_prediction),reference=as.factor(out$binary),positive="R")
tn<-cm$table[2,2]
fn<-cm$table[2,1]
tp<-cm$table[1,1]
fp<-cm$table[1,2]
#concordance
prop.test(x=(tp+tn),n=(tp+tn+fp+fn))
tp+tn
tp+tn+fp+fn
#sens
prop.test(x=(tp),n=(tp+fn))
tp
tp+fn
#spec
prop.test(x=(tn),n=(tn+fp))
tn
tn+fp

#Maj
prop.test(x=(fp),n=(tn+fp))
fp
tn+fp
#Vmj
prop.test(x=(fn),n=(tp+fn))
fn
tp + fn

#PPV
prop.test(x=(tp),n=(tp+fp))
tp
tp+fp

#NPV
prop.test(x=(tn),n=(tn+fn))
tn
tn+fn



colors<-c("#09820B" ,"#03FF25" ,"#FFA703", "#ED3A37")
names(colors)<-c("Correct R","Correct S","M","VM")
ggplot(out2) +
  aes(x=margin, fill=which) +
  geom_bar(stat="count") + theme_minimal() +
  scale_fill_manual(values=colors) + 
  scale_x_continuous(breaks = seq(min(out2$margin),max(out2$margin),1)) + 
  xlab("Residual") +
  labs(fill="") +
  ggtitle("Gentamicin")


############co-amox no blaTEM ecoli################
ec_match2<-select(ec_match,guuid,Coamox_lower,Coamox_upper) %>% filter(!is.na(Coamox_upper))
amrfinder<-read_tsv('amrfinder.tsv')
ec_match2<-filter(ec_match2,guuid %in% amrfinder$Name)
amrfinder<-amrfinder%>% filter(Class=='BETA-LACTAM' ) %>% filter(Scope=='core')  %>%  select(guuid=Name,gene=`Gene symbol`)
length(unique(amrfinder$gene))
one_gene<-amrfinder %>% group_by(guuid) %>% count() %>% filter(n==1)
blaTEM<-filter(amrfinder,gene=='blaTEM-1')
exclude<-filter(one_gene,guuid %in% blaTEM$guuid)
ec_match2<-filter(ec_match2,!guuid %in% exclude$guuid)
amrfinder$gene<-str_replace_all(amrfinder$gene,'-','_')
coamox<-read_tsv('amrfinder.tsv')  %>% filter(Name %in% ec_match2$guuid) %>% filter(grepl('BETA-LACTAM',Class )) %>%  select(guuid=Name,gene=`Gene symbol`) 
coamox$gene<-str_replace_all(coamox$gene,'-','_')

count<-amrfinder %>% filter(str_length(guuid) ==36) %>%  group_by(gene) %>% count() %>% filter(n>=10)


deal_with_rare_genes <- function(gene_name) {
  if (startsWith(gene_name, "ampC")) {
    gsub("^ampC_.*", "ampC_other", gene_name)
  } else if (!grepl("_", gene_name)) {
    paste0(gene_name, "_other")
  } else {
    gsub("_[^_]+$", "_other", gene_name)
  }
}

amrfinder$gene <- ifelse(amrfinder$gene %in% count$gene,
                         amrfinder$gene,
                         sapply(amrfinder$gene, deal_with_rare_genes))

singletons<-amrfinder %>%  group_by(gene) %>% count() %>% filter(n<5)

amrfinder$gene<-ifelse(amrfinder$gene %in% singletons$gene,"other gene",amrfinder$gene)


common_coamox<-filter(coamox,gene %in% count$gene)
length(unique(coamox$gene)) # 49 amp genes
length(unique(common_coamox$gene)) # 14 occur >=10
x<-data.frame(ec_match2$guuid)
x<-left_join(x,amrfinder,by=c("ec_match2.guuid"="guuid"))
x<-select(x,guuid=ec_match2.guuid,gene)
amrfinder<-data.frame(table(x)) %>% spread(gene,Freq)
names(amrfinder)<-str_replace_all(names(amrfinder),'-','_')
amrfinder[,2:ncol(amrfinder)][amrfinder[,2:ncol(amrfinder)]>1]<-1

genes<-names(amrfinder)
genes<-genes[!genes =='guuid']

amrfinder<-left_join(amrfinder,ec_match2,by=c("guuid"="guuid"))

mlst<-read_tsv('./mlst.tsv',col_names = F) %>% select(guuid=X1,spec=X2,mlst=X3)


big<-c('131','95','69','73')
mlst$mlst<-ifelse(mlst$mlst %in% big, mlst$mlst, 'other')

phylogroups<-read_tsv('./phylogroups.tsv')
mlst<-left_join(mlst,phylogroups,by=c("guuid"))
mlst$mlst<-ifelse(mlst$mlst =='other',mlst$phylogroup,mlst$mlst)
other<-c('U','U/cryptic','cryptic','E','G')
mlst$mlst<-ifelse(mlst$mlst %in% other,'other',mlst$mlst)

mlst<-select(mlst,guuid,mlst)

amrfinder<-left_join(amrfinder,mlst,by=c("guuid"))

out=NULL

names(amrfinder)<-janitor::make_clean_names(names(amrfinder))
genes<-janitor::make_clean_names(genes)


genes<-append(genes,"mlst")
amrfinder$mlst<-factor(amrfinder$mlst,levels=c('B2','131','95','73','69','A','B1','C','D','F','other'))
(Y <- with(amrfinder, Surv(coamox_lower, coamox_upper, type = "interval2")))
x<-amrfinder$guuid
amrfinder<-select(amrfinder,-guuid,-coamox_lower,-coamox_lower)

a<-amrfinder[,grep("^V",colnames(amrfinder))]
for(gene in genes){
  print(which(genes ==gene))
  
  print(gene)
  
  b<-select(amrfinder,gene)
  
  #c<-cbind(b,a)
  c<-b
  t<-as.formula(paste("Y", paste(names(c), collapse = " + "), sep=' ~ '))
  
  
  
  m<-broom::tidy(survreg(t ,data=c,dist="gaussian"),conf.int=T)
  
  m<-if(gene=='mlst'){m[2:13,]}else{m[2,]}    
  out=rbind(out,data.frame(m))
}
include<-out

amrfinder_save<-amrfinder
include<-filter(include,!gene %in% mlst$mlst)

amrfinder<-amrfinder[names(amrfinder) %in% include$term | names(amrfinder) %in% names(a) ]
t<-as.formula(paste("Y", paste(names(amrfinder), collapse = " + "), sep=' ~ '))
model<-survreg(t, data=amrfinder,dist="gaussian")
out2<-broom::tidy(survreg(t, data=amrfinder,dist="gaussian"),conf.int=T)
out2$which<-"multivariable no population"

amrfinder<-amrfinder_save
include<-out
amrfinder<-amrfinder[names(amrfinder) %in% include$term | names(amrfinder) %in% names(a) | names(amrfinder) =='mlst' ]
t<-as.formula(paste("Y", paste(names(amrfinder), collapse = " + "), sep=' ~ '))
model<-survreg(t, data=amrfinder,dist="gaussian")
out3<-broom::tidy(survreg(t, data=amrfinder,dist="gaussian"),conf.int=T)
out3$which<-"multivariable population"

coamox_r_genes<-filter(out3,round(estimate,1) >0) %>% filter(!grepl('mlst',term))
coamox_line<-log2(8) - model$coefficients[1]

out$which<-"univariable"
out<-rbind(out,out2,out3)

count$gene<-janitor::make_clean_names(count$gene)

out<-filter(out,term %in% count$gene| grepl('mlst',term) | term=='other_gene')

coamox_genes<-read_tsv('amrfinder.tsv') %>% filter(grepl('BETA-LACTAM',Subclass))
coamox_genes<-unique(coamox_genes$`Gene symbol`)
coamox_genes<-janitor::make_clean_names(coamox_genes)

a<-ifelse(out$term %in% coamox_genes,'red','black')

out<-filter(out,!grepl('blaEC',gene))
out$estimate<-as.numeric(out$estimate)
out$conf.high<-as.numeric(out$conf.high)
out$conf.low<-as.numeric(out$conf.low)


coamox_nosingle_blatem<-ggplot(out) +
  aes(x=term,y=estimate,color=which,group=which) +
  geom_point(position = position_dodge(0.9)) + geom_errorbar(data=out,aes(ymin=conf.low,ymax=conf.high),position="dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90,color=a)) +
  xlab('Gene') + ylab('Effect on log2 MIC') + 
  labs(color="Model") +
  ggtitle("Co-amoxiclav") + geom_hline(yintercept =log2(8) - model$coefficients[1],linetype="dashed")



interaction_model<-broom::tidy(survreg(Y ~ (amp_c_c_11t + amp_c_c_42t + amp_c_other + amp_c_t_32a + bla_carb_other + 
                                              bla_cmy_other + bla_ctx_m_14 + bla_ctx_m_15 + bla_ctx_m_27 + 
                                              bla_ctx_m_other + bla_dha_other + bla_oxa_1 + bla_oxa_other + 
                                              bla_shv_1 + bla_shv_other + bla_tem + bla_tem_1 + bla_tem_30 + 
                                              bla_tem_40 + bla_tem_other + fts_i_other + other_gene + mlst)^2,data=amrfinder,dist="gaussian"))

interaction_terms<-filter(interaction_model,p.value < 0.01) %>% filter(grepl(':',term))
interaction_terms$term<-ifelse(grepl('mlst',interaction_terms$term),sub("mlst.*","mlst",interaction_terms$term),interaction_terms$term)
interaction_terms$term
interaction_terms<-unique(interaction_terms$term)

model<-survreg(Y ~ amp_c_c_11t + amp_c_c_42t + amp_c_other + amp_c_t_32a + bla_carb_other + 
                 bla_cmy_other + bla_ctx_m_14 + bla_ctx_m_15 + bla_ctx_m_27 + 
                 bla_ctx_m_other + bla_dha_other + bla_oxa_1 + bla_oxa_other + 
                 bla_shv_1 + bla_shv_other + bla_tem + bla_tem_1 + bla_tem_30 + 
                 bla_tem_40 + bla_tem_other + fts_i_other + other_gene + mlst,data=amrfinder,dist="gaussian")

elimination<-  step(model, direction="backward")

selected_vars <- names(coef(elimination))
selected_vars<-selected_vars[!selected_vars=="(Intercept)"]
selected_vars<-gsub('mlst.*','mlst',selected_vars)
selected_vars<-unique(selected_vars)

interaction_string <- paste(interaction_terms, collapse = " + ")
new_formula <- as.formula(paste("Y ~ ", paste(selected_vars, collapse = " + "), " + ", interaction_string))
model<-survreg(new_formula,data=amrfinder,dist="gaussian")
elimination<-step(model,direction="backward")

###leave one out cross validation


coamox_gene<-read_tsv('amrfinder.tsv') %>% filter(Class=='BETA-LACTAM') %>% filter(Scope=='core')


out<-NULL
amrfinder$guuid<-x
for (i in 1:nrow(amrfinder)) {
  print(i)
  model<-survreg( Y[-i,] ~ amp_c_c_11t + amp_c_c_42t + amp_c_other + 
                    amp_c_t_32a + bla_carb_other + bla_cmy_other + bla_ctx_m_14 + 
                    bla_ctx_m_27 + bla_ctx_m_other + bla_dha_other + bla_oxa_1 + 
                    bla_oxa_other + bla_shv_1 + bla_shv_other + bla_tem + bla_tem_1 + 
                    bla_tem_30 + bla_tem_40 + bla_tem_other + fts_i_other + other_gene + 
                    mlst + amp_c_t_32a:bla_tem_30 + bla_cmy_other:bla_oxa_1 + 
                    bla_oxa_1:bla_ctx_m_15 + bla_tem:bla_ctx_m_15 + mlst:bla_ctx_m_15 + 
                    bla_shv_other:bla_tem_1 + bla_tem:mlst + bla_tem_1:other_gene + 
                    bla_tem_1:mlst, data = amrfinder[-i,], 
                  dist = "gaussian")
  
  prediction <- predict(model, newdata = amrfinder[i, ])
  correct<-ifelse(prediction >= as.numeric(Y[i,1]) & prediction <= as.numeric(Y[i,2]),1,0)
  
  out<-rbind(out,data.frame(correct,prediction,as.numeric(Y[i,1]), as.numeric(Y[i,2]),amrfinder$guuid[i]))
}


out$amrfinder_prediction<-ifelse(out$amrfinder.guuid.i. %in% coamox_gene$Name,'R','S')
out$amrfinder_prediction<-as.factor(out$amrfinder_prediction)
#breakpoint >8 = log2(8) = 3
out$binary<-ifelse(out$as.numeric.Y.i..2.. >3,'R','S')
out$binary_prediction<-ifelse(out$prediction>3,'R','S')
out$binary<-as.factor(out$binary)
out$binary_prediction<-as.factor(out$binary_prediction)

extend_margin <- function(out, n){
  out$margin <- NA # Initialize margin column
  for(i in 0:n){
    condition <- out$prediction >= (out$as.numeric.Y.i..1.. - i) & 
      out$prediction <= (out$as.numeric.Y.i..2.. + i)
    out$margin[is.na(out$margin) & condition] <- i
  }
  return(out)
}


out2<-extend_margin(out,10)
out2$margin<-ifelse(out2$prediction < out2$as.numeric.Y.i..1..,out2$margin *-1,out2$margin)

out2$which<- ifelse(out2$binary_prediction==out2$binary & out2$binary=='S', 'Correct S',
                    ifelse(out2$binary_prediction==out2$binary & out2$binary =='R','Correct R',
                           ifelse(out2$binary_prediction=='S' & out2$binary =='R','VM',
                                  ifelse(out2$binary_prediction =='R' & out2$binary =='S','M',NA))))

out$within_1<-ifelse(out$prediction >= (out$as.numeric.Y.i..1.. -1) & out$prediction <= (out$as.numeric.Y.i..2.. +1),1,0)
out$within_2<-ifelse(out$prediction >= (out$as.numeric.Y.i..1.. -2) & out$prediction <= (out$as.numeric.Y.i..2.. +2),1,0)


caret::confusionMatrix(out$binary_prediction,out$binary)
table(out$correct)
prop.test(1466,(1466+1404))
1466+1404


table(out$within_1)
prop.test(2447,(2447+423))
2447+423

cm<-caret::confusionMatrix(as.factor(out$binary_prediction),reference=as.factor(out$binary),positive="R")
tn<-cm$table[2,2]
fn<-cm$table[2,1]
tp<-cm$table[1,1]
fp<-cm$table[1,2]
#concordance
prop.test(x=(tp+tn),n=(tp+tn+fp+fn))
tp+tn
tp+tn+fp+fn
#sens
prop.test(x=(tp),n=(tp+fn))
tp
tp+fn
#spec
prop.test(x=(tn),n=(tn+fp))
tn
tn+fp

#Maj
prop.test(x=(fp),n=(tn+fp))
fp
tn+fp
#Vmj
prop.test(x=(fn),n=(tp+fn))
fn
tp + fn

#PPV
prop.test(x=(tp),n=(tp+fp))
tp
tp+fp

#NPV
prop.test(x=(tn),n=(tn+fn))
tn
tn+fn

colors<-c("#09820B" ,"#03FF25" ,"#FFA703", "#ED3A37")
names(colors)<-c("Correct R","Correct S","M","VM")
coamox_residual<-ggplot(out2) +
  aes(x=margin, fill=which) +
  geom_bar(stat="count") + theme_minimal() +
  scale_fill_manual(values=colors) + 
  scale_x_continuous(breaks = seq(min(out2$margin),max(out2$margin),1)) + 
  xlab("Residual") +
  labs(fill="") +
  ggtitle("Co-amoxiclav")



#########line genes/R genes###########
amp_r_genes$confirmed<-ifelse(amp_r_genes$conf.low>0,1,0)
amp_r_genes$which<-"amp"
amp_line_genes$confirmed<-ifelse(amp_line_genes$conf.low > (log2(amp_line)),1,0)

coamox_r_genes$confirmed<-ifelse(coamox_r_genes$conf.low>0,1,0)
coamox_r_genes$which<-"coamox"
coamox_line_genes$confirmed<-ifelse(coamox_line_genes$conf.low > (log2(coamox_line)),1,0)

ceftriaxone_r_genes$confirmed<-ifelse(ceftriaxone_r_genes$conf.low>0,1,0)
ceftriaxone_r_genes$which<-"ceftriaxone"
ceftriaxone_line_genes$confirmed<-ifelse(ceftriaxone_line_genes$conf.low > (log2(ceftriaxone_line)),1,0)

cefuroxime_r_genes$confirmed<-ifelse(cefuroxime_r_genes$conf.low>0,1,0)
cefuroxime_r_genes$which<-"cefuroxime"
cefuroxime_line_genes$confirmed<-ifelse(cefuroxime_line_genes$conf.low > (log2(cefuroxime_line)),1,0)

cipro_r_genes$confirmed<-ifelse(cipro_r_genes$conf.low>0,1,0)
cipro_r_genes$which<-"cipro"
cipro_line_genes$confirmed<-ifelse(cipro_line_genes$conf.low > (log2(cipro_line)),1,0)

gent_r_genes$confirmed<-ifelse(gent_r_genes$conf.low>0,1,0)
gent_r_genes$which<-"gent"
gent_line_genes$confirmed<-ifelse(gent_line_genes$conf.low > (log2(gent_line)),1,0)

piptaz_r_genes$confirmed<-ifelse(piptaz_r_genes$conf.low>0,1,0)
piptaz_r_genes$which<-"piptaz"
piptaz_line_genes$confirmed<-ifelse(piptaz_line_genes$conf.low > (log2(piptaz_line)),1,0)

cotrim_r_genes$confirmed<-ifelse(cotrim_r_genes$conf.low>0,1,0)
cotrim_r_genes$which<-"cotrim"
cotrim_line_genes$confirmed<-ifelse(cotrim_line_genes$conf.low > (log2(cotrim_line)),1,0)

rgenes<-rbind(amp_r_genes,coamox_r_genes,ceftriaxone_r_genes,cefuroxime_r_genes,cipro_r_genes,gent_r_genes,
              piptaz_r_genes,cotrim_r_genes)

rgenes<- filter(rgenes,
                !term=='(Intercept)' &
                  !term=='Log(scale)' &
                  !grepl('other',term))# &
                  #!grepl('_mp_',term))
table(rgenes$confirmed)

confirmed_rgenes<-filter(rgenes,confirmed==1)

linegenes<-rbind(amp_line_genes,coamox_line_genes,ceftriaxone_line_genes,cefuroxime_line_genes,cipro_line_genes,gent_line_genes,
                 piptaz_line_genes,cotrim_line_genes)

linegenes<- filter(linegenes,
                !term=='(Intercept)' &
                  !term=='Log(scale)' &
                  !grepl('other',term)) #&
                  #!grepl('_mp_',term))
confirmed_linegenes<-filter(linegenes,confirmed==1)

amrfinder<-read_tsv('./amrfinder.tsv') 
amrfinder<-filter(amrfinder,grepl('BETA-LACTAM',Class) | grepl('QUINOLONE',Class) | grepl('AMINOGLYCOSIDE',Class) |
                    grepl('TRIMETHOPRIM',Class) | grepl('SULFONAMIDE',Class))
amrfinder_unique<-distinct(amrfinder,`Gene symbol`,.keep_all = T)
table(amrfinder_unique$`Element subtype`)
amrfinder<-filter(amrfinder,Name %in% ec_match$guuid)
count<-amrfinder %>% filter(str_length(Name) ==36) %>%  group_by(`Gene symbol`) %>% count() %>% filter(n>=10)
amrfinder_common<-filter(amrfinder_unique,`Gene symbol` %in% count$`Gene symbol`)
table(amrfinder_common$`Element subtype`)

amp_all_genes$which<-"amp"
coamox_all_genes$which<-"coamox"
ceftriaxone_all_genes$which<-"ceftriaxone"
cefuroxime_all_genes$which<-"cefuroxime"
cipro_all_genes$which<-"cipro"
gent_all_genes$which<-"gent"
piptaz_all_genes$which<-"piptaz"
cotrim_all_genes$which<-"cotrim"


allgenes<-rbind(amp_all_genes,coamox_all_genes,ceftriaxone_all_genes,cefuroxime_all_genes,cipro_all_genes,gent_all_genes,
                piptaz_all_genes,cotrim_all_genes)

allgenes<- filter(allgenes,
                !term=='(Intercept)' &
                  !term=='Log(scale)' &
                  !grepl('other',term)) #&
                  #!grepl('_mp_',term))

classify<-read_tsv('amrfinder.tsv') %>% filter(`Element type`=='AMR') %>%  select(gene=`Gene symbol`,type=`Element subtype`) %>% distinct()
classify$gene<-janitor::make_clean_names(classify$gene)
allgenes<-left_join(allgenes,classify,by=c("term"="gene"))
table(allgenes$type)

rgenes<-left_join(rgenes,classify,by=c("term"="gene"))
table(rgenes$type)

linegenes<-left_join(linegenes,classify,by=c("term"="gene"))
table(linegenes$type)

##########fig 1#################

amp_amrfinder<-read_tsv('amrfinder.tsv') %>% filter(grepl('BETA-LACTAM',Class)) %>% filter(Scope=="core") %>% 
  filter(!grepl('blaTEMp_',`Gene symbol`))
amp_amrfinder<-amp_amrfinder %>% 
  group_by(Name) %>% 
  count()

amp<-select(ec_match,guuid,Ampicillin)
amp<-left_join(amp,amp_amrfinder,by=c("guuid"="Name"))
amp$n<-ifelse(is.na(amp$n),0,amp$n)
amp<-filter(amp,!is.na(Ampicillin))

amp<- amp %>% 
  group_by(Ampicillin,n) %>% 
  summarise(count=n())

amp$Ampicillin<-factor(amp$Ampicillin,levels = c('<=2','4','8','>8'))
amp_plot<-ggplot(amp) +
  aes(x=Ampicillin,y=count,fill=n) +
  geom_bar(stat="identity") +
  theme_minimal() +
  scale_fill_viridis_c() +
  xlab("MIC") +
  labs(fill="Number of ARGs/\nmutations") + 
  geom_vline(xintercept = as.numeric(which(levels(amp$Ampicillin) == "8")) + 0.5, linetype = "dashed") +
  ggtitle("Ampicillin")

cip_amrfinder<-read_tsv('amrfinder.tsv') %>% filter(grepl('QUINOLONE',Class)) %>% filter(Scope=="core")
cip_amrfinder<-cip_amrfinder %>% 
  group_by(Name) %>% 
  count()

cip<-select(ec_match,guuid,Ciprofloxacin)
cip<-left_join(cip,cip_amrfinder,by=c("guuid"="Name"))
cip$n<-ifelse(is.na(cip$n),0,cip$n)
cip<-filter(cip,!is.na(Ciprofloxacin))

cip$Ciprofloxacin<-case_when(cip$Ciprofloxacin == '<=0.125' ~ '<=0.25',
                             cip$Ciprofloxacin == '<=0.25' ~ '<=0.25',
                             cip$Ciprofloxacin == '0.25' ~ '<=0.25',
                             cip$Ciprofloxacin == '0.5' ~ '0.5',
                             cip$Ciprofloxacin == '1' ~ '1',
                             cip$Ciprofloxacin == '>1' ~ '>1')
cip<- cip %>% 
  group_by(Ciprofloxacin,n) %>% 
  summarise(count=n())


cip$Ciprofloxacin<-factor(cip$Ciprofloxacin,levels = c('<=0.25','0.5','1','>1'))
cip_plot<-ggplot(cip) +
  aes(x=Ciprofloxacin,y=count,fill=n) +
  geom_bar(stat="identity") +
  theme_minimal() +
  scale_fill_viridis_c() +
  xlab("MIC") +
  labs(fill="Number of ARGs/\nmutations")+
  geom_vline(xintercept = as.numeric(which(levels(cip$Ciprofloxacin) == "0.5")) + 0.5, linetype = "dashed") +
  ggtitle("Ciprofloxacin")
  

piptaz_amrfinder<-read_tsv('amrfinder.tsv') %>% filter(grepl('BETA-LACTAM',Class)) %>% filter(Scope=="core") %>% 
  filter(!grepl("blaTEMp_",`Gene symbol`))
piptaz_amrfinder<-piptaz_amrfinder %>% 
  group_by(Name) %>% 
  count()

piptaz<-select(ec_match,guuid,Piptaz)
piptaz<-left_join(piptaz,piptaz_amrfinder,by=c("guuid"="Name"))
piptaz$n<-ifelse(is.na(piptaz$n),0,piptaz$n)
piptaz<-filter(piptaz,!is.na(Piptaz))

piptaz<- piptaz %>% 
  group_by(Piptaz,n) %>% 
  summarise(count=n())

piptaz$Piptaz<-factor(piptaz$Piptaz,levels = c('<=4/4','8/4','16/4','32/4','64/4','>64/4'))
piptaz_plot<-ggplot(piptaz) +
  aes(x=Piptaz,y=count,fill=n) +
  geom_bar(stat="identity") +
  theme_minimal() +
  scale_fill_viridis_c() +
  xlab("MIC") +
  labs(fill="Number of ARGs/\nmutations") +
  geom_vline(xintercept = as.numeric(which(levels(piptaz$Piptaz) == "8/4")) + 0.5, linetype = "dashed") +
  ggtitle("Piperacillin-Tazobactam")
  
gent_amrfinder<-read_tsv('amrfinder.tsv') %>% filter(grepl('GENTAMICIN',Subclass)) %>% filter(Scope=="core")
gent_amrfinder<-gent_amrfinder %>% 
  group_by(Name) %>% 
  count()

gent<-select(ec_match,guuid,Gentamicin)
gent<-left_join(gent,gent_amrfinder,by=c("guuid"="Name"))
gent$n<-ifelse(is.na(gent$n),0,gent$n)
gent<-filter(gent,!is.na(Gentamicin))

gent<- gent %>% 
  group_by(Gentamicin,n) %>% 
  summarise(count=n())

gent$Gentamicin<-factor(gent$Gentamicin,levels = c('<=1','2','4','>4'))
gent_plot<-ggplot(gent) +
  aes(x=Gentamicin,y=count,fill=n) +
  geom_bar(stat="identity") +
  theme_minimal() +
  scale_fill_viridis_c(breaks=seq(0,2,1)) +
  xlab("MIC") +
  labs(fill="Number of ARGs/\nmutations") +
  geom_vline(xintercept = as.numeric(which(levels(gent$Gentamicin) == "2")) + 0.5, linetype = "dashed") +
  ggtitle("Gentamicin") 


ceft_amrfinder<-read_tsv('amrfinder.tsv') %>% filter(grepl('CEPHALOSPORIN',Subclass)) %>% filter(Scope=="core")
ceft_amrfinder<-ceft_amrfinder %>% 
  group_by(Name) %>% 
  count()

ceft<-select(ec_match,guuid,Ceftriaxone)

ceft$Ceftriaxone<- case_when(ceft$Ceftriaxone == "<=0.5" ~ "<=1",
                             ceft$Ceftriaxone == "<=1" ~ "<=1",
                             ceft$Ceftriaxone == "1" ~ "<=1",
                             ceft$Ceftriaxone == "2" ~ "2",
                             ceft$Ceftriaxone == "4" ~ "4",
                             ceft$Ceftriaxone == ">4" ~ ">4",
                             TRUE ~ ceft$Ceftriaxone)

ceft<-left_join(ceft,ceft_amrfinder,by=c("guuid"="Name"))
ceft$n<-ifelse(is.na(ceft$n),0,ceft$n)
ceft<-filter(ceft,!is.na(Ceftriaxone))

ceft<- ceft %>% 
  group_by(Ceftriaxone,n) %>% 
  summarise(count=n())

ceft$Ceftriaxone<-factor(ceft$Ceftriaxone,levels = c('<=1','2','4','>4'))
ceft_plot<-ggplot(ceft) +
  aes(x=Ceftriaxone,y=count,fill=n) +
  geom_bar(stat="identity") +
  theme_minimal() +
  scale_fill_viridis_c() +
  xlab("MIC") +
  labs(fill="Number of ARGs/\nmutations") +
  geom_vline(xintercept = as.numeric(which(levels(ceft$Ceftriaxone) == "2")) + 0.5, linetype = "dashed") +
  ggtitle("Ceftriaxone")


cef_amrfinder<-read_tsv('amrfinder.tsv') %>% filter(grepl('CEPHALOSPORIN',Subclass)) %>% filter(Scope=="core")
cef_amrfinder<-cef_amrfinder %>% 
  group_by(Name) %>% 
  count()

cef<-select(ec_match,guuid,Cefuroxime)
cef<-left_join(cef,cef_amrfinder,by=c("guuid"="Name"))
cef$n<-ifelse(is.na(cef$n),0,cef$n)
cef<-filter(cef,!is.na(Cefuroxime))

cef<- cef %>% 
  group_by(Cefuroxime,n) %>% 
  summarise(count=n())

cef$Cefuroxime<-factor(cef$Cefuroxime,levels = c('<=2','4','8','>8'))
cef_plot<-ggplot(cef) +
  aes(x=Cefuroxime,y=count,fill=n) +
  geom_bar(stat="identity") +
  theme_minimal() +
  scale_fill_viridis_c() +
  xlab("MIC") +
  labs(fill="Number of ARGs/\nmutations") +
  geom_vline(xintercept = as.numeric(which(levels(cef$Cefuroxime) == "8")) + 0.5, linetype = "dashed") +
  ggtitle("Cefuroxime")


coamox_amrfinder<-read_tsv('amrfinder.tsv') %>% filter(grepl('BETA-LACTAM',Class)) %>% filter(Scope=="core") %>% 
  filter(!grepl("blaTEMp_",`Gene symbol`))
coamox_amrfinder<-coamox_amrfinder %>% 
  group_by(Name) %>% 
  count()

coamox<-select(ec_match,guuid,Coamox) %>% filter(!Coamox == '>8/2')
coamox<-left_join(coamox,coamox_amrfinder,by=c("guuid"="Name"))
coamox$n<-ifelse(is.na(coamox$n),0,coamox$n)
coamox<-filter(coamox,!is.na(Coamox))

coamox<- coamox %>% 
  group_by(Coamox,n) %>% 
  summarise(count=n())

coamox$Coamox<-factor(coamox$Coamox,levels = c('<=2/2','4/2','8/2','16/2','32/2','>32/2'))
coamox_plot<-ggplot(coamox) +
  aes(x=Coamox,y=count,fill=n) +
  geom_bar(stat="identity") +
  theme_minimal() +
  scale_fill_viridis_c() +
  xlab("MIC") +
  labs(fill="Number of ARGs/\nmutations") +
  geom_vline(xintercept = as.numeric(which(levels(coamox$Coamox) == "8/2")) + 0.5, linetype = "dashed") +
  ggtitle("Co-amoxiclav")


cotrim_amrfinder<-read_tsv('amrfinder.tsv') %>% filter(grepl('TRIMETHOPRIM',Class)|grepl('SULFONAMIDE',Class)) %>% filter(Scope=="core")
cotrim_amrfinder<-cotrim_amrfinder %>% 
  group_by(Name) %>% 
  count()

cotrim<-select(ec_match,guuid,Cotrim)
cotrim<-left_join(cotrim,cotrim_amrfinder,by=c("guuid"="Name"))
cotrim$n<-ifelse(is.na(cotrim$n),0,cotrim$n)
cotrim<-filter(cotrim,!is.na(Cotrim))

cotrim<- cotrim %>% 
  group_by(Cotrim,n) %>% 
  summarise(count=n())

cotrim$Cotrim<-factor(cotrim$Cotrim,levels = c('<=1/19','2/38','4/76','>4/76'))
cotrim_plot<-ggplot(cotrim) +
  aes(x=Cotrim,y=count,fill=n) +
  geom_bar(stat="identity") +
  theme_minimal() +
  scale_fill_viridis_c() +
  xlab("MIC") +
  labs(fill="Number of ARGs/\nmutations") +
  geom_vline(xintercept = as.numeric(which(levels(cotrim$Cotrim) == "4/76")) + 0.5, linetype = "dashed") +
  ggtitle("Co-trimoxazole")

fig1<-amp_plot + coamox_plot + piptaz_plot + cef_plot + ceft_plot + gent_plot + cip_plot + cotrim_plot 

#############fig sup dist##############

p<-coamox_sup_dist / piptaz_sup_dist / cefur_sup_dist



out3<-survreg(coamox_t, data=amrfinder,dist="gaussian")
test<-predict(out3,new_data=amrfinder)

  