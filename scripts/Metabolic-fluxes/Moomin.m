%%% change gene identifiers in model
 load('Recon2.v04.mat');
 dict=tdfread('ensembl_gene_id_2_model_gene_id.csv',',');
 dict.model_gene=num2str(dict.model_gene);
 dict.model_gene=cellstr(dict.model_gene);
 dict.model_gene=strtrim(dict.model_gene);
 ensembGenes=cellstr(dict.ensembl_gene_id);
 [~, Imod, Idict] = intersect(modelR204.genes, dict.model_gene);
 modelR204.genes(Imod)=ensembGenes(Idict);
 save('modelR204_fixed', 'modelR204');

%%% Run for each Series (example for Series 5):

load('modelR204_fixed.mat');
changeCobraSolver('ibm_cplex','milp');

data=tdfread('EBSeq_Series5.csv', ',');
inactiveGenes=tdfread('low_counts_Series5.csv', ',');

inactiveGenes=inactiveGenes.Gene;
inactiveRxns=findRxnsInActiveWithGenes(modelR204, inactiveGenes);
model=removeRxns(modelR204, inactiveRxns);

model5 = moomin(model, data, 'enumerate', 500, 'solverTimeLimit', 50000, 'stoichiometry', 0);
save('model5', 'model5');