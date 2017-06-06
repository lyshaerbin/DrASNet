# DrASNet
Identify the mutations mediated alternative splicing perturbations in cancer

The DrASNet.R is the main function to identify the mutation-AS pairs in cancer

The DrASNet.R script will need the other three scirpts: Perturbed_AS_personalized.R, Mutation_AS.R and Trans_pri.R;

Perturbed_AS_personalized.R is used to identify the patient specific alternative splicing perturbations;

Mutation_AS.R is used to identify all the mutation-AS pairs in each patient;

Trans_pri.R is the greed method to identify the driver mutations.

Example input files included:

CHOL_mean.txt: The PSI profile for a specific cancer type;

CHOL_sample.txt: The sample information for a specific cancer type;

CHOL_sample_mut_code.txt: The mutations in each cancer patient;

network.txt: Two columns represent the interacting genes in protein-protein interaction network;

Output files:

Person_AS_event.txt: The perturbed alternative splicing events in each cancer patient;

Trans_mut_AS.txt: The trans-mutation-AS pairs;

Trans_gene_proi.txt: The trans driver genes;

Cis_mut_AS.txt: The cis-mutation-AS pairs;

Cis_gene_proi.txt: The cis driver genes;

DrAS_mut_AS_trans.txt: The trans mutation-alternative splicing events identified by DrAS-Net;

DrAS_mut_AS_cis.txt: The cis mutation-alternative splicing events identified by DrAS-Net;

If you have any question, please contact:

Nidhi Sahni: nsahni@mdanderson.org;

M. Madan Babu: madanm@mrc-lmb.cam.ac.uk

Song Yi: syi2@mdanderson.org

