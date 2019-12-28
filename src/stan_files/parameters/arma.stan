vector[nt] phi0; 
// gqs() cant not deal with ? yet - as long as that's not fixed
// estimate phi and theta anyways
//matrix[meanstructure ? nt : 0, meanstructure ? nt : 0 ] phi;
//matrix[meanstructure ? nt : 0, meanstructure ? nt : 0 ] theta;
matrix<lower = -1, upper = 1>[nt,nt] phi;
matrix<lower = -1, upper = 1>[nt,nt] theta;
