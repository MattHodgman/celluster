#requires Rtools, Rcpp.
#input: data(z x a array, za), R from a QR of the
#first xx obs of data, a < ¼xx z (a x a array), a, z.
#output: R wrt QR of data (a x a array).
library(Rcpp)
cppFunction(’NumericMatrix xp(const DataFrame&MMG, NumericMatrix &IC, int &a, int &z, int &xx) {
    NumericMatrix GIJ ¼ IC;
    double u, v, x, ghi, ghj;
    NumericVector JC(a), JCT(z);

    for(int i ¼ xx; i < z; þþi) {

        for(int k ¼ 0; k < a; þþk) {
            JCT ¼ MMG[k];
            JC[k]¼JCT[i];
        }

        for(int j ¼ 0; j < a; þþj) {
            u ¼ GIJ(j,j);
            v ¼ JC[j];
            x ¼ sqrt(u*u þ v*v);
            u/¼x;
            v/¼x;
            
            for(int k ¼ j; k < a; þþk) {
                ghi ¼ JC[k];
                ghj ¼ GIJ(j,k);
                JC[k]¼-v*ghj þ u*ghi;
                GIJ(j,k)¼u*ghj þ v*ghi;
            }
        }
    }

    return(GIJ);
}’)

LDRseq <- function(dat, blockSize ¼ 100000) {
    n_c <- ncol(dat)
    n_r <- nrow(dat)
    n_blocks <- ceiling(n_r/blockSize)
    n_last <- n_r
    if(n_last ¼=0) {n_last <- blockSize}
    
    #create indices for each block
    ep <- (1:(n_blocks-1))*blockSize
    ep <- c(ep, ep[n_blocks-1]þn_last)
    sp <- (0:(n_blocks - 1)) * blockSize þ 1

    #initialize R
    R_ <- qr.R(qr(as.matrix(dat[sp[1]:ep[1],])))
    for (i in 2:n_blocks)   {
        block <- dat[sp[i]:ep[i]]
        nRB <- nrow(block)
        R_ <- xp(block, R_, n_c, nRB, 0)
    }
    R_
}