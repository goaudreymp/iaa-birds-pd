          seed = -1   * random if -1
       seqfile = fullpass_concat.fasta
      treefile = fullcalib_oliveros.tre
      mcmcfile = mcmc_pos1.txt
       outfile = out_mcmcpos1.txt

         ndata = 1    * number of partitions
       seqtype = 0    * 0: nucleotides; 1:codons; 2:AAs
       usedata = 2 ./in.BV   * 0: no data; 1:seq like; 2:normal approximation; 3:out.BV (in.BV)
         clock = 3    * 1: global clock; 2: independent rates; 3: correlated rates
       TipDate = 0 100 
       RootAge = 'B(0.535,0.665)'  * safe constraint on root age, used if no fossil for root

         model = 4    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
         alpha = 0.3124    * alpha for gamma rates at sites, from output of beast runs
         ncatG = 5    * No. categories in discrete gamma

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

       BDparas = 0.175 0.086 0.022   * birth, death, sampling
   kappa_gamma = 6 2      * gamma prior for kappa
   alpha_gamma = 1 1      * gamma prior for alpha

   rgene_gamma = 2 42 1    * gammaDir prior for rate for genes  2/42 = 0.0476 overall substitution rate
  sigma2_gamma = 1 10 1   * gammaDir prior for sigma^2     (for clock=2 or 3)

      finetune = 1: .1 .1 .1 .1 .1 .1 * auto (0 or 1): times, musigma2, rates, mixing, paras, FossilErr

         print = 1
        burnin = 80000
      sampfreq = 500
       nsample = 20000
