      seqfile = Med15_phylogeny_Os_outgroup_b90.phy   * sequence data filename
     treefile = RAxML_bestTree.Med15_phylogeny_Os_outgroup_b90_H2b      * tree structure file name
      outfile = results.H2b.txt   * main result file name

        noisy = 1      * 0,1,2,3,9: how much rubbish on the screen
      verbose = 1      * 1:detailed output
      runmode = 0      * 0:user defined tree

      seqtype = 1      * 1:codons
    CodonFreq = 2      * 0:equal, 1:F1X4, 2:F3X4, 3:F61

        model = 2      * 0:one omega ratio for all branches
                       * 1:separate omega for each branch
                       * 2:user specified dN/dS ratios for branches

      NSsites = 0      * 

        icode = 0      * 0:universal code

    fix_kappa = 0      * 1:kappa fixed, 0:kappa to be estimated
        kappa = 2      * initial or fixed kappa

    fix_omega = 0      * 1:omega fixed, 0:omega to be estimated 
        omega = 0.2    * initial omega

    cleandata = 1  * remove sites with ambiguity data (1:yes, 0:no)?
