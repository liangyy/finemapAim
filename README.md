# run-aim
Run AIM [here](https://github.com/jzou1115/aim).

Here we wrap it up as an R function taking genotype, ASC data, eQTL results and returning the PIP and CS.
I mostly follow the [paper]((https://github.com/jzou1115/aim)) "Methods" section "Application to GTEx data" but instead of using CAVIAR as fine-mapping solver we deploy susieR since it is much faster when we include all SNPs in the cis-window (thousands of SNPs).
I use both of the meta-analysis scheme S^M1 and S^M2 as suggested by the paper. 
And we report the results of the scheme with smaller CS size on average.

