# Calculate GDVIS parameters

This function calculates GDVIS parameters. Input is based on LDSC
estimates. h2 needs to be on the observed scale with 50:50
ascertainment. For 2D the input is: h2_sub1.con: the heritability of the
GWAS of the subtype1-cases versus controls. h2_se_sub1.con: the standard
error of the heritability name_sub1: the name of subtype1-cases N_sub1:
the number of individuals in the subtype1-cases group h2_sub2.con: the
heritability of the GWAS of the subtype2-cases versus controls.
h2_se_sub2.con: the standard error of the heritability name_sub2: the
name of subtype2-cases N_sub1: the number of individuals in the
subtype2-cases group name_allcases: the name of the pool of the
subtype1-cases and subtype2-cases group name_con: the name of the
controls rg_sub1.con_sub2.con: the genetic correlation between the GWAS
subtype1-cases versus controls and the GWAS subtype2-cases versus
controls rg_se_sub1.con_sub2.con: the standard error of the genetic
correlation folder_location: the folder that you want the data to be
stored at filename: the name of the file for saving pop.prev_case: the
population prevalence of the cases There are optional variables that you
can also add, this will allow GDVIS to check its calculated values
against LDSC values, the result of which can be found in the logfile.
However, this is not necessary and you can create the input list without
these variables. optional_LDSC_rg_allcases.con_sub1.con: the genetic
correlation between the GWAS all cases versus controls and the GWAS
subtype1-cases versus controls
optional_LDSC_rg_se_allcases.con_sub1.con: the standard error of the
genetic correlation optional_LDSC_rg_allcases.con_sub2.con: the genetic
correlation between the GWAS all cases versus controls and the GWAS
subtype2-cases versus controls
optional_LDSC_rg_se_allcases.con_sub2.con: the standard error of the
genetic correlation optional_LDSC_h2_sub1.sub2: the heritability of the
GWAS of subtype1-cases versus subtype2-cases
optional_LDSC_h2_se_sub1.sub2: the standard error of the heritability
optional_LDSC_rg_sub1.con_sub1.sub2: the genetic correlation of the GWAS
subtype1-cases versus controls and the GWAS subtype1-cases versus
subtype2-cases optional_LDSC_rg_se_sub1.con_sub1.sub2: the standard
error of the genetic correlation optional_LDSC_rg_sub2.con_sub1.sub2:
the genetic correlation of the GWAS subtype2-cases versus controls and
the GWAS subtype1-cases versus subtype2-cases
optional_LDSC_rg_se_sub2.con_sub1.sub2: the standard error of the
genetic correlation optional_LDSC_h2_allcases.con: the heritability of
the GWAS of allversus controls optional_LDSC_h2_se_allcases.con: the
standard error of the heritability

## Usage

``` r
GDVIS_calc(triangle.input.list, log_fun = message, webversion = FALSE)
```

## Arguments

- triangle.input.list:

  input list with all the parameters

- log_fun:

  needs to stay on message, important for webversion

- webversion:

  needs to stay on FALSE

## Details

If you want to plot a subtype with an external trait, you need this
additional input: plot_3D: this tells GDVIS that you want to calculate a
subtype with an external trait and should be set to TRUE h2_ext: the
heritabitlity of the external trait h2_se_ext: the standard error of the
heritabilty rg_sub1.con_ext: the genetic correlation of the GWAS
subtype1-cases versus controls and the GWAS of the external trait
rg_se_sub1.con_ext: the standard error of the genetic correlation
rg_sub2.con_ext: the genetic correlation of the GWAS subtype2-cases
versus controls and the GWAS of the external trait rg_se_sub2.con_ext:
the standard error of the genetic correlation name_ext: the name of the
external trait pop.prev_ext: the population prevalence of the external
trait, if the trait is continuous, set the prevalence to 0.5 There is
some optional input here as well, will allow GDVIS to check its
calculated values against LDSC values: optional_LDSC_rg_sub1.sub2_ext:
the genetic correlation of the GWAS subtype1-cases versus subtype2-cases
and the GWAS of the external trait optional_LDSC_rg_se_sub1.sub2_ext:
the standard error of the genetic correlation
optional_LDSC_rg_allcases.con_ext: the genetic correlation of the GWAS
all cases and the GWAS of the external trait
optional_LDSC_rg_se_allcases.con_ext: the standard error of the genetic
correlation If you want to compare two different subtypes, you need the
following: all input for both subtypes as described above, but staarting
with the triangle1./triangle2. e.g. triangle1.h2_sub1.con and
triangle2.h2_sub2.con plot_2D.2D: this tells GDVIS to run in the subtype
vs subtype mode rg_triangle1.sub1.sub2_triangle2.sub1.sub2: the genetic
correlation of the GWAS subtype1-cases versus subtype2-cases from
subtype definition A and the GWAS subtype1-cases versus subtype2-cases
from subtype definition B rg_se_triangle1.sub1.sub2_triangle2.sub1.sub2:
the standard error the the genetic correlation
rg_triangle1.sub1.con_triangle2.sub1.con: the genetic correlation of the
GWAS subtype1-cases versus controls from subtype definition A and the
GWAS subtype1-cases versus controls from subtype definition B
rg_se_triangle1.sub1.con_triangle2.sub1.con: the standard error the the
genetic correlation rg_triangle1.sub2.con_triangle2.sub1.con: the
genetic correlation of the GWAS subtype12-cases versus controls from
subtype definition A and the GWAS subtype1-cases versus controls from
subtype definition B rg_se_triangle1.sub2.con_triangle2.sub1.con: the
standard error the the genetic correlation
rg_triangle1.sub1.con_triangle2.sub2.con: the genetic correlation of the
GWAS subtype1-cases versus controls from subtype definition A and the
GWAS subtype2-cases versus controls from subtype definition B
rg_se_triangle1.sub1.con_triangle2.sub2.con: the standard error the the
genetic correlation rg_triangle1.sub2.con_triangle2.sub2.con: the
genetic correlation of the GWAS subtype2-cases versus controls from
subtype definition A and the GWAS subtype2-cases versus controls from
subtype definition B rg_se_triangle1.sub2.con_triangle2.sub2.con: the
standard error the the genetic correlation
