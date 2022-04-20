# Changelog

## 0.1.0 (2022-04-20)


### Features

* add compatibility test against alignment module." ([7237bf0](https://www.github.com/hydra-genetics/snv_indels/commit/7237bf09705c7507ebc06fd149727f866746e04d))
* add conventional-prs workflow ([ceba918](https://www.github.com/hydra-genetics/snv_indels/commit/ceba918d775b9ffe8e81a03e2ff12108aa0e0b35))
* Add logging, refactor and fix some bugs ([8c5c09f](https://www.github.com/hydra-genetics/snv_indels/commit/8c5c09f8363d0cfc9af3e76e2b5daf379f8395d7))
* add pull-request template ([14250d1](https://www.github.com/hydra-genetics/snv_indels/commit/14250d1ac84ee7e2afc514f770082fce9bd87bbb))
* add resources to rules. ([d1a0efd](https://www.github.com/hydra-genetics/snv_indels/commit/d1a0efdd97ec02391def69066e7a4fb651900ea4))
* Add tabix indexing rule. Also adapted merg_vcf to take indexed vcf as input. ([a27cb03](https://www.github.com/hydra-genetics/snv_indels/commit/a27cb03ac1b864a8a9fc62d974f903ce1cbb0b68))
* added annotate rule ([601c698](https://www.github.com/hydra-genetics/snv_indels/commit/601c698829f1fd124244c0e1ec163c597428d4ef))
* Added decompose, normalize and sort of vcf files ([69a14b7](https://www.github.com/hydra-genetics/snv_indels/commit/69a14b7bb323bbbd857c4b9232861a9b16ad1219))
* added release-please workflow ([8e171e1](https://www.github.com/hydra-genetics/snv_indels/commit/8e171e1d87b2e600e7679381412acb3e8301d9eb))
* Added rule ensemble_vcf (bcbio recall). Also made sort_vcf more general. ([4d6464c](https://www.github.com/hydra-genetics/snv_indels/commit/4d6464c3a8fd205802a0a0ca1770358d318b3cbc))
* Match rules to standard, optimize, add functions ([c3185f9](https://www.github.com/hydra-genetics/snv_indels/commit/c3185f9f3b4adda639771c29373c6b62782d3ec4))
* moved rule: filter_vcf_on_format move to filtering module ([ebd1f41](https://www.github.com/hydra-genetics/snv_indels/commit/ebd1f41fd54e91e8e27a8759d560993d3f668c52))
* moved rule: filter_vcf_on_format move to filtering module ([3e3207e](https://www.github.com/hydra-genetics/snv_indels/commit/3e3207e4a8cc7793578cd646dd62216afe120f04))
* mutect2 rule ([202e7c2](https://www.github.com/hydra-genetics/snv_indels/commit/202e7c2303764da41c336d6e372a6705a456d663))
* new rule bed_split. Changed reference and bam files ([cf9dcd7](https://www.github.com/hydra-genetics/snv_indels/commit/cf9dcd72e5f69fa6609ef86038d077989e55172d))
* New rule bgzip for vcf-files ([3ab1abd](https://www.github.com/hydra-genetics/snv_indels/commit/3ab1abd14f4aa23e9596dc5c588379afe5dd68be))
* New rule merge vcf ([0f275ba](https://www.github.com/hydra-genetics/snv_indels/commit/0f275ba331cb975a2c44864533d53d238948ba02))
* New rule mutect2_gvcf that generates a gvcf instead of vcf ([0e1abcc](https://www.github.com/hydra-genetics/snv_indels/commit/0e1abccf58250e4267f0b792c707af5cc5bce6ed))
* New rule: fixAF ([bdd4efd](https://www.github.com/hydra-genetics/snv_indels/commit/bdd4efded783d0d5cdc758e774eb5eb480c4282d))
* New rule: merge_gvcf ([cf250cd](https://www.github.com/hydra-genetics/snv_indels/commit/cf250cda317010e19873302fa016679e0dbd405c))
* New vardict rule ([517129f](https://www.github.com/hydra-genetics/snv_indels/commit/517129f36059fdd78ea02b004abc8111d10f5e58))
* Remove snver, fix error message ([0ed6b16](https://www.github.com/hydra-genetics/snv_indels/commit/0ed6b16f59e082232887d88f0b315cb40be49be1))
* Rule for filtering of vcf files. Also added env for conda. ([1566924](https://www.github.com/hydra-genetics/snv_indels/commit/15669247e622438e3362011522fde4a70d149375))
* update version of common container ([52b7b6b](https://www.github.com/hydra-genetics/snv_indels/commit/52b7b6b86d5b60dd776d6486e7b975a6fdf6671b))
* Use barcode as part of units index and set wildcard constr ([d52c2cd](https://www.github.com/hydra-genetics/snv_indels/commit/d52c2cd4b62499d593d07a291571d3fa7a26d027))
* Use g.vcf extension for gvcf files ([eacc3a3](https://www.github.com/hydra-genetics/snv_indels/commit/eacc3a3606b846b6e91dd81b29ec81649511ad1c))


### Bug Fixes

* Added bam-index as input to vardict ([4c6d72b](https://www.github.com/hydra-genetics/snv_indels/commit/4c6d72b61100cabc73f8cfca523ef60c7665a131))
* change run column name to flowcell in units. ([5385f22](https://www.github.com/hydra-genetics/snv_indels/commit/5385f226d758de6c6fa773caa90b1164ef7162f7))
* Downgrade conda package versions to avoid crash ([061365a](https://www.github.com/hydra-genetics/snv_indels/commit/061365a9746dac670d4b085576950537bb1f93a2))
* fixed bam location based on alignment prep-release ([d08edb3](https://www.github.com/hydra-genetics/snv_indels/commit/d08edb374b8eeca48d0870cce212e28328bf69fb))
* lock singularity to version that works with our images from dockerhub ([5115f39](https://www.github.com/hydra-genetics/snv_indels/commit/5115f3905da86e48ba9fb0f760a0693328b38c27))
* Vardict skips reads with mapq 0. Also, correct columns in bed file as default ([72be218](https://www.github.com/hydra-genetics/snv_indels/commit/72be2184629e4d5266714540145a25c929dbb1d3))


### Documentation

* Add latest rule graph ([037d108](https://www.github.com/hydra-genetics/snv_indels/commit/037d108866d5224b30fbd97c0388666cfb626624))
* update pull-request template ([79f7cfb](https://www.github.com/hydra-genetics/snv_indels/commit/79f7cfb4bf2fff00d020b113eae8f717dee07061))
* Update README ([4a237f0](https://www.github.com/hydra-genetics/snv_indels/commit/4a237f0f2aee97995dfc782d608794e24ac494cc))
