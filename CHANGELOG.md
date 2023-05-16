# Changelog

## [0.5.0](https://www.github.com/hydra-genetics/snv_indels/compare/v0.4.0...v0.5.0) (2023-05-16)


### Features

* update snakemake version, allow range up to version 8 ([fdb57fe](https://www.github.com/hydra-genetics/snv_indels/commit/fdb57fe4ba61dd6469502f8876fca0f8ef2b4087))


### Bug Fixes

* **deeptrio:** specify bam index files in the input for make examples ([eddf9f5](https://www.github.com/hydra-genetics/snv_indels/commit/eddf9f5a6a78bf2357aa1e89b8a215a8f90537e1))

## [0.4.0](https://www.github.com/hydra-genetics/snv_indels/compare/v0.3.0...v0.4.0) (2023-05-02)


### Features

* Add deeptrio ([6fa6952](https://www.github.com/hydra-genetics/snv_indels/commit/6fa6952625c46aa0148af68b2555303ade657a0a))
* add deepvariant ([ca464cf](https://www.github.com/hydra-genetics/snv_indels/commit/ca464cf18b794d0d306e36516bc3bef74964c4bc))
* make gatk_pass_filter configurable ([7f65825](https://www.github.com/hydra-genetics/snv_indels/commit/7f65825d9f7a61a4c3bdf6282ad8886762d951f0))
* remove conda env/testing ([bc72d83](https://www.github.com/hydra-genetics/snv_indels/commit/bc72d83fa5a77c120d641e396a98f8a50c6d2e26))


### Bug Fixes

* fix_af.py ([cd0a0d2](https://www.github.com/hydra-genetics/snv_indels/commit/cd0a0d2afae59d9abf4a09a6f9f431012e83c767))
* update config ([9145272](https://www.github.com/hydra-genetics/snv_indels/commit/9145272c760afef96e67de4749315d7a37b63177))
* Update workflow/schemas/config.schema.yaml ([4a6a939](https://www.github.com/hydra-genetics/snv_indels/commit/4a6a9395ac485b45d73326b8e4725da58036916d))


### Documentation

* update CODEOWNERS ([d7f8fa1](https://www.github.com/hydra-genetics/snv_indels/commit/d7f8fa1f3ecb181b3daed293e2551bf8f4a7ac95))

## [0.3.0](https://www.github.com/hydra-genetics/snv_indels/compare/v0.2.0...v0.3.0) (2022-11-09)


### Features

* make config.yaml location more flexible ([d0a5296](https://www.github.com/hydra-genetics/snv_indels/commit/d0a52962b071cad5450cdaf7e55ad4918ba9414d))
* make configfile/confgilefiles argument mandatory ([a860631](https://www.github.com/hydra-genetics/snv_indels/commit/a8606314bf3ea63af71c54ca3232a54c152faebf))
* update snakemake-version ([3d5ed33](https://www.github.com/hydra-genetics/snv_indels/commit/3d5ed33228032923f30bb2c399693f0a9aa2d04f))


### Bug Fixes

* added multiallelic variants to mutect2_pass_filter ([94d1a42](https://www.github.com/hydra-genetics/snv_indels/commit/94d1a4287a76e08d0022b9ebb29ea4196c38fd8c))
* make sure tabulate version is under 0.9.0 ([9bb1fbe](https://www.github.com/hydra-genetics/snv_indels/commit/9bb1fbe8aeb745998f894fa5e76393104e883892))
* set strict mode for conda ([47a7581](https://www.github.com/hydra-genetics/snv_indels/commit/47a7581df62ac26680ab57d92a21859e89b36e5d))


### Documentation

* update compatibility ([00805a9](https://www.github.com/hydra-genetics/snv_indels/commit/00805a96f891a4e628887c95d464202b5fa8fd2e))
* update readme information ([fffecea](https://www.github.com/hydra-genetics/snv_indels/commit/fffeceadff77936412df35ee03b12c350dc1ef29))

## [0.2.0](https://www.github.com/hydra-genetics/snv_indels/compare/v0.1.1...v0.2.0) (2022-06-01)


### Features

* added mutect2 filtering rule ([6a831f2](https://www.github.com/hydra-genetics/snv_indels/commit/6a831f2faa94dde69f41f0f6114634d580fda2d7))
* added mutect2 filtering schemas ([0e09f6d](https://www.github.com/hydra-genetics/snv_indels/commit/0e09f6dc94ec48bee5c67c7c9574c9befdddee7f))
* added rule for hard filtering of mutect2 ([983f889](https://www.github.com/hydra-genetics/snv_indels/commit/983f889dfaec1685f5d6516ac7d5a9063232e30c))
* update wrapper version ([669460f](https://www.github.com/hydra-genetics/snv_indels/commit/669460fa97f203e336307369b1ed1e1cda8c5380))


### Bug Fixes

* adapted tests after new rulenames ([6f51241](https://www.github.com/hydra-genetics/snv_indels/commit/6f512414ea74770c35d3d5c9e56f1211d047f479))
* added temp on output ([1b9c950](https://www.github.com/hydra-genetics/snv_indels/commit/1b9c950a37dbc09f5cd4d071a6d710f4e448addd))
* bugfixes ([25ad129](https://www.github.com/hydra-genetics/snv_indels/commit/25ad1294c4dfc044a69b4cfaffd74860705c0a60))
* bugfixes and testing ([fd71a2c](https://www.github.com/hydra-genetics/snv_indels/commit/fd71a2c7befa4943978b9de3b56b3f672c255384))
* correct singularity ([2867d9d](https://www.github.com/hydra-genetics/snv_indels/commit/2867d9d773a8c5a8a12d10056098d2e5e1f5d748))
* update compatibility config-file ([79d76e6](https://www.github.com/hydra-genetics/snv_indels/commit/79d76e6a9f4b9b6d93e011b5afc819bccf900a58))
* update config, mutect2 to gatk_mutect2 ([883b579](https://www.github.com/hydra-genetics/snv_indels/commit/883b5797e4b5ad1785344647d84dd57657247f59))
* update path, mutect2_gvcf to gatk_mutect2_gvcf ([111cd4e](https://www.github.com/hydra-genetics/snv_indels/commit/111cd4ebf7475031fd3d8c5cb4bb2911ad8e0895))
* wrong parameter name ([3630413](https://www.github.com/hydra-genetics/snv_indels/commit/363041330235702741e34d2d9c747b05f980c008))


### Documentation

* new rulegraph and updated readme ([58bcf0d](https://www.github.com/hydra-genetics/snv_indels/commit/58bcf0dd414677dee9bcdf9ecd42681430277e93))
* update compatibility list ([53cf137](https://www.github.com/hydra-genetics/snv_indels/commit/53cf13718e6e030f37cd1d05a3360779dd7bb808))

### [0.1.1](https://www.github.com/hydra-genetics/snv_indels/compare/v0.1.0...v0.1.1) (2022-04-21)


### Documentation

* update README to include haplotypecaller ([14e1367](https://www.github.com/hydra-genetics/snv_indels/commit/14e136766d0ceee89eb2e5a6191e004fc0bca21f))

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
