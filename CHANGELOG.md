# Changelog

## [1.3.1](https://github.com/hydra-genetics/snv_indels/compare/v1.3.0...v1.3.1) (2025-11-06)


### Bug Fixes

* **mosaicforecast:** Update path parameter to use input directory ([e10fd45](https://github.com/hydra-genetics/snv_indels/commit/e10fd45df3001c5f7675437a6fb144923ee4bde5))

## [1.3.0](https://github.com/hydra-genetics/snv_indels/compare/v1.2.0...v1.3.0) (2025-09-24)


### Features

* add functions for compiling paths to input bam files & output files for whatshap rules ([e8d96ef](https://github.com/hydra-genetics/snv_indels/commit/e8d96ef2d6159523ecb10f68582d2b319769d2c3))


### Bug Fixes

* add umap to config_mosaic.yaml ([2db78f5](https://github.com/hydra-genetics/snv_indels/commit/2db78f5d9336070b55d0b363f526a8dbcfd9a0be))
* change callers name to clairs_to in fix_af.py ([e68a3af](https://github.com/hydra-genetics/snv_indels/commit/e68a3af69cbb57ae9f48164fa2e14d593191149e))
* correct alignment path in get_input_bam() ([e2c36ee](https://github.com/hydra-genetics/snv_indels/commit/e2c36ee9fd4e83f04e6b4c05472e26973c69e5b0))
* flexible bam input to deepsomatic_t_only ([579a52a](https://github.com/hydra-genetics/snv_indels/commit/579a52aecacfbff4978663ec74a9f9535dbfb03e))
* input.bam in message not input.aln ([262b227](https://github.com/hydra-genetics/snv_indels/commit/262b2275cf89f14461d619071c50ae9e129bac89))
* remove redundant line ([d41441e](https://github.com/hydra-genetics/snv_indels/commit/d41441e36de8c17f58c2b863e0fe2b25f689e47c))
* update fix_af.py ([08d33ba](https://github.com/hydra-genetics/snv_indels/commit/08d33ba9301b14115445c2f26b41cb584ca3a2bc))
* use aln not bam in input to whatshap_haplotag ([9aac4c8](https://github.com/hydra-genetics/snv_indels/commit/9aac4c8b484249200eb7f6f78cc7d5a14cc09f5e))


### Documentation

* add whatshap rules to DAG ([124bd59](https://github.com/hydra-genetics/snv_indels/commit/124bd598fa6e9fd1190105f46daefd43d76700e7))
* list all inputs alogn with their names in whatshap rules ([fdd67da](https://github.com/hydra-genetics/snv_indels/commit/fdd67dab77a3dd63c5c28ac500f62bbd357ab6bf))
* update pull_request_template.md ([71d82e1](https://github.com/hydra-genetics/snv_indels/commit/71d82e1cef1960d958a3d38bdd972f0950b9f625))
* update whatshap_haplotag in rules.schema.yaml ([9d20025](https://github.com/hydra-genetics/snv_indels/commit/9d20025efdca62fae71dfdc24fbba401269191b3))

## [1.2.0](https://www.github.com/hydra-genetics/snv_indels/compare/v1.1.0...v1.2.0) (2025-08-29)


### Features

* add bam with type in SM field ([35e0a8a](https://www.github.com/hydra-genetics/snv_indels/commit/35e0a8af21e12e4ed6bdfd858a6c5d017dc7893d))
* add ClairS-TO ([2c35701](https://www.github.com/hydra-genetics/snv_indels/commit/2c35701e73bc6d5267701d16a5f140080aaf334f))
* add DeepSomatic ([ec6e9b7](https://www.github.com/hydra-genetics/snv_indels/commit/ec6e9b735b916933ff8a46cbaecd06450ed7194d))
* add hiphase ([c8f45a6](https://www.github.com/hydra-genetics/snv_indels/commit/c8f45a677bd7ef3dcc8f0899ecdcdee0dc0e827f))
* add script for deepmosaic input and information around deepmosaic ([b1e9f0d](https://www.github.com/hydra-genetics/snv_indels/commit/b1e9f0d54b9d0268d6b1b0c2e5a662f6b7161311))
* create a pacbio deepvariant rule ([c358d0b](https://www.github.com/hydra-genetics/snv_indels/commit/c358d0b092fbb4b4ea2530c7be226dc22e37e4ab))
* create hydra-rule for deepmosaic ([0a85b35](https://www.github.com/hydra-genetics/snv_indels/commit/0a85b354173cb0f6e339e61a7399bebfe0248e74))
* mosaicforecast ([f504e4d](https://www.github.com/hydra-genetics/snv_indels/commit/f504e4da1bfa89a3a9751ba1895118013788683f))


### Bug Fixes

* add integration stuff and fix pycodestyle ([51318fb](https://www.github.com/hydra-genetics/snv_indels/commit/51318fbbc907ed24e1a6eda8d55a9c298038506c))
* build mkdocs github workflow ([71106e3](https://www.github.com/hydra-genetics/snv_indels/commit/71106e30b275bbf8b1a8d576a8bcb6cac8a3b432))
* build mkdocs github workflow ([e48c91e](https://www.github.com/hydra-genetics/snv_indels/commit/e48c91e7431c3c53b32793715ea2bc3f0ecf8743))
* build mkdocs github workflow ([367144d](https://www.github.com/hydra-genetics/snv_indels/commit/367144dea8f8a35b1a4b3ccf271cc8d9c5da7232))
* compliant name pattern ([5a1937c](https://www.github.com/hydra-genetics/snv_indels/commit/5a1937cbcad4d3eecdd31cad54c1ea9ec6bbbff8))
* compliant name pattern ([6108f49](https://www.github.com/hydra-genetics/snv_indels/commit/6108f499fd36334e773c9ee9b396fad582394553))
* credentials for pulling containers ([29665b2](https://www.github.com/hydra-genetics/snv_indels/commit/29665b291804efe7fc93c76689fd9541ff09b112))
* **deepvariant:** add missing extra params to shell ([d015a8e](https://www.github.com/hydra-genetics/snv_indels/commit/d015a8e890c3b8faca0cc85896da372804c58136))
* give info in pon.vcf ([1e68010](https://www.github.com/hydra-genetics/snv_indels/commit/1e680105d3da36f557b57e5ffa0896a27566bcc3))
* **hiphase:** add missing bai outputfile ([9ced9d9](https://www.github.com/hydra-genetics/snv_indels/commit/9ced9d9d73b685edbee98d6565cb2373e8f7ca88))
* integration and snakefmt ([dd1a7a5](https://www.github.com/hydra-genetics/snv_indels/commit/dd1a7a5fe1825982c1c79e869fcae85b715c8a00))
* linting ([faa4bcd](https://www.github.com/hydra-genetics/snv_indels/commit/faa4bcd2ee650e5040ff1acc7e049bc18a174821))
* linting ([3a86154](https://www.github.com/hydra-genetics/snv_indels/commit/3a86154f1e1ec6ab34795b1c8db7c6d525bc03c2))
* linting ([b9179f8](https://www.github.com/hydra-genetics/snv_indels/commit/b9179f8973c32e16e636c5f6d0f3d36a7bc7de82))
* make deepsomatic panel of normal postprocess more flexible ([2b4d1e8](https://www.github.com/hydra-genetics/snv_indels/commit/2b4d1e824fa937dc5193368b47f9c22fbe967338))
* phrasing ([6bd5cc7](https://www.github.com/hydra-genetics/snv_indels/commit/6bd5cc7e47af86f1995f1d73a9d7b11c8e30cf32))
* remove unecessary output ([b88a9b2](https://www.github.com/hydra-genetics/snv_indels/commit/b88a9b2d0b5e91965adc98aef615640aa690599b))
* remove unecessary output ([132ddac](https://www.github.com/hydra-genetics/snv_indels/commit/132ddac5582c1d477fdce5b5350e3ec45d9fedb2))
* remove unecessary output ([9b1ac9d](https://www.github.com/hydra-genetics/snv_indels/commit/9b1ac9dba0b9bc9471072bef1bff70dcc02ac29f))
* separate rule definition for deepvariant pacbio ([7736b3d](https://www.github.com/hydra-genetics/snv_indels/commit/7736b3d34e31e74350898b75a461f0a6455faa32))
* typo ([b5361aa](https://www.github.com/hydra-genetics/snv_indels/commit/b5361aa859609af939c3312610b5e9a2521fa332))
* typo ([d4cd5a9](https://www.github.com/hydra-genetics/snv_indels/commit/d4cd5a97f15d64ae624659ef2045f78b4315877c))
* wrong model used ([f3e4f2b](https://www.github.com/hydra-genetics/snv_indels/commit/f3e4f2ba3c113116b40788a1b73578a546ae855b))


### Documentation

* infos for deep-learning based models ([640df56](https://www.github.com/hydra-genetics/snv_indels/commit/640df5607fd5645daad86983f196fc49f2d1addc))

## [1.1.0](https://www.github.com/hydra-genetics/snv_indels/compare/v1.0.0...v1.1.0) (2024-10-30)


### Features

* Add option to merge af for complex variants ([b8129c3](https://www.github.com/hydra-genetics/snv_indels/commit/b8129c305099e9a7e909b4eebb48eefc0f64c86d))
* New rule to handle af from complex variants ([b881c57](https://www.github.com/hydra-genetics/snv_indels/commit/b881c57915f1754dd595535856fe767f292f59a3))


### Bug Fixes

* Added noqa-tag on import to pass pycodestyle ([b88b397](https://www.github.com/hydra-genetics/snv_indels/commit/b88b397c6663fe279e07a5984c5ef542b15c9485))
* Fixed snakefmt linting error ([e9ae6b8](https://www.github.com/hydra-genetics/snv_indels/commit/e9ae6b8a2fc92eeabcf1e76580fa5519bf14ac1b))
* Handling of invalid input parameter ([223785c](https://www.github.com/hydra-genetics/snv_indels/commit/223785cb4278fe4eee7d680dd4559de93da6ef88))
* KeyError handling and additional unittest ([1604d9c](https://www.github.com/hydra-genetics/snv_indels/commit/1604d9cf3b398755eb33d1a50ff15cbf24330e1d))
* pycodestyle-fix ([9a5cee6](https://www.github.com/hydra-genetics/snv_indels/commit/9a5cee6d8526404f7c75c2e86adfea1ee1ad095b))
* Removed deprecated Mambaforge added Miniforge ([629e1b5](https://www.github.com/hydra-genetics/snv_indels/commit/629e1b553d7d49d523700d0695f69a98025f37e3))
* Removed non used format- and info-keys ([3186f38](https://www.github.com/hydra-genetics/snv_indels/commit/3186f389116f155209b2af9fabb3670298789913))
* Rulegraph as png ([a7e6b48](https://www.github.com/hydra-genetics/snv_indels/commit/a7e6b48cb2fbe5806ebd62c31e79489d2351e6b5))
* Update config/config.yaml ([1c2cacf](https://www.github.com/hydra-genetics/snv_indels/commit/1c2cacfa9efacfb42be4e37fa2441b249e36d823))
* Update to Miniforge3 ([b21d8f1](https://www.github.com/hydra-genetics/snv_indels/commit/b21d8f1c5168cd009f21941ea52a5ed59ce805f1))
* Updated rulegraph with merge_af_complex_var. ([2419bff](https://www.github.com/hydra-genetics/snv_indels/commit/2419bfff6db13939ad4ed191729b8349a48295f9))
* Use apptainer instead of singularity ([fcfa01f](https://www.github.com/hydra-genetics/snv_indels/commit/fcfa01f7dc004721fa1e436e6918f4a221fcdb5f))

## [1.0.0](https://www.github.com/hydra-genetics/snv_indels/compare/v0.6.0...v1.0.0) (2024-04-18)


### ⚠ BREAKING CHANGES

* constrain bgzip and tabix to only work within module

### Features

* add intermediate results directory ([fe5226f](https://www.github.com/hydra-genetics/snv_indels/commit/fe5226f91f260807da479ce0383bbcc42761f560))
* added ruleorders ([be2d4f6](https://www.github.com/hydra-genetics/snv_indels/commit/be2d4f6d505ce29817eabb6890aafa4515764c82))
* constrain bgzip and tabix to only work within module ([e37c9ed](https://www.github.com/hydra-genetics/snv_indels/commit/e37c9ed2a3e4f9be4952588c66ada4a8bb98467f))
* use the run_deepvariant script ([1465ad0](https://www.github.com/hydra-genetics/snv_indels/commit/1465ad0b89d518a1494b1378eed8f5a812a21c75))


### Bug Fixes

* pulp ([6c11be8](https://www.github.com/hydra-genetics/snv_indels/commit/6c11be8f326f9445ee27df05defecd472148557c))
* remove extra pulp in requirements ([82477de](https://www.github.com/hydra-genetics/snv_indels/commit/82477deeeba01a5bd4d82715e1d29f3a54fb45cb))
* **tabix:** fix log filename ([a64e1fd](https://www.github.com/hydra-genetics/snv_indels/commit/a64e1fd42d0572e05e201aed3a814ee548fe013a))


### Documentation

* add deepvariant to softwares.md ([87988a0](https://www.github.com/hydra-genetics/snv_indels/commit/87988a0c32c5e217e94ac54633554659a6d36279))
* explain the optional gVCF output for deepvariant ([de20ce7](https://www.github.com/hydra-genetics/snv_indels/commit/de20ce7f7dc9675f0883c72304acf1f8c9bd5ed1))
* update mkdocs snakemake plugin ([298983b](https://www.github.com/hydra-genetics/snv_indels/commit/298983b3bdf733e1f6202c48b82b3b3987c65a21))

## [0.6.0](https://www.github.com/hydra-genetics/snv_indels/compare/v0.5.0...v0.6.0) (2023-10-24)


### Features

* **fix_af:** add deepvariant ([08489b8](https://www.github.com/hydra-genetics/snv_indels/commit/08489b80d705c289ebe1da3aee96321610face87))


### Bug Fixes

* sort shouldn't require a tbi file ([bee2ec0](https://www.github.com/hydra-genetics/snv_indels/commit/bee2ec0324b805d9f44f111fd60844502e9f1306))
* update schema with correct key value for gatk_mutect2_filter ([1f147fb](https://www.github.com/hydra-genetics/snv_indels/commit/1f147fb79ff7ba939e3bd59350f467018add4ca8))


### Documentation

* add file used to build docs at readthedocs ([0710507](https://www.github.com/hydra-genetics/snv_indels/commit/0710507da577871917e04778d3101b0e00b5d4d4))
* added all remaining rules ([4992b47](https://www.github.com/hydra-genetics/snv_indels/commit/4992b47fd7ee97f7c2590697d445a51eb7df1b15))
* added description of rule graph generation ([f077c95](https://www.github.com/hydra-genetics/snv_indels/commit/f077c958f0c9fc3548c5a9e8fb5b9cabea06d2b2))
* added documentation for rules up to gatk ([bb81083](https://www.github.com/hydra-genetics/snv_indels/commit/bb810835096421129fbad3dff5f39920a8c01990))
* added rule graph ([c5b04bb](https://www.github.com/hydra-genetics/snv_indels/commit/c5b04bb5c4aceb0d3442a0e5b31480f05565fe05))
* fix table in intro.md ([1d62e06](https://www.github.com/hydra-genetics/snv_indels/commit/1d62e06460160641bcdde9b1f17029b3f4da537c))
* include bugfixed rule ([1199c11](https://www.github.com/hydra-genetics/snv_indels/commit/1199c11f6c60afab28cf531202d52d4cc22935bf))
* update docs ([f2da493](https://www.github.com/hydra-genetics/snv_indels/commit/f2da49318544ccf763a7c0b5ad2aeda4c90fa90e))
* update mkdocs plugin version ([decd144](https://www.github.com/hydra-genetics/snv_indels/commit/decd144698f4f1031b0e95fb71a03db455cbd4e9))
* update plugin version ([8c2eee8](https://www.github.com/hydra-genetics/snv_indels/commit/8c2eee87e5d57696b250313203d473a6df6d9419))
* update plugin version ([ae5291c](https://www.github.com/hydra-genetics/snv_indels/commit/ae5291c2ec2698d11d3c8a5c59e441df7fc989af))
* update simplified rulegraph for readthedocs ([3e1d751](https://www.github.com/hydra-genetics/snv_indels/commit/3e1d7510a04fd0969aac9be1b9b36c5faca7cfad))

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
