import numpy as np
import real_data_analyzer

distro = np.load('distro.npy')

cluster_count = 5000
fp_rate = 0.05
fn_rate = 0.25
preval_mean = 0.01
preval_std = 0.01
cls_size_list = np.random.choice(np.arange(1, 727), cluster_count, p=distro)
ibd_generator = real_data_analyzer.LocalIBDGenerator(cls_size_list)
ibd_generator.generate_file(fp_rate,fn_rate,'graph_file_address')
ibd_generator.create_phenotype(preval_mean,preval_std)
ibd_generator.write_to_file('pheno_file_address')
analyzer = real_data_analyzer.LocalIBDGraph('graph_file_address',simulated=True,
                                 pheno_path='pheno_file_address')
analyzer.partition()
stats = analyzer.calculate_stats()