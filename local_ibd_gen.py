import numpy as np

class LocalIBDGenerator(object):
    def __init__(self,cluster_distro):
        self.cluster_distro = cluster_distro    
        self.create_graphs()
    def create_graphs(self):
        inds = np.sum(self.cluster_distro)
        cls_count = self.cluster_distro.shape[0]
        self.membership_array = np.zeros((inds))
        index = 0
        sample_counter = 0
        self.cls_dict = []
        for graph_size in self.cluster_distro:
            temp_array = set()
            for i in range(graph_size):
                temp_array.add(sample_counter)
                self.membership_array[sample_counter] = index
                sample_counter += 1
            self.cls_dict.append(temp_array)
            index += 1
        self.total_count = inds
        return self.membership_array, self.cls_dict
    def generate_file(self, fp_rate, fn_rate,path):
        with open(path+'.abc','w') as output_file:
            people_count = self.membership_array.shape[0]
        
            counter = 0
            for cls_set in self.cls_dict:
                cls_list = list(cls_set)
                for i in range(len(cls_list)):
                    for j in range(i + 1, len(cls_list)):
                        if not bernoulli.rvs(fn_rate):
                            output_file.write('\t'.join([str(item) for item in [cls_list[i],cls_list[j]]])+'\n')
                        
                            counter += 1
            counter = int(fp_rate * counter)
            while counter > 0:
                chosen_ones = np.random.choice(people_count, 2, replace=False)
                if self.membership_array[chosen_ones[0]] != self.membership_array[chosen_ones[1]]:
                    output_file.write('\t'.join([str(item) for item in [chosen_ones[0],chosen_ones[1]]])+'\n')
                    counter -= 1
    def create_phenotype(self,general_mean,general_std):
        phenotype_array = np.zeros(self.membership_array.shape[0],dtype=object)
        for i in range(self.membership_array.shape[0]):
            phenotype_array[i] = []
        prob_array = np.zeros(self.membership_array.shape[0],dtype=float)
        self.pheno_preval = []
        pheno_count = 0 
        for i in range(len(self.cls_dict)):
            cls_listed = list(self.cls_dict[i])
            if len(self.cls_dict[i]) >= THRESHOLD:
                #temp_pheno_array = np.zeros(self.membership_array.shape[0],dtype=bool)
                #temp_pheno_array[:] = False
                preval = np.abs(norm.rvs(loc=general_mean,scale=general_std))
                general_count = self.membership_array.shape[0]-len(self.cls_dict[i])
                general_affected_count = binom.rvs(general_count,preval)
                prob_array[:] = 1
                prob_array[cls_listed] = 0
                prob_array = prob_array/general_count
                general_index_list = np.random.choice(self.total_count,general_affected_count,replace=False,p=prob_array)
                for item in general_index_list:
                    phenotype_array[item].append(pheno_count)
                
                affected_count = int(binom.isf(PVAL,len(self.cls_dict[i]),preval))
                index_list = np.random.choice(len(self.cls_dict[i]),affected_count,replace=False)
                for index in index_list:
                    phenotype_array[cls_listed[index]].append(pheno_count)
                    
                pheno_count += 1
                self.pheno_preval.append(preval)
                
        
        
        
        
        self.phenotype_array = phenotype_array
        
        return self.phenotype_array,self.pheno_preval
            

    def write_to_file(self,path):
        result = {'total_count':self.total_count,'cls_dict':self.cls_dict,
                    'membership_array':self.membership_array,'pheno_preval':self.pheno_preval,
                    'phenotyep_array':self.phenotype_array}
        pkl.dump(result,open(path+'pkl','wb'))
            
                
