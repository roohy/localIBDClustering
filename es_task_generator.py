import sys,os
import subprocess

PREFIX = 'exome_'
PYTHON_SCRIPT_ADDRESS = '/sc/arion/projects/ipm2/roohy/ukbb_ibd/ibdmapping/cls_stats/extract_stats.py'


if __name__ == '__main__':
    directory = sys.argv[1]
    startChr = int(sys.argv[2])
    endChr = int(sys.argv[3])
    subDir = sys.argv[4]
    outputAddr = sys.argv[5]
    for i in range(startChr,endChr+1):
        basePath = os.path.join(directory,str(i),'graphs')
        with open(os.path.join(basePath,'info_file')) as infoFile:
            for line in infoFile:
                data = line.strip().split()
                #data = [int(item) for item in data]
                subfileAddr = os.path.join(subDir,'subfile_'+str(i)+'_'+data[0])
                with open(subfileAddr,'w') as subfile:
                    subfile.write('#BSUB -J libd_cls_'+str(i)+'_'+data[0]+'\n')
                    subfile.write('#BSUB -P acc_ipm2\n#BSUB -q premium\n#BSUB -n 1\n'+
                         '#BSUB -R "span[hosts=1] affinity[core(2, same=socket, exclusive=(socket, injob))]"\n'+
                         '#BSUB -R rusage[mem=50000]\n#BSUB -W 0:30\n')
                    subfile.write('#BSUB -o '+os.path.join(subDir,str(i)+'_'+data[0]+'.stdout')+'\n')
                    subfile.write('#BSUB -e '+os.path.join(subDir,str(i)+'_'+data[0]+'.stderr')+'\n')
                    subfile.write('#BSUB -L /bin/bash\n')
                    subfile.write('module load anaconda3\n')
                    subfile.write('source activate /sc/arion/projects/ipm2/roohy/conda/envs/ukbb\n')
                    args = ['python', PYTHON_SCRIPT_ADDRESS,os.path.join(basePath,data[0],'_abc'),
                        os.path.join(basePath,data[0],'_abc.cls'),
                        os.path.join(outputAddr,str(i)+'_'+data[0])]
                    subfile.write(' '.join(args))
                subprocess.run('bsub < '+subfileAddr,shell=True)
                print(subfileAddr)


