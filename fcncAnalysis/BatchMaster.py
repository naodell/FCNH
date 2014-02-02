import sys, os, subprocess, fileinput, math, tempfile, datetime


def get_current_time():

    now = datetime.datetime.now()
    currentTime = '{0:02d}{1:02d}{2:02d}_{3:02d}{4:02d}{5:02d}'.format(now.year, now.month, now.day, now.hour, now.minute, now.second)
    return currentTime


def make_directory(filePath, clear = True):

    if not os.path.exists(filePath):
        os.system('mkdir -p '+filePath)

    if clear and len(os.listdir(filePath)) != 0:
        os.system('rm '+filePath+'/*')


class JobConfig():
    '''Class for storing configuration for each dataset'''
    def __init__(self, dataName = 'TEST', inDir = '/uscms/home/naodell/nobackup/TEST', nJobs = 1, arguments = 'TEST 0 muon'):
        self._dataName  = dataName
        self._inDir     = inDir
        self._nJobs     = nJobs
        self._args      = arguments


class BatchMaster():
    '''A tool for submitting batch jobs'''
    def __init__(self, configList, shortQueue = False, stageDir = 'outputBatch', executable = 'execBatch.csh', selection = 'test'):
        self._selection  = selection
        self._current    = os.path.abspath('.')
        self._stageDir   = stageDir
        self._configList = configList
        self._executable = executable
        self._shortQueue = shortQueue
    

    def split_jobs_by_dataset(self, directory, nJobs):

        fileList = os.listdir(directory)
        nFiles = len(fileList)
        
        # Split files to requested number.  Cannot exceed
        # the number of files being run over.
        if nJobs > nFiles:
            nJobs = nFiles

        nFilesPerJob = int(math.ceil(float(nFiles)/float(nJobs)))
        fileSplit = [fileList[i:i+nFilesPerJob] for i in range(0, len(fileList), nFilesPerJob)]


        return fileSplit


    def make_batch_lpc(self, cfg, sources):
        '''
        Prepares for submission to lpc.  Does the following:

        1. Generates input_files.txt with files to run over
        2. Write batch configuration file
        '''

        ## Writing the batch config file
        batch_tmp = open('.batch_tmp_{0}'.format(cfg._dataName,), 'w')
        batch_tmp.write('Universe              = vanilla\n')
        batch_tmp.write('Should_Transfer_Files = YES\n')
        batch_tmp.write('WhenToTransferOutput  = ON_EXIT\n')
        batch_tmp.write('Notification          = Never\n')
        batch_tmp.write('Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000\n')
        batch_tmp.write('notify_user           = brian.lee.pollack@cern.ch\n')
        if self._shortQueue:
            batch_tmp.write('+LENGTH               = "SHORT"\n')
        batch_tmp.write('\n')

        for i, source in enumerate(sources):

            ## make file with list of inputs ntuples for the analyzer
            input = open('input_{1}_{2}.txt'.format(self._stageDir, cfg._dataName, str(i+1)), 'w')

            path = cfg._inDir
            if cfg._inDir[:5] == '/pnfs':
                path = 'dcap://cmsdca1.fnal.gov:24140/pnfs/fnal.gov/usr' + cfg._inDir[5:]

            for file in source:
                input.write(path+'/'+file+'\n')
            input.close()

            batch_tmp.write('Arguments             = {0} {1} {2}\n'.format(cfg._dataName, i+1, cfg._args))
            batch_tmp.write('Executable            = {0}\n'.format(self._executable))
            batch_tmp.write('Transfer_Input_Files  = source.tar.gz, input_{0}_{1}.txt\n'.format(cfg._dataName, i+1))
            batch_tmp.write('Output                = reports/{0}_{1}_$(Cluster)_$(Process).stdout\n'.format(cfg._dataName, i+1))
            batch_tmp.write('Error                 = reports/{0}_{1}_$(Cluster)_$(Process).stderr\n'.format(cfg._dataName, i+1))
            batch_tmp.write('Log                   = reports/{0}_{1}_$(Cluster)_$(Process).log   \n'.format(cfg._dataName, i+1))
            batch_tmp.write('Queue\n\n')

        batch_tmp.close()
        

    def submit_to_batch(self, bSystem = 'lpc'):
        '''
        Submits batch jobs to batch.  Currently only works
        for lpc batch system, but should be updated for more 
        general use
        '''

        self._stageDir  = self._stageDir + '/' + get_current_time()
        make_directory(self._stageDir, clear=False)

        print 'Creating tarball of current workspace...'
        os.system('tar czf {0}/source.tar.gz ../../Analysis_CMS 2> /dev/null'.format(self._stageDir))
        print 'Done!'

        subprocess.call('cp {0} {1}'.format(self._executable, self._stageDir), shell=True)
        os.chdir(self._stageDir)
        make_directory('reports', clear=False)
        
        if bSystem is 'lpc':

            for cfg in self._configList:
                sourceFiles = self.split_jobs_by_dataset(cfg._inDir, cfg._nJobs)

                self.make_batch_lpc(cfg, sourceFiles)
                subprocess.call('condor_submit .batch_tmp_{0}'.format(cfg._dataName), shell=True)


